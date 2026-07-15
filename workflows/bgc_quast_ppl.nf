/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ANNOTATION          } from '../subworkflows/local/annotation'
include { BGC_PREDICTION      } from '../subworkflows/local/bgc_prediction'
include { BGCQUAST_COMPARISON } from '../subworkflows/local/bgcquast'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GUNZIP as GUNZIP_INPUT_PREP     } from '../modules/nf-core/gunzip/main'
include { SEQKIT_SEQ as SEQKIT_SEQ_LENGTH } from '../modules/nf-core/seqkit/seq/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BGC_QUAST_PPL {
    take:
    ch_samplesheet // channel: samplesheet read in from --input ; [ meta, fasta, faa, gbk ]

    main:

    ch_versions = Channel.empty()
    ch_bgcquast_run_count = Channel.value(0)

    // compare-samples: note when the run has just one sample (bgc-quast accepts n=1).
    if (params.bgc_quast_mode == 'compare-samples') {
        ch_samplesheet.count().subscribe { n ->
            if (n == 1) {
                log.info("[bgc_quast_ppl] Running compare-samples with a single sample.")
            }
        }
    }

    /*
        REFERENCE RESOLUTION (compare-to-reference only)
        There is exactly ONE reference per run, marked in the samplesheet with type r/R.
        It compares against one or many query genomes (type q/Q). (A different reference = a
        separate run.) The reference row's sample name becomes the bgc-quast --ref-name.
        The reference is a genome too, so it runs through the SAME prep -> annotation ->
        prediction lane as queries (tagged is_reference), predicted ONCE, then fanned out
        to every query.
    */
    ch_reference_rows = Channel.empty() // [ meta(is_reference), [ ref_fasta, [], [] ] ]
    ch_query_ref_link = Channel.empty() // [ rid, query_meta ] : links every query to the single reference
    ch_ref_name       = Channel.empty() // val: reference display name (bgc-quast --ref-name)
    ch_query_samples  = ch_samplesheet  // default: every row is a query

    if (params.bgc_quast_mode == 'compare-to-reference') {
        // Split the reference row (type r/R) from the query rows (type q/Q).
        // Samplesheet shape is already validated in PIPELINE_INITIALISATION.
        def ch_split = ch_samplesheet.branch { meta, fasta, faa, gbk ->
            reference: (meta.type ?: '').toLowerCase() == 'r'
            query: true
        }

        ch_query_samples = ch_split.query

        // Pin the single reference row to a value channel so it can be reused.
        def ch_ref  = ch_split.reference.first()

        // The reference's sample name is its display name (bgc-quast --ref-name).
        ch_ref_name = ch_ref.map { meta, fasta, faa, gbk -> meta.id }

        // one reference row -> predicted once, tagged is_reference
        ch_reference_rows = ch_ref.map { meta, fasta, faa, gbk ->
            [[id: meta.id, category: 'all', is_reference: true], [fasta, faa, gbk]]
        }

        // link every query to the reference id so the fan-out lines up per query
        ch_query_ref_link = ch_split.query
            .combine(ch_ref_name)
            .map { meta, fasta, faa, gbk, rid -> [rid, meta + [category: 'long', is_reference: false]] }
    }

    /*
        INPUT PREP
        Queries and (when present) the reference enter the same lane. is_reference tags them
        so the prediction outputs can be split again afterwards.
    */
    ch_query_rows = ch_query_samples
        .map { meta, fasta, faa, gbk -> [meta + [category: 'all', is_reference: false], [fasta, faa, gbk]] }

    ch_input_prep = ch_query_rows
        .mix(ch_reference_rows)
        .transpose()
        .branch {
            compressed: it[1].toString().endsWith('.gz')
            uncompressed: it[1]
        }

    GUNZIP_INPUT_PREP(ch_input_prep.compressed)
    ch_versions = ch_versions.mix(GUNZIP_INPUT_PREP.out.versions)

    // Merge uncompressed and newly-decompressed files into one input channel
    ch_intermediate_input = GUNZIP_INPUT_PREP.out.gunzip
        .mix(ch_input_prep.uncompressed)
        .groupTuple()
        .map { meta, files ->
            def fasta_found = files.find { it.toString().tokenize('.').last().matches('fasta|fas|fna|fa') }
            def faa_found = files.find { it.toString().endsWith('.faa') }
            def gbk_found = files.find { it.toString().tokenize('.').last().matches('gbk|gbff') }
            def fasta = fasta_found != null ? fasta_found : []
            def faa = faa_found != null ? faa_found : []
            def gbk = gbk_found != null ? gbk_found : []

            [meta, fasta, faa, gbk]
        }
        .branch { meta, fasta, faa, gbk ->
            preannotated: gbk != []
            fastas: true
        }

    // Duplicate and length-filter contigs for BGC screening (speeds up BGC, avoids 'no hits' fails)
    if (params.run_bgc_screening) {
        SEQKIT_SEQ_LENGTH(ch_intermediate_input.fastas.map { meta, fasta, faa, gbk -> [meta, fasta] })
        ch_input_for_annotation = ch_intermediate_input.fastas
            .map { meta, fasta, protein, gbk -> [meta, fasta] }
            .mix(SEQKIT_SEQ_LENGTH.out.fastx.map { meta, fasta -> [meta + [category: 'long'], fasta] })
            .filter { meta, fasta ->
                if (fasta != [] && fasta.isEmpty()) {
                    log.warn("[bgc_quast_ppl] Sample ${meta.id} has no contigs longer than ${params.bgc_mincontiglength} bp. Will not be screened for BGCs.")
                }
                !fasta.isEmpty()
            }
        ch_versions = ch_versions.mix(SEQKIT_SEQ_LENGTH.out.versions)
    }
    else {
        ch_input_for_annotation = ch_intermediate_input.fastas.map { meta, fasta, protein, gbk -> [meta, fasta] }
    }

    /*
        ANNOTATION
    */
    if (params.run_bgc_screening) {
        ANNOTATION(ch_input_for_annotation)
        ch_versions = ch_versions.mix(ANNOTATION.out.versions)

        ch_new_annotation = ch_input_for_annotation
            .join(ANNOTATION.out.faa)
            .join(ANNOTATION.out.gbk)
    }
    else {
        ch_new_annotation = ch_intermediate_input.fastas
    }

    // Mix preannotated samples back with the newly annotated ones
    ch_prepped_input = ch_new_annotation
        .filter { meta, fasta, faa, gbk -> meta.category != 'long' }
        .mix(ch_intermediate_input.preannotated)
        .multiMap { meta, fasta, faa, gbk ->
            fastas: [meta, fasta]
            faas: [meta, faa]
            gbks: [meta, gbk]
        }

    if (params.run_bgc_screening) {
        ch_prepped_input_long = ch_new_annotation
            .filter { meta, fasta, faa, gbk -> meta.category == 'long' }
            .mix(ch_intermediate_input.preannotated)
            .multiMap { meta, fasta, faa, gbk ->
                fastas: [meta, fasta]
                faas: [meta, faa]
                gbks: [meta, gbk]
            }
    }

    /*
        BGC PREDICTION + bgc-quast
    */
    if (params.run_bgc_screening) {
        BGC_PREDICTION(
            ch_prepped_input_long.fastas,
            ch_prepped_input_long.faas.filter { meta, file ->
                if (file != [] && file.isEmpty()) {
                    log.warn("[bgc_quast_ppl] Annotation of sample ${meta.id} produced an empty FAA file. BGC tools needing it will be skipped.")
                }
                !file.isEmpty()
            },
            ch_prepped_input_long.gbks.filter { meta, file ->
                if (file != [] && file.isEmpty()) {
                    log.warn("[bgc_quast_ppl] Annotation of sample ${meta.id} produced an empty GBK file. BGC tools needing it will be skipped.")
                }
                !file.isEmpty()
            },
        )
        ch_versions = ch_versions.mix(BGC_PREDICTION.out.versions)

        // Split prediction outputs and the prepped genomes into query vs reference lanes.
        ch_pred_as     = BGC_PREDICTION.out.antismash_json.branch { meta, f -> reference: meta.is_reference; query: true }
        ch_pred_db     = BGC_PREDICTION.out.deepbgc_tsv.branch    { meta, f -> reference: meta.is_reference; query: true }
        ch_pred_ge     = BGC_PREDICTION.out.gecco_clusters.branch { meta, f -> reference: meta.is_reference; query: true }
        ch_long_fastas = ch_prepped_input_long.fastas.branch      { meta, f -> reference: meta.is_reference; query: true }

        // Fan the single reference's results out to EVERY query, re-keyed to the query meta
        // so the joins in BGCQUAST_COMPARISON line up per query. Empty in non-reference modes.
        def fan_to_queries = { ref_ch ->
            ref_ch
                .map { meta, f -> [meta.id, f] }
                .combine(ch_query_ref_link, by: 0) // [ rid, f, query_meta ]
                .map { rid, f, qmeta -> [qmeta, f] }
        }

        BGCQUAST_COMPARISON(
            ch_pred_as.query,
            ch_pred_db.query,
            ch_pred_ge.query,
            ch_long_fastas.query,
            fan_to_queries(ch_pred_as.reference),
            fan_to_queries(ch_pred_db.reference),
            fan_to_queries(ch_pred_ge.reference),
            fan_to_queries(ch_long_fastas.reference),
            ch_ref_name,
        )
        ch_versions = ch_versions.mix(BGCQUAST_COMPARISON.out.versions)
        ch_bgcquast_run_count = BGCQUAST_COMPARISON.out.results.count()
    }

    //
    // Collate and save software versions
    //
    // NOTE: QUAST reports its version through a `topic: versions` channel, not .out.versions.
    // Collecting topic versions into this YAML is a separate nf-core wiring step (TODO).
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'bgc_quast_ppl_software_versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    emit:
    versions     = ch_versions                          // channel: [ path(versions.yml) ]
    bgcquast_runs = ch_bgcquast_run_count                // channel: val(Integer) number of bgc-quast runs
}
