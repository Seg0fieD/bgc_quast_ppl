/*
    Main workflow: input preparation, contig length filter, annotation, BGC
    prediction, optional BiG-SCAPE and the bgc-quast comparison.
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { ANNOTATION          } from '../subworkflows/local/annotation'
include { BGC_PREDICTION      } from '../subworkflows/local/bgc_prediction'
include { BGCQUAST_COMPARISON } from '../subworkflows/local/bgcquast'
include { BIGSCAPE_ANALYSIS   } from '../subworkflows/local/bigscape_analysis'

include { GUNZIP as GUNZIP_INPUT_PREP     } from '../modules/nf-core/gunzip/main'
include { SEQKIT_SEQ as SEQKIT_SEQ_LENGTH } from '../modules/nf-core/seqkit/seq/main'

// Main workflow
workflow BGC_QUAST_PPL {
    take:
    ch_samplesheet // channel: [ meta, fasta, faa, gbk ] from --input

    main:

    ch_versions = Channel.empty()
    def yellow = params.monochrome_logs ? '' : "\033[1;93m"
    def orange = params.monochrome_logs ? '' : "\033[1;38;5;208m"
    def reset  = params.monochrome_logs ? '' : "\033[0m"
    def white  = params.monochrome_logs ? '' : "\033[97m"
    def banner = "=".multiply(100)

    // Printed at end of run; the live progress display erases log.info notes.
    def sample_warnings = []
    workflow.onComplete {
        sample_warnings.each { println it }
    }

    ch_bgcquast_run_count = Channel.value(0)

    if (params.bgc_quast_mode == 'compare-samples') {
        ch_samplesheet.count().subscribe { n ->
            if (n == 1) {
                log.info("${white}${banner}${reset}\n" +
                    "${orange}[bgc_quast_ppl] NOTE:    Running " +
                    "compare-samples mode with a single sample.${reset}\n" +
                    "${white}${banner}${reset}")
            }
        }
    }

    // The reference runs the same lane as the queries and is predicted once.
    ch_reference_rows = Channel.empty() // [ meta, [ ref_fasta, [], [] ] ]
    ch_ref_name       = Channel.empty() // reference name, --ref-name
    ch_query_samples  = ch_samplesheet  // all rows; reference removed below

    if (params.bgc_quast_mode == 'compare-to-reference') {
        // The samplesheet holds exactly one reference; checked at start-up.
        def ch_split = ch_samplesheet.branch { meta, fasta, faa, gbk ->
            reference: (meta.type ?: '').toLowerCase() == 'r'
            query: true
        }

        ch_query_samples = ch_split.query

        // Pin the single reference row to reuse it.
        def ch_ref  = ch_split.reference.first()

        ch_ref_name = ch_ref.map { meta, fasta, faa, gbk -> meta.id }

        ch_reference_rows = ch_ref.map { meta, fasta, faa, gbk ->
            [[id: meta.id, category: 'all', is_reference: true], [fasta, faa, gbk]]
        }

    }

    // One lane for queries and the reference; meta.is_reference splits the
    // prediction outputs back into query and reference channels.
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

    // Merge uncompressed and decompressed files into one input channel.
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

    // Length-filter contigs for BGC screening.
    if (params.run_bgc_screening) {
        SEQKIT_SEQ_LENGTH(ch_intermediate_input.fastas.map { meta, fasta, faa, gbk -> [meta, fasta] })
        ch_input_for_annotation = SEQKIT_SEQ_LENGTH.out.fastx
            .map { meta, fasta -> [meta + [category: 'long'], fasta] }
            .filter { meta, fasta ->
                if (fasta != [] && fasta.isEmpty()) {
                    def msg = "${white}${banner}${reset}\n" +
                        "${yellow}[bgc_quast_ppl] WARNING: Sample ${meta.id} " +
                        "has no contigs longer than " +
                        "${params.bgc_mincontiglength} bp.\n" +
                        "                         It will not be screened " +
                        "for BGCs.${reset}\n${white}${banner}${reset}"
                    sample_warnings << msg.toString()
                }
                !fasta.isEmpty()
            }
        ch_versions = ch_versions.mix(SEQKIT_SEQ_LENGTH.out.versions)
    }
    else {
        ch_input_for_annotation = ch_intermediate_input.fastas.map { meta, fasta, protein, gbk -> [meta, fasta] }
    }

    // Annotation
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

    // BGC Prediction and  bgc-quast
    if (params.run_bgc_screening) {
        BGC_PREDICTION(
            ch_prepped_input_long.fastas,
            ch_prepped_input_long.faas.filter { meta, file ->
                if (file != [] && file.isEmpty()) {
                    def msg = "${white}${banner}${reset}\n" +
                        "${yellow}[bgc_quast_ppl] WARNING: Annotation " +
                        "of sample ${meta.id} produced an empty FAA file.\n" +
                        "                         BGC tools that need it " +
                        "will be skipped.${reset}\n${white}${banner}${reset}"
                    sample_warnings << msg.toString()
                }
                !file.isEmpty()
            },
            ch_prepped_input_long.gbks.filter { meta, file ->
                if (file != [] && file.isEmpty()) {
                    def msg = "${white}${banner}${reset}\n" +
                        "${yellow}[bgc_quast_ppl] WARNING: Annotation " +
                        "of sample ${meta.id} produced an empty GBK file.\n" +
                        "                         BGC tools that need it " +
                        "will be skipped.${reset}\n${white}${banner}${reset}"
                    sample_warnings << msg.toString()
                }
                !file.isEmpty()
            },
        )
        ch_versions = ch_versions.mix(BGC_PREDICTION.out.versions)

        ch_pred_as     = BGC_PREDICTION.out.antismash_json.branch { meta, f -> reference: meta.is_reference; query: true }
        ch_pred_as_gbk = BGC_PREDICTION.out.antismash_gbk.branch  { meta, f -> reference: meta.is_reference; query: true }
        ch_pred_db     = BGC_PREDICTION.out.deepbgc_tsv.branch    { meta, f -> reference: meta.is_reference; query: true }
        ch_pred_ge     = BGC_PREDICTION.out.gecco_clusters.branch { meta, f -> reference: meta.is_reference; query: true }
        ch_long_fastas = ch_prepped_input_long.fastas.branch      { meta, f -> reference: meta.is_reference; query: true }
        ch_long_gbks   = ch_prepped_input_long.gbks.branch        { meta, f -> reference: meta.is_reference; query: true }

        // BiG-SCAPE runs on queries only, reference branch not used.
        ch_pred_ge_gbk = BGC_PREDICTION.out.gecco_gbk.branch   { meta, f -> reference: meta.is_reference; query: true }
        ch_pred_db_gbk = BGC_PREDICTION.out.deepbgc_gbk.branch { meta, f -> reference: meta.is_reference; query: true }

        ch_bigscape_dir = Channel.value([[:]])

        if (params.run_bigscape && params.bgc_quast_mode == 'compare-samples') {
            BIGSCAPE_ANALYSIS(
                ch_pred_as_gbk.query,
                ch_pred_ge_gbk.query,
                ch_pred_db_gbk.query,
                ch_pred_db.query,
            )
            ch_versions     = ch_versions.mix(BIGSCAPE_ANALYSIS.out.versions)
            ch_bigscape_dir = BIGSCAPE_ANALYSIS.out.dirs
        }

        BGCQUAST_COMPARISON(
            ch_pred_as.query,
            ch_pred_db.query,
            ch_pred_ge.query,
            ch_long_fastas.query,
            ch_long_gbks.query,
            ch_pred_as.reference,
            ch_pred_db.reference,
            ch_pred_ge.reference,
            ch_long_fastas.reference,
            ch_long_gbks.reference,
            ch_ref_name,
            ch_pred_as_gbk.query,
            ch_bigscape_dir,
            ch_pred_as_gbk.reference,
        )

        ch_versions = ch_versions.mix(BGCQUAST_COMPARISON.out.versions)
        ch_bgcquast_run_count = BGCQUAST_COMPARISON.out.results.count()
    }

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'bgc_quast_ppl_software_versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    emit:
    versions      = ch_versions           // channel: [ path(versions.yml) ]
    bgcquast_runs = ch_bgcquast_run_count // channel: bgc-quast run count
}

