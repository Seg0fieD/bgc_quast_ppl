/*
    Assemble per-mode inputs and run bgc-quast.
    Modes: compare-tools, compare-samples, compare-to-reference.
    QUAST runs only in compare-to-reference, unless --bgc_quast_quastdir is supplied.
*/

include { QUAST    } from '../../modules/nf-core/quast/main'
include { BGCQUAST } from '../../modules/local/bgcquast'
include { BIGSCAPE } from '../../modules/local/bigscape'

workflow BGCQUAST_COMPARISON {
    take:
    antismash_json     // [ meta, json ]  query
    deepbgc_tsv        // [ meta, tsv  ]  query (optional per sample)
    gecco_clusters     // [ meta, tsv  ]  query (optional per sample)
    genomes            // [ meta, fasta ] query contigs (--genome and QUAST consensus)
    ref_antismash_json // [ meta, json ]  reference, keyed by query id
    ref_deepbgc_tsv    // [ meta, tsv  ]  reference
    ref_gecco_clusters // [ meta, tsv  ]  reference
    ref_genome         // [ meta, fasta ] reference genome, keyed by query id
    ref_name           // val: reference display name (--ref-name)
    antismash_gbk      // [ meta, [ gbk ] ] query antiSMASH region GBKs (BiG-SCAPE input)
    
    main:
    ch_versions    = Channel.empty()
    ch_bgcquast_in = Channel.empty()
    def mode       = params.bgc_quast_mode

    def proper = [antismash: 'antiSMASH', deepbgc: 'DeepBGC', gecco: 'GECCO']

    /*
        BiG-SCAPE side branch. Runs once per pipeline execution over every sample's
        antiSMASH region GBKs together. Off unless --run_bigscape. Empty list means
        "no folder", which BGCQUAST reads as "do not pass --bigscape-output-dir".
    */
    ch_bigscape_dir = Channel.value([[]])

    if (params.run_bigscape) {
        if (params.bgc_bigscape_dir) {
            // User supplied a finished folder. Skip the run, same as bgc_quast_quastdir.
            ch_bigscape_dir = Channel.value(file(params.bgc_bigscape_dir, checkIfExists: true))
        }
        else {
            def pfam_hmm = file(params.bgc_bigscape_pfam, checkIfExists: true)

            // Pair "<sample_id>_<original_filename>" with its file, then sort so the two
            // lists the module receives stay index-aligned. The prefix is the join key
            // bgc-quast reverses; ".region" must survive for --include-gbk to accept it.
            ch_bigscape_stage = antismash_gbk
                .flatMap { meta, gbks ->
                    (gbks instanceof List ? gbks : [gbks]).collect { g ->
                        ["${meta.id}_${g.name}".toString(), g]
                    }
                }
                .toSortedList { a, b -> a[0] <=> b[0] }
                .filter { rows -> rows.size() > 0 }

            BIGSCAPE(
                ch_bigscape_stage.map { rows -> rows.collect { it[0] } },
                ch_bigscape_stage.map { rows -> rows.collect { it[1] } },
                Channel.value(pfam_hmm.parent),
                Channel.value(pfam_hmm.name),
            )
            ch_versions     = ch_versions.mix(BIGSCAPE.out.versions)
            ch_bigscape_dir = BIGSCAPE.out.results.ifEmpty { [[]] }
        }
    }

    if (mode == 'compare-tools') {
        // One run per sample
        def tool_order = ['antismash', 'deepbgc', 'gecco']

        ch_bgcquast_in = antismash_json.map { meta, f -> [meta, 'antismash', f] }
            .mix(deepbgc_tsv.map    { meta, f -> [meta, 'deepbgc', f] })
            .mix(gecco_clusters.map { meta, f -> [meta, 'gecco', f] })
            .groupTuple(by: 0)
            .map { meta, tools, files ->
                def idx           = (0..<tools.size()).toList().sort { tool_order.indexOf(tools[it]) }
                def ordered_files = idx.collect { files[it] }
                [meta, ordered_files]
            }
            .join(genomes)
            .map { meta, files, genome ->
                // No --names: bgc-quast auto-labels columns by detected tool.
                [meta + [leaf: "${meta.id}"], files, genome, [], [], [], []]
            }
    }
    else if (mode == 'compare-samples') {
        // One run per tool 
        def by_tool = { ch, tool ->
            ch.join(genomes).map { meta, f, g -> [tool, meta.id, f, g] }
        }

        ch_bgcquast_in = by_tool(antismash_json, 'antismash')
            .mix(by_tool(deepbgc_tsv, 'deepbgc'))
            .mix(by_tool(gecco_clusters, 'gecco'))
            .groupTuple(by: 0)
            .combine(ch_bigscape_dir)
            .map { tool, ids, files, gens, bsdir ->
                [
                    [id: "compare_samples_${tool}", bgcquast_names: ids.join(','), leaf: proper[tool]],
                    files, gens, [], [], [],
                    tool == 'antismash' ? bsdir : [],
                ]
            }
    }
    else if (mode == 'compare-to-reference') {
        // Single reference genome, reused by QUAST
        ch_ref_genome_file = ref_genome.map { meta, g -> g }.first()

        ch_query_ordered = genomes
            .map { meta, g -> [meta.id, g] }
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { rows -> [rows.collect { it[0] }, rows.collect { it[1] }] }

        // One QUAST run over all queries vs the reference, unless a dir is supplied.
        if (params.bgc_quast_quastdir) {
            ch_quast_dir = Channel.value(file(params.bgc_quast_quastdir, checkIfExists: true))
        }
        else {
            def ch_quast_in = ch_query_ordered.combine(ch_ref_genome_file)

            QUAST(
                ch_quast_in.map { ids, gs, _r -> [[id: 'quast', labels: ids.join(',')], gs] },
                ch_quast_in.map { ids, _gs, r -> [[id: 'quast'], r] },
                Channel.value([[id: 'quast'], []]),
            )
            ch_versions  = ch_versions.mix(QUAST.out.versions)
            ch_quast_dir = QUAST.out.results.map { meta, dir -> dir }.first()
        }

     // One run per tool: ordered query predictions + reference prediction, reference genome, and QUAST dir.
        def per_tool_ref = { qch, rch, tool ->
            qch.join(genomes)
                .map { meta, qfile, genome -> [meta.id, qfile, genome] }
                .toSortedList { a, b -> a[0] <=> b[0] }
                .filter { rows -> rows.size() > 0 }
                .map { rows ->
                    [rows.collect { it[0] }, rows.collect { it[1] }, rows.collect { it[2] }]
                }
                .combine(rch.map { meta, f -> f })
                .combine(ch_ref_genome_file)
                .combine(ch_quast_dir)
                .combine(ref_name)
                .map { names, files, gens, rfile, rgen, qdir, rid ->
                    [
                        [id: "compare_to_reference_${tool}", bgcquast_names: names.join(','), ref_name: rid, leaf: proper[tool]],
                        files, gens, qdir, rfile, rgen, [],
                    ]
                }
        }

        ch_bgcquast_in = per_tool_ref(antismash_json, ref_antismash_json, 'antismash')
            .mix(per_tool_ref(deepbgc_tsv,    ref_deepbgc_tsv,    'deepbgc'))
            .mix(per_tool_ref(gecco_clusters, ref_gecco_clusters, 'gecco'))

        // Empty means no predictions from the reference or any query.
        ch_bgcquast_in = ch_bgcquast_in.ifEmpty {
            error(
                "[bgc_quast_ppl] compare-to-reference produced no comparisons.\n" +
                "                The reference or all query samples yielded no usable BGC predictions,\n" +
                "                so QUAST and bgc-quast never ran.\n" +
                "                Check that the reference genome passes the contig-length filter, \n" +
                "                and its annotated, and produces antiSMASH/DeepBGC/GECCO output(s)."
            )
        }
    }
    else {
        // auto-mode infer from input
        error("[bgc_quast_ppl] bgc_quast_mode='${mode}' is not supported yet. Use compare-tools, compare-samples, or compare-to-reference.")
    }

    BGCQUAST(ch_bgcquast_in)
    ch_versions = ch_versions.mix(BGCQUAST.out.versions)

    emit:
    results  = BGCQUAST.out.results // [ meta, files ]
    tsv      = BGCQUAST.out.tsv     // [ meta, report.tsv ]
    versions = ch_versions          // [ path(versions.yml) ]
}