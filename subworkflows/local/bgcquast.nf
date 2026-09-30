/*
    Per-mode bgc-quast inputs and runs: compare-tools, compare-samples,
    compare-to-reference; QUAST in compare-to-reference only, unless
    --bgc_quast_quastdir supplies its output.
*/

include { QUAST    } from '../../modules/nf-core/quast/main'
include { BGCQUAST } from '../../modules/local/bgcquast'

workflow BGCQUAST_COMPARISON {
    take:
    // DeepBGC and GECCO query channels omit samples with no BGCs.
    antismash_json     // [ meta, json ]     query antiSMASH results
    deepbgc_tsv        // [ meta, tsv ]      query DeepBGC results
    gecco_clusters     // [ meta, tsv ]      query GECCO results
    genomes            // [ meta, fasta ]    query contigs, QUAST input
    genome_gbks        // [ meta, gbk ]      query annotation, --genome
    ref_antismash_json // [ meta, json ]     reference antiSMASH results
    ref_deepbgc_tsv    // [ meta, tsv ]      reference DeepBGC results
    ref_gecco_clusters // [ meta, tsv ]      reference GECCO results
    ref_genome         // [ meta, fasta ]    reference contigs, QUAST input
    ref_genome_gbk     // [ meta, gbk ]      reference annotation
    ref_name           // val                reference name, --ref-name
    antismash_gbk      // [ meta, [ gbk ] ]  query antiSMASH region GBKs
    bigscape_dir       // [ [ tool: dir ] ]  BiG-SCAPE folders, or [ [:] ]
    ref_antismash_gbk  // [ meta, [ gbk ] ]  reference antiSMASH region GBKs

    main:
    ch_versions    = Channel.empty()
    ch_bgcquast_in = Channel.empty()
    def mode       = params.bgc_quast_mode

    def proper = [antismash: 'antiSMASH', deepbgc: 'DeepBGC', gecco: 'GECCO']

    // Log colours, dropped when --monochrome_logs is set.
    def orange      = params.monochrome_logs ? '' : "\033[38;5;208m"
    def orange_bold = params.monochrome_logs ? '' : "\033[1;38;5;208m"
    def pink        = params.monochrome_logs ? '' : "\033[1;38;5;197m"
    def red         = params.monochrome_logs ? '' : "\033[1;31m"
    def yellow      = params.monochrome_logs ? '' : "\033[1;93m"
    def white       = params.monochrome_logs ? '' : "\033[97m"
    def creset      = params.monochrome_logs ? '' : "\033[0m"
    def banner      = "=".multiply(100)

    // antiSMASH writes a JSON even with zero regions; GECCO and DeepBGC
    // write no file at all.
    def ch_antismash_json = antismash_json.join(antismash_gbk).map { meta, json, _gbks -> [meta, json] }
    def ch_ref_antismash_json = ref_antismash_json.join(ref_antismash_gbk).map { meta, json, _gbks -> [meta, json] }

    // A skipped tool never ran, so it must not be reported.
    def active = []
    if (!params.bgc_skip_antismash) active << 'antiSMASH'
    if (!params.bgc_skip_deepbgc)   active << 'DeepBGC'
    if (!params.bgc_skip_gecco)     active << 'GECCO'

    def ch_found_ids = ch_antismash_json.map { meta, _f -> ['antiSMASH', meta.id] }
        .mix(deepbgc_tsv.map    { meta, _f -> ['DeepBGC', meta.id] })
        .mix(gecco_clusters.map { meta, _f -> ['GECCO', meta.id] })
        .toList()
        .map { rows -> [rows] }

    // Samples with no BGC from any active tool, listed again at run end.
    def no_bgc_notes = []

    genomes.map { meta, _g -> meta.id }
        .toSortedList()
        .map { ids -> [ids] }
        .combine(ch_found_ids)
        .subscribe { ids, rows ->
            def have = rows.groupBy { it[0] }.collectEntries { t, v -> [(t): v.collect { it[1] } as Set] }
            def tools_missing_for = [:]
            active.each { t ->
                ids.findAll { !(have[t] ?: [] as Set).contains(it) }
                    .each { id -> tools_missing_for.get(id, []) << t }
            }
            tools_missing_for.each { id, tools ->
                def what = tools.size() > 1 ? 'those reports' : 'that report'
                no_bgc_notes << [
                    "'${id}': no BGCs predicted, so ${tools.join(', ')} " +
                        "produced no result;",
                    "this sample has no column in ${what}.",
                ]
            }
            active.each { t ->
                if (ids.every { !(have[t] ?: [] as Set).contains(it) }) {
                    no_bgc_notes << [
                        "${t} found no BGCs in any sample, so no ${t} " +
                            "report was produced.",
                    ]
                }
            }
        }

    def ch_ref_found = ch_ref_antismash_json.map { _m, _f -> 'antiSMASH' }
        .mix(ref_deepbgc_tsv.map    { _m, _f -> 'DeepBGC' })
        .mix(ref_gecco_clusters.map { _m, _f -> 'GECCO' })
        .toList()
        .map { found -> [found] }

    ref_name.combine(ch_ref_found).subscribe { rid, found ->
        def tools = active.findAll { !found.contains(it) }
        if (tools) {
            def what = tools.size() > 1 ?
                'those reports were' : 'that report was'
            no_bgc_notes << [
                "reference '${rid}': no BGCs predicted, so " +
                    "${tools.join(', ')} produced no result;",
                "${what} not produced.",
            ]
        }
    }
    def ref_id = null

    ref_name.subscribe { rid -> ref_id = rid }

    def run_info = workflow   

    // workflow is null inside onComplete, so it is held here.
    workflow.onComplete {
        if (no_bgc_notes && run_info.success) {
            def pad   = ' '.multiply(16)
            def lines = []
            no_bgc_notes.eachWithIndex { note, i ->
                if (i > 0) { lines << '&' }
                lines.addAll(note)
            }
            println("${white}${banner}${creset}\n" +
                "${orange_bold}[bgc_quast_ppl] " +
                lines.join('\n' + pad) + "${creset}\n" +
                "${white}${banner}${creset}")
        }
    }

    if (mode == 'compare-tools') {
        // One bgc-quast run per sample.
        def tool_order = ['antismash', 'deepbgc', 'gecco']

        ch_bgcquast_in = ch_antismash_json.map { meta, f -> [meta, 'antismash', f] }
            .mix(deepbgc_tsv.map    { meta, f -> [meta, 'deepbgc', f] })
            .mix(gecco_clusters.map { meta, f -> [meta, 'gecco', f] })
            .groupTuple(by: 0)
            .map { meta, tools, files ->
                def idx           = (0..<tools.size()).toList().sort { tool_order.indexOf(tools[it]) }
                def ordered_files = idx.collect { files[it] }
                [meta, ordered_files]
            }
            .join(genome_gbks)
            .map { meta, files, genome ->
                // No --names: bgc-quast labels columns by the tool it detects.
                [meta + [leaf: "${meta.id}"], files, genome, [], [], [], []]
            }
    }
    else if (mode == 'compare-samples') {
        // One bgc-quast run per tool.
        def by_tool = { ch, tool ->
            ch.join(genome_gbks).map { meta, f, g -> [tool, meta.id, f, g] }
        }

        ch_bgcquast_in = by_tool(ch_antismash_json, 'antismash')
            .mix(by_tool(deepbgc_tsv, 'deepbgc'))
            .mix(by_tool(gecco_clusters, 'gecco'))
            .groupTuple(by: 0)
            // Lists sorted by sample id, so every report has the same column
            // order; groupTuple keeps arrival order, which differs per run.
            .map { tool, ids, files, gens ->
                def idx = (0..<ids.size()).toList().sort { ids[it] }
                [tool, idx.collect { ids[it] }, idx.collect { files[it] }, idx.collect { gens[it] }]
            }
            .combine(bigscape_dir)
            .map { tool, ids, files, gens, bsmap ->
                [
                    [id: "compare_samples_${tool}", bgcquast_names: ids.join(','), leaf: proper[tool]],
                    files, gens, [], [], [],
                    bsmap[tool] ?: [],
                ]
            }
    }
    else if (mode == 'compare-to-reference') {
        // One reference genome: contigs for QUAST, GenBank for bgc-quast.
        ch_ref_genome_file = ref_genome.map { meta, g -> g }.first()
        ch_ref_genome_gbk  = ref_genome_gbk.map { meta, g -> g }.first()

        ch_query_ordered = genomes
            .map { meta, g -> [meta.id, g] }
            .toSortedList { a, b -> a[0] <=> b[0] }
            .map { rows -> [rows.collect { it[0] }, rows.collect { it[1] }] }

        // One QUAST run of all queries against the reference, unless a folder
        // is supplied.
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

        // One run per tool: sorted query predictions, the reference prediction
        // and genome, and the QUAST folder.
        def per_tool_ref = { qch, rch, tool ->
            qch.join(genome_gbks)
                .map { meta, qfile, genome -> [meta.id, qfile, genome] }
                .toSortedList { a, b -> a[0] <=> b[0] }
                .filter { rows -> rows.size() > 0 }
                .map { rows ->
                    [rows.collect { it[0] }, rows.collect { it[1] }, rows.collect { it[2] }]
                }
                .combine(rch.map { meta, f -> f })
                .combine(ch_ref_genome_gbk)
                .combine(ch_quast_dir)
                .combine(ref_name)
                .map { names, files, gens, rfile, rgen, qdir, rid ->
                    [
                        [id: "compare_to_reference_${tool}", bgcquast_names: names.join(','), ref_name: rid, leaf: proper[tool]],
                        files, gens, qdir, rfile, rgen, [],
                    ]
                }
        }

        ch_bgcquast_in = per_tool_ref(ch_antismash_json, ch_ref_antismash_json, 'antismash')
            .mix(per_tool_ref(deepbgc_tsv,    ref_deepbgc_tsv,    'deepbgc'))
            .mix(per_tool_ref(gecco_clusters, ref_gecco_clusters, 'gecco'))

        // Empty when neither the reference nor any query has a usable
        // prediction.
        ch_bgcquast_in = ch_bgcquast_in.ifEmpty {
            error(
                "\n${white}${banner}${creset}\n" +
                "${red}[bgc_quast_ppl] The reference '${ref_id}' has no " +
                "predicted BGCs,\n" +
                "                so bgc-quast did not run.\n" +
                "                Check whether ${active.join(', ')} " +
                "predicted any BGC\n" +
                "                in that genome. The same message " +
                "appears if no query\n" +
                "                sample has a predicted BGC " +
                "either.${creset}\n" +
                "${white}${banner}${creset}"
            )
        }
    }
    else {
        error(
            "\n${white}${banner}${creset}\n" +
            "${pink}[bgc_quast_ppl] --bgc_quast_mode '${mode}' " +
            "is not supported.${creset}\n" +
            "${pink}                Please use compare-tools, " +
            "compare-samples or compare-to-reference.${creset}\n" +
            "${white}${banner}${creset}"
        )
    }

    BGCQUAST(ch_bgcquast_in)
    ch_versions = ch_versions.mix(BGCQUAST.out.versions)

    emit:
    results  = BGCQUAST.out.results // [ meta, files ]
    tsv      = BGCQUAST.out.tsv     // [ meta, report.tsv ]
    versions = ch_versions          // [ path(versions.yml) ]
}
