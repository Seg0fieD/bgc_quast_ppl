/*
    Assemble per-mode inputs and run bgc-quast.
    Modes: compare-tools, compare-samples, compare-to-reference.
    QUAST runs only in compare-to-reference, unless --bgc_quast_quastdir is supplied.
*/

include { QUAST    } from '../../modules/nf-core/quast/main'
include { BGCQUAST } from '../../modules/local/bgcquast'
include { BIGSCAPE as BIGSCAPE_ANTISMASH } from '../../modules/local/bigscape'
include { BIGSCAPE as BIGSCAPE_GECCO     } from '../../modules/local/bigscape'
include { BIGSCAPE as BIGSCAPE_DEEPBGC   } from '../../modules/local/bigscape'
include { BIGSCAPE_DOWNLOAD_DB  } from '../../modules/local/bigscape_download_db'
include { DEEPBGC_SPLIT_GBK     } from '../../modules/local/deepbgc_split_gbk'

workflow BGCQUAST_COMPARISON {
    take:
    antismash_json     // [ meta, json ]  query
    deepbgc_tsv        // [ meta, tsv  ]  query (optional per sample)
    gecco_clusters     // [ meta, tsv  ]  query (optional per sample)
    genomes            // [ meta, fasta ] query contigs, QUAST input
    genome_gbks        // [ meta, gbk ]   query annotation (--genome)
    ref_antismash_json // [ meta, json ]  reference
    ref_deepbgc_tsv    // [ meta, tsv  ]  reference
    ref_gecco_clusters // [ meta, tsv  ]  reference
    ref_genome         // [ meta, fasta ] reference contigs, QUAST input
    ref_genome_gbk     // [ meta, gbk ]   reference annotation
    ref_name           // val: reference display name (--ref-name)
    antismash_gbk      // [ meta, [ gbk ] ] query antiSMASH region GBKs (BiG-SCAPE input)
    gecco_gbk          // [ meta, [ gbk ] ] query GECCO cluster GBKs (BiG-SCAPE input)
    deepbgc_gbk        // [ meta, gbk ]     query DeepBGC multi-record GBK (split first)
    ref_antismash_gbk  // [ meta, [ gbk ] ] reference antiSMASH region GBKs

    main:
    ch_versions    = Channel.empty()
    ch_bgcquast_in = Channel.empty()
    def mode       = params.bgc_quast_mode

    def proper = [antismash: 'antiSMASH', deepbgc: 'DeepBGC', gecco: 'GECCO']

    // Log colours, dropped when --monochrome_logs is set.
    def orange      = params.monochrome_logs ? '' : "\033[38;5;208m"
    def orange_bold = params.monochrome_logs ? '' : "\033[1;38;5;208m"
    def pink        = params.monochrome_logs ? '' : "\033[1;38;5;197m"
    def yellow      = params.monochrome_logs ? '' : "\033[1;93m"
    def hi          = params.monochrome_logs ? '' : "\033[4m"
    def noh         = params.monochrome_logs ? '' : "\033[24m"
    def white       = params.monochrome_logs ? '' : "\033[97m"
    def creset      = params.monochrome_logs ? '' : "\033[0m"

    // antiSMASH writes a JSON even with zero regions; GECCO and DeepBGC write no file at all.
    def ch_antismash_json = antismash_json.join(antismash_gbk).map { meta, json, _gbks -> [meta, json] }
    def ch_ref_antismash_json = ref_antismash_json.join(ref_antismash_gbk).map { meta, json, _gbks -> [meta, json] }

    // A skipped tool never ran, so it must not be reported as having found nothing.
    def active = []
    if (!params.bgc_skip_antismash) active << 'antiSMASH'
    if (!params.bgc_skip_deepbgc)   active << 'DeepBGC'
    if (!params.bgc_skip_gecco)     active << 'GECCO'

    def ch_found_ids = ch_antismash_json.map { meta, _f -> ['antiSMASH', meta.id] }
        .mix(deepbgc_tsv.map    { meta, _f -> ['DeepBGC', meta.id] })
        .mix(gecco_clusters.map { meta, _f -> ['GECCO', meta.id] })
        .toList()
        .map { rows -> [rows] }

    // Samples that no active tool could predict a BGC for, repeated when the run ends.
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
                no_bgc_notes << "${orange_bold}[bgc_quast_ppl] '${id}': no BGCs predicted, so ${tools.join(', ')} produced no result -- this sample has no column in ${tools.size() > 1 ? 'those reports' : 'that report'}.${creset}".toString()
            }
            active.each { t ->
                if (ids.every { !(have[t] ?: [] as Set).contains(it) }) {
                    no_bgc_notes << "${orange_bold}[bgc_quast_ppl] ${t} found no BGCs in any sample, so no ${t} report was produced.${creset}".toString()
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
            no_bgc_notes << "${orange_bold}[bgc_quast_ppl] reference '${rid}': no BGCs predicted, so ${tools.join(', ')} produced no result -- ${tools.size() > 1 ? 'those reports were' : 'that report was'} not produced.${creset}".toString()
        }
    }
    def ref_id = null

    ref_name.subscribe { rid -> ref_id = rid }

    workflow.onComplete {
        if (no_bgc_notes) {
            println ''
            no_bgc_notes.each { println it }
        }
    }

    /*
        BiG-SCAPE side branch. One run per tool over every sample's GBKs together.
        compare-samples only, and off unless --run_bigscape. The channel carries a
        [tool: dir] map; a tool missing from the map gets [], which BGCQUAST reads as
        "no --bigscape-output-dir". The map is always wrapped in a list, because
        combine() spreads one level.
    */
    def bs_tools = []
    if (!params.bgc_skip_antismash) { bs_tools << 'antismash' }
    if (!params.bgc_skip_gecco)     { bs_tools << 'gecco' }
    if (!params.bgc_skip_deepbgc)   { bs_tools << 'deepbgc' }

    ch_bigscape_dir = Channel.value([[:]])

    // The flag is never an error in another mode; the run just continues without it.
    if (params.run_bigscape && mode != 'compare-samples') {
        log.warn(
            "${yellow}[bgc_quast_ppl] BiG-SCAPE does not run in ${mode} mode, \n" +
            "       so BiG-SCAPE and its related steps are skipped. \n" +
            "       For the gene cluster family analysis, run the pipeline in \n" +
            "       ${hi}compare-samples${noh} mode with ${hi}--run_bigscape${noh}.${creset}"
        )
    }

    if (params.run_bigscape && mode == 'compare-samples') {
        // --bgc_bigscape_dir is a parent holding per-tool subfolders, mirroring the published
        // bgc_quast/bigscape/ layout. Tools without a subfolder still run normally.
        def given = [:]

        if (params.bgc_bigscape_dir) {
            file(params.bgc_bigscape_dir, checkIfExists: true)

            bs_tools.each { t ->
                def sub = file("${params.bgc_bigscape_dir}/${t}")
                if (sub.exists() && sub.isDirectory()) {
                    given[t] = sub
                }
            }

            if (!given) {
                error(
                    "[bgc_quast_ppl] --bgc_bigscape_dir contains no per-tool subfolder.\n" +
                    "                Expected at least one of: ${bs_tools.join(', ')}\n" +
                    "                Looked in: ${params.bgc_bigscape_dir}\n" +
                    "                Point it at a previous run's bgc_quast/bigscape/ folder."
                )
            }

            log.info("${orange}            [bgc_quast_ppl] BiG-SCAPE folder supplied for: ${given.keySet().join(', ')}${creset}")
        }

        def to_run = bs_tools.findAll { !given.containsKey(it) }

        if (to_run) {
            log.info("${orange}            [bgc_quast_ppl] BiG-SCAPE will run for: ${to_run.join(', ')}${creset}")
        }

        // Bare map on purpose. It is wrapped once, at the end, before combine() sees it.
        ch_bigscape_run = Channel.value([:])

        if (to_run) {
            // Pfam: use the supplied pressed copy, otherwise download and press one.
            // Resolved once and shared by every run.
            def ch_pfam_dir
            def ch_pfam_name

            if (params.bgc_bigscape_pfam) {
                def pfam_hmm = file(params.bgc_bigscape_pfam, checkIfExists: true)
                ch_pfam_dir  = Channel.value(pfam_hmm.parent)
                ch_pfam_name = Channel.value(pfam_hmm.name)
            }
            else {
                BIGSCAPE_DOWNLOAD_DB()
                ch_versions  = ch_versions.mix(BIGSCAPE_DOWNLOAD_DB.out.versions)
                ch_pfam_dir  = BIGSCAPE_DOWNLOAD_DB.out.db
                ch_pfam_name = Channel.value('Pfam-A.hmm')
            }

            // Stage each GBK as "<sample_id>_<original_filename>" and sort, so the two lists
            // the module receives stay index-aligned. bgc-quast reverses that prefix later,
            // and ".region" or "_cluster_" must survive for BiG-SCAPE's --include-gbk.
            def stage_gbks = { ch ->
                ch.flatMap { meta, gbks ->
                        (gbks instanceof List ? gbks : [gbks]).collect { g ->
                            ["${meta.id}_${g.name}".toString(), g]
                        }
                    }
                    .toSortedList { a, b -> a[0] <=> b[0] }
                    .filter { rows -> rows.size() > 0 }
            }

            ch_bigscape_results = Channel.empty()

            if ('antismash' in to_run) {
                def st = stage_gbks(antismash_gbk)
                BIGSCAPE_ANTISMASH(
                    'antismash',
                    st.map { rows -> rows.collect { it[0] } },
                    st.map { rows -> rows.collect { it[1] } },
                    ch_pfam_dir,
                    ch_pfam_name,
                )
                ch_versions         = ch_versions.mix(BIGSCAPE_ANTISMASH.out.versions)
                ch_bigscape_results = ch_bigscape_results.mix(BIGSCAPE_ANTISMASH.out.results)
            }

            if ('gecco' in to_run) {
                def st = stage_gbks(gecco_gbk)
                BIGSCAPE_GECCO(
                    'gecco',
                    st.map { rows -> rows.collect { it[0] } },
                    st.map { rows -> rows.collect { it[1] } },
                    ch_pfam_dir,
                    ch_pfam_name,
                )
                ch_versions         = ch_versions.mix(BIGSCAPE_GECCO.out.versions)
                ch_bigscape_results = ch_bigscape_results.mix(BIGSCAPE_GECCO.out.results)
            }

            if ('deepbgc' in to_run) {
                // The .bgc.tsv rides along because the BGC numbering comes from it, not the GBK.
                DEEPBGC_SPLIT_GBK(deepbgc_gbk.join(deepbgc_tsv, failOnDuplicate: true))
                ch_versions = ch_versions.mix(DEEPBGC_SPLIT_GBK.out.versions)

                def st = stage_gbks(DEEPBGC_SPLIT_GBK.out.gbk)
                BIGSCAPE_DEEPBGC(
                    'deepbgc',
                    st.map { rows -> rows.collect { it[0] } },
                    st.map { rows -> rows.collect { it[1] } },
                    ch_pfam_dir,
                    ch_pfam_name,
                )
                ch_versions         = ch_versions.mix(BIGSCAPE_DEEPBGC.out.versions)
                ch_bigscape_results = ch_bigscape_results.mix(BIGSCAPE_DEEPBGC.out.results)
            }

            // toList() always emits once, so a tool with an empty GBK channel simply leaves
            // its key out of the map.
            ch_bigscape_run = ch_bigscape_results
                .toList()
                .map { rows -> rows.collectEntries { t, d -> [(t): d] } }
        }

        // Supplied folders and freshly run ones cannot overlap: to_run excludes the supplied.
        // Wrapped in a list here, once, because combine() spreads one level.
        ch_bigscape_dir = ch_bigscape_run.map { ran -> [given + ran] }
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
            // groupTuple keeps arrival order, which varies between tools and between runs.
            // Reorder all three lists by sample id so every report has the same columns.
            .map { tool, ids, files, gens ->
                def idx = (0..<ids.size()).toList().sort { ids[it] }
                [tool, idx.collect { ids[it] }, idx.collect { files[it] }, idx.collect { gens[it] }]
            }
            .combine(ch_bigscape_dir)
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

        // One QUAST run over all queries against the reference, unless a folder is supplied.
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

        // One run per tool: ordered query predictions, the reference prediction and genome,
        // and the QUAST folder.
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

        // Empty means neither the reference nor any query produced a usable prediction.
        def banner = "=".multiply(100)

        ch_bgcquast_in = ch_bgcquast_in.ifEmpty {
            error(
                "\n${white}${banner}${creset}\n" +
                "${pink}[bgc_quast_ppl] The reference '${ref_id}' has no predicted BGCs, so QUAST and bgc-quast did not run.${creset}\n\n" +
                "${pink}  Check whether ${active.join(', ')} predicted any BGC in that genome.${creset}\n" +
                "${pink}  The same message appears if no query sample has a predicted BGC either.${creset}\n" +
                "${white}${banner}${creset}"
            )
        }
    }
    else {
        error("[bgc_quast_ppl] bgc_quast_mode='${mode}' is not supported yet. Use compare-tools, compare-samples, or compare-to-reference.")
    }

    BGCQUAST(ch_bgcquast_in)
    ch_versions = ch_versions.mix(BGCQUAST.out.versions)

    emit:
    results  = BGCQUAST.out.results // [ meta, files ]
    tsv      = BGCQUAST.out.tsv     // [ meta, report.tsv ]
    versions = ch_versions          // [ path(versions.yml) ]
}
