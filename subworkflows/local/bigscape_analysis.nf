/*
    Gene cluster family analysis with BiG-SCAPE, compare-samples only.
    One BiG-SCAPE run per active BGC tool over every sample's GBKs together.
*/

include { BIGSCAPE as BIGSCAPE_ANTISMASH } from '../../modules/local/bigscape'
include { BIGSCAPE as BIGSCAPE_GECCO     } from '../../modules/local/bigscape'
include { BIGSCAPE as BIGSCAPE_DEEPBGC   } from '../../modules/local/bigscape'
include { BIGSCAPE_DOWNLOAD_DB           } from '../../modules/local/bigscape_download_db'
include { DEEPBGC_SPLIT_GBK              } from '../../modules/local/deepbgc_split_gbk'

workflow BIGSCAPE_ANALYSIS {
    take:
    antismash_gbk // [ meta, [ gbk ] ] query antiSMASH region GBKs
    gecco_gbk     // [ meta, [ gbk ] ] query GECCO cluster GBKs
    deepbgc_gbk   // [ meta, gbk ]     query DeepBGC multi-record GBK (split first)
    deepbgc_tsv   // [ meta, tsv ]     query DeepBGC BGC table (numbers the split)

    main:
    ch_versions = Channel.empty()

    def orange = params.monochrome_logs ? '' : "\033[38;5;208m"
    def creset = params.monochrome_logs ? '' : "\033[0m"
    def hi     = params.monochrome_logs ? '' : "\033[4m"
    def noh    = params.monochrome_logs ? '' : "\033[24m"

    def bs_tools = []
    if (!params.bgc_skip_antismash) { bs_tools << 'antismash' }
    if (!params.bgc_skip_gecco)     { bs_tools << 'gecco' }
    if (!params.bgc_skip_deepbgc)   { bs_tools << 'deepbgc' }

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
                "\n${white}${banner}${creset}\n" +
                "${pink}[bgc_quast_ppl] ${hi}--bgc_bigscape_dir${noh} " +
                "contains no per-tool subfolder.\n" +
                "                Expected at least one of: " +
                "${bs_tools.join(', ')}\n" +
                "                Looked in: ${params.bgc_bigscape_dir}\n" +
                "                Point it at a previous run's " +
                "bgc_quast/bigscape/ folder.${creset}\n" +
                "${white}${banner}${creset}"
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

    emit:
    dirs     = ch_bigscape_dir // val: [ [ tool: dir ] ] BiG-SCAPE folder per tool
    versions = ch_versions     // [ path(versions.yml) ]
}
