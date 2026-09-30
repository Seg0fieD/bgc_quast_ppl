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
    deepbgc_gbk   // [ meta, gbk ]     query DeepBGC multi-record GBK to split
    deepbgc_tsv   // [ meta, tsv ]     query DeepBGC BGC table for the split

    main:
    ch_versions = Channel.empty()

    def orange = params.monochrome_logs ? '' : "\033[1;38;5;208m"
    def pink   = params.monochrome_logs ? '' : "\033[1;38;5;197m"
    def red    = params.monochrome_logs ? '' : "\033[1;31m"
    def white  = params.monochrome_logs ? '' : "\033[97m"
    def creset = params.monochrome_logs ? '' : "\033[0m"
    def banner = "=".multiply(100)

    def bs_tools = []
    if (!params.bgc_skip_antismash) { bs_tools << 'antismash' }
    if (!params.bgc_skip_gecco)     { bs_tools << 'gecco' }
    if (!params.bgc_skip_deepbgc)   { bs_tools << 'deepbgc' }

    // BiG-SCAPE runs afresh for any tool without a subfolder.
    def given = [:]
    def notes = []

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
                "${pink}[bgc_quast_ppl] --bgc_bigscape_dir " +
                "contains no per-tool subfolder.\n" +
                "                Expected at least one of: " +
                "${bs_tools.join(', ')}\n" +
                "                Looked in: ${params.bgc_bigscape_dir}\n" +
                "                Point it at a previous run's " +
                "bgc_quast/bigscape/ folder.${creset}\n" +
                "${white}${banner}${creset}"
            )
        }

        notes << "BiG-SCAPE folder supplied for: ${given.keySet().join(', ')}"
    }

    def to_run = bs_tools.findAll { !given.containsKey(it) }

    if (to_run) {
        notes << " BiG-SCAPE will run for tool(s): ${to_run.join(', ')}"
    }

    if (notes) {
        def pad = ' '.multiply(16)
        log.info("\n${white}${banner}${creset}\n" +
            "${orange}[bgc_quast_ppl] ${notes.join('\n' + pad)}${creset}\n" +
            "${white}${banner}${creset}")
    }

    // Bare map; wrapped once at the end, since combine() unwraps one level.
    ch_bigscape_run = Channel.value([:])

    if (to_run) {
        // Pfam: the supplied pressed copy, else downloaded and pressed once,
        // then shared by every BiG-SCAPE run.
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

        // Named <sample_id>_<file>, sorted so names and files pair by index;
        // bgc-quast reads the prefix, BiG-SCAPE the '.region' or '_cluster_'.
        def stage_gbks = { ch, tool ->
            ch.flatMap { meta, gbks ->
                    (gbks instanceof List ? gbks : [gbks]).collect { g ->
                        ["${meta.id}_${g.name}".toString(), g, meta.id]
                    }
                }
                .toSortedList { a, b -> a[0] <=> b[0] }
                .map { rows ->
                      def clashes = rows.groupBy { it[0] }
                          .findAll { _n, r -> r.size() > 1 }
                      if (clashes) {
                          def pad   = ' '.multiply(16)
                          def lines = clashes.collect { n, r ->
                              def ids = r.collect { it[2] }.unique().join(', ')
                              "${pad}${n} (samples: ${ids})"
                          }
                          error("\n${white}${banner}${creset}\n" +
                              "${red}[bgc_quast_ppl] BiG-SCAPE input file " +
                              "names clash for ${tool}:\n" +
                              lines.join('\n') + "\n" +
                              "${pad}Rename one of these samples in the " +
                              "samplesheet and run again.${creset}\n" +
                              "${white}${banner}${creset}")
                    }
                    rows
                }                
                .filter { rows -> rows.size() > 0 }
        }

        ch_bigscape_results = Channel.empty()

        if ('antismash' in to_run) {
            def st = stage_gbks(antismash_gbk, 'antiSMASH')
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
            def st = stage_gbks(gecco_gbk, 'GECCO')
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
            // The .bgc.tsv is joined in; the BGC numbering comes from it.
            DEEPBGC_SPLIT_GBK(deepbgc_gbk.join(deepbgc_tsv, failOnDuplicate: true))
            ch_versions = ch_versions.mix(DEEPBGC_SPLIT_GBK.out.versions)

            def st = stage_gbks(DEEPBGC_SPLIT_GBK.out.gbk, 'DeepBGC')
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

        // toList() always emits once; a tool with no GBKs gets no map entry.
        ch_bigscape_run = ch_bigscape_results
            .toList()
            .map { rows -> rows.collectEntries { t, d -> [(t): d] } }
    }

    // Supplied and freshly run folders cannot overlap: to_run excludes the
    // supplied ones. Wrapped in a list, since combine() unwraps one level.
    ch_bigscape_dir = ch_bigscape_run.map { ran -> [given + ran] }

    emit:
    dirs     = ch_bigscape_dir // val: [ [ tool: dir ] ] folder per tool
    versions = ch_versions     // [ path(versions.yml) ]
}
