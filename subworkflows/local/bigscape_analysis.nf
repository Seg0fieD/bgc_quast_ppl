/*
    Gene cluster family analysis with BiG-SCAPE, compare-samples only.
    One BiG-SCAPE run per active BGC tool over every sample's GBKs together.
*/

include { BIGSCAPE as BIGSCAPE_ANTISMASH } from '../../modules/local/bigscape'
include { BIGSCAPE as BIGSCAPE_GECCO     } from '../../modules/local/bigscape'
include { BIGSCAPE as BIGSCAPE_DEEPBGC   } from '../../modules/local/bigscape'
include { BIGSCAPE_DOWNLOAD_DB           } from '../../modules/local/bigscape_download_db'
include { DEEPBGC_SPLIT_GBK              } from '../../modules/local/deepbgc_split_gbk'

// Input GBKs absent from the BiG-SCAPE clustering: sample, region, length.
def bigscapeLeftOut(rows, tsvs) {
    def all_tsv = tsvs instanceof List ? tsvs : [tsvs]
    def mix_tsv = all_tsv.findAll { t -> t.toString().contains('/mix/') }
    def tsv     = (mix_tsv ?: all_tsv).sort { t -> t.toString() }[0]
    def lines   = tsv.readLines()
    def col     = lines ? (lines[0].split('\t') as List).indexOf('GBK') : -1
    if (col < 0) {
        return []
    }

    def kept = lines.drop(1).collect { line ->
        def cells = line.split('\t')
        cells.size() > col ? cells[col] : ''
    } as Set

    return rows
        .findAll { row -> !kept.contains(row[0].replaceFirst(/\.gbk$/, '')) }
        .collect { row ->
            def locus  = row[1].withReader { r -> r.readLine() } ?: ''
            def length = locus =~ /\s(\d+)\s+bp\b/
            [row[2], row[1].baseName, length.find() ? length.group(1) as Long : null]
        }
}

// Run-end note lines for BGCs absent from the BiG-SCAPE clustering, by cause.
def leftOutNote(items) {
    def bigscape_max_len = 500000
    def order  = ['antiSMASH', 'GECCO', 'DeepBGC']
    def sorted = items.sort(false) { a, b ->
        order.indexOf(a[0]) <=> order.indexOf(b[0]) ?:
            a[1] <=> b[1] ?: a[2] <=> b[2]
    }
    def too_long = { r -> r[3] != null && r[3] >= bigscape_max_len }
    def big      = sorted.findAll(too_long)
    def other    = sorted.findAll { r -> !too_long(r) }
    def entry    = { r ->
        def size = r[3] == null ? 'length unknown' :
            String.format(Locale.US, '%,d bp', r[3]) + ', which is ' +
            (too_long(r) ? 'over' : 'under') + ' the 500 kb limit'
        ["  ${r[0]}, ${r[1]}: ${r[2]}", "    ${size}"]
    }

    def one   = sorted.size() == 1
    def lines = [
        "BiG-SCAPE skipped ${sorted.size()} ${one ? 'BGC' : 'BGCs'}, so " +
            "${one ? 'it has' : 'they have'} no gene cluster family.",
        "The report counts ${one ? 'it' : 'them'} in " +
            '"# BGCs not clustered".',
        'BiG-SCAPE accepts only BGCs shorter than 500 kb (500,000 bp).',
    ]
    if (big) {
        lines << 'Skipped for length:'
        big.each { r -> lines.addAll(entry(r)) }
    }
    if (other) {
        lines << 'Skipped for another reason, such as being identical to ' +
            'another file:'
        other.each { r -> lines.addAll(entry(r)) }
    }
    return lines*.toString()
}

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
                "${pink}[bgc_quast_ppl] ERROR:   --bgc_bigscape_dir " +
                "contains no per-tool subfolder.\n" +
                "                         Expected at least one of: " +
                "${bs_tools.join(', ')}\n" +
                "                         Looked in: " +
                "${params.bgc_bigscape_dir}\n" +
                "                         Point it at a previous run's " +
                "bgc_quast/bigscape/ folder.${creset}\n" +
                "${white}${banner}${creset}"
            )
        }

        notes << "BiG-SCAPE folder supplied for: ${given.keySet().join(', ')}"
    }

    def to_run = bs_tools.findAll { !given.containsKey(it) }

    if (to_run) {
        notes << "BiG-SCAPE will run for tool(s): ${to_run.join(', ')}"
    }

    if (notes) {
        def pad = ' '.multiply(25)
        log.info("\n${white}${banner}${creset}\n" +
            "${orange}[bgc_quast_ppl] NOTE:    " +
            "${notes.join('\n' + pad)}${creset}\n" +
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
                        def pad   = ' '.multiply(25)
                        def lines = clashes.collect { n, r ->
                            def ids = r.collect { it[2] }.unique().join(', ')
                            "${pad}${n} (samples: ${ids})"
                        }
                        error("\n${white}${banner}${creset}\n" +
                            "${red}[bgc_quast_ppl] ERROR:   BiG-SCAPE input " +
                            "file names clash for ${tool}:\n" +
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

        def left_out = []
        def note_left_out = { st, ch_tsv, tool ->
            st.map { rows -> [rows] }
                .combine(ch_tsv.map { tsvs -> [tsvs] })
                .subscribe { rows, tsvs ->
                    bigscapeLeftOut(rows, tsvs).each { r ->
                        left_out << [tool] + r
                    }
                }
        }

        def run_info = workflow

        workflow.onComplete {
            if (left_out && run_info.success) {
                def pad = ' '.multiply(25)
                println("${white}${banner}${creset}\n" +
                    "${orange}[bgc_quast_ppl] NOTE:    " +
                    leftOutNote(left_out).join('\n' + pad) + "${creset}\n" +
                    "${white}${banner}${creset}")
            }
        }

        if ('antismash' in to_run) {
            def st = stage_gbks(antismash_gbk, 'antiSMASH')
            BIGSCAPE_ANTISMASH(
                'antismash',
                st.map { rows -> rows.collect { it[0] } },
                st.map { rows -> rows.collect { it[1] } },
                ch_pfam_dir,
                ch_pfam_name,
            )
            note_left_out(st, BIGSCAPE_ANTISMASH.out.clustering, 'antiSMASH')
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
            note_left_out(st, BIGSCAPE_GECCO.out.clustering, 'GECCO')
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
            note_left_out(st, BIGSCAPE_DEEPBGC.out.clustering, 'DeepBGC')
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
