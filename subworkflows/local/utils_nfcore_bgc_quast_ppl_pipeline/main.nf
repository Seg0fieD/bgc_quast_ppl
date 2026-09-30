/*
    Pipeline initialisation and completion for bgc_quast_ppl: pre-run checks,
    samplesheet parsing, failure explainer and completion notices.
*/

include { UTILS_NFSCHEMA_PLUGIN   } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { samplesheetToList       } from 'plugin/nf-schema'
include { completionEmail         } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary       } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification          } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'

workflow PIPELINE_INITIALISATION {
    take:
    version           // boolean: print the version and exit
    validate_params   // boolean: validate parameters against the schema
    monochrome_logs   // boolean: plain log output, no colour codes
    nextflow_cli_args // list: positional Nextflow command-line arguments
    outdir            // string: output directory
    input             // string: samplesheet path

    main:
    ch_versions = Channel.empty()

    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    // Pipeline checks run before nf-schema, so their messages print first.
    validateAntismashMode()
    def sheet = validateSamplesheetContent(input)
    if (params.bgc_quast_mode == 'compare-to-reference') {
        validateReferenceSamplesheet(sheet)
    }
    validatePreRunEnvironment(input)

    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
    )

    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    Channel.fromList(samplesheetToList(sheet, "${projectDir}/assets/schema_input.json"))
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

workflow PIPELINE_COMPLETION {
    take:
    email           // string: address for the completion email
    email_on_fail   // string: address for the email on failure
    plaintext_email // boolean: plain-text email instead of HTML
    outdir          // path: output directory
    monochrome_logs // boolean: plain log output, no colour codes
    hook_url        // string: webhook URL for notifications
    bgcquast_runs   // channel: val(Integer), bgc-quast runs with results

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    // workflow is null inside onComplete; a channel cannot be read there.
    def wf                 = workflow
    def bgcquast_run_total = 0
    bgcquast_runs.subscribe { n -> bgcquast_run_total = n }

    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                [],
            )
        }

        if (wf.errorMessage || bgcquast_run_total > 0) {
            completionSummary(monochrome_logs)
        }
        else {
            reportNoComparison()
        }

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        explainPipelineError()
    }
}

/*
    Terminal escape codes: colours and reset; all empty under
    --monochrome_logs.
*/
def style() {
    def on = !params.monochrome_logs
    return [
        pink  : on ? "\033[1;38;5;197m" : '',
        red   : on ? "\033[1;31m"       : '',
        yellow: on ? "\033[1;93m"       : '',
        white : on ? "\033[97m"         : '',
        reset : on ? "\033[0m"          : '',
    ]
}

/*
    Message block: the first line prefixed with [bgc_quast_ppl], every later
    line indented by the given number of spaces.
*/
def block(lines, indent) {
    def pad = ' '.multiply(indent)
    def out = ["[bgc_quast_ppl] ${lines[0]}"]
    lines.drop(1).each { line -> out << "${pad}${line}" }
    return out.join('\n')
}

/*
    Start-up error text: message lines in pink between two white banners.
*/
def framed(lines) {
    def styl   = style()
    def banner = "=".multiply(100)
    return "\n${styl.white}${banner}${styl.reset}\n" +
        "${styl.pink}${block(lines, 16)}${styl.reset}\n" +
        "${styl.white}${banner}${styl.reset}"
}

/*
    Samplesheet contents: trimmed header cells and trimmed cells of each
    non-blank data row; both empty for a sheet without a data row.
*/
def readSheet(path) {
    def lines = file(path).readLines().findAll { it.trim() }
    def cells = lines.collect { line ->
        line.split(',', -1).collect { it.trim() }
    }
    if (cells.size() < 2) {
        return [header: cells ? cells[0] : [], rows: []]
    }
    return [header: cells[0], rows: cells.drop(1)]
}

/*
    Samplesheet check for every mode: column order, missing paths, duplicate
    names, and names filled in from the file name. Returns the sheet to parse,
    a rewritten temporary copy when any name was filled in.
*/
def validateSamplesheetContent(input) {
    def sheet = readSheet(input)
    if (!sheet.rows) {
        error(framed([
            "The input samplesheet is empty.",
            "Provide at least one sample.",
        ]))
    }

    def header   = sheet.header
    def ref_mode = params.bgc_quast_mode == 'compare-to-reference'

    if (header.size() < 2 || header[0..1] != ['sample', 'fasta']) {
        error(framed([
            "Samplesheet columns are out of order or missing.",
            "Use: sample,fasta,type",
            "The type column belongs to " +
                "compare-to-reference only.",
        ]))
    }

    def si = header.indexOf('sample')
    def fi = header.indexOf('fasta')
    def ti = header.indexOf('type')

    def rewritten = false
    def out_lines = [header.join(',')]
    def seen      = [] as Set

    sheet.rows.each { cells ->
        def name = si < cells.size() ? cells[si] : ''
        def path = fi < cells.size() ? cells[fi] : ''
        def type = (ti >= 0 && ti < cells.size()) ? cells[ti].toLowerCase() : ''

        if (name && !path) {
            error(framed([
                "Sample or reference '${name}' has no path.",
                "Add the file path or directory for that row.",
            ]))
        }

        if (name && seen.contains(name)) {
            error(framed([
                "Duplicate sample name '${name}' in the samplesheet.",
                "Sample names must be unique.",
            ]))
        }

        if (!name && path) {
            def base   = file(path).name.replaceFirst(/\.(fasta|fas|fna|fa)(\.gz)?$/, '')
            def suffix = ref_mode ? (type == 'r' ? '_ref' : '_query') : ''
            def cand   = "${base}${suffix}"
            def n      = 2
            while (seen.contains(cand)) {
                cand = "${base}${suffix}_${n}"
                n++
            }
            name      = cand
            cells[si] = name
            rewritten = true
            log.info("[bgc_quast_ppl] No sample name given for ${path}; " +
                "using '${name}' from the file name.")
        }

        seen << name

        if (path && !(path ==~ /^(https?|ftp):\/\/.*/)) {
            def expanded = path.startsWith('~')
                ? path.replaceFirst('~', System.getProperty('user.home'))
                : path
            if (!file(expanded).exists()) {
                def role = (ref_mode && type == 'r') ? 'reference' : 'sample'
                error(framed([
                    "The path for ${role} '${name}' does not exist:",
                    "${path}",
                    "Check the file path or directory.",
                ]))
            }
        }

        out_lines << cells.join(',')
    }

    if (rewritten) {
        def tmp = File.createTempFile('bgc_quast_ppl_samplesheet_', '.csv')
        tmp.deleteOnExit()
        tmp.text = out_lines.join('\n') + '\n'
        return tmp.absolutePath
    }
    return input
}

/*
    compare-to-reference samplesheet check: a type column, q/r values only,
    and exactly one reference row.
*/
def validateReferenceSamplesheet(input) {
    def mode  = "compare-to-reference"
    def sheet = readSheet(input)

    if (!sheet.header.contains('type')) {
        error(framed([
            "${mode} needs a 'type' column in the",
            "samplesheet. Add it and run again.",
        ]))
    }

    def ti        = sheet.header.indexOf('type')
    def ref_count = 0

    sheet.rows.eachWithIndex { cells, idx ->
        def raw = ti < cells.size() ? cells[ti] : ''
        def t   = raw.toLowerCase()
        if (!(t in ['q', 'r'])) {
            error(framed([
                "${mode}: row ${idx + 2} has type='${raw}', " +
                    "which is not valid.",
                "Use q/Q for a query or r/R for the reference.",
            ]))
        }
        if (t == 'r') {
            ref_count++
        }
    }

    if (ref_count != 1) {
        error(framed([
            "${mode} needs exactly one reference row",
            "(type r/R). Found ${ref_count}.",
        ]))
    }
}

/*
    antiSMASH mode check: --bgc_antismash_minimal and --bgc_antismash_full
    are mutually exclusive.
*/
def validateAntismashMode() {
    def mini = "--bgc_antismash_minimal"
    def full = "--bgc_antismash_full"

    if (params.bgc_antismash_minimal && params.bgc_antismash_full) {
        error(framed([
            "${mini} and ${full} cannot be set together.",
            "Minimal is the default; pass ${full} only for the full analysis.",
        ]))
    }
}

/*
    BiG-SCAPE cutoffs present in a BiG-SCAPE output_files folder, read from
    the _c<cutoff> suffix of each subfolder name.
*/
def cutoffsIn(dir) {
    def found = []
    dir.listFiles().each { d ->
        def m = (d.name =~ /_c([0-9]*\.?[0-9]+)$/)
        if (d.isDirectory() && m) {
            found << (m[0][1] as Double)
        }
    }
    return found.unique().sort()
}

/*
    Problem lines for a --bgc_bigscape_cutoff missing from the given cutoffs;
    null when the cutoff is present.
*/
def cutoffProblem(cuts, whose) {
    def want = params.bgc_bigscape_cutoff as Double
    if (cuts.any { Math.abs(it - want) < 1e-9 }) {
        return null
    }
    return [
        "Cutoff ${params.bgc_bigscape_cutoff} is not available for ${whose}.",
        "Available cutoffs: ${cuts.join(', ')}",
        "Choose one of those with --bgc_bigscape_cutoff.",
    ]
}

/*
    BiG-SCAPE pre-run check, compare-samples with --run_bigscape only: active
    tools, Pfam database, supplied results folder, cutoff, sample names.
    Appends to the given problem and warning lists.
*/
def checkBigscape(input, problems, warnings) {
    if (params.bgc_skip_antismash && params.bgc_skip_deepbgc && params.bgc_skip_gecco) {
        problems << [
            "--run_bigscape is set but every " +
                "BGC tool is skipped.",
            "BiG-SCAPE clusters the BGCs predicted from the tools, so it has",
            "nothing to work on. Enable at least one of tools: " +
                "antiSMASH, DeepBGC",
            "or GECCO.",
        ]
    }

    def pfam_flag = "--bgc_bigscape_pfam"

    if (!params.bgc_bigscape_pfam && !params.bgc_bigscape_dir) {
        warnings << [
            "No ${pfam_flag} given. Pfam will be downloaded and pressed",
            "automatically (about 400 MB, one-off). " +
                "Pass --save_db to keep it, or",
            "${pfam_flag} to use a copy you already have.",
        ]
    }

    if (params.bgc_bigscape_pfam) {
        def hmm = file(params.bgc_bigscape_pfam)
        if (!hmm.exists()) {
            problems << [
                "Pfam file not found: ${params.bgc_bigscape_pfam}",
                "This must be the .hmm file, not the folder holding it.",
            ]
        }
        else {
            def missing = ['h3f', 'h3i', 'h3m', 'h3p'].findAll {
                !file("${hmm}.${it}").exists()
            }
            if (missing) {
                problems << [
                    "Pfam is not pressed. Missing beside ${hmm.name}: " +
                        "${missing.collect { '.' + it }.join(' ')}",
                    "Fix: run  hmmpress ${hmm}",
                ]
            }
        }
    }

    def dir_flag = "--bgc_bigscape_dir"
    def bs_dir   = params.bgc_bigscape_dir ? file(params.bgc_bigscape_dir) : null

    if (bs_dir && !bs_dir.exists()) {
        problems << ["${dir_flag} path not found: ${params.bgc_bigscape_dir}"]
    }

    def bs_tools = []
    if (!params.bgc_skip_antismash) { bs_tools << 'antismash' }
    if (!params.bgc_skip_gecco)     { bs_tools << 'gecco' }
    if (!params.bgc_skip_deepbgc)   { bs_tools << 'deepbgc' }

    def listed = params.bgc_bigscape_cutoffs.toString()
        .split(',').collect { it.trim() as Double }

    // Each supplied tool folder is checked on its own; tools still to run
    // are checked against --bgc_bigscape_cutoffs.
    if (!bs_dir) {
        def p = cutoffProblem(listed, 'this run')
        if (p) { problems << p }
    }
    else if (bs_dir.exists()) {
        def supplied = [:]
        bs_tools.each { t ->
            def d = file("${params.bgc_bigscape_dir}/${t}")
            if (d.exists() && d.isDirectory()) { supplied[t] = d }
        }

        if (!supplied) {
            problems << [
                "${dir_flag} has no per-tool subfolder.",
                "Looked in: ${params.bgc_bigscape_dir}",
                "Expected at least one of: ${bs_tools.join(', ')}",
                "Point it at a previous run's 'bgc_quast/bigscape/' folder.",
            ]
        }
        else {
            def not_bigscape = [
                "Each per-tool subfolder must be a folder BiG-SCAPE",
                "wrote, holding 'output_files'.",
            ]
            supplied.each { tool, dir ->
                def of    = file("${dir}/output_files")
                def found = of.exists() ? cutoffsIn(of) : []
                if (!of.exists()) {
                    problems << [
                        "No BiG-SCAPE results for ${tool}.",
                        "Looked for: ${of}",
                    ] + not_bigscape
                }
                else if (!found) {
                    problems << [
                        "No BiG-SCAPE results for ${tool}.",
                        "${of} holds no cutoff folders.",
                    ] + not_bigscape
                }
                else {
                    def p = cutoffProblem(found, "the supplied ${tool} folder")
                    if (p) { problems << p }
                }
            }

            def to_run = bs_tools.findAll { !supplied.containsKey(it) }
            if (to_run) {
                def whose = "the tools still to run (${to_run.join(', ')})"
                def p     = cutoffProblem(listed, whose)
                if (p) { problems << p }
            }
        }
    }

    // A BGC file-name marker inside a sample name breaks the report join.
    def markers = []
    if (!params.bgc_skip_antismash) { markers << '.region' }
    if (!params.bgc_skip_gecco || !params.bgc_skip_deepbgc) { markers << '_cluster_' }

    def sheet = (input && file(input).exists()) ? readSheet(input) : null
    def si    = sheet ? sheet.header.indexOf('sample') : -1

    if (markers && si >= 0) {
        sheet.rows.eachWithIndex { cells, idx ->
            def name = si < cells.size() ? cells[si] : ''
            def hit  = markers.find { name.contains(it) }
            if (hit) {
                problems << [
                    "Sample name contains '${hit}' (row ${idx + 2}): ${name}",
                    "BiG-SCAPE results are joined on the file name, and",
                    "'${hit}' in a sample id breaks that. Rename the sample.",
                ]
            }
        }
    }
}

/*
    Pre-run environment check: samplesheet, tool databases, BiG-SCAPE
    inputs, Docker and output folder. Problems halt the run together;
    warnings print and the run continues.
*/
def validatePreRunEnvironment(input) {
    def styl     = style()
    def problems = []
    def warnings = []

    def sheet_ok = input && file(input).exists()
    if (!sheet_ok) {
        problems << ["Samplesheet not found: ${input}"]
    }

    if (!params.bgc_skip_antismash) {
        if (!params.bgc_antismash_db) {
            problems << ["antiSMASH is on but " +
                "--bgc_antismash_db is not set."]
        }
        else if (!file(params.bgc_antismash_db).exists()) {
            problems << ["antiSMASH database folder not found: " +
                "${params.bgc_antismash_db}"]
        }
    }

    if (!params.bgc_skip_deepbgc) {
        if (!params.bgc_deepbgc_db) {
            problems << ["DeepBGC is on but " +
                "--bgc_deepbgc_db is not set."]
        }
        else if (!file(params.bgc_deepbgc_db).exists()) {
            problems << ["DeepBGC database folder not found: " +
                "${params.bgc_deepbgc_db}"]
        }
    }

    if (params.bgc_quast_quastdir && !file(params.bgc_quast_quastdir).exists()) {
        problems << ["--bgc_quast_quastdir " +
            "path not found: ${params.bgc_quast_quastdir}"]
    }

    // bgc-quast reads only the mix bin, so no other binning mode is useful.
    if (params.bgc_bigscape_classify != 'none') {
        problems << [
            "--bgc_bigscape_classify must be 'none', " +
                "not '${params.bgc_bigscape_classify}'.",
            "bgc-quast reads the single mixed bin, so another binning",
            "mode would only add bins that the report never reads.",
        ]
    }

    if (params.run_bigscape && params.bgc_quast_mode != 'compare-samples') {
        problems << [
            "BiG-SCAPE does not run in " +
                "${params.bgc_quast_mode}.",
            "'--run_bigscape' works only in " +
                "compare-samples mode.",
            "Remove --run_bigscape, or switch to " +
                "--bgc_quast_mode " +
                "compare-samples.",
        ]
    }

    if (params.run_bigscape && params.bgc_quast_mode == 'compare-samples') {
        checkBigscape(input, problems, warnings)
    }

    if (sheet_ok) {
        def sheet = readSheet(input)
        def fi    = sheet.header.indexOf('fasta')
        sheet.rows.eachWithIndex { cells, idx ->
            def fp = (fi >= 0 && fi < cells.size()) ? cells[fi] : ''
            if (!fp) {
                return
            }
            if (!file(fp).exists()) {
                problems << ["FASTA not found (row ${idx + 2}): ${fp}"]
            }
            else if (!(fp ==~ /(?i).*\.(fa|fasta|fna)(\.gz)?$/)) {
                warnings << ["Row ${idx + 2} file may not be FASTA: ${fp}"]
            }
        }
    }

    if (workflow.containerEngine == 'docker') {
        try {
            def p = ['docker', 'info'].execute()
            p.waitForOrKill(8000)
            if (p.exitValue() != 0) {
                problems << ["Docker does not seem to be running. " +
                    "Start Docker Desktop and retry."]
            }
        }
        catch (Exception ignored) {
            warnings << ["Could not check Docker status. " +
                "Make sure Docker Desktop is running."]
        }
    }

    if (params.outdir) {
        try {
            def od = file(params.outdir)
            if (od.exists() && !od.canWrite()) {
                warnings << ["Output folder may not be writable: " +
                    "${params.outdir}"]
            }
        }
        catch (Exception ignored) {
        }
    }

    // Nextflow prints "WARN: " before each warning, hence 22 spaces.
    warnings.each { lines ->
        log.warn("${styl.yellow}${block(lines, 22)}${styl.reset}")
    }

    if (problems) {
        def fix = [" Run cannot start. Please fix:"]
        problems.each { lines -> fix.addAll(lines.collect { "  ${it}" }) }
        error(framed(fix))
    }
}

/*
    Failure catalogue: per step, the process name to match, a display name,
    hints keyed by known error signature, and a fallback hint.
*/
def failureCatalogue() {
    def antismash_db = [
        'antiSMASH could not load its database. The folder given in',
        '--bgc_antismash_db is incomplete or is not a version 8',
        'database. This pipeline runs antiSMASH v8 and needs a',
        'matching v8 database.',
    ]
    def antismash_short = [
        'No contig in this sample was long enough for antiSMASH to',
        'scan. Use a longer or better assembly, or set',
        '--bgc_mincontiglength lower so shorter contigs pass the',
        'length filter.',
    ]
    def antismash_output = [
        'antiSMASH stopped before writing its results. This is',
        'usually too little memory or disk space. Raise',
        '--max_memory, free some disk, then run again with',
        '-resume. The task folder\'s .command.err has the cause.',
    ]
    def antismash_generic = [
        'antiSMASH failed. Check that --bgc_antismash_db points at an',
        'antiSMASH v8 database and that the input contigs are long',
        'enough to scan.',
    ]
    def split_path = [
        'A contig name in the DeepBGC .bgc.tsv holds a slash or a',
        'space, so it cannot become a file name. Rename the contigs',
        'in the input assembly.',
    ]
    def split_number = [
        'A contig name in the DeepBGC .bgc.tsv reads as a number, so',
        'the split and bgc-quast would build different identifiers.',
        'Rename that contig.',
    ]
    def split_record = [
        'A record in the DeepBGC GenBank file does not hold exactly',
        'one cluster feature, so it is not a normal DeepBGC result.',
        'Re-run that sample.',
    ]
    def split_generic = [
        'Splitting the DeepBGC GenBank file failed. The message above',
        'names the record or TSV row that could not be matched. The',
        'GenBank file and the .bgc.tsv must come from one DeepBGC run.',
    ]
    def deepbgc_db = [
        'DeepBGC could not find its model files. Set --bgc_deepbgc_db',
        'to the folder holding the downloaded DeepBGC database.',
    ]
    def deepbgc_generic = [
        'DeepBGC failed. Check that --bgc_deepbgc_db points at the',
        'downloaded DeepBGC database folder.',
    ]
    def gecco_generic = [
        'GECCO failed. Check that the sample was annotated and has',
        'predicted genes to scan.',
    ]
    def quast_generic = [
        'QUAST failed. Check the query contigs and the reference genome',
        'given in the samplesheet.',
    ]
    def pfam_network = [
        'Could not reach the Pfam FTP server. Check the network, or',
        'download Pfam-A.hmm yourself and pass it with',
        '--bgc_bigscape_pfam.',
    ]
    def pfam_http = [
        'The Pfam download URL returned an error, so the pinned release',
        'may have moved. Check --bgc_bigscape_pfam_url, or download',
        'Pfam-A.hmm yourself and pass it with --bgc_bigscape_pfam.',
    ]
    def pfam_disk = [
        'Not enough disk for Pfam. It needs roughly 4 GB free in the',
        'Nextflow work directory once unpacked and pressed.',
    ]
    def pfam_press = [
        'The Pfam file downloaded but could not be pressed, so it is',
        'probably truncated. Delete the work directory and run again.',
    ]
    def pfam_generic = [
        'Downloading Pfam failed. Download Pfam-A.hmm yourself, run',
        'hmmpress on it, and pass it with --bgc_bigscape_pfam.',
    ]
    def bigscape_domains = [
        'BiG-SCAPE found no protein domains, so every distance came out',
        '1.0 and no families were built. The Pfam database given in',
        '--bgc_bigscape_pfam is wrong or empty. Check it is a real',
        'Pfam-A.hmm with its .h3f/.h3i/.h3m/.h3p files beside it.',
    ]
    def bigscape_press = [
        'BiG-SCAPE tried to press the Pfam database and could not write',
        'to that folder. Run  hmmpress /path/to/Pfam-A.hmm  once by',
        'hand, then run the pipeline again.',
    ]
    def bigscape_output = [
        'BiG-SCAPE produced no output_files/ folder, which means no BGCs',
        'reached it. Check that the tool it ran for found clusters,',
        'a run where every sample has zero BGCs gives it nothing to',
        'cluster.',
    ]
    def bigscape_input = [
        'BiG-SCAPE read zero GBK files. Its --include-gbk filter needs',
        '"region" or "cluster" in each file name. antiSMASH writes',
        '"region", GECCO and DeepBGC write "_cluster_". Check the',
        'staged names in gbk_input/ inside the failed task folder.',
    ]
    def bigscape_generic = [
        'BiG-SCAPE failed. Check --bgc_bigscape_pfam points at a pressed',
        'Pfam-A.hmm file and that the tool it ran for produced BGC GBKs.',
    ]
    def bgcquast_cutoff = [
        'The BiG-SCAPE cutoff you asked for is not in the results. The',
        'message above lists the cutoffs that are there. Pick one of',
        'them with --bgc_bigscape_cutoff.',
    ]
    def bgcquast_bigscape = [
        'The BiG-SCAPE folder is empty or is the wrong folder.',
        '--bgc_bigscape_dir should point at a parent holding per-tool',
        'subfolders (antismash/, gecco/, deepbgc/), each containing',
        '"output_files".',
    ]
    def bgcquast_generic = [
        'bgc-quast failed. Check that the prediction files, the query',
        'FASTA and the QUAST output folder all reached this step.',
    ]

    // First match wins, so a longer process name precedes its prefix.
    return [
        [
            process   : 'ANTISMASH_ANTISMASH',
            name      : 'antiSMASH',
            signatures: [
                'Modules failing prerequisites'   : antismash_db,
                'No matching database in location': antismash_db,
                'too short'                       : antismash_short,
                'Missing output file'             : antismash_output,
            ],
            generic   : antismash_generic,
        ],
        [
            process   : 'DEEPBGC_SPLIT_GBK',
            name      : 'DeepBGC GBK split',
            signatures: [
                'cannot be used in a file name'    : split_path,
                'changes type when bgc-quast reads': split_number,
                'cluster features, expected 1'     : split_record,
            ],
            generic   : split_generic,
        ],
        [
            process   : 'DEEPBGC',
            name      : 'DeepBGC',
            signatures: [
                'DEEPBGC_DOWNLOADS_DIR'                  : deepbgc_db,
                'DeepBGC models directory does not exist': deepbgc_db,
            ],
            generic   : deepbgc_generic,
        ],
        [
            process   : 'GECCO',
            name      : 'GECCO',
            signatures: [:],
            generic   : gecco_generic,
        ],
        [
            process   : 'QUAST',
            name      : 'QUAST',
            signatures: [:],
            generic   : quast_generic,
        ],
        [
            process   : 'BIGSCAPE_DOWNLOAD_DB',
            name      : 'Pfam download',
            signatures: [
                'ConnectionError'         : pfam_network,
                'HTTPError'               : pfam_http,
                'No space left on device' : pfam_disk,
                'hmmpress did not produce': pfam_press,
            ],
            generic   : pfam_generic,
        ],
        [
            process   : 'BIGSCAPE',
            name      : 'BiG-SCAPE',
            signatures: [
                '0 hsps found in this run': bigscape_domains,
                'hmmpress'                : bigscape_press,
                'Missing output file'     : bigscape_output,
                'No files found'          : bigscape_input,
            ],
            generic   : bigscape_generic,
        ],
        [
            process   : 'BGCQUAST',
            name      : 'bgc-quast',
            signatures: [
                'was not found in the output'  : bgcquast_cutoff,
                'No BiG-SCAPE clustering files': bgcquast_bigscape,
            ],
            generic   : bgcquast_generic,
        ],
    ]
}

/*
    Failure message: the failed step and a hint matched on its error report,
    in red between white banners; the raw report only with --bgc_quast_debug.
*/
def explainPipelineError() {
    def styl = style()
    try {
        def report = (workflow.errorReport ?: '') + '\n' + (workflow.errorMessage ?: '')

        // Failed step: last ':' segment of the process, "(sample)" removed.
        def leaf = ''
        def pm   = (report =~ /Process `([^`]+)`/)
        if (pm.find()) {
            leaf = pm.group(1).replaceAll(/\s*\(.*\)$/, '').tokenize(':')[-1]
        }

        // No process name: a start-up check stopped the run and printed.
        if (!leaf) {
            return
        }

        def hit = failureCatalogue().find {
            leaf == it.process || leaf.startsWith(it.process)
        }

        def title = hit
            ? "[bgc_quast_ppl] The ${hit.name} step failed."
            : "[bgc_quast_ppl] The run stopped and the failing step " +
                "could not be identified."

        def detail = hit
            ? (hit.signatures.find { report.contains(it.key) }?.value ?: hit.generic)
            : [
                "Read the error printed above this banner.",
                "It names the failing task and its work folder; open " +
                    ".command.err there for the whole " +
                    "error output.",
            ]

        def pad  = ' '.multiply(16)
        def body = ["${styl.red}${title}${styl.reset}"]
        detail.each { line -> body << "${styl.red}${pad}${line}${styl.reset}" }
        body << "${styl.red}${pad}Troubleshooting: " +
            "https://github.com/Seg0fieD/bgc_quast_ppl" +
            "#14-troubleshooting${styl.reset}"

        // log.error, not println, so it prints after Nextflow's error report.
        def banner = "=".multiply(100)
        log.error("\n${styl.white}${banner}${styl.reset}\n" + body.join('\n') +
            "\n${styl.white}${banner}${styl.reset}")

        if (params.bgc_quast_debug && report.trim()) {
            log.error("${styl.pink}[bgc_quast_ppl] " +
                "--bgc_quast_debug: " +
                "full error report below:\n${report.trim()}${styl.reset}")
        }
    }
    catch (Exception e) {
        log.error("${styl.pink}[bgc_quast_ppl] error handler failed: " +
            "${e}${styl.reset}")
    }
}

/*
    No-comparison notice: a run that ended without error but with no
    bgc-quast result, which would otherwise look like success.
*/
def reportNoComparison() {
    def styl   = style()
    def pad    = ' '.multiply(16)
    def banner = "=".multiply(100)
    def lines  = [
        "[bgc_quast_ppl] Pipeline did NOT complete successfully.",
        "${pad}No BGC comparison was produced. " +
            "bgc-quast never ran, usually",
        "${pad}because every sample was dropped " +
            "before prediction (for example",
        "${pad}all contigs were shorter than " +
            "${params.bgc_mincontiglength} bp, or annotation",
        "${pad}produced no genes). Use longer or " +
            "better assemblies, or lower",
        "${pad}--bgc_mincontiglength, then run again.",
    ]
    println ''
    println "${styl.white}${banner}${styl.reset}"
    lines.each { line -> println "${styl.red}${line}${styl.reset}" }
    println "${styl.white}${banner}${styl.reset}"
}
