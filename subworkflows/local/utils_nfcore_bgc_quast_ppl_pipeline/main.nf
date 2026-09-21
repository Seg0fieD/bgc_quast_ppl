//
// Subworkflow with functionality specific to the bgc_quast_ppl pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN   } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { samplesheetToList       } from 'plugin/nf-schema'
include { completionEmail         } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary       } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification          } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'

// ANSI pink, used for every message that stops the run before any task starts.
def pink(msg) {
    def esc = "\033"
    return "${esc}[95m${msg}${esc}[0m"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Version flag, and the run's parameter dump written to the output folder.
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Pipeline-specific parameter checks.
    //
    validateInputParameters()

    //
    // antiSMASH minimal and full are mutually exclusive.
    //
    validateAntismashMode()

    //
    // Samplesheet content check for every mode. Runs before nf-schema so the
    // per-sample messages appear first. Returns the samplesheet to parse.
    //
    def sheet = validateSamplesheetContent(input)

    //
    // compare-to-reference: the type column and a single reference row.
    //
    if (params.bgc_quast_mode == 'compare-to-reference') {
        validateReferenceSamplesheet(sheet)
    }

    //
    // Pre-run environment checks: paths, databases, Docker.
    //
    validatePreRunEnvironment(input)

    //
    // Validate parameters against the schema and print the parameter summary.
    //
    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
    )

    //
    // Check the config provided to the pipeline.
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Samplesheet channel built from --input.
    //
    Channel.fromList(samplesheetToList(sheet, "${projectDir}/assets/schema_input.json"))
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    bgcquast_runs   // channel: val(Integer) number of bgc-quast runs that produced results

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    // The workflow handle is null inside the onComplete closure, so capture it here.
    def wf = workflow

    //
    // Completion email and summary
    //
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

        // The standard summary on error or when bgc-quast ran; otherwise a notice that
        // the run succeeded without producing anything.
        def comparison_ran = comparisonProduced(outdir)
        if (wf.errorMessage || comparison_ran) {
            completionSummary(monochrome_logs)
        }
        else {
            reportNoComparison(monochrome_logs)
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Pipeline-specific parameter checks; none at present.
//
def validateInputParameters() {
}

//
// Samplesheet content check for all modes: column order, duplicate names, missing
// paths, and auto-naming of blank sample cells. Returns the samplesheet to parse,
// which is a normalised temporary file when any name was filled in.
//
def validateSamplesheetContent(input) {
    def lines = file(input).readLines().findAll { it.trim() }
    if (lines.size() < 2) {
        error(pink("[bgc_quast_ppl] The input samplesheet is empty.\n" +
            "                Provide at least one sample."))
    }

    def header   = lines[0].split(',', -1).collect { it.trim() }
    def ref_mode = params.bgc_quast_mode == 'compare-to-reference'

    // The first two columns must be sample,fasta in that order. Anything else usually
    // means the columns are jumbled. The type column is checked separately.
    def expected = ['sample', 'fasta']
    if (header.size() < expected.size() || header[0..1] != expected) {
        error(pink("[bgc_quast_ppl] Samplesheet columns are out of order or missing.\n" +
            "                Use: sample,fasta,type\n" +
            "                The type column belongs to compare-to-reference only."))
    }

    def si = header.indexOf('sample')
    def fi = header.indexOf('fasta')
    def ti = header.indexOf('type')

    def rewritten = false
    def out_lines = [lines[0]]
    def seen      = [] as Set

    lines[1..-1].eachWithIndex { line, idx ->
        def cells = line.split(',', -1).collect { it.trim() }
        def name  = si < cells.size() ? cells[si] : ''
        def path  = fi < cells.size() ? cells[fi] : ''
        def type  = (ti >= 0 && ti < cells.size()) ? cells[ti].toLowerCase() : ''

        if (name && !path) {
            error(pink("[bgc_quast_ppl] Sample or reference '${name}' has no path.\n" +
                "                Add the file path or directory for that row."))
        }

        // Duplicate names would silently cross-wire two samples downstream.
        if (name && seen.contains(name)) {
            error(pink("[bgc_quast_ppl] Duplicate sample name '${name}' in the samplesheet.\n" +
                "                Sample names must be unique."))
        }

        // A blank sample cell is filled from the file name, with a _query or _ref
        // suffix in compare-to-reference mode and a number if that name is taken.
        if (!name && path) {
            def base   = file(path).name.replaceFirst(/\.(fasta|fas|fna|fa)(\.gz)?$/, '')
            def suffix = ref_mode ? (type == 'r' ? '_ref' : '_query') : ''
            def cand   = "${base}${suffix}"
            def n      = 2
            while (seen.contains(cand)) {
                cand = "${base}${suffix}_${n}"
                n++
            }
            name = cand
            cells[si] = name
            log.info("[bgc_quast_ppl] No sample name given for ${path}; " +
                "using '${name}' from the file name.")
            rewritten = true
        }

        seen << name

        if (path && !(path ==~ /^(https?|ftp):\/\/.*/)) {
            def expanded = path.startsWith('~')
                ? path.replaceFirst('~', System.getProperty('user.home'))
                : path
            if (!file(expanded).exists()) {
                def role = (ref_mode && type == 'r') ? 'reference' : 'sample'
                error(pink("[bgc_quast_ppl] The path for ${role} '${name}' does not exist:\n" +
                    "                ${path}\n" +
                    "                Check the file path or directory."))
            }
        }

        out_lines << cells.join(',')
    }

    // A new file is written only when a name was filled in.
    if (rewritten) {
        def tmp = File.createTempFile('bgc_quast_ppl_samplesheet_', '.csv')
        tmp.deleteOnExit()
        tmp.text = out_lines.join('\n') + '\n'
        return tmp.absolutePath
    }
    return input
}

//
// compare-to-reference samplesheet check: a type column, valid q/r values, and
// exactly one reference row. Content and empty-cell checks run upstream.
//
def validateReferenceSamplesheet(input) {
    def lines  = file(input).readLines().findAll { it.trim() }
    def header = lines[0].split(',', -1).collect { it.trim() }

    if (!header.contains('type')) {
        error(pink("[bgc_quast_ppl] compare-to-reference needs a 'type' column in the\n" +
            "                samplesheet. Add it and run again."))
    }

    def ti        = header.indexOf('type')
    def ref_count = 0

    lines[1..-1].eachWithIndex { line, idx ->
        def cells  = line.split(',', -1)
        def rownum = idx + 2
        def t      = cells[ti].trim().toLowerCase()
        if (!(t in ['q', 'r'])) {
            error(pink("[bgc_quast_ppl] compare-to-reference: row ${rownum} has " +
                "type='${cells[ti].trim()}', which is not valid.\n" +
                "                Use q/Q for a query or r/R for the reference."))
        }
        if (t == 'r') { ref_count++ }
    }

    if (ref_count != 1) {
        error(pink("[bgc_quast_ppl] compare-to-reference needs exactly one reference row\n" +
            "                (type r/R). Found ${ref_count}."))
    }
}

//
// antiSMASH mode check: minimal is the default, and --bgc_antismash_full cannot be
// combined with --bgc_antismash_minimal.
//
def validateAntismashMode() {
    def cli = workflow.commandLine ?: ''
    def minimal_typed = cli.contains('--bgc_antismash_minimal')
    def full_typed    = cli.contains('--bgc_antismash_full')

    if (minimal_typed && full_typed) {
        error(pink("[bgc_quast_ppl] --bgc_antismash_minimal and --bgc_antismash_full cannot\n" +
            "                both be set. Minimal is the default; pass\n" +
            "                --bgc_antismash_full only for the full analysis."))
    }
}

//
// Pre-run environment check: samplesheet, tool databases, BiG-SCAPE inputs, Docker
// and the output folder. Problems are collected and halt the run; warnings print
// and the run continues.
//
def validatePreRunEnvironment(input) {
    def problems = []
    def warnings = []

    // Samplesheet file exists
    def sheet = input ? file(input) : null
    if (!sheet || !sheet.exists()) {
        problems << "Samplesheet not found: ${input}"
    }

    // antiSMASH database, only when antiSMASH runs
    if (!params.bgc_skip_antismash) {
        if (!params.bgc_antismash_db) {
            problems << "antiSMASH is on but --bgc_antismash_db is not set."
        }
        else if (!file(params.bgc_antismash_db).exists()) {
            problems << "antiSMASH database folder not found: ${params.bgc_antismash_db}"
        }
    }

    // DeepBGC database, only when DeepBGC runs
    if (!params.bgc_skip_deepbgc) {
        if (!params.bgc_deepbgc_db) {
            problems << "DeepBGC is on but --bgc_deepbgc_db is not set."
        }
        else if (!file(params.bgc_deepbgc_db).exists()) {
            problems << "DeepBGC database folder not found: ${params.bgc_deepbgc_db}"
        }
    }

    // QUAST folder override, if given, must exist
    if (params.bgc_quast_quastdir && !file(params.bgc_quast_quastdir).exists()) {
        problems << "--bgc_quast_quastdir path not found: ${params.bgc_quast_quastdir}"
    }

    // BiG-SCAPE, only when it is switched on and the mode actually runs it
    if (params.run_bigscape && params.bgc_quast_mode == 'compare-samples') {
        if (params.bgc_skip_antismash && params.bgc_skip_deepbgc && params.bgc_skip_gecco) {
            problems << "--run_bigscape is set but every BGC tool is skipped.\n" +
                "     BiG-SCAPE clusters the BGCs those tools predict, so it has\n" +
                "     nothing to work on. Enable at least one of antiSMASH, DeepBGC\n" +
                "     or GECCO."
        }

        if (!params.bgc_bigscape_pfam && !params.bgc_bigscape_dir) {
            warnings << "No --bgc_bigscape_pfam given. Pfam will be downloaded and pressed\n" +
                "     automatically (about 400 MB, one-off). Pass --save_db to keep it, or\n" +
                "     --bgc_bigscape_pfam to use a copy you already have."
        }

        if (params.bgc_bigscape_pfam) {
            def hmm = file(params.bgc_bigscape_pfam)
            if (!hmm.exists()) {
                problems << "Pfam file not found: ${params.bgc_bigscape_pfam}\n" +
                    "     This must be the .hmm file, not the folder holding it."
            }
            else {
                def missing = ['h3f', 'h3i', 'h3m', 'h3p'].findAll {
                    !file("${hmm}.${it}").exists()
                }
                if (missing) {
                    problems << "Pfam is not pressed. Missing beside ${hmm.name}: " +
                        "${missing.collect { '.' + it }.join(' ')}\n" +
                        "     Fix: run  hmmpress ${hmm}"
                }
            }
        }

        if (params.bgc_bigscape_dir && !file(params.bgc_bigscape_dir).exists()) {
            problems << "--bgc_bigscape_dir path not found: ${params.bgc_bigscape_dir}"
        }

        // The report cutoff must be one that will exist. --bgc_bigscape_dir is a parent
        // of per-tool subfolders, so every supplied folder is checked on its own, and
        // any tool still due to run is checked against the cutoff list instead.
        def want     = params.bgc_bigscape_cutoff as Double
        def bs_tools = []
        if (!params.bgc_skip_antismash) { bs_tools << 'antismash' }
        if (!params.bgc_skip_gecco)     { bs_tools << 'gecco' }
        if (!params.bgc_skip_deepbgc)   { bs_tools << 'deepbgc' }

        def cutoffsIn = { dir ->
            def found = []
            dir.listFiles().each { d ->
                if (d.isDirectory()) {
                    def m = (d.name =~ /_c([0-9]*\.?[0-9]+)$/)
                    if (m) { found << (m[0][1] as Double) }
                }
            }
            found.unique().sort()
        }

        def checkCutoffList = { cuts, whose ->
            if (!cuts.any { Math.abs(it - want) < 1e-9 }) {
                problems << "Cutoff ${params.bgc_bigscape_cutoff} is not available for ${whose}.\n" +
                    "     Available cutoffs: ${cuts.join(', ')}\n" +
                    "     Choose one of those with --bgc_bigscape_cutoff."
            }
        }

        if (params.bgc_bigscape_dir && file(params.bgc_bigscape_dir).exists()) {
            def supplied = [:]
            bs_tools.each { t ->
                def d = file("${params.bgc_bigscape_dir}/${t}")
                if (d.exists() && d.isDirectory()) { supplied[t] = d }
            }

            if (!supplied) {
                problems << "--bgc_bigscape_dir has no per-tool subfolder.\n" +
                    "     Looked in: ${params.bgc_bigscape_dir}\n" +
                    "     Expected at least one of: ${bs_tools.join(', ')}\n" +
                    "     Point it at a previous run's bgc_quast/bigscape/ folder."
            }
            else {
                supplied.each { tool, dir ->
                    def of = file("${dir}/output_files")
                    if (!of.exists()) {
                        problems << "No BiG-SCAPE results for ${tool}.\n" +
                            "     Looked for: ${of}\n" +
                            "     Each per-tool subfolder must be a folder BiG-SCAPE\n" +
                            "     wrote, holding 'output_files'."
                    }
                    else {
                        def found = cutoffsIn(of)
                        if (!found) {
                            problems << "No BiG-SCAPE results for ${tool}.\n" +
                                "     ${of} holds no cutoff folders.\n" +
                                "     Each per-tool subfolder must be a folder BiG-SCAPE\n" +
                                "     wrote, holding 'output_files'."
                        }
                        else {
                            checkCutoffList(found, "the supplied ${tool} folder")
                        }
                    }
                }

                // Tools without a supplied folder still run, so they are checked against
                // the cutoff list rather than a folder.
                def willRun = bs_tools.findAll { !supplied.containsKey(it) }
                if (willRun) {
                    def cuts = params.bgc_bigscape_cutoffs.toString()
                        .split(',').collect { it.trim() as Double }
                    checkCutoffList(cuts, "the tools still to run (${willRun.join(', ')})")
                }
            }
        }
        else if (!params.bgc_bigscape_dir) {
            def cuts = params.bgc_bigscape_cutoffs.toString()
                .split(',').collect { it.trim() as Double }
            checkCutoffList(cuts, 'this run')
        }

        // A BGC file-name marker inside a sample id breaks the id the report joins on.
        // antiSMASH writes ".region", GECCO and DeepBGC write "_cluster_".
        def markers = []
        if (!params.bgc_skip_antismash) { markers << '.region' }
        if (!params.bgc_skip_gecco || !params.bgc_skip_deepbgc) { markers << '_cluster_' }

        if (markers && sheet && sheet.exists()) {
            def blines = sheet.readLines().findAll { it.trim() }
            if (blines.size() >= 2) {
                def bheader = blines[0].split(',', -1).collect { it.trim() }
                def bsi     = bheader.indexOf('sample')
                if (bsi >= 0) {
                    blines[1..-1].eachWithIndex { line, idx ->
                        def cells = line.split(',', -1)
                        if (bsi < cells.size()) {
                            def name = cells[bsi].trim()
                            def hit  = markers.find { name.contains(it) }
                            if (hit) {
                                problems << "Sample name contains '${hit}' (row ${idx + 2}): ${name}\n" +
                                    "     BiG-SCAPE results are joined on the file name, and\n" +
                                    "     '${hit}' in a sample id breaks that. Rename the sample."
                            }
                        }
                    }
                }
            }
        }

    }

    // FASTA files listed in the samplesheet
    if (sheet && sheet.exists()) {
        def lines = sheet.readLines().findAll { it.trim() }
        if (lines.size() >= 2) {
            def header = lines[0].split(',', -1).collect { it.trim() }
            def fi = header.indexOf('fasta')
            if (fi >= 0) {
                lines[1..-1].eachWithIndex { line, idx ->
                    def cells  = line.split(',', -1)
                    def rownum = idx + 2
                    if (fi < cells.size()) {
                        def fp = cells[fi].trim()
                        if (fp) {
                            if (!file(fp).exists()) {
                                problems << "FASTA not found (row ${rownum}): ${fp}"
                            }
                            else if (!(fp ==~ /(?i).*\.(fa|fasta|fna)(\.gz)?$/)) {
                                warnings << "Row ${rownum} file may not be FASTA: ${fp}"
                            }
                        }
                    }
                }
            }
        }
    }

    // Docker, only when the docker engine is active
    if (workflow.containerEngine == 'docker') {
        try {
            def p = ['docker', 'info'].execute()
            p.waitForOrKill(8000)
            if (p.exitValue() != 0) {
                problems << "Docker does not seem to be running. Start Docker Desktop and retry."
            }
        }
        catch (Exception e) {
            warnings << "Could not check Docker status. Make sure Docker Desktop is running."
        }
    }

    // Output folder writable, warning only
    if (params.outdir) {
        try {
            def od = file(params.outdir)
            if (od.exists() && !od.canWrite()) {
                warnings << "Output folder may not be writable: ${params.outdir}"
            }
        }
        catch (Exception e) {
            // ignore
        }
    }

    warnings.each { log.warn("[bgc_quast_ppl] ${it}") }

    // Every blocking problem is printed together, then the run halts.
    if (problems) {
        def msg = problems.collect { " - ${it}" }.join('\n')
        error(pink("[bgc_quast_ppl] Cannot start. Please fix:\n${msg}"))
    }
}

//
// Failure explainer: the step that failed and what to do about it, matched on the
// process name and on known error signatures. The raw report is printed only with
// --bgc_quast_debug.
//
def explainPipelineError() {
    try {
        def report = (workflow.errorReport ?: '') + '\n' + (workflow.errorMessage ?: '')

        // Failed step name: last ':' segment, trailing "(sample)" removed.
        def leaf = ''
        def pm = (report =~ /Process `([^`]+)`/)
        if (pm.find()) {
            def full = pm.group(1).replaceAll(/\s*\(.*\)$/, '')
            leaf = full.tokenize(':')[-1]
        }

        // No process name means the run stopped on one of the checks above, which
        // already printed its own message.
        if (!leaf) {
            return
        }

        // Per step: process to match, display name, known error signatures, and a
        // fallback. Matching is exact or by prefix, and the first hit wins, so a more
        // specific process name must be listed before a shorter one it starts with.
        def tools = [
            [
                process   : 'ANTISMASH_ANTISMASH',
                name      : 'antiSMASH',
                signatures: [
                    [ match: 'Modules failing prerequisites',
                      hint : 'antiSMASH could not load its database. The folder given in\n' +
                             '  --bgc_antismash_db is incomplete or is not a version 8\n' +
                             '  database. This pipeline runs antiSMASH v8 and needs a\n' +
                             '  matching v8 database.' ],
                    [ match: 'No matching database in location',
                      hint : 'antiSMASH could not load its database. The folder given in\n' +
                             '  --bgc_antismash_db is incomplete or is not a version 8\n' +
                             '  database. This pipeline runs antiSMASH v8 and needs a\n' +
                             '  matching v8 database.' ],
                    [ match: 'too short',
                      hint : 'No contig in this sample was long enough for antiSMASH to\n' +
                             '  scan. Use a longer or better assembly, or set\n' +
                             '  --bgc_mincontiglength lower so shorter contigs pass the\n' +
                             '  length filter.' ],
                    [ match: 'Missing output file',
                      hint : 'antiSMASH finished but found no BGCs in this sample, so it\n' +
                             '  wrote no HTML result files while the module still requires\n' +
                             '  them. Mark the antiSMASH HTML outputs as optional so a\n' +
                             '  no-cluster result is allowed.' ],
                ],
                generic   : 'antiSMASH failed. Check that --bgc_antismash_db points at an\n' +
                            '  antiSMASH v8 database and that the input contigs are long\n' +
                            '  enough to scan.',
            ],
            [
                process   : 'DEEPBGC_SPLIT_GBK',
                name      : 'DeepBGC GBK split',
                signatures: [
                    [ match: 'cannot be used in a file name',
                      hint : 'A contig name in the DeepBGC .bgc.tsv holds a slash or a\n' +
                             '  space, so it cannot become a file name. Rename the contigs\n' +
                             '  in the input assembly.' ],
                    [ match: 'changes type when bgc-quast reads',
                      hint : 'A contig name in the DeepBGC .bgc.tsv reads as a number, so\n' +
                             '  the split and bgc-quast would build different identifiers.\n' +
                             '  Rename that contig.' ],
                    [ match: 'cluster features, expected 1',
                      hint : 'A record in the DeepBGC GenBank file does not hold exactly\n' +
                             '  one cluster feature, so it is not a normal DeepBGC result.\n' +
                             '  Re-run that sample.' ],
                ],
                generic   : 'Splitting the DeepBGC GenBank file failed. The message above\n' +
                            '  names the record or TSV row that could not be matched. The\n' +
                            '  GenBank file and the .bgc.tsv must come from one DeepBGC run.',
            ],
            [
                process   : 'DEEPBGC',
                name      : 'DeepBGC',
                signatures: [
                    [ match: 'DEEPBGC_DOWNLOADS_DIR',
                      hint : 'DeepBGC could not find its model files. Set --bgc_deepbgc_db\n' +
                             '  to the folder holding the downloaded DeepBGC database.' ],
                    [ match: 'DeepBGC models directory does not exist',
                      hint : 'DeepBGC could not find its model files. Set --bgc_deepbgc_db\n' +
                             '  to the folder holding the downloaded DeepBGC database.' ],
                ],
                generic   : 'DeepBGC failed. Check that --bgc_deepbgc_db points at the\n' +
                            '  downloaded DeepBGC database folder.',
            ],
            [
                process   : 'GECCO',
                name      : 'GECCO',
                signatures: [],
                generic   : 'GECCO failed. Check that the sample was annotated and has\n' +
                            '  predicted genes to scan.',
            ],
            [
                process   : 'QUAST',
                name      : 'QUAST',
                signatures: [],
                generic   : 'QUAST failed. Check the query contigs and the reference genome\n' +
                            '  given in the samplesheet.',
            ],
            [
                process   : 'BIGSCAPE_DOWNLOAD_DB',
                name      : 'Pfam download',
                signatures: [
                    [ match: 'ConnectionError',
                      hint : 'Could not reach the Pfam FTP server. Check the network, or\n' +
                             '  download Pfam-A.hmm yourself and pass it with\n' +
                             '  --bgc_bigscape_pfam.' ],
                    [ match: 'HTTPError',
                      hint : 'The Pfam download URL returned an error, so the pinned release\n' +
                             '  may have moved. Check --bgc_bigscape_pfam_url, or download\n' +
                             '  Pfam-A.hmm yourself and pass it with --bgc_bigscape_pfam.' ],
                    [ match: 'No space left on device',
                      hint : 'Not enough disk for Pfam. It needs roughly 4 GB free in the\n' +
                             '  Nextflow work directory once unpacked and pressed.' ],
                    [ match: 'hmmpress did not produce',
                      hint : 'The Pfam file downloaded but could not be pressed, so it is\n' +
                             '  probably truncated. Delete the work directory and run again.' ],
                ],
                generic   : 'Downloading Pfam failed. Download Pfam-A.hmm yourself, run\n' +
                            '  hmmpress on it, and pass it with --bgc_bigscape_pfam.',
            ],
            [
                process   : 'BIGSCAPE',
                name      : 'BiG-SCAPE',
                signatures: [
                    [ match: '0 hsps found in this run',
                      hint : 'BiG-SCAPE found no protein domains, so every distance came out\n' +
                             '  1.0 and no families were built. The Pfam database given in\n' +
                             '  --bgc_bigscape_pfam is wrong or empty. Check it is a real\n' +
                             '  Pfam-A.hmm with its .h3f/.h3i/.h3m/.h3p files beside it.' ],
                    [ match: 'hmmpress',
                      hint : 'BiG-SCAPE tried to press the Pfam database and could not write\n' +
                             '  to that folder. Run  hmmpress /path/to/Pfam-A.hmm  once by\n' +
                             '  hand, then run the pipeline again.' ],
                    [ match: 'Missing output file',
                      hint : 'BiG-SCAPE produced no output_files/ folder, which means no BGCs\n' +
                             '  reached it. Check that the tool it ran for found clusters, \n' +
                             '  a run where every sample has zero BGCs gives it nothing to\n' +
                             '  cluster.' ],
                    [ match: 'No files found',
                      hint : 'BiG-SCAPE read zero GBK files. Its --include-gbk filter needs\n' +
                             '  "region" or "cluster" in each file name. antiSMASH writes\n' +
                             '  "region", GECCO and DeepBGC write "_cluster_". Check the\n' +
                             '  staged names in gbk_input/ inside the failed task folder.' ],
                ],
                generic   : 'BiG-SCAPE failed. Check --bgc_bigscape_pfam points at a pressed\n' +
                            '  Pfam-A.hmm file and that the tool it ran for produced BGC GBKs.',
            ],
            [
                process   : 'BGCQUAST',
                name      : 'bgc-quast',
                signatures: [
                    [ match: 'was not found in the output',
                      hint : 'The BiG-SCAPE cutoff you asked for is not in the results. The\n' +
                             '  message above lists the cutoffs that are there. Pick one of\n' +
                             '  them with --bgc_bigscape_cutoff.' ],
                    [ match: 'No BiG-SCAPE clustering files',
                      hint : 'The BiG-SCAPE folder is empty or is the wrong folder.\n' +
                             '  --bgc_bigscape_dir should point at a parent holding per-tool\n' +
                             '  subfolders (antismash/, gecco/, deepbgc/), each containing\n' +
                             '  "output_files".' ],
                ],
                generic   : 'bgc-quast failed. Check that the prediction files, the query\n' +
                            '  FASTA and the QUAST output folder all reached this step.',
            ],
        ]

        def hit = tools.find { leaf == it.process || leaf.startsWith(it.process) }

        def banner = "=".multiply(100)
        def red    = params.monochrome_logs ? '' : "\033[1;31m"
        def white  = params.monochrome_logs ? '' : "\033[97m"
        def reset  = params.monochrome_logs ? '' : "\033[0m"
        def hi     = params.monochrome_logs ? '' : "\033[4m"
        def noh    = params.monochrome_logs ? '' : "\033[24m"

        def title = hit
            ? "[bgc_quast_ppl] The ${hit.name} step failed."
            : "[bgc_quast_ppl] The run stopped and the failing step could not be identified."

        def detail = hit
            ? (hit.signatures.find { report.contains(it.match) }?.hint ?: hit.generic)
            : "Read the error printed above this banner.\n" +
              "It names the failing task and its work folder; open " +
              "${hi}.command.err${noh} there for the whole error output."

        // println, not log.error, so the banners stay uncoloured.
        println ''
        println "${white}${banner}${reset}"
        println "${red}${title}${reset}"
        println ''
        detail.readLines().each { println "${red}  ${it}${reset}" }
        println "${white}${banner}${reset}"

        if (params.bgc_quast_debug && report.trim()) {
            log.error(pink("[bgc_quast_ppl] --bgc_quast_debug: full error report below:\n" +
                "${report.trim()}"))
        }

        log.error("Please refer to troubleshooting docs: " +
            "https://nf-co.re/docs/usage/troubleshooting")
    }
    catch (Exception e) {
        log.error(pink("[bgc_quast_ppl] error handler failed: ${e}"))
    }
}

//
// True when the mode's output folder exists and is not empty. The check reads the
// output folder rather than a channel, because a channel cannot be read on completion.
//
def comparisonProduced(outdir) {
    try {
        def mode_dir = params.bgc_quast_mode.replaceAll('-', '_')
        def out_dir  = file("${outdir}/bgc_quast/${mode_dir}")
        return out_dir.exists() && out_dir.list() && out_dir.list().size() > 0
    }
    catch (Exception e) {
        log.warn("[bgc_quast_ppl] completion check failed: ${e}")
        return true
    }
}

//
// Notice printed when the run ended without error but bgc-quast produced no
// comparison, which otherwise looks like success.
//
def reportNoComparison(monochrome_logs) {
    def red    = monochrome_logs ? '' : "\033[1;31m"
    def white  = monochrome_logs ? '' : "\033[97m"
    def reset  = monochrome_logs ? '' : "\033[0m"
    def banner = "=".multiply(100)
    println ''
    println "${white}${banner}${reset}"
    println "${red}[bgc_quast_ppl] Pipeline did NOT complete successfully.${reset}"
    println ''
    println "${red}  No BGC comparison was produced. bgc-quast never ran, usually because every${reset}"
    println "${red}  sample was dropped before prediction (for example all contigs were shorter than${reset}"
    println "${red}  ${params.bgc_mincontiglength} bp, or annotation produced no genes).${reset}"
    println ''
    println "${red}  Use longer or better assemblies, or lower --bgc_mincontiglength, then run again.${reset}"
    println "${white}${banner}${reset}"
}

//
// Samplesheet channel check: every run of one sample shares a datatype.
//
def validateInputSamplesheet(input) {
    def (metas, fastas) = input[1..2]

    def endedness_ok = metas.collect { meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of " +
            "the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [metas[0], fastas]
}
