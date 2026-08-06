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

// ANSI pink for validation error messages.
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
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // antiSMASH: minimal and full are mutually exclusive.
    //
    validateAntismashMode()

    //
    // All-modes samplesheet content check. Runs before nf-schema so our per-sample
    // messages fire first. Returns the sheet to parse (normalised if any sample name
    // was auto-filled from its filename).
    //
    def sheet = validateSamplesheetContent(input)

    //
    // Compare-to-reference: enforce the type column and a single reference row.
    //
    if (params.bgc_quast_mode == 'compare-to-reference') {
        validateReferenceSamplesheet(sheet)
    }

    //
    // Pre-run environment checks: paths, databases, Docker.
    //
    validatePreRunEnvironment(input)

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
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

    // Capture the workflow handle; it is null inside the onComplete closure.
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

        // Show the standard summary on error or when bgc-quast ran; otherwise flag the false success.
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
// Check and validate pipeline parameters
//
def validateInputParameters() {
}

//
// All-modes samplesheet content check. Runs before schema parsing so it can report
// clear, per-sample messages. Auto-names rows with a blank sample column and returns
// the samplesheet path to use for parsing (a normalised temp file if any name was added).
//
def validateSamplesheetContent(input) {
    def lines = file(input).readLines().findAll { it.trim() }
    if (lines.size() < 2) {
        error(pink("[bgc_quast_ppl] The input samplesheet is empty. Please provide at least one sample."))
    }

    def header    = lines[0].split(',', -1).collect { it.trim() }
    def ref_mode  = params.bgc_quast_mode == 'compare-to-reference'

    // The first two columns must be sample,fasta in order. A mismatch usually means the
    // columns are jumbled. The type column (compare-to-reference) is checked separately.
    def expected = ['sample', 'fasta']
    if (header.size() < expected.size() || header[0..1] != expected) {
        error(pink("[bgc_quast_ppl] Samplesheet columns are out of order or missing. Use: sample,fasta,type (type only in compare-to-reference)."))
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

        // Sample name given but no path.
        if (name && !path) {
            error(pink("[bgc_quast_ppl] Sample or reference '${name}' has no path mentioned. Please add the correct file path or directory location for that sample or reference."))
        }

        // User-typed name: reject duplicates so samples are not silently cross-wired.
        if (name && seen.contains(name)) {
            error(pink("[bgc_quast_ppl] Duplicate sample name '${name}' in the samplesheet. Sample names must be unique."))
        }

        // No sample name: auto-name from the filename, with a _query/_ref suffix by type
        // in compare-to-reference mode (no suffix in the other modes). A numeric suffix is
        // added if the generated name clashes with one already seen.
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
            log.info("[bgc_quast_ppl] No sample name given for ${path}; using '${name}' from the filename.")
            rewritten = true
        }

        seen << name

        // Path present: check the file exists, with a per-sample/reference message.
        if (path && !(path ==~ /^(https?|ftp):\/\/.*/)) {
            def expanded = path.startsWith('~') ? path.replaceFirst('~', System.getProperty('user.home')) : path
            if (!file(expanded).exists()) {
                def role = (ref_mode && type == 'r') ? 'reference' : 'sample'
                error(pink("[bgc_quast_ppl] The path for ${role} '${name}' does not exist: ${path}. Please check the file path or directory location."))
            }
        }

        out_lines << cells.join(',')
    }

    // Only write a new file if a name was auto-filled; otherwise use the original.
    if (rewritten) {
        def tmp = File.createTempFile('bgc_quast_ppl_samplesheet_', '.csv')
        tmp.deleteOnExit()
        tmp.text = out_lines.join('\n') + '\n'
        return tmp.absolutePath
    }
    return input
}

//
// compare-to-reference type-column check: type present, valid values, exactly one reference.
// Content and empty-cell checks are handled upstream in validateSamplesheetContent.
//
def validateReferenceSamplesheet(input) {
    def lines  = file(input).readLines().findAll { it.trim() }
    def header = lines[0].split(',', -1).collect { it.trim() }

    if (!header.contains('type')) {
        error(pink("[bgc_quast_ppl] compare-to-reference needs a 'type' column in the samplesheet. Please add it and restart the run."))
    }

    def ti        = header.indexOf('type')
    def ref_count = 0

    lines[1..-1].eachWithIndex { line, idx ->
        def cells  = line.split(',', -1)
        def rownum = idx + 2
        def t      = cells[ti].trim().toLowerCase()
        if (!(t in ['q', 'r'])) {
            error(pink("[bgc_quast_ppl] compare-to-reference: row ${rownum} type='${cells[ti].trim()}' is invalid. Use q/Q (query) or r/R (reference)."))
        }
        if (t == 'r') { ref_count++ }
    }

    if (ref_count != 1) {
        error(pink("[bgc_quast_ppl] compare-to-reference needs exactly one reference row (type r/R). Found ${ref_count}."))
    }
}

//
// Minimal mode is the default; --bgc_antismash_full switches to full analysis mode. Passing both throws errors.
//
def validateAntismashMode() {
    def cli = workflow.commandLine ?: ''
    def minimal_typed = cli.contains('--bgc_antismash_minimal')
    def full_typed    = cli.contains('--bgc_antismash_full')

    if (minimal_typed && full_typed) {
        error(pink("[bgc_quast_ppl] --bgc_antismash_minimal and --bgc_antismash_full cannot both be set. Minimal is the default; pass --bgc_antismash_full only if you want the full analysis."))
    }
}

//
// Pre-run environment checks. problems -> collected, then halt. warnings -> printed, continue.
//
def validatePreRunEnvironment(input) {
    def problems = []
    def warnings = []

    // Samplesheet file exists
    def sheet = input ? file(input) : null
    if (!sheet || !sheet.exists()) {
        problems << "Samplesheet not found: ${input}"
    }

    // antiSMASH database (only if antiSMASH is on)
    if (!params.bgc_skip_antismash) {
        if (!params.bgc_antismash_db) {
            problems << "antiSMASH is on but --bgc_antismash_db is not set."
        }
        else if (!file(params.bgc_antismash_db).exists()) {
            problems << "antiSMASH database folder not found: ${params.bgc_antismash_db}"
        }
    }

    // DeepBGC database (only if DeepBGC is on)
    if (!params.bgc_skip_deepbgc) {
        if (!params.bgc_deepbgc_db) {
            problems << "DeepBGC is on but --bgc_deepbgc_db is not set."
        }
        else if (!file(params.bgc_deepbgc_db).exists()) {
            problems << "DeepBGC database folder not found: ${params.bgc_deepbgc_db}"
        }
    }

    // QUAST directory override, if given, must exist
    if (params.bgc_quast_quastdir && !file(params.bgc_quast_quastdir).exists()) {
        problems << "--bgc_quast_quastdir path not found: ${params.bgc_quast_quastdir}"
    }

    // BiG-SCAPE (only if it is switched on)
    if (params.run_bigscape) {
        if (params.bgc_skip_antismash) {
            problems << "--run_bigscape needs antiSMASH, but --bgc_skip_antismash is set. Nothing would feed BiG-SCAPE."
        }

       if (!params.bgc_bigscape_pfam && !params.bgc_bigscape_dir) {
            warnings << "No --bgc_bigscape_pfam given. Pfam will be downloaded and pressed automatically (about 400 MB, one-off).\n     Pass --save_db to keep it, or --bgc_bigscape_pfam to use a copy you already have."
        }

        if (params.bgc_bigscape_pfam) {
            def hmm = file(params.bgc_bigscape_pfam)
            if (!hmm.exists()) {
                problems << "Pfam file not found: ${params.bgc_bigscape_pfam}  (this must be the .hmm file, not its folder)"
            }
            else {
                def missing = ['h3f', 'h3i', 'h3m', 'h3p'].findAll { !file("${hmm}.${it}").exists() }
                if (missing) {
                    problems << "Pfam is not pressed. Missing beside ${hmm.name}: ${missing.collect { '.' + it }.join(' ')}\n     Fix: run  hmmpress ${hmm}"
                }
            }
        }

        if (params.bgc_bigscape_dir && !file(params.bgc_bigscape_dir).exists()) {
            problems << "--bgc_bigscape_dir path not found: ${params.bgc_bigscape_dir}"
        }

        // The report cutoff must be one BiG-SCAPE actually computes
        def cuts = params.bgc_bigscape_cutoffs.toString().split(',').collect { it.trim() as Double }
        if (!cuts.any { Math.abs(it - (params.bgc_bigscape_cutoff as Double)) < 1e-9 }) {
            problems << "--bgc_bigscape_cutoff ${params.bgc_bigscape_cutoff} is not in --bgc_bigscape_cutoffs '${params.bgc_bigscape_cutoffs}'."
        }

        // ".region" in a sample id breaks the BGC id the report joins on
        if (sheet && sheet.exists()) {
            def blines = sheet.readLines().findAll { it.trim() }
            if (blines.size() >= 2) {
                def bheader = blines[0].split(',', -1).collect { it.trim() }
                def bsi     = bheader.indexOf('sample')
                if (bsi >= 0) {
                    blines[1..-1].eachWithIndex { line, idx ->
                        def cells = line.split(',', -1)
                        if (bsi < cells.size() && cells[bsi].trim().contains('.region')) {
                            problems << "Sample name contains '.region' (row ${idx + 2}): ${cells[bsi].trim()}\n     BiG-SCAPE results are joined on the file name, and '.region' in a sample id breaks that. Rename the sample."
                        }
                    }
                }
            }
        }

        if (params.bgc_quast_mode != 'compare-samples') {
            warnings << "BiG-SCAPE is only wired into compare-samples so far. In '${params.bgc_quast_mode}' it will run and publish its folder, but no GCF rows will appear in the report."
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

    // Docker running (only when the docker engine is active)
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

    // Output folder writable — warn only
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

    // Print warnings; run continues
    warnings.each { log.warn("[bgc_quast_ppl] ${it}") }

    // Print all blocking problems together, then halt
    if (problems) {
        def msg = problems.collect { " - ${it}" }.join('\n')
        error(pink("[bgc_quast_ppl] Cannot start. Please fix:\n${msg}"))
    }
}

//
// On failure, report which step failed and how to fix it. Raw error only with --bgc_quast_debug.
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

        // Per step: process to match, display name, known error signatures, generic fallback.
        def tools = [
            [
                process   : 'ANTISMASH_ANTISMASH',
                name      : 'antiSMASH',
                signatures: [
                    [ match: 'Modules failing prerequisites',
                    hint : 'antiSMASH could not load its database. The path/directory in --bgc_antismash_db is missing files or is not a version 8 database. \n  This pipeline runs antiSMASH v8, which needs a matching antiSMASH v8 database. Set --bgc_antismash_db to a v8 database folder.' ],
                    [ match: 'No matching database in location',
                    hint : 'antiSMASH could not load its database. The path/directory in --bgc_antismash_db is missing files or is not a version 8 database. \n  This pipeline runs antiSMASH v8, which needs a matching antiSMASH v8 database. Set --bgc_antismash_db to a v8 database folder.' ],
                    [ match: 'too short',
                    hint : 'No contig in this sample was long enough for antiSMASH to scan. Use a longer or better assembly, or \n  set --bgc_mincontiglength lower so shorter contigs pass the length filter.' ],
                    [ match: 'Missing output file',
                    hint : 'antiSMASH finished but found no BGCs in this sample. With no clusters to show, it did not write its HTML result files, \n  but the pipeline still requires them. Fix: mark the antiSMASH HTML outputs as optional in the antiSMASH module so a no-cluster result is allowed.' ],
                ],
                generic   : 'antiSMASH failed. Check that --bgc_antismash_db points to an appropriate antiSMASH v8 database and that the input contigs are long enough to scan.',
            ],
            [
                process   : 'DEEPBGC',
                name      : 'DeepBGC',
                signatures: [
                    [ match: 'DEEPBGC_DOWNLOADS_DIR',
                    hint : 'DeepBGC could not find its model files. Set --bgc_deepbgc_db to the folder or path that holds the downloaded DeepBGC database.' ],
                    [ match: 'DeepBGC models directory does not exist',
                    hint : 'DeepBGC could not find its model files. Set --bgc_deepbgc_db to the folder or path that holds the downloaded DeepBGC database.' ],
                ],
                generic   : 'DeepBGC failed. Check that --bgc_deepbgc_db points to the downloaded DeepBGC database folder.',
            ],
            [
                process   : 'GECCO',
                name      : 'GECCO',
                signatures: [],
                generic   : 'GECCO failed. Check that the sample was annotated and has predicted genes to scan.',
            ],
            [
                process   : 'QUAST',
                name      : 'QUAST',
                signatures: [],
                generic   : 'QUAST failed. Check the query contigs and the reference genome given in the samplesheet.',
            ],
            [
                process   : 'BIGSCAPE_DOWNLOAD_DB',
                name      : 'Pfam download',
                signatures: [
                    [ match: 'ConnectionError',
                    hint : 'Could not reach the Pfam FTP server. Check the network, or download Pfam-A.hmm yourself and pass it with --bgc_bigscape_pfam.' ],
                    [ match: 'HTTPError',
                    hint : 'The Pfam download URL returned an error. The pinned release may have moved. \n  Check --bgc_bigscape_pfam_url, or download Pfam-A.hmm yourself and pass it with --bgc_bigscape_pfam.' ],
                    [ match: 'No space left on device',
                    hint : 'Not enough disk for Pfam. It needs roughly 4 GB free in the Nextflow work directory once unpacked and pressed.' ],
                    [ match: 'hmmpress did not produce',
                    hint : 'The Pfam file downloaded but could not be pressed, so it is probably truncated or corrupt. Delete the work directory and run again.' ],
                ],
                generic   : 'Downloading Pfam failed. Download Pfam-A.hmm yourself, run hmmpress on it, and pass it with --bgc_bigscape_pfam.',
            ],
            [
                process   : 'BIGSCAPE',
                name      : 'BiG-SCAPE',
                signatures: [
                    [ match: '0 hsps found in this run',
                    hint : 'BiG-SCAPE found no protein domains, so every distance came out 1.0 and no families were built. \n  The Pfam database in --bgc_bigscape_pfam is wrong or empty. Check it is a real Pfam-A.hmm and that its .h3f/.h3i/.h3m/.h3p files are beside it.' ],
                    [ match: 'hmmpress',
                    hint : 'BiG-SCAPE tried to press the Pfam database and could not write to that folder. \n  Run  hmmpress /path/to/Pfam-A.hmm  once by hand, then re-run the pipeline.' ],
                    [ match: 'Missing output file',
                    hint : 'BiG-SCAPE produced no output_files/ folder, which means no BGCs reached it. \n  Check that antiSMASH found clusters — a run where every sample has zero regions gives BiG-SCAPE nothing to cluster.' ],
                    [ match: 'No files found',
                    hint : 'BiG-SCAPE read zero GBK files. Its --include-gbk filter needs "region" in each file name. \n  Check the staged names in gbk_input/ inside the failed task folder.' ],
                ],
                generic   : 'BiG-SCAPE failed. Check --bgc_bigscape_pfam points at a pressed Pfam-A.hmm file and that antiSMASH produced region GBKs.',
            ],
            [
                process   : 'BGCQUAST',
                name      : 'bgc-quast',
                signatures: [],
                generic   : 'bgc-quast failed. Check that the prediction files, the query FASTA, and the QUAST output folder all reached this step.',
            ],
        ]

        def hit = tools.find { leaf == it.process || leaf.startsWith(it.process) }

        def banner = "=".multiply(100)

        if (hit) {
            def sig    = hit.signatures.find { report.contains(it.match) }
            def detail = sig ? sig.hint : hit.generic
            log.error(
                "\n${banner}\n" +
                "[bgc_quast_ppl] The ${hit.name} step failed.\n\n" +
                "  ${detail}\n" +
                "${banner}"
            )
        }
        else {
            log.error(
                "\n${banner}\n" +
                "[bgc_quast_ppl] The pipeline stopped with an error.\n\n" +
                "  See the message above, and open the failing task's .command.err\n" +
                "  file for the full details.\n" +
                "${banner}"
            )
        }

        if (params.bgc_quast_debug && report.trim()) {
            log.error(pink("[bgc_quast_ppl] --bgc_quast_debug: full error report below:\n${report.trim()}"))
        }

        log.error("Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
    }
    catch (Exception e) {
        log.error(pink("[bgc_quast_ppl] error handler failed: ${e}"))
    }
}

//
// True if bgc-quast produced a comparison folder. Checks the output dir, not a channel.
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
// Run ended clean but bgc-quast never ran. Print a red failure notice.
//
def reportNoComparison(monochrome_logs) {
    def red    = monochrome_logs ? '' : "\033[1;31m"
    def reset  = monochrome_logs ? '' : "\033[0m"
    def banner = "=".multiply(100)
    log.error(
        "${red}\n${banner}\n" +
        "[bgc_quast_ppl] Pipeline did NOT complete successfully.\n\n" +
        "  No BGC comparison was produced. bgc-quast never ran, usually because every\n" +
        "  sample was dropped before prediction (for example all contigs were shorter than\n" +
        "  ${params.bgc_mincontiglength} bp, or annotation produced no genes).\n\n" +
        "  Use longer or better assemblies, or lower --bgc_mincontiglength, then run again.\n" +
        "${banner}${reset}"
    )
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastas) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype
    def endedness_ok = metas.collect { meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [metas[0], fastas]
}
