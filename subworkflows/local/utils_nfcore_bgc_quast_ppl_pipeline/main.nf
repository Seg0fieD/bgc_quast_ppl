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
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // antiSMASH: minimal and full are mutually exclusive.
    //
    validateAntismashMode()

    //
    // Compare-to-reference: enforce samplesheet shape correction
    //
    if (params.bgc_quast_mode == 'compare-to-reference') {
        validateReferenceSamplesheet(input)
    }

    //
    // Pre-run environment checks: paths, databases, Docker.
    // Collects all blocking problems and reports them together.
    //
    validatePreRunEnvironment(input)

    //
    // Create channel from input file provided through params.input
    //
    Channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
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

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

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

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
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
// compare-to-reference samplesheet check: sample/fasta/type present, no empty cells,
// exactly one reference row (type r/R). Simple comma split (no quoted-comma support).
//
def validateReferenceSamplesheet(input) {
    def lines = file(input).readLines().findAll { it.trim() }
    if (lines.size() < 2) {
        error("[bgc_quast_ppl] compare-to-reference: samplesheet has no data rows.")
    }

    def header = lines[0].split(',', -1).collect { it.trim() }
    ['sample', 'fasta', 'type'].each { col ->
        if (!header.contains(col)) {
            error("[bgc_quast_ppl] compare-to-reference needs a '${col}' column. Found: ${header.join(', ')}")
        }
    }

    def si = header.indexOf('sample')
    def fi = header.indexOf('fasta')
    def ti = header.indexOf('type')
    def ref_count = 0

    lines[1..-1].eachWithIndex { line, idx ->
        def cells  = line.split(',', -1)
        def rownum = idx + 2
        [si, fi, ti].each { ci ->
            if (ci >= cells.size() || cells[ci].trim() == '') {
                error("[bgc_quast_ppl] compare-to-reference: empty cell in row ${rownum}. sample, fasta and type must all be filled.")
            }
        }
        def t = cells[ti].trim().toLowerCase()
        if (!(t in ['q', 'r'])) {
            error("[bgc_quast_ppl] compare-to-reference: row ${rownum} type='${cells[ti].trim()}' is invalid. Use q/Q (query) or r/R (reference).")
        }
        if (t == 'r') { ref_count++ }
    }

    if (ref_count != 1) {
        error("[bgc_quast_ppl] compare-to-reference needs exactly one reference row (type r/R). Found ${ref_count}.")
    }
}

//
// antiSMASH minimal vs full. Minimal is the default (runs when neither flag is given).
// Full runs only when explicitly requested. Passing BOTH flags is an error.
//
def validateAntismashMode() {
    def cli = workflow.commandLine ?: ''
    def minimal_typed = cli.contains('--bgc_antismash_minimal')
    def full_typed    = cli.contains('--bgc_antismash_full')

    if (minimal_typed && full_typed) {
        error("[bgc_quast_ppl] --bgc_antismash_minimal and --bgc_antismash_full cannot both be set. Minimal is the default; pass --bgc_antismash_full only if you want the full analysis.")
    }
}

//
// Error Handling: Pre-run environment checks. 
//   problems -> block the run (collected > display together > then halt)
//   warnings -> printed, run continues
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
        error("[bgc_quast_ppl] Cannot start. Please fix:\n${msg}")
    }
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
