#!/usr/bin/env nextflow
/*
    bgc_quast_ppl: BGC prediction and comparison, from input contigs through
    annotation and antiSMASH, DeepBGC and GECCO to the bgc-quast report.
*/

include { BGC_QUAST_PPL           } from './workflows/bgc_quast_ppl'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_bgc_quast_ppl_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_bgc_quast_ppl_pipeline'

workflow NFCORE_BGC_QUAST_PPL {
    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    BGC_QUAST_PPL(
        samplesheet
    )

    emit:
    bgcquast_runs = BGC_QUAST_PPL.out.bgcquast_runs
}

workflow {
    main:
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    NFCORE_BGC_QUAST_PPL(
        PIPELINE_INITIALISATION.out.samplesheet
    )

    // Printed here, after the workflow is built, so it follows the step list.
    def pre_run_warning = PIPELINE_INITIALISATION.out.pre_run_warning.val
    if (pre_run_warning) {
        log.info(pre_run_warning)
    }

    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_BGC_QUAST_PPL.out.bgcquast_runs
    )
}
