include { paramsSummaryLog   } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'

workflow UTILS_NFSCHEMA_PLUGIN {

    ttake:
    input_workflow      // workflow: object nf-schema reads metadata from
    validate_params     // boolean:  validate the parameters
    parameters_schema   // string:   path to the params JSON schema; must match validation.parametersSchema.
                        //           empty = use the configured schema or "${projectDir}/nextflow_schema.json".
                        //           should not be empty for meta pipelines
    main:

    //
    // Print parameter summary that differ from the default given in the JSON schema
    //
    if(parameters_schema) {
        log.info paramsSummaryLog(input_workflow, parameters_schema:parameters_schema)
    } else {
        log.info paramsSummaryLog(input_workflow)
    }

    //
    // Validate params against nextflow_schema.json (or validation.parametersSchema).
    //
    if(validate_params) {
        if(parameters_schema) {
            validateParameters(parameters_schema:parameters_schema)
        } else {
            validateParameters()
        }
    }

    emit:
    dummy_emit = true
}

