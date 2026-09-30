/*
    Parameter summary and validation through the nf-schema plugin.
*/

include { paramsSummaryLog   } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'

workflow UTILS_NFSCHEMA_PLUGIN {

    take:
    input_workflow      // workflow: object nf-schema reads metadata from
    validate_params     // boolean:  validate the parameters
    parameters_schema   // string:   params schema path matching
                        //           validation.parametersSchema, or empty
    main:
    // Summary of the parameters that differ from the schema defaults.
    if(parameters_schema) {
        log.info paramsSummaryLog(input_workflow, parameters_schema:parameters_schema)
    } else {
        log.info paramsSummaryLog(input_workflow)
    }

    // Validation against parameters_schema, else the configured schema.
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

