/*
    Nextflow pipeline utilities: version print, parameter dump to JSON and
    the conda channel check.
*/
workflow UTILS_NEXTFLOW_PIPELINE {
    take:
    print_version        // boolean: print version
    dump_parameters      // boolean: dump parameters
    outdir               // path   : results directory
    check_conda_channels // boolean: check conda channels

    main:

    // Pipeline name and version on --version, then exit.
    if (print_version) {
        log.info("${workflow.manifest.name} ${getWorkflowVersion()}")
        System.exit(0)
    }

    // Parameter dump to JSON in the output directory.
    if (dump_parameters && outdir) {
        dumpParametersToJSON(outdir)
    }

    // Conda channel order check, warning on a missing or misordered channel.
    if (check_conda_channels) {
        checkCondaChannels()
    }

    emit:
    dummy_emit = true
}


// Workflow version string: manifest version and short commit id.
def getWorkflowVersion() {
    def version_string = "" as String
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }

    if (workflow.commitId) {
        def git_shortsha = workflow.commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }

    return version_string
}

// Run parameters as a timestamped JSON file in <outdir>/pipeline_info.
def dumpParametersToJSON(outdir) {
    def timestamp = new java.util.Date().format('yyyy-MM-dd_HH-mm-ss')
    def filename  = "params_${timestamp}.json"
    def temp_pf   = new File(workflow.launchDir.toString(), ".${filename}")
    def jsonStr   = groovy.json.JsonOutput.toJson(params)
    temp_pf.text  = groovy.json.JsonOutput.prettyPrint(jsonStr)

    nextflow.extension.FilesEx.copyTo(temp_pf.toPath(), "${outdir}/pipeline_info/params_${timestamp}.json")
    temp_pf.delete()
}

// Conda channel check: conda-forge then bioconda, else a warning.
def checkCondaChannels() {
    def parser = new org.yaml.snakeyaml.Yaml()
    def channels = []
    try {
        def config = parser.load("conda config --show channels".execute().text)
        channels = config.channels
    }
    catch (NullPointerException e) {
        log.debug(e)
        log.warn("Could not verify conda channel configuration.")
        return null
    }
    catch (IOException e) {
        log.debug(e)
        log.warn("Could not verify conda channel configuration.")
        return null
    }

    // Required channels, highest priority first; all must be present.
    def required_channels_in_order = ['conda-forge', 'bioconda']
    def channels_missing = ((required_channels_in_order as Set) - (channels as Set)) as Boolean

    // Channel order check.
    def channel_priority_violation = required_channels_in_order != channels.findAll { ch -> ch in required_channels_in_order }

    if (channels_missing | channel_priority_violation) {
        log.warn """\
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            There is a problem with your Conda configuration!
            You will need to set-up the conda-forge and bioconda channels correctly.
            Please refer to https://bioconda.github.io/
            The observed channel order is
            ${channels}
            but the following channel order is required:
            ${required_channels_in_order}
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        """.stripIndent(true)
    }
}
