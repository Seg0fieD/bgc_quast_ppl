process BGCQUAST {
    tag "${meta.id}"
    label 'process_low'

    // bgc-quast has no biocontainer and is not a packaged tool. It is run in place from
    // bin/bgc-quast/ (option 1). conda provides its Python deps (see bin/bgc-quast/environment.yml).
    conda "${projectDir}/bin/bgc-quast/environment.yml"

    input:
    tuple val(meta), path(mining_results), path(genome), path(quast_dir), path(reference_mining), path(reference_genome)

    output:
    // bgc-quast writes report.html/tsv/txt and bgc-quast.log straight into the output dir.
    // We give it a fixed internal dir and publish its CONTENTS, so the final folder name is
    // decided only by publishDir (conf/modules.config) -> bgc_quast/<mode>/<meta.leaf>.
    tuple val(meta), path("bgcquast_out/*")          , emit: results
    tuple val(meta), path("bgcquast_out/report.tsv") , optional: true, emit: tsv
    tuple val(meta), path("bgcquast_out/report.html"), optional: true, emit: html
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // --names is wired from the subworkflow (which tools/samples are present), carried on meta.
    def names_arg            = meta.bgcquast_names    ? "--names ${meta.bgcquast_names}"                    : ''
    // File-dependent flags are built here from staged inputs; flat tuning flags come via ext.args (modules.config).
    def genome_list = genome ? (genome instanceof List ? genome : [genome]) : []
    // bgc-quast matches each --genome file to its mining result BY FILENAME LABEL.
    // The prepped query FASTA is named "<id>_long.fasta", so strip "_long" and symlink to
    // "<id>.fasta" so the genome label matches the mining-result label.
    def genome_renames = genome_list.collect { "ln -sf \$WORKDIR/${it.name} \$WORKDIR/renamed/" + it.name.replaceFirst(/_long\./, '.') }.join('\n    ')
    def genome_arg     = genome_list ? "--genome " + genome_list.collect { "\$WORKDIR/renamed/" + it.name.replaceFirst(/_long\./, '.') }.join(' ') : ''
    def quast_arg            = quast_dir        ? "--quast-output-dir \$WORKDIR/${quast_dir}"               : ''
    def reference_mining_arg = reference_mining ? "--reference-mining-result \$WORKDIR/${reference_mining}" : ''
    def reference_genome_arg = reference_genome ? "--reference-genome \$WORKDIR/${reference_genome}"        : ''

    """
    # Capture the task dir, then run the tool from its own folder so its `from src.*`
    # imports and any bundled-config reads resolve. All staged paths and the output are
    # given as absolute paths because we change directory.
    WORKDIR=\$PWD
    mkdir -p \$WORKDIR/bgcquast_out
    mkdir -p \$WORKDIR/renamed
    ${genome_renames}

    cd ${projectDir}/bin/bgc-quast

    python3 bgc-quast.py \\
        ${mining_results.collect { "\$WORKDIR/${it}" }.join(' \\\n        ')} \\
        ${genome_arg} \\
        ${quast_arg} \\
        ${reference_mining_arg} \\
        ${reference_genome_arg} \\
        ${names_arg} \\
        --threads ${task.cpus} \\
        --output-dir \$WORKDIR/bgcquast_out \\
        ${args}

    cat <<-END_VERSIONS > \$WORKDIR/versions.yml
    "${task.process}":
        bgc-quast: \$(cat ${projectDir}/bin/bgc-quast/VERSION.txt)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p bgcquast_out
    touch bgcquast_out/report.tsv
    touch bgcquast_out/report.html
    touch bgcquast_out/report.txt
    touch bgcquast_out/bgc-quast.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bgc-quast: \$(cat ${projectDir}/bin/bgc-quast/VERSION.txt)
    END_VERSIONS
    """
}
