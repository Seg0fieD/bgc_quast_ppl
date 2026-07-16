process BGCQUAST {
    tag "${meta.id}"
    label 'process_low'

    // Run in place from bin/bgc-quast/; no biocontainer exists.
    conda "${projectDir}/bin/bgc-quast/environment.yml"

    input:
    tuple val(meta), path(mining_results), path(genome), path(quast_dir), path(reference_mining), path(reference_genome)

    output:
    // Final published folder is set by publishDir (conf/modules.config).
    tuple val(meta), path("bgcquast_out/*")          , emit: results
    tuple val(meta), path("bgcquast_out/report.tsv") , optional: true, emit: tsv
    tuple val(meta), path("bgcquast_out/report.html"), optional: true, emit: html
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def names_arg = meta.bgcquast_names ? "--names ${meta.bgcquast_names}" : ''
    def genome_list = genome ? (genome instanceof List ? genome : [genome]) : []
    // Symlink "<id>_long.fasta" to "<id>.fasta" so --genome matches the mining-result label.
    def genome_renames = genome_list.collect { "ln -sf \$WORKDIR/${it.name} \$WORKDIR/renamed/" + it.name.replaceFirst(/_long\./, '.') }.join('\n    ')
    def genome_arg     = genome_list ? "--genome " + genome_list.collect { "\$WORKDIR/renamed/" + it.name.replaceFirst(/_long\./, '.') }.join(' ') : ''
    def quast_arg            = quast_dir        ? "--quast-output-dir \$WORKDIR/${quast_dir}"               : ''
    def reference_mining_arg = reference_mining ? "--reference-mining-result \$WORKDIR/${reference_mining}" : ''
    def reference_genome_arg = reference_genome ? "--reference-genome \$WORKDIR/${reference_genome}"        : ''

    """
    # Run from bin/bgc-quast/ so `from src.*` imports resolve; staged paths passed as absolute.
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