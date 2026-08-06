process BIGSCAPE {
    tag "bigscape"
    label 'process_high'

    // Container is the tested route; conda is here so a conda user is not blocked.
    conda "bioconda::bigscape=2.0.3"
    container "quay.io/biocontainers/bigscape:2.0.3--pyhdfd78af_0"

    input:
    val  names                              // <sample_label>_<original_gbk_filename>, one per gbk
    path gbks, stageAs: 'raw*/*'            // antiSMASH region GBKs, all samples
    path pfam_dir                           // folder holding the .hmm plus .h3f .h3i .h3m .h3p
    val  pfam_name                          // basename of the .hmm inside pfam_dir

    output:
    path "bigscape"                                              , emit: results
    path "bigscape/output_files/**/*_clustering_c*.tsv"           , emit: clustering, optional: true
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def name_list = names instanceof List ? names : [names]
    def gbk_list  = gbks  instanceof List ? gbks  : [gbks]

    if (name_list.size() != gbk_list.size()) {
        error(
            "[bgc_quast_ppl] BIGSCAPE: got ${name_list.size()} names for ${gbk_list.size()} GBK files.\n" +
            "                These two lists must be built from the same collect() and pair by index."
        )
    }

    // Nextflow has no per-file rename. Pair names to staged paths by index and symlink
    // into a flat gbk_input/. The <sample>_ prefix is the join key bgc-quast reverses,
    // and ".region" must survive so BiG-SCAPE's --include-gbk filter accepts the file.
    def stage_cmds = (0..<gbk_list.size())
        .collect { i -> "ln -s \"\$WORKDIR/${gbk_list[i]}\" \"\$WORKDIR/gbk_input/${name_list[i]}\"" }
        .join('\n    ')

    """
    WORKDIR=\$PWD
    mkdir -p \$WORKDIR/gbk_input
    ${stage_cmds}

    bigscape cluster \\
        -i \$WORKDIR/gbk_input \\
        -o \$WORKDIR/bigscape \\
        -p \$WORKDIR/${pfam_dir}/${pfam_name} \\
        -c ${task.cpus} \\
        -l bigscape \\
        ${args}

    BIGSCAPE_VERSION=\$( { bigscape --version 2>&1 || true; } | tail -n 1 )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigscape: \${BIGSCAPE_VERSION:-2.0.3}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p bigscape/output_files/bigscape_2026-01-01_00-00-00_c0.3/mix
    touch bigscape/output_files/bigscape_2026-01-01_00-00-00_c0.3/mix/mix_clustering_c0.3.tsv
    touch bigscape/index.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigscape: 2.0.3
    END_VERSIONS
    """
}