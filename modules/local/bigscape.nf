process BIGSCAPE {
    tag "${prefix}"
    label 'process_high'

    // Tested with the container; conda is kept for users without Docker.
    conda "bioconda::bigscape=2.0.3"
    container "quay.io/biocontainers/bigscape:2.0.3--pyhdfd78af_0"

    input:
    val  prefix                  // tool name: antismash, gecco or deepbgc
    val  names                   // <sample>_<GBK file name>, one per GBK
    path gbks, stageAs: 'raw*/*' // region or cluster GBKs of all samples
    path pfam_dir                // .hmm file with its .h3f .h3i .h3m .h3p
    val  pfam_name               // file name of the .hmm inside pfam_dir

    output:
    // One folder per tool, published under bigscape/ by conf/modules.config.
    tuple val(prefix), path("${prefix}")                             , emit: results
    path "${prefix}/output_files/**/*_clustering_c*.tsv"             , emit: clustering, optional: true
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def name_list = names instanceof List ? names : [names]
    def gbk_list  = gbks  instanceof List ? gbks  : [gbks]


    // GBKs symlinked under their paired names, sorted for a stable task hash.
    // Names keep <sample>_ for bgc-quast, .region/_cluster_ for BiG-SCAPE.
    def stage_cmds = (0..<gbk_list.size()).toList()
        .sort { i -> name_list[i] }
        .collect { i -> "ln -s \"\$WORKDIR/${gbk_list[i]}\" \"\$WORKDIR/gbk_input/${name_list[i]}\"" }
        .join('\n    ')

    """
    WORKDIR=\$PWD
    mkdir -p \$WORKDIR/gbk_input
    ${stage_cmds}

    # Fixed hash seed, so BiG-SCAPE loads its inputs in one order every run.
    export PYTHONHASHSEED=0

    # Seeds numpy in every Python process, so affinity propagation is stable.
    echo 'import numpy' > \$WORKDIR/sitecustomize.py
    echo 'numpy.random.seed(0)' >> \$WORKDIR/sitecustomize.py
    export PYTHONPATH="\$WORKDIR\${PYTHONPATH:+:\$PYTHONPATH}"

    bigscape cluster \\
        -i \$WORKDIR/gbk_input \\
        -o \$WORKDIR/${prefix} \\
        -p \$WORKDIR/${pfam_dir}/${pfam_name} \\
        -c ${task.cpus} \\
        -l ${prefix} \\
        ${args}

    BIGSCAPE_VERSION=\$( { bigscape --version 2>&1 || true; } | tail -n 1 )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigscape: \${BIGSCAPE_VERSION:-2.0.3}
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${prefix}/output_files/${prefix}_2026-01-01_00-00-00_c0.3/mix
    touch ${prefix}/output_files/${prefix}_2026-01-01_00-00-00_c0.3/mix/mix_clustering_c0.3.tsv
    touch ${prefix}/index.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bigscape: 2.0.3
    END_VERSIONS
    """
}