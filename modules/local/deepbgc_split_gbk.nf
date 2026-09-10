process DEEPBGC_SPLIT_GBK {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::deepbgc=0.1.31"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/deepbgc:0.1.31--pyhca03a8a_0'
        : 'biocontainers/deepbgc:0.1.31--pyhca03a8a_0'}"

    input:
    tuple val(meta), path(bgc_gbk)

    output:
    tuple val(meta), path("split/*.gbk"), optional: true, emit: gbk
    path "versions.yml"                 ,                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p split

    python3 <<'PYEOF'
    from collections import defaultdict
    from pathlib import Path

    import Bio
    from Bio import SeqIO

    src = Path("${bgc_gbk}")
    out_dir = Path("split")

    # ORDER IS LOAD-BEARING: bgc-quast numbers DeepBGC BGCs by record order, so never sort,
    # reorder or parallelise -- a slip silently puts every family on the wrong BGC.
    counters = defaultdict(int)
    written = 0

    for record in SeqIO.parse(src, "genbank"):
        sequence_id = record.name
        counters[sequence_id] += 1
        number = counters[sequence_id]

        # BiG-SCAPE rejects a cluster feature without this note; --force-gbk cannot rescue it.
        for feature in record.features:
            if feature.type == "cluster":
                feature.qualifiers["note"] = ["Cluster number: %d" % number]

        # "cluster" must stay in the name for BiG-SCAPE's --include-gbk default.
        target = out_dir / ("%s_cluster_%d.gbk" % (sequence_id, number))
        SeqIO.write(record, target, "genbank")
        written += 1

    # No BGCs is a result, not a failure, so the emit is optional.
    if written:
        print("Split %d record(s) into %s/" % (written, out_dir), flush=True)
    else:
        print("No BGCs in %s -- nothing to split for sample '${meta.id}'." % src, flush=True)

    Path("versions.yml").write_text(
        '"${task.process}":\\n'
        '    biopython: ' + Bio.__version__ + '\\n'
    )
    PYEOF
    """

    stub:
    """
    mkdir -p split
    touch split/CONTIG_1_cluster_1.gbk
    touch split/CONTIG_1_cluster_2.gbk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: 1.83
    END_VERSIONS
    """
}