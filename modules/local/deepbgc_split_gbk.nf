process DEEPBGC_SPLIT_GBK {
    tag "${meta.id}"
    label 'process_single'

    // Same image as DEEPBGC_PIPELINE. biopython is a deepbgc dependency, so nothing extra.
    conda "bioconda::deepbgc=0.1.31"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/deepbgc:0.1.31--pyhca03a8a_0'
        : 'biocontainers/deepbgc:0.1.31--pyhca03a8a_0'}"

    input:
    tuple val(meta), path(bgc_gbk)

    output:
    tuple val(meta), path("split/*.gbk"), emit: gbk
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p split

    python3 <<'PYEOF'
    import sys
    from collections import defaultdict
    from pathlib import Path

    import Bio
    from Bio import SeqIO

    src = Path("${bgc_gbk}")
    out_dir = Path("split")

    # DeepBGC writes one GBK holding every BGC. BiG-SCAPE reads only the first record of a
    # multi-record file (big_scape/genbank/gbk.py:511-517), so each record becomes its own file.
    #
    # ORDER IS LOAD-BEARING. bgc-quast numbers DeepBGC BGCs with a running per-sequence counter
    # over the rows of .bgc.tsv (genome_mining_parser.py:266-268). Nothing in the GBK carries
    # that number, so it is rebuilt here from record order, which was verified to match TSV row
    # order on all three test samples. Do not sort, reorder or parallelise this loop -- a slip
    # puts every gene cluster family on the wrong BGC, with no error and no warning.
    counters = defaultdict(int)
    written = 0

    for record in SeqIO.parse(src, "genbank"):
        sequence_id = record.name
        counters[sequence_id] += 1
        number = counters[sequence_id]

        # BiG-SCAPE's AS4 path requires this note on the cluster feature
        # (big_scape/genbank/region.py:223-231). DeepBGC does not write one, so every file
        # is rejected with InvalidGBKError without it. --force-gbk cannot help: that fallback
        # only fires when there is no cluster feature at all.
        for feature in record.features:
            if feature.type == "cluster":
                feature.qualifiers["note"] = ["Cluster number: %d" % number]

        # "cluster" must stay in the name: BiG-SCAPE's --include-gbk defaults to cluster,region.
        # The <sample_label>_ prefix is added later, by the BIGSCAPE staging loop.
        target = out_dir / ("%s_cluster_%d.gbk" % (sequence_id, number))
        SeqIO.write(record, target, "genbank")
        written += 1

    if not written:
        sys.exit("No records found in " + str(src))
    print("Split %d record(s) into %s/" % (written, out_dir), flush=True)

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