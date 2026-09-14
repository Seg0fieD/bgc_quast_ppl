process DEEPBGC_SPLIT_GBK {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::deepbgc=0.1.31"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/deepbgc:0.1.31--pyhca03a8a_0'
        : 'biocontainers/deepbgc:0.1.31--pyhca03a8a_0'}"

    input:
    tuple val(meta), path(bgc_gbk), path(bgc_tsv)

    output:
    tuple val(meta), path("split/*.gbk"), optional: true, emit: gbk
    path "versions.yml"                 ,                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p split

    python3 <<'PYEOF'
    import csv
    import sys
    from collections import defaultdict
    from pathlib import Path

    import Bio
    from Bio import SeqIO

    src = Path("${bgc_gbk}")
    tsv = Path("${bgc_tsv}")
    out_dir = Path("split")

    # Values pandas turns into NaN when bgc-quast reads this same TSV.
    PANDAS_NA = {
        "", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan",
        "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a", "nan", "null",
    }

    def fail(message):
        sys.exit("DEEPBGC_SPLIT_GBK [${meta.id}]: " + message)

    def survives_pandas(name):
        if name in PANDAS_NA:
            return False
        for cast in (int, float):
            try:
                return str(cast(name)) == name
            except ValueError:
                continue
        return True

    with tsv.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\\t"))

    if rows:
        for column in ("sequence_id", "bgc_candidate_id"):
            if column not in rows[0]:
                fail("%s has no '%s' column." % (tsv.name, column))

    # bgc-quast numbers each BGC by its position among the rows of this TSV, so the
    # number is taken from here and never from the GBK, whose LOCUS can be truncated.
    counter = defaultdict(int)
    expected = {}
    order = []
    for line, row in enumerate(rows, start=2):
        sequence_id = (row["sequence_id"] or "").strip()
        candidate_id = (row["bgc_candidate_id"] or "").strip()
        if not sequence_id or not candidate_id:
            fail("%s line %d has a blank sequence_id or bgc_candidate_id." % (tsv.name, line))
        if any(c in sequence_id for c in "/\\\\") or sequence_id != "".join(sequence_id.split()):
            fail("sequence_id %r cannot be used in a file name." % sequence_id)
        if not survives_pandas(sequence_id):
            fail("sequence_id %r changes type when bgc-quast reads the TSV." % sequence_id)
        if candidate_id in expected:
            fail("bgc_candidate_id %r appears twice in %s." % (candidate_id, tsv.name))
        counter[sequence_id] += 1
        expected[candidate_id] = (sequence_id, counter[sequence_id])
        order.append(candidate_id)

    seen = []
    for position, record in enumerate(SeqIO.parse(src, "genbank"), start=1):
        clusters = [f for f in record.features if f.type == "cluster"]
        if len(clusters) != 1:
            fail("record %d of %s has %d cluster features, expected 1." % (position, src.name, len(clusters)))
        candidate_ids = clusters[0].qualifiers.get("bgc_candidate_id", [])
        if len(candidate_ids) != 1:
            fail("record %d of %s carries no single bgc_candidate_id qualifier." % (position, src.name))
        candidate_id = candidate_ids[0].strip()
        if candidate_id not in expected:
            fail("%s is in %s but not in %s." % (candidate_id, src.name, tsv.name))
        if candidate_id in seen:
            fail("%s appears in %s more than once." % (candidate_id, src.name))
        seen.append(candidate_id)

        sequence_id, number = expected[candidate_id]
        # BiG-SCAPE rejects a cluster feature without this note; --force-gbk cannot rescue it.
        clusters[0].qualifiers["note"] = ["Cluster number: %d" % number]
        target = out_dir / ("%s_cluster_%d.gbk" % (sequence_id, number))
        if target.exists():
            fail("%s would be written twice." % target.name)
        SeqIO.write(record, target, "genbank")

    written = set(seen)
    missing = [c for c in order if c not in written]
    if missing:
        fail("%d TSV row(s) have no GBK record, first is %s." % (len(missing), missing[0]))
    if len(seen) != len(expected):
        fail("wrote %d file(s) for %d TSV row(s)." % (len(seen), len(expected)))

    if expected:
        print("Split %d record(s) into %s/ (GBK order matches TSV order: %s)"
              % (len(seen), out_dir, "yes" if seen == order else "no"), flush=True)
    else:
        print("No BGCs for sample '${meta.id}' -- nothing to split.", flush=True)

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