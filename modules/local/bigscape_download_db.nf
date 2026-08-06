process BIGSCAPE_DOWNLOAD_DB {
    tag "pfam"
    label 'process_single'

    // Same image as BIGSCAPE. The hmmpress *binary* is not in this container, but
    // pyhmmer and requests both are (bigscape depends on them), so the download and
    // the press need nothing extra. This is the same press call BiG-SCAPE itself
    // makes in big_scape/hmm/hmmer.py:78-82.
    conda "bioconda::bigscape=2.0.3"
    container "quay.io/biocontainers/bigscape:2.0.3--pyhdfd78af_0"

    output:
    path "pfam_db"     , emit: db
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def url          = params.bgc_bigscape_pfam_url
    def pfam_release = (url =~ /Pfam(\d+\.\d+)/) ? (url =~ /Pfam(\d+\.\d+)/)[0][1] : 'current_release'

    """
    mkdir -p pfam_db

    python3 <<'PYEOF'
    import gzip
    import shutil
    import sys
    from pathlib import Path

    import pyhmmer
    import requests
    from pyhmmer.hmmer import hmmpress
    from pyhmmer.plan7 import HMMFile

    url = "${url}"
    out_dir = Path("pfam_db")
    gz_path = out_dir / "Pfam-A.hmm.gz"
    hmm_path = out_dir / "Pfam-A.hmm"

    print("Downloading " + url, flush=True)
    with requests.get(url, stream=True, timeout=(30, 900)) as response:
        response.raise_for_status()
        with open(gz_path, "wb") as handle:
            for chunk in response.iter_content(chunk_size=1 << 20):
                handle.write(chunk)
    print("Downloaded {:.1f} MiB".format(gz_path.stat().st_size / 1048576), flush=True)

    print("Unpacking", flush=True)
    with gzip.open(gz_path, "rb") as src, open(hmm_path, "wb") as dst:
        shutil.copyfileobj(src, dst, length=1 << 20)
    gz_path.unlink()
    print("Unpacked {:.1f} MiB".format(hmm_path.stat().st_size / 1048576), flush=True)

    print("Pressing", flush=True)
    with HMMFile(hmm_path) as hmm_file:
        hmmpress(hmm_file, hmm_path)

    missing = [e for e in (".h3f", ".h3i", ".h3m", ".h3p")
               if not Path(str(hmm_path) + e).exists()]
    if missing:
        sys.exit("hmmpress did not produce: " + " ".join(missing))
    print("Pressed. Pfam is ready.", flush=True)

    Path("versions.yml").write_text(
        '"${task.process}":\\n'
        '    pfam: ${pfam_release}\\n'
        '    pyhmmer: ' + pyhmmer.__version__ + '\\n'
    )
    PYEOF
    """

    stub:
    def url          = params.bgc_bigscape_pfam_url
    def pfam_release = (url =~ /Pfam(\d+\.\d+)/) ? (url =~ /Pfam(\d+\.\d+)/)[0][1] : 'current_release'

    """
    mkdir -p pfam_db
    touch pfam_db/Pfam-A.hmm
    touch pfam_db/Pfam-A.hmm.h3f
    touch pfam_db/Pfam-A.hmm.h3i
    touch pfam_db/Pfam-A.hmm.h3m
    touch pfam_db/Pfam-A.hmm.h3p

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pfam: ${pfam_release}
        pyhmmer: 0.12.0
    END_VERSIONS
    """
}