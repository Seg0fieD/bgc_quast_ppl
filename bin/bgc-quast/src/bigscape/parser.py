"""Read BiG-SCAPE cluster output and map each BGC to its gene cluster family (GCF).

BiG-SCAPE names a BGC by its input file name. The pipeline stages every antiSMASH
region GBK as `<sample_label>_<original_filename>`, so this module reverses that to
rebuild bgc-quast's own BGC id and keys results on (sample_label, bgc_id).

Nothing here raises into the main run. A missing or unreadable folder yields {}.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

# {cutoff: {(sample_label, bgc_id): family_id}}
FamilyMap = Dict[str, Dict[Tuple[str, str], str]]

# Column in the clustering TSV that holds the input GBK file name.
GBK_COLUMN = "GBK"
FAMILY_COLUMN = "Family"

REGION_MARKER = ".region"


def normalize_cutoff(value) -> str:
    """Normalise a cutoff to a canonical string: 0.30 -> '0.3', '0.700' -> '0.7'.

    Used so `--bigscape-cutoff 0.3` matches a folder named `..._c0.30`.
    Returns the input as a stripped string if it is not a number.
    """
    try:
        return f"{float(value):g}"
    except (TypeError, ValueError):
        return str(value).strip()


def strip_sample_prefix(stem: str, sample_labels: List[str]) -> Optional[Tuple[str, str]]:
    """Split `<sample_label>_<rest>` into (sample_label, rest).

    Matches the LONGEST known label, so labels like `assembly` and `assembly_10`
    do not collide. Returns None if no known label matches.
    """
    for label in sorted(sample_labels, key=len, reverse=True):
        prefix = f"{label}_"
        if stem.startswith(prefix):
            return label, stem[len(prefix):]
    return None


def bgc_id_from_record(record_name: str) -> Optional[str]:
    """Turn `CONTIG_2.region001` into `CONTIG_2.1`.

    bgc-quast builds its id as sequence_id + "." + the raw region_number qualifier,
    which is not zero padded (genome_mining_parser.py:144). So the padding is stripped.
    """
    if REGION_MARKER not in record_name:
        return None
    sequence_id, _, region_part = record_name.rpartition(REGION_MARKER)
    if not sequence_id:
        return None
    digits = "".join(c for c in region_part if c.isdigit())
    if not digits:
        return None
    return f"{sequence_id}.{int(digits)}"


def find_clustering_files(bigscape_output_dir: Path) -> Dict[str, Path]:
    """Locate one `mix` clustering TSV per cutoff.

    Folder names carry a run label and a timestamp, e.g.
    `output_files/bigscape_2026-08-04_12-52-26_c0.3/mix/mix_clustering_c0.30.tsv`,
    so nothing is built from a label or a timestamp. The cutoff is read off the
    folder name, after the last `_c`.

    Returns {normalised_cutoff: path}. Empty if nothing is found.
    """
    root = Path(bigscape_output_dir) / "output_files"
    if not root.is_dir():
        return {}

    found: Dict[str, Path] = {}
    for cutoff_dir in sorted(root.glob("*_c*")):
        if not cutoff_dir.is_dir():
            continue  # skips *_full.network and other stray files

        _, marker, raw_cutoff = cutoff_dir.name.rpartition("_c")
        if not marker or not raw_cutoff:
            continue

        matches = sorted(cutoff_dir.glob("mix/*_clustering_c*.tsv"))
        if not matches:
            continue  # this cutoff produced no mix bin; skip it, do not fail

        found[normalize_cutoff(raw_cutoff)] = matches[0]

    return found


def parse_clustering_file(
    tsv_path: Path,
    sample_labels: List[str],
    log=None,
) -> Dict[Tuple[str, str], str]:
    """Parse one clustering TSV into {(sample_label, bgc_id): family_id}."""
    table = pd.read_csv(tsv_path, sep="\t", dtype=str)

    missing = [c for c in (GBK_COLUMN, FAMILY_COLUMN) if c not in table.columns]
    if missing:
        if log:
            log.warning(
                f"BiG-SCAPE clustering file {tsv_path} is missing column(s): "
                f"{', '.join(missing)}. Skipping it."
            )
        return {}

    families: Dict[Tuple[str, str], str] = {}
    unmatched = 0

    for gbk_name, family_id in zip(table[GBK_COLUMN], table[FAMILY_COLUMN]):
        if not isinstance(gbk_name, str) or not isinstance(family_id, str):
            continue

        stem = gbk_name.strip()
        if stem.endswith(".gbk"):
            stem = stem[: -len(".gbk")]

        split = strip_sample_prefix(stem, sample_labels)
        if split is None:
            unmatched += 1
            continue

        sample_label, record_name = split
        bgc_id = bgc_id_from_record(record_name)
        if bgc_id is None:
            unmatched += 1
            continue

        families[(sample_label, bgc_id)] = family_id.strip()

    if unmatched and log:
        log.warning(
            f"{unmatched} row(s) in {tsv_path.name} did not match a known sample "
            f"label or had no region number. Known labels: {', '.join(sample_labels)}"
        )

    return families


def parse_bigscape(
    bigscape_output_dir,
    sample_labels: List[str],
    log=None,
) -> FamilyMap:
    """Read a whole BiG-SCAPE output folder.

    Args:
        bigscape_output_dir: the folder passed to BiG-SCAPE as `-o`.
        sample_labels: bgc-quast display labels, used to strip filename prefixes.
        log: optional bgc-quast Logger.

    Returns:
        {cutoff: {(sample_label, bgc_id): family_id}}. Empty on any failure.
    """
    try:
        directory = Path(bigscape_output_dir)
        if not directory.is_dir():
            if log:
                log.warning(f"BiG-SCAPE output directory not found: {directory}")
            return {}

        clustering_files = find_clustering_files(directory)
        if not clustering_files:
            if log:
                log.warning(
                    f"No BiG-SCAPE clustering files under {directory}/output_files. "
                    "Expected output_files/*_c*/mix/*_clustering_c*.tsv"
                )
            return {}

        result: FamilyMap = {}
        for cutoff, tsv_path in clustering_files.items():
            try:
                parsed = parse_clustering_file(tsv_path, sample_labels, log=log)
            except Exception as e:
                if log:
                    log.warning(f"Failed to read {tsv_path}: {e}")
                continue
            if parsed:
                result[cutoff] = parsed

        if log and result:
            log.info(
                "BiG-SCAPE results loaded for cutoff(s): "
                f"{', '.join(sorted(result))}"
            )
        return result

    except Exception as e:
        if log:
            log.warning(f"Failed to read BiG-SCAPE output: {e}")
        return {}


def select_cutoff(families: FamilyMap, requested) -> Dict[Tuple[str, str], str]:
    """Pick one cutoff's family map.

    Unlike the rest of this module this DOES raise, because a cutoff the user asked
    for that was never run is a real mistake, and silently falling back would hide it.
    """
    if not families:
        return {}

    key = normalize_cutoff(requested)
    if key not in families:
        available = ", ".join(sorted(families))
        raise ValueError(
            f"BiG-SCAPE cutoff '{requested}' was not found in the output. "
            f"Available cutoff(s): {available}"
        )
    return families[key]
