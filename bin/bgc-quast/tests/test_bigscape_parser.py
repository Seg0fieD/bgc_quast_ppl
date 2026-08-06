"""Tests for src/bigscape/parser.py.

These build a fake BiG-SCAPE output tree, so they need no real run, no Pfam
and no BiG-SCAPE install.
"""

from pathlib import Path

import pytest

from src.bigscape.parser import (
    bgc_id_from_record,
    find_clustering_files,
    normalize_cutoff,
    parse_bigscape,
    parse_clustering_file,
    select_cutoff,
    strip_sample_prefix,
)

HEADER = "Record\tGBK\tRecord_Type\tRecord_Number\tCC\tFamily"

# Three samples, same genome, six BGCs each -> six families of three.
ROWS = [
    ("reference_CONTIG_2.region001.gbk", "FAM_00001"),
    ("assembly_10_CONTIG_2.region001.gbk", "FAM_00001"),
    ("assembly_20_CONTIG_2.region001.gbk", "FAM_00001"),
    ("reference_CONTIG_6.region001.gbk", "FAM_00002"),
    ("assembly_10_CONTIG_6.region001.gbk", "FAM_00002"),
    ("assembly_20_CONTIG_6.region001.gbk", "FAM_00002"),
]

LABELS = ["reference", "assembly_10", "assembly_20"]


def write_clustering_tsv(path: Path, rows=ROWS) -> None:
    """Write a minimal clustering TSV with the real column names."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [HEADER]
    for gbk, family in rows:
        record = gbk[: -len(".gbk")]
        lines.append(f"{record}\t{gbk}\tregion\t1\tCC_1\t{family}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.fixture
def bigscape_out(tmp_path: Path) -> Path:
    """A realistic BiG-SCAPE output folder with three cutoffs."""
    root = tmp_path / "bigscape_out" / "output_files"
    stamp = "bigscape_2026-08-04_12-52-26"

    for cutoff in ("0.3", "0.5", "0.7"):
        cutoff_dir = root / f"{stamp}_c{cutoff}"
        # BiG-SCAPE zero-pads the file name but not the folder name.
        write_clustering_tsv(cutoff_dir / "mix" / f"mix_clustering_c{cutoff}0.tsv")
        # Non-mix bins must be ignored.
        (cutoff_dir / "NRPS").mkdir(parents=True, exist_ok=True)
        (cutoff_dir / "record_annotations.tsv").write_text("x\n", encoding="utf-8")

    # Stray files at the same level must be ignored.
    (root / f"{stamp}_full.network").write_text("x\n", encoding="utf-8")
    (root / "stray_c0.9.tsv").write_text("x\n", encoding="utf-8")

    return tmp_path / "bigscape_out"


# --- small helpers ---------------------------------------------------------


def test_normalize_cutoff_drops_padding():
    assert normalize_cutoff("0.30") == "0.3"
    assert normalize_cutoff(0.7) == "0.7"
    assert normalize_cutoff("0.700") == "0.7"
    assert normalize_cutoff(0.3) == normalize_cutoff("0.30")


def test_normalize_cutoff_passes_through_non_numbers():
    assert normalize_cutoff(" abc ") == "abc"


def test_strip_sample_prefix_takes_the_longest_match():
    labels = ["assembly", "assembly_10"]
    assert strip_sample_prefix("assembly_10_CONTIG_2.region001", labels) == (
        "assembly_10",
        "CONTIG_2.region001",
    )
    assert strip_sample_prefix("assembly_CONTIG_2.region001", labels) == (
        "assembly",
        "CONTIG_2.region001",
    )


def test_strip_sample_prefix_returns_none_for_unknown_label():
    assert strip_sample_prefix("mystery_CONTIG_2.region001", LABELS) is None


def test_bgc_id_strips_zero_padding():
    assert bgc_id_from_record("CONTIG_2.region001") == "CONTIG_2.1"
    assert bgc_id_from_record("CONTIG_2.region012") == "CONTIG_2.12"
    assert bgc_id_from_record("NZ_CP069563.1.region003") == "NZ_CP069563.1.3"


def test_bgc_id_returns_none_without_a_region():
    assert bgc_id_from_record("CONTIG_2") is None
    assert bgc_id_from_record("CONTIG_2.regionXYZ") is None


# --- discovery -------------------------------------------------------------


def test_find_clustering_files_reads_cutoff_off_the_folder_name(bigscape_out: Path):
    found = find_clustering_files(bigscape_out)
    assert sorted(found) == ["0.3", "0.5", "0.7"]
    assert found["0.3"].name == "mix_clustering_c0.30.tsv"


def test_find_clustering_files_ignores_network_and_stray_files(bigscape_out: Path):
    found = find_clustering_files(bigscape_out)
    for path in found.values():
        assert path.parent.name == "mix"
    assert "0.9" not in found


def test_find_clustering_files_skips_a_cutoff_with_no_mix_bin(bigscape_out: Path):
    mix = bigscape_out / "output_files" / "bigscape_2026-08-04_12-52-26_c0.5" / "mix"
    for tsv in mix.glob("*.tsv"):
        tsv.unlink()

    found = find_clustering_files(bigscape_out)
    assert sorted(found) == ["0.3", "0.7"]


def test_find_clustering_files_returns_empty_when_dir_is_missing(tmp_path: Path):
    assert find_clustering_files(tmp_path / "nope") == {}


def test_find_clustering_files_works_with_a_user_supplied_label(tmp_path: Path):
    """A user's own folder has a different run label. Globbing must still work."""
    root = tmp_path / "their_out" / "output_files"
    write_clustering_tsv(root / "myrun_2025-01-01_00-00-00_c0.4" / "mix" / "mix_clustering_c0.40.tsv")

    assert sorted(find_clustering_files(tmp_path / "their_out")) == ["0.4"]


# --- parsing ---------------------------------------------------------------


def test_parse_clustering_file_keys_on_sample_and_bgc_id(bigscape_out: Path):
    tsv = find_clustering_files(bigscape_out)["0.3"]
    families = parse_clustering_file(tsv, LABELS)

    assert families[("reference", "CONTIG_2.1")] == "FAM_00001"
    assert families[("assembly_10", "CONTIG_2.1")] == "FAM_00001"
    assert len(families) == len(ROWS)


def test_same_bgc_id_in_three_samples_does_not_overwrite(bigscape_out: Path):
    tsv = find_clustering_files(bigscape_out)["0.3"]
    families = parse_clustering_file(tsv, LABELS)

    shared = [k for k in families if k[1] == "CONTIG_6.1"]
    assert sorted(label for label, _ in shared) == sorted(LABELS)


def test_parse_clustering_file_returns_empty_on_missing_columns(tmp_path: Path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("Record\tRecord_Type\n" + "x\tregion\n", encoding="utf-8")

    assert parse_clustering_file(bad, LABELS) == {}


def test_unknown_labels_are_dropped_not_crashed(bigscape_out: Path):
    tsv = find_clustering_files(bigscape_out)["0.3"]
    families = parse_clustering_file(tsv, ["reference"])

    assert set(label for label, _ in families) == {"reference"}


# --- top level -------------------------------------------------------------


def test_parse_bigscape_returns_every_cutoff(bigscape_out: Path):
    families = parse_bigscape(bigscape_out, LABELS)

    assert sorted(families) == ["0.3", "0.5", "0.7"]
    assert len(families["0.3"]) == len(ROWS)
    assert len(set(families["0.3"].values())) == 2  # two families in the fixture


def test_parse_bigscape_never_raises_on_a_missing_folder(tmp_path: Path):
    assert parse_bigscape(tmp_path / "gone", LABELS) == {}


def test_parse_bigscape_never_raises_on_a_truncated_file(bigscape_out: Path):
    tsv = find_clustering_files(bigscape_out)["0.3"]
    tsv.write_text(HEADER + "\nreference_CONTIG_2.region001.gbk", encoding="utf-8")

    families = parse_bigscape(bigscape_out, LABELS)
    assert sorted(families) == ["0.5", "0.7"]


def test_parse_bigscape_returns_empty_when_output_files_is_absent(tmp_path: Path):
    empty = tmp_path / "bigscape_out"
    empty.mkdir()
    assert parse_bigscape(empty, LABELS) == {}


# --- cutoff selection ------------------------------------------------------


def test_select_cutoff_matches_across_padding(bigscape_out: Path):
    families = parse_bigscape(bigscape_out, LABELS)

    assert select_cutoff(families, 0.3) == families["0.3"]
    assert select_cutoff(families, "0.30") == families["0.3"]


def test_select_cutoff_lists_what_is_available(bigscape_out: Path):
    families = parse_bigscape(bigscape_out, LABELS)

    with pytest.raises(ValueError) as excinfo:
        select_cutoff(families, 0.9)

    message = str(excinfo.value)
    assert "0.9" in message
    assert "0.3, 0.5, 0.7" in message


def test_select_cutoff_returns_empty_when_there_is_no_data():
    assert select_cutoff({}, 0.3) == {}