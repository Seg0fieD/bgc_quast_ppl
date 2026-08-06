"""Tests for src/bigscape/metrics.py.

The fixture is deliberately messier than the real test data: three columns, one
family in all three, one in two of three, and three singletons. The real 18-GBK run
gives six identical families of three, which cannot tell `shared` from `core` or
exercise the singleton path at all.
"""

from pathlib import Path

import pytest

from src.bigscape.metrics import (
    ANTISMASH_TOOL,
    GCF_METRIC_NAMES,
    BigscapeMetricsCalculator,
    annotate_gcf_flags,
    bgcs_in_gcf_count,
    build_bigscape_metadata,
    gcf_count,
    mean_bgcs_per_gcf,
    shared_gcf_count,
    singleton_bgc_count,
    summarise_cutoff,
    unique_gcf_count,
)
from src.genome_mining_result import Bgc, GenomeMiningResult
from src.reporting.report_config import ReportConfigManager

LABELS = ["reference", "assembly_10", "assembly_20"]

# family -> which columns hold it
#   F1  all three          shared, and core
#   F2  reference + a10    shared, NOT core
#   F3  reference only     unique, singleton
#   F4  assembly_10 only   unique, singleton
#   F5  assembly_20 only   unique, singleton
LAYOUT = {
    "reference": [("CONTIG_1.1", "F1"), ("CONTIG_2.1", "F2"), ("CONTIG_3.1", "F3")],
    "assembly_10": [("CONTIG_1.1", "F1"), ("CONTIG_4.1", "F2"), ("CONTIG_5.1", "F4")],
    "assembly_20": [("CONTIG_1.1", "F1"), ("CONTIG_6.1", "F5")],
}


def make_bgc(bgc_id: str, gcf_id=None, product: str = "NRPS") -> Bgc:
    return Bgc(
        bgc_id=bgc_id,
        sequence_id=bgc_id.split(".")[0],
        start=0,
        end=1000,
        product_types=[product],
        gcf_id=gcf_id,
    )


def make_result(label: str, bgcs, tool: str = ANTISMASH_TOOL) -> GenomeMiningResult:
    return GenomeMiningResult(
        input_file=Path(f"/tmp/{label}.json"),
        input_file_label=label,
        mining_tool=tool,
        bgcs=bgcs,
        display_label=label,
    )


@pytest.fixture
def results():
    """Three antiSMASH results with gcf_id already attached."""
    return [
        make_result(label, [make_bgc(bid, fam) for bid, fam in rows])
        for label, rows in ((lbl, LAYOUT[lbl]) for lbl in LABELS)
    ]


@pytest.fixture
def families():
    """The same layout as the parser's {(label, bgc_id): family} map."""
    flat = {}
    for label, rows in LAYOUT.items():
        for bgc_id, fam in rows:
            flat[(label, bgc_id)] = fam
    return flat


# --- flags -----------------------------------------------------------------


def test_annotate_marks_singletons_and_shared(results):
    assert annotate_gcf_flags(results) is True

    by_fam = {b.gcf_id: b for r in results for b in r.bgcs}
    assert by_fam["F1"].gcf_is_shared is True
    assert by_fam["F2"].gcf_is_shared is True   # 2 of 3 columns still counts as shared
    assert by_fam["F3"].gcf_is_shared is False
    assert by_fam["F1"].gcf_is_singleton is False
    assert by_fam["F3"].gcf_is_singleton is True


def test_annotate_returns_false_when_nothing_has_a_family():
    empty = [make_result("a", [make_bgc("CONTIG_1.1")])]
    assert annotate_gcf_flags(empty) is False
    assert empty[0].bgcs[0].gcf_is_shared is False


# --- metric functions ------------------------------------------------------


@pytest.mark.parametrize(
    "func", [gcf_count, bgcs_in_gcf_count, singleton_bgc_count,
             shared_gcf_count, unique_gcf_count, mean_bgcs_per_gcf]
)
def test_every_metric_returns_none_without_families(func):
    """None drops the row instead of printing a misleading 0."""
    assert func([make_bgc("CONTIG_1.1"), make_bgc("CONTIG_2.1")]) is None
    assert func([]) is None


def test_reference_column_values(results):
    annotate_gcf_flags(results)
    bgcs = results[0].bgcs  # F1, F2, F3

    assert gcf_count(bgcs) == 3
    assert singleton_bgc_count(bgcs) == 1      # F3
    assert bgcs_in_gcf_count(bgcs) == 2        # F1, F2
    assert shared_gcf_count(bgcs) == 2         # F1, F2
    assert unique_gcf_count(bgcs) == 1         # F3
    assert mean_bgcs_per_gcf(bgcs) == 1.0


def test_smallest_column_values(results):
    annotate_gcf_flags(results)
    bgcs = results[2].bgcs  # F1, F5

    assert gcf_count(bgcs) == 2
    assert singleton_bgc_count(bgcs) == 1      # F5
    assert shared_gcf_count(bgcs) == 1         # F1
    assert unique_gcf_count(bgcs) == 1         # F5


def test_shared_and_core_are_different(results):
    """F2 sits in 2 of 3 columns: shared, but not present everywhere."""
    annotate_gcf_flags(results)
    ref, a10, a20 = (r.bgcs for r in results)

    assert shared_gcf_count(ref) == 2
    assert "F2" in {b.gcf_id for b in ref}
    assert "F2" not in {b.gcf_id for b in a20}


def test_mean_counts_bgcs_not_families(results):
    annotate_gcf_flags(results)
    doubled = results[0].bgcs + [make_bgc("CONTIG_9.1", "F1")]
    assert mean_bgcs_per_gcf(doubled) == 4 / 3


# --- calculator ------------------------------------------------------------


@pytest.fixture
def compare_samples_config():
    """The real YAML block, so a renamed metric breaks this test."""
    return ReportConfigManager().get_config("compare_samples")


def test_calculator_emits_metrics(results, compare_samples_config):
    values = BigscapeMetricsCalculator(results, compare_samples_config).calculate_metrics()

    assert values
    names = {v.metric_name for v in values}
    assert names == set(GCF_METRIC_NAMES)


def test_calculator_returns_nothing_without_families(compare_samples_config):
    plain = [make_result(l, [make_bgc("CONTIG_1.1")]) for l in LABELS]
    assert BigscapeMetricsCalculator(plain, compare_samples_config).calculate_metrics() == []


def test_calculator_skips_non_antismash_tools(results, compare_samples_config):
    gecco = make_result("gecco_run", [make_bgc("CONTIG_1.1", "F1")], tool="GECCO")
    values = BigscapeMetricsCalculator(results + [gecco], compare_samples_config).calculate_metrics()

    assert all(v.mining_tool == ANTISMASH_TOOL for v in values)
    assert str(gecco.input_file) not in {str(v.file_path) for v in values}


def test_calculator_returns_nothing_when_only_other_tools(compare_samples_config):
    gecco = [make_result("g", [make_bgc("CONTIG_1.1", "F1")], tool="GECCO")]
    assert BigscapeMetricsCalculator(gecco, compare_samples_config).calculate_metrics() == []


def test_calculator_groups_by_product_type(results, compare_samples_config):
    results[0].bgcs[0].product_types = ["PKS"]
    values = BigscapeMetricsCalculator(results, compare_samples_config).calculate_metrics()

    groups = {tuple(sorted(v.grouping.items())) for v in values}
    assert () in groups                                    # the ungrouped Total row
    assert (("product_type", "PKS"),) in groups


# --- payload ---------------------------------------------------------------


DISPLAY = {
    "gcf_count": "# GCFs",
    "bgcs_in_gcf_count": "# BGCs in GCFs",
    "singleton_bgc_count": "# singleton BGCs",
    "shared_gcf_count": "# shared GCFs",
    "unique_gcf_count": "# unique GCFs",
    "mean_bgcs_per_gcf": "Mean BGCs per GCF",
}


def test_summarise_matches_the_metric_functions(results, families):
    """The table path and the dropdown path must not drift apart."""
    annotate_gcf_flags(results)
    summary = summarise_cutoff(families, LABELS, DISPLAY)

    rows = {row[0]: row[1:] for row in summary["rows"]}
    assert rows["# GCFs"] == ["3", "3", "2"]
    assert rows["# singleton BGCs"] == ["1", "1", "1"]
    assert rows["# BGCs in GCFs"] == ["2", "2", "1"]
    assert rows["# shared GCFs"] == ["2", "2", "1"]
    assert rows["# unique GCFs"] == ["1", "1", "1"]
    assert rows["Mean BGCs per GCF"] == ["1.0", "1.0", "1.0"]

    # same numbers the registered metrics produce for the first column
    assert int(rows["# GCFs"][0]) == gcf_count(results[0].bgcs)
    assert int(rows["# shared GCFs"][0]) == shared_gcf_count(results[0].bgcs)


def test_summarise_keeps_the_given_column_order(families):
    reversed_order = list(reversed(LABELS))
    summary = summarise_cutoff(families, reversed_order, DISPLAY)

    assert list(summary["sets"]) == reversed_order
    gcf_row = next(r for r in summary["rows"] if r[0] == "# GCFs")
    assert gcf_row[1:] == ["2", "3", "3"]


def test_summarise_sets_hold_family_ids(families):
    summary = summarise_cutoff(families, LABELS, DISPLAY)

    assert summary["sets"]["reference"] == ["F1", "F2", "F3"]
    assert summary["sets"]["assembly_20"] == ["F1", "F5"]


def test_summarise_handles_a_column_with_no_families(families):
    summary = summarise_cutoff(families, LABELS + ["empty_sample"], DISPLAY)

    assert summary["sets"]["empty_sample"] == []
    gcf_row = next(r for r in summary["rows"] if r[0] == "# GCFs")
    assert gcf_row[-1] == "0"


def test_metadata_covers_every_cutoff(families):
    payload = build_bigscape_metadata(
        {"0.3": families, "0.5": families, "0.7": families},
        LABELS, "0.5", DISPLAY,
    )

    assert payload["cutoffs"] == ["0.3", "0.5", "0.7"]
    assert payload["default_cutoff"] == "0.5"
    assert payload["columns"] == LABELS
    assert set(payload["rows"]) == {"0.3", "0.5", "0.7"}
    assert payload["report_url"].endswith("bigscape/index.html")


def test_metadata_sorts_cutoffs_numerically(families):
    payload = build_bigscape_metadata(
        {"0.7": families, "0.3": families, "0.15": families},
        LABELS, "0.3", DISPLAY,
    )
    assert payload["cutoffs"] == ["0.15", "0.3", "0.7"]


def test_metadata_falls_back_when_default_cutoff_is_absent(families):
    payload = build_bigscape_metadata({"0.5": families}, LABELS, "0.3", DISPLAY)
    assert payload["default_cutoff"] == "0.5"


def test_metadata_is_none_when_there_is_nothing_to_show(families):
    assert build_bigscape_metadata({}, LABELS, "0.3", DISPLAY) is None
    assert build_bigscape_metadata({"0.3": families}, [], "0.3", DISPLAY) is None