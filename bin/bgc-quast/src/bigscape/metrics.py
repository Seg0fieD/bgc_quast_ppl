"""Gene cluster family (GCF) metrics built from BiG-SCAPE results.

Six metrics are registered here. They read `bgc.gcf_id`, which
`pipeline_helper.parse_input()` fills in for the cutoff the user selected.

Three of them (shared, unique, singleton) need to see every column at once, so the
calculator marks each BGC with two transient flags first. This copies the pattern
`CompareToolsMetricsCalculator` uses for `is_unique` (metrics_calculators.py:360-363).

Row labels for these metrics live in configs/report_config.yaml under `compare_samples`.
A metric with no YAML entry produces no row.
"""

from typing import Dict, Iterable, List, Optional, Set, Tuple

from src.genome_mining_result import Bgc, GenomeMiningResult
from src.reporting.metrics import METRIC_REGISTRY, metric
from src.reporting.metrics_calculators import BasicMetricsCalculator
from src.reporting.report_config import ReportConfig
from src.reporting.report_data import MetricValue

# Only antiSMASH results can carry GCF data. BiG-SCAPE cannot read GECCO or DeepBGC output.
ANTISMASH_TOOL = "antiSMASH"

# The metrics this module registers, in report row order.
GCF_METRIC_NAMES = [
    "gcf_count",
    "bgcs_in_gcf_count",
    "singleton_bgc_count",
    "shared_gcf_count",
    "unique_gcf_count",
    "mean_bgcs_per_gcf",
]


# ---------------------------------------------------------------------------
# Transient flags
# ---------------------------------------------------------------------------


def annotate_gcf_flags(results: List[GenomeMiningResult]) -> bool:
    """Mark every BGC with `gcf_is_singleton` and `gcf_is_shared`.

    Both facts are global, so they are computed once across all results and then
    stashed on the BGC objects, where the per-group metric functions can read them.

    - gcf_is_singleton: this BGC's family has exactly one member in the whole run.
    - gcf_is_shared:    this BGC's family appears in more than one column.

    Returns True if any BGC carried a gcf_id at all.
    """
    fam_sizes: Dict[str, int] = {}
    fam_columns: Dict[str, Set[str]] = {}

    for result in results:
        label = result.display_label or result.input_file_label
        for bgc in result.bgcs:
            if not bgc.gcf_id:
                continue
            fam_sizes[bgc.gcf_id] = fam_sizes.get(bgc.gcf_id, 0) + 1
            fam_columns.setdefault(bgc.gcf_id, set()).add(label)

    for result in results:
        for bgc in result.bgcs:
            if bgc.gcf_id:
                bgc.gcf_is_singleton = fam_sizes[bgc.gcf_id] == 1
                bgc.gcf_is_shared = len(fam_columns[bgc.gcf_id]) > 1
            else:
                bgc.gcf_is_singleton = False
                bgc.gcf_is_shared = False

    return bool(fam_sizes)


def _families_in_group(bgcs: Iterable[Bgc]) -> Dict[str, int]:
    """Count BGCs per family within one group. Ignores BGCs with no family."""
    counts: Dict[str, int] = {}
    for bgc in bgcs:
        if bgc.gcf_id:
            counts[bgc.gcf_id] = counts.get(bgc.gcf_id, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# Metric functions
#
# Each returns None when no BGC in the group has a family, which drops the row
# instead of printing a misleading 0 (metrics_calculators.py:124). This is the
# same pattern as mean_gene_per_bgc (metrics.py:121-132).
# ---------------------------------------------------------------------------


@metric("gcf_count")
def gcf_count(bgcs: Iterable[Bgc]) -> Optional[int]:
    """Number of distinct gene cluster families in this column."""
    counts = _families_in_group(bgcs)
    if not counts:
        return None
    return len(counts)


@metric("bgcs_in_gcf_count")
def bgcs_in_gcf_count(bgcs: Iterable[Bgc]) -> Optional[int]:
    """BGCs whose family has at least two members across the whole run."""
    items = [b for b in bgcs if b.gcf_id]
    if not items:
        return None
    return sum(1 for b in items if not getattr(b, "gcf_is_singleton", False))


@metric("singleton_bgc_count")
def singleton_bgc_count(bgcs: Iterable[Bgc]) -> Optional[int]:
    """BGCs that are alone in their family."""
    items = [b for b in bgcs if b.gcf_id]
    if not items:
        return None
    return sum(1 for b in items if getattr(b, "gcf_is_singleton", False))


@metric("shared_gcf_count")
def shared_gcf_count(bgcs: Iterable[Bgc]) -> Optional[int]:
    """Families in this column that also appear in at least one other column."""
    items = [b for b in bgcs if b.gcf_id]
    if not items:
        return None
    return len({b.gcf_id for b in items if getattr(b, "gcf_is_shared", False)})


@metric("unique_gcf_count")
def unique_gcf_count(bgcs: Iterable[Bgc]) -> Optional[int]:
    """Families found in this column only."""
    items = [b for b in bgcs if b.gcf_id]
    if not items:
        return None
    return len({b.gcf_id for b in items if not getattr(b, "gcf_is_shared", False)})


@metric("mean_bgcs_per_gcf")
def mean_bgcs_per_gcf(bgcs: Iterable[Bgc]) -> Optional[float]:
    """Average number of this column's BGCs per family it contains."""
    counts = _families_in_group(bgcs)
    if not counts:
        return None
    return sum(counts.values()) / len(counts)


# ---------------------------------------------------------------------------
# Calculator
# ---------------------------------------------------------------------------


class BigscapeMetricsCalculator(BasicMetricsCalculator):
    """Adds the GCF metrics to a compare-samples report.

    Only antiSMASH results are considered. Returns an empty list when no BGC carries
    a family, so a run without BiG-SCAPE adds nothing at all to the report.
    """

    def __init__(self, results: List[GenomeMiningResult], config: ReportConfig):
        self.results = results
        self.config = config

    def calculate_metrics(self) -> List[MetricValue]:
        antismash_results = [r for r in self.results if r.mining_tool == ANTISMASH_TOOL]
        if not antismash_results:
            return []

        if not annotate_gcf_flags(antismash_results):
            return []

        # The compare_samples config also lists four metrics upstream never implemented
        # (bgc_diversity, sample_similarity, core_bgcs, accessory_bgcs). They are not in
        # METRIC_REGISTRY, so asking for them would raise. Take only ours.
        configured = {m.name for m in self.config.metrics}
        metric_names = [n for n in GCF_METRIC_NAMES if n in configured]
        if not metric_names:
            return []

        all_metrics: List[MetricValue] = []

        for grouping_dims in self._generate_grouping_combinations(self.config):
            for result in antismash_results:
                try:
                    all_metrics.extend(
                        self._calculate_all_metrics_for_bgcs(
                            result.bgcs,
                            result.input_file,
                            result.mining_tool,
                            metric_names,
                            grouping_dims,
                        )
                    )
                except Exception as e:
                    print(f"Warning: Error calculating GCF metrics for {result.input_file}: {e}")

        return all_metrics


# ---------------------------------------------------------------------------
# Payload for the HTML report
#
# The table shows one cutoff. The dropdown needs every cutoff, and those numbers
# cannot come from bgc.gcf_id, which only holds the selected one. So these work
# straight off the parser's {cutoff: {(label, bgc_id): family}} map.
# ---------------------------------------------------------------------------


def _index_cutoff(
    families_at_cutoff: Dict[Tuple[str, str], str],
) -> Tuple[Dict[str, Dict[str, int]], Dict[str, Set[str]], Dict[str, int]]:
    """Return (per-column family counts, family -> columns, family -> total size)."""
    by_label: Dict[str, Dict[str, int]] = {}
    fam_columns: Dict[str, Set[str]] = {}
    fam_sizes: Dict[str, int] = {}

    for (label, _bgc_id), family in families_at_cutoff.items():
        column = by_label.setdefault(label, {})
        column[family] = column.get(family, 0) + 1
        fam_columns.setdefault(family, set()).add(label)
        fam_sizes[family] = fam_sizes.get(family, 0) + 1

    return by_label, fam_columns, fam_sizes


def summarise_cutoff(
    families_at_cutoff: Dict[Tuple[str, str], str],
    column_labels: List[str],
    display_names: Dict[str, str],
) -> Dict[str, object]:
    """Build the `rows` and `sets` for one cutoff.

    Args:
        families_at_cutoff: {(sample_label, bgc_id): family_id} for this cutoff.
        column_labels: the report's own column order. Never sort this.
        display_names: {metric_name: row label}, taken from the report config.

    Returns:
        {"rows": [[row_label, cell, cell, ...], ...], "sets": {label: [family, ...]}}
    """
    by_label, fam_columns, fam_sizes = _index_cutoff(families_at_cutoff)

    stats: Dict[str, Dict[str, float]] = {}
    for label in column_labels:
        counts = by_label.get(label, {})
        n_bgcs = sum(counts.values())
        n_gcfs = len(counts)
        singletons = sum(c for f, c in counts.items() if fam_sizes.get(f, 0) == 1)
        shared = sum(1 for f in counts if len(fam_columns.get(f, ())) > 1)

        stats[label] = {
            "gcf_count": n_gcfs,
            "bgcs_in_gcf_count": n_bgcs - singletons,
            "singleton_bgc_count": singletons,
            "shared_gcf_count": shared,
            "unique_gcf_count": n_gcfs - shared,
            "mean_bgcs_per_gcf": (n_bgcs / n_gcfs) if n_gcfs else 0.0,
        }

    rows = []
    for name in GCF_METRIC_NAMES:
        label = display_names.get(name, name)
        cells = []
        for column in column_labels:
            value = stats[column][name]
            cells.append(f"{value:.1f}" if isinstance(value, float) else str(value))
        rows.append([label] + cells)

    sets = {label: sorted(by_label.get(label, {})) for label in column_labels}

    return {"rows": rows, "sets": sets}


def build_bigscape_metadata(
    families: Dict[str, Dict[Tuple[str, str], str]],
    column_labels: List[str],
    default_cutoff: str,
    display_names: Dict[str, str],
    report_url: str = "../../bigscape/index.html",
) -> Optional[Dict[str, object]]:
    """Build the whole `metadata["bigscape"]` payload, one entry per cutoff.

    Returns None when there is nothing to show, so the key is simply absent and the
    HTML report renders exactly as it does today.
    """
    if not families or not column_labels:
        return None

    cutoffs = sorted(families, key=lambda c: float(c))
    rows: Dict[str, object] = {}
    sets: Dict[str, object] = {}

    for cutoff in cutoffs:
        summary = summarise_cutoff(families[cutoff], column_labels, display_names)
        rows[cutoff] = summary["rows"]
        sets[cutoff] = summary["sets"]

    return {
        "cutoffs": cutoffs,
        "default_cutoff": default_cutoff if default_cutoff in families else cutoffs[0],
        "columns": column_labels,
        "rows": rows,
        "sets": sets,
        "report_url": report_url,
    }