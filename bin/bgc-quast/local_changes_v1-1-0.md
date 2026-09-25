# Local changes to the vendored bgc-quast

This folder is upstream `gurevichlab/bgc-quast` at tag `bgc-quast-1.1.0`
(commit `45f50f2`, released 2026-09-10), with the changes below applied for
the `bgc_quast_ppl` pipeline. Every difference from that tag is listed here.

All changes serve one feature: an optional gene cluster family (GCF) layer
from BiG-SCAPE, shown in the compare-samples report. Nothing is active
unless `--bigscape-output-dir` is passed. Without that flag this copy
behaves as upstream v1.1.0.

## Re-measuring the change surface

The upstream clone's git directory is kept outside the working tree, so the
surface can be re-measured at any time:

    git --git-dir=v1.1.0/bgc-quast.git --work-tree=bin/bgc-quast status --short

The result must be the modified files in section 2 plus the new files in
section 1. Anything else is an unrecorded change.

## 1. New files

Three source files and two test files. None of them are upstream code.

- `bgc_quast/bigscape/__init__.py` - empty, marks the package.
- `bgc_quast/bigscape/parser.py` (240 lines) - reads BiG-SCAPE clustering
  files and maps `(sample label, BGC id)` to a family identifier, per
  cutoff. `bgc_id_from_record` rebuilds the bgc-quast BGC identifier from a
  staged region GenBank file name; for antiSMASH it emits
  `<contig>.reg.<number>`, which is the form v1.1.0's
  `get_antismash_bgc_id` produces for the region level.
- `bgc_quast/bigscape/metrics.py` (308 lines) - six GCF metrics, the
  calculator that produces them, and the metadata payload the HTML panel
  reads.
- `tests/test_bigscape_parser.py` (322 lines).
- `tests/test_bigscape_metrics.py` (289 lines).

## 2. Modified upstream files

### `bgc_quast/option_parser.py`

A new `BiG-SCAPE` argument group with two options: `--bigscape-output-dir`
(`-b`), the folder BiG-SCAPE was given as `-o`, and `--bigscape-cutoff`,
which cutoff the report table shows. Validation rejects a cutoff outside
`(0, 1]` and rejects a cutoff given without an output directory.

Effect on upstream behaviour: two added options. No existing option
changes.

### `bgc_quast/config.py`

`bigscape_cutoff` added to the `Config` dataclass, read from `config.yaml`
with a fallback of `0.3`, and overridden by `--bigscape-cutoff` when given.

Effect: none unless the flag is used.

### `bgc_quast/configs/config.yaml`

One key appended: `bigscape_cutoff: 0.3`.

Effect: none on its own; it supplies the default above.

### `bgc_quast/configs/report_config.yaml`

Two changes.

Six GCF metric entries added under `compare_samples`: `gcf_count`,
`bgcs_in_gcf_count`, `singleton_bgc_count`, `shared_gcf_count`,
`unique_gcf_count` and `mean_bgcs_per_gcf`. A comment marks the metrics
upstream declares but has not implemented, so the two groups are not
confused.

The `sample_group` grouping dimension is commented out rather than deleted.
Upstream ships it with placeholder group names (`group_A`, `group_B`,
`group_C`) that no pipeline input supplies. Commenting has the same effect
as deleting, keeps the diff small, and leaves the local change visible to a
reader.

Effect: the GCF rows appear only when BiG-SCAPE data is present.
`sample_group` is inert.

### `bgc_quast/genome_mining_result.py`

One field on the `Bgc` dataclass: `gcf_id`, defaulting to `None`, holding
the family a BGC belongs to at the selected cutoff.

Effect: none when BiG-SCAPE is not run; the field stays `None`.

### `bgc_quast/pipeline_helper.py`

Imports `parse_bigscape` and `select_cutoff`. Two attributes added:
`bigscape_families` and `bigscape_report_url`.

A block runs after display labels are assigned and deduplicated, because
the join key is `(display label, BGC id)` and both parts must be final. It
parses the BiG-SCAPE folder, selects the cutoff, writes `gcf_id` onto every
matching BGC, and logs how many BGCs were assigned to a family. It then
builds a relative link to BiG-SCAPE's own HTML report.

The block refuses to run when the results hold more than one mining tool. A
BiG-SCAPE folder describes one tool, so in compare-tools mode the families
would be attached to BGCs they did not come from. A warning is logged and
the GCF rows are skipped.

Two keyword arguments are passed on to `build_report`.

Effect: none without `--bigscape-output-dir`.

Note: the v1.0.0 copy had a wildcard import of the metrics module here, to
force the `@metric` decorators to register. It was deliberately not
re-applied. `report_builder.py` imports that module anyway, so the
decorators still run, and a wildcard import is the first thing a reviewer
objects to.

### `bgc_quast/reporting/report_builder.py`

Imports the calculator, the metadata builder and `normalize_cutoff`. Two
keyword arguments added to `build_report`. The compare-samples branch,
which upstream leaves empty, now builds the GCF metrics and the panel
metadata when family data is present.

Column labels for the panel follow the order of `results`, which is the
report's own column order. Sorting them here would break the match with the
table.

Effect: the compare-samples branch is a no-op without family data, as
upstream.

### `bgc_quast/html_report/build_report.js`

293 lines added; 5 lines replaced.

Five new functions: `gcfVennRegions`, `drawVennGcf`, `buildVennPicker`,
`buildGcfSummaryTable` and `initGcfPanel`.

The five replaced lines are the panel dispatch. Upstream shows its Python
plots panel in compare-tools only. It now also opens in compare-samples
when family data is present, showing the GCF panel, with the tab relabelled
`GCF overlap`. The compare-tools path is unchanged.

Effect: no change to any existing report without family data.

### `bgc_quast/html_report/report.css`

80 lines appended at the end, under a `BiG-SCAPE GCF PANEL` comment.
Nothing removed and no upstream rule edited. 639 lines upstream, 719 here.

Effect: new rules only, all scoped to the GCF panel.

### `pyproject.toml`

`bgc_quast.bigscape` added to the package list.

Effect: none on the pipeline, which runs the tool in place from this
folder. Without it an installed copy would be missing the package.

### `tests/test_pipeline_helper.py`

An upstream test. The asserted `build_report` argument list gained
`bigscape_families` and `bigscape_report_url`, because the signature
changed.

Effect: the test now matches the local signature and would fail against
clean upstream.

## 3. Changes from the v1.0.0 copy that were deliberately dropped

Upstream fixed both of these itself in v1.1.0, so the local patches were
removed rather than re-applied.

- A widened regex for fuzzy BGC coordinates. Upstream's v1.1.0 pattern is
  wider than the local one was.
- A one-line DOI correction in `README.md`. Upstream corrected the same
  DOI, so the vendored README now differs from upstream in no way.

## 4. Incidental differences, to be reverted

These are not intentional changes and carry no behaviour.

- `reporting/report_builder.py`: upstream's `# TODO:` marker reads
  `# TODO_`, and a whitespace-only line was added at the end of the file.
- `config.py` and `option_parser.py`: two whitespace-only lines.
- `dev/prism_validation/*/bgc-quast.log`: five upstream example log files
  are absent from this copy.

## 5. What to re-check after any future upstream pull

- `configs/report_config.yaml`: the `sample_group` grouping dimension is
  restored by a pull and must be commented out again.
- `tests/test_pipeline_helper.py`: a pull restores the upstream argument
  list and the two keywords must be added back.
- The antiSMASH BGC identifier form. v1.1.0 changed it from
  `<contig>.<number>` to `<contig>.reg.<number>` and nothing errored - the
  family rows simply vanished from one report. After any upgrade, compare
  the per-tool assignment count in `bgc-quast.log` against that tool's
  total BGC count in the report. They must be equal.
