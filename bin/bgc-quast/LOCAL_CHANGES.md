# Local changes to the vendored bgc-quast

This folder is a vendored copy of [gurevichlab/bgc-quast](https://github.com/gurevichlab/bgc-quast)
v1.0.0, used by the `bgc_quast_ppl` Nextflow pipeline. The upstream repo is **not** modified.
Every change made here is listed below so it can be reviewed and, later, offered upstream.

**Purpose of the changes:** add BiG-SCAPE gene cluster family (GCF) results to the
compare-samples report. bgc-quast reads a finished BiG-SCAPE output folder, exactly the way
it already reads a QUAST output folder. It does not run BiG-SCAPE.

**Guarantee:** without `--bigscape-output-dir`, every report is byte-identical to upstream.
Verified by running compare-samples with and without the flag and diffing `report.txt`,
`report.tsv` and `report.html`.

Last updated: 2026-08-12. The antiSMASH round (Phases 0–6) is complete and verified end to end
on pipeline output. The GECCO/DeepBGC round is in progress: the parser, the tool filters and the
report-link derivation are done and verified; the pipeline side is partly built.

---

## New files

### `src/bigscape/__init__.py`
**What:** Empty file that makes `src.bigscape` a package.
**Why:** All new Python lives in one folder, so the change surface is one directory instead
of edits scattered through `src/`.
**Upstream impact:** None.

### `src/bigscape/parser.py`
**What:** Reads a BiG-SCAPE `cluster` output folder and returns
`{cutoff: {(sample_label, bgc_id): family_id}}`. Also `normalize_cutoff` and `select_cutoff`.
**Why:** BiG-SCAPE names each BGC by its input file name. The pipeline stages every antiSMASH
region GBK as `<sample_label>_<original_name>`, and this reverses that to rebuild bgc-quast's
own `sequence_id.region_number` id. The key includes the sample because the same `bgc_id`
(e.g. `CONTIG_6.1`) occurs in more than one sample.
Output folder names carry a run label and a timestamp, so paths are found by globbing
`output_files/*_c*/mix/*_clustering_c*.tsv`, never built from a fixed string.
**Upstream impact:** None. Nothing calls it unless the new flag is given. A missing or broken
folder returns `{}` and never raises into the run.
**Update 2026-08-11 — two naming schemes.** `bgc_id_from_record` now auto-detects which marker
a staged file name carries: `.region` first (antiSMASH), then `_cluster_` (GECCO, and DeepBGC
once its multi-record GBK is split). `.region` is tried first, so the antiSMASH path is
byte-for-byte the same code it was. The two branches differ on purpose: antiSMASH pads the
region number in the file name while bgc-quast reads the unpadded qualifier, so the padding is
stripped; GECCO and DeepBGC pad nothing and bgc-quast keeps the raw text, so the number is kept
verbatim. `rpartition` takes the last marker, which handles both a sequence id containing a dot
(`NC_003888.3_cluster_2`) and one containing `_cluster_`.
**Upstream impact:** Still none. A name with neither marker returns `None` and the row is
dropped, exactly as before.

### `src/bigscape/metrics.py`
**What:** Six registered metrics (`gcf_count`, `bgcs_in_gcf_count`, `singleton_bgc_count`,
`shared_gcf_count`, `unique_gcf_count`, `mean_bgcs_per_gcf`), a `BigscapeMetricsCalculator`,
and the builders for the HTML payload.
**Why:** Three of the metrics need to see every column at once, so the calculator marks each
BGC with two transient flags first. This copies the `is_unique` pattern in
`CompareToolsMetricsCalculator`. Each metric returns `None` when no BGC in the group has a
family, so the row is dropped instead of printing a misleading `0` — same as
`mean_gene_per_bgc`.
**Upstream impact:** None. Importing the module registers the metrics, but a metric only
produces a row when BGCs actually carry a family.
**Update 2026-08-12 — the tool filter is gone.** `BigscapeMetricsCalculator` used to keep only
`mining_tool == "antiSMASH"` results. It now uses every result it is given. BiG-SCAPE reads GECCO
and DeepBGC GBKs too, and in `compare-samples` each bgc-quast run already holds exactly one tool,
so the filter was a no-op there. A run holding more than one tool is stopped earlier, in
`pipeline_helper.py`. The `ANTISMASH_TOOL` constant stays — the tests use it — but it now records
the exact spelling `genome_mining_parser.py` writes, not a rule.
**Upstream impact:** Still none without `--bigscape-output-dir`.

### `tests/test_bigscape_parser.py`
**What:** 31 tests for the parser. Builds a fake output tree, so no BiG-SCAPE, Pfam or Docker
needed.
**Why:** Covers path discovery, cutoff normalisation, longest-prefix label matching, and the
"never raise" contract. Nine tests added 2026-08-11 for the `_cluster_` marker: GECCO and split
DeepBGC names, a sequence id carrying a dot, the number kept verbatim, the last marker winning,
broken names returning `None`, `.region` taking priority, and a sample label that itself
contains `_cluster_`.
**Upstream impact:** None.

### `tests/test_bigscape_metrics.py`
**What:** 25 tests for the metrics, the calculator and the HTML payload.
**Why:** The fixture is deliberately messier than the real test data — one family in three
columns, one in two of three, three singletons — because the real data gives six identical
families and cannot tell "shared" from "present everywhere".
One test checks that the table numbers and the dropdown numbers agree; they come from two
different code paths and must not drift.
**Upstream impact:** None.
**Update 2026-08-12 — two tests replaced, count unchanged at 25.** Both asserted the tool filter
that was removed from `metrics.py`, so they described behaviour that no longer exists. They were
replaced rather than repaired: `test_calculator_skips_non_antismash_tools` became
`test_calculator_uses_every_tool`, and `test_calculator_returns_nothing_when_only_other_tools`
became `test_calculator_works_for_a_non_antismash_tool`. The three payload tests that assert on
`report_url` were deliberately left untouched — see the `report_builder.py` note below.

---

## Edited files

### `src/genome_mining_result.py`
**What:** One optional field on `Bgc`: `gcf_id: Optional[str] = None`, appended last.
**Why:** Somewhere to hold the family for the selected cutoff. Appended because `bgc_id` and
`sequence_id` have no defaults, so every existing positional construction keeps working.
**Upstream impact:** None. Defaults to `None`.

### `src/genome_mining_parser.py`

**What:** One regex in `parse_antismash_json`. The region location pattern
`\[(\d+):(\d+)\](?:\((\+|-)\))?` became `\[<?(\d+):>?(\d+)\](?:\((\+|-)\))?`.

**Why:** An upstream bug, not ours. A BGC that runs off the end of a contig gets a fuzzy
GenBank bound — antiSMASH writes `[0:>24677](+)`. The old pattern does not match it, the
`ValueError` becomes an `InvalidInputException`, and `parse_input_mining_result_files` then
falls through to the next parser. `parse_deepbgc_json` accepts any JSON and returns an empty
list, so **the whole file** is silently reported as a DeepBGC run with zero BGCs. Found on a
draft isolate assembly with 261 contigs: one fuzzy region cost all 63 of that sample's regions.
In compare-samples the mismatched tool aborts the run; in compare-tools it would have produced
a wrong report with no warning.

The coordinate itself is correct — `<` and `>` mark uncertainty about the boundary, not about
the number — so ignoring the marker is the right reading.

**Upstream impact:** Files that parsed before parse identically. Files that were previously
misrouted to DeepBGC now parse as antiSMASH, which is what they are.

**Verified:** pytest failure list unchanged at 38 lines.

### `configs/config.yaml`
**What:** One new setting, `bigscape_cutoff: 0.3`.
**Why:** Which cutoff the report table opens on. Kept beside the other tuning defaults.
**Upstream impact:** None.

### `src/config.py`
**What:** `bigscape_cutoff: float = 0.3` on `Config`, read in `load_config` with
`cfg.get(...)`, plus a CLI override next to the existing four.
**Why:** Mirrors `compare_tools_overlap_threshold` exactly.
The field has a default (unlike its neighbours) so any `Config(...)` call that omits it still
works. `cfg.get` rather than `cfg[...]` so an older `config.yaml` still loads.
**Upstream impact:** None.

### `src/option_parser.py`
**What:** A "BiG-SCAPE" argument group with `--bigscape-output-dir` / `-b` and
`--bigscape-cutoff`. One check in `validate_arguments`: a cutoff without a directory is an
error, and the cutoff must be in (0, 1].
**Why:** `-b` was the only sensible short flag still free (`-o -g -t -r -q -R -h` are taken).
`--bigscape-cutoff` is long-only, matching `--merge-distance`, `--min-bgc-length` and
`--edge-distance`. Its argparse default is `None`, not `0.3`, so the YAML default wins unless
the user actually types the flag.
**Upstream impact:** None. Two new optional flags.
**Update 2026-08-12 — a duplicate line removed.** `add_mode_specific_arguments` created the
`Compare-tools` argument group twice in a row. The first call was dead — argparse suppresses a
section with no actions, so `--help` never showed it — but it was left over from the BiG-SCAPE
group edit. Deleted.
**Upstream impact:** None. `--help` output is unchanged.

### `configs/report_config.yaml`
**What:** Two changes.
1. Added six metric entries to the `compare_samples` block, with display names and one
   `precision: 1`.
2. **Deleted the `sample_group` grouping dimension** from that block.
**Why (1):** Row labels come from this file, not from Python. A registered metric with no
entry here produces no row. The four unimplemented upstream metrics (`bgc_diversity`,
`sample_similarity`, `core_bgcs`, `accessory_bgcs`) were left exactly as they were.
**Why (2):** See "Upstream bugs found" below.
**Upstream impact:** The six new rows appear only when BiG-SCAPE data is present. The
`sample_group` deletion is a bug fix; verified it changes no output (see below).
**Note 2026-08-12:** the six rows needed **no change** to work for GECCO and DeepBGC. They were
already worded generically ("in this sample"), which was confirmed by hand runs against both
tools' BiG-SCAPE output.

### `src/pipeline_helper.py`
**What:** Two imports, one attribute (`self.bigscape_families`), an attach block at the end
of `parse_input()`, and one extra argument in the `build_report(...)` call.
**Why:** The attach block sets `bgc.gcf_id` for antiSMASH results. It sits **after**
`assign_and_deduplicate_display_labels`, because the join key is
`(display_label, bgc_id)` and `display_label` does not exist before that call.
`parse_bigscape` swallows errors and returns `{}`; `select_cutoff` raises. A broken folder
should not kill a run, but asking for a cutoff that was never computed should say so rather
than quietly showing a different one.
**Upstream impact:** None. The whole block is inside `if self.args.bigscape_output_dir`.
**Update 2026-08-12 — mixed-tool guard, and the report link.** Two changes. First, the attach
block no longer skips non-antiSMASH results, and a guard replaces that skip: if a BiG-SCAPE
folder is given and the run holds more than one `mining_tool`, it warns and skips the GCF layer
entirely. One folder describes one tool, so matching its families against three tools' BGCs would
produce numbers that look plausible and are wrong. Second, `self.bigscape_report_url` is derived
here from the folder passed to `--bigscape-output-dir` and handed to `build_report()`. The
`../../` in that link describes the pipeline's published tree, not anything about bgc-quast, so
it is built here rather than inside `src/bigscape/metrics.py`, which stays free of path handling.
**Upstream impact:** None. Both live inside `if self.args.bigscape_output_dir`.
**Update 2026-08-12 (second) — the derivation must not resolve the path.** The first version
called `Path(...).resolve()` on the folder it was given. Nextflow stages an input directory as a
**symlink**, so resolving followed the link into the work directory, the parent stopped being
`bigscape`, and the link came out one level short. Dropping `.resolve()` fixed it, and the path
is already absolute so nothing else was needed. Verified end to end: the report now carries
`"report_url": "../../bigscape/antismash/index.html"`.
**Upstream impact:** None. Upstream has no staged-symlink layout, so the behaviour is the same
either way there.

### `src/reporting/report_builder.py`
**What:** Three imports, a new defaulted keyword argument `bigscape_families=None` on
`build_report()`, and the previously empty `COMPARE_SAMPLES` branch filled in.
**Why:** `build_report()` never sees the CLI args, and `bgc.gcf_id` holds only the selected
cutoff, so the full multi-cutoff map has to arrive as an argument for the HTML dropdown to
have anything to switch between. The branch mirrors the `COMPARE_TOOLS` branch directly
above it, and ships its payload through the existing `metadata` dict — the same channel
compare-tools uses for `pairwise_by_run` — so no new HTML placeholder was needed.
The whole branch is wrapped in `if bigscape_families:`, so with BiG-SCAPE off it does exactly
what the old bare `...` did.
**Upstream impact:** None. The new argument is optional and defaults to `None`.
**Update 2026-08-12.** The `column_labels` list no longer filters on antiSMASH, and a second
defaulted keyword argument `bigscape_report_url=None` was added beside `bigscape_families`. It is
forwarded to `build_bigscape_metadata` **only when set**, so that function keeps its own default
and its existing tests are untouched.
**Upstream impact:** None. Both arguments are optional and default to `None`.

### `src/html_report/build_report.js`
**What:** Added `drawVennGcf` (3-set Venn), `gcfVennRegions` (computes the seven
region counts from family membership sets), `buildGcfSummaryTable`, and
`initGcfPanel`. Changed the tab gate in `DOMContentLoaded` so `'pyplots'` is added
to `METRIC_TABS_BY_MODE.compare_samples` only when `reportMetadata.bigscape`
exists, and added a `compare_samples` branch beside the existing
`mode === 'compare_tools'` one that calls `initGcfPanel`.
**Why:** The existing `drawVenn` is 2-set, tool-labelled and driven by precomputed
counts in `metadata.pairwise_by_run`. A GCF Venn needs three circles and real set
membership, so it is new code rather than a change to `drawVenn`. Writing the
function alone was not enough: `METRIC_TABS_BY_MODE.compare_samples` was `['bgcs']`
and `initVennPanel` ran only under `compare_tools`, so the panel never rendered.
**Upstream impact:** None without BiG-SCAPE. `drawVenn` and `initVennPanel` are
untouched, so the compare-tools Venn is unchanged. The tab array is computed at
runtime, so a compare-samples run with no `reportMetadata.bigscape` shows only the
"All BGCs" tab, exactly as before. The `pyplots` button and `pythonPlotsPanel` div
already existed in `report_template.html`; the template was not modified.
**Note 2026-08-12:** this file needed **no change** for GECCO and DeepBGC. `drawVennGcf` was
already tool-agnostic, and compare-samples has three columns whatever the tool. Confirmed by
hand runs: both new reports render a correct 3-circle Venn, and for DeepBGC the seven region
counts sum to 40, which is exactly the number of families counted independently from the
BiG-SCAPE clustering file.

### `src/html_report/report.css`
**What:** Appended `.gcf-*` classes at the end of the file for the GCF panel: the
cutoff dropdown, the summary table, the note line, and the link.
**Why:** The new panel needed styling and the existing classes did not fit.
**Upstream impact:** None. Additions only, all namespaced under `.gcf-`. No
existing selector was changed, so nothing already in the report can be affected.

---

## Not touched

- `src/reporting/report_formatter.py` — the `metadata` dict already reaches the browser.
- `src/html_report/report_template.html` — no new placeholder needed; the `pyplots` tab
  button and `pythonPlotsPanel` div already exist.
- `src/reporting/metrics.py` and `src/reporting/metrics_calculators.py` — the GCF versions
  live in `src/bigscape/` instead.
- `dev/html_report_experiments/` — a stale duplicate of the HTML assets. Left alone.
- The four unimplemented `compare_samples` metrics in `configs/report_config.yaml`.

---

## Pipeline-side work — not in this folder

The Nextflow side of both rounds touched **pipeline files only**, none of which are vendored
bgc-quast code, so nothing was added to the per-file log above:

`modules/local/bigscape.nf`, `modules/local/bigscape_download_db.nf`,
`modules/local/deepbgc_split_gbk.nf`, `modules/local/bgcquast.nf`,
`subworkflows/local/bgc_prediction.nf`, `subworkflows/local/bgcquast.nf`,
`workflows/bgc_quast_ppl.nf`,
`subworkflows/local/utils_nfcore_bgc_quast_ppl_pipeline/main.nf`, `conf/modules.config`,
`nextflow.config`, `nextflow_schema.json`, `README.md`.

Those files are ours and are tracked by git normally. The pytest failure list was unchanged
throughout.

**One exception worth naming.** The GECCO/DeepBGC round's Phase 4 restructured the published
BiG-SCAPE folder into `bgc_quast/bigscape/<tool>/`, and that **did** require a vendored change —
the `.resolve()` fix logged under `src/pipeline_helper.py` above. A pipeline layout change
reaching into vendored code is exactly the kind of coupling worth keeping small, which is why the
link is derived in one place and passed down as a plain string.

---

## The GECCO and DeepBGC round — where it stands

The goal is to let BiG-SCAPE run on GECCO and DeepBGC predictions too, not only antiSMASH, so
all three compare-samples reports grow GCF rows. Neither tool's TSV is readable by BiG-SCAPE,
but both also emit GBKs, which are.

**Done in this folder:**

1. The parser learned a second naming scheme, `_cluster_` — see `src/bigscape/parser.py`.
2. The places that filtered on antiSMASH were opened and a mixed-tool guard put in their place —
   see `src/bigscape/metrics.py`, `src/pipeline_helper.py`, `src/reporting/report_builder.py`.
3. The BiG-SCAPE report link is now derived from the folder bgc-quast was given, rather than
   hardcoded — see `src/pipeline_helper.py`.

**Verified by hand, no pipeline**, by pointing bgc-quast at BiG-SCAPE folders produced outside
Nextflow. Every BGC joined in both runs — 18 of 18 for GECCO, 115 of 115 for DeepBGC — and the
table, the Venn and BiG-SCAPE's own clustering file agree on every number.

**Done outside this folder:** a module that splits DeepBGC's multi-record GBK and adds the
`/note="Cluster number: N"` qualifier BiG-SCAPE requires, and the per-tool `BIGSCAPE` aliasing
with nested output. Neither touches vendored code.

**Not done:** the channel wiring that actually feeds GECCO and DeepBGC GBKs to BiG-SCAPE, and the
params to switch tools on. Neither will appear in this file.

---

## Upstream bugs found

Two problems in the vendored copy that are not ours, reported here so they can go upstream.

### 1. `report_config.yaml` declared a grouping dimension that does not exist

The `compare_samples` block listed a `sample_group` grouping dimension. Only two grouping
keys are registered in `metrics.py`: `completeness` and `product_type`. Any calculator built
on that config hits `ValueError: Unknown grouping key`, which is swallowed by the
`try/except Exception` in `metrics_calculators.py`. The result is a
`Warning: Error calculating metrics for <file>` line and silently dropped rows.

It never fired before because the `COMPARE_SAMPLES` branch was empty, but
`report_writer.py` already loads that block on every compare-samples run.

**Fix applied:** deleted the three lines.
**Verified safe:** `sample_group` appears nowhere else in the repository. compare-samples was
run before and after the deletion; `report.txt`, `report.tsv` and `report.html` were all
byte-identical, and the test-suite result was unchanged.

### 2. `tests/` is out of sync with `src/`

On a clean checkout the vendored test suite does not pass — roughly 20 failures and 18 errors
out of 114. The tests were written against an older version of the source:

- they use `RunningMode.UNKNOWN`; the enum has only `COMPARE_TO_REFERENCE`, `COMPARE_TOOLS`
  and `COMPARE_SAMPLES`
- `determine_running_mode()` is called with an older signature
- they expect `"Unknown"` where the code now produces `"Unknown completeness"` and
  `"Unknown product"`

**Not fixed here.** Repairing it means rewriting most of the suite or changing behaviour, and
either would swamp this diff. It is a separate piece of work.

**What we do instead:** record the failure list once as a baseline, and require that it does
not change after each edit. Every change above was checked that way, and the new
`tests/test_bigscape_*.py` files pass in full.

Run from `bin/bgc-quast`; `../../tmp_test/` is the repo's scratch folder.

```bash
python -m pytest tests/ -q 2>&1 | grep -E "^(FAILED|ERROR)" | sort > ../../tmp_test/pytest_baseline.txt
```

The baseline is **38 `FAILED`/`ERROR` lines**. The passing count grows as tests are added — it
went from `20 failed / 123 passed / 18 errors` to `20 failed / 132 passed / 18 errors` when the
nine parser tests landed. **The guard is that the failure list is identical, not that the suite
is green.**

---

## One thing to know before testing a change here

**Nothing in this folder is part of a Nextflow task hash.** The `BGCQUAST` process only calls
`python3 bgc-quast.py`, so editing any file here leaves the hash unchanged and `-resume` will
serve the cached report. After any change in this folder, re-run **without** `-resume`, or the
result is stale.a a 