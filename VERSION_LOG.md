# bgc_quast_ppl: Version log

All notable changes to this pipeline are listed here, newest first.

---

## v2.1.0 - 2026-09-30

### Changed behaviour

- `--run_bigscape` with compare-tools or compare-to-reference now stops the
  run at start-up with an error. Before, the run was accepted. BiG-SCAPE
  works only in compare-samples mode.
- `--bgc_quast_mode auto`, or any other unsupported mode, now stops the run
  at start-up, before the parameter summary.
- `--bgc_bigscape_classify` now accepts only `none`. Any other value stops
  the run at start-up. The report reads only that one bin.
- Only Nextflow 25.10.x runs (25.10.4 or newer). Other versions stop at
  start-up.

### Added

- A start-up header with the pipeline name and version.
- Every error, warning and note now sits in a white frame. The colour shows
  the kind: pink for start-up errors, red for errors during the run, yellow
  for warnings, orange for notes.
- Start-up warnings are shown together in one frame, under the list of
  steps, with `&` between them.
- A note when compare-samples mode runs with a single sample.
- A new check stops the run when two BiG-SCAPE input files would get the
  same name, for example samples `a` and `a_b`. The error names the samples
  to rename.

### Fixed

- A blank sample name in the samplesheet is now filled from the file name.
  Before, the input check rejected it.
- The "no BGCs predicted" notes no longer appear when a tool crashed. They
  are shown only when the run succeeds.
- Warnings about samples with no long contigs or an empty annotation no
  longer vanish from the screen. They are shown at the end of the run.
- The error for a reference with no BGCs now says that bgc-quast did not
  run. It wrongly said QUAST did not run either.
- Options in messages print as plain text, without underline.
- One dashed line under the parameter summary instead of two.
- Start-up errors now respect `--monochrome_logs`.

### Other

- Code comments and descriptions rewritten throughout. Results are
  unchanged.

---

## v2.0.1 - 2026-09-27

### Fixed

- HTML reports were blank in v2.0.0. A syntax error in the report script
  stopped the page from loading. The tables, charts and notes now show
  again.

### Changed

- The BiG-SCAPE steps now run as their own stage, named
  `BIGSCAPE_ANALYSIS` in the progress display. Results are unchanged.

---

## v2.0.0 - 2026-09-27

This version moves to bgc-quast 1.1.0 and removes two parameters, so some
older commands no longer work. See [Parameters](#parameters) below.

### Added

- Reports now show the **mean BGC length in genes**. The pipeline gives
  bgc-quast the Pyrodigal gene annotation instead of the plain genome.
- Compare-tools mode writes four new files per sample: a table of every BGC
  from every tool, a table and an HTML page of overlapping BGCs, and a
  GenBank file with all BGCs marked on the genome.
- A short note under the report table explains how the Venn diagram numbers
  differ from the table numbers. It appears only when there are more than
  three samples.
- A new README with a table of contents, every parameter, notes and a
  troubleshooting guide.
- This version log.

### Changed

- bgc-quast updated from version 1.0.0 to **1.1.0**. The BiG-SCAPE changes
  were added again on top of the new version.
- **Total BGC length** is replaced by **total BGC span**, which bgc-quast
  1.1.0 measures differently. Numbers in this row are not comparable with
  reports from earlier versions.
- Each genome is now annotated **once**. Before, every genome was annotated
  twice and one copy was thrown away. `--save_annotations` now saves one
  folder instead of two.

### Fixed

- **BiG-SCAPE family numbers are now the same on every run.** Before, the
  same input could give slightly different gene cluster families. BiG-SCAPE
  loaded its input files in a random order, which changed some distances.
  The pipeline now fixes that order and seeds BiG-SCAPE's random numbers.
- A run that produced no comparison could still report success if the
  output folder held results from an earlier run. The pipeline now checks
  what the current run produced.

### Removed

- `--bgc_quast_output_bgcs`. bgc-quast 1.1.0 no longer has this option. The
  BGC GenBank file is now always written.
- The `prodigal` choice for `--annotation_tool`, and its four
  `--annotation_prodigal_*` options. Prodigal's output cannot be read by
  antiSMASH, DeepBGC or GECCO, so it could never produce a report.

### Parameters

| Old parameter | New parameter |
|---|---|
| `--bgc_quast_output_bgcs` | removed |
| `--annotation_tool prodigal` | removed, use `pyrodigal` |
| `--annotation_prodigal_singlemode` | removed |
| `--annotation_prodigal_closed` | removed |
| `--annotation_prodigal_transtable` | removed |
| `--annotation_prodigal_forcenonsd` | removed |

---

## v1.1.0 - 2026-09-23

This version adds gene cluster families with BiG-SCAPE. All older commands
still work.

### Added

- Optional **BiG-SCAPE** step, turned on with `--run_bigscape`. It groups
  similar BGCs into gene cluster families and adds family rows and a Venn
  diagram to the reports. It runs once for each prediction tool
  (antiSMASH, DeepBGC and GECCO), in compare-samples mode only.
- Automatic download of Pfam 38.2 when no Pfam file is given. The release is
  fixed on purpose, because a different Pfam release can change the
  families.
- `--bgc_bigscape_dir` to reuse BiG-SCAPE results from an earlier run.
- Several family cut-offs per run, with a cut-off menu in the HTML report.
- More checks before the run starts, with clear messages that list every
  problem at once.

### Changed

- A sample in which a tool finds no BGCs is left out of that tool's report,
  and the run says so. antiSMASH now behaves like GECCO and DeepBGC here.
- With an empty reference in compare-to-reference mode, the pipeline stops
  with a clear message instead of making an empty report.
- BiG-SCAPE now runs only in compare-samples mode. In the other modes
  `--run_bigscape` is skipped with a warning, and the run continues.
- The pipeline diagram was updated.

### Fixed

- antiSMASH BGCs cut off by a contig edge are now read correctly.
- GECCO now reads the gene annotation made by the pipeline.

### Removed

- An unused step that unzipped Pyrodigal's gene nucleotide file.

### Parameters

| Old parameter | New parameter |
|---|---|
| - | `--run_bigscape` |
| - | `--bgc_bigscape_pfam` |
| - | `--bgc_bigscape_pfam_url` |
| - | `--bgc_bigscape_dir` |
| - | `--bgc_bigscape_cutoffs` |
| - | `--bgc_bigscape_cutoff` |
| - | `--bgc_bigscape_classify` |

---

## v1.0.0 - 2026-08-05

First version.

### Added

- Input preparation: unzip, and removal of contigs shorter than
  `--bgc_mincontiglength` (3,000 bp by default).
- Gene annotation with Pyrodigal.
- BGC prediction with antiSMASH, DeepBGC and GECCO. Each tool can be
  skipped. antiSMASH runs in minimal mode by default; `--bgc_antismash_full`
  runs the full analysis.
- Comparison with bgc-quast 1.0.0, kept inside the pipeline, in three modes:
  compare-samples, compare-tools and compare-to-reference.
- QUAST alignment, run automatically in compare-to-reference mode. The
  reference is marked in the samplesheet with a `type` column.
- Checks on the samplesheet and settings before the run starts, and plain
  explanations when a step fails.
- A clear notice when a run ends without producing any comparison.
- Example input data and example results.

The input preparation, annotation and prediction steps are based on
[nf-core/funcscan](https://github.com/nf-core/funcscan) 3.0.0. Its AMP and
ARG screening, MultiQC, comBGC and taxonomic classification are not included.
