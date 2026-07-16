# bgc_quast_ppl

A Nextflow pipeline for **biosynthetic gene cluster (BGC) prediction and comparison**.

It takes one or more genome assemblies, predicts BGCs with three tools (**antiSMASH**, **DeepBGC**, **GECCO**), and then compares the results with **bgc-quast**. Built on the [nf-core/funcscan](https://github.com/nf-core/funcscan) 3.0.0 framework.

```
samplesheet  →  prep  →  annotation  →  BGC prediction  →  bgc-quast  →  results
                                                  (+ QUAST, in compare-to-reference mode)
```

---

## 1. What you need before you start

| Requirement | Notes |
|---|---|
| **Nextflow 25.10.5** | This exact version. Do **not** use 26.04 (it breaks this pipeline). |
| **Docker** | Must be installed and running. All tools run inside containers, so you don't install antiSMASH/DeepBGC/GECCO yourself. |
| **antiSMASH database (v8)** | Downloaded once, stored locally. See section 3. |
| **DeepBGC database** | Downloaded once, stored locally. See section 3. |
| Genome assemblies | Nucleotide FASTA files (`.fasta` or `.fasta.gz`). |

GECCO needs no external database — its model ships inside the container.

---

## 2. Get the pipeline

```bash
git clone <YOUR_REPO_URL> bgc_quast_ppl
cd bgc_quast_ppl
```

Pin the Nextflow version in **every terminal** you run the pipeline from:

```bash
export NXF_VER=25.10.5
```

Tip: run that line once per terminal session, before any `nextflow` command.

---

## 3. Databases

You need two databases on disk **before** running. Point the pipeline at them with **absolute paths** (full paths starting from `/`), not relative ones.

- **antiSMASH v8 database** → passed with `--bgc_antismash_db`
- **DeepBGC database** → passed with `--bgc_deepbgc_db`

These are downloaded using each tool's own standard download command. The exact commands depend on the tool version, so follow the official instructions:

- antiSMASH: <https://docs.antismash.secondarymetabolites.org/> (the database must match **antiSMASH 8**, which is why the folder is named `antismash_db_v8`).
- DeepBGC: <https://github.com/Merck/deepbgc> (`deepbgc download`).

Once downloaded, note the full path to each, for example:

```
~/bgc_quast_ppl/db/antismash_db_v8
~/bgc_quast_ppl/db/deepbgc_db
```

You will pass these paths on the command line (shown below).

---

## 4. The input samplesheet

The input is a CSV file. The columns are:

| Column | Meaning |
|---|---|
| `sample` | A unique name for the genome (no spaces). |
| `fasta` | Path to the genome FASTA (`.fasta` or `.fasta.gz`). |
| `type` | **Only needed in compare-to-reference mode.** `q` = query, `r` = reference. Case does not matter (`q`/`Q`, `r`/`R`). |

### For compare-samples and compare-tools

The `type` column is not needed. A two-column sheet is enough:

```csv
sample,fasta
assembly_10,~/bgc_quast_ppl/data/assembly_10.fasta.gz
assembly_20,~/bgc_quast_ppl/data/assembly_20.fasta.gz
```

### For compare-to-reference

You must include the `type` column, with **exactly one** reference row (`r`) and one or more query rows (`q`):

```csv
sample,fasta,type
assembly_10,~/bgc_quast_ppl/data/assembly_10.fasta.gz,q
assembly_20,~/bgc_quast_ppl/data/assembly_20.fasta.gz,q
reference,~/bgc_quast_ppl/data/reference.fasta.gz,r
```

The reference row's `sample` name (here `reference`) is used as the reference label in the bgc-quast report.

---

## 5. The three modes

You pick the mode with `--bgc_quast_mode`. Choose **one** per run.

| Mode | What it compares | When to use it |
|---|---|---|
| `compare-samples` *(default)* | The same tool across all your samples. You get one report per tool. | "How do my genomes compare to each other?" |
| `compare-tools` | The three tools against each other, per sample. You get one report per sample. | "For this genome, how do antiSMASH, DeepBGC, and GECCO differ?" |
| `compare-to-reference` | Your query genomes against one reference genome. QUAST runs automatically. | "How do my genomes' BGCs compare to a known reference?" |

---

## 6. How to run

Run from inside the `bgc_quast_ppl` directory. Replace the paths with your own.

### Compare-samples (default)

```bash
export NXF_VER=25.10.5

nextflow run . \
  -profile docker \
  --input data/samplesheet.csv \
  --outdir results \
  --bgc_quast_mode compare-samples \
  --bgc_antismash_db ~/bgc_quast_ppl/db/antismash_db_v8 \
  --bgc_deepbgc_db ~/bgc_quast_ppl/db/deepbgc_db \
  --max_cpus 4 \
  --max_memory 24.GB
```

### Compare-tools

Same command, change the mode:

```bash
  --bgc_quast_mode compare-tools \
```

### Compare-to-reference

Use a samplesheet that has the `type` column and one reference row (section 4):

```bash
export NXF_VER=25.10.5

nextflow run . \
  -profile docker \
  --input data/samplesheet_ref.csv \
  --outdir results \
  --bgc_quast_mode compare-to-reference \
  --bgc_antismash_db ~/bgc_quast_ppl/db/antismash_db_v8 \
  --bgc_deepbgc_db ~/bgc_quast_ppl/db/deepbgc_db \
  --max_cpus 4 \
  --max_memory 24.GB
```

QUAST runs on its own in this mode — you do not set it up.

---

## 7. Common options

All optional. Defaults are shown.

### Resources

| Option | Default | Meaning |
|---|---|---|
| `--max_cpus` | `2` | Most CPUs any single step may use. |
| `--max_memory` | `8.GB` | Most memory any single step may use. |
| `--max_time` | `24.h` | Time limit per step. |

Set these to fit your machine (e.g. `--max_cpus 4 --max_memory 24.GB`).

### Choosing tools

| Option | Effect |
|---|---|
| `--bgc_skip_antismash` | Skip antiSMASH. |
| `--bgc_skip_deepbgc` | Skip DeepBGC. |
| `--bgc_skip_gecco` | Skip GECCO. |

By default all three run.

### antiSMASH depth

antiSMASH runs in **minimal** mode by default (faster, core BGC detection only).

| Option | Effect |
|---|---|
| `--bgc_antismash_full` | Run the full analysis instead of minimal. |
| `--bgc_antismash_minimal` | Force minimal (this is already the default). |

Do **not** pass both at once — the pipeline stops with an error if you do.

### bgc-quast tuning

| Option | Default | Meaning |
|---|---|---|
| `--bgc_quast_edge_distance` | `100` | Distance (bp) from a contig edge used to call a BGC "incomplete". |
| `--bgc_quast_min_bgc_length` | `0` | Ignore BGCs shorter than this. `0` = no minimum. |
| `--bgc_quast_merge_distance` | `0` | Merge BGCs closer than this. |
| `--bgc_quast_overlap_fraction` | `0.9` | Overlap fraction used when matching BGCs. |
| `--bgc_quast_output_bgcs` | `false` | Also write the individual BGC sequences. |
| `--bgc_quast_quastdir` | — | Supply your own QUAST output directory (compare-to-reference only). QUAST is then skipped. |
| `--bgc_quast_debug` | `false` | Print the raw error output on failure (for troubleshooting). |

### Other

| Option | Default | Meaning |
|---|---|---|
| `--bgc_mincontiglength` | `3000` | Contigs shorter than this are filtered out before prediction. |
| `--save_annotations` | `false` | Keep the intermediate gene-annotation files. |

---

## 8. Where the results go

Everything lands under the folder you gave to `--outdir`.

```
results/
├── bgc_quast/
│   ├── compare_samples/            # compare-samples mode
│   │   ├── antiSMASH/report.tsv
│   │   ├── DeepBGC/report.tsv
│   │   └── GECCO/report.tsv
│   ├── compare_tools/              # compare-tools mode
│   │   └── <sample>/report.tsv
│   ├── compare_to_reference/       # compare-to-reference mode
│   │   ├── antiSMASH/report.tsv
│   │   ├── DeepBGC/report.tsv
│   │   └── GECCO/report.tsv
│   └── quast/                      # only in compare-to-reference mode
└── pipeline_info/                  # run reports, timeline, and DAG diagram
```

Only the folder for the mode you ran will be present. The main results are the `report.tsv` files.

---

## 9. Resuming a run

If a run stops partway, you can continue from where it left off by adding `-resume`:

```bash
nextflow run . -profile docker --input ... --outdir results -resume
```

Keep the same `--outdir`. Do not change database options on a `-resume`.

---

## 10. Troubleshooting

- **"command not found: nextflow"** — install Nextflow, then re-run `export NXF_VER=25.10.5`.
- **Docker errors / nothing runs** — make sure Docker is installed and running before you start.
- **Wrong Nextflow version** — this pipeline needs `25.10.5`. Newer versions (26.04+) will fail. Always set `export NXF_VER=25.10.5`.
- **Database errors from antiSMASH or DeepBGC** — check the paths you passed are **absolute** (start with `/`) and point to the correct folders. The antiSMASH database must be the **v8** database.
- **A run failed and you want detail** — add `--bgc_quast_debug` to see the raw error message.
- **"Pipeline did NOT complete successfully" (red box)** — the run finished but produced no comparison (for example, all contigs were shorter than `--bgc_mincontiglength`, so nothing reached the tools). Check your input genomes.

---

## 11. Quick reference

```bash
# 1. one-time setup
git clone <YOUR_REPO_URL> bgc_quast_ppl
cd bgc_quast_ppl
# download antiSMASH v8 and DeepBGC databases (see section 3)

# 2. every run
export NXF_VER=25.10.5

nextflow run . \
  -profile docker \
  --input data/samplesheet.csv \
  --outdir results \
  --bgc_quast_mode compare-samples \
  --bgc_antismash_db /ABS/PATH/db/antismash_db_v8 \
  --bgc_deepbgc_db /ABS/PATH/db/deepbgc_db \
  --max_cpus 4 \
  --max_memory 24.GB
```
