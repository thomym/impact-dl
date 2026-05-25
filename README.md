# Epigenomic Enrichment Pipeline

A pipeline for computing tissue/cell-type enrichment of GWAS signals using
sequence-based epigenomic predictions ([Sei](https://github.com/FunctionLab/sei-framework)),
matched background null sets, and ontology-based profile annotation.

> **Status.** Steps 1 – 4 are stable and configured from a single `paths.yaml`.
> Reference data (gnomAD, 1KG, OWL files) is downloaded by the user following
> the instructions below.

---

## Pipeline overview

```
GWAS Summary Statistics (hg19)
        │
        ▼
┌──────────────────────┐
│  gwas_preprocess/    │  QC: MAF filtering (gnomAD), strand ambiguity
│                      │  removal, imputation quality filtering
└────────┬─────────────┘
         ▼
┌──────────────────────┐
│  pre_prediction/     │  1. LD clumping (PLINK, r²≥0.69, 250 kb, p<5e-8)
│  (Snakemake)         │  2. VEP annotation → exclude coding clumps
│                      │  3. VCF export (PLINK2, hg19 ref alleles)
│                      │  4. Sample N matched background clump sets
│                      │     (default N=300; matched on chr, MAF, clump
│                      │      size, length)
│                      │  5. Convert mock datasets → VCFs
└────────┬─────────────┘
         ▼
┌──────────────────────┐
│  prediction_pipeline/│  Sei inference on GPU (real + N mocks)
│  (SLURM jobs)        │  Per profile: max |effect| per clump →
│                      │  average top-K clump scores (default K=50) →
│                      │  GWAS-level score
│                      │  Empirical p-value: fraction of mock scores ≥ real
└────────┬─────────────┘
         ▼
┌──────────────────────┐
│  ontology_matching/  │  Map Sei profiles → EFO/BTO/CL/CLO terms via
│                      │  BioPortal Annotator API, expand to ancestors,
│                      │  test enrichment (NES, MWU, Fisher) per term
└──────────────────────┘
```

All knobs in the diagram (r², kb window, p-value, N mocks, K, etc.) are
configurable in `paths.yaml`; the defaults shown match the values used in the
companion paper.

---

## Quick start

```bash
# 1. Clone (this repo + sei-framework)
git clone <your-fork-of-impact-dl>          # ← replace with the real URL
git clone https://github.com/FunctionLab/sei-framework

# 2. Create environments
micromamba env create -f impact-dl/m-env.yaml
micromamba env create -f impact-dl/sei-env.yaml

# 3. Download reference data — see "Reference data" below

# 4. Configure paths
cd impact-dl
cp paths.yaml.example paths.yaml
$EDITOR paths.yaml                          # set work_dir + gwas_names

# 5. Run step 1 (per GWAS), then step 2 (Snakemake)
micromamba activate m-env
python gwas_preprocess/prepare_gwas.py     --discovery <gwas_name>
bash   gwas_preprocess/qc_discovery_data.sh --gwas_name <gwas_name>
cd pre_prediction && snakemake --cores 8 --snakefile Snakefile
```

Full details for each step are below.

---

## Requirements

### Software

| Tool | Version tested | Purpose |
|------|---------------|---------|
| PLINK | 1.9 | LD clumping, bfile filtering |
| PLINK2 | 2.0 | VCF export with reference allele correction |
| bcftools | ≥ 1.14 | VCF processing |
| Ensembl VEP | `ensemblorg/ensembl-vep:release_113.0` | Coding variant annotation (containerized) |
| Snakemake | ≥ 7.0 | `pre_prediction/` workflow |
| SLURM | any | `prediction_pipeline/` job submission |
| udocker or docker | recent | Running the VEP container |

### Python environments

Two micromamba environments are provided at the repo root:

- **`m-env.yaml`** (env name: `m-env`) — steps 1, 2, 3c, 3d, 4
- **`sei-env.yaml`** (env name: `sei-env`) — steps 3a, 3b (Sei inference, requires GPU)

```bash
micromamba env create -f m-env.yaml
micromamba env create -f sei-env.yaml

micromamba activate m-env       # steps 1, 2, 3c, 3d, 4
micromamba activate sei-env     # steps 3a, 3b
```

---

## Repository structure

```
impact-dl/
├── paths.yaml.example          ← template (copy to paths.yaml and edit)
├── paths.py                    ← shared yaml loader with ${var} interpolation
├── m-env.yaml, sei-env.yaml    ← conda environments
│
├── gwas_preprocess/            ← Step 1
│   ├── prepare_gwas.py
│   ├── qc_discovery_data.sh
│   └── parse_args.sh
│
├── pre_prediction/             ← Step 2 (Snakemake)
│   ├── Snakefile
│   ├── gen_pre_clumping.sh
│   ├── main_gwas_to_vcf.sh
│   ├── mock_dataset_to_vcf.sh
│   ├── sample_null_blocks.py
│   └── rename_chrs.txt
│
├── prediction_pipeline/        ← Step 3
├── ontology_matching/          ← Step 4
│   ├── 1-four_ontologies.py
│   ├── 2-ancestor_and_collapse.py
│   ├── 3-enrichment_per_gwas.py
│   └── functions_for_snp_gsealike.py
└── README.md
```

---

## Directory layout (user-side)

This pipeline assumes one top-level **`work_dir`**. The default layout is:

```
<work_dir>/
├── impact-dl/                       ← this repo
├── sei-framework/                   ← cloned from FunctionLab/sei-framework
│   ├── resources/
│   │   └── hg19_UCSC.fa             ← from sei-framework's download_data.sh
│   └── model/
│       └── target.names
│
├── resources/                       ← reference data (see "Reference data")
│   ├── gnomad/<version>/
│   ├── 1kg/
│   ├── vep_cache/
│   └── ontology/
│
├── gwas/                            ← one subdir per GWAS (inputs)
│   └── <gwas_name>/
│       └── gwas_raw.tsv
│
└── results/                         ← one subdir per GWAS (outputs)
    └── <gwas_name>/
```

You don't have to put everything inside `work_dir`. Every path is overridable
in `paths.yaml` — e.g. if `sei-framework` already lives at
`/vol/scratch/shared/sei-framework`, point `sei_framework:` at it directly.

---

## Reference data

| What | Source | Goes to | Size |
|------|--------|---------|------|
| Sei model weights + hg19 FASTA + `target.names` | sei-framework's `download_data.sh` | `<work_dir>/sei-framework/` | ~5 GB |
| gnomAD MAF tables | built from a gnomAD release (see below) | `<resources>/gnomad/<version>/` | varies |
| 1KG PLINK bfiles + EUR panel | 1000 Genomes FTP | `<resources>/1kg/` | ~30 GB |
| VEP offline cache (GRCh37) | Ensembl VEP `INSTALL.pl` | `<resources>/vep_cache/` | ~25 GB |
| Ontology OWL files | OLS / OBO Foundry | `<resources>/ontology/` | a few MB each |

### 1. Sei model + hg19 FASTA

```bash
cd <work_dir>/sei-framework/resources
bash download_data.sh
```

This populates `sei-framework/resources/hg19_UCSC.fa` and `sei-framework/model/target.names`.

### 2. gnomAD MAF tables

The pipeline reads per-chromosome TSVs with columns `CHROM POS ID REF ALT MAF`,
plus one file concatenated across all chromosomes. The `MAF` column is the
European allele frequency from gnomAD. Produce these from a gnomAD release
using any tool that writes those columns (it gets folded to the true
minor-allele frequency downstream).

Set `gnomad_version` in `paths.yaml` to match your gnomAD directory name. The
per-chromosome files feed `gwas_preprocess/prepare_gwas.py`; the concatenated
file feeds `pre_prediction/gen_pre_clumping.sh`.

### 3. 1KG PLINK bfiles + EUR panel

```bash
# Phased 1KG VCFs → PLINK bfiles (adjust to your preferred 1KG release)
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/...
plink --vcf ALL.chr*.vcf.gz --make-bed --out <resources>/1kg/ds

# Population panel (used to restrict LD to EUR samples)
wget -P <resources>/1kg/ \
  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```

### 4. VEP offline cache (GRCh37)

Follow the [standard VEP cache install](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html),
but install into `<resources>/vep_cache/` so the pipeline scripts can mount it
into the container at `/opt/vep/.vep`:

```bash
mkdir -p <resources>/vep_cache
udocker pull ensemblorg/ensembl-vep:release_113.0
udocker run \
  -v <resources>/vep_cache:/opt/vep/.vep \
  ensemblorg/ensembl-vep:release_113.0 \
  INSTALL.pl --AUTO cf -s homo_sapiens -y GRCh37 --CACHEDIR /opt/vep/.vep
```

Substitute `docker` for `udocker` if you use Docker.

### 5. Ontology OWL files *(needed only for step 4)*

Download `efo.owl`, `bto.owl`, `cl.owl`, `clo.owl` from
[OLS](https://www.ebi.ac.uk/ols4) or [OBO Foundry](https://obofoundry.org/)
into `<resources>/ontology/` (paths overridable via `efo_owl` / `bto_owl` /
`cl_owl` / `clo_owl` in `paths.yaml`).

Step 4 also reads two files that live in the repo under
`ontology_matching/data/` (resolved via `${repo}` in `paths.yaml`):

- `sup_table_2_sei_article.csv` — standardized cell-type names (`${sup2_file}`).
- `bioportal_cache/` — pre-computed BioPortal annotation caches
  (`${bioportal_cache_dir}`). With these present, `1-four_ontologies.py` reads
  cached annotations instead of querying the BioPortal API.

---

## Configuration

All filesystem paths and pipeline knobs live in a single file at the repo
root: **`paths.yaml`**.

```bash
cd <work_dir>/impact-dl
cp paths.yaml.example paths.yaml
$EDITOR paths.yaml          # at minimum: set `work_dir` and `gwas_names`
```

How it works:

- Any value can reference an earlier one with `${var}` — e.g.
  `gwas_root: ${work_dir}/gwas`.
- Any line can be replaced with an absolute path to point outside
  `work_dir` — e.g. `sei_framework: /vol/scratch/shared/sei-framework`.
- All scripts that need paths import the shared loader at
  [`paths.py`](paths.py), so behavior is consistent across stages.
- Override the yaml file location with:
  - `--paths_yaml /other/paths.yaml` (Python scripts), or
  - `snakemake --config paths_yaml=/other/paths.yaml` (Snakemake).

`paths.yaml` is gitignored; commit only `paths.yaml.example`.

You can also debug the resolved config:

```bash
python paths.py --paths_yaml paths.yaml
```

---

## Step-by-step usage

### Step 1 — Prepare GWAS (`gwas_preprocess/`) — *optional*

> **Optional.** This step is a convenience preprocessor that turns a raw GWAS
> summary-statistics file into the `gwas.QC.Transformed` format the rest of the
> pipeline consumes. If you already have QC'd summary statistics, **skip this
> step** and provide a file matching the input contract below.

Input:  `${gwas_root}/${gwas_name}/gwas_raw.tsv`
Output: `${gwas_root}/${gwas_name}/gwas.QC.Transformed`

#### Input contract — `gwas.QC.Transformed`

If you bypass Step 1, your file must be a tab-separated table with this column
order (the rest of the pipeline reads columns positionally):

```
SNP  CHR  BP  A1  A2  MAF  SE  P  N  INFO  BETA  OR
```

- `SNP` = rsID, `CHR` = 1–22, `BP` = hg19 position, `A1` = effect allele,
  `A2` = other allele, `MAF` = minor-allele frequency, `BETA`/`OR` = effect.
- The QC filters Step 1 applies (you should apply equivalents): `MAF > 0.01`,
  `INFO > 0.8`, no strand-ambiguous SNPs (A/T, G/C), no duplicate SNP IDs.

**1a. Reformat headers and merge MAF from gnomAD**

```bash
micromamba activate m-env

python gwas_preprocess/prepare_gwas.py --discovery <gwas_name>
```

`work_dir`, `auxiliary_path` and `version` default to `gwas_root`, `resources`
and `gnomad_version` from `paths.yaml`. Override any of them with the
corresponding flag (`--work_dir`, `--auxiliary_path`, `--version`) or point
at a different yaml with `--paths_yaml`. Other knobs (`--discovery_population
EUR`, `--default_maf 0.05`, `--maf_th 0`, `--N <sample_size>`, `--chrs 1,2,...`)
have sensible defaults; see `--help`.

`gwas_raw.tsv` is the original downloaded GWAS summary statistics file.
The script auto-detects and renames non-standard column headers (see
[`prepare_gwas.py`](gwas_preprocess/prepare_gwas.py) for the full mapping)
and merges MAF from the per-chromosome gnomAD tables under
`<resources>/gnomad/<version>/`.
Intermediate output: `${gwas_root}/${gwas_name}/gwas.tsv`.

**1b. Quality control filtering**

```bash
bash gwas_preprocess/qc_discovery_data.sh --gwas_name <gwas_name>
```

`work_dir` defaults to `gwas_root` from `paths.yaml`; override with
`--work_dir <path>` or `--paths_yaml <other.yaml>`.

Filters:
- MAF > 0.01
- INFO (imputation quality) > 0.8
- Remove strand-ambiguous SNPs (A/T, G/C)
- Remove duplicate SNP IDs

Final output: `${gwas_root}/${gwas_name}/gwas.QC.Transformed`.

---

### Step 2 — Pre-prediction (`pre_prediction/`)

Reads `<repo>/paths.yaml` for all paths and pipeline knobs.

```bash
micromamba activate m-env

cd pre_prediction
snakemake --cores 8 --snakefile Snakefile -n    # dry-run (recommended first)
snakemake --cores 8 --snakefile Snakefile       # real run
```

Use a different config file:

```bash
snakemake --cores 8 --snakefile Snakefile \
          --config paths_yaml=/path/to/other.yaml
```

The workflow runs four rules in order, once per GWAS listed under
`gwas_names` in `paths.yaml`:

| Rule | Output |
|------|--------|
| `pre_clumping` | Filtered 1KG bfiles + `gwas_for_clumping.txt` |
| `main_gwas_to_vcf` | `${results_root}/${gwas}/clumping_and_vcfs_outputs/main_gwas_snps_correct_ref_final.vcf` |
| `sample_null_blocks` | N mock datasets at `${results_root}/${gwas}/mock_datasets/mock_dataset{1..N}.tsv` |
| `mock_to_vcf` | `${results_root}/${gwas}/vcfs_mock_datasets/mock_dataset{idx}/mock_dataset{idx}_snps_correct_ref_final.vcf` |

`sample_null_blocks` is memory-heavy; the Snakefile caps it at one concurrent
invocation via a `sample_null_blocks_concurrent` resource. Raise it only if
you have enough RAM.

VEP runs via the container declared in `paths.yaml`
(`vep_image: ensemblorg/ensembl-vep:release_113.0`,
`container_runtime: udocker|docker`). The host's `vep_cache` is mounted into
the container as `/opt/vep/.vep`.

Run on a SLURM cluster:

```bash
snakemake --cores 8 --snakefile Snakefile \
    --cluster "sbatch --mem={resources.mem_mb}M -c {threads}" \
    --jobs 50
```

(Snakemake ≥ 8.0 also supports `--executor slurm` via the
`snakemake-executor-plugin-slurm` plugin.)

---

### Step 3 — Sei predictions (`prediction_pipeline/`)

Runs Sei on the real (clumped) GWAS VCF and on every mock VCF, then
aggregates the per-SNP effects into per-profile scores and computes empirical
p-values.

**Disk note:** peak usage is ~600 MB × `n_mock_datasets` (≈ 180 GB at the
default of 300 mocks). Stored under `${results_root}/${gwas_name}/sei_outputs_*/`.
Ref/alt prediction h5s are deleted after each Sei run; only the diffs h5 is kept.

**Two SLURM jobs (GPU), then two Python scripts (CPU, run locally — NOT on cluster):**

> **Submit pattern**: cd into `impact_dl/` before `sbatch` so that
> `$SLURM_SUBMIT_DIR` points at the repo (this is how the slurm scripts find
> `paths.py`). If you'd rather sbatch from anywhere, pass
> `--repo_root /path/to/impact_dl` as a script arg.

#### 3a. Sei on the real GWAS

```bash
cd /path/to/impact_dl
sbatch --partition=<your_partition> \
       prediction_pipeline/1-sei_main.slurm \
       --gwas_name <gwas_name>
```

Output: `${results_root}/${gwas_name}/sei_outputs_main_gwas/chromatin-profiles-hdf5/main_gwas_snps_correct_ref_final_diffs.h5`.
Re-submitting is safe: if the diffs h5 already exists, the script exits early.
Logs land in `${results_root}/${gwas_name}/slurm_logs/`.

#### 3b. Sei on the mocks (SLURM array, batched)

Most clusters limit the total *queued* (not running) jobs per user — typically
~100 via the `normal` QOS. Submitting all 300 mock array tasks at once usually
hits `QOSMaxSubmitJobPerUserLimit`. We provide a thin batching wrapper that
submits chunks of `BATCH_SIZE` (default 80, safely under the limit) one at a
time, waiting for each to drain.

**Run inside `tmux` or `screen`** — this blocks until all 300 mocks finish
(several hours typically), with progress every 30 s:

```bash
cd /path/to/impact_dl
tmux            # or screen
bash prediction_pipeline/2-sei_mocks.sh \
    --gwas_name <gwas_name> \
    --partition <your_partition>
```

Override defaults via flags:
- `--batch_size 80` — array chunk submitted per batch (raise if your QOS allows)
- `--array_concurrency 20` — the `%N` cap on running tasks per array
- `--total <N>` — defaults to `n_mock_datasets` from `paths.yaml`

Each per-task script is resumable: completed mocks are skipped via the
diffs.h5 existence check, so killing and re-running `submit_mocks.sh` is safe.

**Bare sbatch** (single batch, no wrapper) is fine for small N:
```bash
sbatch --partition=<your_partition> --array=1-50%20 \
       prediction_pipeline/_sei_mocks_array.slurm --gwas_name <gwas_name>
```

#### 3c. Aggregate mock scores

```bash
micromamba activate m-env
python impact-dl/prediction_pipeline/3-aggregate_mocks.py --gwas_name <gwas_name>
```

Reads all `n_mock_datasets` Sei h5 files, computes per-clump max → top-K → mean
to produce one 21,907-vector per mock, saved as
`${results_root}/${gwas_name}/sei_score_aggregations_all/mock_scores_aggregation_50.csv`.

By default uses up to 8 parallel workers (each loads one h5 — ~600 MB);
override with `--n_workers <N>` if you have more memory headroom.

#### 3d. Empirical p-values (real vs. mock)

```bash
micromamba activate m-env
python impact-dl/prediction_pipeline/4-empirical_pvalues.py --gwas_name <gwas_name>
```

Produces:
- `empirical_pvalues50.csv` — per-profile `Empirical_P_Upper`, `Empirical_P_Lower`, `Empirical_P_Two_Sided`, `Effect_Size`.
- `top_50_snps_by_annotation_main.json` — top-50 SNPs per Sei profile.
- `all_snp_scores_.json` — per-SNP-per-profile scores (large).
- `empirical_pvalues50.versions.txt` — sidecar recording the numpy/scipy/conda-env that produced the CSV.

> **Environment matters.** Run steps 3c and 3d in the **same** Python
> environment for reproducibility. Numerical results can drift slightly
> across numpy/scipy versions; the `.versions.txt` sidecar records what was
> used.

#### Debugging helper

`view_diffs_h5.py` loads a Sei diffs h5 + row labels + `target_names` into a
DataFrame for inspection. Default `target_names` comes from `paths.yaml`.

```bash
python impact-dl/prediction_pipeline/view_diffs_h5.py \
    ${results_root}/${gwas_name}/sei_outputs_main_gwas/chromatin-profiles-hdf5/main_gwas_snps_correct_ref_final_diffs.h5
```

---

### Step 4 — Ontology matching & enrichment (`ontology_matching/`)

Maps Sei sequence profiles to biological ontology terms (EFO, BTO, CL, CLO)
via the BioPortal Annotator API, expands to ancestor terms, then runs
GSEA-like enrichment (NES + MWU + Fisher) per term.

```bash
micromamba activate m-env
cd /path/to/impact_dl/ontology_matching
```

The first two scripts are **GWAS-independent** — run once for the whole
pipeline. The third is **per-GWAS**.

#### 4.1 `1-four_ontologies.py` — sup2 → BioPortal terms (run once)

Annotates each cell-type name in `sup_table_2` against EFO, CLO, CL, and BTO
via the BioPortal Annotator API.

```bash
export BIOPORTAL_API_KEY=...        # only needed if a cache miss forces an API call
python 1-four_ontologies.py
```

The Annotator API is rate-limited (queries can take hours for a fresh sup2).
Pre-computed caches (`annotation_cache_*.pkl` + `bioportal_class_cache.json`)
are in the repo at `${bioportal_cache_dir}`; the script reads them instead of
querying the API. To regenerate, set `BIOPORTAL_API_KEY` and delete the cache
files.

Output: `${ontology_matched_dir}/matched_sup2_<ont>.csv` for each ontology.

#### 4.2 `2-ancestor_and_collapse.py` — expand + collapse (run once)

Loads the OWL files, expands each matched term to its ancestors (minus
diseases / generic terms), then collapses any chain whose descendants share
the same annotation coverage (keeps the leaf).

```bash
python 2-ancestor_and_collapse.py
```

OWL files default to `${ontology_dir}/{efo,bto,cl,clo}.owl` — download from
[OLS](https://www.ebi.ac.uk/ols4) (a few MB each).

Output: `${ontology_collapsed_dir}/<ont>_term_to_annotations_collapsed.csv`
plus a per-ontology `*_collapse_log.csv`.

#### 4.3 `3-enrichment_per_gwas.py` — per-GWAS enrichment

For one GWAS, tests each collapsed term for enrichment using:

- one-sided Mann-Whitney U (`MWU_PValue`) on `Empirical_P_Upper`
- GSEA-like running ES + NES + permutation p-value (`NES`, `ES_PValue`) —
  seeded for reproducibility (`--seed`, default 42)

```bash
python 3-enrichment_per_gwas.py --gwas_name <gwas_name>
```

Pre-flight assertions (strict — refuses to run otherwise):
- `${results_root}/${gwas_name}/sei_score_aggregations_all/empirical_pvalues50.csv` exists
- `vcfs_mock_datasets/` contains exactly `n_mock_datasets` subdirs
- `clumps_cleaned_noncoding.clumped` has ≥ 50 clumps

Outputs (under `${results_root}/${gwas_name}/`):
- `ontology_enrichment_<ont>.csv` — one row per testable term
- `term_to_annotation_mappings_<ont>.csv` — flat term→annotation map

Useful flags: `--ontologies EFO CL` (subset), `--min_cells N` (override
`min_cells_per_term`), `--seed N` (RNG seed).

A free BioPortal account is required only to refresh the matched CSVs from a
new sup2 file: https://bioportal.bioontology.org
