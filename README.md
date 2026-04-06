# IMPACT-DL

This is a preliminary code version

>

## Overview

IMPACT-DL identifies cell types and tissues enriched for predicted noncoding regulatory effects in GWAS data. It couples the [Sei](https://github.com/FunctionLab/sei-framework) deep learning model (which predicts variant effects across 21,907 chromatin profiles) with LD-aware aggregation and ontology-guided organization to produce interpretable maps of trait-relevant cellular contexts.

## Architecture & Data Flow

```
GWAS Summary Statistics (33 traits, European ancestry, hg19)
        │
        ▼
┌──────────────────────┐
│  gwas_preprocess/     │  QC: MAF filtering (gnomAD v4.1), strand ambiguity
│                       │  removal, imputation quality filtering
└────────┬─────────────┘
         ▼
┌──────────────────────┐
│  pre_prediction/      │  1. LD clumping (PLINK, r²≥0.69, 250kb, p<5e-8)
│  (Snakemake)          │  2. VEP annotation → exclude coding clumps
│                       │  3. VCF export (PLINK2, hg19 ref alleles)
│                       │  4. Sample 300 matched background clump sets
│                       │     (matched on chr, MAF±0.1, clump size±15, length±20%)
│                       │  5. Convert 300 mock datasets → VCFs
└────────┬─────────────┘
         ▼
┌──────────────────────┐
│  prediction_pipeline/ │  Sei inference on GPU (real + 300 mocks)
│  (SLURM jobs)         │  Per profile: max |effect| per clump →
│                       │  average top 50 clump scores → GWAS-level score
│                       │  Empirical p-value: fraction of mock scores ≥ real
│                       │  
└────────┬─────────────┘
         ▼
┌──────────────────────┐
│  ontology_matching/   │  Map 21,907 Sei profiles → EFO/BTO/CL/CLO terms
│                       │  via BioPortal Annotator API (top 4 matches)
│                       │  
└────────┬─────────────┘
         ▼
┌──────────────────────┐
│  for_figures/         │  Per term: partition profiles into in-set vs. rest
│                       │  GSEA-like ES + MWU test + Fisher's exact test
│                       │  
└──────────────────────┘
```

## Directory Structure

| Directory | Purpose |
|---|---|
| `gwas_preprocess/` | GWAS QC, gnomAD MAF annotation, allele harmonization |
| `pre_prediction/` | Snakemake pipeline: LD clumping, VEP filtering, null set generation, VCF creation |
| `prediction_pipeline/` | Sei GPU inference + streaming score aggregation via SLURM |
| `ontology_matching/` | BioPortal-based ontology matching, ancestor expansion, enrichment testing |
| `for_figures/` | GSEA-like enrichment functions, top SNP extraction, visualization |

## Setup

### System Requirements

- **Linux** 
- **SLURM** cluster with GPU nodes (for Sei inference)
- **CUDA**-capable GPU

### External Tools

| Tool | Version | Purpose |
|---|---|---|
| [PLINK 1.9](https://www.cog-genomics.org/plink/) | 1.9 | LD clumping (r²≥0.69), LD table computation |
| [PLINK 2](https://www.cog-genomics.org/plink/2.0/) | 2.0 | VCF export with `--ref-from-fa` (hg19 ref alleles) |
| [bcftools](http://www.htslib.org/) | ≥1.10 | VCF manipulation, ref allele correction |
| [Snakemake](https://snakemake.readthedocs.io/) | ≥7.0 | Pipeline orchestration |
| [udocker](https://indigo-dc.github.io/udocker/) | ≥1.3 | Rootless container execution (for VEP + bcftools on HPC) |
| [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/) | — | Variant consequence annotation (coding clump exclusion) | called in the code vep_container, run with udocker (in the pre_processing files, installation specified below)
| [Sei framework](https://github.com/FunctionLab/sei-framework) | — | Sequence-based variant effect prediction (21,907 chromatin profiles) |

### udocker + VEP Setup

VEP is used to identify protein-altering variants so that LD clumps containing coding variants can be excluded. It runs via [udocker](https://indigo-dc.github.io/udocker/) on HPC nodes without root access.

```bash
pip install udocker

# VEP container
udocker pull ensemblorg/ensembl-vep
udocker create --name=vep_container ensemblorg/ensembl-vep
```

VEP requires a local cache directory (expected at `/opt/vep/.vep` inside the container). Download from [Ensembl](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html).

### Python Environments

The pipeline requires **two separate conda/micromamba environments** due to Sei's dependency constraints.

**1. `sei-env`** — Sei GPU inference (Python 3.6)

Used by the SLURM jobs in `prediction_pipeline/` to run Sei variant effect predictions.

```bash
micromamba create -f sei-env.yaml
micromamba activate sei-env
```

**2. `m-env`** — Everything else (Python 3.9)

Used for GWAS preprocessing, Snakemake orchestration, score aggregation, ontology matching, enrichment analysis, and figure generation. Includes `bcftools`, `snakemake`, and all Python dependencies.

```bash
micromamba create -f m-env.yaml
micromamba activate m-env
```

> ⚠️ The two environments use different Python versions and incompatible pickle protocols. Sei inference must run in `sei-env`; all other scripts must run in `m-env`.

### Reference Data

- **1000 Genomes Phase 3** — PLINK binary files (`.bed/.bim/.fam`) + European population panel
- **gnomAD v4.1** — Per-chromosome European allele frequencies. Download from [gnomAD downloads](https://gnomad.broadinstitute.org/downloads). The pipeline expects a combined file with columns `CHROM, POS, REF, ALT, MAF`.
- **hg19 FASTA** — Included with the Sei framework (`resources/hg19_UCSC.fa`)
- **Ontology files** — EFO, BTO, CL, CLO in OWL format from their respective sources
- **[BioPortal API key](https://bioportal.bioontology.org/accounts/new)** — For ontology annotation matching

### Configuration

Edit `pre_prediction/config.yaml` with your local paths:

```yaml
gwas_root: "/path/to/gwas/summary/stats"
results_root: "/path/to/output"
kg_bfile: "/path/to/1kg/plink_prefix"
gnomad_file: "/path/to/gnomad/v4.1/chr_all_eur_maf"
ref_fasta: "/path/to/sei-framework/resources/hg19_UCSC.fa"
pvalue_threshold: 5e-8
n_mock_datasets: 300
```

Also update hardcoded paths in `ontology_matching/` and `prediction_pipeline/` scripts, and set your BioPortal API key:

```bash
export BIOPORTAL_API_KEY="your_key_here"
```

