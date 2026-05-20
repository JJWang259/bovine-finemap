# bovine-finemap

A pipeline for GWAS, Bayesian fine-mapping, and epigenomic enrichment analyses in cattle.

---

## Table of Contents
- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Input Data](#input-data)
- [Pipeline Steps](#pipeline-steps)
  - [Step 1: GWAS with SLEMM](#step-1-gwas-with-slemm)
  - [Step 2: Fine-Mapping with BFMAP](#step-2-fine-mapping-with-bfmap)
  - [Step 3: Functional Enrichment](#step-3-functional-enrichment)
- [Output Files](#output-files)
- [Publications](#publications)
- [Contact](#contact)

---

## Overview

This repository implements a three-step pipeline for identifying and functionally interpreting causal genetic variants underlying complex traits in cattle:

1. **GWAS** using SLEMM-GWA, a reliability-weighted linear mixed model approach for large-scale sequence-level association testing
2. **Bayesian fine-mapping** using BFMAP and SuSiE-adj, with optional regulatory magnitude (RM) prior weights derived from epigenomic annotations across 55 bovine tissue contexts
3. **Functional enrichment** using GEMRICH to test whether fine-mapped variants are enriched in tissue-specific open chromatin regions (OCRs)

---

## Repository Structure

```
bovine-finemap/
├── README.md
└── scripts/
    ├── 01a_lmm_fit.sh                     # Step 1a: LMM model fitting (SLEMM)
    ├── 01b_gwas_association.sh            # Step 1b: GWAS (SLEMM-GWA)
    ├── 02a_identify_candidate_regions.R   # Step 2a: Define candidate regions
    ├── 02b_construct_bfmap_grm.sh         # Step 2b: Construct BFMAP GRM
    ├── 02c_estimate_heritability.sh       # Step 2c: Estimate heritability with BFMAP
    ├── 02d_finemapping_bfmap.sh           # Step 2d: BFMAP forward selection fine-mapping
    ├── 02e_summarise_finemapping.R        # Step 2e: Aggregate fine-mapping results
    ├── 03_susie_adj.R                     # Step 3: SuSiE-adj with optional RM prior weights
    └── 04_gemrich_enrichment.R            # Step 4: OCR enrichment analysis
```

---

## Dependencies


### Software

| Tool | Version | Purpose | Link |
|------|---------|---------|------|
| SLEMM | any | GWAS | [jiang18/slemm](https://github.com/jiang18/slemm) |
| BFMAP | ≥ 0.65 | Bayesian fine-mapping | [jiang18/bfmap](https://github.com/jiang18/bfmap) |
| GEMRICH | any | Large-effect enrichment analysis | [jiang18/gemrich](https://github.com/jiang18/gemrich) |
| PLINK | any | Genotype data processing | [cog-genomics.org/plink](https://www.cog-genomics.org/plink/) |
| R | ≥ 4.1 | Data processing and scripting | [r-project.org](https://www.r-project.org) |

---

## Input Data

| Data | Source | Notes |
|------|--------|-------|
| Genotypes (50,309 bulls) | Not publicly released | Contact authors |
| Phenotypes (30 traits, DRP) | Not publicly released | Contact authors |
| Fine-mapping summary stats | [Dryad](https://datadryad.org/dataset/doi:10.5061/dryad.vmcvdnd3q) | PCIP for all 30 traits |
| Tissue-specific OCR BEDs | [Zenodo](https://zenodo.org/10.5281/zenodo.12216791) | 47 tissues |
| Regulatory magnitude scores | [Zenodo](https://zenodo.org/10.5281/zenodo.12216791) | 28M SNPs |
| gkm-SVM weights | [Zenodo](https://zenodo.org/10.5281/zenodo.12216791) | 206 tissue-mark contexts |
| Reference genome (ARS-UCD1.2) | [Ensembl v105](https://ensembl.org) | — |

---

## Pipeline Steps

### Step 1: GWAS with SLEMM

```bash
bash scripts/gwas/run_slemm_gwas.sh \
  --geno  data/genotypes \
  --pheno data/phenotypes.txt \
  --trait MilkYield \
  --out   results/gwas/MilkYield

Rscript scripts/gwas/identify_peaks.R \
  --gwas    results/gwas/MilkYield_gwas.txt \
  --pthresh 5e-8 \
  --out     results/gwas/MilkYield
```

### Step 2: Fine-Mapping with BFMAP

```bash
# Option A: Run BFMAP from scratch
bash scripts/finemapping/run_bfmap.sh \
  --gwas    results/gwas/MilkYield_gwas.txt \
  --regions results/gwas/MilkYield_regions.txt \
  --geno    data/genotypes \
  --pheno   data/phenotypes.txt \
  --trait   MilkYield \
  --out     results/finemapping/MilkYield

# Option B: Use published fine-mapping summary statistics
Rscript scripts/utils/parse_dryad_finemapping.R \
  --input finemapping_summary_stats_dairy_2025.csv \
  --out   results/finemapping/MilkYield_parsed \
  --trait MilkYield

# SuSiE-adj with regulatory magnitude prior weights
Rscript scripts/finemapping/run_susie_adj.R \
  --bfmap_out results/finemapping/MilkYield_bfmap.txt \
  --rm_scores data/regulatory_magnitude_scores.txt \
  --out       results/finemapping/MilkYield_susie_adj
```

### Step 3: Functional Enrichment

```bash
Rscript scripts/enrichment/run_gemrich_enrichment.R \
  --finemapping_stats results/finemapping/MilkYield_bfmap.txt \
  --ocr_dir           data/tissue_specific_OCRs/ \
  --trait             MilkYield \
  --out               results/enrichment/MilkYield_OCR
```

Or run all steps at once:
```bash
cp config_template.sh config.sh   # edit paths and traits
bash run_pipeline.sh --config config.sh
```

---

## Output Files

| File | Description |
|------|-------------|
| `results/gwas/<trait>_gwas.txt` | GWAS summary statistics |
| `results/gwas/<trait>_peaks.txt` | Lead SNPs per independent peak |
| `results/gwas/<trait>_regions.txt` | 2-Mb candidate windows for fine-mapping |
| `results/finemapping/<trait>_all_signals.txt` | BFMAP posterior inclusion probabilities (PCIP) |
| `results/finemapping/<trait>_susie_adj_RM_weighted.txt` | SuSiE-adj PIPs with RM prior weights |
| `results/finemapping/<trait>_high_confidence_variants.txt` | Variants with PIP > 0.5 |
| `results/enrichment/<trait>_OCR_enrichment_results.txt` | Per-tissue enrichment statistics |
| `results/figures/` | Locus plots and credible set comparison figures |

---

## Publications

If you use these scripts, please cite:

```
Guan D*, Bruscadin JJ*, Yang W*, Prowse-Wilkins C*, Wang J*, et al. (2025).
An integrated multi-tissue atlas of epigenomic landscapes and regulatory
elements in the bovine genome. Nature Genetics *(under review)*

Wang J, Gao Y, Toghiani S, Cole JB, Maltecca C, Ma L, Jiang J. (2025).
Genome-wide association study and fine-mapping using imputed sequences to
prioritize candidate genes for 30 complex traits in 50,309 Holstein bulls.
Journal of Dairy Science. https://doi.org/10.3168/jds.2025-27058
```

---

## Contact

**Junjian Wang**  
Department of Animal Science, North Carolina State University  
📧 jwang259@ncsu.edu  
🔗 [jjwang259.github.io](https://jjwang259.github.io/)

**Jicai Jiang**  
Department of Animal Science, North Carolina State University  
📧 jicai_jiang@ncsu.edu  
🔗 [cals.ncsu.edu/animal-science/people/jicai-jiang](https://cals.ncsu.edu/animal-science/people/jicai-jiang/)

For bug reports or questions, please open an [Issue](https://github.com/JJWang259/bovine-finemap/issues).
