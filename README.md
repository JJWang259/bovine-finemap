# bovine-finemap

A pipeline for GWAS, Bayesian fine-mapping, and epigenomic enrichment analyses in cattle.
It implements four steps: (1) GWAS using SLEMM-GWA, a reliability-weighted linear mixed model approach for large-scale sequence-level association testing; (2) Bayesian fine-mapping using BFMAP forward selection; (3) SuSiE-adj fine-mapping with optional regulatory magnitude (RM) prior weights derived from epigenomic annotations; and (4) functional enrichment using GEMRICH to test whether BFMAP fine-mapped variants are enriched in tissue-specific open chromatin regions.


---

## Table of Contents
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [Input Data](#input-data)
- [Pipeline Steps](#pipeline-steps)
  - [Step 1: GWAS with SLEMM](#step-1-gwas-with-slemm)
  - [Step 2: Fine-Mapping with BFMAP](#step-2-fine-mapping-with-bfmap)
  - [Step 3: SuSiE-adj with Regulatory Magnitude Priors](#step-3-susie-adj-with-regulatory-magnitude-priors)
  - [Step 4: Functional Enrichment with GEMRICH](#step-4-functional-enrichment-with-gemrich)
- [Output Files](#output-files)
- [Publications](#publications)
- [Contact](#contact)

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
| MPH | any | Polygenic variance partitioning | [jiang18/mph](https://github.com/jiang18/mph) |
| ld_adjuster | any | LD adjustment for fine-mapping | [jiang18/ld_adjuster](https://github.com/jiang18/ld_adjuster) |
| PLINK | any | Genotype data processing | [cog-genomics.org/plink](https://www.cog-genomics.org/plink/) |
| R | ≥ 4.1 | Data processing and scripting | [r-project.org](https://www.r-project.org) |


### R Packages

```r
install.packages(c("data.table", "parallel","susieR"))

# GEMRICH (install from GitHub)
devtools::install_github("jiang18/gemrich")
```

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

### Required files

| File | Format | Description |
|------|--------|-------------|
| `geno_model.*` | PLINK binary | Model SNP genotypes for LMM fitting and BFMAP GRM construction |
| `geno_test.*` | PLINK binary | Genotypes for association testing and fine-mapping |
| `<trait>.csv` | CSV | Per-trait phenotype file (columns: IID, trait[, reliability]) |
| `snp.info.csv` | CSV | SNP name list for SLEMM variance component estimation and BFMAP GRM construction |
| `annotation.bed` | BED | Functional annotations for GEMRICH (columns: chr, start, end, category) |
| `snplist.csv` | CSV | SNP universe for GEMRICH background proportion calculation (columns: chr, pos) |

### Phenotype file format
```
IID,<trait>[,<reliability>]
HO123456,0.312,0.85
HO123457,1.120,0.91
```

The reliability column is optional. If included, specify its column name in
`ERROR_WEIGHT_NAME` in `01a_lmm_fit.sh`, `02c_estimate_heritability.sh`,
`02d_finemapping_bfmap.sh`, `03d_run_mph.sh`, and `04c_estimate_variant_effect.sh`.



### SNP info file format (`snp.info.csv`)

A single-column CSV listing SNP IDs used for variance component estimation
in SLEMM and GRM construction in BFMAP. No header required.
```
1_112_C
1_370_G
1_689_ACCTG
```

### SNP list file for GEMRICH (`snplist.csv`)

A two-column CSV containing all genotyped SNPs, used to compute background
annotation category proportions for GEMRICH enrichment analysis.

```
chr,pos
1,112
1,370
```

Generate from PLINK1 BIM: `awk '{print $1","$4}' geno.bim | sed '1s/^/chr,pos\n/' > snplist.csv`

Generate from PLINK2 PVAR: `awk 'NR>1 {print $1","$2}' geno.pvar | sed '1s/^/chr,pos\n/' > snplist.csv`

---

## Pipeline Steps

### Step 1: GWAS with SLEMM

GWAS is run in two substeps. First, a linear mixed model (LMM) is fit per trait to estimate variance components (`01a`). The fitted model is then used to run genome-wide association tests (`01b`), with results concatenated into a single per-trait file.

Edit paths at the top of each script before running.

```bash
bash scripts/01a_lmm_fit.sh
bash scripts/01b_gwas_association.sh
```

`01b` supports both merged and per-chromosome genotype files via a `{CHR}` placeholder in `GENO_TEST_TEMPLATE`.

For a full description of SLEMM parameters, see the [SLEMM documentation](https://github.com/jiang18/slemm).

**Output:** `<trait>_GWAS_All.txt`

---

### Step 2: Fine-Mapping with BFMAP

Fine-mapping is run in five substeps. Candidate regions are first identified from significant GWAS peaks (2a), a GRM is constructed for BFMAP (2b), heritability is estimated per trait (2c), fine-mapping is performed per region (2d), and results are aggregated into a single summary table (2e). Edit paths at the top of each script before running.

**2a: Identify candidate regions**
```bash
Rscript scripts/02a_identify_candidate_regions.R
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `p1` | 5e-7 | Primary significance threshold for lead variant |
| `p2` | 5e-6 | Secondary threshold for counting supporting markers |
| `scan` | 5,000,000 bp | Minimum gap to define a new region |
| `min_sig_variants` | 1 | Minimum variants passing p1 to retain a region |
| `min_secondary_variants` | 3 | Minimum variants passing p2 to retain a region |

Output: `<trait>_candidate_regions.csv`

**2b: Construct BFMAP GRM**

Note: the BFMAP GRM format is not compatible with MPH's GRM format. This step only needs to be run once per genotype dataset.

```bash
bash scripts/02b_construct_bfmap_grm.sh
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `GRM_TYPE` | 2 | GRM type: 1 = centered, 2 = standardized |

Output: `bfmap_grm.*`

**2c: Estimate heritability**

```bash
bash scripts/02c_estimate_heritability.sh
```

Output: `<trait>.varcomp.csv`

**2d: Fine-mapping**

```bash
bash scripts/02d_finemapping_bfmap.sh
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `WINDOW` | 1,000,000 bp | Defines the candidate region as ±1 Mb around the top associated position for each significant GWAS peak |

Output: per-region BFMAP results (`<region_id>.csv`) and region index (`<trait>.range.csv`) in `<trait>_bfmap/`.

For a full description of BFMAP parameters, see the [BFMAP documentation](https://github.com/jiang18/bfmap).

**2e: Aggregate fine-mapping results**

```bash
Rscript scripts/02e_summarise_finemapping.R
```
Aggregates BFMAP outputs across one or more traits into a single summary table.

Output: `fmap_all.csv`

### Step 3: SuSiE-adj with Regulatory Magnitude Priors

```bash
Rscript scripts/03_susie_adj.R
```
RM scores are available at [Zenodo](https://zenodo.org/10.5281/zenodo.12216791). Running without `--rm_scores` produces unweighted credible sets for comparison.

---

**Step 4: Functional enrichment with GEMRICH**

GEMRICH applies an MLE-based model to estimate enrichment of fine-mapped
large-effect signals across functional annotation categories.

```bash
Rscript scripts/04_gemrich_enrichment.R
```

Output: `<group>.<pv>.enrichment.csv`, `.mle.csv`, `.var.csv`, `.logll.csv`, `.counts.csv`

For a full description of GEMRICH functions and usage, see the [GEMRICH documentation](https://github.com/jiang18/gemrich).

Tissue-specific OCR BED files are available at [Zenodo](https://zenodo.org/10.5281/zenodo.12216791).
BFMAP Fine-mapping summary statistics are avaiable at [Dryad](https://datadryad.org/dataset/doi:10.5061/dryad.vmcvdnd3q).

---

## Output Files

| File | Description |
|------|-------------|
| `<trait>_GWAS_All.txt` | GWAS summary statistics across all chromosomes (01b) |
| `<trait>_candidate_regions.csv` | Candidate regions for fine-mapping (02a) |
| `<trait>_bfmap/<region_id>.csv` | Per-region BFMAP fine-mapping results (02d) |
| `fmap_all.csv` | Aggregated fine-mapping summary across all traits (02e) |
| `susie_adj_RM_weighted.txt` | SuSiE-adj PIPs with RM prior weights (03) |
| `susie_adj_cs_comparison.txt` | Credible set size comparison (03) |
| `<group>.<pv>.enrichment.csv` | GEMRICH enrichment estimates across annotation categories (04) |

---

## Publications

This pipeline was developed for and used in the following studies. If you find it helpful, please consider citing the relevant paper(s):
 
- **Dailu Guan†, Jennifer Jessica Bruscadin†, Wenjing Yang†, Claire Prowse-Wilkins†, Junjian Wang†, Bruna Petry, Houcheng Li, Huicong Zhang, Xiaoqin Xu, Ying Wang, Zhangyuan Pan, Yuri Utsunomiya, Gordon K. Murdoch, Kimberly M. Davenport, Yue Xing, Guosong Wang, Christian Maltecca, Li Ma, Gonzalo Rincon, Sarah Corum, Pengcheng Lyu, Nader Deeb, Virgínia Mara Pereira Ribeiro, Jennifer J. Michal, Zhiping Weng, Jicai Jiang, Lingzhao Fang, Honglin Jiang, Brenda M. Murdoch, Monique Rijnkels, Clare A. Gill, Timothy P. L. Smith, Zhihua Jiang, Wansheng Liu, James Reecy, Juan F. Medrano, James E. Koltes, Pablo J. Ross, Huaijun Zhou\*.** An integrated multi-tissue atlas of epigenomic landscapes and regulatory elements in the bovine genome. *Nature Communications* (under review)
- **Junjian Wang, Yahui Gao, Sajjad Toghiani, John B. Cole, Christian Maltecca, Li Ma, Jicai Jiang.** Genome-wide association study and fine-mapping using imputed sequences to prioritize candidate genes for 30 complex traits in 50,309 Holstein bulls. *Journal of Dairy Science*, 2025. https://doi.org/10.3168/jds.2025-27058

Please also cite the underlying tools:

- **SLEMM:** Cheng J, et al. (2023). *Bioinformatics*, 39(4). https://doi.org/10.1093/bioinformatics/btad127
- **BFMAP:** Jiang J, et al. (2019). *Communications Biology*, 2, 212. https://doi.org/10.1038/s42003-019-0454-y; Wang J, et al. (2025). *Briefings in Bioinformatics*, 26(6). https://doi.org/10.1093/bib/bbaf614
- **SuSiE-adj:** Wang J, et al. (2025). *Briefings in Bioinformatics*, 26(6). https://doi.org/10.1093/bib/bbaf614; Zou Y, et al. (2022). *PLoS Genetics*, 18(7), e1010299. https://doi.org/10.1371/journal.pgen.1010299; Wang G, et al. (2020). *Journal of the Royal Statistical Society Series B*, 82(5), 1273–1300. https://doi.org/10.1111/rssb.12388
- **MPH:** Jiang J. (2024). *Bioinformatics*, 40(5). https://doi.org/10.1093/bioinformatics/btae298
- **GEMRICH:** Jiang J. https://github.com/jiang18/gemrich
- **PLINK:** Chang CC, et al. (2015). *GigaScience*, 4(1). https://doi.org/10.1186/s13742-015-0047-8

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
