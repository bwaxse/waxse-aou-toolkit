# Genomics Analysis - Ancestry PCA & Genotype Preparation

**Research Focus**: Genomics infrastructure and population structure analysis using All of Us genomic data

**Platform**: All of Us Research Workbench (OMOP CDM v8) with Hail genomics framework
**Primary Methods**: Ancestry inference via PCA, genotype quality control, VCF preparation
**Language**: Python (Jupyter notebooks) with R (optional visualization)

---

## Getting Started

**BEFORE running any notebooks in this project:**

Run `_reference/verily/00_setup_workspace.ipynb` once per JupyterLab session
- This initializes `WORKSPACE_CDR`, `WORKSPACE_BUCKET`, etc.
- Required for all subsequent notebooks to work
- See main `CLAUDE.md` for full setup details

Once environment variables are set, follow the analysis pipeline in order (see workflow below).

---

## Project Workflow Overview

This project provides shared genomics infrastructure used by downstream GWAS analyses (HPV, Sarcoid, etc.).

The pipeline focuses on:
1. **Ancestry inference** via principal component analysis (PCA)
2. **Genotype quality control** and format conversion
3. **Reference data preparation** for use in multi-ancestry GWAS

```
Z.001: Prepare Ancestry PGen → Z.003: Prepare Pruned BED →
Z.005: Merge & Run PCA → Z.007: Assess & Format PCA →
Z.009: Prepare Ancestry Genotypes & VCF
```

Additionally:
- **001**: Cohort definition with Hail framework
- **G.11**: GPR15 gene-specific cohort analysis (GI inflammation case study)

---

## File Map by Analysis Stage

### Z.001: Prepare Ancestry PGen
Converts All of Us genomic data to Plink format for ancestry analysis.

**Files:**
- `Z.001_prepare_ancestry_pgen.py` / `.ipynb` - Primary conversion
- `Z.001_prepare_ancestry_pgen_for_pcs.py` / `.ipynb` - PCA-specific variant extraction

**Workflow:**
1. Load All of Us genomic files (`.mt` Hail MatrixTable format)
2. Filter to ancestry-informative variants
3. Export to Plink format (`.bed`, `.bim`, `.fam`)

**Outputs:**
- Plink genotype files for downstream analysis
- Ancestry-specific variant sets

**Key Functions:**
- Hail MatrixTable operations
- Plink export

---

### Z.003: Prepare Pruned BED
LD-prunes genotypes to reduce multicollinearity in PCA.

**Files:**
- `Z.003_prepare_pruned_bed.py` / `.ipynb`

**Workflow:**
1. Read Plink files from Z.001
2. Apply LD-pruning (removes correlated variants)
3. Export pruned genotypes

**Outputs:**
- Pruned `.bed`, `.bim`, `.fam` files for PCA

**Why Pruning?**
- Reduces computational cost of PCA
- Removes redundant information (variants in high LD)
- Improves numerical stability of eigendecomposition

---

### Z.005: Merge & Run PCA
Merges ancestry-specific genotypes and performs principal component analysis.

**Files:**
- `Z.005_merge_bed_run_pca.py` / `.ipynb`

**Workflow:**
1. Merge pruned genotypes across ancestry groups
2. Calculate ancestry-specific principal components
3. Extract top 20 PCs (used as covariates in GWAS)

**Outputs:**
- PCA scores (person × PC matrix)
- Loadings and variance explained
- Ancestry-specific PC sets

**Machine Requirements:**
- Sufficient memory for full genotype matrix
- Typical: n2-standard-32 or larger

**Tools Used:**
- Plink (`--pca` command)
- Python pandas/numpy for post-processing

---

### Z.007: Assess & Format PCA
Quality control and formatting of PCA results for downstream use.

**Files:**
- `Z.007_assess_and_format_pca.py` / `.ipynb`
- `Z.008_pca_gt_correlation (R).py` / `.ipynb` - Optional R visualization

**QC Checks:**
- Variance explained by each PC
- Scree plot (cumulative variance)
- PC correlation with known ancestry labels
- Outlier detection

**Outputs:**
- Formatted PC file for GWAS null models
- Summary statistics and plots
- R script for advanced visualization (optional)

---

### Z.009: Prepare Ancestry Genotypes & VCF
Exports full genotypes in formats suitable for downstream analyses.

**Files:**
- `Z.009_prepare_ancestry_genotypes.py` / `.ipynb` - Ancestry-specific genotypes
- `Z.009_prepare_vcf.py` / `.ipynb` - VCF format export

**Workflow:**
1. Load original All of Us genomic data
2. Filter to specific ancestry group(s)
3. Export to analysis formats:
   - Plink (`.bed`, `.bim`, `.fam`) - for GWAS testing
   - VCF - for other genomics tools

**Outputs:**
- Per-ancestry Plink files
- VCF files (chromosome-wise or genome-wide)

**Uses:**
- SAIGE GWAS (requires ancestry-specific genotypes)
- Population genetics analysis
- Integration with external tools

---

### 001: Cohort Definition (Hail)
General cohort extraction and characterization using Hail framework.

**Files:**
- `001 Cohort (Hail).py`

**Features:**
- Loads All of Us participant metadata
- Extracts persons with whole genome sequencing available
- Merges with OMOP clinical data
- Prepares for genomic analyses

---

### G.11: GPR15 Cohort - GI Inflammation
Case study: Inflammatory bowel disease (IBD) cohort focusing on GPR15 locus.

**Files:**
- `G.11 GPR15 Cohort - GI inflammation.py` / `.ipynb`

**Purpose:**
- Demonstrates detailed cohort definition workflow
- Focuses on gastrointestinal inflammation phenotype
- Shows integration of genomic + clinical data

---

## Data Flow to Downstream Projects

This genomics project produces outputs used by:

### HPV GWAS (`hpv/`)
- Uses ancestry-specific PCs from Z.007
- Uses pruned genotypes from Z.003 for null model
- Uses full genotypes from Z.009 for GWAS testing

### Sarcoid GWAS (`sarcoid/`)
- Same data flow as HPV project
- Includes ancestry-specific stratification

### Other Analyses
- Custom PCA for new phenotypes
- Population structure characterization
- Ancestry classification for non-GWAS projects

---

## Key Concepts

### Ancestry-Informative Variants
Genetic variants that differ in frequency across continental ancestry groups. Used for:
- Inferring ancestry of study participants
- Ancestry stratification in GWAS
- Population genetics analysis

### LD-Pruning
Removes variants in high linkage disequilibrium (LD > threshold). Benefits:
- Reduces redundancy (highly correlated variants excluded)
- Faster computation
- More stable statistical results
- Common threshold: r² > 0.1 at 500kb windows

### Principal Component Analysis (PCA)
Dimensionality reduction technique that:
- Identifies major axes of genetic variation
- Maps to continental ancestry (PC1, PC2, etc.)
- Generates ancestry-specific covariates for GWAS
- Standard: use top 20 PCs in GWAS null models

---

## Technical Details

### Hail Framework
Scalable genomics analysis platform built on Apache Spark.

**Key Data Structure: MatrixTable**
- Genomic data organized as matrix (variants × samples)
- Row annotations: variant metadata (position, allele, consequences)
- Column annotations: sample metadata (ancestry, phenotype)
- Entry data: genotypes (0/0, 0/1, 1/1)

**Common Operations:**
```python
mt = hl.read_matrix_table("path/to/data.mt")
mt = mt.filter_rows(...)  # Filter variants
mt = mt.filter_cols(...)  # Filter samples
hl.export_plink(mt, "output_prefix")  # Export to Plink format
```

### Plink Format
Standard format for genetic data:
- `.bed` - Binary genotype data
- `.bim` - Variant information (position, alleles, etc.)
- `.fam` - Sample information (family, phenotype, covariates)

### VCF Format
Variant Call Format: standardized text format for genomic variants. Advantages:
- Human-readable
- Standard across genomics tools
- Includes variant annotations
- Larger file size than Plink

---

## Machine Type Recommendations

| Task | Machine Type | Cost | Notes |
|------|-------------|------|-------|
| PCA (sample size < 100k) | n2-standard-32 | $$ | 32 vCPU, 128GB RAM |
| PCA (sample size > 100k) | n2-standard-64 | $$$ | 64 vCPU, 256GB RAM |
| Plink operations | c2-standard-8 | $ | 8 vCPU, 32GB RAM (usually sufficient) |
| Hail processing | n2-highmem-16+ | $$ | High memory for large MatrixTables |

---

## File Naming Convention

- `Z.###` - Infrastructure/utility analyses
- `G.##` - Gene-specific analyses
- `001` - Cohort definitions
- Numerical prefix indicates analysis order

---

## Getting Help

- **All of Us Support**: [Researcher Help Desk](https://support.researchallofus.org/hc/en-us/)
- **Hail Documentation**: https://hail.is/docs/0.2/
- **GWAS Projects**: See `hpv/CLAUDE.md` or `sarcoid/CLAUDE.md` for downstream usage
- **Main Repository**: See main `CLAUDE.md` for overall toolkit guidance

---

**Last Updated**: January 2026
**Status**: Active - Used by downstream GWAS analyses
