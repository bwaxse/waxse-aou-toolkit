# HPV Genetic Risk Association Study

**Research Question**: Identify genetic variants associated with HPV-related anal and vulvovaginal/cervical disease in the All of Us cohort

**Platform**: All of Us Research Workbench (OMOP CDM v8)
**Analysis Type**: GWAS (Genome-Wide Association Study) with multi-ancestry meta-analysis
**Primary Methods**: SAIGE (binary traits), METAL (meta-analysis), HLA typing analysis

---

## Getting Started

**BEFORE running any notebooks in this project:**

1. Run `_reference/verily/00_setup_workspace.ipynb` once per JupyterLab session
   - This initializes `WORKSPACE_CDR`, `WORKSPACE_BUCKET`, etc.
   - Required for all subsequent notebooks to work
   - See main `CLAUDE.md` for full setup details

2. Run notebooks in order: `B01 → B02 → B03 → B04 → B05 → B06 → B07`
   - Each stage depends on outputs from previous stage
   - See "Project Workflow Overview" below for detailed pipeline description

---

## Project Workflow Overview

This project follows a systematic 7-stage pipeline (B01-B07) for conducting GWAS on HPV-associated conditions:

```
B01: Cohort Definition → B02: SAIGE GWAS → B03: Meta-analysis →
B04: Post-processing → B05: Download → B06: HLA Visualization → B07: HLA QC
```

Each stage produces outputs used by subsequent stages. Files exist in both `.py` and `.ipynb` formats with identical code.

---

## File Organization

**Source files (tracked in git):**
- `.py` files in project root - Read these for understanding code and development

**Notebooks (generated on demand):**
- `.ipynb` files in `notebooks/` subdirectory - Generated for Verily Workbench uploads
- Request: "convert B01.a HPV Cohort v0.py to notebook"
- Both tracked in git for collaboration, masked from Claude context to save tokens

---

## File Map by Analysis Stage

### B01: Cohort Definition (3 cohorts)
Extracts HPV case/control cohorts from All of Us OMOP data and prepares for GWAS.

#### Files:
- **B01.a HPV Cohort v0.py** - General HPV cohort (all anatomic sites combined)
  - Extracts persons with HPV-related diagnoses from OMOP `condition_occurrence`
  - Applies inclusion criteria: whole genome sequencing available, specific HPV concept IDs
  - Adds demographics (race, ethnicity, sex, age) and socioeconomic covariates (employment, income, education, insurance)
  - Creates Table 1 descriptive statistics
  - Runs PheWAS analysis using PheTK
  - **Key Functions**: `polars_gbq()`, `update_covariates()`, `describe_group()`
  - **Inputs**: OMOP tables (person, condition_occurrence, observation)
  - **Outputs**: `hpv_gwas_cohort.csv`, `hpv_phewas_cohort.csv`, `hpv_phewas_results.csv`

- **B01.b Anal HPV Cohort v1 (5.20).py** - Anal HPV-specific cohort
  - Focuses on anal canal HPV diagnoses (ICD-10: D01.3, etc.)
  - Stratifies by HIV status for subgroup analysis
  - **Outputs**: `hpv_anal_gwas_cohort.csv`

- **B01.b Vulvovaginal HPV Cohort v2 (12.12).py** - Vulvovaginal/cervical HPV cohort
  - Female-only cohort with cervical/vulvar/vaginal HPV
  - **Outputs**: `hpv_vv_gwas_cohort.csv`

#### Key Patterns:
- BigQuery SQL with OMOP CDM structure
- Covariate engineering: binary indicators for demographics/SES categories
- Hierarchical concept ID searching (ICD9CM, ICD10CM, SNOMED)

---

### B02: SAIGE GWAS Analysis (4 files)
Runs genome-wide association tests using SAIGE across multiple ancestries.

#### Files:
- **B02.a Anal HPV SAIGE GWAS v0.py** - Initial anal HPV GWAS (all ancestries)
  - Sets up SAIGE null model and hypothesis testing for 5 ancestries
  - **Key Functions**: `dsub_script()`, `run_saige_null()`, `run_saige_test()`, `check_dsub_status()`
  - **Null Model Covariates**: `imputed_sex`, `age_at_last_ehr`, `ancPC1-20`
  - **Machine Type**: c4-standard-8 (8 vCPUs)
  - **Containerized**: Uses `wzhou88/saige:1.3.6` Docker image

- **B02.b HIV- Anal HPV SAIGE GWAS v1.py** - HIV-negative stratified analysis
  - Subset analysis excluding HIV+ individuals
  - Same SAIGE pipeline as v0

- **B02.b HIV- Anal HPV SAIGE GWAS v1 redo.py** - Redo/validation run

- **B02.c HIV- Vulvovaginal HPV SAIGE GWAS v2.py** - Female-specific GWAS
  - Vulvovaginal/cervical HPV analysis

#### SAIGE Workflow:
```bash
# Step 1: Fit null GLMM model using pruned genotypes
step1_fitNULLGLMM.R --plinkFile=pruned_genotypes \
    --phenoFile=gwas_metadata.tsv \
    --phenoCol=condition__hpv \
    --covarColList=imputed_sex,age_at_last_ehr,ancPC1-20 \
    --traitType=binary \
    --outputPrefix=saige_null_model

# Step 2: Run genome-wide tests per chromosome (parallelized via dsub)
step2_SPAtests.R --vcfFile=chr${CHR}.vcf.gz \
    --GMMATmodelFile=saige_null_model.rda \
    --varianceRatioFile=saige_null_model.varianceRatio.txt \
    --minMAC=20 \
    --SAIGEOutputFile=gwas_results_chr${CHR}.txt
```

#### Key Components:
- **dsub submission helper**: Manages distributed computing on Google Cloud Batch
- **Reference data paths**: Points to pre-computed pruned genotypes, imputed VCFs
- **Ancestry-specific processing**: EUR, AFR, AMR, EAS, SAS analyzed separately
- **Outputs**: Per-chromosome GWAS results in `swarm_gwas/` directories

---

### B03: Meta-Analysis (2 files)
Combines ancestry-specific GWAS results using METAL.

#### Files:
- **B03 METAL Meta-analysis (EUR AFR AMR).py**
  - Merges chromosome-wise results into single file per ancestry
  - Validates sample sizes and file completeness
  - Runs METAL inverse-variance weighted meta-analysis
  - **Key Functions**: `combine_gwas_results()`, `validate_saige_inputs()`, `run_metal()`
  - **METAL Configuration**:
    - Scheme: STDERR (inverse-variance weighting)
    - Weight: `n_total` (total sample size per ancestry)
    - Auto-flip alleles based on frequency
  - **Machine Type**: n2d-standard-8
  - **Container**: Custom `bwaxse/metal` image

- **B03 Summarize GWAS.py**
  - Quality control and summarization of GWAS results
  - Checks for expected files, missing chromosomes
  - Generates summary statistics tables

#### METAL Script Template:
```metal
MARKER vid
WEIGHT n_total
ALLELE effect_allele other_allele
FREQ effect_allele_frequency
PVAL p_value
EFFECT beta
STDERR standard_error

SCHEME STDERR
AVERAGEFREQ ON
MINMAXFREQ ON

PROCESS {ancestry1_file}
PROCESS {ancestry2_file}
PROCESS {ancestry3_file}

ANALYZE HETEROGENEITY
```

---

### B04: Post-Processing (2 files)
Cleans and formats GWAS results for downstream analysis.

#### Files:
- **B04.a Anal GWAS Post-Processing.py** - Post-process anal HPV results
- **B04.b Vulvovaginal GWAS Post-Processing.py** - Post-process vulvovaginal results

#### Typical Operations:
- Filter by MAF, INFO score, missingness thresholds
- Annotate variants with rsIDs, genes, functional annotations
- Calculate lambda_GC (genomic inflation)
- Identify genome-wide significant hits (p < 5×10⁻⁸)

---

### B05: Data Download (1 file)
Downloads final summary statistics for sharing/publication.

- **B05 Download Sumstats.py**
  - Exports cleaned GWAS results
  - Formats for standard repositories (GWAS Catalog, dbGaP)
  - Creates README/metadata files

---

### B06: HLA Visualization (2 files)
Specialized visualization for HLA region associations.

#### Files:
- **B06.1 Plot HLA Results 4-Digit.py** - 4-digit HLA resolution plots
  - Creates Manhattan plots overlaying HLA alleles on genomic backbone
  - Forest plots comparing effect sizes across ancestries
  - **Key Functions**: `create_hla_manhattan_plot()`, `create_hla_forest_plot()`
  - **Visualization Library**: gwaslab, matplotlib, seaborn
  - **Data Sources**:
    - HLA typing results from `HLA/four_digit_HPV_{ancestry}_HLA_results_freq.csv`
    - GWAS summary statistics from B02/B03
  - **Plot Features**:
    - Gene track with HLA-A, HLA-B, HLA-C, HLA-DR, HLA-DQ, HLA-DP loci
    - Class I (red) vs Class II (blue) color coding
    - Bonferroni significance thresholds
    - Multi-panel comparison across ancestries

- **B06.1 Plot HLA Results 6-Digit.py** - 6-digit HLA resolution (higher resolution)
  - Same structure as 4-digit but with more specific alleles

#### HLA Analysis Workflow:
1. Load HLA typing results and GWAS summary stats
2. Extract HLA region variants (chr6:29,602,238-33,409,896)
3. Map HLA alleles to genomic positions using gene annotations
4. Create visualizations:
   - **Manhattan**: Shows -log10(p) vs genomic position
   - **Forest**: Shows odds ratios with 95% CI across ancestries
5. Identify protective vs risk alleles

---

### B07: HLA Quality Control (1 file)
Quality checks for HLA imputation and association results.

- **B07 HLA DMA DRA Check.py**
  - Checks allele counts for HLA-DMA and HLA-DRA genes
  - These genes often have AC=1 (monomorphic) and need to be excluded
  - Validates HLA imputation quality metrics

---

## Key Reusable Components

### 1. BigQuery Helper Functions
```python
def polars_gbq(query: str) -> pl.DataFrame:
    """Execute BigQuery SQL and return polars DataFrame"""
    client = bigquery.Client()
    query_job = client.query(query)
    rows = query_job.result()
    return pl.from_arrow(rows.to_arrow())
```

### 2. dsub Distributed Computing Helper
```python
def dsub_script(machine_type, envs, in_params, out_params, script, ...):
    """Submit dsub job to Google Cloud Batch"""
    # Constructs dsub command with:
    # - Machine specs (CPU, memory, disk)
    # - Environment variables
    # - Input/output file mappings
    # - Docker container image
    # - Preemptible instances for cost savings
```

### 3. SAIGE Pipeline Functions
```python
def run_saige_null(anc, trait, trait_type, covs, covs_discrete, script):
    """Fit SAIGE null GLMM model"""

def run_saige_test(anc, trait, script, chroms=range(1,23)):
    """Run SAIGE hypothesis tests per chromosome"""
```

### 4. METAL Meta-Analysis Functions
```python
def combine_gwas_results(my_bucket, in_dir, trait_type, output_path):
    """Merge chromosome-wise GWAS results into single file"""

def validate_saige_inputs(trait, ancestries, base_output_folder):
    """Check SAIGE outputs exist and extract sample sizes"""

def run_metal(trait, ancestries, base_output_folder):
    """Execute METAL meta-analysis across ancestries"""
```

### 5. Covariate Engineering
```python
def update_covariates(df, sarcoid=True):
    """Create binary indicator variables for demographics/SES"""
    # Age bins: 17-24, 25-34, 35-44, 45-54, 55-64, 65-74, 75+
    # Race: black, white, asian, hispanic, multip, mena, aian, nhpi
    # SES: employment, income, education, insurance
```

### 6. Table 1 Generation
```python
def describe_group(df, sarcoid=True):
    """Generate descriptive statistics table"""
    # Median (IQR) for continuous variables
    # n (%) for categorical variables
    # Suppresses counts < 20 per All of Us policy
```

---

## Common Patterns and Architecture

### Data Flow
```
OMOP BigQuery Tables
  ↓
Cohort Definition (polars/pandas)
  ↓
GWAS Metadata TSV (person_id, phenotype, covariates)
  ↓
SAIGE (distributed via dsub)
  ↓
Per-chromosome results
  ↓
METAL meta-analysis
  ↓
Post-processing & Visualization
```

### File Naming Conventions
- **B##** prefix indicates analysis stage sequence
- **v#** suffix indicates iteration/version
- Dates in parentheses indicate last major update (5.20 = May 20)
- Ancestry codes: EUR (European), AFR (African), AMR (Admixed American), EAS (East Asian), SAS (South Asian)

### Storage Structure
```
{workspace_bucket}/
├── saige_gwas/{version}/{ancestry}/{trait}/
│   ├── gwas_metadata.tsv           # Phenotype + covariates
│   ├── saige_null_model.rda        # Fitted null model
│   ├── saige_null_model.varianceRatio.txt
│   └── swarm_gwas/
│       └── gwas_results_chr{1-22}.txt  # Per-chromosome results
├── HLA/
│   ├── four_digit_HPV_{ancestry}_HLA_results_freq.csv
│   └── meta_HPV_4D_1.txt           # Meta-analysis HLA results
└── dsub/logs/                      # Job execution logs
```

### Ancestry-Specific Processing
All GWAS analyses stratify by genetic ancestry:
- **EUR**: European ancestry (largest N)
- **AFR**: African ancestry
- **AMR**: Admixed American (Hispanic/Latino)
- **EAS**: East Asian ancestry (smaller N, sometimes excluded)
- **SAS**: South Asian ancestry (smaller N, sometimes excluded)

Ancestry assignment uses predicted ancestry from `cohort_metadata__pcs__phenotypes.tsv` (column: `ancestry_pred_other`).

### Covariates Standard Set
```python
covariates = ['imputed_sex', 'age_at_last_ehr'] + ['ancPC{}'.format(str(x)) for x in range(1, 21)]
# imputed_sex: M/F from genotype data
# age_at_last_ehr: Age at end of study period
# ancPC1-PC20: Genetic principal components (within-ancestry)
```

### Statistical Thresholds
- **Genome-wide significance**: p < 5×10⁻⁸ (red line on Manhattan plots)
- **Bonferroni correction**: p < 0.05 / N_tests (purple line, varies by analysis)
- **Minimum MAC**: 20 (minimum minor allele count for inclusion)

### Docker Images Used
- **SAIGE**: `wzhou88/saige:1.3.6` - GWAS analysis
- **METAL**: Custom `bwaxse/metal` image - Meta-analysis
- **Base environment**: `us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.2.14`

### Cloud Computing Patterns
- **Null model fitting**: Single job, c4-standard-8, non-preemptible
- **Chromosome-wise tests**: 22 parallel jobs, c4-standard-8, preemptible
- **Meta-analysis**: Single job, n2d-standard-8, preemptible
- **Storage**: Google Cloud Storage buckets
- **Orchestration**: dsub with Google Cloud Batch backend

---

## Dependencies

### Python Libraries
```python
# Data manipulation
pandas, polars, numpy, scipy, pyarrow

# Cloud integration
google-cloud-bigquery, google-cloud-storage

# Visualization
matplotlib, seaborn, gwaslab

# Specialized analysis
PheTK  # PheWAS analysis toolkit for All of Us

# Utilities
jinja2  # Template rendering
subprocess  # Shell command execution
```

### External Tools (via Docker)
- **SAIGE** (1.3.6): GWAS for binary/quantitative traits with mixed models
- **METAL**: Inverse-variance weighted meta-analysis
- **PLINK2**: Genotype file manipulation

### Reference Data Sources
Located in `{src_bucket}/data/`:
- `cohort_metadata__pcs__phenotypes.tsv.gz` - Pre-computed demographics, PCs, example phenotypes
- `stg105/{ancestry}/targimp_genotypes_chr{1-22}.vcf.gz` - Imputed genotypes (per chromosome)
- `stg105/{ancestry}/reference_samples.txt` - Samples to exclude from analysis
- `stg303/{ancestry}/pruned_genotypes.{bed,bim,fam}` - LD-pruned variants for null model
- Ensembl GTF (hg38) for gene annotations

---

## Study-Specific Notes

### HPV Concept IDs (OMOP)
The cohort definition uses multiple OMOP concept IDs to capture HPV diagnoses:
- **45572221**: HPV-related concept
- **1567478**: Hierarchical concept for HPV
- **35206462, 1572289**: Additional HPV-related ICD codes
- **44826393, 44826394, 44829052, 44820669**: Site-specific HPV codes
- **44828175, 44835099, 44835101, 44835098, 44835100, 44831613**: Cervical/anal/vulvar HPV
- **35205611, 44828697**: Additional anatomic site codes

Exclusion criteria include specific conditions (concept IDs: 44829737, 35205776, 44824945, 35225089).

### HIV Stratification
Some analyses stratify by HIV status:
- HIV+ individuals have different HPV risk profiles
- Separate GWAS conducted for HIV- subgroup
- Allows identification of HIV-specific genetic modifiers

### HLA Region Special Handling
- **HLA-DMA and HLA-DRA** are excluded (monomorphic, AC=1)
- HLA region: chr6:29,602,238-33,409,896 (hg38)
- Separate imputation and association testing for HLA alleles (4-digit and 6-digit resolution)
- Visualization uses special gene tracks and color coding (Class I = red, Class II = blue)

---

## Usage Examples

### Run Complete GWAS Pipeline
```python
# 1. Define cohort
# Run B01.a HPV Cohort v0.py to generate hpv_gwas_cohort.csv

# 2. Prepare GWAS metadata
# Run B02.a to merge cohort with genomic PCs and run SAIGE

# 3. Meta-analyze across ancestries
# Run B03 METAL Meta-analysis

# 4. Visualize HLA region
# Run B06.1 Plot HLA Results 4-Digit.py
```

### Check Job Status
```python
check_dsub_status(user='bwaxse', full=False, age='3d')
```

### Access Results
```python
# Load GWAS results
eur_results = pl.read_csv(
    f'{my_bucket}/saige_gwas/v1_redo/eur/condition__hpv/gwas_results.tsv.gz',
    separator='\t'
)

# Load meta-analysis results
metal_results = pd.read_csv(
    f'{my_bucket}/saige_gwas/v1_redo/metal/condition__hpv/condition__hpv_afr_amr_eur1.tsv',
    sep='\t'
)
```

---

## Tips for Claude Code Reuse

1. **For new binary phenotype GWAS**: Adapt B01 cohort definition with new concept IDs, then run B02-B03 pipeline unchanged

2. **For quantitative trait GWAS**: Modify `trait_type='binary'` to `'quantitative'` in B02, adjust covariates as needed

3. **For single-ancestry analysis**: Skip B03 METAL step, proceed directly to B04 post-processing

4. **For HLA-focused study**: Use B06 visualization templates, swap in new HLA typing results

5. **Common function extraction opportunities**:
   - BigQuery cohort extraction pattern (B01)
   - SAIGE null model + test workflow (B02)
   - METAL meta-analysis wrapper (B03)
   - HLA visualization functions (B06)

6. **Reusable bash scripts**:
   - `run_saige_null_model.sh`
   - `run_saige_chrom.sh`
   - `run_metal.sh`

---

## Limitations and Known Issues

1. **Small sample sizes in EAS/SAS**: Often insufficient power, excluded from meta-analysis
2. **HLA-DMA/DRA monomorphic**: Must be removed before analysis
3. **Dsub c4 machine types**: Don't work with google-batch provider (use c2 or n2)
4. **OMOP concept coverage**: May miss HPV diagnoses coded with non-standard vocabularies
5. **Ancestry prediction accuracy**: Admixed individuals may be misclassified
6. **Missing covariate handling**: SAIGE handles missing data but may reduce effective sample size

---

## Metadata
- **Author**: Bennett Waxse
- **Platform**: All of Us Research Workbench
- **Data Version**: CDR v8
- **Last Updated**: December 2025
- **License**: MIT
