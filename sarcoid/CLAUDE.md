# Sarcoidosis Genetic Risk Association Study

**Research Question**: Identify genetic variants associated with sarcoidosis in the All of Us cohort

**Platform**: All of Us Research Workbench (OMOP CDM v8)
**Analysis Type**: GWAS (Genome-Wide Association Study) with multi-ancestry meta-analysis
**Primary Methods**: SAIGE (binary trait), METAL (meta-analysis), dual case definition sensitivity analysis

---

## Getting Started

**BEFORE running any notebooks in this project:**

1. Run `_reference/verily/00_setup_workspace.ipynb` once per JupyterLab session
   - This initializes `WORKSPACE_CDR`, `WORKSPACE_BUCKET`, etc.
   - Required for all subsequent notebooks to work
   - See main `CLAUDE.md` for full setup details

2. Run notebooks in order: `01 → 02 → 03 → 04 → 05`
   - Each stage depends on outputs from previous stage
   - See "Project Workflow Overview" below for detailed pipeline description

---

## Project Workflow Overview

This project follows a systematic 5-stage pipeline (01-05) for conducting GWAS on sarcoidosis:

```
01: Cohort Definition → 02: SAIGE GWAS → 03: Meta-analysis →
04: Post-processing & Visualization → 05: Download
```

**Unique Feature**: The pipeline includes TWO parallel analyses with different case definitions for sensitivity analysis:
- **≥ 1 code**: More inclusive (captures all diagnosed cases)
- **≥ 2 codes**: More stringent (requires confirmation on separate dates)

Each stage produces outputs used by subsequent stages. Files exist in both `.py` and `.ipynb` formats with identical code.

---

## File Organization

**Source files (tracked in git):**
- `.py` files in project root - Read these for understanding code and development

**Notebooks (generated on demand):**
- `.ipynb` files in `notebooks/` subdirectory - Generated for Verily Workbench uploads
- Request: "convert 01 Sarcoid Cohort.py to notebook"
- Both tracked in git for collaboration, masked from Claude context to save tokens

---

## File Map by Analysis Stage

### 01: Cohort Definition (1 file)
Extracts sarcoidosis case/control cohorts from All of Us OMOP data with sophisticated exclusion logic.

#### File:
- **01 Sarcoid Cohort.py** - Comprehensive cohort extraction and characterization
  - Uses `AouQueries` class for flexible OMOP concept searching
  - Searches for "sarcoid" across ICD9CM, ICD10CM, SNOMED vocabularies
  - Extracts from both `condition_occurrence` and `observation` tables
  - Implements special V-code handling for ICD9CM codes
  - Applies exclusion criteria to remove competing granulomatous diagnoses from controls
  - Adds demographics, procedures, and covariates
  - Creates descriptive statistics and visualizations
  - **Key Classes**: `AouQueries`
  - **Key Methods**:
    - `find_diagnosis_codes()` - Search OMOP concept tables with flexible text/code matching
    - `person_code_df()` - Create person-level diagnosis dataframe
    - `find_procedure_codes()` - Search procedure tables
    - `person_procedure_df()` - Create person-level procedure dataframe
  - **Inputs**: OMOP tables (condition_occurrence, observation, person, procedure_occurrence), genomic_metrics.tsv, ancestry_specific_pcs.parquet
  - **Outputs**:
    - `all_sarcoid_cases.csv` - Cases with exclusion diagnoses removed from controls
    - `sarcoid_cases_without_exclusions.csv` - All sarcoid cases before applying exclusions

#### Cohort Definition Strategy:
- **Inclusion**: Text search for "sarcoid" across standard vocabularies
- **Exclusion Logic** (applied to controls only):
  - Competing granulomatous disease diagnoses
  - Exclusion drug usage within 1 year of sarcoid diagnosis
  - Creates `has_exclusion_condition` and `exclusion_drug` flags
- **Case Counting**: Tracks `sarcoid_n` (number of sarcoid diagnosis codes per person)

#### Key Patterns:
- **V-code Special Handling**: ICD9CM V codes require concept_relationship table joins for correct vocabulary attribution
- **Dual Path Extraction**: Searches both `source_value` and `source_concept_id` paths to maximize sensitivity
- **Flexible Querying**: Jinja2 templating for dynamic SQL generation based on search criteria

---

### 02: SAIGE GWAS Analysis (2 files)
Runs genome-wide association tests using SAIGE with two different case definitions.

#### Files:
- **02 SAIGE GWAS (≥ 1 code).py** - Analysis using all cases with at least 1 diagnosis
  - Uses `all_sarcoid_cases.csv` cohort
  - Higher sample size, more inclusive
  - **Output folder**: `saige_gwas/min_1/`

- **02 SAIGE GWAS (≥ 2 codes).py** - More stringent analysis requiring 2+ diagnoses on separate dates
  - Filters to cases with `sarcoid_n >= 2`
  - Lower sample size, higher diagnostic confidence
  - Reduces false positives from coding errors or rule-out diagnoses
  - **Output folder**: `saige_gwas/min_2/`

#### Shared Configuration:
- **Ancestries**: EUR (European), AFR (African), plus transancestry analysis
- **Covariates**: `sex_binary`, `end_of_study_age`, `ancPC1-20` (ancestry-specific PCs)
- **Discrete Covariates**: `sex_binary` (for sex chromosome handling)
- **Trait Type**: Binary
- **SAIGE Parameters**:
  - `minMAC=20` - Minimum minor allele count
  - `LOCO=TRUE` - Leave one chromosome out (reduces proximal contamination)
  - Firth correction for rare variants
- **Machine Type**: c4-standard-8 (8 vCPUs)
- **Container**: `wzhou88/saige:1.3.6` Docker image

#### Key Functions:
```python
def dsub_script() # Submit distributed computing jobs to Google Cloud Batch
def run_saige_null() # Fit null GLMM model using pruned genotypes
def run_saige_step2() # Run hypothesis tests per chromosome (parallelized)
def check_dsub_status() # Monitor job execution status
def parse_pcs() # Extract and format principal components
def binarize_sex_chromosomes() # Create binary sex variable (XX=0, XY=1)
def print_ancestry_summary() # Display case/control counts by ancestry
```

#### SAIGE Workflow:
```bash
# Step 1: Fit null GLMM model
step1_fitNULLGLMM.R --plinkFile={anc}_merged_genotypes \
    --phenoFile=gwas_metadata.tsv \
    --phenoCol=condition__sarcoid \
    --covarColList=sex_binary,end_of_study_age,ancPC1-20 \
    --qCovarColList=sex_binary \
    --traitType=binary \
    --outputPrefix=saige_null_model

# Step 2: Run genome-wide tests (distributed, 1 job per chromosome)
step2_SPAtests.R --bedFile=genotypes_chr${CHR}.bed \
    --GMMATmodelFile=saige_null_model.rda \
    --varianceRatioFile=saige_null_model.varianceRatio.txt \
    --minMAC=20 \
    --LOCO=TRUE \
    --is_Firth_beta=TRUE \
    --SAIGEOutputFile=gwas_results_chr${CHR}.txt
```

#### Key Components:
- **Sex handling**: Binary encoding (0=female, 1=male) from dragen_sex_ploidy
- **Ancestry assignment**: Uses pre-computed ancestry predictions and ancestry-specific PCs
- **Reference data**:
  - Pruned genotypes: `{anc}_merged_genotypes.{bed,bim,fam}` for null model
  - Full genotypes: `genotypes_chr{1-22}.{bed,bim,fam}` for testing
- **Outputs**: Per-chromosome GWAS results in `{output_folder}/{anc}/condition__sarcoid/swarm_gwas/`

---

### 03: Meta-Analysis (1 file)
Combines ancestry-specific GWAS results using METAL inverse-variance weighting.

#### File:
- **03 METAL Meta-analysis (EUR AFR).py**
  - Merges chromosome-wise results into single file per ancestry
  - Validates sample sizes and file completeness
  - Runs METAL inverse-variance weighted meta-analysis for EUR + AFR
  - **Key Functions**: `combine_gwas_results()`, `validate_saige_inputs()`, `extract_sample_sizes()`, `run_metal()`
  - **METAL Configuration**:
    - Scheme: STDERR (inverse-variance weighting)
    - Weight: `n_total` (total sample size per ancestry)
    - Auto-flip alleles based on frequency
    - Heterogeneity testing enabled
  - **Machine Type**: n2d-standard-8
  - **Container**: Custom `bwaxse/metal` image

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

PROCESS {eur_file}
PROCESS {afr_file}

ANALYZE HETEROGENEITY
```

#### Outputs:
- Meta-analysis results: `condition__sarcoid_afr_eur1.tsv`
- Heterogeneity info: `condition__sarcoid_afr_eur1.tsv.info`

---

### 04: Post-Processing (2 files)
Quality control, lead SNP identification, and visualization.

#### Files:
- **04 GWAS Post-Processing.py** - Comprehensive QC and visualization using gwaslab
  - Format conversion for gwaslab compatibility
  - Manhattan and QQ plot generation
  - Multi-ancestry stacked plots
  - Regional association plots for lead SNPs
  - Lead SNP identification with LD-based clumping (500kb windows)
  - Annotation with GWAS Catalog known loci
  - **Key Functions**: `convert_snp_format()`, `add_chr_pos_from_snpid()`, `annotate_with_known_loci()`, `process_gwas()`

- **04 Summarize GWAS.py** - Summary statistics and quality checks
  - Validates SAIGE output completeness
  - Checks for missing chromosomes
  - Generates summary tables by ancestry
  - Lambda GC calculation (genomic inflation)
  - Counts genome-wide significant hits (p < 5×10⁻⁸)

#### Typical Operations:
- Filter by MAF, INFO score thresholds
- Extract lead SNPs (p < 5×10⁻⁸)
- LD-based clumping to identify independent signals
- Annotate with rsIDs, nearest genes, known GWAS loci
- Generate publication-quality plots

---

### 05: Data Download (1 file)
Prepares summary statistics for sharing and publication.

- **05 Download Sumstats.py**
  - Splits large GWAS result files into <100MB chunks for download
  - Creates chunked files: `{ancestry}_gwas_results_chunk_{001-NNN}.tsv.gz`
  - Exports cleaned GWAS results
  - Prepares metadata files for repositories (GWAS Catalog, dbGaP)
  - **Key Function**: `split_gwas_file()`

---

## Key Reusable Components

### 1. AouQueries Class - Flexible OMOP Searching
```python
class AouQueries:
    def __init__(version, bucket):
        """Initialize with All of Us CDR version and workspace bucket"""

    def find_diagnosis_codes(
        vocabularies=['ICD9CM', 'ICD10CM', 'SNOMED'],
        search_terms=['sarcoid'],
        exclude_terms=None,
        exact_codes={'ICD10CM': ['D86.0', 'D86.9']},
        pattern_codes={'ICD10CM': ['D86.%']},
        exclude_codes={'ICD10CM': ['D86.1']},
        person_ids=None
    ):
        """
        Flexible diagnosis code search with multiple matching strategies:
        - search_terms: Text search in concept_name (e.g., 'sarcoid')
        - exact_codes: Exact code matching by vocabulary
        - pattern_codes: SQL LIKE pattern matching (e.g., 'D86.%')
        - exclude_codes: Remove specific codes
        - exclude_terms: Remove concepts with certain text

        Special handling for ICD9CM V codes (requires concept_relationship mapping)
        Searches both source_value and source_concept_id paths
        """

    def person_code_df(count_threshold=1):
        """
        Create person-level summary with:
        - case: 0/1 binary indicator (≥ threshold codes)
        - {condition}_n: Total diagnosis count per person
        - earliest_date, latest_date: First and last diagnosis dates
        """

    def find_procedure_codes(
        vocabularies=['CPT4', 'HCPCS'],
        search_terms=None,
        exact_codes=None,
        pattern_codes=None
    ):
        """Search procedure_occurrence table with similar flexibility"""

    def person_procedure_df():
        """Create person-level procedure summary"""
```

### 2. SAIGE Pipeline Functions
```python
def run_saige_null(anc, trait, trait_type, covs, covs_discrete, script):
    """
    Fit SAIGE null GLMM model

    Args:
        anc: Ancestry code (eur, afr, amr, etc.)
        trait: Phenotype name (e.g., 'condition__sarcoid')
        trait_type: 'binary' or 'quantitative'
        covs: List of covariate column names
        covs_discrete: List of categorical covariates (for sex chromosome handling)
        script: dsub submission script
    """

def run_saige_step2(anc, trait, script, chroms=range(1,23)):
    """
    Run SAIGE hypothesis tests per chromosome (parallelized)

    Submits 22 separate dsub jobs (one per chromosome)
    Uses preemptible instances for cost savings
    """

def check_dsub_status(user='bwaxse', full=False, age='3d'):
    """Monitor dsub job execution status"""

def parse_pcs(df, pc_prefix='ancPC'):
    """
    Extract principal components from wide-format dataframe

    Returns: List of PC column names (ancPC1, ancPC2, ..., ancPC20)
    """

def binarize_sex_chromosomes(df):
    """
    Create binary sex variable from dragen_sex_ploidy

    Encoding:
    - XX → 0 (female)
    - XY → 1 (male)
    - Other → excluded
    """

def print_ancestry_summary(df, trait):
    """Display case/control counts stratified by ancestry"""
```

### 3. METAL Meta-Analysis Functions
```python
def combine_gwas_results(my_bucket, in_dir, trait_type, output_path):
    """
    Merge chromosome-wise GWAS results into single file

    Concatenates gwas_results_chr1.txt through chr22.txt
    Compresses output with gzip
    """

def validate_saige_inputs(trait, ancestries, base_output_folder):
    """
    Check SAIGE outputs exist and are complete

    Returns: Dict with sample size info (n_cases, n_controls per ancestry)
    """

def extract_sample_sizes(filepath):
    """
    Parse n_cases and n_controls from GWAS results header

    SAIGE outputs include sample size in result file metadata
    """

def run_metal(trait, ancestries, base_output_folder):
    """
    Execute METAL meta-analysis across ancestries

    Creates METAL script file
    Submits dsub job with bwaxse/metal container
    Returns output path
    """
```

### 4. Post-Processing Functions
```python
def convert_snp_format(df):
    """
    Convert SAIGE variant IDs to gwaslab-compatible format

    SAIGE: chr:pos:ref:alt
    gwaslab: chr_pos_ref_alt
    """

def add_chr_pos_from_snpid(df):
    """
    Extract CHR and POS columns from variant ID string

    Required for gwaslab plotting functions
    """

def annotate_with_known_loci(lead_snps_df, gwas_catalog_df, window_kb=500):
    """
    Compare lead SNPs to GWAS Catalog known loci

    Adds 'KNOWN' column indicating overlap with published findings
    Uses genomic window matching (default 500kb)
    """

def process_gwas(sumstats_path, output_prefix, known_loci_df=None):
    """
    Complete gwaslab processing pipeline:
    1. Load summary statistics
    2. QC filtering (MAF, INFO thresholds)
    3. Manhattan plot
    4. QQ plot
    5. Regional plots for lead SNPs
    6. Export lead SNP list
    """

def split_gwas_file(input_path, output_prefix, max_mb=100):
    """
    Split large GWAS files into chunks for download

    Args:
        max_mb: Maximum file size in megabytes

    Creates numbered chunks: {prefix}_chunk_001.tsv.gz, etc.
    """
```

---

## Common Patterns and Architecture

### Data Flow
```
OMOP BigQuery Tables
  ↓
AouQueries.find_diagnosis_codes() → Jinja2 SQL Template → Polars DataFrame
  ↓
Person-level aggregation (≥1 vs ≥2 code thresholds)
  ↓
CSV Cohort Files + Exclusion Logic
  ↓
Merge with genomic metadata (PCs, sex, age)
  ↓
GWAS Metadata TSV (person_id, phenotype, covariates)
  ↓
SAIGE (distributed via dsub, 22 jobs per ancestry)
  ↓
Per-chromosome results → Combined TSV.GZ
  ↓
METAL meta-analysis (EUR + AFR)
  ↓
Post-processing (gwaslab) & Visualization
  ↓
Chunked download files
```

### File Naming Conventions
- **##** prefix indicates analysis stage (01-05)
- **min_1** vs **min_2** output folders distinguish case definitions:
  - `min_1/` - ≥ 1 diagnostic code
  - `min_2/` - ≥ 2 diagnostic codes
- Ancestry codes: EUR (European), AFR (African), AMR (Admixed American), EAS (East Asian), SAS (South Asian)

### Storage Structure
```
{workspace_bucket}/
├── data/cohorts/
│   ├── all_sarcoid_cases.csv              # Cases with exclusions applied
│   └── sarcoid_cases_without_exclusions.csv  # All cases before exclusions
├── saige_gwas/min_1/                      # ≥1 code analysis
│   ├── {ancestry}/condition__sarcoid/
│   │   ├── gwas_metadata.tsv              # Phenotype + covariates
│   │   ├── saige_null_model.rda           # Fitted null model
│   │   ├── saige_null_model.varianceRatio.txt
│   │   ├── swarm_gwas/
│   │   │   └── gwas_results_chr{1-22}.txt # Per-chromosome results
│   │   └── gwas_results.tsv.gz            # Combined results
│   └── metal/condition__sarcoid/
│       ├── condition__sarcoid_afr_eur1.tsv     # Meta-analysis results
│       └── condition__sarcoid_afr_eur1.tsv.info # Heterogeneity info
├── saige_gwas/min_2/                      # ≥2 codes analysis (same structure)
└── dsub/logs/                             # Job execution logs
```

### Ancestry-Specific Processing
Primary focus on EUR and AFR ancestries:
- **EUR**: European ancestry (largest sample size)
- **AFR**: African ancestry
- **Transancestry**: Attempted but may be underpowered for AMR/EAS/SAS

Ancestry assignment uses predicted ancestry with ancestry-specific principal components (ancPC1-20).

### Covariates Standard Set
```python
# Ancestry-specific (EUR, AFR)
covariates = ['sex_binary', 'end_of_study_age'] + ['ancPC{}'.format(str(x)) for x in range(1, 21)]
qCovarColList = ['sex_binary']  # Categorical covariate for sex chromosome handling

# Trans-ancestry
all_covariates = ['sex_binary', 'end_of_study_age'] + ['PC{}'.format(str(x)) for x in range(1, 17)]
all_qCovarColList = ['sex_binary']
```

### Statistical Thresholds
- **Genome-wide significance**: p < 5×10⁻⁸ (standard GWAS threshold)
- **Minimum MAC**: 20 (minimum minor allele count for inclusion)
- **LD-based clumping**: 500kb windows for lead SNP identification

### Docker Images Used
- **SAIGE**: `wzhou88/saige:1.3.6` - GWAS analysis
- **METAL**: Custom `bwaxse/metal` image - Meta-analysis
- **Base environment**: `us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.2.14`

### Cloud Computing Patterns
- **Null model fitting**: Single job per ancestry, c4-standard-8, non-preemptible
- **Chromosome-wise tests**: 22 parallel jobs per ancestry, c4-standard-8, preemptible
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

# Specialized tools
omop_unifier  # Custom OMOP data extraction and harmonization

# Utilities
jinja2  # SQL template rendering
subprocess  # Shell command execution
```

### External Tools (via Docker)
- **SAIGE** (1.3.6): GWAS for binary/quantitative traits with mixed models
- **METAL**: Inverse-variance weighted meta-analysis
- **PLINK2**: Genotype file manipulation

### Reference Data Sources
Located in `{src_bucket}/data/`:
- `ancestry_specific_pcs.parquet` - Pre-computed ancestry-specific principal components
- `genomic_metrics.tsv` - Sex ploidy and QC metrics
- `{ancestry}_merged_genotypes.{bed,bim,fam}` - LD-pruned variants for null model
- `genotypes_chr{1-22}.{bed,bim,fam}` - Full genotype data per chromosome

---

## Study-Specific Notes

### Sarcoidosis Case Definition

The study uses a **text-based search** strategy rather than hardcoded concept IDs:

```python
# Search across standard vocabularies
vocabularies = ["ICD9CM", "ICD10CM", "SNOMED"]
search_terms = ["sarcoid"]

# Captures variants like:
# - Sarcoidosis
# - Sarcoid granuloma
# - Pulmonary sarcoidosis
# - Cardiac sarcoidosis
# etc.
```

**Rationale**: Flexible text search is more robust to vocabulary updates and captures all sarcoid-related concepts.

### Dual Case Definition Strategy

**Purpose**: Sensitivity analysis to assess phenotype robustness

| Definition | Sample Size | Misclassification Risk | Use Case |
|------------|-------------|------------------------|----------|
| **≥ 1 code** | Higher | Moderate (includes rule-outs, single errors) | Discovery, maximum power |
| **≥ 2 codes** | Lower | Lower (requires confirmation) | Validation, high specificity |

**Expected Result**: If genetic associations are consistent across both definitions, this strengthens confidence in findings.

**Implementation Details**:
- Both analyses run identical SAIGE pipelines
- Only difference: `count_threshold` in `person_code_df()` method
- Separate output folders: `min_1/` vs `min_2/`

### Exclusion Criteria Logic

Competing diagnoses are removed from **controls only** (not from cases):

```python
# Pseudo-code
if has_exclusion_condition OR (exclusion_drug_date within 1yr of sarcoid_dx):
    if is_sarcoid_case:
        status = 'case'  # Keep as case
    else:
        status = 'excluded'  # Remove from control group
```

**Rationale**:
- Avoids misclassifying cases with similar granulomatous conditions
- Maintains case definition integrity
- Improves specificity without sacrificing sensitivity

### V-Code Special Handling

ICD9CM V codes (e.g., V42.0 for "Kidney replaced by transplant") require special processing:

1. Extract V codes from source tables
2. Join to `concept_relationship` table
3. Map to correct vocabulary via related concepts
4. Merge back with non-V code results

**Why Needed**: OMOP CDM maps V codes through relationships rather than direct vocabulary attribution.

---

## Usage Examples

### Run Complete ≥1 Code GWAS Pipeline
```python
# 1. Define cohort
# Run 01 Sarcoid Cohort.py to generate all_sarcoid_cases.csv

# 2. Run GWAS for EUR and AFR
# Run 02 SAIGE GWAS (≥ 1 code).py
# This creates null models and runs chromosome-wise tests

# 3. Meta-analyze across ancestries
# Run 03 METAL Meta-analysis (EUR AFR).py

# 4. Post-process and visualize
# Run 04 GWAS Post-Processing.py and 04 Summarize GWAS.py

# 5. Prepare for download
# Run 05 Download Sumstats.py
```

### Run ≥2 Codes Analysis for Sensitivity
```python
# Same workflow, use 02 SAIGE GWAS (≥ 2 codes).py instead
# Results go to saige_gwas/min_2/ folder
# Compare lead SNPs between min_1/ and min_2/ for robustness check
```

### Check Job Status
```python
check_dsub_status(user='bwaxse', full=False, age='3d')
```

### Extract Cases with Custom Definition
```python
# Initialize query builder
aou = AouQueries()

# Search for sarcoid diagnoses
aou.find_diagnosis_codes(
    search_terms=["sarcoid"],
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"]
)

# Execute query
results_df = aou.polars_gbq(aou.query_text)

# Create person-level summary (≥3 codes for very high specificity)
cases_df = aou.person_code_df(count_threshold=3)
```

### Access Results
```python
# Load EUR GWAS results (≥1 code)
eur_results = pl.read_csv(
    f'{my_bucket}/saige_gwas/min_1/eur/condition__sarcoid/gwas_results.tsv.gz',
    separator='\t'
)

# Load meta-analysis results
metal_results = pd.read_csv(
    f'{my_bucket}/saige_gwas/min_1/metal/condition__sarcoid/condition__sarcoid_afr_eur1.tsv',
    sep='\t'
)

# Compare ≥1 vs ≥2 code lead SNPs
min1_leads = pd.read_csv(f'{my_bucket}/saige_gwas/min_1/lead_snps.tsv', sep='\t')
min2_leads = pd.read_csv(f'{my_bucket}/saige_gwas/min_2/lead_snps.tsv', sep='\t')

# Check concordance
shared_loci = pd.merge(min1_leads, min2_leads, on='nearest_gene', suffixes=('_min1', '_min2'))
```

---

## Tips for Claude Code Reuse

1. **For new binary phenotype GWAS with variable case definitions**:
   - Adapt 01 cohort definition with new search terms/codes
   - Modify `count_threshold` parameter for sensitivity analyses
   - Run 02-05 pipeline unchanged

2. **For continuous trait GWAS**:
   - Keep 01 cohort extraction pattern
   - Modify `trait_type='binary'` to `'quantitative'` in 02
   - Adjust covariates as needed (may not need qCovarColList)

3. **For single-ancestry analysis**:
   - Skip 03 METAL step
   - Proceed directly to 04 post-processing on ancestry-specific results

4. **For flexible OMOP concept searching**:
   - Use `AouQueries` class with custom search terms
   - Mix text search, exact codes, and pattern matching as needed
   - Handle V codes automatically with special logic

5. **Common function extraction opportunities**:
   - `AouQueries` class for any OMOP data extraction (highly reusable!)
   - SAIGE null model + test workflow (02)
   - METAL meta-analysis wrapper (03)
   - gwaslab processing pipeline (04)

6. **Reusable bash scripts**:
   - SAIGE null model fitting
   - SAIGE chromosome-wise testing
   - METAL meta-analysis

---

## Comparison to HPV Pipeline

### Similarities:
- Same SAIGE/METAL workflow
- Identical covariate structure and statistical thresholds
- Distributed computing via dsub
- Multi-ancestry meta-analysis focus

### Key Differences:

| Feature | Sarcoid | HPV |
|---------|---------|-----|
| **Case definition** | ≥1 vs ≥2 codes (dual analysis) | Single definition with HIV exclusion |
| **Phenotypes** | Single condition (sarcoidosis) | Multiple anatomic sites (anal, vulvovaginal, cervical) |
| **Concept search** | Text-based ("sarcoid") | Specific concept IDs |
| **Exclusions** | Granulomatous diseases (from controls) | HIV+ participants (separate analysis) |
| **Ancestries** | EUR + AFR focus | EUR + AFR + AMR |
| **HLA analysis** | Not performed | Extensive 4-digit and 6-digit HLA imputation |
| **Output folders** | `min_1/` and `min_2/` | `v1_redo/` |
| **Class-based design** | `AouQueries` class for querying | Function-based approach |

---

## Limitations and Known Issues

1. **Text-based search**: May capture non-sarcoidosis granulomatous conditions despite exclusion logic
2. **Small sample sizes in AMR/EAS/SAS**: Transancestry analysis may be underpowered
3. **Exclusion criteria**: Specific excluded diagnoses not fully documented in code comments
4. **V-code complexity**: Special handling adds query complexity and potential failure points
5. **Case definition sensitivity**: ≥1 code may include rule-out diagnoses or coding errors
6. **No HLA analysis**: Unlike HPV pipeline, no specialized HLA imputation or visualization
7. **Missing covariate handling**: SAIGE handles missing data but may reduce effective sample size

---

## Metadata
- **Author**: Bennett Waxse
- **Platform**: All of Us Research Workbench
- **Data Version**: CDR v8
- **Last Updated**: December 2025
- **License**: MIT
