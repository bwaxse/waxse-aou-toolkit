# All of Us Data Science Toolkit - Architecture

**Purpose**: Cross-project patterns, reusable components, and architectural guidelines for All of Us Research

**Owner**: Bennett Waxse
**Platform**: All of Us Research Workbench (OMOP CDM v8)
**Last Updated**: December 2025

---

## Overview

This toolkit contains production-ready code for genomic and epidemiological analyses using All of Us data. Projects share common architectural patterns that enable code reuse and maintain consistency.

### Current Projects
1. **HPV** (`/hpv/`): GWAS of HPV-associated cancers with multi-ancestry meta-analysis
2. **BMI** (`/bmi/`): BMI data harmonization with temporal matching

### Planned Projects
- Additional GWAS phenotypes borrowing from HPV pipeline
- Other covariate harmonization modules borrowing from BMI approach
- PheWAS and association analyses

---

## Cross-Project Patterns

### 1. Dual File Format Strategy
**Pattern**: Every analysis has both `.py` and `.ipynb` files with identical code

**Rationale**:
- `.ipynb`: Interactive development, visualization, teaching
- `.py`: Version control, code review, imports, automation

**Implementation**:
```bash
# Generate from notebook
jupyter nbconvert --to python my_analysis.ipynb

# Or export from Jupyter UI
```

**Example**:
```
hpv/
├── B01.a HPV Cohort v0.ipynb    # Notebook for interactive use
└── B01.a HPV Cohort v0.py       # Python script (same code)
```

---

### 2. BigQuery → Polars Pipeline
**Pattern**: All OMOP data extraction uses BigQuery client → Arrow → Polars

**Code Template**:
```python
from google.cloud import bigquery
import polars as pl

def polars_gbq(query: str) -> pl.DataFrame:
    """
    Execute BigQuery SQL and return polars DataFrame
    Standard pattern used across all projects
    """
    client = bigquery.Client()
    query_job = client.query(query)
    rows = query_job.result()
    df = pl.from_arrow(rows.to_arrow())
    return df
```

**Why Polars?**
- Faster than pandas for large datasets
- Native BigQuery integration via Arrow (zero-copy)
- Better memory efficiency
- Type safety

**When to convert to pandas?**
- Complex datetime operations (polars datetime still maturing)
- Integration with scikit-learn, statsmodels
- Legacy code that expects pandas

**Example Usage**:
```python
# Extract cohort (polars)
cohort_query = f"""
SELECT person_id, birth_datetime, race_concept_id
FROM {version}.person
WHERE person_id IN (SELECT person_id FROM my_cohort)
"""
cohort_df = polars_gbq(cohort_query)

# Clean and transform (polars - fast)
cohort_df = cohort_df.with_columns([
    pl.col('birth_datetime').cast(pl.Date).alias('dob'),
    pl.col('race_concept_id').replace(race_mapping).alias('race')
])

# Convert for analysis (pandas)
cohort_pd = cohort_df.to_pandas()
```

---

### 3. OMOP Data Extraction Pattern
**Pattern**: Standardized approach to querying OMOP CDM

**Core OMOP Tables Used**:
```sql
-- Person table: Demographics
{version}.person
  person_id, birth_datetime, race_concept_id, ethnicity_concept_id,
  gender_concept_id, sex_at_birth_concept_id

-- Condition occurrence: Diagnoses
{version}.condition_occurrence
  person_id, condition_concept_id, condition_start_date,
  condition_source_value, condition_source_concept_id

-- Observation: Survey responses (PPI)
{version}.observation
  person_id, observation_concept_id, value_source_concept_id,
  observation_date

-- Measurement: Labs, vitals, anthropometrics
{version}.measurement
  person_id, measurement_concept_id, value_as_number,
  unit_concept_id, measurement_date

-- Concept: Lookup table for codes
{version}.concept
  concept_id, concept_name, concept_code, vocabulary_id

-- CB search tables: All of Us cohort builder
{version}.cb_search_person
  person_id, has_whole_genome_variant, has_lr_whole_genome_variant

{version}.cb_search_all_events
  person_id, concept_id, entry_date, is_standard
```

**Standard Query Pattern**:
```sql
-- 1. Define eligible participants (WITH CTE)
WITH distinct_participants AS (
    SELECT DISTINCT person_id
    FROM {version}.observation
    WHERE observation_concept_id IN (...)
    -- Apply inclusion criteria
)

-- 2. Join to get demographics
SELECT
    p.person_id,
    p.birth_datetime,
    race.concept_name AS race,
    ethnicity.concept_name AS ethnicity
FROM {version}.person p
INNER JOIN distinct_participants dp ON p.person_id = dp.person_id
LEFT JOIN {version}.concept race ON p.race_concept_id = race.concept_id
LEFT JOIN {version}.concept ethnicity ON p.ethnicity_concept_id = ethnicity.concept_id
```

**Concept ID Lookup**:
```python
# Always join concept table to get human-readable names
concept_query = f"""
SELECT concept_id, concept_name, concept_code, vocabulary_id
FROM {version}.concept
WHERE vocabulary_id IN ('ICD10CM', 'ICD9CM', 'SNOMED', 'LOINC')
  AND concept_name LIKE '%diabetes%'
"""
```

---

### 4. Covariate Engineering Pattern
**Pattern**: Transform categorical OMOP data into analysis-ready binary indicators

**Template** (from HPV and BMI projects):
```python
def engineer_covariates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create binary indicator variables for categorical covariates

    Common categories:
    - Demographics: age bins, race/ethnicity, sex
    - SES: employment, income, education, insurance
    - Clinical: comorbidities, medications
    """

    # Age binning
    age_bins = [17, 24, 34, 44, 54, 64, 74, float('inf')]
    age_labels = ['17_24', '25_34', '35_44', '45_54', '55_64', '65_74', '75_plus']
    for i, label in enumerate(age_labels):
        df[label] = (
            (df['age'] >= age_bins[i]) & (df['age'] < age_bins[i+1])
        ).astype(int)

    # Race indicators (multi-hot encoding)
    df['black'] = (df['race'] == 'Black').astype(int)
    df['white'] = (df['race'] == 'White').astype(int)
    df['asian'] = (df['race'] == 'Asian').astype(int)
    df['hispanic'] = (df['ethnicity'] == 'Hispanic or Latino').astype(int)

    # Sex (binary encoding with reference category)
    df['female'] = (df['sex_at_birth'] == 'Female').astype(int)

    # SES indicators
    df['employed'] = (df['employment_status'] == 'employed').astype(int)
    df['college_grad'] = (df['education'] == 'coll_grad').astype(int)

    return df
```

**Usage in GWAS**:
```python
# Prepare phenotype file for SAIGE
covariates = ['female', 'black', 'white', 'asian', 'hispanic'] + \
             [f'ancPC{i}' for i in range(1, 21)] + \
             ['17_24', '25_34', '35_44', '45_54', '55_64', '65_74', '75_plus']

gwas_metadata = engineer_covariates(demographics_df)
gwas_metadata[['person_id', 'phenotype'] + covariates].to_csv('gwas_metadata.tsv', sep='\t')
```

---

### 5. Distributed Computing with dsub
**Pattern**: Google Cloud Batch job submission for parallelizable tasks

**Core Helper Function**:
```python
def dsub_script(
    machine_type: str,      # 'c4-standard-8', 'n2d-standard-16'
    envs: Dict[str, str],   # Environment variables
    in_params: Dict[str, str],   # Input file mappings (local name → GCS path)
    out_params: Dict[str, str],  # Output file mappings
    boot_disk: int = 100,   # GB
    disk_size: int = 150,   # GB
    image: str = '...',     # Docker image
    script: str = '...',    # Bash script to run
    preemptible: bool = True  # Use preemptible instances?
) -> None:
    """
    Submit job to Google Cloud Batch

    Standard pattern:
    1. Define machine specs
    2. Map environment variables (parameters)
    3. Map input files (downloads from GCS to /mnt/data/input/)
    4. Map output files (uploads from container to GCS)
    5. Specify Docker image with analysis tools
    6. Provide bash script that orchestrates the analysis
    """
    # Build dsub command
    dsub_cmd = build_dsub_command(...)
    os.system(dsub_cmd)
```

**Example - SAIGE GWAS** (per-chromosome parallelization):
```python
# Run 22 jobs in parallel (one per chromosome)
for chrom in range(1, 23):
    dsub_script(
        machine_type='c4-standard-8',
        envs={'CHR': chrom},
        in_params={
            'INPUT_VCF': f'{bucket}/genotypes_chr{chrom}.vcf.gz',
            'INPUT_NULL_MODEL': f'{bucket}/saige_null_model.rda'
        },
        out_params={
            'OUTPUT_FILE': f'{bucket}/gwas_chr{chrom}.txt'
        },
        image='wzhou88/saige:1.3.6',
        script='run_saige_chrom.sh'
    )
```

**Job Monitoring**:
```python
def check_dsub_status(user=None, age='3d'):
    """Check status of running/recent jobs"""
    cmd = f"dstat --provider google-batch --user {user} --age {age} --status '*'"
    subprocess.run(cmd, shell=True)
```

**Common Machine Types**:
- `c4-standard-8`: CPU-intensive (GWAS, meta-analysis) - 8 vCPUs
- `n2d-standard-16`: Memory-intensive (large data processing) - 16 vCPUs
- `c2-standard-60`: High-throughput (many parallel tasks) - 60 vCPUs

**Cost Optimization**:
- Use `preemptible=True` for interruptible workloads (3-4x cheaper)
- Right-size machine types (don't over-provision)
- Use pd-ssd disks for I/O-intensive tasks, pd-standard otherwise

---

### 6. SAIGE GWAS Workflow
**Pattern**: Two-step GWAS using SAIGE (used across phenotypes)

**Step 1: Fit Null Model** (once per phenotype per ancestry)
```bash
#!/bin/bash
# run_saige_null_model.sh

step1_fitNULLGLMM.R \
    --plinkFile=${INPUT_BED/.bed/} \    # Pruned genotypes for GRM
    --phenoFile=${INPUT_METADATA} \     # Phenotype + covariates TSV
    --phenoCol=${TRAIT} \               # Column name for phenotype
    --covarColList=${COVARIATES} \      # Comma-separated covariate list
    --sampleIDColinphenoFile=person_id \
    --traitType=${TRAIT_TYPE} \         # 'binary' or 'quantitative'
    --nThreads=8 \
    --outputPrefix=${OUT_PREFIX}        # Outputs: .rda, .varianceRatio.txt
```

**Step 2: Hypothesis Testing** (parallelized per chromosome)
```bash
#!/bin/bash
# run_saige_chrom.sh

step2_SPAtests.R \
    --vcfFile=${INPUT_VCF} \            # Imputed genotypes for chromosome
    --vcfFileIndex=${INPUT_VCF_IX} \
    --vcfField=DS \                     # Dosage field
    --chrom=${CHR} \
    --GMMATmodelFile=${INPUT_NULL_RDA} \
    --varianceRatioFile=${INPUT_NULL_VARRAT} \
    --minMAC=20 \                       # Minimum minor allele count
    --SAIGEOutputFile=${OUTPUT_FILE}    # GWAS results
```

**Reusable Python Wrapper**:
```python
def run_saige_gwas_pipeline(
    trait: str,
    trait_type: str,  # 'binary' or 'quantitative'
    ancestries: List[str],  # ['eur', 'afr', 'amr']
    cohort_dfs: List[pd.DataFrame],
    covariates: List[str],
    output_folder: str
):
    """
    Complete SAIGE GWAS pipeline

    1. Prepare metadata (merge cohorts with genomic PCs)
    2. Fit null models (one job per ancestry)
    3. Run tests (22 jobs per ancestry, per chromosome)
    4. Combine results
    """
    for anc in ancestries:
        # Prepare metadata
        prepare_gwas_metadata(cohort_dfs, anc, output_folder)

        # Fit null model (blocks until complete)
        run_saige_null(anc, trait, trait_type, covariates, ...)

        # Run tests (parallelized)
        run_saige_test(anc, trait, chroms=range(1,23), ...)

    # Combine chromosome-wise results
    for anc in ancestries:
        combine_gwas_results(anc, trait, output_folder)
```

---

### 7. Meta-Analysis with METAL
**Pattern**: Inverse-variance weighted meta-analysis across ancestries

**METAL Configuration Script**:
```metal
# metal_script.txt
MARKER vid                     # Variant ID
WEIGHT n_total                 # Sample size for weighting
ALLELE effect_allele other_allele
FREQ effect_allele_frequency
PVAL p_value
EFFECT beta
STDERR standard_error
SEPARATOR TAB

SCHEME STDERR                  # Inverse-variance weighting
AVERAGEFREQ ON                 # Average allele frequencies
MINMAXFREQ ON                  # Report min/max frequencies

PROCESS ancestry1_results.tsv
PROCESS ancestry2_results.tsv
PROCESS ancestry3_results.tsv

ANALYZE HETEROGENEITY          # Test for effect heterogeneity
```

**Python Wrapper**:
```python
def run_metal_meta_analysis(
    trait: str,
    ancestries: List[str],
    input_files: Dict[str, str],  # {ancestry: file_path}
    output_dir: str
):
    """
    Run METAL meta-analysis

    1. Validate inputs (check files exist, extract sample sizes)
    2. Create METAL script
    3. Submit dsub job with METAL container
    4. Parse results
    """
    # Validate
    sample_sizes = validate_inputs(input_files)

    # Create METAL script
    metal_script = create_metal_script(input_files, sample_sizes)

    # Run via dsub
    dsub_script(
        machine_type='n2d-standard-8',
        in_params=input_files,
        out_dirs={'OUTPUT_PATH': output_dir},
        envs={'TRAIT': trait, 'SAMPLE_SIZES': ','.join(map(str, sample_sizes))},
        image='bwaxse/metal',
        script='run_metal.sh'
    )
```

---

### 8. Temporal Data Matching
**Pattern**: Match time-stamped measurements to cohort reference times (from BMI project, generalizable)

**Core Algorithm**:
```python
def hierarchical_temporal_matching(
    cohort_df: pd.DataFrame,      # Has person_id, time_zero
    measurement_df: pd.DataFrame,  # Has person_id, date, value
    time_col: str = 'time_zero',
    include_post: bool = True,
    max_prior_days: int = 365,
    max_post_days: int = 90
) -> pd.DataFrame:
    """
    Match measurements to cohort with hierarchical preference

    Hierarchy:
    1. Prior measurements (closer to time_zero = higher priority)
    2. Post measurements (heavily penalized, fallback only)

    Returns: cohort_df with matched measurements and quality flags
    """
    # Merge all possible matches
    merged = pd.merge(cohort_df, measurement_df, on='person_id', how='left')

    # Calculate temporal offset
    merged['days_diff'] = (merged['date'] - merged[time_col]).dt.days

    # Filter to valid time windows
    valid = (
        ((merged['days_diff'] <= 0) & (abs(merged['days_diff']) <= max_prior_days)) |
        ((merged['days_diff'] > 0) & (merged['days_diff'] <= max_post_days) & include_post)
    )
    merged = merged[valid]

    # Hierarchical priority (lower = better)
    merged['priority'] = np.where(
        merged['days_diff'] <= 0,
        abs(merged['days_diff']),      # Prior: prefer closest
        1000 + abs(merged['days_diff'])  # Post: heavily penalized
    )

    # Keep best match per person
    merged = merged.sort_values(['person_id', 'priority'])
    result = merged.drop_duplicates('person_id', keep='first')

    # Add quality flags
    result['timing'] = np.where(result['days_diff'] <= 0, 'prior', 'post')
    result['quality'] = pd.cut(
        abs(result['days_diff']),
        bins=[0, 30, 90, 180, float('inf')],
        labels=['excellent', 'good', 'acceptable', 'poor']
    )

    return result
```

**Reusable for**:
- Lab values (HbA1c, cholesterol, etc.)
- Vitals (blood pressure, heart rate)
- Anthropometrics (BMI, waist circumference)
- Any timestamped covariate

**Strategies for Different Study Types**:
```python
# GWAS: Strict prior only (avoid reverse causation)
gwas_match = hierarchical_temporal_matching(
    cohort, measurements,
    include_post=False,
    max_prior_days=365
)

# Drug study: Recent baseline (close to treatment start)
drug_match = hierarchical_temporal_matching(
    cohort, measurements,
    include_post=False,
    max_prior_days=90
)

# Longitudinal: Flexible windows
longitudinal_match = hierarchical_temporal_matching(
    cohort, measurements,
    include_post=True,
    max_prior_days=730,
    max_post_days=180
)
```

---

### 9. Table 1 Generation
**Pattern**: Descriptive statistics tables following epidemiology standards

**Template**:
```python
def generate_table1(
    df: pd.DataFrame,
    group_col: str = None,  # Optional stratification variable
    suppress_small_counts: bool = True  # All of Us data use policy
) -> pd.DataFrame:
    """
    Generate Table 1 with descriptive statistics

    Continuous variables: Median (IQR)
    Categorical variables: n (%)
    Small counts: Suppress if < 20 (All of Us policy)
    """

    def describe_continuous(series):
        """Format: median (Q1, Q3)"""
        q1, median, q3 = series.quantile([0.25, 0.5, 0.75])
        return f"{median:.1f} ({q1:.1f}, {q3:.1f})"

    def describe_categorical(series, total):
        """Format: n (%) or <20 (<%)"""
        count = series.sum()
        pct = 100 * count / total

        if suppress_small_counts and 1 <= count < 20:
            return f"<20 (<{20/total*100:.1f})"
        else:
            return f"{count} ({pct:.1f})"

    # Build table
    results = {}

    # Continuous
    results['Age, median (IQR)'] = describe_continuous(df['age'])

    # Categorical
    for var in ['female', 'white', 'black', 'asian', 'hispanic', 'employed']:
        results[f'{var.capitalize()}, n (%)'] = describe_categorical(
            df[var], len(df)
        )

    return pd.DataFrame(results, index=[0])
```

---

### 10. HLA Visualization
**Pattern**: Specialized Manhattan and forest plots for HLA region

**Manhattan Plot with Gene Tracks**:
```python
def create_hla_manhattan(
    gwas_results: pl.DataFrame,  # Genome-wide results
    hla_results: pl.DataFrame,    # HLA allele results
    gene_annotations: pl.DataFrame,  # HLA gene positions
    ancestry: str
):
    """
    Manhattan plot overlaying HLA alleles on genomic background

    Features:
    - Gray background: GWAS SNPs in HLA region
    - Colored points: HLA alleles (Class I = red, Class II = blue)
    - Gene track: HLA gene positions below x-axis
    - Significance lines: Bonferroni, genome-wide (5e-8)
    """
    fig, ax = plt.subplots()

    # Background SNPs
    ax.scatter(gwas_results['POS'], -np.log10(gwas_results['P']),
               c='lightgray', s=5, alpha=0.4)

    # HLA alleles
    for hla_gene in hla_class_i_genes:
        gene_data = hla_results.filter(pl.col('gene') == hla_gene)
        ax.scatter(gene_data['pos'], -np.log10(gene_data['P']),
                  c='red', s=20, label=hla_gene)

    # Gene track (below x-axis)
    for gene in gene_annotations.iter_rows(named=True):
        rect = Rectangle((gene['start'], -2.5),
                        gene['end'] - gene['start'], 0.6,
                        facecolor=gene_color[gene['name']])
        ax.add_patch(rect)

    return fig, ax
```

**Forest Plot for Multi-Ancestry Comparison**:
```python
def create_forest_plot(
    results_dict: Dict[str, pl.DataFrame],  # {ancestry: results}
    bonf_thresholds: Dict[str, float]
):
    """
    Forest plot showing odds ratios across ancestries

    Features:
    - Rows: HLA alleles (significant in ≥1 ancestry)
    - Columns: Ancestries (offset vertically)
    - Points: OR with 95% CI error bars
    - Colors: By ancestry
    - Shapes: EUR=circle, AFR=square, AMR=triangle
    - Alpha: 1.0 if significant, 0.2 if not
    """
    fig, ax = plt.subplots()

    for ancestry, results in results_dict.items():
        y_offset = ancestry_offsets[ancestry]
        significant = results['P'] < bonf_thresholds[ancestry]

        for _, row in results.iterrows():
            y_pos = marker_positions[row['allele']] + y_offset
            alpha = 1.0 if row['P'] < bonf_thresholds[ancestry] else 0.2

            # Plot CI
            ax.plot([row['OR_lower'], row['OR_upper']], [y_pos, y_pos],
                   color=ancestry_colors[ancestry], alpha=alpha, linewidth=2)

            # Plot point estimate
            ax.scatter(row['OR'], y_pos,
                      marker=ancestry_markers[ancestry],
                      color=ancestry_colors[ancestry],
                      s=30, alpha=alpha)

    ax.axvline(x=1, color='black', linewidth=2)  # Null line
    ax.set_xscale('log')
    return fig, ax
```

---

## Reusable Function Library

### Core Utilities (use across all projects)

#### 1. Data Extraction
```python
# polars_gbq() - BigQuery → polars DataFrame
from bmi_functions import polars_gbq

# Usage
df = polars_gbq(f"""
SELECT * FROM {version}.person LIMIT 1000
""")
```

#### 2. Distributed Computing
```python
# dsub_script() - Submit Google Cloud Batch job
from hpv.B02... import dsub_script, check_dsub_status

# Usage
dsub_script(
    machine_type='c4-standard-8',
    envs={'PARAM': 'value'},
    in_params={'INPUT': 'gs://bucket/file.txt'},
    out_params={'OUTPUT': 'gs://bucket/output.txt'},
    script='my_analysis.sh'
)

# Monitor
check_dsub_status(age='1d')
```

#### 3. Temporal Matching
```python
# hierarchical_temporal_matching() - Match measurements to reference time
from bmi_functions import hierarchical_temporal_matching

# Usage
matched_cohort = hierarchical_temporal_matching(
    cohort_df,
    measurement_df,
    time_col='enrollment_date',
    include_post=False,
    max_prior_days=365
)
```

#### 4. Covariate Engineering
```python
# update_covariates() - Create binary indicators
from hpv.B01... import update_covariates

# Usage
df_with_covariates = update_covariates(demographics_df, sarcoid=False)
```

#### 5. SAIGE Pipeline
```python
# run_saige_null(), run_saige_test() - GWAS workflow
from hpv.B02... import run_saige_null, run_saige_test

# Usage
run_saige_null('eur', 'condition__diabetes', 'binary', covariates, ...)
run_saige_test('eur', 'condition__diabetes', script='run_saige_chrom.sh')
```

#### 6. Meta-Analysis
```python
# run_metal() - Combine ancestry-specific GWAS
from hpv.B03... import run_metal, combine_gwas_results

# Usage
combine_gwas_results('eur', 'condition__diabetes', trait_type='binary', ...)
run_metal('condition__diabetes', ['eur', 'afr', 'amr'], output_folder)
```

#### 7. Table 1
```python
# describe_group() - Descriptive statistics
from hpv.B01... import describe_group

# Usage
table1 = describe_group(cohort_df, sarcoid=False)
```

---

## Common Workflows

### Workflow 1: New Binary Phenotype GWAS
**Borrow from**: HPV project (B01-B03)

**Steps**:
1. **Define cohort** (adapt `B01.a HPV Cohort v0.py`):
   - Replace HPV concept IDs with your phenotype concept IDs
   - Keep demographics/covariate extraction code unchanged
   - Use `update_covariates()` function as-is

2. **Run SAIGE** (adapt `B02.a Anal HPV SAIGE GWAS v0.py`):
   - Update `traits` dict: `{'condition__my_phenotype': 'binary'}`
   - Keep `covariates` list same (or customize)
   - Run `run_saige_null()` and `run_saige_test()` unchanged

3. **Meta-analyze** (use `B03 METAL Meta-analysis` unchanged):
   - Just update trait name in function calls
   - All logic is trait-agnostic

4. **Visualize** (adapt `B06.1 Plot HLA Results` if HLA-relevant):
   - Swap in your GWAS results
   - Customizegene regions of interest

**Time saved**: ~80% (mostly just concept ID swapping)

---

### Workflow 2: Quantitative Trait GWAS
**Borrow from**: HPV project with modifications

**Steps**:
1. Extract continuous phenotype (e.g., HbA1c) instead of binary diagnosis
2. Change `trait_type='binary'` → `'quantitative'` in SAIGE calls
3. Update METAL column mappings (quantitative traits have 'N' instead of 'N_cases'/'N_controls')
4. All other code unchanged

---

### Workflow 3: New Covariate Harmonization
**Borrow from**: BMI project

**Steps**:
1. **Identify concept IDs** for your measurement (e.g., blood pressure):
   ```sql
   SELECT concept_id, concept_name
   FROM {version}.concept
   WHERE concept_name LIKE '%blood pressure%'
     AND domain_id = 'Measurement'
   ```

2. **Adapt extraction** (copy `extract_bmi_data()` from `bmi_functions.py`):
   - Replace BMI/weight/height concept IDs with yours
   - Keep unit handling structure

3. **Reuse cleaning pipeline**:
   - `clean_units()` - update problematic unit list
   - `apply_clean_outliers()` - use as-is (change column names)
   - `convert_units()` - update conversion factors
   - `calculate_bmi_with_validation()` - adapt validation logic

4. **Temporal matching** - use unchanged:
   - `hierarchical_temporal_matching()` works for any measurement

---

### Workflow 4: PheWAS Analysis
**Borrow from**: HPV project (B01.a)

**Steps**:
1. Extract cohort with `polars_gbq()`
2. Use `update_covariates()` to engineer covariates
3. Use PheTK package (already integrated):
   ```python
   from PheTK.PheWAS import PheWAS

   phewas = PheWAS(
       phecode_version="X",
       cohort_csv_path="my_cohort.csv",
       covariate_cols=["age", "sex", "race"],
       independent_variable_of_interest="exposure",
       min_cases=50
   )
   phewas.run()
   ```

---

## File Organization Conventions

### Project Structure
```
project_name/
├── CLAUDE.md                    # Project-specific documentation
├── B01_descriptive_name.py      # Stage 1: Cohort definition
├── B01_descriptive_name.ipynb
├── B02_descriptive_name.py      # Stage 2: Primary analysis
├── B02_descriptive_name.ipynb
├── B03_descriptive_name.py      # Stage 3: Secondary analysis
├── ...
└── outputs/                     # (in GCS, not local)
    ├── cohort.csv
    ├── results.tsv
    └── figures/
```

### Naming Conventions
- **B##** prefix: Analysis stage (B01, B02, ..., B99)
  - Defines execution order
  - Later stages depend on earlier stages
- **Descriptive name**: What the analysis does
- **Version suffix** (optional): `v0`, `v1`, `v2`, or date `(5.20)`
- **Iteration suffix**: `redo`, `final`, `updated`

### GCS Storage Structure
```
{workspace_bucket}/
├── cohorts/
│   ├── {project}/
│   │   └── {phenotype}_cohort.csv
├── saige_gwas/
│   ├── {version}/
│   │   ├── {ancestry}/
│   │   │   ├── {trait}/
│   │   │   │   ├── gwas_metadata.tsv
│   │   │   │   ├── saige_null_model.rda
│   │   │   │   ├── gwas/
│   │   │   │   │   └── gwas_results_chr{1-22}.txt
│   │   │   │   └── gwas_results.tsv.gz
├── metal/
│   ├── {trait}/
│   │   └── {trait}_{ancestry1}_{ancestry2}_...1.tsv
├── HLA/
│   └── {resolution}_digit_{phenotype}_{ancestry}_HLA_results.csv
└── dsub/logs/
    └── {job-name}/{user}/{date}/{job-id}.log
```

---

## Docker Images

### Standard Images
- **All of Us base**: `us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.2.14`
  - Pre-installed: pandas, numpy, matplotlib, google-cloud-bigquery
  - Use for: Data extraction, cleaning, visualization

- **SAIGE**: `wzhou88/saige:1.3.6`
  - Pre-installed: SAIGE, PLINK2, R
  - Use for: GWAS analysis

### Custom Images (built by user)
- **METAL**: `{artifact_registry}/bwaxse/metal`
  - Pre-installed: METAL, awk, gzip
  - Use for: Meta-analysis

### Building Custom Images
```dockerfile
# Dockerfile
FROM ubuntu:20.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    wget

# Install METAL
RUN wget http://csg.sph.umich.edu/abecasis/metal/download/Linux-metal.tar.gz && \
    tar -xzf Linux-metal.tar.gz && \
    mv metal /usr/local/bin/

CMD ["/bin/bash"]
```

```bash
# Build and push
docker build -t {artifact_registry}/bwaxse/metal .
docker push {artifact_registry}/bwaxse/metal
```

---

## All of Us Data Policies

### Cell Suppression Rules
**Policy**: Suppress counts 1-19 in all outputs

**Implementation**:
```python
def suppress_small_counts(count, threshold=20):
    """Return '<20' if count is 1-19, otherwise actual count"""
    if 1 <= count < threshold:
        return f"<{threshold}"
    else:
        return str(count)

# Usage in Table 1
table1['n'] = table1['n'].apply(suppress_small_counts)
```

### Genomic Data Access
- Whole genome sequencing: `has_whole_genome_variant = 1`
- Long-read sequencing: `has_lr_whole_genome_variant = 1`
- Always filter to participants with genomic data before GWAS

### Workspace Variables
```python
import os

# CDR version (changes over time)
version = os.environ['WORKSPACE_CDR']
# Example: '`nih-nci-dceg-cc-strides-subs-sc-us-central1-1-cdr.cdr_controlled`'

# Workspace bucket (user-specific persistent storage)
my_bucket = os.environ['WORKSPACE_BUCKET']
# Example: 'gs://fc-secure-...'

# Google Cloud project
project = os.environ['GOOGLE_PROJECT']
# Example: 'terra-vpc-sc-...'

# User email
user_email = os.environ['OWNER_EMAIL']
# Example: 'user@researchallofus.org'
```

---

## Dependencies

### Python Libraries (Standard)
```python
# Data manipulation
pandas >= 1.3
polars >= 0.15
numpy >= 1.20
scipy >= 1.7

# Cloud integration
google-cloud-bigquery >= 2.0
google-cloud-storage >= 1.40

# Visualization
matplotlib >= 3.3
seaborn >= 0.11

# Utilities
pyarrow >= 6.0  # BigQuery → Arrow → polars
jinja2 >= 3.0   # Template rendering
```

### Domain-Specific Libraries
```python
# GWAS visualization
gwaslab >= 3.4

# PheWAS
PheTK  # All of Us phenotype toolkit

# Statistical genetics
# (installed in Docker containers, not in Python env)
# SAIGE, PLINK2, METAL
```

---

## Future Enhancements

### Planned Components
1. **Polygenic Risk Score (PRS)** calculation module
   - Borrowing from SAIGE output format
   - Standardized PRS calculation across ancestries

2. **Survival analysis** pipeline
   - Time-to-event outcomes from OMOP
   - Cox proportional hazards with genetic covariates

3. **Multi-trait GWAS**
   - MTAG (multi-trait analysis of GWAS)
   - Genetic correlation estimation

4. **Annotation pipeline**
   - Functional annotation of GWAS hits
   - Gene-based testing
   - Pathway enrichment

### Code Reuse Strategy
As new projects are added, identify:
1. **Common functions** → Extract to shared `utils/` module
2. **Common workflows** → Create template notebooks
3. **Common configs** → Parameterize with YAML/JSON

---

## Tips for Claude Code Generation

### When Generating New Code
1. **Always start with an existing project** as template
2. **Identify the stage** in the pipeline (cohort, GWAS, visualization, etc.)
3. **Find the closest match** in existing files (e.g., binary phenotype → HPV)
4. **Swap specific elements**:
   - Concept IDs for new phenotype
   - File paths for new project name
   - Keep all infrastructure code unchanged

### Code Review Checklist
- [ ] Uses `polars_gbq()` for BigQuery queries?
- [ ] Applies cell suppression (`< 20` for counts 1-19)?
- [ ] Uses dsub for distributed computing (not manual loops)?
- [ ] Follows file naming convention (B## prefix)?
- [ ] Has both .py and .ipynb versions?
- [ ] Documents key functions with docstrings?
- [ ] Stores outputs in GCS (not local filesystem)?

### Common Mistakes to Avoid
1. **Don't** hard-code workspace bucket paths (use `os.environ['WORKSPACE_BUCKET']`)
2. **Don't** use pandas for BigQuery (use polars via `polars_gbq()`)
3. **Don't** run GWAS in serial (use dsub parallelization)
4. **Don't** skip cell suppression (All of Us policy violation)
5. **Don't** assume CDR version (always use `os.environ['WORKSPACE_CDR']`)

---

## Metadata
- **Author**: Bennett Waxse
- **Organization**: All of Us Research Program
- **Last Updated**: December 2025
- **License**: MIT
- **Status**: Living document - updated as new projects are added
