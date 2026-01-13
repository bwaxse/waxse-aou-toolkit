# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Repository Overview

**Purpose**: Production-ready code for genomic and epidemiological analyses using All of Us Research Workbench data

**Platform**: All of Us Research Workbench (OMOP CDM v8) on Google Cloud Platform
**Data Model**: OMOP Common Data Model (CDM)
**Primary Language**: Python (Jupyter notebooks + exported .py scripts)
**License**: MIT

This repository contains multiple completed and active research projects (HPV GWAS, Sarcoid GWAS, wearables analysis, medication extraction) with comprehensive documentation and reusable components.

---

## Getting Started

### Step 1: Initialize Workspace (First Thing in Each Session)

Before running any analysis notebooks, you **must** initialize your workspace environment variables.

**Run this notebook first**: `_reference/verily/00_setup_workspace.ipynb`

This single notebook (~2 seconds to run) will set up:
- `GOOGLE_CLOUD_PROJECT` - Your Google Cloud project ID
- `WORKSPACE_CDR` - BigQuery dataset with All of Us CDR v8
- `WORKSPACE_BUCKET` - Persistent workspace GCS bucket
- `WORKSPACE_TEMP_BUCKET` - Temporary workspace GCS bucket

**Workflow:**
1. Open a new JupyterLab session in Verily Workbench
2. Run `_reference/verily/00_setup_workspace.ipynb` once
3. All other notebooks can now access these variables via `os.environ`
4. Run this again if you start a new JupyterLab session

**In your analysis notebooks**, access variables like:
```python
import os

WORKSPACE_CDR = os.environ['WORKSPACE_CDR']
WORKSPACE_BUCKET = os.environ['WORKSPACE_BUCKET']
```

See the notebook itself for more details and troubleshooting.

---

## Quick Decision Guide

**READ PROJECT-SPECIFIC CLAUDE.md FILES FIRST**

Each project has detailed documentation with full API references, workflows, and reusable functions:
- `hpv/CLAUDE.md` - Multi-ancestry GWAS pipeline (7-stage)
- `sarcoid/CLAUDE.md` - GWAS with dual case definitions (5-stage)
- `lc_wearables/CLAUDE.md` - Fitbit wearables data processing

**When to use what:**
- New GWAS project → Read `hpv/CLAUDE.md` or `sarcoid/CLAUDE.md` first, use as template
- OMOP concept searching → Read `sarcoid/CLAUDE.md`, use `AouQueries` class from `sarcoid/01 Sarcoid Cohort.py`
- Wearables/device data → Read `lc_wearables/CLAUDE.md`
- ICD code extraction → See `trusted_queries/icd_codes.py` for `phetk_icd_query()` function
- Medication extraction → See `meds/EXPORT 01.12_extraction_medication_ingredients_rx_norm.py`

---

## File Structure

### File Format Strategy (.py Source + Generated Notebooks)

Analysis code is maintained in **Python scripts** (`.py` files) for version control and Claude context efficiency:

**Source files (tracked in git):**
- `.py` files - **Source of truth**, in project roots (`hpv/`, `sarcoid/`, etc.)
  - Clean Python code for reading and editing
  - Better git diffs and code review
  - Function imports across projects
  - Searchable and diff-friendly
  - **This is what Claude reads** (saves 60% context tokens)

**Notebook files (generated on demand):**
- `.ipynb` files - **Created for Verily Workbench upload**, stored in `<project>/notebooks/` subdirectories
  - Tracked in git for collaboration
  - Masked from Claude via `.claudeignore` (saves tokens)
  - Generated when needed: request "convert X.py to notebook"
  - Upload to Verily, run interactively, then regenerate as needed

**Workflow:**
1. **Develop**: Claude writes/edits `.py` files in project root
2. **Upload**: Request "create notebook for X.py" when ready to run in Verily
3. **Run**: Upload generated `.ipynb` to Verily Workbench and execute interactively
4. **Iterate**: Edit `.py` file, regenerate `.ipynb` as needed

**Why this approach?**
- `.py` files are token-efficient for Claude context
- Both formats stay in git (good for collaboration)
- Clear separation: code (in .py) vs execution format (.ipynb)
- Flexibility: regenerate notebooks anytime from source

### Project Directories

```
waxse-aou-toolkit/
├── hpv/              # HPV GWAS pipeline (production, multi-ancestry, 7-stage)
├── sarcoid/          # Sarcoid GWAS pipeline (production, dual case definitions, 5-stage)
├── genomics/         # Genomics analysis (ancestry PCA, genotype preparation)
├── lc_wearables/     # Wearables data analysis (Fitbit processing)
├── autoencoder_fever/# Autoencoder analysis
├── meds/             # Medication extraction utilities
├── trusted_queries/  # Reusable query functions (ICD code extraction)
├── lc_kir/           # KIR genetics project (early stage)
└── _reference/       # Reference materials and documentation
    ├── all_of_us_tables/  # OMOP table schemas
    ├── phecode/           # Phecode mappings
    ├── verily/            # Verily platform examples
    └── resource_monitoring_template.md
```

### Converting Python Scripts to Notebooks

When you're ready to upload a `.py` file to Verily Workbench for interactive execution:

**Request format:**
"Convert `hpv/B01.a HPV Cohort v0.py` to a notebook"

**Claude will:**
1. Read the `.py` file
2. Generate a properly formatted `.ipynb` file
3. Save to `hpv/notebooks/B01.a HPV Cohort v0.ipynb`
4. Preserve code structure, comments, markdown cells, and cell order

**Then you can:**
- Navigate to the `notebooks/` subdirectory
- Upload the `.ipynb` file to your Verily Workbench
- Run cells interactively
- Modify as needed in Verily

**Note**: Notebooks are ephemeral for interactive work. If you want to save changes back to code:
1. Edit the `.py` file directly
2. Request a fresh notebook conversion
3. Upload the new version to Verily

This ensures your source `.py` files stay as the authoritative version.

---

## Core Technologies and Patterns

### Data Access

**BigQuery to Polars** (preferred):
```python
def polars_gbq(query: str) -> pl.DataFrame:
    """Execute BigQuery SQL and return polars DataFrame"""
    client = bigquery.Client()
    query_job = client.query(query)
    rows = query_job.result()
    return pl.from_arrow(rows.to_arrow())
```
- Faster than pandas for large datasets
- Better memory efficiency
- Use this pattern consistently

**Environment Variables** (required):
```python
WORKSPACE_CDR = os.environ['WORKSPACE_CDR']  # BigQuery dataset
WORKSPACE_BUCKET = os.environ['WORKSPACE_BUCKET']  # GCS bucket
```
Never hardcode workspace paths - always use environment variables.

### OMOP CDM Querying

**Flexible Concept Searching** - Use `AouQueries` class from `sarcoid/01 Sarcoid Cohort.py`:
```python
aou = AouQueries(version='v8', bucket=WORKSPACE_BUCKET)

# Text-based search across vocabularies
aou.find_diagnosis_codes(
    vocabularies=['ICD9CM', 'ICD10CM', 'SNOMED'],
    search_terms=['sarcoid'],
    exact_codes={'ICD10CM': ['D86.0', 'D86.9']},
    pattern_codes={'ICD10CM': ['D86.%']},
    exclude_codes={'ICD10CM': ['D86.1']},
    exclude_terms=['sarcoidosis cutis']
)

# Create person-level summary
cases_df = aou.person_code_df(count_threshold=2)  # ≥2 codes = case
```

**Special ICD V-Code Handling** - ICD9CM V codes require `concept_relationship` joins (see `trusted_queries/icd_codes.py`).

### Distributed Computing (GWAS)

**dsub** for parallelized jobs on Google Cloud Batch:
```python
def dsub_script(machine_type, envs, in_params, out_params, script, ...):
    """Submit job to Google Cloud Batch

    Args:
        machine_type: e.g., 'c4-standard-8', 'n2d-standard-8'
        envs: dict of environment variables
        in_params: dict of input file paths (GCS)
        out_params: dict of output file paths (GCS)
        script: bash script to execute
    """
```

**SAIGE GWAS Workflow** (see `hpv/B02*.py` or `sarcoid/02*.py`):
1. Fit null GLMM model (`run_saige_null()`) - single job
2. Run chromosome-wise tests (`run_saige_step2()`) - 22 parallel jobs
3. Check status with `check_dsub_status()`

**METAL Meta-Analysis** (see `hpv/B03*.py` or `sarcoid/03*.py`):
```python
run_metal(trait='condition__sarcoid', ancestries=['eur', 'afr'], base_output_folder='saige_gwas/min_1')
```

### Standard Imports and Setup

All notebooks use consistent setup (see `hpv/B01.a HPV Cohort v0.py:1-50`):
```python
from google.cloud import bigquery
import pandas as pd
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime, date
from jinja2 import Template
import os

# Environment
WORKSPACE_CDR = os.environ['WORKSPACE_CDR']
WORKSPACE_BUCKET = os.environ['WORKSPACE_BUCKET']

# Plotting style
sns.set(style="whitegrid", font_scale=0.9)
palette = ['#0173b2', '#de8f05', '#8de5a1', '#d55e00', '#029e73',
           '#cc78bc', '#ece133', '#56b4e9', '#ca9161', '#fbafe4', '#949494']
sns.set_palette(sns.color_palette(palette))
```

---

## Key Reusable Functions

### From HPV/Sarcoid Projects

**Cohort Definition**:
- `polars_gbq(query)` - BigQuery → polars DataFrame (used everywhere)
- `AouQueries` class - Flexible OMOP concept searching (sarcoid)
- `phetk_icd_query(ds)` - ICD code extraction with V-code handling (trusted_queries)
- `update_covariates(df)` - Demographic covariate engineering (hpv)
- `describe_group(df)` - Table 1 descriptive statistics with <20 suppression (hpv)

**GWAS Pipeline**:
- `dsub_script()` - Submit distributed computing jobs (hpv/sarcoid)
- `run_saige_null()` - Fit SAIGE null model (hpv/sarcoid)
- `run_saige_step2()` - Run SAIGE chromosome tests (hpv/sarcoid)
- `check_dsub_status()` - Monitor job status (hpv/sarcoid)
- `combine_gwas_results()` - Merge per-chromosome results (hpv/sarcoid)
- `validate_saige_inputs()` - Check SAIGE output completeness (hpv/sarcoid)
- `run_metal()` - METAL meta-analysis wrapper (hpv/sarcoid)

**Post-Processing**:
- `convert_snp_format()` - SAIGE → gwaslab format (sarcoid)
- `add_chr_pos_from_snpid()` - Parse variant IDs (sarcoid)
- `process_gwas()` - Complete gwaslab pipeline (sarcoid)
- `split_gwas_file()` - Split large files for download (sarcoid)

**Utilities**:
- `parse_pcs(df)` - Extract principal components (sarcoid)
- `binarize_sex_chromosomes(df)` - Sex encoding from ploidy (sarcoid)
- `print_ancestry_summary(df, trait)` - Display case/control counts (sarcoid)

---

## All of Us Platform Specifics

### Data Policies (CRITICAL)

**Small Count Suppression**:
- Display `"<20"` for counts 1-19 in ALL outputs
- Applies to tables, plots, and text descriptions
- Never show exact counts below 20

**Genomic Data**:
- Filter to `has_whole_genome_variant = 1` before GWAS
- Check genomic_metrics.tsv for QC metrics

**Workspace Variables**:
- Always use `os.environ['WORKSPACE_CDR']` (never hardcode dataset)
- Always use `os.environ['WORKSPACE_BUCKET']` (never hardcode bucket)

### Reference Data Locations

Standard reference data in `{src_bucket}/data/`:
- `cohort_metadata__pcs__phenotypes.tsv.gz` - Demographics, PCs, phenotypes
- `ancestry_specific_pcs.parquet` - Ancestry-specific PCs
- `genomic_metrics.tsv` - Sex ploidy, QC metrics
- `{ancestry}_merged_genotypes.{bed,bim,fam}` - LD-pruned genotypes for null models
- `genotypes_chr{1-22}.{bed,bim,fam}` - Full genotypes for GWAS

### Ancestry Codes

Standard ancestry labels (from `ancestry_pred_other` column):
- `EUR` - European ancestry
- `AFR` - African ancestry
- `AMR` - Admixed American (Hispanic/Latino)
- `EAS` - East Asian ancestry
- `SAS` - South Asian ancestry

**IMPORTANT**: Use consistent casing (uppercase) throughout pipeline. Define once:
```python
ANCESTRIES = ['EUR', 'AFR', 'AMR']  # Define at top of file
```

### Machine Types and Costs

**Interactive Notebooks**:
- Light processing: n2-standard-16 (16 vCPU, 64GB RAM)
- Large cohorts: n2-standard-32 or n2-standard-64

**GWAS Jobs (dsub)**:
- SAIGE null model: c4-standard-8 or n2-highmem-16
- SAIGE chr tests: c4-standard-8 (22 parallel jobs, use preemptible)
- METAL meta-analysis: n2d-standard-8

**NOTE**: c4 machine types don't work with google-batch provider. Use c2 or n2 families instead.

See `_reference/resource_monitoring_template.md` for monitoring code and cost tracking.

---

## GWAS Pipeline Architecture

### File Naming Conventions

HPV uses B## prefix:
- `B01` - Cohort definition
- `B02` - SAIGE GWAS
- `B03` - METAL meta-analysis
- `B04` - Post-processing
- `B05` - Download
- `B06` - HLA visualization
- `B07` - HLA QC

Sarcoid uses ## prefix (01-05).

**Version/iteration suffixes**:
- `v0`, `v1`, `v2` - Major iterations
- Dates in parentheses - Last update (e.g., `(5.20)` = May 20)
- `redo` suffix - Rerun/validation

### SAIGE Configuration

**Standard covariates**:
```python
# Ancestry-specific
covariates = ['sex_binary', 'end_of_study_age'] + ['ancPC{}'.format(i) for i in range(1, 21)]
qCovarColList = ['sex_binary']  # Categorical (for sex chromosome)

# HPV uses slightly different names
covariates = ['imputed_sex', 'age_at_last_ehr'] + ['ancPC{}'.format(i) for i in range(1, 21)]
```

**Standard parameters**:
- `minMAC=20` - Minimum minor allele count
- `LOCO=TRUE` - Leave one chromosome out
- `traitType=binary` or `quantitative`
- Firth correction for rare variants

**Docker images**:
- SAIGE: `wzhou88/saige:1.3.6`
- METAL: Custom `bwaxse/metal` image
- Base environment: `us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.2.14`

### Statistical Thresholds

- **Genome-wide significance**: p < 5×10⁻⁸
- **Bonferroni**: p < 0.05 / N_tests (varies by analysis)
- **Minimum MAC**: 20
- **LD clumping window**: 500kb (for lead SNP identification)

---

## Common Workflows

### Starting a New GWAS Project

1. **Choose template**: HPV (multi-site) or Sarcoid (single condition)
2. **Read relevant CLAUDE.md** for full context
3. **Copy template files**:
   ```bash
   cp hpv/B01.a\ HPV\ Cohort\ v0.py my_project/01_cohort.py
   cp hpv/B02.a\ Anal\ HPV\ SAIGE\ GWAS\ v0.py my_project/02_gwas.py
   ```
4. **Modify cohort definition** (swap concept IDs/search terms)
5. **Keep infrastructure code** (dsub, SAIGE, METAL functions)
6. **Run pipeline**: 01 → 02 → 03 → 04 → 05

### Extracting OMOP Data

**Option 1: Hardcoded concept IDs** (HPV pattern):
```python
concept_ids = [45572221, 1567478, 35206462]  # HPV-related concepts
query = f"""
SELECT person_id, condition_start_date, condition_concept_id
FROM {WORKSPACE_CDR}.condition_occurrence
WHERE condition_concept_id IN ({','.join(map(str, concept_ids))})
"""
df = polars_gbq(query)
```

**Option 2: Flexible text search** (Sarcoid pattern):
```python
aou = AouQueries()
aou.find_diagnosis_codes(
    search_terms=['sarcoid'],
    vocabularies=['ICD9CM', 'ICD10CM', 'SNOMED']
)
df = aou.polars_gbq(aou.query_text)
cases = aou.person_code_df(count_threshold=2)
```

**Option 3: ICD codes with V-code handling**:
```python
from trusted_queries.icd_codes import phetk_icd_query
query = phetk_icd_query(WORKSPACE_CDR)
df = polars_gbq(f"SELECT * FROM ({query}) WHERE ICD LIKE 'D86%'")
```

### Running Distributed Jobs

```python
# 1. Submit job
script = """
#!/bin/bash
step1_fitNULLGLMM.R --plinkFile={PLINK} --phenoFile={PHENO} ...
"""
dsub_script(
    machine_type='c4-standard-8',
    envs={'ANCESTRY': 'eur'},
    in_params={'PLINK': f'{bucket}/genotypes', 'PHENO': f'{bucket}/metadata.tsv'},
    out_params={'OUTPUT': f'{bucket}/results/'},
    script=script,
    preemptible=True
)

# 2. Check status
check_dsub_status(user='bwaxse', age='3d')

# 3. Access results
results = pl.read_csv(f'{WORKSPACE_BUCKET}/results/gwas_results.tsv.gz', separator='\t')
```

---

## Claude Code Hooks (Safety Checks)

This repository has custom hooks in `.claude/` that will warn or block certain operations:

**Blocked**:
- `gsutil rm -r` (recursive GCS deletion) - Prevents accidental deletion of months of work

**Warnings**:
- Large dsub jobs (≥100 tasks) - Cost reminder
- Ancestry codes in GWAS files - Consistency check reminder
- `print()` in .py files - Consider using logging module
- Notebook exports - Reminder to sync .ipynb and .py

**Reminders**:
- On session end - Export notebooks to .py before finishing

These hooks are defined in `.claude/hookify.*.local.md` files. If a hook blocks you incorrectly, set `enabled: false` in the relevant file.

---

## Best Practices

### Code Organization

- Keep notebooks focused on single analysis stage
- Extract reusable functions to importable .py files
- Document function parameters and return types
- Use descriptive variable names (not `df1`, `df2`)

### Performance

- Use `polars_gbq()` instead of pandas for BigQuery (faster, less memory)
- Use preemptible instances for parallelizable jobs (60-80% cost savings)
- Store outputs in GCS, not local filesystem
- Monitor memory usage with `psutil` (see `_reference/resource_monitoring_template.md`)

### Data Quality

- Always validate OMOP concept IDs before using
- Check for duplicate person_ids after joins
- Verify ancestry consistency across pipeline stages
- Suppress counts <20 in all outputs
- Apply physiologically plausible ranges for measurements

### Version Control

- Commit both .ipynb and .py files together
- Export notebooks before committing: `jupyter nbconvert --to python *.ipynb`
- Use descriptive commit messages
- Tag major analysis versions

### Documentation

- Each project should have a `CLAUDE.md` file
- Document key functions with docstrings
- Note any deviations from standard patterns
- Track resource usage and costs

---

## Troubleshooting

### Common Issues

**Out of memory errors**:
- Upgrade machine type (n2-standard-32 → n2-standard-64)
- Use `polars` instead of `pandas`
- Process data in chunks
- Free memory between operations: `del df; gc.collect()`

**dsub job failures**:
- Check logs: `gsutil cat gs://{bucket}/dsub/logs/{job_id}.log`
- c4 machines don't work with google-batch (use c2/n2)
- Increase timeout for long-running jobs
- Non-preemptible for null model (can't be interrupted)

**ICD V-code issues**:
- Use `phetk_icd_query()` from `trusted_queries/icd_codes.py`
- V codes require `concept_relationship` table joins
- Don't query V codes directly from condition_occurrence

**Ancestry mismatches**:
- Define `ANCESTRIES = ['EUR', 'AFR']` once at top of file
- Use consistent casing (uppercase)
- Verify filtering in cohort, GWAS, and meta-analysis stages

**SAIGE convergence failures**:
- Check phenotype distribution (case/control balance)
- Try different covariate sets
- Increase memory for large sample sizes
- Check for missing covariates

---

## Dependencies

### Python Libraries
```
pandas>=1.3.0
polars>=0.18.0
numpy>=1.20.0
scipy>=1.7.0
matplotlib>=3.4.0
seaborn>=0.11.0
google-cloud-bigquery>=3.0.0
google-cloud-storage>=2.0.0
pyarrow>=8.0.0
jinja2>=3.0.0
gwaslab  # GWAS visualization
```

### External Tools (via Docker)
- SAIGE 1.3.6 - Binary/quantitative trait GWAS
- METAL - Inverse-variance weighted meta-analysis
- PLINK2 - Genotype file manipulation
- dsub - Google Cloud Batch job submission

### All of Us Platform
- Access to Researcher Workbench required
- Registered user tier for genomic data
- Appropriate Data Use Agreement

---

## Getting Help

- **All of Us Support**: [Researcher Help Desk](https://support.researchallofus.org/hc/en-us)
- **Project Documentation**: See individual `{project}/CLAUDE.md` files
- **Issues**: Open GitHub issue in this repository

---

**Last Updated**: January 2026
**Author**: Bennett Waxse
