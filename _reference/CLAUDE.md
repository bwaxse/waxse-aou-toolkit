# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Directory Overview

**Purpose**: Reference materials, documentation, and example code for All of Us Research Workbench development

**Status**: Reference library (read-only)

This directory contains supporting materials that are referenced by active projects but not directly executed.

---

## Directory Structure

```
_reference/
├── all_of_us_tables/      # OMOP table schemas and metadata
├── phecode/               # Phecode mapping files
├── verily/                # Verily Workbench platform examples
└── resource_monitoring_template.md  # Template for tracking compute resources
```

---

## Subdirectories

### all_of_us_tables/

**Purpose**: OMOP CDM table schemas, reference data, and metadata for All of Us CDR v8

**Contents**:
- `All of Us Controlled Tier Dataset v8 CDR Data Dictionary (C2024Q3R8) - OMOP-Compatible Tables.tsv`
  - Complete data dictionary for OMOP tables
  - Column names, data types, descriptions
  - All 108+ OMOP tables included

- `All of Us Controlled Tier Dataset v8 CDR Data Dictionary (C2024Q3R8) - Wearables.tsv`
  - Wearables-specific table schemas
  - Fitbit device data tables
  - Activity, heart rate, sleep table structures

- `concept_domains.tsv`
  - OMOP concept domains (Condition, Drug, Measurement, etc.)
  - Domain hierarchy and relationships
  - Used for understanding concept organization

- `condition_vocabularies.tsv`
  - Vocabulary IDs used in condition_occurrence table
  - ICD9CM, ICD10CM, SNOMED counts
  - Helps determine which vocabularies to search

- `measurement_concepts.tsv`
  - Common measurement concepts with codes
  - Lab tests, vitals, anthropometrics
  - Includes LOINC codes and units

- `table_row_counts.tsv`
  - Row counts for all OMOP tables in CDR v8
  - Useful for estimating query performance
  - Updated per data release

- `table_schemas.tsv`
  - Simplified table schemas
  - Primary keys, foreign keys, indices
  - Quick reference for joins

- `visit_concepts.tsv`
  - Visit type concepts (Inpatient, Outpatient, ER, etc.)
  - Visit occurrence concepts used in All of Us
  - Useful for filtering by care setting

- `vocabulary_structure.tsv`
  - OMOP vocabulary hierarchy
  - Standard vs classification vs source vocabularies
  - Vocabulary relationships

**Usage**:
```python
# Load table schema reference
schemas = pd.read_csv("_reference/all_of_us_tables/table_schemas.tsv", sep="\t")

# Check which vocabularies to search for conditions
vocabs = pd.read_csv("_reference/all_of_us_tables/condition_vocabularies.tsv", sep="\t")

# Find common measurement concepts
measurements = pd.read_csv("_reference/all_of_us_tables/measurement_concepts.tsv", sep="\t")
```

---

### phecode/

**Purpose**: Phecode (phenotype code) mapping files for converting ICD codes to research phenotypes

**Contents**:
- `phecodeX_info.csv`
  - Phecode definitions and metadata
  - Phecode labels, categories, clinical descriptions
  - Hierarchical structure (phecode → parent category)

- `phecodeX_unrolled_ICD_CM.csv`
  - ICD-CM to phecode mappings
  - Unrolled from hierarchical structure
  - One row per ICD code → phecode mapping
  - Includes exclude ranges for phecode case definitions

**What are Phecodes?**

Phecodes are clinically meaningful phenotype groupings derived from ICD billing codes:
- Aggregate related ICD codes (e.g., all hypertension codes → phecode 401)
- More suitable for genetic research than raw ICD codes
- Hierarchical structure enables phenome-wide analyses (PheWAS)
- Include exclusion criteria to define clean case/control groups

**Example Structure**:
```
phecodeX_info.csv:
phecode | phenotype | category | sex | excl_phecodes | excl_phenotypes
401.1   | Essential hypertension | circulatory system | Both | 401.2 | Secondary hypertension

phecodeX_unrolled_ICD_CM.csv:
ICD  | vocabulary | phecode | excl_phecode | excl_phenotype
I10  | ICD10CM    | 401.1   | 401.2        | Secondary hypertension
I11.9| ICD10CM    | 401.1   | 401.2        | Secondary hypertension
```

**Usage**:
```python
# Load phecode mappings
phecode_info = pd.read_csv("_reference/phecode/phecodeX_info.csv")
phecode_map = pd.read_csv("_reference/phecode/phecodeX_unrolled_ICD_CM.csv")

# Map ICD codes to phecodes
icd_df = # ... your ICD codes
icd_with_phecode = icd_df.merge(
    phecode_map[['ICD', 'vocabulary', 'phecode']],
    on=['ICD', 'vocabulary'],
    how='left'
)

# Get phecode descriptions
phecode_labeled = icd_with_phecode.merge(
    phecode_info[['phecode', 'phenotype', 'category']],
    on='phecode',
    how='left'
)
```

**Integration with Toolkit**:
- Use with `trusted_queries/icd_codes.py` to map ICD extractions to phecodes
- Useful for PheWAS studies (test all phecodes vs genetic variant)
- Exclusion criteria can be applied for case definition

---

### verily/

**Purpose**: Example notebooks demonstrating Verily Workbench platform features

**Contents**:
- `verily_01_1. Setting_Env_Variables.ipynb`
  - How to set and access workspace environment variables
  - WORKSPACE_CDR, WORKSPACE_BUCKET, GOOGLE_PROJECT
  - Best practices for environment configuration

- `verily_01_2. Setting_Env_Variables_p2.ipynb`
  - Additional environment variable examples
  - Custom variables for project-specific paths

- `verily_02_1. Accessing All of Us Data.ipynb`
  - Connecting to BigQuery
  - Querying OMOP tables
  - Loading data to pandas/polars DataFrames
  - Best practices for data access

**Usage**:
These notebooks are **reference examples** showing platform patterns. Copy patterns into your analysis notebooks.

**Key Patterns from These Notebooks**:

```python
# Environment variables
import os
WORKSPACE_CDR = os.environ['WORKSPACE_CDR']
WORKSPACE_BUCKET = os.environ['WORKSPACE_BUCKET']
GOOGLE_PROJECT = os.environ['GOOGLE_PROJECT']

# BigQuery client
from google.cloud import bigquery
client = bigquery.Client(project=GOOGLE_PROJECT)

# Query pattern
query = f"""
SELECT * FROM {WORKSPACE_CDR}.person
LIMIT 10
"""
result = client.query(query).to_dataframe()

# GCS file access
from google.cloud import storage
storage_client = storage.Client()
bucket = storage_client.bucket(WORKSPACE_BUCKET.replace('gs://', ''))
blob = bucket.blob('path/to/file.csv')
blob.download_to_filename('local_file.csv')
```

---

## resource_monitoring_template.md

**Purpose**: Template for tracking compute resource usage and costs during analysis projects

**Contents**:
- Python code for monitoring memory and CPU usage
- Template tables for recording machine types and performance
- Cost estimation tables
- Lessons learned sections
- Machine type decision records

**When to Use**:
Copy this template into project-specific CLAUDE.md files to track:
- Which machine types were used for which operations
- Peak memory/CPU utilization
- Operation durations
- Estimated costs
- Performance bottlenecks
- Recommendations for future runs

**Key Monitoring Function**:
```python
import psutil
import os
from datetime import datetime

def print_resource_usage(label=""):
    """Monitor memory and CPU usage with recommendations"""
    memory = psutil.virtual_memory()
    cpu_percent = psutil.cpu_percent(interval=1)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    print(f"\nResource Usage: {label}")
    print(f"Timestamp: {timestamp}")
    print(f"Memory: {memory.used / 1e9:.1f}GB / {memory.total / 1e9:.1f}GB "
          f"({memory.percent:.1f}% used)")
    print(f"CPU: {cpu_percent:.1f}% used across {os.cpu_count()} cores")

    if memory.percent > 85:
        print("⚠️  HIGH MEMORY - Consider upgrading machine type")
    elif memory.percent < 40 and cpu_percent < 40:
        print("💰 LOW USAGE - Could downsize to save costs")
    else:
        print("✅ Good utilization")
```

**Usage in Notebooks**:
```python
# Add to notebook after loading large datasets
df = polars_gbq(query)
print_resource_usage("After loading diagnosis data")

# Track expensive operations
merged = df1.join(df2, on='person_id')
print_resource_usage("After join operation")

# Before/after comparisons
print_resource_usage("Before GWAS")
run_saige_null()
print_resource_usage("After SAIGE null model")
```

**Benefits**:
- Identify memory bottlenecks
- Right-size machine types (avoid over/under-provisioning)
- Estimate costs for future similar projects
- Document performance for reproducibility

---

## Usage Patterns

### Quick Reference Lookup

```python
# Find measurement concept for BMI
measurements = pd.read_csv("_reference/all_of_us_tables/measurement_concepts.tsv", sep="\t")
bmi_concepts = measurements[measurements['concept_name'].str.contains('BMI', case=False)]

# Check table row counts before query
row_counts = pd.read_csv("_reference/all_of_us_tables/table_row_counts.tsv", sep="\t")
condition_rows = row_counts[row_counts['table_name'] == 'condition_occurrence']['row_count'].values[0]
print(f"Querying {condition_rows:,} condition records")

# Map ICD to phecode
phecode_map = pd.read_csv("_reference/phecode/phecodeX_unrolled_ICD_CM.csv")
d86_phecode = phecode_map[phecode_map['ICD'].str.startswith('D86')]
```

### Project Setup

When starting a new project:
1. Copy `resource_monitoring_template.md` into project CLAUDE.md
2. Reference `all_of_us_tables/` for OMOP table schemas
3. Use `verily/` notebooks as examples for data access patterns
4. Check `phecode/` if doing PheWAS analysis

---

## Integration with Active Projects

### HPV/Sarcoid GWAS Projects

- Reference table schemas to understand OMOP joins
- Use phecode mappings for PheWAS analysis
- Apply resource monitoring during compute-intensive steps
- Check measurement concepts for covariates (BMI, labs)

### Medication Extraction

- Reference `concept_domains.tsv` to understand Drug domain
- Use `vocabulary_structure.tsv` to verify RxNorm hierarchies
- Check `table_row_counts.tsv` to estimate drug_exposure query time

### Wearables Analysis

- Reference `All of Us...Wearables.tsv` for Fitbit table schemas
- Use `visit_concepts.tsv` to understand activity contexts
- Apply resource monitoring for large wearables datasets

---

## Maintenance

### Updating Reference Files

When All of Us releases new CDR versions:

1. **Update data dictionaries**:
   - Download new data dictionary from All of Us
   - Replace files in `all_of_us_tables/`
   - Note CDR version in filename

2. **Update row counts**:
   ```sql
   -- Run in BigQuery
   SELECT table_name, row_count
   FROM `{WORKSPACE_CDR}.__TABLES__`
   ORDER BY table_name
   ```
   - Save as `table_row_counts.tsv`

3. **Update concept references**:
   ```sql
   -- Common measurement concepts
   SELECT concept_id, concept_name, vocabulary_id, concept_code
   FROM `{WORKSPACE_CDR}.concept`
   WHERE domain_id = 'Measurement'
     AND standard_concept = 'S'
   ORDER BY concept_name
   ```
   - Save as `measurement_concepts.tsv`

4. **Update phecode mappings**:
   - Check for new phecode releases at phewascatalog.org
   - Update `phecodeX_info.csv` and `phecodeX_unrolled_ICD_CM.csv`

---

## Notes for Claude Code

### When Working with Reference Files

1. **These are read-only**: Don't modify reference files during analysis
2. **Copy patterns, don't hardcode**: Extract patterns but don't hardcode file paths
3. **Check versions**: Ensure reference files match your CDR version
4. **Document usage**: If using reference data, note which files in analysis docs

### Common Tasks

**Find OMOP table schema**:
```python
schemas = pd.read_csv("_reference/all_of_us_tables/table_schemas.tsv", sep="\t")
person_schema = schemas[schemas['table_name'] == 'person']
```

**Estimate query performance**:
```python
row_counts = pd.read_csv("_reference/all_of_us_tables/table_row_counts.tsv", sep="\t")
drug_rows = row_counts[row_counts['table_name'] == 'drug_exposure']['row_count'].values[0]
print(f"Expect ~{drug_rows / 1e6:.0f}M row scan")
```

**Map to phecodes**:
```python
phecode_map = pd.read_csv("_reference/phecode/phecodeX_unrolled_ICD_CM.csv")
icd_df = icd_df.merge(phecode_map, on=['ICD', 'vocabulary'], how='left')
```

### Troubleshooting

**File not found**:
- Use absolute paths: `os.path.join(repo_root, '_reference', 'phecode', 'file.csv')`
- Or copy to workspace bucket for cloud access

**Outdated reference data**:
- Check CDR version in filename vs your WORKSPACE_CDR
- Update reference files if CDR version mismatch

**Large file load times**:
- Reference files are small (<100MB) and load quickly
- If slow, check network connection to workspace

---

**Last Updated**: January 2026
**Maintainer**: Bennett Waxse
**CDR Version**: C2024Q3R8 (v8)
