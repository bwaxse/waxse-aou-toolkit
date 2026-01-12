# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Directory Overview

**Purpose**: Reusable, production-tested query functions for common All of Us data extraction tasks

**Status**: Utility library
**Primary File**: `icd_codes.py`

This directory contains trusted, well-tested query functions that can be imported and reused across projects. These functions handle complex OMOP CDM patterns correctly.

---

## Files

### icd_codes.py

**Purpose**: Extract ICD-9-CM and ICD-10-CM diagnosis codes from OMOP condition_occurrence and observation tables with special V-code handling

**Key Function**: `phetk_icd_query(ds)`

---

## Main Function: phetk_icd_query()

### Function Signature

```python
def phetk_icd_query(ds: str) -> str:
    """
    Generate SQL query to extract ICD codes from All of Us OMOP data

    This query is optimized for the All of Us platform with special handling
    for ICD9CM V codes, which require concept_relationship table joins for
    correct vocabulary attribution.

    Args:
        ds (str): Google BigQuery dataset ID (e.g., from os.environ['WORKSPACE_CDR'])

    Returns:
        str: SQL query that generates a table with:
            - person_id: Participant identifier
            - date: Diagnosis date (condition_start_date or observation_date)
            - ICD: ICD code (concept_code)
            - vocabulary_id: ICD9CM or ICD10CM
    """
```

### What It Does

Extracts all ICD diagnosis codes from **four different paths** in OMOP:

1. `condition_occurrence` via `condition_source_value` (direct code match)
2. `condition_occurrence` via `condition_source_concept_id` (concept ID match)
3. `observation` via `observation_source_value` (direct code match)
4. `observation` via `observation_source_concept_id` (concept ID match)

**Special V-Code Handling**: ICD9CM V codes (e.g., V42.0 for kidney transplant) overlap between ICD9CM and ICD10CM vocabularies. The function uses `concept_relationship` table to ensure correct vocabulary attribution.

### Why This Query Exists

**Problem**: ICD codes in OMOP can appear in multiple tables and fields:
- `condition_occurrence` - Primary diagnosis table
- `observation` - Sometimes used for diagnoses in source data
- `source_value` - String code directly
- `source_concept_id` - Mapped concept ID

Additionally, V codes are ambiguous and require special handling to determine if they're ICD9CM or ICD10CM.

**Solution**: This query:
1. Searches all four paths
2. Handles V codes specially via concept_relationship
3. Returns UNION DISTINCT of all results
4. Ensures vocabulary_id is correct

---

## Query Structure

### Main Components

#### 1. icd_query (Base Extraction)

Extracts from 4 paths with UNION DISTINCT:

```sql
-- Path 1: condition_occurrence via source_value
SELECT DISTINCT
    co.person_id,
    co.condition_start_date AS date,
    c.vocabulary_id,
    c.concept_code AS ICD,
    co.condition_concept_id AS concept_id
FROM condition_occurrence AS co
INNER JOIN concept AS c
    ON co.condition_source_value = c.concept_code
WHERE c.vocabulary_id IN ("ICD9CM", "ICD10CM")

UNION DISTINCT

-- Path 2: condition_occurrence via source_concept_id
SELECT DISTINCT
    co.person_id,
    co.condition_start_date AS date,
    c.vocabulary_id,
    c.concept_code AS ICD,
    co.condition_concept_id AS concept_id
FROM condition_occurrence AS co
INNER JOIN concept AS c
    ON co.condition_source_concept_id = c.concept_id
WHERE c.vocabulary_id IN ("ICD9CM", "ICD10CM")

UNION DISTINCT

-- Path 3: observation via source_value
-- Path 4: observation via source_concept_id
-- (Similar structure for observation table)
```

**Key Point**: Uses both `source_value` (string) and `source_concept_id` (ID) joins to maximize capture.

#### 2. v_icd_vocab_query (V-Code Correction)

Fixes vocabulary attribution for V codes:

```sql
SELECT DISTINCT
    v_icds.person_id,
    v_icds.date,
    v_icds.ICD,
    c.vocabulary_id  -- Correct vocabulary from concept_relationship
FROM (
    SELECT * FROM icd_query
    WHERE ICD LIKE "V%"
) AS v_icds
INNER JOIN concept_relationship AS cr
    ON v_icds.concept_id = cr.concept_id_1
INNER JOIN concept AS c
    ON cr.concept_id_2 = c.concept_id
WHERE
    c.vocabulary_id IN ("ICD9CM", "ICD10CM")
    AND v_icds.ICD = c.concept_code
    AND NOT v_icds.vocabulary_id != c.vocabulary_id
```

**Purpose**:
- V codes appear in both ICD9CM and ICD10CM
- `concept_relationship` table provides correct mapping
- This sub-query replaces incorrect vocabulary_id with correct one

#### 3. final_query (Combine Results)

Merges non-V codes (correct as-is) with corrected V codes:

```sql
-- Non-V codes (vocabulary_id already correct)
SELECT DISTINCT person_id, date, ICD, vocabulary_id
FROM icd_query
WHERE NOT ICD LIKE "V%"

UNION DISTINCT

-- V codes (vocabulary_id corrected)
SELECT DISTINCT *
FROM v_icd_vocab_query
```

---

## Usage Examples

### Basic Usage

```python
from trusted_queries.icd_codes import phetk_icd_query
import os
from google.cloud import bigquery

# Get workspace CDR
WORKSPACE_CDR = os.environ['WORKSPACE_CDR']

# Generate query
query = phetk_icd_query(WORKSPACE_CDR)

# Execute query
client = bigquery.Client()
result = client.query(f"SELECT * FROM ({query})").to_dataframe()
```

### Filter to Specific ICD Codes

```python
# Sarcoidosis (D86.x)
sarcoid_query = f"""
SELECT DISTINCT person_id, date, ICD, vocabulary_id
FROM ({phetk_icd_query(WORKSPACE_CDR)})
WHERE ICD LIKE 'D86%'
"""

sarcoid_df = client.query(sarcoid_query).to_dataframe()
```

### Person-Level Case Counting

```python
# Count diagnoses per person
query = f"""
SELECT
    person_id,
    COUNT(DISTINCT date) AS n_diagnoses,
    MIN(date) AS first_diagnosis,
    MAX(date) AS last_diagnosis
FROM ({phetk_icd_query(WORKSPACE_CDR)})
WHERE ICD LIKE 'D86%'
GROUP BY person_id
"""

case_counts = client.query(query).to_dataframe()
```

### V-Code Examples

```python
# Extract V codes (these will have correct vocabulary)
v_code_query = f"""
SELECT * FROM ({phetk_icd_query(WORKSPACE_CDR)})
WHERE ICD LIKE 'V%'
ORDER BY person_id, date
"""

v_codes = client.query(v_code_query).to_dataframe()

# Example V codes:
# - V42.0: Kidney replaced by transplant (ICD9CM)
# - V58.69: Long-term use of other medications (ICD9CM)
# - V10.3: Personal history of malignant neoplasm (ICD9CM)
```

---

## Integration with Other Toolkit Components

### With AouQueries Class (Sarcoid Project)

The `AouQueries` class in `sarcoid/01 Sarcoid Cohort.py` implements similar logic but with more flexibility:
- Text-based searching
- Pattern matching
- Exclusion criteria
- Procedure codes

**When to use which**:
- **phetk_icd_query()**: When you need all ICD codes quickly
- **AouQueries**: When you need flexible searching with text/patterns/exclusions

### With PheWAS Analysis

This query was originally developed for PheWAS Toolkit (PheTK) integration:

```python
from trusted_queries.icd_codes import phetk_icd_query

# Get all ICD codes
icd_df = polars_gbq(f"SELECT * FROM ({phetk_icd_query(WORKSPACE_CDR)})")

# Map to phecodes (using phecode mappings from _reference/phecode/)
# Run PheWAS analysis
```

---

## V-Code Special Handling: Why It's Necessary

### The Problem

ICD9CM V codes (V01-V91) overlap with ICD10CM in the OMOP concept table:

| ICD Code | Could Be ICD9CM | Could Be ICD10CM |
|----------|-----------------|------------------|
| V42.0 | Kidney replaced by transplant | — |
| V58.69 | Long-term use of medications | — |
| V10.3 | Personal history of cancer | — |

Without special handling, the `vocabulary_id` field from direct joins may be incorrect.

### The Solution

Use `concept_relationship` table to find the **related standard concept**, which has the correct vocabulary:

```
source_concept_id → concept_relationship → standard_concept_id → vocabulary_id
```

This two-hop approach ensures V codes get the correct ICD9CM/ICD10CM attribution.

### Example

```sql
-- Without V-code handling (WRONG for some V codes)
SELECT c.vocabulary_id, c.concept_code
FROM condition_occurrence co
JOIN concept c ON co.condition_source_value = c.concept_code
WHERE c.concept_code = 'V42.0'
-- May return ICD10CM incorrectly

-- With V-code handling (CORRECT)
SELECT c2.vocabulary_id, c1.concept_code
FROM condition_occurrence co
JOIN concept c1 ON co.condition_source_concept_id = c1.concept_id
JOIN concept_relationship cr ON c1.concept_id = cr.concept_id_1
JOIN concept c2 ON cr.concept_id_2 = c2.concept_id
WHERE c1.concept_code = 'V42.0'
  AND c1.concept_code = c2.concept_code
-- Returns ICD9CM correctly
```

---

## Performance Considerations

### Query Complexity

This query is **complex** with:
- 4 base queries UNION DISTINCT
- Subquery for V code correction
- Multiple joins

**Expected runtime**: 30-60 seconds on full All of Us dataset

### Optimization Tips

1. **Filter early**: Add WHERE clauses to the subquery:
   ```python
   query = f"""
   SELECT * FROM ({phetk_icd_query(WORKSPACE_CDR)})
   WHERE ICD LIKE 'D86%'  -- Filter after query runs
   """
   ```

2. **Materialize results**: For repeated use, save to a table:
   ```python
   query = f"""
   CREATE OR REPLACE TABLE `{project}.{dataset}.all_icd_codes` AS
   SELECT * FROM ({phetk_icd_query(WORKSPACE_CDR)})
   """
   ```

3. **Use specific person_ids**: If you already have a cohort:
   ```python
   query = f"""
   SELECT * FROM ({phetk_icd_query(WORKSPACE_CDR)})
   WHERE person_id IN UNNEST(@person_ids)
   """
   ```

---

## Comparison to AouQueries

| Feature | phetk_icd_query() | AouQueries.find_diagnosis_codes() |
|---------|-------------------|-----------------------------------|
| **ICD extraction** | ✅ All ICD codes | ✅ All ICD codes + SNOMED |
| **V-code handling** | ✅ Automatic | ✅ Automatic |
| **Text search** | ❌ | ✅ (e.g., "sarcoid") |
| **Pattern matching** | ❌ | ✅ (e.g., "D86.%") |
| **Exact codes** | ❌ | ✅ (e.g., ['D86.0']) |
| **Exclusion terms** | ❌ | ✅ (e.g., exclude "cutis") |
| **Return format** | Raw SQL string | SQL string with parameters |
| **Complexity** | Simple function | Class with methods |

**Use phetk_icd_query() when**: You need all ICD codes quickly and simply

**Use AouQueries when**: You need flexible searching with text/patterns/exclusions

---

## Common Modifications

### Add ICD-9-PCS Procedure Codes

```python
def phetk_icd_pcs_query(ds):
    """Extract ICD procedure codes"""
    # Similar structure but use procedure_occurrence table
    # and filter to ICD9Proc vocabulary
```

### Add Date Range Filter

```python
def phetk_icd_query_filtered(ds, start_date, end_date):
    """Extract ICD codes within date range"""
    base_query = phetk_icd_query(ds)
    return f"""
    SELECT * FROM ({base_query})
    WHERE date BETWEEN '{start_date}' AND '{end_date}'
    """
```

### Extract Specific Vocabulary Only

```python
def phetk_icd10_query(ds):
    """Extract only ICD10CM codes"""
    base_query = phetk_icd_query(ds)
    return f"""
    SELECT * FROM ({base_query})
    WHERE vocabulary_id = 'ICD10CM'
    """
```

---

## Limitations and Future Enhancements

### Current Limitations

1. **ICD only**: Doesn't extract SNOMED, LOINC, or other vocabularies
2. **No hierarchy**: Doesn't use concept_ancestor for parent/child relationships
3. **No filtering**: Returns all ICD codes (user must filter)
4. **Condition context**: Doesn't include condition_type_concept_id (primary vs secondary diagnosis)

### Potential Enhancements

1. **Add SNOMED**: Expand to include SNOMED-CT diagnosis codes
2. **Add condition type**: Include primary/secondary/admitting diagnosis flags
3. **Add hierarchy**: Use concept_ancestor for "any sarcoidosis" queries
4. **Parameterize**: Add optional parameters for vocabularies, date ranges
5. **Optimize**: Create materialized view for frequently accessed data

---

## Dependencies

```python
# None - pure SQL generation
```

This function only generates SQL strings, so it has no Python dependencies beyond the standard library.

---

## Notes for Claude Code

### When Working with This File

1. **Don't modify the V-code logic**: It's carefully designed and tested
2. **Test modifications**: V-code handling is complex, test any changes thoroughly
3. **Document changes**: If adding features, document the new SQL patterns

### Common Tasks

**Create variant for different vocabulary**:
```python
# Copy phetk_icd_query and modify WHERE clauses
def phetk_snomed_query(ds):
    # Change: WHERE c.vocabulary_id IN ("ICD9CM", "ICD10CM")
    # To: WHERE c.vocabulary_id = "SNOMED"
```

**Add person_id filter parameter**:
```python
def phetk_icd_query(ds, person_ids=None):
    base_query = # ... existing logic
    if person_ids:
        return f"""
        SELECT * FROM ({base_query})
        WHERE person_id IN UNNEST(@person_ids)
        """
    return base_query
```

### Troubleshooting

**Slow performance**:
- This is expected for full dataset scan
- Filter results after query completion
- Consider materializing to temp table

**Missing diagnoses**:
- Check if codes exist in concept table
- Verify vocabulary_id is ICD9CM or ICD10CM
- Some source systems use local codes not in OMOP

**V-code issues**:
- V-codes are ICD9CM only (ICD-10 uses Z codes)
- If V-codes not appearing, check concept_relationship table exists

---

**Last Updated**: January 2026
**Author**: Bennett Waxse (PheTK integration)
