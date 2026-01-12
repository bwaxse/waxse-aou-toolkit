# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Directory Overview

**Purpose**: Extract and harmonize medication data from All of Us OMOP drug_exposure table

**Status**: Production utility script
**Primary File**: `EXPORT 01.12_extraction_medication_ingredients_rx_norm.py` (exported from .ipynb)

This directory contains a comprehensive medication extraction pipeline that maps drug exposures to RxNorm ingredient-level concepts with sophisticated filtering and cleaning.

---

## Main File: EXPORT 01.12_extraction_medication_ingredients_rx_norm.py

### What It Does

Extracts **systemic outpatient medication exposures** from All of Us data with the following processing steps:

1. **Extraction**: Queries drug_exposure table with RxNorm ingredient mapping via cb_criteria_ancestor
2. **Route Filtering**: Removes topical/non-systemic routes (creams, inhalers, patches, etc.)
3. **Route Imputation**: Infers "oral" route from drug name when route is missing
4. **Setting Filtering**: Removes inpatient/emergency/ICU medications (focus on chronic outpatient use)
5. **Ingredient Mapping**: Maps to root-level RxNorm ingredient concepts (rank 1 in cb_criteria)
6. **Dose Extraction**: Parses strength in mg from drug names where available
7. **Export**: Produces person-level and grouped summaries

### Key Query Pattern

```python
is_descendants_q = f"""
    SELECT
        d_exposure.person_id,
        d_standard_concept.concept_name AS standard_concept_name,
        d_exposure.drug_exposure_start_date,
        d_exposure.drug_exposure_end_date,
        d_exposure.verbatim_end_date,
        d_exposure.refills,
        d_exposure.quantity,
        d_exposure.days_supply,
        d_route.concept_name AS route_concept_name,
        d_visit.concept_name AS visit_occurrence_concept_name,
        cr.concept_id AS cb_criteria_root_ancestor_id,
        cr_concept.concept_name AS cb_criteria_root_ancestor_name  -- Ingredient name
    FROM {CDR}.drug_exposure AS d_exposure
    LEFT JOIN {CDR}.concept AS d_standard_concept
        ON d_exposure.drug_concept_id = d_standard_concept.concept_id
    LEFT JOIN {CDR}.concept AS d_route
        ON d_exposure.route_concept_id = d_route.concept_id
    LEFT JOIN {CDR}.visit_occurrence AS v
        ON d_exposure.visit_occurrence_id = v.visit_occurrence_id
    LEFT JOIN {CDR}.concept AS d_visit
        ON v.visit_concept_id = d_visit.concept_id
    LEFT JOIN {CDR}.cb_criteria_ancestor AS ca
        ON d_exposure.drug_concept_id = ca.descendant_id
    LEFT JOIN (
        SELECT DISTINCT cr.concept_id, CAST(cr.id AS STRING) AS id
        FROM {CDR}.cb_criteria AS cr
        WHERE cr.full_text LIKE '%_rank1]%'  -- Only root-level (rank 1) ancestors
    ) cr
        ON ca.ancestor_id = cr.concept_id
    LEFT JOIN {CDR}.concept AS cr_concept
        ON cr.concept_id = cr_concept.concept_id;
"""
```

**Key Features**:
- Uses `cb_criteria_ancestor` table to map drugs to ingredient-level concepts
- Filters to `_rank1` ancestors (root ingredients, not drug forms or brands)
- Joins visit_occurrence for setting information (inpatient vs outpatient)
- Preserves multiple date fields for temporal analysis

---

## Filtering Strategy

### 1. Route Exclusions (Non-Systemic)

Removes medications that are not systemically absorbed:
```python
words_remove_rows_insensitive = [
    "topical", "spray", "lotion", "cream", "inhaler", "rectal", "suppository",
    "inhalation", "sponge", "patch", "metered", "ophthalmic", "ointment",
    "shampoo", "gel", "pad", "soap", "mucosal", "vaginal", "insert", "vaccine",
    "film", "meningococcal B", "virus", "hpv"
]

medication_df_filtered = medication_df.filter(
    ~pl.col("standard_concept_name").str.to_lowercase().str.contains(remove_pattern)
)
```

**Rationale**: For systemic exposure studies (GWAS, drug-disease associations), only systemically absorbed medications are relevant.

### 2. Route Imputation

When route is missing or coded as "No matching concept", infers oral route from drug name:
```python
oral_words = ["oral", "capsule", "softgel", "tablet", "chewable"]
oral_pattern = "|".join(map(re.escape, oral_words))

medication_df_updated = medication_df_filtered.with_columns(
    pl.when(pl.col("standard_concept_name").str.to_lowercase().str.contains(oral_pattern))
    .then(pl.lit("Oral route"))
    .otherwise(pl.col("route_concept_name"))
    .alias("route_concept_name")
)
```

**Rationale**: Many oral medications have missing route codes but the form (tablet, capsule) is in the name.

### 3. Retained Routes

After cleaning, keeps only systemic routes:
```python
route_keep_rows_insensitive = [
    "sublingual", "subcutaneous", "oral", "no matching concept", "NA", "Null",
    "intravenous", "intramuscular"
]
```

### 4. Setting Exclusions (Inpatient)

Removes medications from acute care settings:
```python
words_remove_rows_insensitive = ["emergency", "inpatient", "intensive care"]

medication_df_oral_outpatient = medication_df_oral.filter(
    pl.col("visit_occurrence_concept_name").is_null() |
    ~pl.col("visit_occurrence_concept_name").str.to_lowercase().str.contains(remove_pattern)
)
```

**Rationale**:
- Focus on chronic/maintenance medications (outpatient prescriptions)
- Avoid acute treatment medications that don't reflect long-term use
- Includes NULL settings (likely outpatient pharmacy fills)

---

## Output Files

### 1. meds_extract_grouped.csv
**Format**: Grouped by drug + ingredient + route
```
standard_concept_name | cb_criteria_root_ancestor_name | route_concept_name | count
```

**Use**:
- Quality check for ingredient mapping
- Identify most common drug forms per ingredient
- Verify route assignments

### 2. meds_extract.csv
**Format**: Unique person counts per ingredient
```
cb_criteria_root_ancestor_name | count
```

**Use**:
- Person-level exposure counts for each ingredient
- Identify prevalent medications in cohort
- Input for medication exposure GWAS or PheWAS

---

## Dose Extraction Feature

Extracts numeric strength from drug names:
```python
medication_df_doses = medication_df_oral_outpatient.with_columns(
    pl.col("standard_concept_name")
    .str.extract(r'(\d+(\.\d+)?)\s*MG', 1)  # Extract numbers before "MG"
    .cast(pl.Float64)
    .alias("med_strength_mg")
)
```

**Captures**:
- "Metformin 500 MG" → 500.0
- "Lisinopril 10 MG Oral Tablet" → 10.0
- "Atorvastatin 20.5 MG" → 20.5

**Limitations**:
- Only captures mg units (not mcg, IU, etc.)
- Misses combination products with multiple active ingredients
- Doesn't capture dosing frequency

---

## Date Handling Strategy

The script considers three date fields:
1. `drug_exposure_start_date` - Official start date
2. `drug_exposure_end_date` - Official end date
3. `verbatim_end_date` - Original end date from source system

**Processing Logic** (see lines 296-299):
- Use latest available date for end date
- If end date missing, use start date
- Verbatim dates take precedence when available

---

## Reusable Patterns

### RxNorm Ingredient Mapping via cb_criteria

The `cb_criteria` and `cb_criteria_ancestor` tables provide hierarchical drug classification:
- **_rank1**: Root ingredient level (e.g., "Metformin")
- **_rank2**: Drug form level (e.g., "Metformin 500 MG Oral Tablet")
- **_rank3+**: Specific branded products

**To map to ingredients**:
```sql
LEFT JOIN {CDR}.cb_criteria_ancestor AS ca
    ON d_exposure.drug_concept_id = ca.descendant_id
LEFT JOIN (
    SELECT DISTINCT cr.concept_id
    FROM {CDR}.cb_criteria AS cr
    WHERE cr.full_text LIKE '%_rank1]%'
) cr
    ON ca.ancestor_id = cr.concept_id
```

### Polars String Filtering Pattern

Flexible case-insensitive filtering with multiple terms:
```python
# Define exclusion terms
exclude_terms = ["topical", "cream", "ointment"]
pattern = "|".join(map(re.escape, exclude_terms))

# Filter (exclude matches)
df_filtered = df.filter(
    ~pl.col("column_name").str.to_lowercase().str.contains(pattern)
)

# Filter (keep matches)
df_filtered = df.filter(
    pl.col("column_name").str.to_lowercase().str.contains(pattern)
)
```

**Why `re.escape()`**: Escapes special regex characters (e.g., "." in "0.5 MG")

---

## Integration with Other Projects

### For GWAS Studies

Use ingredient-level exposures as binary traits:
```python
# Load medication exposures
med_exposures = pl.read_csv("meds_extract.csv")

# Filter to specific ingredient
metformin_users = med_exposures.filter(
    pl.col("cb_criteria_root_ancestor_name") == "Metformin"
)

# Create binary exposure variable for cohort
cohort = cohort.join(
    metformin_users.select(["person_id"]).with_columns(pl.lit(1).alias("exposed_to_metformin")),
    on="person_id",
    how="left"
).with_columns(
    pl.col("exposed_to_metformin").fill_null(0)
)
```

### For PheWAS Studies

Use as exposure variable to identify medication-associated conditions.

### For Temporal Analysis

Preserve date fields for:
- Time-to-event analyses
- Exposure windows relative to outcomes
- Treatment pattern analysis

---

## Data Quality Considerations

### Missing Ingredient Mappings

Not all drug_exposure records map to rank 1 ingredients:
```python
null_count = medication_df_oral_outpatient.filter(
    medication_df_oral_outpatient["cb_criteria_root_ancestor_name"].is_not_null()
).shape[0]
print(f"Number of rows with RxNorm ingredient mapping: {null_count}")
```

**Common reasons for NULL**:
- Non-standard drug concepts
- Missing cb_criteria hierarchy
- Devices/supplies coded as drugs
- Custom/compounded medications

### Route Missingness

Many drug records have `route_concept_id = 0` ("No matching concept"):
- Script attempts to impute from drug name
- Review `route_concept_name` value counts to assess quality

### Visit Occurrence Missingness

NULL visit_occurrence_id:
- Often indicates pharmacy fills (outpatient)
- Kept in final dataset (assumed outpatient)

---

## Common Use Cases

### 1. Extract Exposure to Specific Drug Class

```python
# Load data
med_df = polars_gbq(is_descendants_q)
med_filtered = apply_all_filters(med_df)  # Apply route/setting filters

# Filter to statins (example)
statins = ["Atorvastatin", "Simvastatin", "Rosuvastatin", "Pravastatin"]
statin_pattern = "|".join(statins)

statin_users = med_filtered.filter(
    pl.col("cb_criteria_root_ancestor_name").str.contains(statin_pattern)
).select(["person_id"]).unique()
```

### 2. Calculate Medication Exposure Duration

```python
# Using start and end dates
med_df_with_duration = med_filtered.with_columns(
    (pl.col("drug_exposure_end_date") - pl.col("drug_exposure_start_date")).alias("duration_days")
)

# Aggregate per person per ingredient
person_ingredient_exposure = med_df_with_duration.group_by(
    ["person_id", "cb_criteria_root_ancestor_name"]
).agg([
    pl.col("duration_days").sum().alias("total_exposure_days"),
    pl.col("drug_exposure_start_date").min().alias("first_exposure"),
    pl.col("drug_exposure_end_date").max().alias("last_exposure"),
    pl.len().alias("n_prescriptions")
])
```

### 3. Identify Polypharmacy

```python
# Count distinct ingredients per person
polypharmacy = med_filtered.group_by("person_id").agg(
    pl.col("cb_criteria_root_ancestor_name").n_unique().alias("n_medications")
)

# Flag high polypharmacy (≥5 medications)
polypharmacy = polypharmacy.with_columns(
    (pl.col("n_medications") >= 5).alias("high_polypharmacy")
)
```

---

## Limitations and Future Enhancements

### Current Limitations

1. **Dose units**: Only extracts mg, not mcg/IU/g/mL
2. **Combination products**: Single strength extracted (doesn't handle multiple active ingredients)
3. **Dosing frequency**: Not captured (would need days_supply and quantity)
4. **ATC codes**: Commented out in script (could add for therapeutic class grouping)
5. **Date logic**: End date calculation logic mentioned but incomplete in script

### Potential Enhancements

1. **ATC classification**: Enable ancestor_atc_code join for therapeutic groupings
2. **Dose normalization**: Convert all units to standard (e.g., mg equivalents)
3. **Defined Daily Dose (DDD)**: Calculate cumulative exposure in DDDs
4. **Treatment episodes**: Cluster prescriptions into treatment episodes using gaps
5. **Adherence metrics**: Calculate proportion of days covered (PDC) or medication possession ratio (MPR)

---

## Dependencies

```python
import os
import polars as pl
from google.cloud import bigquery
import re
```

**Required**:
- `polars>=0.18.0` - DataFrame library (faster than pandas)
- `google-cloud-bigquery>=3.0.0` - BigQuery client
- `re` - Regular expressions (standard library)

---

## Notes for Claude Code

### When Working with This File

1. **Filtering order matters**: Route filtering → Route imputation → Setting filtering
2. **Polars syntax**: Uses method chaining, not pandas-style operations
3. **Case sensitivity**: All string matching is case-insensitive via `.str.to_lowercase()`
4. **NULL handling**: Polars distinguishes between NULL and empty string

### Common Tasks

**Add new route exclusion**:
```python
words_remove_rows_insensitive.append("new_route_term")
```

**Add new ingredient of interest**:
```python
target_ingredients = ["Metformin", "Insulin"]
pattern = "|".join(target_ingredients)
target_df = med_filtered.filter(
    pl.col("cb_criteria_root_ancestor_name").str.contains(pattern)
)
```

**Extract different dose unit**:
```python
# For mcg instead of mg
.str.extract(r'(\d+(\.\d+)?)\s*MCG', 1)
```

### Troubleshooting

**Memory issues**:
- drug_exposure table is large (~20M rows)
- Filter early in pipeline
- Use polars (more efficient than pandas)

**Missing ingredients**:
- Check cb_criteria table has `_rank1` entries
- Verify drug_concept_id is standard concept
- Some drugs may not have ingredient-level mapping

**Route quality**:
- High proportion of "No matching concept" is expected
- Review standard_concept_name to verify imputation logic catches oral forms

---

**Last Updated**: January 2026
**Author**: Bennett Waxse
