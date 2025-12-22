# BMI Harmonization for All of Us Research

**Purpose**: Extract, clean, and harmonize BMI (Body Mass Index) data with advanced quality control and temporal matching for All of Us cohort studies

**Platform**: All of Us Research Workbench (OMOP CDM v8)
**Analysis Type**: Data harmonization and temporal matching
**Primary Output**: Clean BMI measurements matched to cohort time points

---

## Project Overview

This project provides battle-tested code for BMI data extraction, harmonization, and temporal matching from All of Us OMOP data. It handles:
- Multiple measurement concept IDs (weight, height, BMI)
- Unit conversion and validation (pounds→kg, inches→cm)
- Outlier detection with unit-aware thresholds
- Height carry-forward for adult populations
- Sophisticated temporal matching strategies

The code is designed as **reusable functions** that can be imported into any All of Us analysis requiring BMI as a covariate or outcome.

---

## File Map

### bmi_functions.py
**Standalone function library** for BMI harmonization. All functions are self-contained and can be imported individually.

**Key exports:**
```python
# Complete pipeline
harmonize_bmi_complete(cohort_dfs, version, matching_strategy, time_col, ...)

# Individual components
extract_bmi_data(person_ids, version)
clean_units(df)
apply_clean_outliers(df, column, unit, sigma_threshold)
convert_units(df)
enhanced_height_carryforward(df, max_carryforward_days)
calculate_bmi_with_validation(df)
hierarchical_temporal_matching(cohort_df, bmi_df, time_col, ...)
analyze_temporal_matching_quality(merged_cohorts)
plot_temporal_matching_quality(merged_cohorts)

# Convenience wrappers for specific study types
bmi_for_gwas(cohort_dfs, version, time_col)
bmi_for_drug_study(cohort_dfs, version, time_col)
bmi_flexible(cohort_dfs, version, time_col)
```

**Function dependencies**: None (imports: pandas, polars, numpy, matplotlib, seaborn, google-cloud-bigquery)

### bmi_harmonization.ipynb
**Notebook demonstrating** the complete BMI harmonization workflow with:
- Step-by-step examples
- Quality control visualizations
- Multiple temporal matching strategies
- Usage patterns for different research scenarios

**Purpose**: Teaching/documentation - shows how to use `bmi_functions.py` in practice

---

## Workflow Overview

```
1. Extract BMI data from OMOP
   ↓
2. Clean units and remove invalid measurements
   ↓
3. Remove outliers (4-sigma, unit-aware)
   ↓
4. Convert to standard units (kg, cm)
   ↓
5. Enhanced height carry-forward (3 years for adults)
   ↓
6. Calculate/validate BMI
   ↓
7. Temporal matching with cohorts
   ↓
8. Quality assessment and visualization
```

---

## Key Functions

### 1. Complete Pipeline
```python
def harmonize_bmi_complete(
    cohort_dfs: List[pd.DataFrame],
    version: str,
    matching_strategy: str = 'standard',
    time_col: str = 'time_zero',
    height_carryforward_days: int = 1095,
    **matching_kwargs
) -> Tuple[pl.DataFrame, List[pd.DataFrame], Dict]:
    """
    End-to-end BMI harmonization with flexible temporal matching

    Args:
        cohort_dfs: List of DataFrames with person_id and time reference column
        version: All of Us CDR version (e.g., os.environ['WORKSPACE_CDR'])
        matching_strategy: 'standard', 'strict_prior', 'flexible', 'drug_study'
        time_col: Column name for reference time in cohort DataFrames
        height_carryforward_days: Days to carry forward height (default: 1095 = 3 years)

    Returns:
        bmi_df: Harmonized BMI measurements (polars DataFrame)
        merged_cohorts: List of cohort DataFrames with matched BMI
        quality_summary: Dict with matching quality metrics
    """
```

**Matching strategies:**
- **'standard'**: Most common use case
  - Prefer prior measurements (up to 365 days before time_zero)
  - Allow post measurements as backup (up to 90 days after)

- **'strict_prior'**: GWAS/genetic studies
  - Only prior measurements (no post-time_zero BMI)
  - Max 365 days before time_zero

- **'flexible'**: Longitudinal studies
  - Longer windows: 730 days prior, 180 days post

- **'drug_study'**: Baseline characteristics
  - Recent prior only: 90 days before, no post
  - Ensures BMI is measured close to drug initiation

---

### 2. Data Extraction

```python
def extract_bmi_data(
    person_ids: Union[List, pd.Series],
    version: str
) -> pl.DataFrame:
    """
    Extract BMI-related measurements from OMOP measurement table

    Handles:
    - Multiple input formats (list, Series, or list of DataFrames)
    - Weight concept IDs: 3010220, 3013762, 3025315, 3027492, 3023166
    - Height concept IDs: 3036798, 3023540, 3019171, 3036277, 40655804
    - BMI concept IDs: 3038553, 4245997

    Returns: polars DataFrame with columns:
        person_id, date, wt, wt_units, ht, ht_units, bmi, bmi_units
    """
```

**Validated concept IDs** (tested against All of Us data):
- **Weight**: Covers measurements in pounds, kilograms
- **Height**: Covers measurements in inches, centimeters
- **BMI**: Direct BMI recordings

---

### 3. Unit Cleaning

```python
def clean_units(df: pl.DataFrame) -> pl.DataFrame:
    """
    Remove measurements with problematic or ambiguous units

    Cleaning rules:
    - Weight: Remove 'No matching concept'
    - Height: Remove 'percent' (ambiguous)
    - BMI: Remove 'ratio', 'no value', 'No matching concept'
    """
```

---

### 4. Outlier Detection

```python
def apply_clean_outliers(
    df: pl.DataFrame,
    column: str,
    unit: Optional[str] = None,
    sigma_threshold: float = 4
) -> pl.DataFrame:
    """
    Remove outliers using sigma-based thresholds BY UNIT TYPE

    Key innovation: Calculate mean/SD separately for each unit
    - Prevents removing valid measurements just because they're in different units
    - Default 4-sigma threshold (very conservative)

    Example:
    - Weight in kg: mean=75, sd=15  →  bounds: [15, 135]
    - Weight in lbs: mean=165, sd=33 →  bounds: [33, 297]
    """
```

---

### 5. Unit Conversion

```python
def convert_units(df: pl.DataFrame) -> pl.DataFrame:
    """
    Standardize measurements to kg and cm

    Conversions:
    - inch (US) → cm: multiply by 2.54
    - pound (US) → kg: multiply by 0.45359237
    """
```

---

### 6. Enhanced Height Carry-Forward

```python
def enhanced_height_carryforward(
    df: pl.DataFrame,
    max_carryforward_days: int = 1095
) -> pl.DataFrame:
    """
    Carry forward height measurements for adults

    Rationale: Adult height changes slowly, so measurements within
    3 years can be used to fill gaps

    Adds 'height_source' column:
    - 'measured': Direct measurement
    - 'carried_forward': Filled from prior measurement
    - 'missing': No measurement available within time window

    Default: 1095 days (3 years) - appropriate for adults
    """
```

---

### 7. BMI Calculation & Validation

```python
def calculate_bmi_with_validation(df: pl.DataFrame) -> pl.DataFrame:
    """
    Calculate BMI and cross-validate with recorded BMI

    Process:
    1. Calculate: BMI = weight_kg / (height_m)²
    2. Compare calculated vs recorded BMI
    3. Use calculated BMI when recorded is missing
    4. Calculate difference for quality assessment

    Returns DataFrame with:
    - bmi: Final BMI (recorded or calculated)
    - bmi_calc: Calculated BMI
    - bmi_diff: Recorded - calculated (for QC)
    """
```

---

### 8. Hierarchical Temporal Matching

```python
def hierarchical_temporal_matching(
    cohort_df: pd.DataFrame,
    bmi_df: pd.DataFrame,
    time_col: str = 'time_zero',
    include_post: bool = True,
    max_prior_days: int = 365,
    max_post_days: int = 90
) -> pd.DataFrame:
    """
    Match BMI measurements to cohort with preference for PRIOR measurements

    Hierarchy (lower priority score = better):
    1. Prior measurements: priority = abs(days_diff)
       - Closer to time_zero is better
    2. Post measurements: priority = 1000 + abs(days_diff)
       - Heavily penalized, only used if no prior available

    Returns cohort DataFrame with:
    - bmi, wt_kg, ht_cm: Matched measurements
    - date: BMI measurement date
    - days_diff: Days from time_zero (negative = prior, positive = post)
    - bmi_timing: 'prior' or 'post'
    - bmi_quality: 'excellent' (<30 days), 'good' (<90 days),
                    'acceptable' (<180 days), 'poor' (>180 days)
    """
```

**Key advantage over simple "closest match"**: Ensures baseline BMI truly represents status *before* the event of interest, critical for causal inference.

---

### 9. Quality Assessment

```python
def analyze_temporal_matching_quality(
    merged_cohorts: List[pd.DataFrame]
) -> Dict:
    """
    Summarize matching quality across cohorts

    Returns dict with:
    - total_participants: N participants across all cohorts
    - participants_with_bmi: N with successful BMI match
    - matching_rate: Proportion with matched BMI
    - timing_distribution: Count of prior vs post matches
    - timing_percentages: Proportion prior vs post
    - quality_distribution: Count by quality category
    - quality_percentages: Proportion by quality category
    - days_diff_stats: Mean, median, std, min, max of temporal offset
    """

def plot_temporal_matching_quality(
    merged_cohorts: List[pd.DataFrame],
    figsize: Tuple[int, int] = (15, 5)
) -> matplotlib.figure.Figure:
    """
    Visualize matching quality

    Creates 3-panel figure:
    1. Histogram of days_diff (temporal offset distribution)
    2. Pie chart of prior vs post matches
    3. Bar chart of quality categories
    """
```

---

## Usage Examples

### Example 1: GWAS Study (Strict Prior Only)
```python
from bmi_functions import bmi_for_gwas

# Load your cohort
cohort_df = pd.read_csv('my_gwas_cohort.csv')
cohort_df['enrollment_date'] = pd.to_datetime(cohort_df['enrollment_date'])

# Get BMI data (strictly prior to enrollment)
bmi_data, [matched_cohort], quality = bmi_for_gwas(
    cohort_dfs=[cohort_df],
    version=os.environ['WORKSPACE_CDR'],
    time_col='enrollment_date'
)

# Check quality
print(f"Matching rate: {quality['matching_rate']:.1%}")
print(f"Prior matches: {quality['timing_percentages']['prior']:.1%}")

# Use matched_cohort for GWAS (has BMI column now)
```

### Example 2: Clinical Outcomes Study (Flexible)
```python
from bmi_functions import harmonize_bmi_complete, plot_temporal_matching_quality

# Multiple cohorts (cases and controls)
cases_df = pd.read_csv('cases.csv')
controls_df = pd.read_csv('controls.csv')

# Harmonize with flexible matching
bmi_data, [cases_matched, controls_matched], quality = harmonize_bmi_complete(
    cohort_dfs=[cases_df, controls_df],
    version=os.environ['WORKSPACE_CDR'],
    matching_strategy='flexible',  # Allow longer windows
    time_col='diagnosis_date'
)

# Visualize quality
plot_temporal_matching_quality([cases_matched, controls_matched])

# Export
cases_matched.to_csv('cases_with_bmi.csv', index=False)
controls_matched.to_csv('controls_with_bmi.csv', index=False)
```

### Example 3: Drug Response Study (Baseline BMI)
```python
from bmi_functions import bmi_for_drug_study

# Need BMI just before drug start (within 90 days)
drug_cohort_df = pd.read_csv('drug_cohort.csv')

bmi_data, [matched], quality = bmi_for_drug_study(
    cohort_dfs=[drug_cohort_df],
    version=os.environ['WORKSPACE_CDR'],
    time_col='drug_start_date'
)

# Check how many had recent BMI
print(f"Excellent quality (< 30 days): {quality['quality_percentages']['excellent']:.1%}")
```

### Example 4: Custom Matching Strategy
```python
from bmi_functions import harmonize_bmi_complete

# Custom: Prior only, but very long window (2 years)
bmi_data, merged_cohorts, quality = harmonize_bmi_complete(
    cohort_dfs=[my_cohort],
    version=os.environ['WORKSPACE_CDR'],
    matching_strategy='standard',  # Start with standard
    time_col='surgery_date',
    include_post=False,           # Override: no post-surgery BMI
    max_prior_days=730,          # Override: 2 years prior
    height_carryforward_days=1825  # Override: 5 years for height
)
```

---

## Architecture Patterns

### 1. Dual DataFrame Approach
- **Polars** for data extraction and transformation (faster, better BigQuery integration)
- **Pandas** for temporal matching and final output (easier datetime handling, scikit-learn compatibility)

```python
# Extract with polars (fast BigQuery → Arrow → polars)
bmi_df = polars_gbq(query)

# Clean and transform with polars (fast columnar operations)
bmi_df = clean_units(bmi_df)
bmi_df = convert_units(bmi_df)

# Convert to pandas for temporal matching (easier date arithmetic)
bmi_pandas = bmi_df.to_pandas()
merged = hierarchical_temporal_matching(cohort_df, bmi_pandas, ...)
```

### 2. Pipeline Composition
Functions are designed to be composable - you can use the complete pipeline OR individual components:

```python
# Option A: Complete pipeline
bmi_data, merged, quality = harmonize_bmi_complete(...)

# Option B: Custom pipeline with individual functions
raw_bmi = extract_bmi_data(person_ids, version)
clean_bmi = clean_units(raw_bmi)
outlier_free = apply_clean_outliers(clean_bmi, 'wt', 'wt_units')
outlier_free = apply_clean_outliers(outlier_free, 'ht', 'ht_units')
# ... continue as needed
```

### 3. Strategy Pattern for Temporal Matching
Pre-defined strategies encode best practices for common research scenarios:

```python
matching_params = {
    'standard': {
        'include_post': True,
        'max_prior_days': 365,
        'max_post_days': 90
    },
    'strict_prior': {...},
    'flexible': {...},
    'drug_study': {...}
}

# User selects strategy, function applies appropriate params
params = matching_params[matching_strategy]
params.update(matching_kwargs)  # Allow overrides
```

---

## Quality Control Features

### 1. Unit-Aware Outlier Detection
**Problem**: Naive outlier detection flags valid measurements in minority units.

**Solution**: Calculate statistics separately for each unit:
```python
# Bad approach (naive)
mean_weight = df['weight'].mean()  # Mixes kg and lbs!
sd_weight = df['weight'].std()

# Good approach (unit-aware)
for unit in df['weight_units'].unique():
    subset = df.filter(pl.col('weight_units') == unit)
    mean = subset['weight'].mean()
    sd = subset['weight'].std()
    # Apply thresholds within unit group
```

### 2. BMI Cross-Validation
Compare calculated BMI (from height/weight) with recorded BMI to detect:
- Data entry errors
- Unit conversion issues
- Implausible measurements

```python
df['bmi_diff'] = df['bmi_recorded'] - df['bmi_calculated']
# Large differences (> 5 BMI units) flag potential issues
```

### 3. Height Source Tracking
Transparency about which height measurements are carried forward vs directly measured:

```python
# height_source column allows sensitivity analysis:
primary_analysis = df  # All participants
sensitivity_analysis = df[df['height_source'] == 'measured']  # Measured only
```

### 4. Temporal Matching Quality Flags
Built-in quality assessment helps decide whether to exclude low-quality matches:

```python
# Filter by quality
high_quality = merged_df[merged_df['bmi_quality'].isin(['excellent', 'good'])]

# Or stratify analysis
by_quality = merged_df.groupby('bmi_quality').apply(analysis_function)
```

---

## OMOP Data Sources

### Measurement Table Structure
```sql
SELECT
    person_id,
    measurement_date,
    measurement_concept_id,  -- Which type of measurement?
    value_as_number,        -- Numeric value
    unit_concept_id         -- Unit of measurement
FROM measurement
WHERE measurement_concept_id IN (
    -- Weight concept IDs
    -- Height concept IDs
    -- BMI concept IDs
)
```

### Concept Table (Units)
```sql
SELECT
    concept_id,
    concept_name  -- e.g., 'kilogram', 'pound (US)', 'centimeter'
FROM concept
WHERE concept_id = measurement.unit_concept_id
```

---

## Convenience Wrappers

### For GWAS
```python
def bmi_for_gwas(cohort_dfs, version, time_col='enrollment_date'):
    """Strict prior-only matching (365 days)"""
    return harmonize_bmi_complete(
        cohort_dfs, version,
        matching_strategy='strict_prior',
        time_col=time_col
    )
```

### For Drug Studies
```python
def bmi_for_drug_study(cohort_dfs, version, time_col='drug_start_date'):
    """Recent baseline only (90 days prior)"""
    return harmonize_bmi_complete(
        cohort_dfs, version,
        matching_strategy='drug_study',
        time_col=time_col
    )
```

### For Longitudinal Studies
```python
def bmi_flexible(cohort_dfs, version, time_col='time_zero'):
    """Longer windows (730 days prior, 180 days post)"""
    return harmonize_bmi_complete(
        cohort_dfs, version,
        matching_strategy='flexible',
        time_col=time_col
    )
```

---

## Tips for Claude Code Reuse

### Adapting for New Phenotypes
If you need similar harmonization for other measurements (e.g., blood pressure, cholesterol):

1. **Copy the extraction pattern** from `extract_bmi_data()`:
   - Replace concept IDs with target measurement concepts
   - Keep the unit handling logic

2. **Reuse cleaning functions** as-is:
   - `clean_units()` - just update the problematic unit list
   - `apply_clean_outliers()` - no changes needed
   - `convert_units()` - update conversion factors

3. **Keep temporal matching unchanged**:
   - `hierarchical_temporal_matching()` works for any measurement

### Borrowing Components
```python
# Import just what you need
from bmi_functions import (
    polars_gbq,                        # BigQuery helper
    hierarchical_temporal_matching,     # Temporal matching
    analyze_temporal_matching_quality,  # QC metrics
    plot_temporal_matching_quality      # QC visualization
)

# Use with your own data extraction
my_data = extract_my_phenotype(person_ids, version)
my_data_cleaned = my_cleaning_pipeline(my_data)
my_data_pandas = my_data_cleaned.to_pandas()

# Temporal matching
matched_cohort = hierarchical_temporal_matching(
    cohort_df, my_data_pandas,
    time_col='diagnosis_date',
    include_post=False
)
```

### Extending for Longitudinal Data
If you need **multiple** BMI measurements per person (not just closest):

```python
def get_all_bmi_in_window(cohort_df, bmi_df, time_col, window_before, window_after):
    """Get ALL BMI measurements within time window"""
    merged = pd.merge(cohort_df, bmi_df, on='person_id')
    merged['days_diff'] = (merged['date'] - merged[time_col]).dt.days
    in_window = (
        (merged['days_diff'] >= -window_before) &
        (merged['days_diff'] <= window_after)
    )
    return merged[in_window].sort_values(['person_id', 'date'])
```

---

## Limitations and Assumptions

1. **Adult populations only**: 3-year height carry-forward inappropriate for children
2. **Missing data mechanism**: Assumes BMI missing at random (MAR)
3. **Temporal bias**: Prefer prior measurements assumes causality direction
4. **Unit trust**: Assumes recorded units are accurate (not always true in EHR)
5. **Extreme outliers**: 4-sigma may still include some implausible values
6. **Pregnancy**: Does not handle pregnancy-related weight changes specially

---

## Dependencies

```python
# Core data manipulation
pandas >= 1.3
polars >= 0.15
numpy >= 1.20

# Visualization
matplotlib >= 3.3
seaborn >= 0.11

# Cloud integration
google-cloud-bigquery >= 2.0

# Utilities
pyarrow >= 6.0  # For BigQuery → polars
```

---

## Metadata
- **Author**: Bennett Waxse
- **Created**: June 2025
- **Platform**: All of Us Research Workbench
- **Data Version**: CDR v8
- **License**: MIT
- **Status**: Production-ready, validated on All of Us data

---

## References

### Validated Concept IDs
All concept IDs have been validated against All of Us Controlled Tier Dataset v8 to ensure they capture relevant measurements without introducing noise.

### Temporal Matching Best Practices
The hierarchical matching approach is based on epidemiological best practices for:
- Case-control studies (prefer prior exposure assessment)
- Cohort studies (baseline covariate measurement)
- Clinical trials (pre-treatment characteristics)

### Height Carry-Forward Rationale
3-year window based on:
- Adult height stability (minimal change after growth plate closure)
- Balance between data availability and measurement validity
- Clinical practice guidelines for using historical height measurements
