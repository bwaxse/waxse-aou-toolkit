# Applying the rcox_phewas_37_z_stat PheRS to a New Cohort

## Quick Start

This notebook (`C.01 Apply rcox_phewas_37_z_stat to new cohort.py`) allows you to apply the validated residualized PheRS (Phenotype Risk Score) to your own cohort in the All of Us Research Program.

**Estimated runtime**: 5-15 minutes depending on cohort size

**Requirements**:
- Access to All of Us Research Workbench with Verily JupyterLab
- A cohort file with `person_id` and `time_zero` columns
- Access to shared reference data files in your workspace bucket

---

## What is This PheRS?

**PheRS** = Phenotype Risk Score

A composite score measuring the burden of NEW (incident) medical diagnoses appearing after a specific index event.

### This Specific Score: `rcox_phewas_37_z_stat_phers`

**Components**:
- **37 phecodes** (phenotype codes): Selected via Cox proportional hazards PheWAS analysis
- **z-statistic weights**: Each phecode weighted by effect size AND statistical precision (HR / SE)
- **Residualized** ("r" prefix): Adjusted for age, sex, race/ethnicity, and EHR characteristics
- **NEW diagnoses only**: Incident ICD codes 28-365 days post-index-event

### Interpretation

| Score | Interpretation |
|-------|-----------------|
| 0 | No qualifying diagnoses detected in post-index period |
| Negative | Fewer diagnoses than expected for this demographic group |
| Positive | More diagnoses than expected for this demographic group |
| Magnitude | Larger absolute values = stronger departure from demographic predictions |

**Example**:
- A person in their 50s with good EHR coverage might be expected to have an rPheRS near 0
- A positive rPheRS of +1.5 means they have more post-index diagnoses than typical for their age/sex/race
- A negative rPheRS of -0.8 means they have fewer diagnoses than expected

---

## How to Use This Notebook

### Step 1: Prepare Your Cohort File

Create a CSV or TSV file with your cohort data. **Required columns**:

- `person_id` (integer): All of Us participant ID
- `time_zero` (datetime): Date of the index event (e.g., COVID diagnosis, enrollment date)

**Example**:
```
person_id,time_zero
1000001,2020-03-15
1000002,2020-03-18
1000003,2020-04-02
```

**Notes**:
- Use ISO format for dates: `YYYY-MM-DD` or `YYYY-MM-DDTHH:MM:SS`
- Keep the file simple—other columns are okay but won't be used
- Save as `.csv` or `.tsv`

### Step 2: Upload to Verily Workbench

1. Log into your All of Us Research Workbench
2. Open JupyterLab
3. Upload this notebook file to your workspace
4. Upload your cohort file to your workspace

### Step 3: Run Setup Cells

Before running this notebook, initialize your workspace:

```python
# Run this in a separate cell FIRST
%env WORKSPACE_CDR=<your-workspace-cdr>
%env WORKSPACE_BUCKET=<your-workspace-bucket>
```

Or run `_reference/verily/00_setup_workspace.ipynb` once per session.

### Step 4: Execute the Notebook

Run through all cells in order. The notebook will:

1. **Load reference data** from your bucket:
   - 37 phecode weights with z-statistics
   - Phecode-to-ICD mappings

2. **Ask for your cohort file path** (interactive prompt)

3. **Extract all ICD codes** for your cohort (by querying condition_occurrence)

4. **Identify NEW post-index diagnoses** (28-365 days after time_zero, not present pre-index)

5. **Calculate raw PheRS** (sum of z-statistic weights)

6. **Retrieve demographics** for your cohort

7. **Residualize PheRS** (fit regression model on your cohort to remove demographic confounding)

8. **Export results** to CSV/TSV

### Step 5: Download Results

The notebook will ask for an output file path and save your cohort with PheRS scores:

```
person_id,time_zero,cox_phewas_37_z_stat_phers,rcox_phewas_37_z_stat_phers
1000001,2020-03-15,2.45,0.82
1000002,2020-03-18,0.00,-1.23
1000003,2020-04-02,3.15,1.45
```

**Output columns**:
- `cox_phewas_37_z_stat_phers`: Raw PheRS (sum of weights)
- `rcox_phewas_37_z_stat_phers`: Residualized PheRS (adjusted for demographics)

---

## Data Flow and Methodology

### 1. ICD Event Extraction

The notebook queries both:
- **condition_occurrence table**: Primary diagnoses
- **observation table**: Secondary diagnoses or codes stored as observations

For each ICD code, it captures:
- ICD9CM codes
- ICD10CM codes
- Special handling for ICD9CM V codes (which overlap with ICD10CM)

**SQL approach**: Uses the `phetk_icd_query()` function which searches four different join paths to maximize code capture.

### 2. Identifying NEW Post-Index Events

An ICD code is "NEW" if:
```
Present in post-index period (28-365 days after time_zero) AND
NOT present in pre-index period (any time before index)
```

**Timeline**:
```
Time ←─────────────── Pre-index ──────────────→ | time_zero | ←─ 28 days → | 28-365 days window | ←─ Post window
                     (All history)              |           |              (Target period)
                                                0 days      28              365 days
```

**Why this approach?**
- Captures incident conditions (new diagnoses), not existing comorbidities
- 28-day lag excludes acute symptoms coinciding with index event
- 365-day window focuses on chronic sequelae

### 3. PheRS Calculation

Raw PheRS = Σ (z-statistic weight for each phecode present)

```python
# Example with 3 phecodes
Person has:
  - Phecode A (weight = 1.5)
  - Phecode C (weight = 0.8)

Raw PheRS = 1.5 + 0.8 = 2.3
```

**Weighting**: Uses z-statistics (HR / SE) from Cox models, which:
- Emphasizes effect sizes (larger HR = larger weight)
- Accounts for precision (lower SE = more confident estimate)
- Balances both statistical and clinical significance

### 4. Residualization

Fits a linear regression model:
```
raw_phers ~ age_spline + sex + race_ethnicity + has_icd9 + has_icd10 + ehr_length
```

Then calculates studentized residuals:
```
residual = observed - predicted
residualized = residual / sqrt(MSE)
```

**Why residualize?**
- Removes confounding from demographics
- Allows fair comparison across different age/sex/racial groups
- Makes rPheRS more interpretable as "departure from demographic expectations"
- Better for association analyses (avoids collider bias)

**Important**: Residualization is fit on YOUR cohort, not transferred from original analysis. This accounts for demographic and EHR differences specific to your population.

---

## Key Statistics and Validation

The notebook automatically reports:

### ICD Extraction
- Number of ICD events extracted
- Distribution across ICD9CM vs ICD10CM
- Unique participants with data

### Event Detection
- Pre-index events (baseline history)
- Post-index events (28-365 days)
- NEW events (in post but not pre)
- Participants with ≥1 new ICD code

### Phecode Coverage
- Phecodes with ≥1 person
- Median number of phecodes per person

### Raw PheRS Distribution
- Mean, median, SD
- % with PheRS = 0 (no qualifying diagnoses)
- Visualization of distribution

### Residualization Quality
- R² of regression model
- Mean/SD of raw vs residualized
- Check for extreme outliers
- Visualization of before/after

---

## Common Issues and Troubleshooting

### Problem: "Could not find phecode weights file"

**Solution**: Ensure file exists in your bucket at:
```
gs://<your-bucket>/data/long_covid/phers/variable_hr/long_covid_min_1_pct_37_phewas_train_phecode_weights.tsv
```

Contact your analysis coordinator to confirm the correct path.

### Problem: Very small numbers of participants have PheRS > 0

**Possible causes**:
1. Your cohort has very recent time_zero dates (not enough time for diagnoses to accumulate)
2. Your cohort has good health status (fewer diagnoses expected)
3. Your cohort has short EHR history (codes may exist but predate index)

**Check**:
- Confirm time_zero dates are correct
- Look at the raw statistics from ICD extraction
- Check "Participants with ≥1 new ICD" count

### Problem: All residuals are zero

**Cause**: Likely all participants have missing demographic data (especially sex)

**Solution**:
- Check your demographics query
- Verify sex_at_birth field is populated in the person table
- The notebook will report how many were excluded due to missing sex

### Problem: Query times out

**Solutions**:
1. Start with a smaller cohort to test
2. Run during off-peak hours
3. Contact your workspace administrator if persists

### Problem: "Cohort file must contain columns: ['person_id', 'time_zero']"

**Solution**: Verify your CSV/TSV has these exact column names. Column names are case-sensitive.

Correct:
```csv
person_id,time_zero
1000001,2020-03-15
```

Incorrect:
```csv
Person_ID,Time_Zero     ← Case mismatch
PersonID,time_zero       ← Wrong names
```

---

## Using Results in Analysis

### For Descriptive Statistics

Report both raw and residualized:
```
Mean (SD) raw PheRS: 2.1 (1.8)
Mean (SD) rPheRS: 0.1 (0.9)
N with rPheRS = 0: 342 (15%)
```

### For Association Analysis (e.g., GWAS, phenotype association)

Use `rcox_phewas_37_z_stat_phers` (residualized) as:
- **Continuous outcome**: Standard linear/logistic regression
- **Binary outcome**: Dichotomize at 0 or median
- **Stratification**: Stratify by quartiles

**Why not raw PheRS?**
- Raw scores are confounded by age, sex, race/ethnicity
- Using raw in association studies leads to spurious findings
- Residualized is more interpretable and valid

### For Visualization

Example visualization code:
```python
import matplotlib.pyplot as plt
import seaborn as sns

# Compare by demographic group
fig, ax = plt.subplots()
sns.violinplot(data=results, x='race', y='rcox_phewas_37_z_stat_phers', ax=ax)
ax.axhline(0, color='red', linestyle='--', alpha=0.5)
plt.title('Residualized PheRS by Race/Ethnicity')
plt.show()
```

---

## Technical Details

### Dependencies

Standard scientific Python stack (all pre-installed in Verily):
- `pandas`, `polars`, `numpy`, `scipy`
- `scikit-learn` (for SplineTransformer)
- `google-cloud-bigquery` (for OMOP queries)
- `matplotlib`, `seaborn` (for visualization)

**No external packages required** — this notebook is completely self-contained.

### Database Queries

The notebook uses the following OMOP tables:
- `condition_occurrence`: Primary diagnoses
- `observation`: Secondary diagnoses/codes
- `concept`: ICD code definitions
- `concept_relationship`: For V-code vocabulary disambiguation
- `person`: Demographics (DOB, sex, race, ethnicity)

### Computational Requirements

- **Machine type**: n2-standard-16 (16 vCPU, 64GB RAM) works for most cohorts
- **Disk space**: ~1-2 GB for intermediate files
- **Network**: ~5-20 minutes depending on cohort size and network speed

### Reproducibility

This notebook:
1. Uses deterministic SQL queries (no randomization)
2. Specifies sklearn seed for reproducibility (if needed for sklearn operations)
3. Saves all parameters (phecode list, weights, spline configuration)
4. Documents exact methodology for residualization

To reproduce results: Run notebook on same CDR version with same cohort file.

---

## Citation and Acknowledgments

**Original Methodology**:
These PheRS were developed in analysis of post-COVID outcomes using Cox proportional hazards PheWAS.

**When publishing results using this score, please cite**:
- The original analysis publication (provided separately)
- This implementation notebook (include link/version)
- Specify CDR version (e.g., "CDR v8")

**Example citation**:
> Phenotype risk scores were calculated using the residualized Cox PheWAS 37-phecode score with z-statistic weighting, residualized for age, sex, race/ethnicity, and EHR characteristics within this cohort.

---

## Questions or Issues?

1. **Technical issues**: Contact your Verily workspace administrator
2. **Methodology questions**: Contact the analysis team (email/contact info provided separately)
3. **Data access issues**: Check All of Us Help Desk: https://support.researchallofus.org/
4. **Notebook bugs**: Report with complete error message and cohort size

---

## Additional Resources

### Useful All of Us Documentation
- [OMOP CDM Overview](https://www.ohdsi.org/web/wiki/doku.php?id=documentation:cdm:single-page)
- [Condition Occurrence Table](https://www.ohdsi.org/web/wiki/doku.php?id=documentation:cdm:condition_occurrence)
- [ICD Code Extraction Best Practices](https://www.ohdsi.org/web/wiki/doku.php?id=documentation:sql_101)

### PheWAS and PheRS Literature
- PheWAS methodology: [PheWAS Wikipedia](https://en.wikipedia.org/wiki/Phenome-wide_association_study)
- Phecode system: [PheCode website](https://phewascatalog.org/)

### Related Analyses
- HPV GWAS pipeline: See `hpv/CLAUDE.md`
- Sarcoid GWAS pipeline: See `sarcoid/CLAUDE.md`
- Wearables analysis: See `lc_wearables/CLAUDE.md`

---

**Last Updated**: January 2026
**Notebook Version**: 1.0
**Status**: Verified and ready for use
