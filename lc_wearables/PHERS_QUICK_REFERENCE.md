# PheRS Application - Quick Reference Guide

## Pre-Flight Checklist

- [ ] I have access to All of Us Research Workbench with Verily JupyterLab
- [ ] I have a CSV/TSV file with `person_id` and `time_zero` columns
- [ ] I have uploaded this file to my Verily workspace
- [ ] My workspace bucket path is: `gs://fc-secure-XXXXXXXXX/`
- [ ] I've run `_reference/verily/00_setup_workspace.ipynb` (or will run setup in notebook)
- [ ] I have the phecode weight file in my bucket (ask coordinator for location)
- [ ] I have the phecode-ICD mapping file in my bucket

## Cohort File Template

```csv
person_id,time_zero
1000001,2020-03-15
1000002,2020-03-18
1000003,2020-04-02
```

**Column requirements**:
- `person_id`: integer (All of Us ID)
- `time_zero`: ISO date format (YYYY-MM-DD)

## Step-by-Step Execution

### 1. Upload Files
- [ ] Notebook: `C.01 Apply rcox_phewas_37_z_stat to new cohort.py`
- [ ] Your cohort file (*.csv or *.tsv)

### 2. Initialize Environment
```python
# Run first (if not already setup)
import os
os.environ['WORKSPACE_CDR'] = 'project.dataset'  # From your setup
os.environ['WORKSPACE_BUCKET'] = 'gs://fc-secure-xyz'  # From your setup
```

### 3. Run Notebook
- Execute each cell in order
- When prompted, enter your cohort file path
- When prompted, enter desired output file path

### 4. Download Results
- Find output file in Verily Files
- Download to your computer
- Use for analysis!

## Output Columns Explained

| Column | Type | Range | Interpretation |
|--------|------|-------|-----------------|
| `person_id` | integer | varies | All of Us participant ID |
| `time_zero` | date | varies | Your index event date |
| `cox_phewas_37_z_stat_phers` | float | 0 to ~20 | **Raw PheRS** - sum of weights |
| `rcox_phewas_37_z_stat_phers` | float | ~-3 to +3 | **Residualized PheRS** - adjusted |

## Score Interpretation Cheat Sheet

```
rPheRS < -2     : Much less disease burden than expected
rPheRS -1 to 0  : About as expected, maybe slightly less
rPheRS 0        : Exactly as expected
rPheRS 0 to +1  : About as expected, maybe slightly more
rPheRS > +2     : Much more disease burden than expected
```

## What to Report in Methods Section

> Residualized Cox PheWAS 37-phecode PheRS was calculated for [N] participants. Raw PheRS values (mean ± SD: X ± Y) were computed as the sum of z-statistic-weighted phecode burdens for incident diagnoses occurring 28-365 days post-index-event. Residualization was performed using linear regression adjusted for age (B-spline, 3 knots), sex, race/ethnicity categories, ICD9/10 vocabulary presence, and EHR length. Resulting studentized residuals (mean ± SD: A ± B) were used for primary analyses.

## What to Report in Results Section

### Baseline Characteristics (Example)
```
N = XXXX
Age, mean (SD): XX (XX) years
Female, N (%): XXX (XX%)
Raw PheRS, mean (SD): X.X (X.X)
Residualized PheRS, mean (SD): 0.1 (0.9)
  Median (IQR): 0.0 (−0.6, 0.8)
  Range: −2.8 to 3.2
Participants with PheRS = 0, N (%): XXX (XX%)
```

## Visualization Examples

### Plot 1: Distribution
```python
import matplotlib.pyplot as plt
plt.hist(results['rcox_phewas_37_z_stat_phers'], bins=50)
plt.axvline(0, color='red', linestyle='--')
plt.xlabel('Residualized PheRS')
plt.ylabel('Number of Participants')
plt.title('Distribution of Residualized PheRS (N=XXXX)')
```

### Plot 2: By Sex
```python
import seaborn as sns
sns.violinplot(data=results, x='sex', y='rcox_phewas_37_z_stat_phers')
plt.axhline(0, color='red', linestyle='--', alpha=0.5)
plt.ylabel('Residualized PheRS')
```

### Plot 3: By Age Group
```python
results['age_group'] = pd.cut(results['age_years'], bins=[0, 30, 50, 70, 100])
sns.boxplot(data=results, x='age_group', y='rcox_phewas_37_z_stat_phers')
plt.axhline(0, color='red', linestyle='--', alpha=0.5)
```

## Key Residualization Statistics (From Notebook)

Look for these in the output:
- **R²**: Should be ~0.3-0.7 (higher = more variance explained by demographics)
- **Mean residual**: Should be ≈ 0 (check: usually ±0.001 or smaller)
- **SD residual**: Should be ≈ 1 (standardized scale)
- **N samples**: Should match your cohort size

If R² is very low (< 0.1):
- May indicate your cohort has unusual demographics
- May indicate poor EHR coverage
- Consider noting in limitations

## Troubleshooting Quick Hits

| Error | Likely Cause | Fix |
|-------|-------------|-----|
| FileNotFoundError for weights | Weights not in bucket | Check path with coordinator |
| "Required columns not found" | CSV has wrong column names | Use exact names: `person_id`, `time_zero` |
| Query timeout | Too large cohort or slow network | Run at off-peak time or smaller subset |
| All PheRS = 0 | No incident diagnoses detected | Check time_zero dates and cohort definition |
| Missing residuals | Missing sex in demographics | Rare; notebook will report how many excluded |

## File Naming Suggestions

When saving your output:
```
cohort_NAME_phers_YYYYMMDD.csv
cohort_STUDYNAME_N1000_phers_pheno_20260115.tsv
```

Include:
- Cohort identifier
- N or size indicator
- Date (for versioning)

## Important Reminders

1. **Use residualized scores** for analyses (not raw)
2. **Report both** raw and residualized in tables
3. **Zero = no qualifying diagnoses**, not missing data
4. **Don't compare across CDR versions** (results may differ)
5. **Save your code** - make notebook reproducible
6. **Check diagnostics** - review R², residual plots, outliers
7. **Document assumptions** - what is your cohort? how defined?

## Reference Data Locations

Ask your coordinator to confirm these paths in your bucket:

**Phecode weights** (required):
```
gs://fc-secure-XXXXXXX/data/long_covid/phers/variable_hr/
  long_covid_min_1_pct_37_phewas_train_phecode_weights.tsv
```

**Phecode-ICD mapping** (required):
```
gs://fc-secure-XXXXXXX/data/pasc/phers/
  phecodeX_unrolled_ICD_CM.csv
```

## Contact/Support

- **Technical issues**: Verily workspace support
- **Data questions**: Contact analysis coordinator
- **Methodology**: See README_PHERS_APPLICATION.md
- **Citation**: See full README for how to cite

## Notes / To Do

```
My cohort details:
- N = _____
- Index event: _________________
- Date range: ___ to ___

Files saved:
- Input cohort: _____________________
- Output results: _____________________

Results to report:
- [ ] N and sample characteristics
- [ ] Raw PheRS: mean (SD), range
- [ ] Residualized PheRS: mean (SD), range
- [ ] % with PheRS = 0
- [ ] Key descriptive plots
- [ ] R² from residualization model
```

---

**Print this page and keep it handy while running the analysis!**
