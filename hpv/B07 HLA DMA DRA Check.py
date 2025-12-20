#!/usr/bin/env python
# coding: utf-8

# In[ ]:


### Sept 15, 2025
### Bennett Waxse


# In[ ]:


#PheTK installation
#!pip install PheTK --upgrade


# In[ ]:


### load all apps
import pandas as pd
import polars as pl
import numpy as np
import os
import subprocess
import statsmodels.api as sm
from scipy.stats import fisher_exact

version = get_ipython().run_line_magic('env', 'WORKSPACE_CDR')
my_bucket = os.getenv('WORKSPACE_BUCKET')


# # Anal

# In[ ]:


## HLA results
hla_cohort_eur = pl.read_csv(f'{my_bucket}/HLA/HVP_EUR_Phenotype_HLA.csv')
hla_cohort_afr = pl.read_csv(f'{my_bucket}/HLA/HVP_AFR_Phenotype_HLA.csv')
hla_cohort_amr = pl.read_csv(f'{my_bucket}/HLA/HVP_AMR_Phenotype_HLA.csv')


# ## DMA

# In[ ]:


hla_cohort_eur = hla_cohort_eur.with_columns(
    pl.when(pl.col('DMA.01.01.01G').is_null())
    .then(1.0)
    .otherwise(2.0)
    .alias('DMA.01.01.01G')
)


# In[ ]:


hla_cohort_afr = hla_cohort_afr.with_columns(
    pl.when(pl.col('DMA.01.01.01G').is_null())
    .then(1.0)
    .otherwise(2.0)
    .alias('DMA.01.01.01G')
)


# In[ ]:


hla_cohort_amr = hla_cohort_amr.with_columns(
    pl.when(pl.col('DMA.01.01.01G').is_null())
    .then(1.0)
    .otherwise(2.0)
    .alias('DMA.01.01.01G')
)


# In[ ]:


hla_cohort_eur.group_by(pl.col('DMA.01.01.01G', 'condition__hpv')).len()


# In[ ]:


hla_cohort_afr.group_by(pl.col('DMA.01.01.01G', 'condition__hpv')).len()


# In[ ]:


hla_cohort_amr.group_by(pl.col('DMA.01.01.01G', 'condition__hpv')).len()


# ### Fisher's Exact

# In[ ]:


All_No_Missing=hla_cohort_eur.to_pandas()

# Set outcome and covariates
outcome = "condition__hpv"
covariates = [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]


# Identify predictor columns (starting from column index 30)
predictor_cols = list(All_No_Missing.columns[29:])
predictor_cols = ['DMA.01.01.01G']

# Initialize list for results
results = []

# Subset and drop rows with NA
outcome = 'condition__hpv'
predictor = 'DMA.01.01.01G'

cols_to_use = [outcome, predictor] + [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]

data_subset = All_No_Missing[cols_to_use].dropna()

# Define design matrix X and target y
X = data_subset[[predictor] + covariates]
X = sm.add_constant(X)
y = data_subset[outcome]

# Count cases and controls
cases = (y == 1).sum()
controls = (y == 0).sum()

# Check if predictor is binary/categorical
unique_vals = X[predictor].nunique()

if unique_vals == 2:
    # Binary predictor - use Fisher's exact test
    contingency_table = pd.crosstab(X[predictor], y)
    print(f"Contingency table for {predictor}:")
    print(contingency_table)

    # Perform Fisher's exact test
    oddsratio, pval = fisher_exact(contingency_table)

    # Convert to log-odds (coefficient) and approximate SE
    coef = np.log(oddsratio) if oddsratio > 0 and np.isfinite(oddsratio) else np.nan

    # Approximate SE using Woolf's method for 2x2 tables
    if not contingency_table.values.min() == 0:  # No zero cells
        a, b, c, d = contingency_table.values.flatten()
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    else:
        # Add continuity correction for zero cells
        a, b, c, d = contingency_table.values.flatten() + 0.5
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

    print(f"Fisher's exact: OR={oddsratio:.3f}, p={pval:.2e}")


# In[ ]:


All_No_Missing=hla_cohort_afr.to_pandas()

# Set outcome and covariates
outcome = "condition__hpv"
covariates = [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]


# Identify predictor columns (starting from column index 30)
predictor_cols = list(All_No_Missing.columns[29:])
predictor_cols = ['DMA.01.01.01G']

# Initialize list for results
results = []

# Subset and drop rows with NA
outcome = 'condition__hpv'
predictor = 'DMA.01.01.01G'

cols_to_use = [outcome, predictor] + [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]

data_subset = All_No_Missing[cols_to_use].dropna()

# Define design matrix X and target y
X = data_subset[[predictor] + covariates]
X = sm.add_constant(X)
y = data_subset[outcome]

# Count cases and controls
cases = (y == 1).sum()
controls = (y == 0).sum()

# Check if predictor is binary/categorical
unique_vals = X[predictor].nunique()

if unique_vals == 2:
    # Binary predictor - use Fisher's exact test
    contingency_table = pd.crosstab(X[predictor], y)
    print(f"Contingency table for {predictor}:")
    print(contingency_table)

    # Perform Fisher's exact test
    oddsratio, pval = fisher_exact(contingency_table)

    # Convert to log-odds (coefficient) and approximate SE
    coef = np.log(oddsratio) if oddsratio > 0 and np.isfinite(oddsratio) else np.nan

    # Approximate SE using Woolf's method for 2x2 tables
    if not contingency_table.values.min() == 0:  # No zero cells
        a, b, c, d = contingency_table.values.flatten()
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    else:
        # Add continuity correction for zero cells
        a, b, c, d = contingency_table.values.flatten() + 0.5
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

    print(f"Fisher's exact: OR={oddsratio:.3f}, p={pval:.2e}")


# In[ ]:


All_No_Missing=hla_cohort_amr.to_pandas()

# Set outcome and covariates
outcome = "condition__hpv"
covariates = [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]


# Identify predictor columns (starting from column index 30)
predictor_cols = list(All_No_Missing.columns[29:])
predictor_cols = ['DMA.01.01.01G']

# Initialize list for results
results = []

# Subset and drop rows with NA
outcome = 'condition__hpv'
predictor = 'DMA.01.01.01G'

cols_to_use = [outcome, predictor] + [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]

data_subset = All_No_Missing[cols_to_use].dropna()

# Define design matrix X and target y
X = data_subset[[predictor] + covariates]
X = sm.add_constant(X)
y = data_subset[outcome]

# Count cases and controls
cases = (y == 1).sum()
controls = (y == 0).sum()

# Check if predictor is binary/categorical
unique_vals = X[predictor].nunique()

if unique_vals == 2:
    # Binary predictor - use Fisher's exact test
    contingency_table = pd.crosstab(X[predictor], y)
    print(f"Contingency table for {predictor}:")
    print(contingency_table)

    # Perform Fisher's exact test
    oddsratio, pval = fisher_exact(contingency_table)

    # Convert to log-odds (coefficient) and approximate SE
    coef = np.log(oddsratio) if oddsratio > 0 and np.isfinite(oddsratio) else np.nan

    # Approximate SE using Woolf's method for 2x2 tables
    if not contingency_table.values.min() == 0:  # No zero cells
        a, b, c, d = contingency_table.values.flatten()
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    else:
        # Add continuity correction for zero cells
        a, b, c, d = contingency_table.values.flatten() + 0.5
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

    print(f"Fisher's exact: OR={oddsratio:.3f}, p={pval:.2e}")


# ## DRA

# In[ ]:


hla_cohort_eur = hla_cohort_eur.with_columns(
    pl.when(pl.col('DRA.01.05').is_null())
    .then(1.0)
    .otherwise(2.0)
    .alias('DRA.01.05')
)


# In[ ]:


hla_cohort_afr = hla_cohort_afr.with_columns(
    pl.when(pl.col('DRA.01.05').is_null())
    .then(1.0)
    .otherwise(2.0)
    .alias('DRA.01.05')
)


# In[ ]:


hla_cohort_amr = hla_cohort_amr.with_columns(
    pl.when(pl.col('DRA.01.05').is_null())
    .then(1.0)
    .otherwise(2.0)
    .alias('DRA.01.05')
)


# In[ ]:


hla_cohort_eur.group_by(pl.col('DRA.01.05', 'condition__hpv')).len()


# In[ ]:


hla_cohort_afr.group_by(pl.col('DRA.01.05', 'condition__hpv')).len()


# In[ ]:


hla_cohort_amr.group_by(pl.col('DRA.01.05', 'condition__hpv')).len()


# ### Fisher's Exact

# In[ ]:


All_No_Missing=hla_cohort_eur.to_pandas()

# Set outcome and covariates
outcome = "condition__hpv"
covariates = [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]


# Identify predictor columns (starting from column index 30)
predictor_cols = list(All_No_Missing.columns[29:])
predictor_cols = ['DRA.01.05']

# Initialize list for results
results = []

# Subset and drop rows with NA
outcome = 'condition__hpv'
predictor = 'DRA.01.05'

cols_to_use = [outcome, predictor] + [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]

data_subset = All_No_Missing[cols_to_use].dropna()

# Define design matrix X and target y
X = data_subset[[predictor] + covariates]
X = sm.add_constant(X)
y = data_subset[outcome]

# Count cases and controls
cases = (y == 1).sum()
controls = (y == 0).sum()

# Check if predictor is binary/categorical
unique_vals = X[predictor].nunique()

if unique_vals == 2:
    # Binary predictor - use Fisher's exact test
    contingency_table = pd.crosstab(X[predictor], y)
    print(f"Contingency table for {predictor}:")
    print(contingency_table)

    # Perform Fisher's exact test
    oddsratio, pval = fisher_exact(contingency_table)

    # Convert to log-odds (coefficient) and approximate SE
    coef = np.log(oddsratio) if oddsratio > 0 and np.isfinite(oddsratio) else np.nan

    # Approximate SE using Woolf's method for 2x2 tables
    if not contingency_table.values.min() == 0:  # No zero cells
        a, b, c, d = contingency_table.values.flatten()
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    else:
        # Add continuity correction for zero cells
        a, b, c, d = contingency_table.values.flatten() + 0.5
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

    print(f"Fisher's exact: OR={oddsratio:.3f}, p={pval:.2e}")


# In[ ]:


All_No_Missing=hla_cohort_afr.to_pandas()

# Set outcome and covariates
outcome = "condition__hpv"
covariates = [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]


# Identify predictor columns (starting from column index 30)
predictor_cols = list(All_No_Missing.columns[29:])
predictor_cols = ['DRA.01.05']

# Initialize list for results
results = []

# Subset and drop rows with NA
outcome = 'condition__hpv'
predictor = 'DRA.01.05'

cols_to_use = [outcome, predictor] + [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]

data_subset = All_No_Missing[cols_to_use].dropna()

# Define design matrix X and target y
X = data_subset[[predictor] + covariates]
X = sm.add_constant(X)
y = data_subset[outcome]

# Count cases and controls
cases = (y == 1).sum()
controls = (y == 0).sum()

# Check if predictor is binary/categorical
unique_vals = X[predictor].nunique()

if unique_vals == 2:
    # Binary predictor - use Fisher's exact test
    contingency_table = pd.crosstab(X[predictor], y)
    print(f"Contingency table for {predictor}:")
    print(contingency_table)

    # Perform Fisher's exact test
    oddsratio, pval = fisher_exact(contingency_table)

    # Convert to log-odds (coefficient) and approximate SE
    coef = np.log(oddsratio) if oddsratio > 0 and np.isfinite(oddsratio) else np.nan

    # Approximate SE using Woolf's method for 2x2 tables
    if not contingency_table.values.min() == 0:  # No zero cells
        a, b, c, d = contingency_table.values.flatten()
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    else:
        # Add continuity correction for zero cells
        a, b, c, d = contingency_table.values.flatten() + 0.5
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

    print(f"Fisher's exact: OR={oddsratio:.3f}, p={pval:.2e}")


# In[ ]:


All_No_Missing=hla_cohort_amr.to_pandas()

# Set outcome and covariates
outcome = "condition__hpv"
covariates = [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]


# Identify predictor columns (starting from column index 30)
predictor_cols = list(All_No_Missing.columns[29:])
predictor_cols = ['DMA.01.01.01G']

# Initialize list for results
results = []

# Subset and drop rows with NA
outcome = 'condition__hpv'
predictor = 'DMA.01.01.01G'

cols_to_use = [outcome, predictor] + [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]

data_subset = All_No_Missing[cols_to_use].dropna()

# Define design matrix X and target y
X = data_subset[[predictor] + covariates]
X = sm.add_constant(X)
y = data_subset[outcome]

# Count cases and controls
cases = (y == 1).sum()
controls = (y == 0).sum()

# Check if predictor is binary/categorical
unique_vals = X[predictor].nunique()

if unique_vals == 2:
    # Binary predictor - use Fisher's exact test
    contingency_table = pd.crosstab(X[predictor], y)
    print(f"Contingency table for {predictor}:")
    print(contingency_table)

    # Perform Fisher's exact test
    oddsratio, pval = fisher_exact(contingency_table)

    # Convert to log-odds (coefficient) and approximate SE
    coef = np.log(oddsratio) if oddsratio > 0 and np.isfinite(oddsratio) else np.nan

    # Approximate SE using Woolf's method for 2x2 tables
    if not contingency_table.values.min() == 0:  # No zero cells
        a, b, c, d = contingency_table.values.flatten()
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    else:
        # Add continuity correction for zero cells
        a, b, c, d = contingency_table.values.flatten() + 0.5
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

    print(f"Fisher's exact: OR={oddsratio:.3f}, p={pval:.2e}")


# # Vulvovaginal Cervical

# In[ ]:


## HLA results
hla_cohort_eur = pl.read_csv(f'{my_bucket}/HLA/HVP_EUR_VV_Phenotype_HLA.csv')
hla_cohort_afr = pl.read_csv(f'{my_bucket}/HLA/HVP_AFR_VV_Phenotype_HLA.csv')
hla_cohort_amr = pl.read_csv(f'{my_bucket}/HLA/HVP_AMR_VV_Phenotype_HLA.csv')


# ## DMA

# In[ ]:


hla_cohort_eur = hla_cohort_eur.with_columns(
    pl.when(pl.col('DMA.01.01.01G').is_null())
    .then(1.0)
    .otherwise(2.0)
    .alias('DMA.01.01.01G')
)


# In[ ]:


hla_cohort_afr = hla_cohort_afr.with_columns(
    pl.when(pl.col('DMA.01.01.01G').is_null())
    .then(1.0)
    .otherwise(2.0)
    .alias('DMA.01.01.01G')
)


# In[ ]:


hla_cohort_amr = hla_cohort_amr.with_columns(
    pl.when(pl.col('DMA.01.01.01G').is_null())
    .then(1.0)
    .otherwise(2.0)
    .alias('DMA.01.01.01G')
)


# In[ ]:


hla_cohort_eur.group_by(pl.col('DMA.01.01.01G', 'condition__hpv')).len()


# In[ ]:


hla_cohort_afr.group_by(pl.col('DMA.01.01.01G', 'condition__hpv')).len()


# In[ ]:


hla_cohort_amr.group_by(pl.col('DMA.01.01.01G', 'condition__hpv')).len()


# ### Fisher's Exact

# In[ ]:


All_No_Missing=hla_cohort_eur.to_pandas()

# Set outcome and covariates
outcome = "condition__hpv"
covariates = [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]


# Identify predictor columns (starting from column index 30)
predictor_cols = list(All_No_Missing.columns[29:])
predictor_cols = ['DMA.01.01.01G']

# Initialize list for results
results = []

# Subset and drop rows with NA
outcome = 'condition__hpv'
predictor = 'DMA.01.01.01G'

cols_to_use = [outcome, predictor] + [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]

data_subset = All_No_Missing[cols_to_use].dropna()

# Define design matrix X and target y
X = data_subset[[predictor] + covariates]
X = sm.add_constant(X)
y = data_subset[outcome]

# Count cases and controls
cases = (y == 1).sum()
controls = (y == 0).sum()

# Check if predictor is binary/categorical
unique_vals = X[predictor].nunique()

if unique_vals == 2:
    # Binary predictor - use Fisher's exact test
    contingency_table = pd.crosstab(X[predictor], y)
    print(f"Contingency table for {predictor}:")
    print(contingency_table)

    # Perform Fisher's exact test
    oddsratio, pval = fisher_exact(contingency_table)

    # Convert to log-odds (coefficient) and approximate SE
    coef = np.log(oddsratio) if oddsratio > 0 and np.isfinite(oddsratio) else np.nan

    # Approximate SE using Woolf's method for 2x2 tables
    if not contingency_table.values.min() == 0:  # No zero cells
        a, b, c, d = contingency_table.values.flatten()
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    else:
        # Add continuity correction for zero cells
        a, b, c, d = contingency_table.values.flatten() + 0.5
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

    print(f"Fisher's exact: OR={oddsratio:.3f}, p={pval:.2e}")


# In[ ]:


All_No_Missing=hla_cohort_afr.to_pandas()

# Set outcome and covariates
outcome = "condition__hpv"
covariates = [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]


# Identify predictor columns (starting from column index 30)
predictor_cols = list(All_No_Missing.columns[29:])
predictor_cols = ['DMA.01.01.01G']

# Initialize list for results
results = []

# Subset and drop rows with NA
outcome = 'condition__hpv'
predictor = 'DMA.01.01.01G'

cols_to_use = [outcome, predictor] + [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]

data_subset = All_No_Missing[cols_to_use].dropna()

# Define design matrix X and target y
X = data_subset[[predictor] + covariates]
X = sm.add_constant(X)
y = data_subset[outcome]

# Count cases and controls
cases = (y == 1).sum()
controls = (y == 0).sum()

# Check if predictor is binary/categorical
unique_vals = X[predictor].nunique()

if unique_vals == 2:
    # Binary predictor - use Fisher's exact test
    contingency_table = pd.crosstab(X[predictor], y)
    print(f"Contingency table for {predictor}:")
    print(contingency_table)

    # Perform Fisher's exact test
    oddsratio, pval = fisher_exact(contingency_table)

    # Convert to log-odds (coefficient) and approximate SE
    coef = np.log(oddsratio) if oddsratio > 0 and np.isfinite(oddsratio) else np.nan

    # Approximate SE using Woolf's method for 2x2 tables
    if not contingency_table.values.min() == 0:  # No zero cells
        a, b, c, d = contingency_table.values.flatten()
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    else:
        # Add continuity correction for zero cells
        a, b, c, d = contingency_table.values.flatten() + 0.5
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

    print(f"Fisher's exact: OR={oddsratio:.3f}, p={pval:.2e}")


# In[ ]:


All_No_Missing=hla_cohort_amr.to_pandas()

# Set outcome and covariates
outcome = "condition__hpv"
covariates = [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]


# Identify predictor columns (starting from column index 30)
predictor_cols = list(All_No_Missing.columns[29:])
predictor_cols = ['DMA.01.01.01G']

# Initialize list for results
results = []

# Subset and drop rows with NA
outcome = 'condition__hpv'
predictor = 'DMA.01.01.01G'

cols_to_use = [outcome, predictor] + [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]

data_subset = All_No_Missing[cols_to_use].dropna()

# Define design matrix X and target y
X = data_subset[[predictor] + covariates]
X = sm.add_constant(X)
y = data_subset[outcome]

# Count cases and controls
cases = (y == 1).sum()
controls = (y == 0).sum()

# Check if predictor is binary/categorical
unique_vals = X[predictor].nunique()

if unique_vals == 2:
    # Binary predictor - use Fisher's exact test
    contingency_table = pd.crosstab(X[predictor], y)
    print(f"Contingency table for {predictor}:")
    print(contingency_table)

    # Perform Fisher's exact test
    oddsratio, pval = fisher_exact(contingency_table)

    # Convert to log-odds (coefficient) and approximate SE
    coef = np.log(oddsratio) if oddsratio > 0 and np.isfinite(oddsratio) else np.nan

    # Approximate SE using Woolf's method for 2x2 tables
    if not contingency_table.values.min() == 0:  # No zero cells
        a, b, c, d = contingency_table.values.flatten()
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    else:
        # Add continuity correction for zero cells
        a, b, c, d = contingency_table.values.flatten() + 0.5
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

    print(f"Fisher's exact: OR={oddsratio:.3f}, p={pval:.2e}")


# ## DRA

# In[ ]:


hla_cohort_eur = hla_cohort_eur.with_columns(
    pl.when(pl.col('DRA.01.05').is_null())
    .then(1.0)
    .otherwise(2.0)
    .alias('DRA.01.05')
)


# In[ ]:


hla_cohort_afr = hla_cohort_afr.with_columns(
    pl.when(pl.col('DRA.01.05').is_null())
    .then(1.0)
    .otherwise(2.0)
    .alias('DRA.01.05')
)


# In[ ]:


hla_cohort_amr = hla_cohort_amr.with_columns(
    pl.when(pl.col('DRA.01.05').is_null())
    .then(1.0)
    .otherwise(2.0)
    .alias('DRA.01.05')
)


# In[ ]:


hla_cohort_eur.group_by(pl.col('DRA.01.05', 'condition__hpv')).len()


# In[ ]:


hla_cohort_afr.group_by(pl.col('DRA.01.05', 'condition__hpv')).len()


# In[ ]:


hla_cohort_amr.group_by(pl.col('DRA.01.05', 'condition__hpv')).len()


# ### Fisher's Exact

# In[ ]:


All_No_Missing=hla_cohort_eur.to_pandas()

# Set outcome and covariates
outcome = "condition__hpv"
covariates = [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]


# Identify predictor columns (starting from column index 30)
predictor_cols = list(All_No_Missing.columns[29:])
predictor_cols = ['DRA.01.05']

# Initialize list for results
results = []

# Subset and drop rows with NA
outcome = 'condition__hpv'
predictor = 'DRA.01.05'

cols_to_use = [outcome, predictor] + [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]

data_subset = All_No_Missing[cols_to_use].dropna()

# Define design matrix X and target y
X = data_subset[[predictor] + covariates]
X = sm.add_constant(X)
y = data_subset[outcome]

# Count cases and controls
cases = (y == 1).sum()
controls = (y == 0).sum()

# Check if predictor is binary/categorical
unique_vals = X[predictor].nunique()

if unique_vals == 2:
    # Binary predictor - use Fisher's exact test
    contingency_table = pd.crosstab(X[predictor], y)
    print(f"Contingency table for {predictor}:")
    print(contingency_table)

    # Perform Fisher's exact test
    oddsratio, pval = fisher_exact(contingency_table)

    # Convert to log-odds (coefficient) and approximate SE
    coef = np.log(oddsratio) if oddsratio > 0 and np.isfinite(oddsratio) else np.nan

    # Approximate SE using Woolf's method for 2x2 tables
    if not contingency_table.values.min() == 0:  # No zero cells
        a, b, c, d = contingency_table.values.flatten()
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    else:
        # Add continuity correction for zero cells
        a, b, c, d = contingency_table.values.flatten() + 0.5
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

    print(f"Fisher's exact: OR={oddsratio:.3f}, p={pval:.2e}")


# In[ ]:


All_No_Missing=hla_cohort_afr.to_pandas()

# Set outcome and covariates
outcome = "condition__hpv"
covariates = [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]


# Identify predictor columns (starting from column index 30)
predictor_cols = list(All_No_Missing.columns[29:])
predictor_cols = ['DRA.01.05']

# Initialize list for results
results = []

# Subset and drop rows with NA
outcome = 'condition__hpv'
predictor = 'DRA.01.05'

cols_to_use = [outcome, predictor] + [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]

data_subset = All_No_Missing[cols_to_use].dropna()

# Define design matrix X and target y
X = data_subset[[predictor] + covariates]
X = sm.add_constant(X)
y = data_subset[outcome]

# Count cases and controls
cases = (y == 1).sum()
controls = (y == 0).sum()

# Check if predictor is binary/categorical
unique_vals = X[predictor].nunique()

if unique_vals == 2:
    # Binary predictor - use Fisher's exact test
    contingency_table = pd.crosstab(X[predictor], y)
    print(f"Contingency table for {predictor}:")
    print(contingency_table)

    # Perform Fisher's exact test
    oddsratio, pval = fisher_exact(contingency_table)

    # Convert to log-odds (coefficient) and approximate SE
    coef = np.log(oddsratio) if oddsratio > 0 and np.isfinite(oddsratio) else np.nan

    # Approximate SE using Woolf's method for 2x2 tables
    if not contingency_table.values.min() == 0:  # No zero cells
        a, b, c, d = contingency_table.values.flatten()
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    else:
        # Add continuity correction for zero cells
        a, b, c, d = contingency_table.values.flatten() + 0.5
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

    print(f"Fisher's exact: OR={oddsratio:.3f}, p={pval:.2e}")


# In[ ]:


All_No_Missing=hla_cohort_amr.to_pandas()

# Set outcome and covariates
outcome = "condition__hpv"
covariates = [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]


# Identify predictor columns (starting from column index 30)
predictor_cols = list(All_No_Missing.columns[29:])
predictor_cols = ['DMA.01.01.01G']

# Initialize list for results
results = []

# Subset and drop rows with NA
outcome = 'condition__hpv'
predictor = 'DMA.01.01.01G'

cols_to_use = [outcome, predictor] + [
    "end_of_study_age", "sex_binary", "ancPC1", "ancPC2", "ancPC3", "ancPC4", "ancPC5",
    "ancPC6", "ancPC7", "ancPC8", "ancPC9", "ancPC10", "ancPC11", "ancPC12", "ancPC13", "ancPC14", "ancPC15",
    "ancPC16", "ancPC17", "ancPC18", "ancPC19", "ancPC20"
]

data_subset = All_No_Missing[cols_to_use].dropna()

# Define design matrix X and target y
X = data_subset[[predictor] + covariates]
X = sm.add_constant(X)
y = data_subset[outcome]

# Count cases and controls
cases = (y == 1).sum()
controls = (y == 0).sum()

# Check if predictor is binary/categorical
unique_vals = X[predictor].nunique()

if unique_vals == 2:
    # Binary predictor - use Fisher's exact test
    contingency_table = pd.crosstab(X[predictor], y)
    print(f"Contingency table for {predictor}:")
    print(contingency_table)

    # Perform Fisher's exact test
    oddsratio, pval = fisher_exact(contingency_table)

    # Convert to log-odds (coefficient) and approximate SE
    coef = np.log(oddsratio) if oddsratio > 0 and np.isfinite(oddsratio) else np.nan

    # Approximate SE using Woolf's method for 2x2 tables
    if not contingency_table.values.min() == 0:  # No zero cells
        a, b, c, d = contingency_table.values.flatten()
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    else:
        # Add continuity correction for zero cells
        a, b, c, d = contingency_table.values.flatten() + 0.5
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)

    print(f"Fisher's exact: OR={oddsratio:.3f}, p={pval:.2e}")

