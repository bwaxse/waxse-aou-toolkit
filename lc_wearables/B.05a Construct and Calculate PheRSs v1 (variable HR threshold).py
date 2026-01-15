#!/usr/bin/env python
# coding: utf-8

# # Construct and Calculate PheRSs

# ## Import Packages and Codes

# In[ ]:


import sys
print(sys.version)


# In[ ]:


from IPython.display import display, HTML
from IPython.display import clear_output
import gc
from google.cloud import bigquery
import pandas as pd
import polars as pl
import pyarrow as pa
import os
import subprocess
import numpy as np, scipy as sps
from scipy import stats
from sklearn.metrics import roc_curve, auc
import matplotlib, matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
from datetime import timedelta, datetime as dt
# import pytz
from typing import Dict, List, Optional, Union
from jinja2 import Template, Environment


# In[ ]:


from phetk.phecode import Phecode


# In[ ]:


# This line allows for the plots to be displayed inline in the Jupyter notebook
get_ipython().run_line_magic('matplotlib', 'inline')

sns.set(style="ticks",font_scale=1)


# In[ ]:


# Define the version of the Curated Data Repository (CDR).
version = get_ipython().run_line_magic('env', 'WORKSPACE_CDR')
print("version: " + version)

# get the bucket name
bucket = os.getenv('WORKSPACE_BUCKET')
print("bucket: " + bucket)


# In[ ]:


# show all columns in pandas
pd.set_option("display.max_columns", None)

# show full column width
pd.set_option('display.max_colwidth', 100)

# Polars string length to 100
pl.Config.set_fmt_str_lengths(100)

# Set the row limit to a higher value
pl.Config.set_tbl_rows(50)


# In[ ]:


palette = ['#0173b2', '#de8f05', '#8de5a1', '#d55e00', '#029e73', '#cc78bc', '#ece133', 
           '#56b4e9', '#949494', '#fbafe4', '#ca9161']

color_dict = {
    'SARS-CoV-2': '#0173b2',
    'Flu': '#de8f05',
    'RSV': '#8de5a1',
    'RV': '#d55e00',
    'hCoV': '#029e73',
    'hMPV': '#cc78bc',
    'PIV': '#ece133',
    'ADV': '#56b4e9',
    'pert': '#949494',
    'M_pna': '#fbafe4',
    'C_pna': '#ca9161'
}

sns.set_palette(sns.color_palette(palette))
sns.color_palette(palette)


# ## Helper Functions

# In[ ]:


def polars_gbq(query):
    """
    Take a SQL query and return result as polars dataframe
    :param query: BigQuery SQL query
    :return: polars dataframe
    """
    client = bigquery.Client()
    query_job = client.query(query)
    rows = query_job.result()
    df = pl.from_arrow(rows.to_arrow())

    return df


# # Load Data

# ## Load PASC Code Lists

# In[ ]:


# # copy csv file to the bucket
# args = ["gsutil", "cp", f"./Zhang_long_covid_min_1_icd.tsv", f'{bucket}/data/pasc/phers/']
# output = subprocess.run(args, capture_output=True)

# # print output from gsutil
# output.stderr

# # copy csv file to the bucket
# args = ["gsutil", "cp", f"./Zhang_long_covid_min_1_lists.tsv", f'{bucket}/data/pasc/phers/']
# output = subprocess.run(args, capture_output=True)

# # print output from gsutil
# output.stderr

# # copy csv file to the bucket
# args = ["gsutil", "cp", f"./phecodeX_unrolled_ICD_CM.csv", f'{bucket}/data/pasc/phers/']
# output = subprocess.run(args, capture_output=True)

# # print output from gsutil
# output.stderr


# In[ ]:


phecodeX_unrolled_icd_df = pl.read_csv(f'{bucket}/data/pasc/phers/phecodeX_unrolled_ICD_CM.csv')


# In[ ]:


min_1_phecodes_df = pl.read_csv(f'{bucket}/data/long_covid/phers/cox_sig_phewas_long_covid_min_1_phecodes.tsv', separator='\t')


# In[ ]:


min_1_phecodes_df = min_1_phecodes_df.select(['phecode', 'phecode_string', 
                                                          'phecode_category',  
                                                          'long_covid_min_1_hazard_ratio_low',
                                                          'long_covid_min_1_hazard_ratio_high',
                                                          'long_covid_min_1_standard_error',
                                                          'long_covid_min_1_hazard_ratio',
                                                         ])


# In[ ]:


# Get unique phecodes from each dataframe
min_1 = set(min_1_phecodes_df['phecode'].unique())


# In[ ]:


len(min_1)


# In[ ]:


zhang_pasc_icd_df = pl.read_csv(f'{bucket}/data/pasc/phers/Zhang_PASC_icd.tsv', separator='\t')


# In[ ]:


zhang_pasc_list_df = pl.read_csv(f'{bucket}/data/pasc/phers/Zhang_PASC_lists.tsv', separator='\t')


# ## Load Cohorts

# In[ ]:


combined_df = pl.read_csv(f'{bucket}/data/long_covid/combined_cohort.csv', separator='\t',
                         try_parse_dates=True)
test_df = pl.read_csv(f'{bucket}/data/long_covid/test_cohort.csv', separator='\t',
                         try_parse_dates=True)
train_df = pl.read_csv(f'{bucket}/data/long_covid/train_cohort.csv', separator='\t',
                         try_parse_dates=True)


# # Get Demographics and Residualization Covariates

# In[ ]:


demographics_q = f"""
SELECT DISTINCT
    p.person_id,
    CAST(p.birth_datetime AS DATE) AS dob,
    p_race_concept.concept_name as race,
    p_ethnicity_concept.concept_name as ethnicity,
    p_sex_at_birth_concept.concept_name as sex_at_birth,
FROM
    {version}.person p
LEFT JOIN
    {version}.concept p_race_concept 
        ON p.race_concept_id = p_race_concept.concept_id 
LEFT JOIN
    {version}.concept p_ethnicity_concept 
        ON p.ethnicity_concept_id = p_ethnicity_concept.concept_id 
LEFT JOIN
    {version}.concept p_sex_at_birth_concept 
        ON p.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id  
"""

demographics_df = polars_gbq(demographics_q)


# In[ ]:


mapping = {
    "None of these": "None",
    "I prefer not to answer": "PNA",
    "Black or African American": "Black",
    "PMI: Skip": "Skip",
    "More than one population": "Multip",
    "Asian": "Asian",
    "Native Hawaiian or Other Pacific Islander": "NHPI",
    "None Indicated": "NoInd",
    "White": "White",
    "Middle Eastern or North African": "MENA",
    "American Indian or Alaska Native": "AIAN"
}
demographics_df = demographics_df.with_columns(
    pl.col("race").map_elements(lambda x: mapping.get(x, x), return_dtype=pl.Utf8)
)


# In[ ]:


mapping = {
    "I prefer not to answer": "PNA",
    "No matching concept": "NoMatch",
    "Male": "Male",
    "Intersex": "Intersex",
    "Female": "Female",
    "PMI: Skip": "Skip",
    "None": "None"
}
demographics_df = demographics_df.with_columns(
    pl.col("sex_at_birth").map_elements(lambda x: mapping.get(x, x), return_dtype=pl.Utf8)
)


# In[ ]:


mapping = {
    "PMI: Prefer Not To Answer": "PNA",
    "No matching concept": "NoMatch",
    "What Race Ethnicity: Race Ethnicity None Of These": "None",
    "Not Hispanic or Latino": "NotHisp",
    "PMI: Skip": "Skip",
    "Hispanic or Latino": "HispLat",
}
demographics_df = demographics_df.with_columns(
    pl.col("ethnicity").map_elements(lambda x: mapping.get(x, x), return_dtype=pl.Utf8)
)


# In[ ]:


def create_demographics_dummies(df):
    """Create race/ethnicity categories for PheWAS"""
    
    # Create the race_ethnicity column using nested when-then logic
    df_with_categories = df.with_columns(
        pl.when(pl.col("ethnicity") == "HispLat")
        .then(pl.lit("Hispanic_Latino"))
        .when((pl.col("ethnicity") == "NotHisp") & (pl.col("race") == "White"))
        .then(pl.lit("White_NonHispanic"))
        .when((pl.col("ethnicity") == "NotHisp") & (pl.col("race") == "Black"))
        .then(pl.lit("Black_NonHispanic"))
        .when((pl.col("ethnicity") == "NotHisp") & (pl.col("race") == "Asian"))
        .then(pl.lit("Asian_NonHispanic"))
        .when((pl.col("ethnicity") == "NotHisp") & (pl.col("race") == "Multip"))
        .then(pl.lit("Multiracial_NonHispanic"))
        .when((pl.col("ethnicity") == "NotHisp") & (pl.col("race") == "AIAN"))
        .then(pl.lit("AIAN_NonHispanic"))
        .when((pl.col("ethnicity") == "NotHisp") & (pl.col("race") == "MENA"))
        .then(pl.lit("MENA_NonHispanic"))
        .when((pl.col("ethnicity") == "NotHisp") & (pl.col("race") == "NHPI"))
        .then(pl.lit("NHPI_NonHispanic"))
        .otherwise(pl.lit("Other_Unknown"))
        .alias("race_ethnicity")
    )
    
    # Create dummy variables (omit White_NonHispanic as reference)
    demographics_final = df_with_categories.with_columns([
        pl.when(pl.col("sex_at_birth") == "Male")
          .then(1)
          .when(pl.col("sex_at_birth") == "Female") 
          .then(0)
          .otherwise(None)  # Exclude ambiguous cases
          .cast(pl.Int8)
          .alias("sex_male"),
        (pl.col("race_ethnicity") == "Hispanic_Latino").cast(pl.Int8).alias("Hispanic_Latino"),
        (pl.col("race_ethnicity") == "Black_NonHispanic").cast(pl.Int8).alias("Black_NonHispanic"),
        (pl.col("race_ethnicity") == "Asian_NonHispanic").cast(pl.Int8).alias("Asian_NonHispanic"),
        (pl.col("race_ethnicity") == "Multiracial_NonHispanic").cast(pl.Int8).alias("Multiracial_NonHispanic"),
        (pl.col("race_ethnicity") == "Other_Unknown").cast(pl.Int8).alias("Other_Race_NonHispanic")
        # Note: White_NonHispanic omitted as reference category
    ])

    return demographics_final.drop('race_ethnicity')


# In[ ]:


demographics_with_dummies = create_demographics_dummies(demographics_df)


# # Identify ICD Codes for PheRS

# In[ ]:


icd_10_df = phecodeX_unrolled_icd_df.filter(pl.col('vocabulary_id')=='ICD10CM').select(['ICD']).unique()


# In[ ]:


zhang_137_pasc_df = zhang_pasc_list_df.filter(pl.col('137_pasc_category').is_not_null()).select(['137_pasc_category', '137_pasc_domain'])
zhang_44_pasc_df = zhang_pasc_list_df.filter(pl.col('44_pasc_category').is_not_null()).select(['44_pasc_category', '44_pasc_domain'])


# In[ ]:


zhang_137_pasc_df = zhang_137_pasc_df.join(
    zhang_pasc_icd_df.select(['ccsr_description', 'icd_10', 'icd_10_string', 'cluster_category']),
    left_on='137_pasc_category',
    right_on='ccsr_description',
    how='left'
)


# In[ ]:


zhang_44_pasc_df = zhang_44_pasc_df.join(
    zhang_pasc_icd_df.select(['ccsr_description', 'icd_10', 'icd_10_string', 'cluster_category']),
    left_on='44_pasc_category',
    right_on='ccsr_description',
    how='left'
)


# In[ ]:


def add_period_to_icd10(icd_code):
    """
    Adds a period to an ICD-10 code in the correct position.
    Standard format is a letter followed by 2 digits (or letter-digit-letter), 
    then a period, then additional digits.
    """
    if not icd_code or len(icd_code) < 3:
        return icd_code
    
    # Insert period after the third character (e.g., M16.12)
    return f"{icd_code[:3]}.{icd_code[3:]}"

# Create a mapping from no-period to with-period ICD-10 codes
icd_mapping = {}

# Option 1: Use a regex pattern to identify where periods should be in icd_10_df
with_period_codes = icd_10_df.select("ICD").to_series().to_list()

# Option 2: If icd_10_df already has the periods, create a mapping by removing periods
for code in with_period_codes:
    # Create a version with period removed
    no_period_code = code.replace(".", "")
    icd_mapping[no_period_code] = code


# In[ ]:


# update zhang_137_pasc_df
zhang_137_pasc_df = zhang_137_pasc_df.with_columns([
    pl.col("icd_10").map_elements(
        lambda x: icd_mapping.get(x, add_period_to_icd10(x)),
        return_dtype = pl.String
    ).alias("icd_10_with_period")
])

# replace the original column
zhang_137_pasc_df = zhang_137_pasc_df.drop("icd_10").rename({"icd_10_with_period": "icd_10"})


# In[ ]:


zhang_137_pasc_df = zhang_137_pasc_df.with_columns(
    pl.lit('ICD10CM').alias('vocabulary_id')
)


# In[ ]:


# update zhang_44_pasc_df
zhang_44_pasc_df = zhang_44_pasc_df.with_columns([
    pl.col("icd_10").map_elements(
        lambda x: icd_mapping.get(x, add_period_to_icd10(x)),
        return_dtype = pl.String
    ).alias("icd_10_with_period")
])

# replace the original column
zhang_44_pasc_df = zhang_44_pasc_df.drop("icd_10").rename({"icd_10_with_period": "icd_10"})


# In[ ]:


zhang_44_pasc_df = zhang_44_pasc_df.with_columns(
    pl.lit('ICD10CM').alias('vocabulary_id')
)


# ## Variable HR Threshold for Cox PheWAS

# In[ ]:


sns.histplot(min_1_phecodes_df['long_covid_min_1_hazard_ratio'], bins=100)


# In[ ]:


filtered_min_1_phecodes_df = min_1_phecodes_df.filter(pl.col('long_covid_min_1_hazard_ratio')<10)


# In[ ]:


sns.histplot(filtered_min_1_phecodes_df['long_covid_min_1_hazard_ratio'], bins=100)

def create_serial_hr_threshold(significant_phecodes_df, hr_col):
    """
    Create PheRS with incrementally stricter HR thresholds
    """
    
    # Sort phecodes by HR (ascending)
    sorted_phecodes = significant_phecodes_df.sort(hr_col)
    
    # Define removal percentiles
    removal_percentiles = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95]
    
    results = {}
    
    for percentile in removal_percentiles:
        n_total = len(sorted_phecodes)
        n_keep = int(n_total * (100 - percentile) / 100)
        
        filtered_phecodes = sorted_phecodes.tail(n_keep) if n_keep > 0 else sorted_phecodes.head(0)
        
        if len(filtered_phecodes) == 0:
            print(f"Percentile {percentile}: No phecodes remaining")
            continue
            
        min_hr = filtered_phecodes[hr_col].min()
        max_hr = filtered_phecodes[hr_col].max()
        
        print(f"Percentile {percentile}: Keeping {n_keep}/{n_total} phecodes (HR range: {min_hr:.2f} - {max_hr:.2f})")
        
        results[percentile] = filtered_phecodes
    
    return results
# In[ ]:


def create_serial_hr_thresholds(significant_phecodes_df, hr_col, 
                               focus_range=(24, 46), step_by_one=True):
    """
    Create PheRS with incrementally stricter HR thresholds.
    Focus on fine-grained steps in the critical range where performance peaks.
    
    Parameters:
    - significant_phecodes_df: DataFrame with phecodes and HR values
    - hr_col: Column name containing hazard ratios
    - focus_range: Tuple of (min_keep, max_keep) for fine-grained stepping
    - step_by_one: If True, step by 1 phecode in focus range, else use percentiles
    """
    
    # Sort phecodes by HR (ascending - lowest HR first)
    sorted_phecodes = significant_phecodes_df.sort(hr_col)
    n_total = len(sorted_phecodes)
    
    # Define thresholds
    thresholds_to_keep = []
    
    # Coarse steps outside the focus range
    coarse_steps = [n_total, int(0.9*n_total), int(0.8*n_total), int(0.7*n_total), int(0.6*n_total)]
    
    # Add steps above focus range
    for step in coarse_steps:
        if step > focus_range[1]:
            thresholds_to_keep.append(step)
    
    # Fine-grained steps in focus range
    if step_by_one:
        # Step by 1 phecode in critical range
        for n_keep in range(focus_range[1], focus_range[0] - 1, -1):
            thresholds_to_keep.append(n_keep)
    else:
        # Use percentile-based steps in focus range
        focus_steps = list(range(focus_range[0], focus_range[1] + 1, 2))
        thresholds_to_keep.extend(reversed(focus_steps))
    
    # Add steps below focus range
    coarse_steps_low = [int(0.4*n_total), int(0.3*n_total), int(0.2*n_total), int(0.1*n_total)]
    for step in coarse_steps_low:
        if step < focus_range[0]:
            thresholds_to_keep.append(step)
    
    # Remove duplicates and sort descending (most inclusive first)
    thresholds_to_keep = sorted(list(set(thresholds_to_keep)), reverse=True)
    
    results = {}
    
    for n_keep in thresholds_to_keep:
        if n_keep <= 0 or n_keep > n_total:
            continue
            
        # Take the top n_keep phecodes (highest HRs)
        filtered_phecodes = sorted_phecodes.tail(n_keep)
        
        if len(filtered_phecodes) == 0:
            print(f"Keeping {n_keep}: No phecodes remaining")
            continue
            
        min_hr = filtered_phecodes[hr_col].min()
        max_hr = filtered_phecodes[hr_col].max()
        
        print(f"Keeping {n_keep}/{n_total} phecodes (HR range: {min_hr:.2f} - {max_hr:.2f})")
        
        # Use n_keep as key for easier interpretation
        results[n_keep] = filtered_phecodes
    
    return results


# In[ ]:


phecode_thresholds = create_serial_hr_thresholds(min_1_phecodes_df, 'long_covid_min_1_hazard_ratio')


# In[ ]:


phecode_thresholds[37].sort('long_covid_min_1_hazard_ratio', descending=True)['phecode', 'phecode_string', 'phecode_category', 'long_covid_min_1_hazard_ratio']


# ## Merge in ICD Codes

# In[ ]:


for percentile in phecode_thresholds.keys():
    phecode_thresholds[percentile] = phecode_thresholds[percentile].join(
        phecodeX_unrolled_icd_df,
        on='phecode',
        how='left'
    )


# # Get Phecode Events

# In[ ]:


phecode = Phecode(platform="aou")


# In[ ]:


phecode.icd_events.head()


# ## Add ICD 9 and 10 flags to demographics

# In[ ]:


# Create binary indicators for ICD-9/10 presence and calculate EHR length
icd_indicators = (
    phecode.icd_events
    .group_by('person_id')
    .agg([
        pl.col('flag').eq(9).any().cast(pl.Int8).alias('has_icd9'),
        pl.col('flag').eq(10).any().cast(pl.Int8).alias('has_icd10'),
        (pl.col('date').max() - pl.col('date').min()).dt.total_days().truediv(365.25).alias('ehr_length_years')
    ])
)


# In[ ]:


demographics_with_dummies = demographics_with_dummies.join(
    icd_indicators,
    on='person_id',
    how='left'
)


# # Phecode Prevalence Weights (Inverse Log Normalized)

# In[ ]:


def identify_new_post_covid_icds(icd_events_df, split_df, post_covid_start=28, post_covid_end=365):
    """
    Identify ICD codes that are NEW in the post-COVID period
    Uses ALL available pre-COVID history
    """
    # Join to get time_zero for each person
    joined_df = icd_events_df.join(
        split_df.select(['person_id', 'time_zero']),
        on='person_id',
        how='inner'
    )
    
    # Calculate days since time_zero
    with_days = joined_df.with_columns([
        (pl.col('date') - pl.col('time_zero')).dt.total_days().alias('days_since_time_zero')
    ])
    
    # Split into pre-COVID and post-COVID periods
    pre_covid_events = with_days.filter(
        pl.col('days_since_time_zero') < 0  # all history before COVID
    )
    
    post_covid_events = with_days.filter(
        (pl.col('days_since_time_zero') >= post_covid_start) & 
        (pl.col('days_since_time_zero') <= post_covid_end)  # 28-365 days post-COVID
    )
    
    # Get unique person-ICD combinations from pre-COVID period
    pre_covid_icds = pre_covid_events.select(['person_id', 'ICD', 'vocabulary_id']).unique()
    
    # Mark pre-COVID ICD codes
    pre_covid_person_icds = pre_covid_icds.with_columns(
        pl.lit(True).alias('had_pre_covid')
    )
    
    # Find post-COVID ICD codes that are new (not present pre-COVID)
    new_post_covid_icds = post_covid_events.join(
        pre_covid_person_icds,
        on=['person_id', 'ICD', 'vocabulary_id'],
        how='left'
    ).filter(
        pl.col('had_pre_covid').is_null()  # Only keep ICD codes not present pre-COVID
    ).drop('had_pre_covid')
    
    print(f"Pre-COVID events: {len(pre_covid_events):,}")
    print(f"Post-COVID events: {len(post_covid_events):,}")
    print(f"New post-COVID ICD events: {len(new_post_covid_icds):,}")
    print(f"Patients with new ICD codes: {new_post_covid_icds.select(pl.n_unique('person_id')).item():,}")
    
    return new_post_covid_icds

# Apply to both cohorts
print('Training set:')
train_new_icd_events = identify_new_post_covid_icds(phecode.icd_events, train_df)
print('\nTest set:')
test_new_icd_events = identify_new_post_covid_icds(phecode.icd_events, test_df)


# ## Long COVID PheWAS Events and Weights 

# In[ ]:


cox_phewas_test_phecode_events = {}
cox_phewas_test_N = {}
cox_phewas_test_phecode_counts = {}

# Test set new events
for percentile in phecode_thresholds.keys():
    # Join ICD events with phecode mapping
    cox_phewas_test_phecode_events[percentile] = test_new_icd_events.join(
        phecode_thresholds[percentile].select(['ICD', 'vocabulary_id', 'phecode']),
        on=['ICD', 'vocabulary_id'],
        how='inner'  # Only keep events that map to a phecode
    )

    # Get total number of patients in the cohort
    cox_phewas_test_N[percentile] = test_new_icd_events.select(pl.n_unique('person_id')).item()

    # Count unique patients for each phecode
    cox_phewas_test_phecode_counts[percentile] = (cox_phewas_test_phecode_events[percentile]
        .group_by('phecode')
        .agg(
            pl.n_unique('person_id').alias('patient_count')
        )
    )


# In[ ]:


cox_phewas_train_phecode_events = {}
cox_phewas_train_N = {}
cox_phewas_train_phecode_counts = {}
cox_phewas_train_weights = {}

# Training set new events and weights
for percentile in phecode_thresholds.keys():
    # Join ICD events with phecode mapping
    cox_phewas_train_phecode_events[percentile] = train_new_icd_events.join(
        phecode_thresholds[percentile].select(['ICD', 'vocabulary_id', 'phecode']),
        on=['ICD', 'vocabulary_id'],
        how='inner'  # Only keep events that map to a phecode
    )

    # Get total number of patients in the cohort
    cox_phewas_train_N[percentile] = train_new_icd_events.select(pl.n_unique('person_id')).item()

    # Count unique patients for each phecode
    cox_phewas_train_phecode_counts[percentile] = (cox_phewas_train_phecode_events[percentile]
        .group_by('phecode')
        .agg(
            pl.n_unique('person_id').alias('patient_count')
        )
    )

    # Calculate weights as log(N/np)
    cox_phewas_train_weights[percentile] = cox_phewas_train_phecode_counts[percentile].with_columns([
        pl.lit(cox_phewas_train_N[percentile]).alias('total_patients'),
        (pl.lit(cox_phewas_train_N[percentile]) / pl.col('patient_count')).log10().alias('weight')
    ])

    # Join with phecode descriptions for better readability
    cox_phewas_train_weights[percentile] = cox_phewas_train_weights[percentile].join(
        phecode_thresholds[percentile].select(['phecode', 'phecode_string', 'phecode_category',
                                               'long_covid_min_1_hazard_ratio_low',
                                               'long_covid_min_1_hazard_ratio_high', 
                                               'long_covid_min_1_standard_error',
                                               'long_covid_min_1_hazard_ratio',]).unique(),
        on='phecode',
        how='left'
    )

    cox_phewas_train_weights[percentile].write_csv(f'{bucket}/data/long_covid/phers/variable_hr/long_covid_min_1_pct_{percentile}_phewas_train_phecode_weights.tsv', separator='\t')


# ## Zhang Events and Weights 

# In[ ]:


# Zhang test events
# Join ICD events with phecode mapping
zhang_137_test_phecode_events_df = test_new_icd_events.join(
    zhang_137_pasc_df.select(['icd_10', 'vocabulary_id', '137_pasc_category']),
    left_on=['ICD', 'vocabulary_id'],
    right_on=['icd_10', 'vocabulary_id'],
    how='inner'  # Only keep events that map to a phecode
)

# Get total number of patients in the cohort
zhang_137_test_N = test_new_icd_events.select(pl.n_unique('person_id')).item()

# Count unique patients for each phecode
zhang_137_test_phecode_counts = (zhang_137_test_phecode_events_df
    .group_by('137_pasc_category')
    .agg(
        pl.n_unique('person_id').alias('patient_count')
    )
)


# In[ ]:


# Zhang training events and weights
# Join ICD events with phecode mapping
zhang_137_train_phecode_events_df = train_new_icd_events.join(
    zhang_137_pasc_df.select(['icd_10', 'vocabulary_id', '137_pasc_category']),
    left_on=['ICD', 'vocabulary_id'],
    right_on=['icd_10', 'vocabulary_id'],
    how='inner'  # Only keep events that map to a phecode
)

# Get total number of patients in the cohort
zhang_137_train_N = train_new_icd_events.select(pl.n_unique('person_id')).item()

# Count unique patients for each phecode
zhang_137_train_phecode_counts = (zhang_137_train_phecode_events_df
    .group_by('137_pasc_category')
    .agg(
        pl.n_unique('person_id').alias('patient_count')
    )
)

# Calculate weights as log(N/np)
zhang_137_train_weights_df = zhang_137_train_phecode_counts.with_columns([
    pl.lit(zhang_137_train_N).alias('total_patients'),
    (pl.lit(zhang_137_train_N) / pl.col('patient_count')).log10().alias('weight')
])

# Join with phecode descriptions for better readability
zhang_137_train_weights_df = zhang_137_train_weights_df.join(
    zhang_137_pasc_df.select(['137_pasc_category', '137_pasc_domain', 'cluster_category']).unique(),
    on='137_pasc_category',
    how='left'
)

zhang_137_train_weights_df.write_csv(f'{bucket}/data/long_covid/phers/zhang_137_train_PASC_phecode_weights.tsv', separator='\t')


# In[ ]:


# Zhang test events
# Join ICD events with phecode mapping
zhang_44_test_phecode_events_df = test_new_icd_events.join(
    zhang_44_pasc_df.select(['icd_10', 'vocabulary_id', '44_pasc_category']),
    left_on=['ICD', 'vocabulary_id'],
    right_on=['icd_10', 'vocabulary_id'],
    how='inner'  # Only keep events that map to a phecode
)

# Get total number of patients in the cohort
zhang_44_test_N = test_new_icd_events.select(pl.n_unique('person_id')).item()

# Count unique patients for each phecode
zhang_44_test_phecode_counts = (zhang_44_test_phecode_events_df
    .group_by('44_pasc_category')
    .agg(
        pl.n_unique('person_id').alias('patient_count')
    )
)


# In[ ]:


# Zhang training events and weights
# Join ICD events with phecode mapping
zhang_44_train_phecode_events_df = train_new_icd_events.join(
    zhang_44_pasc_df.select(['icd_10', 'vocabulary_id', '44_pasc_category']),
    left_on=['ICD', 'vocabulary_id'],
    right_on=['icd_10', 'vocabulary_id'],
    how='inner'  # Only keep events that map to a phecode
)

# Get total number of patients in the cohort
zhang_44_train_N = train_new_icd_events.select(pl.n_unique('person_id')).item()

# Count unique patients for each phecode
zhang_44_train_phecode_counts = (zhang_44_train_phecode_events_df
    .group_by('44_pasc_category')
    .agg(
        pl.n_unique('person_id').alias('patient_count')
    )
)

# Calculate weights as log(N/np)
zhang_44_train_weights_df = zhang_44_train_phecode_counts.with_columns([
    pl.lit(zhang_44_train_N).alias('total_patients'),
    (pl.lit(zhang_44_train_N) / pl.col('patient_count')).log10().alias('weight')
])

# Join with phecode descriptions for better readability
zhang_44_train_weights_df = zhang_44_train_weights_df.join(
    zhang_44_pasc_df.select(['44_pasc_category', '44_pasc_domain', 'cluster_category']).unique(),
    on='44_pasc_category',
    how='left'
)

zhang_44_train_weights_df.write_csv(f'{bucket}/data/long_covid/phers/zhang_44_train_PASC_phecode_weights.tsv', separator='\t')


# # Calculate PheRS

# In[ ]:


def calculate_phers_for_patients(phecode_events, phecode_weights, phenotype, label):
    # Create a patient-phenotype presence matrix (1 if patient has phenotype, 0 otherwise)
    # We only need to know if a patient has a phecode at least once
    patient_phecode_matrix = (phecode_events
        .select(['person_id', phenotype])
        .unique()
        .with_columns(pl.lit(1).alias('has_phecode'))
    )
    
    # Join with weights
    weighted_matrix = patient_phecode_matrix.join(
        phecode_weights.select([phenotype, 'weight']),
        on=phenotype,
        how='left'
    )
    
    # Sum weights for each patient to get PheRS
    phers_scores = weighted_matrix.group_by('person_id').agg(
        pl.sum('weight').alias(label)
    )
    
    return phers_scores


# In[ ]:


def create_weight_variants(phewas_train_weights_df):    
    weights_df = phewas_train_weights_df.with_columns([
        # prevalence-based weight
        pl.col('weight').alias('prevalence_weight'),
               
        # z-statistic
        (pl.col('long_covid_min_1_hazard_ratio') / pl.col('long_covid_min_1_standard_error')).alias('z_stat_weight'),
    ])
    
    return weights_df


# In[ ]:


# Create the weight variants for PheRS approach
cox_phewas_train_weights_variants = {}

for percentile in phecode_thresholds.keys():
    cox_phewas_train_weights_variants[percentile] = create_weight_variants(cox_phewas_train_weights[percentile])


# In[ ]:


all_weight_conditions = []

# Calculate PheRS for each weighting approach
for percentile in phecode_thresholds.keys():
    weight_approaches = [
        ('prevalence_weight', f'cox_phewas_{percentile}_prevalence_phers'),
        ('z_stat_weight', f'cox_phewas_{percentile}_z_stat_phers'), 
    ]
    all_weight_conditions.extend([f'cox_phewas_{percentile}_prevalence_phers', f'cox_phewas_{percentile}_z_stat_phers'])

    for weight_col, phers_label in weight_approaches:
        # Create temporary weights dataframe for this approach
        temp_weights = cox_phewas_train_weights_variants[percentile].select(['phecode', weight_col]).rename({weight_col: 'weight'})

        # Calculate training PheRS
        temp_train_phers = calculate_phers_for_patients(cox_phewas_train_phecode_events[percentile], 
                                                 temp_weights, 
                                                 phenotype='phecode',
                                                 label=phers_label)

        # Calculate test PheRS  
        temp_test_phers = calculate_phers_for_patients(cox_phewas_test_phecode_events[percentile], 
                                                temp_weights, 
                                                phenotype='phecode',
                                                label=phers_label)

        # Join to main dataframes
        train_df = train_df.join(temp_train_phers, on='person_id', how='left').with_columns(
            pl.col(phers_label).fill_null(0)
        )

        test_df = test_df.join(temp_test_phers, on='person_id', how='left').with_columns(
            pl.col(phers_label).fill_null(0)
        )


# In[ ]:


# Calculate PheRS for each patient
zhang_137_train_phers_df = calculate_phers_for_patients(zhang_137_train_phecode_events_df, 
                                               zhang_137_train_weights_df, 
                                               phenotype='137_pasc_category',
                                               label='zhang_137_phers')

zhang_137_test_phers_df = calculate_phers_for_patients(zhang_137_test_phecode_events_df, 
                                               zhang_137_train_weights_df, 
                                               phenotype='137_pasc_category',
                                               label='zhang_137_phers')


# In[ ]:


# Calculate PheRS for each patient
zhang_44_train_phers_df = calculate_phers_for_patients(zhang_44_train_phecode_events_df, 
                                               zhang_44_train_weights_df, 
                                               phenotype='44_pasc_category',
                                               label='zhang_44_phers')

zhang_44_test_phers_df = calculate_phers_for_patients(zhang_44_test_phecode_events_df, 
                                               zhang_44_train_weights_df, 
                                               phenotype='44_pasc_category',
                                               label='zhang_44_phers')


# In[ ]:


train_df = train_df.join(
    zhang_137_train_phers_df,
    on='person_id',
    how='left'
).with_columns(
    pl.col('zhang_137_phers').fill_null(0)
)
train_df = train_df.join(
    zhang_44_train_phers_df,
    on='person_id',
    how='left'
).with_columns(
    pl.col('zhang_44_phers').fill_null(0)
)


# In[ ]:


test_df = test_df.join(
    zhang_137_test_phers_df,
    on='person_id',
    how='left'
).with_columns(
    pl.col('zhang_137_phers').fill_null(0)
)
test_df = test_df.join(
    zhang_44_test_phers_df,
    on='person_id',
    how='left'
).with_columns(
    pl.col('zhang_44_phers').fill_null(0)
)


# ## 0 Scores

# In[ ]:


all_weight_conditions.extend(['zhang_137_phers', 'zhang_44_phers'])

for col in all_weight_conditions:
    zero_pct = (train_df.filter(pl.col(col) == 0).height / train_df.height) * 100
    print(f"{col}: {zero_pct:.1f}% have score = 0")


# # Residualize PheRS Using Train Set

# In[ ]:


from sklearn.preprocessing import SplineTransformer
from sklearn.linear_model import LinearRegression


# ## Add Demographics and Residualization Variables

# In[ ]:


train_df = train_df.join(
    demographics_with_dummies.select(['person_id', 'dob', 'sex_male', 'Hispanic_Latino',
                                      'Black_NonHispanic', 'Asian_NonHispanic', 'Multiracial_NonHispanic', 
                                      'Other_Race_NonHispanic', 'has_icd9', 'has_icd10', 'ehr_length_years']),
    on='person_id',
    how='left'
)
test_df = test_df.join(
    demographics_with_dummies.select(['person_id', 'dob', 'sex_male', 'Hispanic_Latino',
                                      'Black_NonHispanic', 'Asian_NonHispanic', 'Multiracial_NonHispanic', 
                                      'Other_Race_NonHispanic', 'has_icd9', 'has_icd10', 'ehr_length_years']),
    on='person_id',
    how='left'
)


# In[ ]:


train_df = train_df.with_columns([
    (pl.col("time_zero").dt.date() - pl.col("dob")).dt.total_days()
    .truediv(365.25)
    .alias("age_years")
])
test_df = test_df.with_columns([
    (pl.col("time_zero").dt.date() - pl.col("dob")).dt.total_days()
    .truediv(365.25)
    .alias("age_years")
])


# In[ ]:


# Check which columns have NaNs
nan_check = train_df.select([
    pl.col('age_years').is_null().sum().alias('age_years_nulls'),
    pl.col('sex_male').is_null().sum().alias('sex_male_nulls'),
    pl.col('Hispanic_Latino').is_null().sum().alias('Hispanic_Latino_nulls'),
    pl.col('Black_NonHispanic').is_null().sum().alias('Black_NonHispanic_nulls'),
    pl.col('Asian_NonHispanic').is_null().sum().alias('Asian_NonHispanic_nulls'),
    pl.col('Multiracial_NonHispanic').is_null().sum().alias('Multiracial_NonHispanic_nulls'),
    pl.col('Other_Race_NonHispanic').is_null().sum().alias('Other_Race_NonHispanic_nulls'),
])

nan_check


# In[ ]:


print("Removing missing values for residualization...")
original_len = len(train_df)

print(f"Train N: {original_len}")

# Filter out people with missing sex
train_df = train_df.filter(pl.col('sex_male').is_not_null())

print(f"After excluding {original_len-len(train_df)} with missing sex: {len(train_df)}")


# In[ ]:


print("Removing missing values for residualization...")
original_len = len(test_df)

print(f"Test N: {original_len}")

# Filter out people with missing sex
test_df = test_df.filter(pl.col('sex_male').is_not_null())

print(f"After excluding {original_len-len(test_df)} with missing sex: {len(test_df)}")


# In[ ]:


all_weight_conditions


# In[ ]:


def residualize_phers_all_cases(train_df, test_df, phers_cols=all_weight_conditions):
    """
    Residualization including cases with PheRS of 0
    """

    residualized_scores = {}
    models = {}

    for phers_col in phers_cols:
        print(f"Residualizing {phers_col}...")

        # 1. FIT MODEL ON TRAINING SET ONLY
        X_train, age_spline, feature_cols = prepare_residualization_features(train_df)
        train_pd = train_df.to_pandas()
        train_mask = train_pd[phers_col].notna()

        y_train = train_pd[phers_col][train_mask].values
        X_train_clean = X_train[train_mask]

        # Fit residualization model on training data
        resid_model = LinearRegression()
        resid_model.fit(X_train_clean, y_train)

        # Calculate training residuals
        y_train_pred = resid_model.predict(X_train_clean)
        train_residuals = y_train - y_train_pred
        train_mse = np.mean(train_residuals**2)  # Use training MSE for standardization

        # 2. APPLY MODEL TO TRAINING SET
        train_studentized = train_residuals / np.sqrt(train_mse)

        # 3. APPLY SAME MODEL TO TEST SET
        X_test, _, _ = prepare_residualization_features_with_fitted_spline(test_df, age_spline, feature_cols)
        test_pd = test_df.to_pandas()
        test_mask = test_pd[phers_col].notna()

        y_test = test_pd[phers_col][test_mask].values
        X_test_clean = X_test[test_mask]

        # Apply training model to test data
        y_test_pred = resid_model.predict(X_test_clean)
        test_residuals = y_test - y_test_pred
        test_studentized = test_residuals / np.sqrt(train_mse)  # Use training MSE for standardization

        # Store results
        residualized_scores[f'train_{phers_col}'] = train_studentized
        residualized_scores[f'test_{phers_col}'] = test_studentized
        models[phers_col] = (resid_model, age_spline, train_mse)

        print(f"  Training - Original: mean={y_train.mean():.3f}, Residualized: mean={train_studentized.mean():.3f}")
        print(f"  Test - Original: mean={y_test.mean():.3f}, Residualized: mean={test_studentized.mean():.3f}")

    return residualized_scores, models


# In[ ]:


def prepare_residualization_features(df):
    """
    Prepare features for the residualization model
    """
    # Handle both pandas and polars input
    if hasattr(df, 'to_pandas'):  # It's a polars DataFrame
        df_pd = df.to_pandas()
    else:  # It's already pandas
        df_pd = df
    
    # Age spline (cubic with 3 knots)
    age_spline = SplineTransformer(n_knots=3, degree=3, include_bias=False)
    age_features = age_spline.fit_transform(df_pd[['age_years']].values)
    
    # Collect all features
    feature_cols = [
        'sex_male',
        'Hispanic_Latino', 'Black_NonHispanic', 'Asian_NonHispanic', 
        'Multiracial_NonHispanic', 'Other_Race_NonHispanic', 
        'has_icd9', 'has_icd10', 'ehr_length_years'
    ]
        
    # Stack features
    demographic_features = df_pd[feature_cols].values
    X = np.hstack([age_features, demographic_features])
    
    return X, age_spline, feature_cols


# In[ ]:


def prepare_residualization_features_with_fitted_spline(df, fitted_age_spline, feature_cols):
    """
    Apply already-fitted age spline to new data
    """
    df_pd = df.to_pandas()
    
    # Use the already-fitted spline
    age_features = fitted_age_spline.transform(df_pd[['age_years']].values)
    
    # Get demographic features
    demographic_features = df_pd[feature_cols].values
    X = np.hstack([age_features, demographic_features])
    
    return X, fitted_age_spline, feature_cols


# In[ ]:


# Add residualized scores back to your dataframes
def add_residualized_scores_to_dfs(train_df, test_df, residualized_scores):
    """
    Add residualized scores back to train and test dataframes
    """
    # Convert to pandas for easier manipulation
    train_pd = train_df.to_pandas()
    test_pd = test_df.to_pandas()
    
    # Add training residualized scores
    for score_name, scores in residualized_scores.items():
        if score_name.startswith('train_'):
            col_name = score_name.replace('train_', 'r')  # e.g., 'train_phewas_phers' -> 'rphewas_phers'
            
            # Create full array with NaNs
            full_scores = np.full(len(train_pd), np.nan)
            
            # Find which rows had valid scores
            original_col = score_name.replace('train_', '')  # e.g., 'train_phewas_phers' -> 'phewas_phers'
            mask = train_pd[original_col].notna()
            full_scores[mask] = scores
            
            train_pd[col_name] = full_scores
    
    # Add test residualized scores  
    for score_name, scores in residualized_scores.items():
        if score_name.startswith('test_'):
            col_name = score_name.replace('test_', 'r')  # e.g., 'test_phewas_phers' -> 'rphewas_phers'
            
            # Create full array with NaNs
            full_scores = np.full(len(test_pd), np.nan)
            
            # Find which rows had valid scores
            original_col = score_name.replace('test_', '')  # e.g., 'test_phewas_phers' -> 'phewas_phers'
            mask = test_pd[original_col].notna()
            full_scores[mask] = scores
            
            test_pd[col_name] = full_scores
    
    # Convert back to polars
    train_final = pl.from_pandas(train_pd)
    test_final = pl.from_pandas(test_pd)
    
    return train_final, test_final


# - Zero PheRS = "No qualifying high-risk phenotypes"
# - Negative rPheRS = "Even healthier than demographics would predict"
# - Positive rPheRS = "Higher burden than expected for demographic group"

# In[ ]:


# Step 1: Fix the function and apply residualization
residualized_scores, resid_models = residualize_phers_all_cases(train_df, test_df)


# In[ ]:


# Step 2: Add scores back to dataframes
train_final, test_final = add_residualized_scores_to_dfs(train_df, test_df, residualized_scores)


# In[ ]:


train_final.write_csv(f'{bucket}/data/long_covid/phers/variable_hr/residualized_phers_train_plus_set.tsv', 
                      separator='\t')
test_final.write_csv(f'{bucket}/data/long_covid/phers/variable_hr/residualized_phers_test_plus_set.tsv', 
                      separator='\t')

