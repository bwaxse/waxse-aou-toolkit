#!/usr/bin/env python
# coding: utf-8

# # Apply Residualized Cox PheWAS 37 PheRS to a New Cohort
#
# **Purpose**: Apply the `r_phers` score to a new participant cohort
#
# Steps:
# 1. Extract all ICD diagnosis (based on All of Us OMOP CDM)
# 2. Identify NEW post-index-event diagnoses (incident phenotypes)
# 3. Calculate raw PheRS using phecode weights
# 4. Residualize PheRS using participant demographics and EHR characteristics
# 5. Export results
#
# ## Section 1: Introduction and Setup

# In[ ]:


import sys
print(sys.version)


# In[ ]:


from IPython.display import display, HTML
from google.cloud import bigquery
import pandas as pd
import polars as pl
import pyarrow as pa
import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple
from sklearn.preprocessing import SplineTransformer
from sklearn.linear_model import LinearRegression


# In[ ]:


# Set up display options
pd.set_option("display.max_columns", None)
pd.set_option('display.max_colwidth', 100)
pl.Config.set_fmt_str_lengths(100)
pl.Config.set_tbl_rows(50)


# In[ ]:


# Get workspace environment variables
WORKSPACE_CDR = os.environ['WORKSPACE_CDR']
WORKSPACE_BUCKET = os.environ['WORKSPACE_BUCKET']


# In[ ]:


print("""
================================================================================
PheRS BACKGROUND AND INTERPRETATION
================================================================================

- A PheRS (Phenotype Risk Score) is a composite measure of disease burden based 
    on ICD diagnostic codes
- This long COVID PheRS captures a specific pattern of incident phenotypes/diagnoses
  following COVID infection

This Specific Residualized Score: r_phers
- 37 phecodes selected from a Cox PheWAS
- Weighted by z-statistics (Effect Size / Std Error) from Cox models
- Adjusted for age, sex, race/ethnicity, and EHR coverage characteristics
- "r" prefix = residualized (adjusted for demographic confounders)

Interpretation:
- 0 = Expected qualifying phenotypes detected in post-index period
- Negative rPheRS = Fewer phenotypes than expected for this demographic group
- Positive rPheRS = More phenotypes than expected for this demographic group

Event Definition:
- Includes ICD diagnoses 28-365 days post-index-event
- 28-day lag excludes acute/immediate period
- 365-day window targets chronic sequelae (long-term conditions)
- Requires "new" diagnoses (not present in pre-index history)

================================================================================
""")


# ## Section 2: Load Reference Data

# In[ ]:


# Load phecode weights
phecode_weights_df = pl.read_csv('long_covid_phecode_weights.tsv', separator='\t') #############################
print(f"✓ Loaded {len(phecode_weights_df)} phecodes with weights")
print(f"  Columns: {phecode_weights_df.columns}")


# In[ ]:


# Load phecode-ICD mappings
phecode_icd_df = pl.read_csv(phecodeX_unrolled_ICD_CM.csv') #############################
print(f"✓ Loaded {len(phecode_icd_df)} phecode-ICD mappings")
print(f"  Columns: {phecode_icd_df.columns}")


# In[ ]:


# Display summary of phecode weights
print("\n--- Phecode Weight Summary ---")
display(phecode_weights_df.select(['phecode', 'phecode_string', 'z_stat_weight']).head(10))


# ## Section 3: Load User Cohort

# In[ ]:


print("\n=== LOADING USER COHORT ===\n")

# IMPORTANT: Users must modify the path below to their cohort csv file
user_cohort = pl.read_csv(COHORT_FILE_PATH, try_parse_dates=True)
print(f"✓ Loaded cohort with {len(user_cohort)} participants")


# In[ ]:


# Validate required columns
required_cols = ['person_id', 'time_zero']
missing_cols = [col for col in required_cols if col not in user_cohort.columns]

if missing_cols:
    raise ValueError(
        f"Cohort file must contain columns: {required_cols}\n"
        f"Missing columns: {missing_cols}\n"
        f"Available columns: {user_cohort.columns}"
    )

print(f"\n✓ Required columns present: {required_cols}")
print(f"\nCohort preview:")
display(user_cohort.head())

print(f"\nCohort summary:")
print(f"  N participants: {len(user_cohort):,}")
print(f"  Date range for time_zero: {user_cohort['time_zero'].min()} to {user_cohort['time_zero'].max()}")


# ## Section 4: Extract and Process ICD Events

# In[ ]:


def phetk_icd_query(ds: str) -> str:
    """
    Generate SQL query to extract ICD codes from All of Us OMOP data.

    Extracts ICD9CM and ICD10CM codes from both condition_occurrence and
    observation tables with special V-code handling via concept_relationship.

    Args:
        ds (str): Google BigQuery dataset ID (WORKSPACE_CDR)

    Returns:
        str: SQL query string that generates table with:
            - person_id: Participant identifier
            - date: Diagnosis date
            - ICD: ICD code (concept_code)
            - vocabulary_id: ICD9CM or ICD10CM
    """

    # Base ICD extraction from 4 paths
    icd_query = f"""
        (
            SELECT DISTINCT
                co.person_id,
                co.condition_start_date AS date,
                c.vocabulary_id AS vocabulary_id,
                c.concept_code AS ICD,
                co.condition_concept_id AS concept_id
            FROM
                {ds}.condition_occurrence AS co
            INNER JOIN
                {ds}.concept AS c
            ON
                co.condition_source_value = c.concept_code
            WHERE
                c.vocabulary_id in ("ICD9CM", "ICD10CM")
        )
        UNION DISTINCT
        (
            SELECT DISTINCT
                co.person_id,
                co.condition_start_date AS date,
                c.vocabulary_id AS vocabulary_id,
                c.concept_code AS ICD,
                co.condition_concept_id AS concept_id
            FROM
                {ds}.condition_occurrence AS co
            INNER JOIN
                {ds}.concept AS c
            ON
                co.condition_source_concept_id = c.concept_id
            WHERE
                c.vocabulary_id in ("ICD9CM", "ICD10CM")
        )
        UNION DISTINCT
        (
            SELECT DISTINCT
                o.person_id,
                o.observation_date AS date,
                c.vocabulary_id AS vocabulary_id,
                c.concept_code AS ICD,
                o.observation_concept_id AS concept_id
            FROM
                {ds}.observation AS o
            INNER JOIN
                {ds}.concept as c
            ON
                o.observation_source_value = c.concept_code
            WHERE
                c.vocabulary_id in ("ICD9CM", "ICD10CM")
        )
        UNION DISTINCT
        (
            SELECT DISTINCT
                o.person_id,
                o.observation_date AS date,
                c.vocabulary_id AS vocabulary_id,
                c.concept_code AS ICD,
                o.observation_concept_id AS concept_id
            FROM
                {ds}.observation AS o
            INNER JOIN
                {ds}.concept as c
            ON
                o.observation_source_concept_id = c.concept_id
            WHERE
                c.vocabulary_id in ("ICD9CM", "ICD10CM")
        )
    """

    # V-code handling (V codes overlap between ICD9CM and ICD10CM)
    v_icd_vocab_query = f"""
        SELECT DISTINCT
            v_icds.person_id,
            v_icds.date,
            v_icds.ICD,
            c.vocabulary_id
        FROM
            (
                SELECT
                    *
                FROM
                    ({icd_query}) AS icd_events
                WHERE
                    icd_events.ICD LIKE "V%"
            ) AS v_icds
        INNER JOIN
            {ds}.concept_relationship AS cr
        ON
            v_icds.concept_id = cr.concept_id_1
        INNER JOIN
            {ds}.concept AS c
        ON
            cr.concept_id_2 = c.concept_id
        WHERE
            c.vocabulary_id IN ("ICD9CM", "ICD10CM")
        AND
            v_icds.ICD = c.concept_code
        AND NOT
            v_icds.vocabulary_id != c.vocabulary_id
    """

    # Final query: combine non-V codes with corrected V codes
    final_query = f"""
        (
            SELECT DISTINCT
                person_id,
                date,
                ICD,
                vocabulary_id
            FROM
                ({icd_query})
            WHERE
                NOT ICD LIKE "V%"
        )
        UNION DISTINCT
        (
            SELECT DISTINCT
                *
            FROM
                ({v_icd_vocab_query})
        )
    """

    return final_query


# In[ ]:


def polars_gbq(query: str) -> pl.DataFrame:
    """
    Execute BigQuery SQL query and return Polars DataFrame.

    Args:
        query (str): SQL query string

    Returns:
        pl.DataFrame: Query results
    """
    client = bigquery.Client()
    query_job = client.query(query)
    rows = query_job.result()
    df = pl.from_arrow(rows.to_arrow())
    return df


# In[ ]:


print("\n=== EXTRACTING ICD EVENTS ===\n")

# Get the ICD extraction query
icd_query = phetk_icd_query(WORKSPACE_CDR)

# Get cohort person IDs as a string for SQL WHERE clause
person_ids_str = ",".join(str(pid) for pid in user_cohort['person_id'].to_list())

# Execute query to get all ICD codes for cohort
full_icd_query = f"""
    SELECT person_id, date, ICD, vocabulary_id
    FROM ({icd_query})
    WHERE person_id IN ({person_ids_str})
    AND date IS NOT NULL
"""

print("Querying ICD events from condition_occurrence and observation tables...")
print("(This may take 1-2 minutes depending on cohort size)")

icd_events = polars_gbq(full_icd_query)

print(f"✓ Extracted {len(icd_events):,} ICD events")
print(f"  Unique participants with ICD codes: {icd_events.select(pl.n_unique('person_id')).item():,}")
print(f"  Date range: {icd_events['date'].min()} to {icd_events['date'].max()}")


# In[ ]:


# Show vocabulary distribution
vocab_counts = icd_events.group_by('vocabulary_id').agg(
    pl.n_unique('person_id').alias('n_persons'),
    pl.count().alias('n_events')
).sort('n_events', descending=True)

print("\nVocabulary distribution:")
display(vocab_counts)


# In[ ]:


def identify_new_post_index_icds(icd_events_df: pl.DataFrame,
                                 cohort_df: pl.DataFrame,
                                 post_start_days: int = 28,
                                 post_end_days: int = 365) -> pl.DataFrame:
    """
    Identify ICD codes that are NEW (incident) in the post-index-event period.

    An ICD code is considered "NEW" if it appears in the post-index period
    but was NOT present in the pre-index period (baseline history).

    Args:
        icd_events_df: DataFrame with person_id, date, ICD, vocabulary_id
        cohort_df: DataFrame with person_id, time_zero
        post_start_days: Days after time_zero to start counting (default 28)
        post_end_days: Days after time_zero to stop counting (default 365)

    Returns:
        DataFrame with NEW post-index ICD events only
    """

    # Join to get time_zero for each person
    joined_df = icd_events_df.join(
        cohort_df.select(['person_id', 'time_zero']),
        on='person_id',
        how='inner'
    )

    # Calculate days since time_zero
    with_days = joined_df.with_columns([
        (pl.col('date') - pl.col('time_zero')).dt.total_days().alias('days_since_time_zero')
    ])

    # Split into pre-index and post-index periods
    pre_index_events = with_days.filter(
        pl.col('days_since_time_zero') < 0  # All history before index
    )

    post_index_events = with_days.filter(
        (pl.col('days_since_time_zero') >= post_start_days) &
        (pl.col('days_since_time_zero') <= post_end_days)
    )

    # Get unique person-ICD combinations from pre-index period
    pre_index_icds = pre_index_events.select(['person_id', 'ICD', 'vocabulary_id']).unique()

    # Mark pre-index ICD codes
    pre_index_person_icds = pre_index_icds.with_columns(
        pl.lit(True).alias('had_pre_index')
    )

    # Find post-index ICD codes that are new (not present pre-index)
    new_post_index_icds = post_index_events.join(
        pre_index_person_icds,
        on=['person_id', 'ICD', 'vocabulary_id'],
        how='left'
    ).filter(
        pl.col('had_pre_index').is_null()  # Only keep ICD codes not present pre-index
    ).drop('had_pre_index')

    n_pre = len(pre_index_events)
    n_post = len(post_index_events)
    n_new = len(new_post_index_icds)
    n_persons_with_new = new_post_index_icds.select(pl.n_unique('person_id')).item()

    print(f"\nIncident ICD Detection Summary:")
    print(f"  Pre-index events: {n_pre:,}")
    print(f"  Post-index events ({post_start_days}-{post_end_days} days): {n_post:,}")
    print(f"  NEW post-index events: {n_new:,}")
    print(f"  Participants with ≥1 new ICD: {n_persons_with_new:,} / {len(cohort_df):,}")

    return new_post_index_icds


# In[ ]:


print("\n=== IDENTIFYING NEW POST-INDEX DIAGNOSES ===\n")

new_icd_events = identify_new_post_index_icds(
    icd_events,
    user_cohort,
    post_start_days=28,
    post_end_days=365
)


# In[ ]:


# Filter phecode-ICD mapping to only the 37 target phecodes
target_phecodes = set(phecode_weights_df['phecode'].to_list())

phecode_icd_filtered = phecode_icd_df.filter(
    pl.col('phecode').is_in(target_phecodes)
)

print(f"\nPhecode-ICD Mapping:")
print(f"  Target phecodes: {len(target_phecodes)}")
print(f"  Unique ICD codes for targets: {phecode_icd_filtered.select(pl.n_unique('ICD')).item():,}")


# In[ ]:


# Map NEW ICD events to phecodes
phecode_icd_mapping = phecode_icd_filtered.select(['ICD', 'vocabulary_id', 'phecode'])

new_phecode_events = new_icd_events.join(
    phecode_icd_mapping,
    on=['ICD', 'vocabulary_id'],
    how='inner'  # Keep only ICD codes that map to target phecodes
)

print(f"\nPhecode Mapping Results:")
print(f"  NEW ICD events mapping to 37 phecodes: {len(new_phecode_events):,}")
print(f"  Unique phecodes present: {new_phecode_events.select(pl.n_unique('phecode')).item()}")
print(f"  Unique persons with mapped phecodes: {new_phecode_events.select(pl.n_unique('person_id')).item():,}")


# ## Section 5: Calculate Raw PheRS

# In[ ]:


def calculate_phers_for_patients(phecode_events: pl.DataFrame,
                                 phecode_weights: pl.DataFrame,
                                 weight_col: str = 'z_stat_weight') -> pl.DataFrame:
    """
    Calculate PheRS by summing weights of phecodes present in each participant.

    For each participant, identifies all unique phecodes they have (presence/absence)
    and sums the corresponding weights to get PheRS.

    Args:
        phecode_events: DataFrame with person_id, phecode columns
        phecode_weights: DataFrame with phecode and weight columns
        weight_col: Column name in phecode_weights containing the weight to use

    Returns:
        DataFrame with person_id and phers_score columns
    """

    # Create binary matrix: unique person-phecode pairs
    patient_phecode_presence = (phecode_events
        .select(['person_id', 'phecode'])
        .unique()
    )

    # Join with weights
    weighted_matrix = patient_phecode_presence.join(
        phecode_weights.select(['phecode', weight_col]),
        on='phecode',
        how='left'
    )

    # Sum weights for each person to get PheRS
    phers_scores = weighted_matrix.group_by('person_id').agg(
        pl.sum(weight_col).alias('phers_score')
    )

    return phers_scores


# In[ ]:


print("\n=== CALCULATING RAW PheRS ===\n")

# Select z-statistic weight for this analysis
raw_phers = calculate_phers_for_patients(
    new_phecode_events,
    phecode_weights.select(['phecode', 'z_stat_weight']),
    weight_col='z_stat_weight'
)

# Rename to match output column name
raw_phers = raw_phers.rename({'phers_score': 'cox_phewas_37_z_stat_phers'})

print(f"✓ Calculated raw PheRS for {len(raw_phers):,} participants")
print(f"\nRaw PheRS Statistics:")
phers_stats = raw_phers['cox_phewas_37_z_stat_phers'].describe()
print(f"  Min: {phers_stats[0]:.3f}")
print(f"  25th percentile: {phers_stats[1]:.3f}")
print(f"  Median: {phers_stats[2]:.3f}")
print(f"  Mean: {phers_stats[3]:.3f}")
print(f"  75th percentile: {phers_stats[4]:.3f}")
print(f"  Max: {phers_stats[5]:.3f}")


# In[ ]:


# Identify participants with zero PheRS
zero_phers_count = (raw_phers
    .filter(pl.col('cox_phewas_37_z_stat_phers') == 0)
    .height
)
zero_pct = (zero_phers_count / len(raw_phers)) * 100

print(f"\nZero PheRS Distribution:")
print(f"  Participants with PheRS = 0: {zero_phers_count:,} ({zero_pct:.1f}%)")
print(f"  Participants with PheRS > 0: {len(raw_phers) - zero_phers_count:,} ({100-zero_pct:.1f}%)")


# In[ ]:


# Plot distribution
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# Histogram
axes[0].hist(raw_phers['cox_phewas_37_z_stat_phers'].to_list(), bins=50, edgecolor='black', alpha=0.7)
axes[0].set_xlabel('Raw PheRS (z-statistic weighted)')
axes[0].set_ylabel('Number of Participants')
axes[0].set_title('Distribution of Raw PheRS')
axes[0].grid(alpha=0.3)

# Box plot
axes[1].boxplot(raw_phers['cox_phewas_37_z_stat_phers'].to_list(), vert=True)
axes[1].set_ylabel('Raw PheRS')
axes[1].set_title('Box Plot of Raw PheRS')
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.show()


# ## Section 6: Residualize PheRS

# In[ ]:


def get_demographics_and_covariates(cohort_df: pl.DataFrame,
                                    icd_events_df: pl.DataFrame) -> pl.DataFrame:
    """
    Query demographics and calculate EHR-based covariates.

    Retrieves from person table:
    - Age (calculated from DOB at time_zero)
    - Sex at birth
    - Race and ethnicity

    Calculates from ICD events:
    - ICD9/ICD10 presence flags
    - EHR length (duration from first to last ICD event)

    Args:
        cohort_df: Cohort DataFrame with time_zero
        icd_events_df: All ICD events

    Returns:
        DataFrame with covariates for residualization
    """

    # Query demographics from person table
    person_ids_str = ",".join(str(pid) for pid in cohort_df['person_id'].to_list())

    demographics_query = f"""
    SELECT DISTINCT
        p.person_id,
        CAST(p.birth_datetime AS DATE) AS dob,
        p_race_concept.concept_name as race,
        p_ethnicity_concept.concept_name as ethnicity,
        p_sex_at_birth_concept.concept_name as sex_at_birth
    FROM
        {WORKSPACE_CDR}.person p
    LEFT JOIN
        {WORKSPACE_CDR}.concept p_race_concept
            ON p.race_concept_id = p_race_concept.concept_id
    LEFT JOIN
        {WORKSPACE_CDR}.concept p_ethnicity_concept
            ON p.ethnicity_concept_id = p_ethnicity_concept.concept_id
    LEFT JOIN
        {WORKSPACE_CDR}.concept p_sex_at_birth_concept
            ON p.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id
    WHERE p.person_id IN ({person_ids_str})
    """

    print("Querying demographics from person table...")
    demographics_df = polars_gbq(demographics_query)

    # Calculate ICD9/10 indicators and EHR length
    icd_indicators = (
        icd_events_df
        .group_by('person_id')
        .agg([
            pl.col('vocabulary_id').eq('ICD9CM').any().cast(pl.Int8).alias('has_icd9'),
            pl.col('vocabulary_id').eq('ICD10CM').any().cast(pl.Int8).alias('has_icd10'),
            (pl.col('date').max() - pl.col('date').min()).dt.total_days().truediv(365.25).alias('ehr_length_years')
        ])
    )

    # Join demographics with ICD indicators
    combined = demographics_df.join(
        icd_indicators,
        on='person_id',
        how='left'
    )

    # Fill missing ICD flags with 0
    combined = combined.with_columns([
        pl.col('has_icd9').fill_null(0),
        pl.col('has_icd10').fill_null(0),
        pl.col('ehr_length_years').fill_null(0)
    ])

    return combined


# In[ ]:


def create_demographic_dummies(df: pl.DataFrame) -> pl.DataFrame:
    """
    Create dummy variables for race/ethnicity and sex.

    Race/Ethnicity Categories:
    - Hispanic_Latino (reference)
    - White_NonHispanic
    - Black_NonHispanic
    - Asian_NonHispanic
    - Multiracial_NonHispanic
    - Other_Unknown (combined all others)

    Sex:
    - sex_male (1 = Male, 0 = Female, NULL = missing)

    Args:
        df: DataFrame with race, ethnicity, sex_at_birth columns

    Returns:
        DataFrame with dummy variables
    """

    # Standardize race names
    race_mapping = {
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

    df = df.with_columns(
        pl.col("race").map_elements(lambda x: race_mapping.get(x, x), return_dtype=pl.Utf8)
    )

    # Standardize ethnicity names
    ethnicity_mapping = {
        "I prefer not to answer": "PNA",
        "No matching concept": "NoMatch",
        "Not Hispanic or Latino": "NotHisp",
        "Hispanic or Latino": "HispLat",
        "PMI: Prefer Not To Answer": "PNA",
        "PMI: Skip": "Skip",
        "What Race Ethnicity: Race Ethnicity None Of These": "None",
    }

    df = df.with_columns(
        pl.col("ethnicity").map_elements(lambda x: ethnicity_mapping.get(x, x), return_dtype=pl.Utf8)
    )

    # Standardize sex names
    sex_mapping = {
        "Male": "Male",
        "Female": "Female",
        "Intersex": "Intersex",
        "I prefer not to answer": "PNA",
        "PMI: Skip": "Skip",
        "No matching concept": "NoMatch",
        "None": "None"
    }

    df = df.with_columns(
        pl.col("sex_at_birth").map_elements(lambda x: sex_mapping.get(x, x), return_dtype=pl.Utf8)
    )

    # Create race_ethnicity combinations
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

    # Create dummy variables
    demographics_final = df_with_categories.with_columns([
        pl.when(pl.col("sex_at_birth") == "Male")
          .then(1)
          .when(pl.col("sex_at_birth") == "Female")
          .then(0)
          .otherwise(None)
          .cast(pl.Int8)
          .alias("sex_male"),
        (pl.col("race_ethnicity") == "Hispanic_Latino").cast(pl.Int8).alias("Hispanic_Latino"),
        (pl.col("race_ethnicity") == "Black_NonHispanic").cast(pl.Int8).alias("Black_NonHispanic"),
        (pl.col("race_ethnicity") == "Asian_NonHispanic").cast(pl.Int8).alias("Asian_NonHispanic"),
        (pl.col("race_ethnicity") == "Multiracial_NonHispanic").cast(pl.Int8).alias("Multiracial_NonHispanic"),
        (pl.col("race_ethnicity") == "Other_Unknown").cast(pl.Int8).alias("Other_Race_NonHispanic")
    ])

    return demographics_final.drop(['race', 'ethnicity', 'sex_at_birth', 'race_ethnicity'])


# In[ ]:


print("\n=== GETTING DEMOGRAPHICS AND COVARIATES ===\n")

demographics = get_demographics_and_covariates(user_cohort, icd_events)
demographics_dummy = create_demographic_dummies(demographics)

print(f"✓ Retrieved demographics for {len(demographics_dummy):,} participants")

# Calculate age at time_zero
cohort_with_age = user_cohort.join(
    demographics_dummy.select(['person_id', 'dob']),
    on='person_id',
    how='left'
).with_columns(
    ((pl.col("time_zero").dt.date() - pl.col("dob")).dt.total_days() / 365.25).alias("age_years")
)

demographics_dummy = demographics_dummy.join(
    cohort_with_age.select(['person_id', 'age_years']),
    on='person_id',
    how='left'
)

print(f"\nDemographics Summary:")
print(f"  Age at time_zero - Mean: {demographics_dummy['age_years'].mean():.1f}, Range: {demographics_dummy['age_years'].min():.1f}-{demographics_dummy['age_years'].max():.1f}")
print(f"  Male: {demographics_dummy.filter(pl.col('sex_male')==1).height} ({demographics_dummy.filter(pl.col('sex_male')==1).height/len(demographics_dummy)*100:.1f}%)")
print(f"  Has EHR length >0: {demographics_dummy.filter(pl.col('ehr_length_years')>0).height}")


# In[ ]:


def prepare_residualization_features(df: pl.DataFrame) -> Tuple[np.ndarray, SplineTransformer, List[str]]:
    """
    Prepare feature matrix for residualization model.

    Creates:
    - Age spline: Cubic B-spline with 3 knots (6 features after transformation)
    - Demographics: Sex, race/ethnicity dummies
    - EHR characteristics: ICD9/10 flags, EHR length

    Args:
        df: Polars DataFrame with age_years and demographic columns

    Returns:
        Tuple of (feature_matrix, fitted_age_spline, feature_column_names)
    """

    # Convert to pandas for sklearn
    df_pd = df.to_pandas()

    # Age spline (cubic, 3 knots, no bias)
    age_spline = SplineTransformer(n_knots=3, degree=3, include_bias=False)
    age_features = age_spline.fit_transform(df_pd[['age_years']].values)

    # Demographic and EHR features
    feature_cols = [
        'sex_male',
        'Hispanic_Latino', 'Black_NonHispanic', 'Asian_NonHispanic',
        'Multiracial_NonHispanic', 'Other_Race_NonHispanic',
        'has_icd9', 'has_icd10', 'ehr_length_years'
    ]

    demographic_features = df_pd[feature_cols].values

    # Stack all features
    X = np.hstack([age_features, demographic_features])

    return X, age_spline, feature_cols


# In[ ]:


def residualize_phers(cohort_df: pl.DataFrame,
                      raw_phers_df: pl.DataFrame,
                      demographics_df: pl.DataFrame) -> Tuple[pl.DataFrame, Dict]:
    """
    Fit residualization model and calculate residualized PheRS.

    Approach:
    1. Merge raw PheRS with demographics
    2. Fit LinearRegression: raw_phers ~ age_spline + sex + race_ethnicity + has_icd9/10 + ehr_length
    3. Calculate residuals and studentize by training MSE

    Args:
        cohort_df: Cohort DataFrame
        raw_phers_df: DataFrame with cox_phewas_37_z_stat_phers
        demographics_df: DataFrame with covariates

    Returns:
        Tuple of (cohort with residualized scores, model info dict)
    """

    print("\nPreparing residualization...")

    # Merge all data
    merged_df = cohort_df.join(raw_phers_df, on='person_id', how='left')
    merged_df = merged_df.join(demographics_df.select(['person_id', 'age_years', 'sex_male',
                                                       'Hispanic_Latino', 'Black_NonHispanic',
                                                       'Asian_NonHispanic', 'Multiracial_NonHispanic',
                                                       'Other_Race_NonHispanic', 'has_icd9',
                                                       'has_icd10', 'ehr_length_years']),
                              on='person_id', how='left')

    # Fill null PheRS with 0
    merged_df = merged_df.with_columns(
        pl.col('cox_phewas_37_z_stat_phers').fill_null(0)
    )

    # Check for missing covariates
    n_original = len(merged_df)
    missing_check = merged_df.select([
        pl.col('age_years').is_null().sum(),
        pl.col('sex_male').is_null().sum(),
    ])

    # Remove rows with missing sex (required for model)
    merged_df = merged_df.filter(pl.col('sex_male').is_not_null())
    n_after = len(merged_df)

    if n_original > n_after:
        print(f"  Removed {n_original - n_after} participants with missing sex")

    print(f"  Residualizing for {n_after:,} participants")

    # Prepare features for ALL data
    X, age_spline, feature_cols = prepare_residualization_features(merged_df)

    # Get target variable
    merged_pd = merged_df.to_pandas()
    y = merged_pd['cox_phewas_37_z_stat_phers'].values

    # Fit residualization model
    print("  Fitting residualization model...")
    resid_model = LinearRegression()
    resid_model.fit(X, y)

    # Calculate residuals and studentize
    y_pred = resid_model.predict(X)
    residuals = y - y_pred
    mse = np.mean(residuals**2)
    studentized_residuals = residuals / np.sqrt(mse)

    # Calculate R-squared
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)

    # Add residualized scores to dataframe
    merged_df = merged_df.with_columns(
        pl.Series('r_phers', studentized_residuals).alias('r_phers')
    )

    print(f"\nResiduation Model Results:")
    print(f"  R²: {r_squared:.4f}")
    print(f"  Mean (raw): {y.mean():.3f}")
    print(f"  Mean (residualized): {studentized_residuals.mean():.3f}")
    print(f"  SD (raw): {y.std():.3f}")
    print(f"  SD (residualized): {studentized_residuals.std():.3f}")

    model_info = {
        'age_spline': age_spline,
        'residual_mse': mse,
        'r_squared': r_squared,
        'n_samples': len(merged_df)
    }

    return merged_df, model_info


# In[ ]:


print("\n=== RESIDUALIZING PheRS ===\n")

cohort_with_phers, model_info = residualize_phers(
    user_cohort,
    raw_phers,
    demographics_dummy
)


# In[ ]:


# Check residualization quality
print("\nResidual Quality Checks:")

# Check mean is near zero
residuals = cohort_with_phers['r_phers'].to_list()
print(f"  Mean residual: {np.mean(residuals):.6f} (expected ≈ 0)")
print(f"  Median residual: {np.median(residuals):.3f}")

# Check for extreme outliers
n_outliers = sum(1 for r in residuals if abs(r) > 5)
print(f"  Outliers (|residual| > 5): {n_outliers} ({n_outliers/len(residuals)*100:.2f}%)")

# Show distribution
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# Raw vs Residualized
axes[0].hist(cohort_with_phers['cox_phewas_37_z_stat_phers'].to_list(),
             bins=50, alpha=0.6, label='Raw PheRS', edgecolor='black')
axes[0].set_xlabel('Score')
axes[0].set_ylabel('Number of Participants')
axes[0].set_title('Raw PheRS Distribution')
axes[0].grid(alpha=0.3)

# Residualized
axes[1].hist(cohort_with_phers['r_phers'].to_list(),
             bins=50, alpha=0.6, label='Residualized PheRS', color='orange', edgecolor='black')
axes[1].set_xlabel('Studentized Residuals')
axes[1].set_ylabel('Number of Participants')
axes[1].set_title('Residualized PheRS Distribution')
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.show()


# ## Section 7: Export Results and Summary

# In[ ]:


print("\n=== FINAL SUMMARY ===\n")

# Select output columns
output_cols = ['person_id', 'time_zero', 'cox_phewas_37_z_stat_phers', 'r_phers']
output_df = cohort_with_phers.select(output_cols)

print("Output Columns:")
for col in output_cols:
    print(f"  - {col}")

print(f"\nFinal Dataset:")
print(f"  N participants: {len(output_df):,}")
display(output_df.head(10))


# In[ ]:


# Save results
output_path = input("\nEnter path for output file (e.g., cohort_with_phers.csv): ").strip()

if output_path.endswith('.tsv'):
    output_df.write_csv(output_path, separator='\t')
    print(f"✓ Saved to {output_path} (TSV format)")
else:
    output_df.write_csv(output_path, separator=',')
    print(f"✓ Saved to {output_path} (CSV format)")


# In[ ]:


print("""
================================================================================
ANALYSIS COMPLETE
================================================================================

Next Steps:
1. Download the output file from your Verily workspace
2. Use cox_phewas_37_z_stat_phers for raw burden scores
3. Use r_phers for demographic-adjusted analysis

Interpretation Tips:
- Zero scores: Participant has no qualifying post-index diagnoses
- Negative rPheRS: Fewer diagnoses than expected for demographic group
- Positive rPheRS: More diagnoses than expected for demographic group
- Magnitude: Larger absolute values = stronger departure from prediction

Statistical Considerations:
- Raw PheRS: Use for observational descriptive purposes
- Residualized PheRS: Preferred for association analyses (adjusts for confounding)
- When reporting: Include N, mean, SD, and range of both raw and residualized

Questions or Issues?
- Review this notebook's documentation sections
- Check validation checks (zero PheRS %, residualization quality)
- Verify your cohort file had required columns (person_id, time_zero)

Citations:
This implementation is based on methodology from analysis of post-COVID outcomes.
Please cite the original publication when using these scores.

================================================================================
""")

