#!/usr/bin/env python
# coding: utf-8

# # Packages and Setup

# In[ ]:


from IPython.display import display, HTML
from google.cloud import bigquery
import pandas as pd
import pyarrow as pa
import polars as pl
import os
import subprocess
import numpy as np, scipy as sps
# from scipy.stats import linregress
# import scipy.stats as stats
import matplotlib, matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime, date
from typing import Dict, List, Optional, Union
from jinja2 import Template, Environment


# In[ ]:


# This line allows for the plots to be displayed inline in the Jupyter notebook
get_ipython().run_line_magic('matplotlib', 'inline')

sns.set(style="whitegrid",font_scale=0.9)

palette = ['#0173b2', '#de8f05', '#8de5a1', '#d55e00', '#029e73', '#cc78bc', '#ece133', 
           '#56b4e9', '#ca9161', '#fbafe4', '#949494']

sns.set_palette(sns.color_palette(palette))
sns.color_palette(palette)


# In[ ]:


# show all columns in pandas
pd.set_option("display.max_columns", None)

# show full column width
pd.set_option('display.max_colwidth', 100)

# Polars string length to 100
pl.Config.set_fmt_str_lengths(100)

pl.Config.set_tbl_rows(50)


# # Helper Functions

# In[ ]:


def polars_gbq(query):
    """
    Execute a BigQuery SQL query and return results as a polars dataframe
    :param query: BigQuery SQL query
    :return: polars.DataFrame
    """
    client = bigquery.Client()
    query_job = client.query(query)
    rows = query_job.result()
    df = pl.from_arrow(rows.to_arrow())

    return df


# # HPV Cohort

# In[ ]:


# This query represents dataset "GWAS Cohort" for domain "person" and was generated for All of Us Controlled Tier Dataset v8
hpv_person_sql = """
    SELECT
        person.person_id 
    FROM
        `""" + os.environ["WORKSPACE_CDR"] + """.person` person   
    WHERE
        person.PERSON_ID IN (SELECT
            distinct person_id  
        FROM
            `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` cb_search_person  
        WHERE
            cb_search_person.person_id IN (SELECT
                person_id 
            FROM
                `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` p 
            WHERE
                has_whole_genome_variant = 1 
            UNION
            DISTINCT SELECT
                person_id 
            FROM
                `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_person` p 
            WHERE
                has_lr_whole_genome_variant = 1 ) 
            AND cb_search_person.person_id IN (SELECT
                criteria.person_id 
            FROM
                (SELECT
                    DISTINCT person_id, entry_date, concept_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_all_events` 
                WHERE
                    (concept_id IN (45572221) 
                    AND is_standard = 0 )) criteria 
            UNION
            DISTINCT SELECT
                criteria.person_id 
            FROM
                (SELECT
                    DISTINCT person_id, entry_date, concept_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_all_events` 
                WHERE
                    (concept_id IN(SELECT
                        DISTINCT c.concept_id 
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` c 
                    JOIN
                        (SELECT
                            CAST(cr.id as string) AS id       
                        FROM
                            `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` cr       
                        WHERE
                            concept_id IN (1567478)       
                            AND full_text LIKE '%_rank1]%'      ) a 
                            ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                            OR c.path LIKE CONCAT('%.', a.id) 
                            OR c.path LIKE CONCAT(a.id, '.%') 
                            OR c.path = a.id) 
                    WHERE
                        is_standard = 0 
                        AND is_selectable = 1) 
                    AND is_standard = 0 )) criteria 
            UNION
            DISTINCT SELECT
                criteria.person_id 
            FROM
                (SELECT
                    DISTINCT person_id, entry_date, concept_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_all_events` 
                WHERE
                    (concept_id IN (35206462) 
                    AND is_standard = 0 )) criteria 
            UNION
            DISTINCT SELECT
                criteria.person_id 
            FROM
                (SELECT
                    DISTINCT person_id, entry_date, concept_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_all_events` 
                WHERE
                    (concept_id IN(SELECT
                        DISTINCT c.concept_id 
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` c 
                    JOIN
                        (SELECT
                            CAST(cr.id as string) AS id       
                        FROM
                            `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` cr       
                        WHERE
                            concept_id IN (1572289)       
                            AND full_text LIKE '%_rank1]%'      ) a 
                            ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                            OR c.path LIKE CONCAT('%.', a.id) 
                            OR c.path LIKE CONCAT(a.id, '.%') 
                            OR c.path = a.id) 
                    WHERE
                        is_standard = 0 
                        AND is_selectable = 1) 
                    AND is_standard = 0 )) criteria 
            UNION
            DISTINCT SELECT
                criteria.person_id 
            FROM
                (SELECT
                    DISTINCT person_id, entry_date, concept_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_all_events` 
                WHERE
                    (concept_id IN (44826393, 44826394) 
                    AND is_standard = 0 )) criteria 
            UNION
            DISTINCT SELECT
                criteria.person_id 
            FROM
                (SELECT
                    DISTINCT person_id, entry_date, concept_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_all_events` 
                WHERE
                    (concept_id IN (44829052) 
                    AND is_standard = 0 )) criteria 
            UNION
            DISTINCT SELECT
                criteria.person_id 
            FROM
                (SELECT
                    DISTINCT person_id, entry_date, concept_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_all_events` 
                WHERE
                    (concept_id IN (44820669) 
                    AND is_standard = 0 )) criteria 
            UNION
            DISTINCT SELECT
                criteria.person_id 
            FROM
                (SELECT
                    DISTINCT person_id, entry_date, concept_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_all_events` 
                WHERE
                    (concept_id IN (44828175, 44835099, 44835101, 44835098, 44835100, 44831613) 
                    AND is_standard = 0 )) criteria 
            UNION
            DISTINCT SELECT
                criteria.person_id 
            FROM
                (SELECT
                    DISTINCT person_id, entry_date, concept_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_all_events` 
                WHERE
                    (concept_id IN (35205611) 
                    AND is_standard = 0 )) criteria 
            UNION
            DISTINCT SELECT
                criteria.person_id 
            FROM
                (SELECT
                    DISTINCT person_id, entry_date, concept_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_all_events` 
                WHERE
                    (concept_id IN (44828697) 
                    AND is_standard = 0 )) criteria ) 
            AND cb_search_person.person_id NOT IN (SELECT
                criteria.person_id 
            FROM
                (SELECT
                    DISTINCT person_id, entry_date, concept_id 
                FROM
                    `""" + os.environ["WORKSPACE_CDR"] + """.cb_search_all_events` 
                WHERE
                    (concept_id IN(SELECT
                        DISTINCT c.concept_id 
                    FROM
                        `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` c 
                    JOIN
                        (SELECT
                            CAST(cr.id as string) AS id       
                        FROM
                            `""" + os.environ["WORKSPACE_CDR"] + """.cb_criteria` cr       
                        WHERE
                            concept_id IN (44829737, 35205776, 44824945, 35225089)       
                            AND full_text LIKE '%_rank1]%'      ) a 
                            ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                            OR c.path LIKE CONCAT('%.', a.id) 
                            OR c.path LIKE CONCAT(a.id, '.%') 
                            OR c.path = a.id) 
                    WHERE
                        is_standard = 0 
                        AND is_selectable = 1) 
                    AND is_standard = 0 )) criteria ) )"""

hpv_person_df = pd.read_gbq(
    hpv_person_sql,
    dialect="standard",
    use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
    progress_bar_type="tqdm_notebook")


# In[ ]:


hpv_person_df['case']=1


# # Add Covariates to Full EHR Cohort

# In[ ]:


# Heart rate rhythm, Heart rate, Blood pressure panel, Adult Waist Circumference Protocol, PhenX - hip circumference protocol 020801
# Body height, Body weight, Diastolic blood pressure, Systolic blood pressure, Body mass index (BMI) [Ratio], Body weight Measured --pre pregnancy
vitals_exclusion = [3022318, 3027018, 3031203, 40759207, 40765148, 3036277, 3025315, 3012888, 3004249, 3038553, 3022281]
vitals_exclusion_str = ', '.join(map(str, vitals_exclusion))

# unique person_id w/ ICD or SNOMED in source_ values in observation/condition_occurrence
distinct_participants_cte = f"""
WITH distinct_participants AS (
    SELECT DISTINCT o.person_id
    FROM {version}.observation AS o
    JOIN {version}.concept AS c 
        ON o.observation_source_value = c.concept_code
        OR o.observation_source_concept_id = c.concept_id
    WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM', 'SNOMED')

    UNION DISTINCT
    
    SELECT DISTINCT co.person_id
    FROM {version}.condition_occurrence AS co
    JOIN {version}.concept AS c 
        ON co.condition_source_value = c.concept_code
        OR co.condition_source_concept_id = c.concept_id
    WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM', 'SNOMED')
)
"""

total_ehr_q = f"""
{distinct_participants_cte}
SELECT person_id
FROM distinct_participants
ORDER BY person_id
"""

total_ehr_df = polars_gbq(total_ehr_q)


# In[ ]:


total_ehr_df.height


# In[ ]:


# Here, end_of_study_age is age at death, or age at end of study period ()
demographics_q = f"""
{distinct_participants_cte}
SELECT DISTINCT
    p.person_id,
    CAST(p.birth_datetime AS DATE) AS dob,
    p_race_concept.concept_name as race,
    p_ethnicity_concept.concept_name as ethnicity,
    p_sex_at_birth_concept.concept_name as sex_at_birth,
    p_gender_concept.concept_name as gender,
    DATETIME_DIFF(
                IF(DATETIME(death_datetime) IS NULL, DATETIME('2023-10-01'), DATETIME(death_datetime)), 
                DATETIME(birth_datetime), 
                DAY
            )/365.2425 AS end_of_study_age
FROM
    {version}.person p
LEFT JOIN
    {version}.concept p_gender_concept 
        ON p.gender_concept_id = p_gender_concept.concept_id 
LEFT JOIN
    {version}.concept p_race_concept 
        ON p.race_concept_id = p_race_concept.concept_id 
LEFT JOIN
    {version}.concept p_ethnicity_concept 
        ON p.ethnicity_concept_id = p_ethnicity_concept.concept_id 
LEFT JOIN
    {version}.concept p_sex_at_birth_concept 
        ON p.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id  
LEFT JOIN
    {version}.death d
        ON p.person_id = d.person_id
WHERE
    p.person_id IN (SELECT person_id FROM distinct_participants)
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


ppi_employment_q = f"""
{distinct_participants_cte},
categorized AS (
    SELECT 
        obs.person_id,
        CASE
            WHEN obs.value_source_concept_id = 1585960 THEN 'disabled'
            WHEN obs.value_source_concept_id IN (1585953, 1585954) THEN 'employed'
            WHEN obs.value_source_concept_id = 1585959 THEN 'retired'
            WHEN obs.value_source_concept_id = 1585958 THEN 'student'
            WHEN obs.value_source_concept_id = 1585957 THEN 'homemaker'
            WHEN obs.value_source_concept_id IN (1585955, 1585956) THEN 'out of work'
            WHEN obs.value_source_concept_id IN (903079, 903096) THEN 'unk_skip'
        END AS employ,
        CASE
            WHEN obs.value_source_concept_id = 1585960 THEN 1
            WHEN obs.value_source_concept_id IN (1585953, 1585954) THEN 2
            WHEN obs.value_source_concept_id = 1585959 THEN 3
            WHEN obs.value_source_concept_id = 1585958 THEN 4
            WHEN obs.value_source_concept_id = 1585957 THEN 5
            WHEN obs.value_source_concept_id IN (1585955, 1585956) THEN 6
            WHEN obs.value_source_concept_id IN (903079, 903096) THEN 7
        END AS priority
    FROM 
        {version}.observation AS obs
    WHERE 
        obs.observation_source_concept_id IN (1585952)
        AND obs.person_id IN (SELECT person_id FROM distinct_participants)
),
ranked AS (
    SELECT
        person_id,
        employ,
        ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY priority) AS rank
    FROM 
        categorized
)
SELECT 
    person_id,
    employ
FROM 
    ranked
WHERE 
    rank = 1
"""

ppi_employment_df = polars_gbq(ppi_employment_q)


# In[ ]:


demographics_df = demographics_df.join(
    ppi_employment_df, 
    on='person_id',
    how='left' 
)


# In[ ]:


ppi_income_q = f"""
{distinct_participants_cte},
categorized AS (
    SELECT 
        obs.person_id,
        CASE
            WHEN obs.value_source_concept_id = 1585960 THEN 'disabled'
            WHEN obs.value_source_concept_id = 1585376 THEN 'less_10'
            WHEN obs.value_source_concept_id = 1585377 THEN '10_25'
            WHEN obs.value_source_concept_id = 1585378 THEN '25_35'
            WHEN obs.value_source_concept_id = 1585379 THEN '35_50'
            WHEN obs.value_source_concept_id = 1585380 THEN '50_75'
            WHEN obs.value_source_concept_id = 1585381 THEN '75_100'
            WHEN obs.value_source_concept_id = 1585382 THEN '100_150'
            WHEN obs.value_source_concept_id = 1585383 THEN '150_200'
            WHEN obs.value_source_concept_id = 1585384 THEN '200_more'
            WHEN obs.value_source_concept_id IN (903079, 903096) THEN 'skip'
        END AS income,
    FROM 
        {version}.observation AS obs
    WHERE 
        obs.observation_source_concept_id IN (1585375)
        AND obs.person_id IN (SELECT person_id FROM distinct_participants)
)
SELECT 
    person_id,
    income
FROM 
    categorized
"""

ppi_income_df = polars_gbq(ppi_income_q)


# In[ ]:


demographics_df = demographics_df.join(
    ppi_income_df, 
    on='person_id',
    how='left' 
)


# In[ ]:


ppi_edu_q = f"""
{distinct_participants_cte},
categorized AS (
    SELECT 
        obs.person_id,
        CASE
        WHEN obs.value_source_concept_id = 1585941 THEN 'none'
        WHEN obs.value_source_concept_id = 1585942 THEN '1_4'
        WHEN obs.value_source_concept_id = 1585943 THEN '5_8'
        WHEN obs.value_source_concept_id = 1585944 THEN '9_11'
        WHEN obs.value_source_concept_id = 1585945 THEN '12_GED'
        WHEN obs.value_source_concept_id = 1585946 THEN 'coll_1_3'
        WHEN obs.value_source_concept_id = 1585947 THEN 'coll_grad'
        WHEN obs.value_source_concept_id = 1585948 THEN 'adv_degree'
        WHEN obs.value_source_concept_id IN (903079, 903096) THEN 'skip'
        END AS edu,
    FROM 
        {version}.observation AS obs
    WHERE 
        obs.observation_source_concept_id IN (1585940)
        AND obs.person_id IN (SELECT person_id FROM distinct_participants)
)
SELECT 
    person_id,
    edu
FROM 
    categorized
"""

ppi_edu_df = polars_gbq(ppi_edu_q)


# In[ ]:


demographics_df = demographics_df.join(
    ppi_edu_df, 
    on='person_id',
    how='left' 
)


# In[ ]:


ppi_ins_q = f"""
{distinct_participants_cte},
categorized AS (
    SELECT 
        obs.person_id,
        CASE
        WHEN obs.value_source_concept_id = 1585387 THEN '1'
        WHEN obs.value_source_concept_id = 1585388 THEN '0'
        WHEN obs.value_source_concept_id IN (903096, 903087, 903079) THEN '-9'
        END AS ins,
    FROM 
        {version}.observation AS obs
    WHERE 
        obs.observation_source_concept_id IN (1585386)
        AND obs.person_id IN (SELECT person_id FROM distinct_participants)
)
SELECT 
    person_id,
    ins
FROM 
    categorized
"""

ppi_ins_df = polars_gbq(ppi_ins_q)


# In[ ]:


ppi_ins_type_q = f"""
{distinct_participants_cte},
categorized AS (
    SELECT 
        obs.person_id,
        CASE
            WHEN obs.value_source_concept_id = 43529209 THEN 'medicaid'
            WHEN obs.value_source_concept_id = 43529926 THEN 'va'
            WHEN obs.value_source_concept_id = 43529920 THEN 'military'
            WHEN obs.value_source_concept_id = 43529111 THEN 'ihs'
            WHEN obs.value_source_concept_id = 43529210 THEN 'medicare'
            WHEN obs.value_source_concept_id = 43529120 THEN 'emp'
            WHEN obs.value_source_concept_id = 43529119 THEN 'purchased'
            WHEN obs.value_source_concept_id = 43528423 THEN 'other'
            WHEN obs.value_source_concept_id = 43529095 THEN 'none'
            WHEN obs.value_source_concept_id IN (903096, 46237613) THEN '-9'
        END AS ins_type,
        CASE
            WHEN obs.value_source_concept_id = 43529209 THEN 1
            WHEN obs.value_source_concept_id = 43529926 THEN 2
            WHEN obs.value_source_concept_id = 43529920 THEN 3
            WHEN obs.value_source_concept_id = 43529111 THEN 4
            WHEN obs.value_source_concept_id = 43529210 THEN 5
            WHEN obs.value_source_concept_id = 43529120 THEN 6
            WHEN obs.value_source_concept_id = 43529119 THEN 7
            WHEN obs.value_source_concept_id = 43528423 THEN 8
            WHEN obs.value_source_concept_id = 43529095 THEN 0
            WHEN obs.value_source_concept_id IN (903096, 46237613) THEN 7
        END AS priority
    FROM 
        {version}.observation AS obs
    WHERE 
        obs.observation_source_concept_id IN (43528428)
        AND obs.person_id IN (SELECT person_id FROM distinct_participants)
),
ranked AS (
    SELECT
        person_id,
        ins_type,
        ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY priority) AS rank
    FROM 
        categorized
)
SELECT 
    person_id,
    ins_type
FROM 
    ranked
WHERE 
    rank = 1
"""

ppi_ins_type_df = polars_gbq(ppi_ins_type_q)


# In[ ]:


ppi_ins_df = ppi_ins_df.join(
    ppi_ins_type_df, 
    on='person_id',
    how='full', coalesce=True 
)


# In[ ]:


ppi_ins_df = ppi_ins_df.with_columns(
    pl.when((pl.col('ins') == '-9') & (pl.col('ins_type').is_not_null()) & (pl.col('ins_type') == 'none'))
    .then(0)
    .when((pl.col('ins') == '-9') & (pl.col('ins_type').is_not_null()) & (pl.col('ins_type') != '-9'))
    .then(1)
    .otherwise(pl.col('ins'))
    .alias('ins')
)


# In[ ]:


demographics_df = demographics_df.join(
    ppi_ins_df, 
    on='person_id',
    how='left' 
)


# In[ ]:


# total: 358,180
demographics_df.null_count()


# In[ ]:


demographics_df = demographics_df.to_pandas()


# In[ ]:


demographics_df[['race', 'ethnicity']].value_counts()


# In[ ]:


demographics_df = pl.from_pandas(demographics_df)


# In[ ]:


demographics_df = demographics_df.join(
    pl.from_pandas(hpv_person_df),
    on='person_id',
    how='left'
)


# In[ ]:


# Fill null values in the 'case' column with 0
demographics_df = demographics_df.with_columns(
    pl.col('case').fill_null(0)
)


# In[ ]:


demographics_pd = demographics_df.to_pandas()


# In[ ]:


def update_covariates(df, sarcoid=True):
    age_bins = [17, 24, 34, 44, 54, 64, 74, float('inf')]
    age_labels = ['17_24', '25_34', '35_44', '45_54', '55_64', '65_74', '75_plus']
    dx_labels = ['dx_17_24', 'dx_25_34', 'dx_35_44', 'dx_45_54', 'dx_55_64', 'dx_65_74', 'dx_75_plus']

    # Create binary columns for each age group
    for i, label in enumerate(age_labels):
        df[label] = ((df['end_of_study_age'] >= age_bins[i]) & (df['end_of_study_age'] < age_bins[i + 1])).astype(int)
    if sarcoid:
        for i, label in enumerate(dx_labels):
            df[label] = ((df['age_at_diagnosis'] >= age_bins[i]) & (df['age_at_diagnosis'] < age_bins[i + 1])).astype(int)

    # Handling sex_at_birth with specified conditions
    conditions_sex = {'Female': 1, 'Male': 0, 'NoMatch': -9, 'Skip': -9, 'PNA': -9, 'Intersex': -9}
    df['female'] = df['sex_at_birth'].map(conditions_sex).fillna(-9)
    
    # Race
    df['black'] = (df['race'] == 'Black').astype(int)
    df['asian'] = (df['race'] == 'Asian').astype(int)
    df['nhpi'] = (df['race'] == 'NHPI').astype(int)
    df['multip'] = (df['race'] == 'Multip').astype(int)
    df['white'] = (df['race'] == 'White').astype(int)
    df['mena'] = (df['race'] == 'MENA').astype(int)
    df['aian'] = (df['race'] == 'AIAN').astype(int)
    df['race_none'] = (df['race'] == 'None').astype(int)
    df['race_NoInd'] = (df['race'] == 'NoInd').astype(int)
    df['race_skip'] = ((df['race'] == 'Skip') | (df['race'] == 'PNA')).astype(int)

    # Ethnicity
    df['hispanic'] = (df['ethnicity'] == 'HispLat').astype(int)
    df['not_hispanic'] = (df['ethnicity'] == 'NotHisp').astype(int)
    df['ethnicity_none'] = ((df['ethnicity'] == 'None') | (df['ethnicity'] == 'NoMatch')).astype(int)
    df['ethnicity_skip'] = ((df['ethnicity'] == 'Skip') | (df['ethnicity'] == 'PNA')).astype(int)

    # Employment
    df['employed'] = (df['employ'] == 'employed').astype(int)
    df['student'] = (df['employ'] == 'student').astype(int)
    df['homemaker'] = (df['employ'] == 'homemaker').astype(int)
    df['retired'] = (df['employ'] == 'retired').astype(int)
    df['out_of_work'] = (df['employ'] == 'out of work').astype(int)
    df['disabled'] = (df['employ'] == 'disabled').astype(int)
    df['emp_skip'] = (df['employ'] == 'unk_skip').astype(int)

    # Income
    df['less_25'] = ((df['income'] == 'less_10') | (df['income'] == '10_25')).astype(int)
    df['25_100'] = ((df['income'] == '25_35') | (df['income'] == '35_50') | (df['income'] == '50_75') | (df['income'] == '75_100')).astype(int)
    df['100_more'] = ((df['income'] == '100_150') | (df['income'] == '150_200') | (df['income'] == '200_more')).astype(int)
    df['inc_skip'] = (df['income'] == 'skip').astype(int)

    # Education
    df['edu_none'] = (df['edu'] == 'none').astype(int)
    df['k_12'] = ((df['edu'] == '1_4') | (df['edu'] == '5_8') | (df['edu'] == '9_11') | (df['edu'] == '12_GED')).astype(int)
    df['coll_1_3'] = (df['edu'] == 'coll_1_3').astype(int)
    df['coll_grad'] = (df['edu'] == 'coll_grad').astype(int)
    df['adv_degree'] = (df['edu'] == 'adv_degree').astype(int)
    df['edu_skip'] = (df['edu'] == 'skip').astype(int)

    # Insurance type
    df['medicaid'] = (df['ins_type'] == 'medicaid').astype(int)
    df['medicare'] = (df['ins_type'] == 'medicare').astype(int)
    df['emp'] = (df['ins_type'] == 'emp').astype(int)
    df['none'] = (df['ins_type'] == 'none').astype(int)
    df['other'] = ((df['ins_type'] == 'ihs') | (df['ins_type'] == 'purchased other') | (df['ins_type'] == 'military') | (df['ins_type'] == 'va')).astype(int)
    df['ins_skip'] = (df['ins_type'] == '-9').astype(int)

    return df

demographics_pd = update_covariates(demographics_pd, sarcoid=False)


# In[ ]:


def describe_group(df, sarcoid=True):
    """
    Generate descriptive statistics for a group within the DataFrame.
    """
    def custom_format(x):
        return f"{x:.1f}".rstrip('0').rstrip('.') if not pd.isna(x) else "NaN"

    # Descriptive statistics for continuous variables
    age_desc = df['end_of_study_age'].agg(['median', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)]).map(custom_format)
    
    if sarcoid == True:
        dx_age_desc = df['age_at_diagnosis'].agg(['median', lambda x: x.quantile(0.25), lambda x: x.quantile(0.75)]).map(custom_format)

    # Count and percentages for categorical variables
    def calc_percentage(col):
        valid_data = df[col][df[col] != -9]
        count = valid_data.sum()
        total = valid_data.count()
        perc = np.round(count / total * 100, 1) if total != 0 else 0
        
        # Handle small counts
        if 1 <= count < 20:
            count_str = "< 20"
            perc_str = f"< {np.round(20/total * 100, 1)}"
        else:
            count_str = str(count)
            perc_str = custom_format(perc)

        # Only include total if different from df length
        output = f"{count_str} ({perc_str})" if total == len(df) else f"{count_str}/{total} ({perc_str})"
        return output

    # Applying calc_percentage to various columns
    demographics_percentages = {col: calc_percentage(col) for col in ['female', 'white', 'race_skip', 'black', 
                                                                      'asian', 'multip', 'mena', 'aian', 'nhpi', 
                                                                      'race_none', 'race_NoInd', 'hispanic']}

    # Applying calc_percentage to various columns
    ses_percentages = {col: calc_percentage(col) for col in ['employed', 'retired', 'disabled', 'student', 'homemaker', 'out_of_work', 
                                                             'emp_skip', 'less_25', '25_100', '100_more', 'inc_skip', 'edu_none', 'k_12',
                                                             'coll_1_3', 'coll_grad', 'adv_degree', 'edu_skip', 'medicare', 
                                                             'medicaid', 'emp', 'none', 'other', 'ins_skip']}
        
    # Building the results dictionary
    if sarcoid == True:
        results = {
            'Total': len(df),
            'End of Study Age, median (IQR)': f"{age_desc['median']} ({age_desc.iloc[1]}, {age_desc.iloc[2]})",
            'Sarcoid Dx Age, median (IQR)': f"{dx_age_desc['median']} ({dx_age_desc.iloc[1]}, {dx_age_desc.iloc[2]})",
            **{f"{key.capitalize()} (n (%))": value for key, value in demographics_percentages.items()},
            **{f"{key.capitalize()} (n (%))": value for key, value in ses_percentages.items()}

        }    
    else:
        results = {
            'Total': len(df),
            'End of Study Age, median (IQR)': f"{age_desc['median']} ({age_desc.iloc[1]}, {age_desc.iloc[2]})",
            'Sarcoid Dx Age, median (IQR)': f"NA",
            **{f"{key.capitalize()} (n (%))": value for key, value in demographics_percentages.items()},
            **{f"{key.capitalize()} (n (%))": value for key, value in ses_percentages.items()}

        }           

    return pd.DataFrame(results, index=[0])

controls = describe_group(demographics_pd[demographics_pd['case']==1], sarcoid=False)
cases = describe_group(demographics_pd[demographics_pd['case']==0], sarcoid=False)

# Then concatenate with total_ehr_cohort using a single-level index
table_1 = pd.concat([
    pd.DataFrame(controls).assign(category='Controls'),
    pd.DataFrame(cases).assign(category='Cases')
]).set_index('category', append=True)


# In[ ]:


transposed_table = table_1.transpose()
transposed_table


# In[ ]:


demographics_pd.to_csv('hpv_gwas_cohort.csv', index=False)


# # PheWAS

# In[ ]:


from PheTK.PheWAS import PheWAS
from PheTK.Phecode import Phecode
from PheTK.Plot import Plot
from PheTK.Cohort import Cohort


# In[ ]:


phecode = Phecode(platform="aou")
phecode.count_phecode(
    phecode_version="X", 
    icd_version="US",
    phecode_map_file_path=None, 
    output_file_name="phecode_X_counts.csv"
)


# In[ ]:


phewas_cohort = demographics_pd[demographics_pd['female'].isin([1, 0])]


# In[ ]:


phewas_cohort.to_csv('hpv_phewas_cohort.csv', index=False)


# In[ ]:


hpv_phewas = PheWAS(
    phecode_version="X",
    phecode_count_csv_path="phecode_X_counts.csv",
    cohort_csv_path="hpv_phewas_cohort.csv",
    sex_at_birth_col="female",
    male_as_one=True,
    covariate_cols=["end_of_study_age", "female", "white", "black", "asian", "hispanic"],
    independent_variable_of_interest="case",
    min_cases=50,
    min_phecode_count=2,
    output_file_name="hpv_phewas_results.csv"
)
hpv_phewas.run()


# In[ ]:


my_bucket = os.getenv('WORKSPACE_BUCKET')

# copy csv file to the bucket
args = ["gsutil", "cp", f"./hpv_phewas_results.csv", f"{my_bucket}/data/phewas/"]
output = subprocess.run(args, capture_output=True)

# print output from gsutil
output.stderr


# In[ ]:


p = Plot("hpv_phewas_results.csv")
p.manhattan(label_values="p_value", label_count=50, save_plot=False)


# In[ ]:


pgrm_bucket = '{bucket}'
acaf_person_id_df = pl.read_csv(f'{pgrm_bucket}/data/hail/pgrm_filtered_acaf_plink.fam', 
                                separator="\t", 
                                new_columns=["A","B","C","D","E","F"])

acaf_person_id_df = acaf_person_id_df.select('B')


# In[ ]:


hpv_cases = pl.read_csv('hpv_gwas_cohort.csv', try_parse_dates=True)


# In[ ]:


# Count people who are both in hpv_cases with case=1 and in acaf_person_id_df
overlap_count = (
    hpv_cases
    .filter(pl.col("case") == 1)
    .join(
        acaf_person_id_df.select(pl.col("B").alias("person_id")),
        on="person_id",
        how="inner"
    )
    .select(pl.len())
    .item()
)

print(f"Number of people with case=1 that are also in acaf_person_id_df: {overlap_count}")

