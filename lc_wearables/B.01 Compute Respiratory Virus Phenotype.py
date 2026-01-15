#!/usr/bin/env python
# coding: utf-8

# # Respiratory Viruses Cohort
# This notebook describes the queries performed to capture a cohort of those with viral respiratory infections.
# 
# NOTE: The current age calculated in this notebook reflects the age as of today (the day the notebook is ran), NOT the age at enrollment, nor the age at diagnosis. In addition, due to the month and day of birth supression in the Controlled Tier data, the age was calculated using the participant year of birth and today's year.
# 
# This notebook differs from the published workflow (<a href="https://www.nature.com/articles/s41598-025-02183-9">Waxse et al., 2025</a>) in the following ways:
# - added SARS-CoV-2 tests to mirror N3C test list

# ## Import Packages and Codes

# In[ ]:


# !pip install polars && pip install "polars[timezone]"


# In[ ]:


from IPython.display import display, HTML
from google.cloud import bigquery
import pandas as pd
import pyarrow as pa
import polars as pl
import os
import subprocess
import numpy as np
import matplotlib, matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
from datetime import timedelta, datetime as dt
import pytz


# In[ ]:


# expected_versions = {
#     "pandas": "2.2.2",
#     "pyarrow": "9.0.0",
#     "polars": "0.20.26",
#     "numpy": "1.23.5",
#     "matplotlib": "3.7.2",
#     "seaborn": "0.13.2",
#     "pytz": "2023.3"
# }
expected_versions = {
    "pandas": "2.0.3",
    "pyarrow": "9.0.0",
    "polars": "1.19.0",
    "numpy": "1.23.5",
    "matplotlib": "3.7.2",
    "seaborn": "0.12.2",
    "pytz": "2023.3"
}

# Check versions
libraries = [pd, pa, pl, np, matplotlib, sns, pytz]
names = ["pandas", "pyarrow", "polars", "numpy", "matplotlib", "seaborn", "pytz"]

print("Checking library versions:")
for lib, name in zip(libraries, names):
    current_version = lib.__version__
    expected_version = expected_versions.get(name)
    
    if current_version != expected_version:
        print(f"{name}: expected {expected_version}, but found {current_version}")
    else:
        print(f"{name}: version matched")


# In[ ]:


# This line allows for the plots to be displayed inline in the Jupyter notebook
get_ipython().run_line_magic('matplotlib', 'inline')

sns.set(style="whitegrid",font_scale=1)


# In[ ]:


version = get_ipython().run_line_magic('env', 'WORKSPACE_CDR')
version


# In[ ]:


# show all columns in pandas
pd.set_option("display.max_columns", None)

# show full column width
pd.set_option('display.max_colwidth', 100)

# Polars string length to 100
pl.Config.set_fmt_str_lengths(100)

# Set the row limit to a higher value
pl.Config.set_tbl_rows(50)


# ## Helper Functions
# 

# In[ ]:


def polars_gbq(query):
    """Execute BigQuery SQL and return polars DataFrame"""
    client = bigquery.Client()
    return pl.from_arrow(client.query(query).result().to_arrow())

# In[ ]:


def construct_concept_query(version, search_terms, exclusion_patterns):
    """
    Construct a SQL query string for concept codes based on a search term 
    and a list of exclusion patterns.
    
    Args:
    version (str): The dataset version to query against.
    search_terms (str or list of str): One or more terms to search for within concept names.
    exclusion_patterns (list[str]): A list of patterns to exclude from the search results.
    
    Returns:
    str: A SQL query string.
    """
    # Ensure search_terms is a list even if a single term is provided
    if isinstance(search_terms, str):
        search_terms = [search_terms]
        
    # Convert the list of search terms into a part of the WHERE clause
    search_terms_conditions = " OR ".join(
        [f"LOWER(c.concept_name) LIKE '%{term.lower()}%'" for term in search_terms]
    )
    
    # Convert the list of exclusion patterns into a format suitable for the UNNEST(ARRAY[..]) in SQL
    exclusion_patterns_formatted = ', '.join(f"'{pattern}'" for pattern in exclusion_patterns)

    # Define the query
    query_string = f"""
    SELECT 
        DISTINCT c.concept_id,
        c.concept_code,
        c.vocabulary_id,
        c.concept_name
    FROM 
        `{version}.concept` c
    WHERE
        c.vocabulary_id LIKE '%ICD%'
        AND c.domain_id = 'Condition'
        AND ({search_terms_conditions})
        AND NOT EXISTS (
            SELECT 1
            FROM UNNEST(ARRAY[{exclusion_patterns_formatted}]) AS excluded_pattern
            WHERE LOWER(c.concept_name) LIKE excluded_pattern
        )
    ORDER BY c.concept_name DESC
    """

    return query_string


# In[ ]:


def fetch_co(version, query, return_count=False):
    """
    Fetch condition occurrence table based on a given query for 
    codes and optionally return the count of unique patients.
    
    Args:
    version (str): The dataset version to query against.
    query (str): The query that fetches the relevant concept IDs.
    return_count (bool, optional): If True, returns the count of 
    unique patients. Otherwise, returns detailed dataframe.
    
    Returns:
    int or polars.DataFrame: Depending on return_count, 
    returns either the count of unique patients or a dataframe 
    with detailed occurrences.
    """
    # Define the query for condition occurrences
    co_codes = f"""
        SELECT DISTINCT
            co.person_id,
            condition_source_concept_id as concept_id,
            c1.concept_name AS concept_name,
            c1.concept_code as code,
            c1.vocabulary_id as vocab,
            condition_start_datetime as datetime,
            CAST(co.visit_occurrence_id AS STRING) as visit_id,
            v.visit_start_datetime as visit_start,
            c2.concept_name AS visit_name,
        FROM
            {version}.condition_occurrence co
        INNER JOIN
            {version}.concept c1 ON co.condition_source_concept_id = c1.concept_id
        LEFT JOIN
            {version}.visit_occurrence v ON co.visit_occurrence_id = v.visit_occurrence_id AND co.person_id = v.person_id
        LEFT JOIN
            {version}.concept c2 ON v.visit_concept_id = c2.concept_id
        WHERE 
            co.condition_source_concept_id IN (SELECT concept_id FROM ({query}))
        ORDER BY person_id, condition_start_datetime
    """
    
    if return_count:
        # If only the count of unique patients is needed
        count_query = f"SELECT COUNT(*) AS patient_count FROM ({co_codes})"
        result = polars_gbq(count_query)
        return result.get_column('patient_count')[0]
    else:
        # If the detailed dataframe is needed
        detailed_result = polars_gbq(co_codes)
        return detailed_result


# In[ ]:


def fetch_o(version, query, return_count=False):
    """
    Fetch observation table based on a given query for 
    codes and optionally return the count of unique patients.
    
    Args:
    version (str): The dataset version to query against.
    query (str): The query that fetches the relevant concept IDs.
    return_count (bool, optional): If True, returns the count of 
    unique patients. Otherwise, returns detailed dataframe.
    
    Returns:
    int or polars.DataFrame: Depending on return_count, 
    returns either the count of unique patients or a dataframe 
    with detailed occurrences.
    """
    # Define the query for condition occurrences
    o_codes = f"""
        SELECT DISTINCT
            o.person_id,
            observation_source_concept_id as concept_id,
            c1.concept_name AS concept_name,
            c1.concept_code as code,
            c1.vocabulary_id as vocab,
            observation_datetime as datetime,
            CAST(o.visit_occurrence_id AS STRING) as visit_id,
            v.visit_start_datetime as visit_start,
            c2.concept_name AS visit_name,
        FROM
            {version}.observation o
        INNER JOIN
            {version}.concept c1 ON o.observation_source_concept_id = c1.concept_id
        LEFT JOIN
            {version}.visit_occurrence v ON o.visit_occurrence_id = v.visit_occurrence_id AND o.person_id = v.person_id
        LEFT JOIN
            {version}.concept c2 ON v.visit_concept_id = c2.concept_id
        WHERE 
            o.observation_source_concept_id IN (SELECT concept_id FROM ({query}))
        ORDER BY person_id, observation_datetime
    """
    
    if return_count:
        # If only the count of unique patients is needed
        count_query = f"SELECT COUNT(*) AS patient_count FROM ({o_codes})"
        result = polars_gbq(count_query)
        return result.get_column('patient_count')[0]
    else:
        # If the detailed dataframe is needed
        detailed_result = polars_gbq(o_codes)
        return detailed_result


# # Cohort Creation
# <b>Diagnosis</b><br>
# <p>With ICD9, and ICD10 codes present in <code>source_concept_id</code> from the <code>condition_occurrence</code> and <code>observation</code> tables, we will use these to identify diagnoses. Prior work compared <code>concept_id</code> and <code>source_concept_id</code>, with minimal loss and preserved datatype in the latter.<br><br>
# Exclusion arrays were manually curated to exclude non-respiratory or ambiguous codes. Ancestor analysis detected adenoviral infection as a missed ICD code on the first iteration.</p><br>
# <b>Test</b><br>
# <p>LOINC concepts pertaining to each virus were queried from the <code>measurement</code> table.<br>
# Initial attempts added individual tests from the Cohort Builder to relevant respiratory panels, but here we instead select all tests for each virus, and then broaden the search to also include other tests included in panels that use these tests.</p><br>
# <b>Treatment</b><br>
# <p>Antiviral entries for <code>domain_id</code> = 'Drug' for influenza and SARS-CoV-2 were retrieved from the <code>drug_exposure</code> table.<br>
# For influenza, oseltamivir (2-6 days of treatment to avoid prophylaxis), zanamivir (2-6 days of treatment to avoid prophylaxis), and baloxavir were used. For SARS-CoV-2, remdesivir, nirmatrelvir/ritonavir, and molnupiravir were used'</p><br>

# ## Total EHR Population 
# with ICD codes, LOINC measurement, or drug exposure

# In[ ]:


# Heart rate rhythm, Heart rate, Blood pressure panel, Adult Waist Circumference Protocol, PhenX - hip circumference protocol 020801
# Body height, Body weight, Diastolic blood pressure, Systolic blood pressure, Body mass index (BMI) [Ratio], Body weight Measured --pre pregnancy
vitals_exclusion = [3022318, 3027018, 3031203, 40759207, 40765148, 3036277, 3025315, 3012888, 3004249, 3038553, 3022281]
vitals_exclusion_str = ', '.join(map(str, vitals_exclusion))

# unique person_id w/ ICD in source_ values in observation/condition_occurrence, LOINC in measurement (vitals excluded), Drug in drug_exposure
total_ehr_q = f"""
WITH distinct_participants AS (
    SELECT DISTINCT o.person_id
    FROM {version}.observation AS o
    JOIN {version}.concept AS c 
        ON o.observation_source_value = c.concept_code
        OR o.observation_source_concept_id = c.concept_id
    WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM')

    UNION DISTINCT
    
    SELECT DISTINCT co.person_id
    FROM {version}.condition_occurrence AS co
    JOIN {version}.concept AS c 
        ON co.condition_source_value = c.concept_code
        OR co.condition_source_concept_id = c.concept_id
    WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM')
    
    UNION DISTINCT
    
    SELECT DISTINCT m.person_id
    FROM {version}.measurement AS m
    JOIN {version}.concept AS c ON m.measurement_concept_id = c.concept_id
    WHERE c.vocabulary_id = 'LOINC' AND c.concept_id NOT IN ({vitals_exclusion_str})
    
    UNION DISTINCT
    
    SELECT DISTINCT de.person_id
    FROM {version}.drug_exposure AS de
    JOIN {version}.concept AS c ON de.drug_concept_id = c.concept_id
    WHERE c.domain_id = 'Drug'
)
SELECT person_id
FROM distinct_participants
ORDER BY person_id
"""

total_ehr_df = polars_gbq(total_ehr_q)


# In[ ]:


print(f"Total EHR Population with ICD codes, LOINC measurement, or drug exposure: {len(total_ehr_df)}")


# In[ ]:


# earliest and latest dates associated with ICD in source_ values in observation/condition_occurrence, LOINC in measurement (vitals excluded), Drug in drug_exposure
date_ehr_q = f"""
WITH date_ranges AS (
    SELECT 
        observation_date AS event_date
    FROM {version}.observation AS o
    JOIN {version}.concept AS c 
        ON o.observation_source_value = c.concept_code 
        OR o.observation_source_concept_id = c.concept_id
    WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM')
    
    UNION ALL
        
    SELECT 
        condition_start_date AS event_date
    FROM {version}.condition_occurrence AS co
    JOIN {version}.concept AS c 
        ON co.condition_source_value = c.concept_code
        OR co.condition_source_concept_id = c.concept_id  
    WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM')
    
    UNION ALL
    
    SELECT 
        measurement_date AS event_date
    FROM {version}.measurement AS m
    JOIN {version}.concept AS c ON m.measurement_concept_id = c.concept_id
    WHERE c.vocabulary_id = 'LOINC' AND c.concept_id NOT IN ({vitals_exclusion_str})
    
    UNION ALL
    
    SELECT 
        drug_exposure_start_date AS event_date
    FROM {version}.drug_exposure AS de
    JOIN {version}.concept AS c ON de.drug_concept_id = c.concept_id
    WHERE c.domain_id = 'Drug'
)
SELECT
    MIN(event_date) AS start_date,
    MAX(event_date) AS end_date
FROM
    date_ranges 
"""
polars_gbq(date_ehr_q)


# In[ ]:


# There are surprisingly early dates, but they are few in number.
# This modifies the GBQ to only include event_date if there are more than 20 in a given year:

# earliest and latest dates associated with ICD in source_ values in observation/condition_occurrence, LOINC in measurement (vitals excluded), Drug in drug_exposure
date_ehr_q = f"""
WITH date_ranges AS (
    SELECT 
        observation_date AS event_date
    FROM {version}.observation AS o
    JOIN {version}.concept AS c 
        ON o.observation_source_value = c.concept_code 
        OR o.observation_source_concept_id = c.concept_id
    WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM')
    
    UNION ALL
        
    SELECT 
        condition_start_date AS event_date
    FROM {version}.condition_occurrence AS co
    JOIN {version}.concept AS c 
        ON co.condition_source_value = c.concept_code
        OR co.condition_source_concept_id = c.concept_id  
    WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM')
    
    UNION ALL
    
    SELECT 
        measurement_date AS event_date
    FROM {version}.measurement AS m
    JOIN {version}.concept AS c ON m.measurement_concept_id = c.concept_id
    WHERE c.vocabulary_id = 'LOINC' AND c.concept_id NOT IN ({vitals_exclusion_str})
    
    UNION ALL
    
    SELECT 
        drug_exposure_start_date AS event_date
    FROM {version}.drug_exposure AS de
    JOIN {version}.concept AS c ON de.drug_concept_id = c.concept_id
    WHERE c.domain_id = 'Drug'
),
yearly_counts AS (
    SELECT 
        EXTRACT(YEAR FROM event_date) AS year,
        COUNT(*) AS event_count
    FROM date_ranges
    GROUP BY EXTRACT(YEAR FROM event_date)
    HAVING COUNT(*) > 20
)
SELECT
    MIN(dr.event_date) AS start_date,
    MAX(dr.event_date) AS end_date
FROM date_ranges dr
JOIN yearly_counts yc
    ON EXTRACT(YEAR FROM dr.event_date) = yc.year
"""
polars_gbq(date_ehr_q)


# ## Parainfluenza

# ### Parainfluenza Diagnoses

# In[ ]:


para_exclusion = [
    '%vaccin%', 
    '%non-respiratory%',
    '%otitis%',
    '%myocarditis%',
    '%bovine%',
    '%adenovirus%',
    '%myelitis%',
    '%gastrointestinal%',
    '%gastroenteritis%',
    '%central nervous system%',
    '%cns disorder%',
    '%encephalopathy%',
    '%encephalitis%',
    '%virus inconclusive%',
    '%not isolated%',
    '%not detected%',
    '%equine%'
]

# Querying the concept codes
para_codes_q = construct_concept_query(version, 'parainfluenza', para_exclusion)


# In[ ]:


para_codes = polars_gbq(para_codes_q)


# In[ ]:


para_codes


# In[ ]:


# fetch_co(version, para_codes_q, return_count=True)
# 263 v7
# 331 v8


# In[ ]:


# fetch_o(version, para_codes_q, return_count=True)
# 0 v7, v8


# In[ ]:


para_code_df = fetch_co(version, para_codes_q, return_count=False)


# In[ ]:


para_code_df['person_id'].n_unique()
# v7 136
# v8 197


# In[ ]:


para_code_counts = para_code_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)


# In[ ]:


para_codes = para_codes.join(
    para_code_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


para_codes.head()


# ### Parainfluenza Tests

# In[ ]:


# Search panels for parainfluenza tests
para_tests_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    c.vocabulary_id = 'LOINC'
    AND LOWER(c.concept_name) LIKE '%parainfluenza%'
    AND (
        LOWER(c.concept_name) LIKE '%[presence]%'
        OR LOWER(c.concept_name) LIKE '%[identifier]%'
        OR LOWER(c.concept_name) LIKE '%[cycle threshold]%'
        OR LOWER(c.concept_name) LIKE '%[susceptibility]%'
        OR LOWER(c.concept_name) LIKE '%viral load%'
        OR LOWER(c.concept_name) LIKE '%culture%'
        OR c.concept_name LIKE '%NAA%'
        OR c.concept_name LIKE '%Measurement of%'
    )
    AND NOT EXISTS (
        SELECT 1
        FROM UNNEST(ARRAY[
            '%porcine%', 
            '%canine%',
            '%bovine%',
            '%haemophilus%',
            '%cornea%',
            '%igg%',
            '%igm%',
            '%cerebral%',
            '% ab %',
            '%coronavirus%'
        ]) AS excluded_pattern
        WHERE LOWER(c.concept_name) LIKE excluded_pattern
    )
ORDER BY c.concept_name
"""

para_tests = polars_gbq(para_tests_q)


# In[ ]:


m_para_tests = f"""
        SELECT
            m.person_id,
            m.measurement_concept_id as concept_id,
            c.concept_name as concept_name,
            c.concept_code as code,
            c.vocabulary_id as vocab,
            m.measurement_datetime as datetime,
            c2.concept_name AS result,
            CAST(m.visit_occurrence_id AS STRING) as visit_id,
            v.visit_start_datetime as visit_start,
            c1.concept_name AS visit_name
        FROM
            {version}.measurement m
        INNER JOIN
            {version}.concept c ON m.measurement_concept_id = c.concept_id
        LEFT JOIN
            {version}.concept c2 ON m.value_as_concept_id = c2.concept_id
        LEFT JOIN
            {version}.visit_occurrence v ON m.visit_occurrence_id = v.visit_occurrence_id AND m.person_id = v.person_id
        LEFT JOIN
            {version}.concept c1 ON v.visit_concept_id = c1.concept_id      
        WHERE
            m.measurement_concept_id IN (SELECT DISTINCT concept_id FROM ({para_tests_q}))
        ORDER BY 
            person_id, measurement_datetime
        """

# polars_gbq(f"""SELECT COUNT(*) FROM ({m_para_tests})""")
# 58930 v7
# 118773 v8


# In[ ]:


para_test_df = polars_gbq(m_para_tests)


# In[ ]:


para_test_df['person_id'].n_unique()
# 7161 v7
# 14681 v8 


# In[ ]:


para_test_df.group_by(['concept_id', 'concept_name', 'result']).agg(pl.len().alias('count')).sort(['count'], descending=True).head()


# In[ ]:


para_test_counts = para_test_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

para_tests = para_tests.join(
    para_test_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


para_tests.sort('concept_code').head(10)


# ## Human Metapnuemovirus

# ### Human Metapnuemovirus Diagnoses

# In[ ]:


hmpv_exclusion = [
    '%vaccin%', 
    '%non-respiratory%',
    '%otitis%',
    '%myocarditis%',
    '%adenovirus%',
    '%bovine%',
    '%myelitis%',
    '%gastrointestinal%',
    '%gastroenteritis%',
    '%central nervous system%',
    '%cns disorder%',
    '%encephalopathy%',
    '%encephalitis%',
    '%virus inconclusive%',
    '%not isolated%',
    '%not detected%',
    '%equine%'
]

# Querying the concept codes
hmpv_codes_q = construct_concept_query(version, 'metapneumovirus', hmpv_exclusion)


# In[ ]:


hmpv_codes = polars_gbq(hmpv_codes_q)


# In[ ]:


hmpv_codes


# In[ ]:


# fetch_co(version, hmpv_codes_q, return_count=True)
# 608 v7
# 910 v8


# In[ ]:


# fetch_o(version, hmpv_codes_q, return_count=True)
# 0 v7, v8


# In[ ]:


hmpv_code_df = fetch_co(version, hmpv_codes_q, return_count=False)


# In[ ]:


hmpv_code_df['person_id'].n_unique()
# 305 v7
# 454 v8


# In[ ]:


hmpv_code_df.head()


# In[ ]:


hmpv_code_counts = hmpv_code_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

hmpv_codes = hmpv_codes.join(
    hmpv_code_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


hmpv_codes


# ### Human Metapneumovirus Tests

# In[ ]:


# Search panels for hMPV tests
hmpv_tests_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    c.vocabulary_id = 'LOINC'
    AND LOWER(c.concept_name) LIKE '%metapneumovirus%'
    AND (
        LOWER(c.concept_name) LIKE '%[presence]%'
        OR LOWER(c.concept_name) LIKE '%[identifier]%'
        OR LOWER(c.concept_name) LIKE '%[cycle threshold]%'
        OR LOWER(c.concept_name) LIKE '%[susceptibility]%'
        OR LOWER(c.concept_name) LIKE '%viral load%'
        OR LOWER(c.concept_name) LIKE '%culture%'
        OR c.concept_name LIKE '%NAA%'
        OR c.concept_name LIKE '%Measurement of%'
    )
    AND NOT EXISTS (
        SELECT 1
        FROM UNNEST(ARRAY[
            '%porcine%', 
            '%canine%',
            '%bovine%',
            '%haemophilus%',
            '%cornea%',
            '%igg%',
            '%igm%',
            '%cerebral%',
            '% ab %',
            '%adenovirus%'
        ]) AS excluded_pattern
        WHERE LOWER(c.concept_name) LIKE excluded_pattern
    )
ORDER BY c.concept_name
"""

hmpv_tests = polars_gbq(hmpv_tests_q)


# In[ ]:


m_hmpv_tests = f"""
        SELECT
            m.person_id,
            m.measurement_concept_id as concept_id,
            c.concept_name as concept_name,
            c.concept_code as code,
            c.vocabulary_id as vocab,
            m.measurement_datetime as datetime,
            c2.concept_name AS result,
            CAST(m.visit_occurrence_id AS STRING) as visit_id,
            v.visit_start_datetime as visit_start,
            c1.concept_name AS visit_name
        FROM
            {version}.measurement m
        INNER JOIN
            {version}.concept c ON m.measurement_concept_id = c.concept_id
        LEFT JOIN
            {version}.concept c2 ON m.value_as_concept_id = c2.concept_id
        LEFT JOIN
            {version}.visit_occurrence v ON m.visit_occurrence_id = v.visit_occurrence_id AND m.person_id = v.person_id
        LEFT JOIN
            {version}.concept c1 ON v.visit_concept_id = c1.concept_id      
        WHERE
            m.measurement_concept_id IN (SELECT DISTINCT concept_id FROM ({hmpv_tests_q}))
        ORDER BY 
            person_id, measurement_datetime
        """

# polars_gbq(f"""SELECT COUNT(*) FROM ({m_hmpv_tests})""")
# 15715 v7
# 30894 v8


# In[ ]:


hmpv_test_df = polars_gbq(m_hmpv_tests)


# In[ ]:


hmpv_test_df['person_id'].n_unique()
# 6934 v7
# 14088 v8


# In[ ]:


hmpv_test_df.group_by(['concept_id', 'concept_name', 'result']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(10)


# In[ ]:


hmpv_test_df.head()


# In[ ]:


hmpv_test_counts = hmpv_test_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

hmpv_tests = hmpv_tests.join(
    hmpv_test_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


hmpv_tests.head()


# ## Rhinovirus

# ### Rhinovirus Diagnoses

# In[ ]:


rv_exclusion = [
    '%vaccin%', 
    '%non-respiratory%',
    '%otitis%',
    '%myocarditis%',
    '%adenovirus%',
    '%bovine%',
    '%myelitis%',
    '%gastrointestinal%',
    '%gastroenteritis%',
    '%central nervous system%',
    '%cns disorder%',
    '%encephalopathy%',
    '%encephalitis%',
    '%virus inconclusive%',
    '%not isolated%',
    '%not detected%',
    '%equine%'
]

# Querying the concept codes
rv_codes_q = construct_concept_query(version, 'rhinovirus', rv_exclusion)


# In[ ]:


rv_codes = polars_gbq(rv_codes_q)


# In[ ]:


rv_codes


# In[ ]:


# fetch_co(version, rv_codes_q, return_count=True)
# 783 v7
# 1075 v8


# In[ ]:


# fetch_o(version, rv_codes_q, return_count=True)
# 0 v7, v8


# In[ ]:


rv_code_df = fetch_co(version, rv_codes_q, return_count=False)


# In[ ]:


rv_code_df['person_id'].n_unique()
# 365 v7
# 531 v8


# In[ ]:


rv_code_df.head()


# In[ ]:


rv_code_counts = rv_code_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

rv_codes = rv_codes.join(
    rv_code_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


rv_codes


# ### Rhinovirus Tests

# In[ ]:


# Search panels for rhinovirus tests
rv_tests_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    c.vocabulary_id = 'LOINC'
    AND LOWER(c.concept_name) LIKE '%rhinovirus%'
    AND (
        LOWER(c.concept_name) LIKE '%[presence]%'
        OR LOWER(c.concept_name) LIKE '%[identifier]%'
        OR LOWER(c.concept_name) LIKE '%[cycle threshold]%'
        OR LOWER(c.concept_name) LIKE '%[susceptibility]%'
        OR LOWER(c.concept_name) LIKE '%viral load%'
        OR LOWER(c.concept_name) LIKE '%culture%'
        OR c.concept_name LIKE '%NAA%'
        OR c.concept_name LIKE '%Measurement of%'
    )
    AND NOT EXISTS (
        SELECT 1
        FROM UNNEST(ARRAY[
            '%porcine%', 
            '%canine%',
            '%bovine%',
            '%haemophilus%',
            '%cornea%',
            '%igg%',
            '%igm%',
            '%cerebral%',
            '% ab %',
            '%adenovirus%'
        ]) AS excluded_pattern
        WHERE LOWER(c.concept_name) LIKE excluded_pattern
    )
ORDER BY c.concept_name
"""

rv_tests = polars_gbq(rv_tests_q)


# In[ ]:


m_rv_tests = f"""
        SELECT
            m.person_id,
            m.measurement_concept_id as concept_id,
            c.concept_name as concept_name,
            c.concept_code as code,
            c.vocabulary_id as vocab,
            m.measurement_datetime as datetime,
            c2.concept_name AS result,
            CAST(m.visit_occurrence_id AS STRING) as visit_id,
            v.visit_start_datetime as visit_start,
            c1.concept_name AS visit_name
        FROM
            {version}.measurement m
        INNER JOIN
            {version}.concept c ON m.measurement_concept_id = c.concept_id
        LEFT JOIN
            {version}.concept c2 ON m.value_as_concept_id = c2.concept_id
        LEFT JOIN
            {version}.visit_occurrence v ON m.visit_occurrence_id = v.visit_occurrence_id AND m.person_id = v.person_id
        LEFT JOIN
            {version}.concept c1 ON v.visit_concept_id = c1.concept_id      
        WHERE
            m.measurement_concept_id IN (SELECT DISTINCT concept_id FROM ({rv_tests_q}))
        ORDER BY 
            person_id, measurement_datetime
        """

rv_test_df = polars_gbq(m_rv_tests)
len(rv_test_df)
# 14873 v7
# 29609 v8


# In[ ]:


rv_test_df['person_id'].n_unique()
# 6418 v7
# 13506 v8


# In[ ]:


rv_test_df.group_by(['concept_id', 'concept_name', 'result']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(10)


# In[ ]:


rv_test_df.head()


# In[ ]:


rv_test_counts = rv_test_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

rv_tests = rv_tests.join(
    rv_test_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


rv_tests.head(10)


# ## Bocavirus
# Not used - no ICD/SNOMED codes and only nine negative tests.

# ### Bocavirus Diagnoses

# In[ ]:


bocavirus_exclusion = [
    '%vaccin%', 
    '%non-respiratory%',
    '%otitis%',
    '%myocarditis%',
    '%bovine%',
    '%central nervous system%',
    '%cns disorder%',
    '%encephalopathy%',
    '%encephalitis%',
    '%virus inconclusive%',
    '%not isolated%',
    '%not detected%',
    '%equine%'
]

# Querying the concept codes
bocavirus_codes_q = construct_concept_query(version, 'bocavirus', bocavirus_exclusion)


# In[ ]:


bocavirus_codes = polars_gbq(bocavirus_codes_q)


# In[ ]:


bocavirus_codes


# ### Bocavirus Tests

# In[ ]:


# Search panels for rhinovirus tests
bv_tests_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    c.vocabulary_id = 'LOINC'
    AND LOWER(c.concept_name) LIKE '%bocavirus%'
    AND (
        LOWER(c.concept_name) LIKE '%[presence]%'
        OR LOWER(c.concept_name) LIKE '%[identifier]%'
        OR LOWER(c.concept_name) LIKE '%[cycle threshold]%'
        OR LOWER(c.concept_name) LIKE '%[susceptibility]%'
        OR LOWER(c.concept_name) LIKE '%viral load%'
        OR LOWER(c.concept_name) LIKE '%culture%'
        OR c.concept_name LIKE '%NAA%'
        OR c.concept_name LIKE '%Measurement of%'
    )
    AND NOT EXISTS (
        SELECT 1
        FROM UNNEST(ARRAY[
            '%porcine%', 
            '%canine%',
            '%bovine%',
            '%cornea%',
            '%igg%',
            '%igm%',
            '%adenovirus%'
        ]) AS excluded_pattern
        WHERE LOWER(c.concept_name) LIKE excluded_pattern
    )
ORDER BY c.concept_name
"""

bv_tests = polars_gbq(bv_tests_q)


# In[ ]:


m_bv_tests = f"""
        SELECT
            m.person_id,
            m.measurement_concept_id as concept_id,
            c.concept_name as concept_name,
            c.concept_code as code,
            c.vocabulary_id as vocab,
            m.measurement_datetime as datetime,
            c2.concept_name AS result,
            CAST(m.visit_occurrence_id AS STRING) as visit_id,
            v.visit_start_datetime as visit_start,
            c1.concept_name AS visit_name
        FROM
            {version}.measurement m
        INNER JOIN
            {version}.concept c ON m.measurement_concept_id = c.concept_id
        LEFT JOIN
            {version}.concept c2 ON m.value_as_concept_id = c2.concept_id
        LEFT JOIN
            {version}.visit_occurrence v ON m.visit_occurrence_id = v.visit_occurrence_id AND m.person_id = v.person_id
        LEFT JOIN
            {version}.concept c1 ON v.visit_concept_id = c1.concept_id      
        WHERE
            m.measurement_concept_id IN (SELECT DISTINCT concept_id FROM ({bv_tests_q}))
        ORDER BY 
            person_id, measurement_datetime
        """

bv_test_df = polars_gbq(m_bv_tests)
len(bv_test_df)
# 12 v8


# In[ ]:


bv_test_df['person_id'].n_unique()
# 12 v8


# In[ ]:


bv_test_df.group_by(['concept_id', 'concept_name', 'result']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(10)


# In[ ]:


bv_test_counts = bv_test_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

bv_tests = bv_tests.join(
    bv_test_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


bv_tests.head(10)


# ## Adenovirus

# ### Adenovirus Diagnoses

# In[ ]:


adv_inclusion = [
    'adenovirus',
    'adenoviral pneumonia'
]

adv_exclusion = [
    '%vaccin%', 
    '%non-respiratory%',
    '%otitis%',
    '%rhinovirus%',
    '%myocarditis%',
    '%porcine%',
    '%caprine%',
    '%ovine%',
    '%simian%',
    '%murine%',
    '%bovine%',
    '%myelitis%',
    '%food%',
    '%gastrointestinal%',
    '%conjunctivitis%',
    '%enteritis%',
    '%central nervous system%',
    '%cns disorder%',
    '%encephalopathy%',
    '%keratitis%',
    '%hepatitis%',
    '%meningitis%',
    '%encephalitis%',
    '%virus inconclusive%',
    '%not isolated%',
    '%not detected%',
    '%equine%'
]

adv_codes_q = construct_concept_query(version, adv_inclusion, adv_exclusion)


# In[ ]:


adv_codes = polars_gbq(adv_codes_q)


# In[ ]:


adv_codes.slice(0,50)


# In[ ]:


# fetch_co(version, adv_codes_q, return_count=True)
# 401 v7
# 606 v8


# In[ ]:


# fetch_o(version, adv_codes_q, return_count=True)
# 0 v7, v8


# In[ ]:


adv_code_df = fetch_co(version, adv_codes_q, return_count=False)


# In[ ]:


adv_code_df['person_id'].n_unique()
# 162 v7
# 221 v8


# In[ ]:


adv_code_df.head()


# In[ ]:


adv_code_counts = adv_code_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

adv_codes = adv_codes.join(
    adv_code_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


adv_codes


# ### Adenovirus Tests

# In[ ]:


# Search panels for adenovirus tests
adv_tests_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    c.vocabulary_id = 'LOINC'
    AND LOWER(c.concept_name) LIKE '%adenovirus%'
    AND (
        LOWER(c.concept_name) LIKE '%[presence]%'
        OR LOWER(c.concept_name) LIKE '%[identifier]%'
        OR LOWER(c.concept_name) LIKE '%[cycle threshold]%'
        OR LOWER(c.concept_name) LIKE '%[susceptibility]%'
        OR LOWER(c.concept_name) LIKE '%viral load%'
        OR LOWER(c.concept_name) LIKE '%culture%'
        OR c.concept_name LIKE '%NAA%'
        OR c.concept_name LIKE '%Measurement of%'
    )
    AND NOT EXISTS (
        SELECT 1
        FROM UNNEST(ARRAY[
            '%porcine%', 
            '%canine%',
            '%bovine%',
            '%equine%',
            '%haemophilus%',
            '%cornea%',
            '%conjunctival%',
            '%stool%',
            '%urine%',
            '%blood%',
            '%serum%',
            '%plasma%',
            '%genital%',
            '%marrow%',
            '%skin%',
            '%tissue%',
            '%body%',
            '%igg%',
            '%igm%',
            '%cerebral%',
            '%neutralization%',
            '% ab %',
            '%metapneumovirus%'
        ]) AS excluded_pattern
        WHERE LOWER(c.concept_name) LIKE excluded_pattern
    )
ORDER BY c.concept_name
"""

adv_tests = polars_gbq(adv_tests_q)


# In[ ]:


m_adv_tests = f"""
        SELECT
            m.person_id,
            m.measurement_concept_id as concept_id,
            c.concept_name as concept_name,
            c.concept_code as code,
            c.vocabulary_id as vocab,
            m.measurement_datetime as datetime,
            c2.concept_name AS result,
            CAST(m.visit_occurrence_id AS STRING) as visit_id,
            v.visit_start_datetime as visit_start,
            c1.concept_name AS visit_name
        FROM
            {version}.measurement m
        INNER JOIN
            {version}.concept c ON m.measurement_concept_id = c.concept_id
        LEFT JOIN
            {version}.concept c2 ON m.value_as_concept_id = c2.concept_id
        LEFT JOIN
            {version}.visit_occurrence v ON m.visit_occurrence_id = v.visit_occurrence_id AND m.person_id = v.person_id
        LEFT JOIN
            {version}.concept c1 ON v.visit_concept_id = c1.concept_id      
        WHERE
            m.measurement_concept_id IN (SELECT DISTINCT concept_id FROM ({adv_tests_q}))
        ORDER BY 
            person_id, measurement_datetime
        """

adv_test_df = polars_gbq(m_adv_tests)
len(adv_test_df)
# 18741 v7
# 37029 v8


# In[ ]:


adv_test_df['person_id'].n_unique()
# 8675 v7
# 17362 v8


# In[ ]:


adv_test_df.group_by(['concept_id', 'concept_name', 'result']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(10)


# In[ ]:


adv_test_counts = adv_test_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

adv_tests = adv_tests.join(
    adv_test_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


adv_tests.slice(0,50)


# ## RSV

# ### RSV Diagnoses

# In[ ]:


rsv_exclusion = [
    '%vaccin%', 
    '%non-respiratory%',
    '%otitis%',
    '%myocarditis%',
    '%adenovirus%',
    '%bovine%',
    '%myelitis%',
    '%gastrointestinal%',
    '%gastroenteritis%',
    '%central nervous system%',
    '%cns disorder%',
    '%encephalopathy%',
    '%encephalitis%',
    '%virus inconclusive%',
    '%not isolated%',
    '%not detected%',
    '%equine%'
]

# Querying the concept codes
rsv_codes_q = construct_concept_query(version, 'syncytial', rsv_exclusion)


# In[ ]:


rsv_codes = polars_gbq(rsv_codes_q)


# In[ ]:


# fetch_co(version, rsv_codes_q, return_count=True)
# 2553 v7
# 3502 v8


# In[ ]:


# fetch_o(version, rsv_codes_q, return_count=True)
# 0 v7, v8


# In[ ]:


rsv_code_df = fetch_co(version, rsv_codes_q, return_count=False)


# In[ ]:


rsv_code_df['person_id'].n_unique()
# 880 v7
# 1312 v8


# In[ ]:


rsv_code_counts = rsv_code_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

rsv_codes = rsv_codes.join(
    rsv_code_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


rsv_codes


# ### RSV Tests

# In[ ]:


# Search panels for RSV tests
rsv_tests_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    c.vocabulary_id = 'LOINC'
    AND LOWER(c.concept_name) LIKE '%syncytial%'
    AND (
        LOWER(c.concept_name) LIKE '%[presence]%'
        OR LOWER(c.concept_name) LIKE '%[identifier]%'
        OR LOWER(c.concept_name) LIKE '%[cycle threshold]%'
        OR LOWER(c.concept_name) LIKE '%[susceptibility]%'
        OR LOWER(c.concept_name) LIKE '%viral load%'
        OR LOWER(c.concept_name) LIKE '%culture%'
        OR c.concept_name LIKE '%NAA%'
        OR c.concept_name LIKE '%Measurement of%'
    )
    AND NOT EXISTS (
        SELECT 1
        FROM UNNEST(ARRAY[
            '%porcine%', 
            '%canine%',
            '%bovine%',
            '%influenza%',
            '%cornea%',
            '%igg%',
            '%igm%',
            '%cerebral%',
            '% ab %',
            '%metapneumovirus%',
            '%coronavirus%'
        ]) AS excluded_pattern
        WHERE LOWER(c.concept_name) LIKE excluded_pattern
    )
ORDER BY c.concept_name
"""

rsv_tests = polars_gbq(rsv_tests_q)


# In[ ]:


m_rsv_tests = f"""
        SELECT
            m.person_id,
            m.measurement_concept_id as concept_id,
            c.concept_name as concept_name,
            c.concept_code as code,
            c.vocabulary_id as vocab,
            m.measurement_datetime as datetime,
            c2.concept_name AS result,
            CAST(m.visit_occurrence_id AS STRING) as visit_id,
            v.visit_start_datetime as visit_start,
            c1.concept_name AS visit_name
        FROM
            {version}.measurement m
        INNER JOIN
            {version}.concept c ON m.measurement_concept_id = c.concept_id
        LEFT JOIN
            {version}.concept c2 ON m.value_as_concept_id = c2.concept_id
        LEFT JOIN
            {version}.visit_occurrence v ON m.visit_occurrence_id = v.visit_occurrence_id AND m.person_id = v.person_id
        LEFT JOIN
            {version}.concept c1 ON v.visit_concept_id = c1.concept_id      
        WHERE
            m.measurement_concept_id IN (SELECT DISTINCT concept_id FROM ({rsv_tests_q}))
        ORDER BY 
            person_id, measurement_datetime
        """

rsv_test_df = polars_gbq(m_rsv_tests)
len(rsv_test_df)
# 30206 v7
# 73417 v8


# In[ ]:


rsv_test_df['person_id'].n_unique()
# 12529 v7
# 29183 v8


# In[ ]:


rsv_test_df.group_by(['concept_id', 'concept_name', 'result']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(10)


# In[ ]:


rsv_test_df.head()


# In[ ]:


rsv_test_counts = rsv_test_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

rsv_tests = rsv_tests.join(
    rsv_test_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


rsv_tests.head(10)


# ## Coronaviruses

# ### Coronavirus Diagnoses

# In[ ]:


hcov_exclusion = [
    '%vaccin%', 
    '%non-respiratory%',
    '%acute respiratory syndrome coronavirus 2%',
    '%acute respiratory syndrome coronavirus%',
    '%sars%',
    '%2019%',
    '%otitis%',
    '%myocarditis%',
    '%myelitis%',
    '%gastrointestinal%',
    '%gastroenteritis%',
    '%central nervous system%',
    '%cns disorder%',
    '%encephalopathy%',
    '%encephalitis%',
    '%virus inconclusive%',
    '%not isolated%',
    '%not detected%',
    '%equine%'
]

# Querying the concept codes
hcov_codes_q = construct_concept_query(version, 'coronavirus', hcov_exclusion)


# In[ ]:


hcov_codes = polars_gbq(hcov_codes_q)


# In[ ]:


# fetch_co(version, hcov_codes_q, return_count=True)
# 2152 v7
# 3124 v8


# In[ ]:


# fetch_o(version, hcov_codes_q, return_count=True)
#0 v7, v8


# In[ ]:


hcov_code_df = fetch_co(version, hcov_codes_q, return_count=False)


# In[ ]:


hcov_code_df['person_id'].n_unique()
# 1091 v7
# 1694 v8


# In[ ]:


hcov_code_counts = hcov_code_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

hcov_codes = hcov_codes.join(
    hcov_code_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


hcov_codes


# ### Coronavirus Tests

# In[ ]:


# Search panels for coronavirus tests
hcov_tests_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    c.vocabulary_id = 'LOINC'
    AND LOWER(c.concept_name) LIKE '%coronavirus%'
    AND (
        LOWER(c.concept_name) LIKE '%[presence]%'
        OR LOWER(c.concept_name) LIKE '%[identifier]%'
        OR LOWER(c.concept_name) LIKE '%[cycle threshold]%'
        OR LOWER(c.concept_name) LIKE '%[susceptibility]%'
        OR LOWER(c.concept_name) LIKE '%viral load%'
        OR LOWER(c.concept_name) LIKE '%culture%'
        OR c.concept_name LIKE '%NAA%'
        OR c.concept_name LIKE '%Measurement of%'
    )
    AND NOT EXISTS (
        SELECT 1
        FROM UNNEST(ARRAY[
            '%turkey%', 
            '%sars%', 
            '%mers%', 
            '%porcine%', 
            '%canine%',
            '%bovine%',
            '%haemophilus%',
            '%cornea%',
            '%igg%',
            '%igm%',
            '%cerebral%',
            '% ab %'
        ]) AS excluded_pattern
        WHERE LOWER(c.concept_name) LIKE excluded_pattern
    )
ORDER BY c.concept_name
"""

hcov_tests = polars_gbq(hcov_tests_q)


# In[ ]:


m_hcov_tests = f"""
        SELECT
            m.person_id,
            m.measurement_concept_id as concept_id,
            c.concept_name as concept_name,
            c.concept_code as code,
            c.vocabulary_id as vocab,
            m.measurement_datetime as datetime,
            c2.concept_name AS result,
            CAST(m.visit_occurrence_id AS STRING) as visit_id,
            v.visit_start_datetime as visit_start,
            c1.concept_name AS visit_name
        FROM
            {version}.measurement m
        INNER JOIN
            {version}.concept c ON m.measurement_concept_id = c.concept_id
        LEFT JOIN
            {version}.concept c2 ON m.value_as_concept_id = c2.concept_id
        LEFT JOIN
            {version}.visit_occurrence v ON m.visit_occurrence_id = v.visit_occurrence_id AND m.person_id = v.person_id
        LEFT JOIN
            {version}.concept c1 ON v.visit_concept_id = c1.concept_id      
        WHERE
            m.measurement_concept_id IN (SELECT DISTINCT concept_id FROM ({hcov_tests_q}))
        ORDER BY 
            person_id, measurement_datetime
        """

hcov_test_df = polars_gbq(m_hcov_tests)
len(hcov_test_df)
# 51934 v7
# 107347 v8


# In[ ]:


hcov_test_df['person_id'].n_unique()
# 6097 v7
# 13368 v8


# In[ ]:


hcov_test_df.group_by(['concept_id', 'concept_name', 'result']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(10)


# In[ ]:


hcov_test_counts = hcov_test_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

hcov_tests = hcov_tests.join(
    hcov_test_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


hcov_tests.head(10)


# ## Influenza

# ### Influenza Diagnoses

# In[ ]:


flu_exclusion = [
    '%haemophilus%', 
    '%hemophilus%', 
    '%parainfluenza%', 
    '%vaccin%', 
    '%at risk of influenza%', 
    '%post-influenza%', 
    '%pneumonia or influenza nos%',
    '%&/or%',
    '%influenza-like%',
    '%influenza like illness%',
    '%non-respiratory%',
    '%otitis%',
    '%myocarditis%',
    '%myelitis%',
    '%gastrointestinal%',
    '%gastroenteritis%',
    '%central nervous system%',
    '%cns disorder%',
    '%encephalopathy%',
    '%encephalitis%',
    '%virus inconclusive%',
    '%not isolated%',
    '%not detected%',
    '%equine%'
]

# flu_exclusion = ['%equine%']

# Querying the concept codes
flu_codes_q = construct_concept_query(version, 'influenza', flu_exclusion)


# In[ ]:


flu_codes = polars_gbq(flu_codes_q)


# In[ ]:


# fetch_co(version, flu_codes_q, return_count=True)
# 27229 v7
# 39856 v8


# In[ ]:


# fetch_o(version, flu_codes_q, return_count=True)
#0 v7, v8


# In[ ]:


flu_code_df = fetch_co(version, flu_codes_q, return_count=False)


# In[ ]:


flu_code_df['person_id'].n_unique()
# 12431 v7
# 17456 v8


# In[ ]:


flu_code_counts = flu_code_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

flu_codes = flu_codes.join(
    flu_code_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


flu_codes.head(10)


# ### Influenza Tests

# In[ ]:


# Search panels for influenza tests
flu_tests_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    c.vocabulary_id = 'LOINC'
    AND LOWER(c.concept_name) LIKE '%influenza%'
    AND (
        LOWER(c.concept_name) LIKE '%[presence]%'
        OR LOWER(c.concept_name) LIKE '%[identifier]%'
        OR LOWER(c.concept_name) LIKE '%[cycle threshold]%'
        OR LOWER(c.concept_name) LIKE '%[susceptibility]%'
        OR LOWER(c.concept_name) LIKE '%viral load%'
        OR LOWER(c.concept_name) LIKE '%culture%'
        OR c.concept_name LIKE '%NAA%'
        OR c.concept_name LIKE '%Measurement of%'
    )
    AND NOT EXISTS (
        SELECT 1
        FROM UNNEST(ARRAY[
            '%canine%',
            '%bovine%',
            '%equine%',
            '%adamantane%',
            '%neuraminidase%',
            '%parainfluenza%',
            '%haemophilus%',
            '%tissue%',
            '%pericardial%',
            '%cornea%',
            '%serum%',
            '%igg%',
            '%igm%',
            '%sars%',
            '%cerebral%',
            '% ab %',
            '%syncytial%'
        ]) AS excluded_pattern
        WHERE LOWER(c.concept_name) LIKE excluded_pattern
    )
ORDER BY c.concept_name
"""

flu_tests = polars_gbq(flu_tests_q)


# In[ ]:


m_flu_tests = f"""
        SELECT
            m.person_id,
            m.measurement_concept_id as concept_id,
            c.concept_name as concept_name,
            c.concept_code as code,
            c.vocabulary_id as vocab,
            m.measurement_datetime as datetime,
            c2.concept_name AS result,
            CAST(m.visit_occurrence_id AS STRING) as visit_id,
            v.visit_start_datetime as visit_start,
            c1.concept_name AS visit_name
        FROM
            {version}.measurement m
        INNER JOIN
            {version}.concept c ON m.measurement_concept_id = c.concept_id
        LEFT JOIN
            {version}.concept c2 ON m.value_as_concept_id = c2.concept_id
        LEFT JOIN
            {version}.visit_occurrence v ON m.visit_occurrence_id = v.visit_occurrence_id AND m.person_id = v.person_id
        LEFT JOIN
            {version}.concept c1 ON v.visit_concept_id = c1.concept_id      
        WHERE
            m.measurement_concept_id IN (SELECT DISTINCT concept_id FROM ({flu_tests_q}))
        ORDER BY 
            person_id, measurement_datetime
        """

flu_test_df = polars_gbq(m_flu_tests)
len(flu_test_df)
# 167270 v7
# 349728 v8


# In[ ]:


flu_test_df['person_id'].n_unique()
# 33655 v7
# 60805 v8


# In[ ]:


flu_test_df.group_by(['concept_id', 'concept_name', 'result']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(10)


# In[ ]:


flu_test_counts = flu_test_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

flu_tests = flu_tests.join(
    flu_test_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


flu_tests.head(10)


# ### Influenza Treatments

# In[ ]:


# Influenza drugs: oseltamivir
ose_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.vocabulary_id,
    c.domain_id,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    LOWER(c.concept_name) LIKE '%oseltamivir%'
    AND c.domain_id = 'Drug'
ORDER BY c.concept_name DESC
"""

polars_gbq(f"""SELECT COUNT(*) FROM ({ose_q})""")
# 1602 v7
# 1742 v8


# In[ ]:


# Influenza drugs: zanamivir
zan_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.vocabulary_id,
    c.domain_id,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    LOWER(c.concept_name) LIKE '%zanamivir%'
    AND c.domain_id = 'Drug'
ORDER BY c.concept_name DESC
"""

polars_gbq(f"""SELECT COUNT(*) FROM ({zan_q})""")
# 162 v7
# 164 v8


# In[ ]:


# Influenza drugs: baloxavir
bal_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.vocabulary_id,
    c.domain_id,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    LOWER(c.concept_name) LIKE '%baloxavir%'
    AND c.domain_id = 'Drug'
ORDER BY c.concept_name DESC
"""

polars_gbq(f"""SELECT COUNT(*) FROM ({bal_q})""")
# 88 v7
# 93 v8


# In[ ]:


de_flu_drugs = f"""
    SELECT DISTINCT
        de.person_id,
        drug_concept_id as concept_id,
        c1.concept_name AS concept_name,
        c1.concept_code as code,
        c1.vocabulary_id as vocab,
        drug_exposure_start_datetime as datetime,
        drug_exposure_end_datetime as end_datetime,
        days_supply,
        CAST (de.visit_occurrence_id AS STRING) as visit_id,
        v.visit_start_datetime as visit_start,
        c2.concept_name AS visit_name    
    FROM
        {version}.drug_exposure de
    INNER JOIN
        {version}.concept c1 ON de.drug_concept_id = c1.concept_id
    LEFT JOIN
        {version}.visit_occurrence v ON de.visit_occurrence_id = v.visit_occurrence_id AND de.person_id = v.person_id
    LEFT JOIN
        {version}.concept c2 ON v.visit_concept_id = c2.concept_id
    WHERE 
        (de.drug_concept_id IN (SELECT concept_id FROM ({ose_q})) AND
        (de.days_supply > 1 
        OR TIMESTAMP_DIFF(de.drug_exposure_end_datetime, de.drug_exposure_start_datetime, DAY) > 1.75
        OR EXISTS (
            SELECT 1
            FROM {version}.drug_exposure de2
            WHERE de2.person_id = de.person_id
                AND de2.drug_concept_id IN (SELECT concept_id FROM ({ose_q}))
                AND (
                    de2.drug_exposure_start_datetime BETWEEN DATE_ADD(de.drug_exposure_start_datetime, INTERVAL 1 DAY) AND DATE_ADD(de.drug_exposure_start_datetime, INTERVAL 5 DAY)
                    OR
                    de2.drug_exposure_start_datetime BETWEEN DATE_SUB(de.drug_exposure_start_datetime, INTERVAL 5 DAY) AND DATE_SUB(de.drug_exposure_start_datetime, INTERVAL 1 DAY)
                )
            )
        )
        ) 
        OR
        (de.drug_concept_id IN (SELECT concept_id FROM ({zan_q})) AND
        (de.days_supply > 1 
        OR TIMESTAMP_DIFF(de.drug_exposure_end_datetime, de.drug_exposure_start_datetime, DAY) > 1.75
        OR EXISTS (
            SELECT 1
            FROM {version}.drug_exposure de2
            WHERE de2.person_id = de.person_id
                AND de2.drug_concept_id IN (SELECT concept_id FROM ({zan_q}))
                AND (
                    de2.drug_exposure_start_datetime BETWEEN DATE_ADD(de.drug_exposure_start_datetime, INTERVAL 1 DAY) AND DATE_ADD(de.drug_exposure_start_datetime, INTERVAL 5 DAY)
                    OR
                    de2.drug_exposure_start_datetime BETWEEN DATE_SUB(de.drug_exposure_start_datetime, INTERVAL 5 DAY) AND DATE_SUB(de.drug_exposure_start_datetime, INTERVAL 1 DAY)
                )
            )
        )
        ) 
        OR
        (de.drug_concept_id IN (SELECT concept_id FROM ({bal_q}))
        )
    ORDER BY person_id, drug_exposure_start_datetime
"""
flu_drug_df = polars_gbq(de_flu_drugs)
len(flu_drug_df)
# 26934 v7
# 35796 v8


# In[ ]:


flu_drug_df['person_id'].n_unique()
# 12989 v7
# 17759 v8


# In[ ]:


# flu_drugs = pl.concat([ose, zan, bal])


# In[ ]:


# flu_drug_counts = flu_drug_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

# flu_drugs = flu_drugs.join(
#     flu_drug_counts,
#     on=['concept_id', 'concept_name'],
#     how='left'
# ).sort('concept_name', descending=False)


# In[ ]:


conditions = ['oseltamivir', 'zanamivir', 'baloxavir']
unique_counts = {}

for condition in conditions:
    # Filter rows where concept_name contains the condition
    condition_df = flu_drug_df.filter(pl.col('concept_name').str.contains(condition))

    # Count unique person_id's
    unique_counts[condition] = condition_df.select('person_id').unique().shape[0]

unique_counts
# {'oseltamivir': 12915, 'zanamivir': 31, 'baloxavir': 85} v7
# {'oseltamivir': 17647, 'zanamivir': 35, 'baloxavir': 144} v8


# In[ ]:


# Calculate the difference in days between end_datetime and start_datetime
flu_drug_df = flu_drug_df.with_columns(
    ((pl.col("end_datetime") - pl.col("datetime")).dt.total_seconds() / (24*3600)).alias("duration_days")
).with_columns(
    pl.col("duration_days").round(decimals=1)
)
# Reorder columns to place duration_days before days_supply
flu_drug_df = flu_drug_df.select([
    "person_id", "concept_id", "concept_name", "code", "vocab",
    "datetime", "end_datetime", "duration_days", "days_supply",
    "visit_id", "visit_start", "visit_name"
])


# In[ ]:


# Filtering out oseltamivir and zanamivir duration_days 7+ or days_supply > 5 when duration_days is not between 3-6.
# This was selected because numerous days_supply = 10 had duration_days = 5, which would be an appropriate treatment.
flu_drug_df = flu_drug_df.filter(
    ~(
        (pl.col('concept_name').str.contains('oseltamivir') | pl.col('concept_name').str.contains('zanamivir')) & 
        (
            (pl.col('duration_days') > 6) | 
            ((pl.col('days_supply') > 5) & ~pl.col('duration_days').is_between(3, 6))
        )
    )
)


# In[ ]:


len(flu_drug_df)
# 11619 v7
# 15716 v8


# In[ ]:


flu_drug_df['person_id'].n_unique()
# 7460 v7
# 10136 v8


# In[ ]:


conditions = ['oseltamivir', 'zanamivir', 'baloxavir']
unique_counts = {}

for condition in conditions:
    # Filter rows where concept_name contains the condition
    condition_df = flu_drug_df.filter(pl.col('concept_name').str.contains(condition))

    # Count unique person_id's
    unique_counts[condition] = condition_df.select('person_id').unique().shape[0]

unique_counts
# {'oseltamivir': 7393, 'zanamivir': 8, 'baloxavir': 85} v7
# {'oseltamivir': 10021, 'zanamivir': 11, 'baloxavir': 144} v8


# In[ ]:


flu_drug_df = flu_drug_df.drop(['end_datetime', 'duration_days', 'days_supply'])


# In[ ]:


flu_drug_df.head()


# In[ ]:


flu_drug_df.group_by(['concept_id', 'code', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)


# ## COVID

# ### COVID Diagnoses

# In[ ]:


#ncov and hcov-19 yielded no additional COVID codes
covid_inclusion = [
    'sars-cov-2',
    'covid',
    'severe acute respiratory syndrome coronavirus'
]

covid_exclusion = [
    '%vaccin%', 
    '%immunity%', 
    '%non-respiratory%',
    '%otitis%',
    '%myocarditis%',
    '%myelitis%',
    '%cardiomyopathy%',
    '%kidney injury%',
    '%conjunctivitis%',
    '%rhabdomyolysis%',
    '%gastrointestinal%',
    '%gastroenteritis%',
    '%cytopenia%',
    '%exposure%',
    '%risk category for developing complication%',
    '%central nervous system%',
    '%cns disorder%',
    '%encephalopathy%',
    '%encephalitis%',
    '%virus inconclusive%',
    '%result unknown%',
    '%result negative%',
    '%not isolated%',
    '%not detected%',
    '%post%',
    '%indeterminate%',
    '%equivocal%',
    '% igm %',
    '% igg %',
    '% iga %',
    '%antibody%',
    '%equine%'
]

# Querying the concept codes
covid_codes_q = construct_concept_query(version, covid_inclusion, covid_exclusion)


# In[ ]:


covid_codes = polars_gbq(covid_codes_q)


# In[ ]:


covid_codes


# In[ ]:


# fetch_co(version, covid_codes_q, return_count=True)
# 40604 v7
# 88855 v8


# In[ ]:


# fetch_o(version, covid_codes_q, return_count=True)
# 13192 v7
# 14735 v8


# In[ ]:


covid_code_co_df = fetch_co(version, covid_codes_q, return_count=False)


# In[ ]:


covid_code_co_df['person_id'].n_unique()
# 17252 v7
# 37272 v8


# In[ ]:


covid_code_o_df = fetch_o(version, covid_codes_q, return_count=False)


# In[ ]:


covid_code_o_df['person_id'].n_unique()
# 2___ v7
# 5773 v8


# In[ ]:


covid_code_df = pl.concat([covid_code_co_df, covid_code_o_df])


# In[ ]:


covid_code_df['person_id'].n_unique()
# 16974 v7
# 37435 v8


# In[ ]:


covid_code_counts = covid_code_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

covid_codes = covid_codes.join(
    covid_code_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


covid_codes


# ### COVID Tests

# In[ ]:


# Search panels for COVID tests
covid_tests_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    c.vocabulary_id = 'LOINC'
    AND (
        LOWER(c.concept_name) LIKE '%covid%'
        OR LOWER(c.concept_name) LIKE '%sars-cov-2%'
        OR LOWER(c.concept_name) LIKE '%sars-related coronavirus%' -- removed ' rna'
        OR LOWER(c.concept_name) LIKE 'severe acute respiratory syndrome coronavirus' -- added
    )
    AND (
        LOWER(c.concept_name) LIKE '%[presence]%'
        OR LOWER(c.concept_name) LIKE '%[identifier]%'
        OR LOWER(c.concept_name) LIKE '%[cycle threshold%'
        OR LOWER(c.concept_name) LIKE '%[susceptibility]%'
        OR LOWER(c.concept_name) LIKE '%clade [type]%' 
        OR LOWER(c.concept_name) LIKE '%viral load%'
        OR LOWER(c.concept_name) LIKE '%culture%'
        OR c.concept_name LIKE '%NAA%'
        OR LOWER(c.concept_name) LIKE '%measurement of%'
        OR c.concept_name LIKE '%2019-ncov%'
        OR LOWER(c.concept_name) LIKE '%antigen%' -- added
    )
    AND NOT EXISTS (
        SELECT 1
        FROM UNNEST(ARRAY[
            '%porcine%', 
            '%canine%',
            '%bovine%',
            '%influenza%',
            '%haemophilus%',
          --  '%blood%',
          --  '%serum%',
            '%cornea%',
            '%igg%',
            '%igm%',
            '%cerebral%',
            '%donor%',
            '%antibody%',
            '% ab %',
            '%volatile%',
            '%stimulated gamma interferon%'
        ]) AS excluded_pattern
        WHERE LOWER(c.concept_name) LIKE excluded_pattern
    )
ORDER BY c.concept_name
"""

covid_tests = polars_gbq(covid_tests_q)


# In[ ]:


# Note, this does not detect 756055 (OMOP4873969 -- OMOP Extension -- Measurement of Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2))<br>
# It captures 100 positive cases, but given no other use of OMOP extensions alone, will not include.
m_covid_tests = f"""
        SELECT
            m.person_id,
            m.measurement_concept_id as concept_id,
            c.concept_name as concept_name,
            c.concept_code as code,
            c.vocabulary_id as vocab,
            m.measurement_datetime as datetime,
            c2.concept_name AS result,
            CAST(m.visit_occurrence_id AS STRING) as visit_id,
            v.visit_start_datetime as visit_start,
            c1.concept_name AS visit_name
        FROM
            {version}.measurement m
        INNER JOIN
            {version}.concept c ON m.measurement_concept_id = c.concept_id
        LEFT JOIN
            {version}.concept c2 ON m.value_as_concept_id = c2.concept_id
        LEFT JOIN
            {version}.visit_occurrence v ON m.visit_occurrence_id = v.visit_occurrence_id AND m.person_id = v.person_id
        LEFT JOIN
            {version}.concept c1 ON v.visit_concept_id = c1.concept_id      
        WHERE
            m.measurement_concept_id IN (SELECT DISTINCT concept_id FROM ({covid_tests_q}))
        ORDER BY 
            person_id, measurement_datetime
        """

covid_test_df = polars_gbq(m_covid_tests)
len(covid_test_df)
# 275422 v7
# 483390 v8
# 483392 v8 N3C concept update


# In[ ]:


covid_test_df['person_id'].n_unique()
# 93903 v7
# 138161 v8
# 138161 v8 N3C concept update


# In[ ]:


covid_test_df.group_by(['concept_id', 'concept_name', 'result']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(10)


# In[ ]:


covid_test_counts = covid_test_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)

covid_tests = covid_tests.join(
    covid_test_counts,
    on=['concept_id', 'concept_name'],
    how='left'
).sort('concept_code', descending=False)


# In[ ]:


covid_tests.slice(0,50)


# ### COVID Treatments

# In[ ]:


# COVID drugs: remdesivir
rem_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.vocabulary_id,
    c.domain_id,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    LOWER(c.concept_name) LIKE '%remdesivir%'
    AND c.domain_id = 'Drug'
ORDER BY c.concept_name DESC
"""

polars_gbq(f"""SELECT COUNT(*) FROM ({rem_q})""")
# 115 v7
# 125 v8


# In[ ]:


# COVID drugs: nirmatrelvir
nir_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.vocabulary_id,
    c.domain_id,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    LOWER(c.concept_name) LIKE '%nirmatrelvir%'
    AND c.domain_id = 'Drug'
ORDER BY c.concept_name DESC
"""

polars_gbq(f"""SELECT COUNT(*) FROM ({nir_q})""")
# 26 v7
# 54 v8


# In[ ]:


# COVID drugs: molnupiravir
mol_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.vocabulary_id,
    c.domain_id,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    LOWER(c.concept_name) LIKE '%molnupiravir%'
    AND c.domain_id = 'Drug'
ORDER BY c.concept_name DESC
"""

polars_gbq(f"""SELECT COUNT(*) FROM ({mol_q})""")
# 13 v7
# 31 v8


# In[ ]:


de_covid_drugs = f"""
    SELECT DISTINCT
        de.person_id,
        drug_concept_id as concept_id,
        c1.concept_name AS concept_name,
        c1.concept_code as code,
        c1.vocabulary_id as vocab,
        drug_exposure_start_datetime as datetime,
        CAST (de.visit_occurrence_id AS STRING) as visit_id,
        v.visit_start_datetime as visit_start,
        c2.concept_name AS visit_name
    FROM
        {version}.drug_exposure de
    INNER JOIN
        {version}.concept c1 ON de.drug_concept_id = c1.concept_id
    LEFT JOIN
        {version}.visit_occurrence v ON de.visit_occurrence_id = v.visit_occurrence_id AND de.person_id = v.person_id
    LEFT JOIN
        {version}.concept c2 ON v.visit_concept_id = c2.concept_id
    WHERE 
        (de.drug_concept_id IN (SELECT concept_id FROM ({rem_q})) AND
        (de.days_supply > 1 
        OR TIMESTAMP_DIFF(de.drug_exposure_end_datetime, de.drug_exposure_start_datetime, DAY) > 1.75
        OR EXISTS (
            SELECT 1
            FROM {version}.drug_exposure de2
            WHERE de2.person_id = de.person_id
                AND de2.drug_concept_id IN (SELECT concept_id FROM ({rem_q}))
                AND (
                    de2.drug_exposure_start_datetime BETWEEN DATE_ADD(de.drug_exposure_start_datetime, INTERVAL 1 DAY) AND DATE_ADD(de.drug_exposure_start_datetime, INTERVAL 5 DAY)
                    OR
                    de2.drug_exposure_start_datetime BETWEEN DATE_SUB(de.drug_exposure_start_datetime, INTERVAL 5 DAY) AND DATE_SUB(de.drug_exposure_start_datetime, INTERVAL 1 DAY)
                )
            )
        )
        ) 
        OR
        (de.drug_concept_id IN (SELECT concept_id FROM ({nir_q})) AND
        (de.days_supply > 1 
        OR TIMESTAMP_DIFF(de.drug_exposure_end_datetime, de.drug_exposure_start_datetime, DAY) > 1.75
        OR EXISTS (
            SELECT 1
            FROM {version}.drug_exposure de2
            WHERE de2.person_id = de.person_id
                AND de2.drug_concept_id IN (SELECT concept_id FROM ({nir_q}))
                AND (
                    de2.drug_exposure_start_datetime BETWEEN DATE_ADD(de.drug_exposure_start_datetime, INTERVAL 1 DAY) AND DATE_ADD(de.drug_exposure_start_datetime, INTERVAL 5 DAY)
                    OR
                    de2.drug_exposure_start_datetime BETWEEN DATE_SUB(de.drug_exposure_start_datetime, INTERVAL 5 DAY) AND DATE_SUB(de.drug_exposure_start_datetime, INTERVAL 1 DAY)
                )
            )
        )
        )
        OR
        (de.drug_concept_id IN (SELECT concept_id FROM ({mol_q})) AND
        (de.days_supply > 1 
        OR TIMESTAMP_DIFF(de.drug_exposure_end_datetime, de.drug_exposure_start_datetime, DAY) > 1.75
        OR EXISTS (
            SELECT 1
            FROM {version}.drug_exposure de2
            WHERE de2.person_id = de.person_id
                AND de2.drug_concept_id IN (SELECT concept_id FROM ({mol_q}))
                AND (
                    de2.drug_exposure_start_datetime BETWEEN DATE_ADD(de.drug_exposure_start_datetime, INTERVAL 1 DAY) AND DATE_ADD(de.drug_exposure_start_datetime, INTERVAL 5 DAY)
                    OR
                    de2.drug_exposure_start_datetime BETWEEN DATE_SUB(de.drug_exposure_start_datetime, INTERVAL 5 DAY) AND DATE_SUB(de.drug_exposure_start_datetime, INTERVAL 1 DAY)
                )
            )
        )
        )
    ORDER BY person_id, drug_exposure_start_datetime
"""

covid_drug_df = polars_gbq(de_covid_drugs)
len(covid_drug_df)
# 2615 v7
# 15117 v8


# In[ ]:


covid_drug_df['person_id'].n_unique()
# 1329 v7
# 10364 v8


# In[ ]:


covid_drug_df.group_by(['concept_id', 'concept_name', 'code']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(50)


# In[ ]:


conditions = ['remdesivir', 'nirmatrelvir', 'molnupiravir']
unique_counts = {}

for condition in conditions:
    # Filter rows where concept_name contains the condition
    condition_df = covid_drug_df.filter(pl.col('concept_name').str.contains(condition))

    # Count unique person_id's
    unique_counts[condition] = condition_df.select('person_id').unique().shape[0]

unique_counts
# {'remdesivir': 873, 'nirmatrelvir': 426, 'molnupiravir': 36} v7
# {'remdesivir': 1679, 'nirmatrelvir': 7855, 'molnupiravir': 952} v8


# In[ ]:


covid_drug_df.group_by(['concept_id', 'code', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True)


# ## Bordetella pertussus/parapertussis

# ### Pertussis Diagnoses

# In[ ]:


# Pertussis codes
pert_codes_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.vocabulary_id,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    c.vocabulary_id LIKE '%ICD%'
    AND c.domain_id = 'Condition'
    AND (LOWER(c.concept_name) LIKE '%pertussis%' OR LOWER(c.concept_name) LIKE '%bordetella%')
    AND NOT EXISTS (
        SELECT 1
        FROM UNNEST(ARRAY[
            '%vaccin%', 
            '%non-respiratory%',
            '%otitis%',
            '%meningitis%',
            '%myocarditis%',
            '%myelitis%',
            '%conjunctivitis%',
            '%gastrointestinal%',
            '%gastroenteritis%',
            '%central nervous system%',
            '%cns disorder%',
            '%encephalopathy%',
            '%encephalitis%',
            '%bronchiseptica%',
            '%immune%',
            '%not isolated%',
            '%not detected%'
        ]) AS excluded_pattern
        WHERE LOWER(c.concept_name) LIKE excluded_pattern
    )
ORDER BY c.concept_name DESC
"""

pert_codes = polars_gbq(f"""({pert_codes_q}) LIMIT 150""")


# In[ ]:


pert_codes


# In[ ]:


# fetch_co(version, pert_codes_q, return_count=True)
# 43 v7
# 396 v8


# In[ ]:


# fetch_o(version, pert_codes_q, return_count=True)
#0 v7, v8


# In[ ]:


pert_code_df = fetch_co(version, pert_codes_q, return_count=False)


# In[ ]:


pert_code_df['person_id'].n_unique()
# 31 v7
# 30 v8


# In[ ]:


pert_code_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(10)


# ### Pertussis Tests

# In[ ]:


# Search panels for pertussis tests
# looking at panel descendants (below), 3033775 [Bordetella sp indetified...] was missed. This is now added manually in the WHERE clause.

pert_tests_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    c.vocabulary_id = 'LOINC'
    AND (LOWER(c.concept_name) LIKE '%pertussis%' OR LOWER(c.concept_name) LIKE '%bordetella sp identified in specimen by organism specific culture%')
    AND (
        LOWER(c.concept_name) LIKE '%[presence]%'
        OR LOWER(c.concept_name) LIKE '%[identifier]%'
        OR LOWER(c.concept_name) LIKE '%[cycle threshold]%'
        OR LOWER(c.concept_name) LIKE '%[susceptibility]%'
        OR LOWER(c.concept_name) LIKE '%viral load%'
        OR LOWER(c.concept_name) LIKE '%culture%'
        OR c.concept_name LIKE '%NAA%'
        OR c.concept_name LIKE '%Measurement of%'
    )
    AND NOT EXISTS (
        SELECT 1
        FROM UNNEST(ARRAY[
            '%porcine%', 
            '%canine%',
            '%bovine%',
            '%cornea%',
            '%bronchiseptica%',
            '%holmesii%',
            '%igg%',
            '%igm%',
            '%cerebral%',
            '% ab %'
        ]) AS excluded_pattern
        WHERE LOWER(c.concept_name) LIKE excluded_pattern
    )
ORDER BY c.concept_name
"""

# polars_gbq(pert_tests_q)


# In[ ]:


m_pert_tests = f"""
        SELECT
            m.person_id,
            m.measurement_concept_id as concept_id,
            c.concept_name as concept_name,
            c.concept_code as code,
            c.vocabulary_id as vocab,
            m.measurement_datetime as datetime,
            c2.concept_name AS result,
            CAST(m.visit_occurrence_id AS STRING) as visit_id,
            v.visit_start_datetime as visit_start,
            c1.concept_name AS visit_name
        FROM
            {version}.measurement m
        INNER JOIN
            {version}.concept c ON m.measurement_concept_id = c.concept_id
        LEFT JOIN
            {version}.concept c2 ON m.value_as_concept_id = c2.concept_id
        LEFT JOIN
            {version}.visit_occurrence v ON m.visit_occurrence_id = v.visit_occurrence_id AND m.person_id = v.person_id
        LEFT JOIN
            {version}.concept c1 ON v.visit_concept_id = c1.concept_id      
        WHERE
            m.measurement_concept_id IN (SELECT DISTINCT concept_id FROM ({pert_tests_q}))
        ORDER BY 
            person_id, measurement_datetime
        """

pert_test_df = polars_gbq(m_pert_tests)
len(pert_test_df)
# 29303 v7
# 50152 v8


# In[ ]:


pert_test_df['person_id'].n_unique()
# 10902 v7
# 17520 v8


# In[ ]:


pert_test_df.group_by(['concept_id', 'concept_name', 'result']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(10)


# ## Mycoplasma pneumoniae

# ### M. pneumoniae Diagnoses

# In[ ]:


mpa_exclusion = [
    '%non-respiratory%',
    '%otitis%',
    '%myelitis%',
    '%gastrointestinal%',
    '%gastroenteritis%',
    '%central nervous system%',
    '%cns disorder%',
    '%encephalopathy%',
    '%encephalitis%',
    '%hominis%',
    '%pyelonephritis%',
    '%postpartum%',
    '%postabortal%',
    '%pelvic%',
    '%arthritis%',
    '%anemia%',
    '%titers%',
    '%ureaplasma%',
    '%meleagridis%',
    '%mastitis%',
    '%unspecified%',
    '%elsewhere%',
    '%balanitis%',
    '%infertility%',
    '%conjunctivitis%',
    '%genital%',
    '%multiforme%',
    '%swine%',
    '%agglutinin%',
    '%gallisepticum%',
    '%synovitis%',
    '%not isolated%',
    '%not detected%',
    '%equine%'   
]

# Querying the concept codes
mpa_codes_q = construct_concept_query(version, 'mycoplasma', mpa_exclusion)


# In[ ]:


mpa_codes = polars_gbq(mpa_codes_q)


# In[ ]:


mpa_codes


# In[ ]:


# fetch_co(version, mpa_codes_q, return_count=True)
# 279 v7
# 444 v8


# In[ ]:


# fetch_o(version, mpa_codes_q, return_count=True)
#0 v7, v8


# In[ ]:


mpa_code_df = fetch_co(version, mpa_codes_q, return_count=False)


# In[ ]:


mpa_code_df['person_id'].n_unique()
# 171 v7
# 271 v8


# In[ ]:


mpa_code_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(10)


# ### M. pneumoniae Tests

# In[ ]:


# Search panels for Mycoplasma tests
mpa_tests_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    c.vocabulary_id = 'LOINC'
    AND LOWER(c.concept_name) LIKE '%mycoplasma%'
    AND (
        LOWER(c.concept_name) LIKE '%[presence]%'
        OR LOWER(c.concept_name) LIKE '%[identifier]%'
        OR LOWER(c.concept_name) LIKE '%[cycle threshold]%'
        OR LOWER(c.concept_name) LIKE '%[susceptibility]%'
        OR LOWER(c.concept_name) LIKE '%viral load%'
        OR LOWER(c.concept_name) LIKE '%culture%'
        OR c.concept_name LIKE '%NAA%'
        OR c.concept_name LIKE '%Measurement of%'
    )
    AND NOT EXISTS (
        SELECT 1
        FROM UNNEST(ARRAY[
            '%porcine%', 
            '%canine%',
            '%bovine%',
            '%cornea%',
            '%chlamydophila%',
            '%genital%',
            '%agalactiae%',
            '%bovis%',
            '%capricolum%',
            '%fermentans%',
            '%flocculare%',
            '%gallisepticum%',
            '%genitalium%',
            '%hyopneumoniae%',
            '%hominis%',
            '%hyorhinis%',
            '%hyosynoviae%',
            '%iowae%',
            '%meleagridis%',
            '%mycoides%',
            '%penetrans%',
            '%body%',
            '%exudate%',
            '%milk%',
            '%poultry%',
            '%ureaplasma%',
            '%urine%',
            '%suis%',
            '%synoviae%',
            '%urine%',
            '%igg%',
            '%igm%',
            '%cerebral%',
            '% ab %'
        ]) AS excluded_pattern
        WHERE LOWER(c.concept_name) LIKE excluded_pattern
    )
ORDER BY c.concept_name
"""

# polars_gbq(mpa_tests_q)


# In[ ]:


m_mpa_tests = f"""
        SELECT
            m.person_id,
            m.measurement_concept_id as concept_id,
            c.concept_name as concept_name,
            c.concept_code as code,
            c.vocabulary_id as vocab,
            m.measurement_datetime as datetime,
            c2.concept_name AS result,
            CAST(m.visit_occurrence_id AS STRING) as visit_id,
            v.visit_start_datetime as visit_start,
            c1.concept_name AS visit_name
        FROM
            {version}.measurement m
        INNER JOIN
            {version}.concept c ON m.measurement_concept_id = c.concept_id
        LEFT JOIN
            {version}.concept c2 ON m.value_as_concept_id = c2.concept_id
        LEFT JOIN
            {version}.visit_occurrence v ON m.visit_occurrence_id = v.visit_occurrence_id AND m.person_id = v.person_id
        LEFT JOIN
            {version}.concept c1 ON v.visit_concept_id = c1.concept_id      
        WHERE
            m.measurement_concept_id IN (SELECT DISTINCT concept_id FROM ({mpa_tests_q}))
        ORDER BY 
            person_id, measurement_datetime
        """

mpa_test_df = polars_gbq(m_mpa_tests)
len(mpa_test_df)
# 13603 v7
# 21819 v8


# In[ ]:


mpa_test_df['person_id'].n_unique()
# 5909 v7
# 10633 v8


# In[ ]:


mpa_test_df.group_by(['concept_id', 'concept_name', 'result']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(10)


# ## Chlamydia pneumoniae

# ### C. pneumoniae Diagnoses

# In[ ]:


# C. pneumoniae codes
cpna_codes_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.vocabulary_id,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    c.vocabulary_id LIKE '%ICD%'
    AND c.domain_id = 'Condition'
    AND LOWER(c.concept_name) LIKE '%chlamydia%'
    AND NOT EXISTS (
        SELECT 1
        FROM UNNEST(ARRAY[
            '%vaccin%', 
            '%non-respiratory%',
            '%otitis%',
            '%carditis%',
            '%bovine%',
            '%myelitis%',
            '%gastrointestinal%',
            '%gastroenteritis%',
            '%central nervous system%',
            '%cns disorder%',
            '%encephalopathy%',
            '%encephalitis%',
            '%not isolated%',
            '%not detected%',
            '%sexually%',
            '%chlamydiae%',
            '%other chlamydial%',
            '%unspecified%',
            '%genitourinary%',
            '%conjunctivitis%',
            '%and chlamydial%',
            '%or chlamydial%',
            '%venereal%',
            '%trachomatis%',
            '%conjunctiva%',
            '%sequela%',
            '%arthritis%',
            '%hepatitis%',
            '%other specified%',
            '%dacryocystitis%',
            '%pelvic%',
            '%cystitis%',
            '%congenital%',
            '%vaginitis%',
            '%urethritis%',
            '%salpingitis%',
            '%prostatitis%',
            '%peritonitis%',
            '%lymphogranuloma%',
            '%epididy%',
            '%anus%',
            '%derma%',
            '%neonatal%',
            '%cervicitis%',
            '%bartholinitis%',
            '%balanitis%',
            '%chlamydia test%',
            '%psittaci%',
            '%psitacci%',
            '%elisa%',
            '%pcr%',
            '%genital%',
            '%aids%',
            '%equine%'
        ]) AS excluded_pattern
        WHERE LOWER(c.concept_name) LIKE excluded_pattern
    )
    AND concept_id != 438066
ORDER BY c.concept_name DESC
"""

cpna_codes = polars_gbq(f"""({cpna_codes_q}) LIMIT 150""")


# In[ ]:


cpna_codes.head()


# In[ ]:


# fetch_co(version, cpna_codes_q, return_count=True)
# 178 v7
# 231 v8


# In[ ]:


# fetch_o(version, cpna_codes_q, return_count=True)
#0 v7 v8


# In[ ]:


cpna_code_df = fetch_co(version, cpna_codes_q, return_count=False)


# In[ ]:


cpna_code_df['person_id'].n_unique()
# 101 v7
# 134 v8


# In[ ]:


cpna_code_df.group_by(['concept_id', 'concept_name']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(10)


# ### C. pneumoniae Tests

# In[ ]:


# Search panels for C. pneumoniae tests
# after looking at panel descendants (below), 36303621 and 21493344 [Chlamydophila pneumoniae DNA...] were missed. These are now added manually in the WHERE clause.

cpna_tests_q = f"""
SELECT 
    DISTINCT c.concept_id,
    c.concept_code,
    c.concept_name
FROM 
    `{version}.concept` c
WHERE
    c.vocabulary_id = 'LOINC'
    AND (LOWER(c.concept_name) LIKE '%chlamydia%' OR LOWER(c.concept_name) LIKE '%chlamydophila pneumoniae%')
    AND (
        LOWER(c.concept_name) LIKE '%[presence]%'
        OR LOWER(c.concept_name) LIKE '%[identifier]%'
        OR LOWER(c.concept_name) LIKE '%[cycle threshold]%'
        OR LOWER(c.concept_name) LIKE '%[susceptibility]%'
        OR LOWER(c.concept_name) LIKE '%viral load%'
        OR LOWER(c.concept_name) LIKE '%culture%'
        OR c.concept_name LIKE '%NAA%'
        OR c.concept_name LIKE '%Measurement of%'
    )
    AND NOT EXISTS (
        SELECT 1
        FROM UNNEST(ARRAY[
            '%porcine%', 
            '%canine%',
            '%bovine%',
            '%cornea%',
            '%trachomatis%',
            '%psittaci%',
            '%mycoplasma%',
            '%igg%',
            '%igm%',
            '%cerebral%',
            '%anal%',
            '%body%',
            '%cervix%',
            '%conjunctival%',
            '%genital%',
            '%peritoneal%',
            '%stool%',
            '%urethra%',
            '%urine%',
            '%vaginal%',
            '%serum%',
            '%penis%',
            '% ab %'
        ]) AS excluded_pattern
        WHERE LOWER(c.concept_name) LIKE excluded_pattern
    )
ORDER BY c.concept_name
"""

# polars_gbq(cpna_tests_q)


# In[ ]:


c_pna_tests = f"""
        SELECT
            m.person_id,
            m.measurement_concept_id as concept_id,
            c.concept_name as concept_name,
            c.concept_code as code,
            c.vocabulary_id as vocab,
            m.measurement_datetime as datetime,
            c2.concept_name AS result,
            CAST(m.visit_occurrence_id AS STRING) as visit_id,
            v.visit_start_datetime as visit_start,
            c1.concept_name AS visit_name
        FROM
            {version}.measurement m
        INNER JOIN
            {version}.concept c ON m.measurement_concept_id = c.concept_id
        LEFT JOIN
            {version}.concept c2 ON m.value_as_concept_id = c2.concept_id
        LEFT JOIN
            {version}.visit_occurrence v ON m.visit_occurrence_id = v.visit_occurrence_id AND m.person_id = v.person_id
        LEFT JOIN
            {version}.concept c1 ON v.visit_concept_id = c1.concept_id      
        WHERE
            m.measurement_concept_id IN (SELECT DISTINCT concept_id FROM ({cpna_tests_q}))
        ORDER BY 
            person_id, measurement_datetime
        """

cpna_test_df = polars_gbq(c_pna_tests)
len(cpna_test_df)
# 15093 v7
# 30289 v8


# In[ ]:


cpna_test_df['person_id'].n_unique()
# 6696 v7
# 14510 v8


# In[ ]:


cpna_test_df.group_by(['concept_id', 'concept_name', 'result']).agg(pl.len().alias('count')).sort(['count'], descending=True).head(10)


# ## Confirm Complete Respiratory Virus Panel Results
# This uses test result concept_ids to identify panels that contain these results. Identified panel descendents are then queried in the measurement table to ensure that no relevated panel test results were missed.  

# In[ ]:


#Combine test dataframes by concept_id
combined_tests_df = pl.concat([para_test_df[['concept_id']].unique(), 
                               hmpv_test_df[['concept_id']].unique(), 
                               rv_test_df[['concept_id']].unique(), 
                               adv_test_df[['concept_id']].unique(), 
                               rsv_test_df[['concept_id']].unique(), 
                               hcov_test_df[['concept_id']].unique(), 
                               flu_test_df[['concept_id']].unique(), 
                               covid_test_df[['concept_id']].unique()])


# In[ ]:


combined_test_list = tuple(combined_tests_df['concept_id'])
len(combined_test_list)
# 239 v7
# 274 v8
# 275 v8 updated N3C variables


# In[ ]:


# identify which panels contain these test concept_ids + panels 21491676 and 1175736 (Bordatella results)
# also, exclude 37021437, which is a mostly bacterial LRTI panel

panels_q = f"""
        SELECT
            DISTINCT 
            concept_id_2,
            c1.concept_name
        FROM
            {version}.concept_relationship cr
        INNER JOIN
            {version}.concept c1 ON cr.concept_id_2 = c1.concept_id
        WHERE
            (concept_id_1 IN {combined_test_list}
            OR concept_id_2 = 21491676
            OR concept_id_2 = 1175736)
            AND relationship_id = 'Contained in panel'
            AND concept_id_2 != 37021437
        ORDER BY concept_id_2
        """

panels_df = polars_gbq(panels_q)


# In[ ]:


panels_df.head()


# In[ ]:


panels_str = tuple(panels_df['concept_id_2'])
len(panels_str)
# 34 v7
# 41 v8


# In[ ]:


# select all descendants from included panels
panel_descendants = f"""
        SELECT
            concept_id_1,
            c.concept_name,
            relationship_id,
            concept_id_2,
            c1.concept_name
        FROM
            {version}.concept_relationship cr
        INNER JOIN
            {version}.concept c ON cr.concept_id_1 = c.concept_id
        INNER JOIN
            {version}.concept c1 ON cr.concept_id_2 = c1.concept_id
        WHERE
            concept_id_1 IN {panels_str}
            AND relationship_id = 'Panel contains'
        ORDER BY concept_id_1
        """

panel_descendants_df = polars_gbq(panel_descendants)


# In[ ]:


len(panel_descendants_df)
# 359 rows
# 393 rows


# In[ ]:


#Combine test dataframes for pertussis, C. pneumoniae, and M. pneumoniae
combined_added_tests_df = pl.concat([pert_test_df[['concept_id']].unique(), 
                                     mpa_test_df[['concept_id']].unique(), 
                                     cpna_test_df[['concept_id']].unique()])


# In[ ]:


#filtering out original viruses (para, hmpv, rv, boca, adv, rsv, cov, flu, covid)
filtered_panel_descendants_df = panel_descendants_df.filter(
    ~pl.col("concept_id_2").is_in(combined_tests_df["concept_id"])
)

filtered_panel_descendants_df.head() 
# 166 total 359 rows v7
# 173 total 393 rows v8


# In[ ]:


#filtering out added tests (pertussis, M. pna, C. pna) 
filtered_panel_descendants_df = filtered_panel_descendants_df.filter(
    ~pl.col("concept_id_2").is_in(combined_added_tests_df['concept_id'])
)

len(filtered_panel_descendants_df)
#13 more rows removed (153 total) v7
# 159 v8


# In[ ]:


filtered_panel_descendants_list = tuple(filtered_panel_descendants_df['concept_id_2'])
len(filtered_panel_descendants_list)
# 153 v7
# 159 v8


# In[ ]:


# Which of these tests have measurements?
filtered_panel_tests_q = f"""
        SELECT
            measurement_concept_id,
            c.concept_name,
            c.vocabulary_id,
            c.concept_code,
            COUNT(*) as count
        FROM
            {version}.measurement
        INNER JOIN
            {version}.concept c ON measurement_concept_id = c.concept_id
        WHERE
            measurement_concept_id IN {filtered_panel_descendants_list}
        GROUP BY
            measurement_concept_id,
            c.concept_name,
            c.vocabulary_id,
            c.concept_code
        ORDER BY
            count DESC
        """

polars_gbq(filtered_panel_tests_q)


# V7: S. pyogenes can be ignored as the only typical bacterial pathogen. Bordetella holmesii is an underrecognized pathogen.<br>
# V8: Specimen source identified only significant result

# # Review Null Results

# In[ ]:


# List of dataframe names
code_df_names = ['para_code_df', 'hmpv_code_df', 'rv_code_df',
                 'adv_code_df', 'rsv_code_df', 'hcov_code_df', 'flu_code_df', 
                 'covid_code_df', 'pert_code_df', 'mpa_code_df', 'cpna_code_df']

# List of dataframes
code_dfs = [para_code_df, hmpv_code_df, rv_code_df,
            adv_code_df, rsv_code_df, hcov_code_df, flu_code_df, 
            covid_code_df, pert_code_df, mpa_code_df, cpna_code_df]

# Create a list of dictionaries, each containing NaN counts for a dataframe
code_nan_counts = [
    {
        "DataFrame": df_name,
        **{col: count.item() if isinstance(count, pl.Series) else count 
           for col, count in df.null_count().to_dict().items()}
    } 
    for df_name, df in zip(code_df_names, code_dfs)
]
code_row_counts = {df_name: df.height for df_name, df in zip(code_df_names, code_dfs)}

# Convert the list of dictionaries to a DataFrame
code_nan_counts_df = pl.DataFrame(code_nan_counts)

# Add the total row count column
code_nan_counts_df = code_nan_counts_df.with_columns(
    pl.Series("Total Rows", [code_row_counts[name] for name in code_nan_counts_df["DataFrame"]])
)

code_nan_counts_df


# In[ ]:


# List of dataframe names
test_df_names = ['para_test_df', 'hmpv_test_df', 'rv_test_df',
            'adv_test_df', 'rsv_test_df', 'hcov_test_df', 'flu_test_df', 
            'covid_test_df', 'pert_test_df', 'mpa_test_df', 'cpna_test_df']

# List of dataframes
test_dfs = [para_test_df, hmpv_test_df, rv_test_df,
            adv_test_df, rsv_test_df, hcov_test_df, flu_test_df, 
            covid_test_df, pert_test_df, mpa_test_df, cpna_test_df]

# Create a list of dictionaries, each containing NaN counts for a dataframe
nan_counts = [
    {
        "DataFrame": df_name,
        **{col: count.item() if isinstance(count, pl.Series) else count 
           for col, count in df.null_count().to_dict().items()}
    } 
    for df_name, df in zip(test_df_names, test_dfs)
]
row_counts = {df_name: df.height for df_name, df in zip(test_df_names, test_dfs)}

# Convert the list of dictionaries to a DataFrame
nan_counts_df = pl.DataFrame(nan_counts)

# Add the total row count column
nan_counts_df = nan_counts_df.with_columns(
    pl.Series("Total Rows", [row_counts[name] for name in nan_counts_df["DataFrame"]])
)

nan_counts_df


# In[ ]:


# List of dataframe names
drug_df_names = ['flu_drug_df', 'covid_drug_df']

# List of dataframes
drug_dfs = [flu_drug_df, covid_drug_df]

# Create a list of dictionaries, each containing NaN counts for a dataframe
drug_nan_counts = [
    {
        "DataFrame": df_name,
        **{col: count.item() if isinstance(count, pl.Series) else count 
           for col, count in df.null_count().to_dict().items()}
    } 
    for df_name, df in zip(drug_df_names, drug_dfs)
]
drug_row_counts = {df_name: df.height for df_name, df in zip(drug_df_names, drug_dfs)}

# Convert the list of dictionaries to a DataFrame
drug_nan_counts_df = pl.DataFrame(drug_nan_counts)

# Add the total row count column
drug_nan_counts_df = drug_nan_counts_df.with_columns(
    pl.Series("Total Rows", [drug_row_counts[name] for name in drug_nan_counts_df["DataFrame"]])
)

drug_nan_counts_df


# 'visit_occurrence_id' notebook demonstrates a value shift error for visit_occurrence_id when importing into Pandas. Querying visit_occurrence using visit_occurrence_id from code or test entries where there is an ID but not a name confirms the join results (i.e. visit_name is missing from concept table).

# # Recategorize Test Results

# ## Map 'visit_name'

# In[ ]:


# Concatenate 'visit_name' columns from each DataFrame in the three lists and get unique values
all_visit_names_df = pl.concat(
    [df.select('visit_name') for df in drug_dfs] + 
    [df.select('visit_name') for df in code_dfs] + 
    [df.select('visit_name') for df in test_dfs]
).unique()


# In[ ]:


all_visit_names_df.head()


# In[ ]:


#Map 'visit_name' according to the specified breakdown
map_ip = ['Intensive Care', 'Inpatient Visit', 'Emergency Room and Inpatient Visit', 'Inpatient Hospital', 'Hospital', 'Inpatient Psychiatric Facility', 'Psychiatric Hospital', 'Inpatient Hospice']
map_er = ['Emergency Room Visit', 'Emergency Room - Hospital', 'Observation Room']
map_uc = ['Urgent Care Facility', ]
map_facility = ['Rehabilitation Hospital', 'Nursing Facility', 'Non-hospital institution Visit', 'Behavioral Disturbances Assisted Living Facility']
map_amb_clinic = ['Ambulatory Clinic / Center', 'Ambulatory Pain Clinic / Center', 'Ambulatory Oncology Clinic / Center', 'Ambulatory Research Clinic / Center', 'Ambulatory Dental Clinic / Center']
map_amb_proc = ['Outpatient Hospital', 'Ambulatory Surgical Center', 'Ambulatory Endoscopy Cinic / Center', 'Ambulatory Infusion Therapy Clinic / Center', 'Ambulatory Oncological Radiation Clinic / Center'] 
map_amb_rad = ['Ambulatory Radiology Clinic / Center', 'Ambulatory Magnetic Resonance Imaging (MRI) Clinic / Center', 'Ambulatory Mammography Clinic / Center']
map_amb_rehab = ['Ambulatory Rehabilitation Visit', 'Comprehensive Inpatient Rehabilitation Facility'] 
map_amulance = ['Ambulance - Land']
map_op = ['Outpatient Visit', 'Office Visit']
map_cm = ['Case Management Visit']
map_home = ['Home Visit']
map_tele = ['Telehealth']
map_pharm = ['Pharmacy visit', 'Pharmacy']
map_exam = ['Health examination'] 
map_lab = ['Laboratory Visit'] 
map_vac = ['Mass Immunization Center']
map_unk = ['Unknown Value (but present in data)']

mapping_acuity = {
    5: map_ip,
    4: map_er + map_amulance,
    3: map_uc,
    2: map_facility,
    1: map_amb_clinic + map_amb_proc + map_amb_rad + map_amb_rehab  + map_op + map_cm + map_home + map_tele + map_pharm + map_exam + map_lab + map_vac,
    0: map_unk
}

mapping_group = {
    'ip': map_ip,
    'er': map_er,
    'uc': map_uc,
    'facility': map_facility,
    'amb_clinic': map_amb_clinic,
    'amb_proc': map_amb_proc,
    'amb_rad': map_amb_rad,
    'amb_rehab': map_amb_rehab,
    'ambulance': map_amulance,
    'op': map_op,
    'cm': map_cm,
    'home': map_home,
    'tele': map_tele,
    'pharm': map_pharm,
    'exam': map_exam,
    'lab': map_lab,
    'vac': map_vac,
    'unk': map_unk
}

# Create a reverse mapping for easier lookup
reverse_mapping_acuity = {visit: acuity for acuity, visits in mapping_acuity.items() for visit in visits}
reverse_mapping_group = {visit: group for group, visits in mapping_group.items() for visit in visits}

def replace_visits(row):
    acuity = reverse_mapping_acuity.get(row, None)
    group = reverse_mapping_group.get(row, None)
    return acuity, group

def apply_mappings(df):
    # Create new columns using the mappings
    new_df = df.with_columns([
        pl.col('visit_name').map_elements(lambda x: replace_visits(x)[0], return_dtype=pl.Int64).alias('visit_acuity'),
        pl.col('visit_name').map_elements(lambda x: replace_visits(x)[1], return_dtype=pl.Utf8).alias('visit_group')
    ])
    return new_df

# Apply the mappings to each DataFrame in the lists
cleaned_visit_test_dfs = [apply_mappings(df) for df in test_dfs]
cleaned_visit_code_dfs = [apply_mappings(df) for df in code_dfs]
cleaned_visit_drug_dfs = [apply_mappings(df) for df in drug_dfs]


# ## Clean 'result'

# In[ ]:


# get unique test values
test_results_df = pl.concat([df.select('result') for df in test_dfs]).unique()


# In[ ]:


test_results_df.head()


# In[ ]:


total_null_count = sum(df['result'].null_count() for df in cleaned_visit_test_dfs)

print(f"Total number of null values across all DataFrames: {total_null_count}")
# 567 v7
# 1378 v8


# In[ ]:


# Function to replace nulls in the 'result' column with 'null' string
def replace_nulls(df):
    return df.with_columns(
        pl.col('result').fill_null('null')
    )

# Apply the function to each DataFrame in the list
cleaned_visit_test_dfs = [replace_nulls(df) for df in cleaned_visit_test_dfs]


# In[ ]:


df = cleaned_visit_test_dfs[7]
# view odd reuslts - most within COVID results
# df.filter(pl.col('result').is_in(['SARS coronavirus AS', 'Not pregnant', 'Osteoarthritis',
#                                  '31.2', '13.8', '42.9', 'Measurement of Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)']))


# In[ ]:


# To replace values in the 'result' columns of the dataframes according to the specified breakdown, 
# we will define mappings for each category (0, 1, and 2) and then apply these mappings to each dataframe.

# Define the mappings
mapping_0 = ['Not detected', 'Negative', 'Normal', 'Not detected/negative', 'No', 'Non-Reactive', 'Nonreactive', 'Influenza A virus not detected', '0.0']
mapping_1 = ['Detected', 'Abnormal', 'Positive', '+', 'High', 'Presumptive positive', 'Yes', 'Present', 'Influenza A virus', 'Influenza A virus antigen']
mapping_2 = ['No matching concept', 'Indeterminate', 'Not tested', 'Pending', 'Refused', 'Invalid', 'DNR', 'Inconclusive', 
             'Inadequate', 'Equivocal', 'N/A', 'NA', 'Quantity insufficient', 'Nasopharyngeal swab', 'Patient refused', 
             'Test not performed', 'Not given', 'Service comment', 'null', 'Null', 'Uncertain', 'Unknown', 'SARS coronavirus AS', 
             '42.9', '27.2', 'Test not done', 'Not pregnant', 'Service comment', 'Specimen from nasopharyngeal structure',
            'Reported', '31.2', '13.8', 'Nasopharyngeal swab', 'Measurement of Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)', 
             'Osteoarthritis', 'Influenza A virus antigen']

# Create a function to apply the mappings
def replace_results(row):
    if row in mapping_0:
        return 2
    elif row in mapping_1:
        return 1
    elif row in mapping_2:
        return 0
    else:
        return row  # or handle unknown values as you see fit

# Apply replace_results to each DataFrame in the list
cleaned_result_visit_test_dfs = [df.with_columns(pl.col('result').map_elements(replace_results, return_dtype=pl.Int64).cast(pl.Int64)) for df in cleaned_visit_test_dfs]


# ## Categorize 'concept_name'

# In[ ]:


# Concatenate 'result' columns from each DataFrame and get unique values
test_concept_names_df = pl.concat([df.select('concept_name') for df in test_dfs]).unique()


# In[ ]:


test_concept_names_df.write_csv('test_concept_names.csv')


# In[ ]:


total_null_count = sum(df['concept_name'].null_count() for df in test_dfs)

print(f"Total number of null values across all DataFrames: {total_null_count}")


# In[ ]:


os.system(f"gsutil cp test_concept_names_categorized_v8.csv {my_bucket}/data/")
# found at: {bucket or my_bucket}/data/test_concept_names_categorized_v8.csv


# In[ ]:


# Load the CSV file with categorizations
concept_name_categorizations_df = pl.read_csv(f'{bucket or my_bucket}/data/test_concept_names_categorized_v8.csv')


# In[ ]:


concept_name_categorizations_df.head()


# In[ ]:


updated_result_visit_test_dfs = [df.join(concept_name_categorizations_df, on='concept_name', how='left') for df in cleaned_result_visit_test_dfs]


# In[ ]:


updated_result_visit_test_dfs[6].head()


# In[ ]:


# DataFrame names
test_df_names = ['para', 'hMPV', 'RV',
                 'ADV', 'RSV', 'hCOV', 'Flu', 
                 'COVID', 'pertussis', 'M. pna', 'C. pna']

# Check the number of unique 'pathogen' values in each DataFrame
for df, name in zip(updated_result_visit_test_dfs, test_df_names):
    unique_pathogens = df['pathogen'].unique().len()
    
    if unique_pathogens != 1:
        print(f"Quality Check Failed: {name} has {unique_pathogens} unique 'pathogen' values.")
    else:
        print(f"Quality Check Passed: {name} has 1 unique 'pathogen' value.")


# ## Review Counts

# In[ ]:


# Define a mapping of old keys to new keys
key_mapping = {
 'para_code_df': 'para',
 'hmpv_code_df': 'hMPV',
 'rv_code_df': 'RV',
 'adv_code_df': 'ADV',
 'rsv_code_df': 'RSV',
 'hcov_code_df': 'hCOV',
 'flu_code_df': 'Flu',
 'covid_code_df': 'COVID',
 'pert_code_df': 'pertussis',
 'mpa_code_df': 'M. pna',
 'cpna_code_df': 'C. pna'
}

# Apply the mapping to rename keys
for old_key, new_key in key_mapping.items():
    if old_key in code_row_counts:
        code_row_counts[new_key] = code_row_counts.pop(old_key)
        
code_row_counts
#{'para': 230, 'hMPV': 598, 'RV': 687, 'ADV': 396, 'RSV': 2421, 'hCOV': 2079, 'Flu': 24661,
# 'COVID': 45723, 'pertussis': 33, 'M. pna': 254, 'C. pna': 168} #v7 
# {'para': 331, 'hMPV': 910, 'RV': 1075, 'ADV': 606, 'RSV': 3502, 'hCOV': 3124, 'Flu': 39856,
# 'COVID': 103590, 'pertussis': 396, 'M. pna': 444, 'C. pna': 231}


# In[ ]:


# DataFrame names
test_df_names = ['para', 'hMPV', 'RV',
                 'ADV', 'RSV', 'hCOV', 'Flu', 
                 'COVID', 'pertussis', 'M. pna', 'C. pna']

unique_person_id_code_counts = {}

for df_name, df in zip(test_df_names, cleaned_visit_code_dfs):
    # Count the number of unique person_id values in each DataFrame
    unique_count = df.select(pl.col("person_id")).n_unique()
    
    # Store the count in a dictionary with the DataFrame name as the key
    unique_person_id_code_counts[df_name] = unique_count

unique_person_id_code_counts
# {'para': 136, 'hMPV': 305, 'RV': 365, 'ADV': 162, 'RSV': 880, 'hCOV': 1091, 'Flu': 12372, 
# 'COVID': 16974, 'pertussis': 23, 'M. pna': 171, 'C. pna': 101}
# {'para': 197, 'hMPV': 454, 'RV': 531, 'ADV': 221, 'RSV': 1312, 'hCOV': 1694, 'Flu': 17456, 
# 'COVID': 37435, 'pertussis': 30, 'M. pna': 271, 'C. pna': 134}


# In[ ]:


# DataFrame names in the desired order
desired_order = ['COVID', 'Flu', 'para', 'hCOV', 'RSV', 'ADV', 'RV', 'hMPV', 'pertussis', 'C. pna', 'M. pna']

# Reorder data according to the desired order
total_counts = [code_row_counts[label] for label in desired_order]
unique_counts = [unique_person_id_code_counts[label] for label in desired_order]

# Position of bars on the x-axis
ind = np.arange(len(desired_order))

# Size of the plot
fig, ax = plt.subplots(figsize=(12, 8), dpi=300)

# Width of a bar
width = 0.35       

# Plotting
rects1 = ax.bar(ind - width/2, total_counts, width, label='Total Row Counts')
rects2 = ax.bar(ind + width/2, unique_counts, width, label='Unique Person ID Counts')

# Add some text for labels, title, and custom x-axis tick labels
ax.set_ylabel('Counts')
ax.set_title('ICD Code Counts by Virus and Unique Person IDs')
ax.set_xticks(ind)
ax.set_xticklabels([])
ax.legend()

# Table data
table_data = [total_counts, unique_counts]
columns = desired_order
rows = ['Total ICDs', 'Unique Person IDs with ICDs']

# Adding table below the plot
table = ax.table(cellText=table_data, rowLabels=rows, colLabels=columns, loc='bottom', cellLoc='center')

# Set the color of table cells to match the tick color
tick_color = 'gray'  # Replace with the desired color
for key, cell in table.get_celld().items():
    cell.set_edgecolor(tick_color)
    cell.set_linewidth(0.5)

    # Increase the row heights of the table
table.set_fontsize(10)  # Adjust fontsize as needed
table.scale(1, 2.8)  # Scale the rows (1 in x-direction, 1.5 in y-direction)

# Adjust layout to make room for the table
plt.subplots_adjust(left=0.2, bottom=0.5)  # Increase bottom parameter to make more room for the table
plt.legend(loc='upper right')
plt.show()


# In[ ]:


# Define a mapping of old keys to new keys
key_mapping = {
 'flu_drug_df': 'Flu',
 'covid_drug_df': 'COVID'
}

# Apply the mapping to rename keys
for old_key, new_key in key_mapping.items():
    if old_key in drug_row_counts:
        drug_row_counts[new_key] = drug_row_counts.pop(old_key)
        
drug_row_counts
# {'Flu': 11619, 'COVID': 2664} v7
# {'Flu': 15716, 'COVID': 15117} v8


# In[ ]:


# DataFrame names
drug_df_names = ['Flu', 'COVID']

unique_person_id_drug_counts = {}

for df_name, df in zip(drug_df_names, cleaned_visit_drug_dfs):
    # Count the number of unique person_id values in each DataFrame
    unique_count = df.select(pl.col("person_id")).n_unique()
    
    # Store the count in a dictionary with the DataFrame name as the key
    unique_person_id_drug_counts[df_name] = unique_count

unique_person_id_drug_counts
# {'Flu': 7460, 'COVID': 1329} v7
# {'Flu': 10136, 'COVID': 10364} v8


# In[ ]:


# DataFrame names in the desired order
desired_order = ['COVID', 'Flu']

# Reorder data according to the desired order
total_counts = [drug_row_counts[label] for label in desired_order]
unique_counts = [unique_person_id_drug_counts[label] for label in desired_order]

# Position of bars on the x-axis
ind = np.arange(len(desired_order))

# Size of the plot
fig, ax = plt.subplots(figsize=(1, 4), dpi=300)

# Width of a bar
width = 0.35       

# Plotting
rects1 = ax.bar(ind - width/2, total_counts, width, label='Total Row Counts')
rects2 = ax.bar(ind + width/2, unique_counts, width, label='Unique Person ID Counts')

# Add some text for labels, title, and custom x-axis tick labels
ax.set_ylabel('Counts')
ax.set_title('Rx Counts by Virus')
ax.set_xticks(ind)
ax.set_xticklabels([])
ax.legend()

# Table data
table_data = [total_counts, unique_counts]
columns = desired_order
rows = ['Total Treatments', 'Unique Person IDs with Treatment']

# Adding table below the plot
table = ax.table(cellText=table_data, rowLabels=rows, colLabels=columns, loc='bottom', cellLoc='center')

# Set the color of table cells to match the tick color
tick_color = 'gray'  # Replace with the desired color
for key, cell in table.get_celld().items():
    cell.set_edgecolor(tick_color)
    cell.set_linewidth(0.5)

    # Increase the row heights of the table
table.set_fontsize(10)  # Adjust fontsize as needed
table.scale(1, 2.8)  # Scale the rows (1 in x-direction, 1.5 in y-direction)

# Adjust layout to make room for the table
plt.subplots_adjust(left=0.2, bottom=0.5)  # Increase bottom parameter to make more room for the table
plt.legend(loc='upper right')
plt.show()


# In[ ]:


# DataFrame names
test_df_names = ['para', 'hMPV', 'RV',
                 'ADV', 'RSV', 'hCOV', 'Flu', 
                 'COVID', 'pertussis', 'M. pna', 'C. pna']

unique_person_id_test_counts = {}

for df_name, df in zip(test_df_names, cleaned_result_visit_test_dfs):
    # Count the number of unique person_id values in each DataFrame
    unique_count = df.select(pl.col("person_id")).n_unique()
    
    # Store the count in a dictionary with the DataFrame name as the key
    unique_person_id_test_counts[df_name] = unique_count

unique_person_id_test_counts
# {'para': 7161, 'hMPV': 6934, 'RV': 6418, 'ADV': 8675, 'RSV': 12529, 'hCOV': 6097, 'Flu': 33655, 
#  'COVID': 93903, 'pertussis': 10902, 'M. pna': 5909, 'C. pna': 6696}
# {'para': 14681, 'hMPV': 14088, 'RV': 13506, 'ADV': 17362, 'RSV': 29183, 'hCOV': 13368, 'Flu': 60805, 
#  'COVID': 138161, 'pertussis': 17520, 'M. pna': 10633, 'C. pna': 14510}


# In[ ]:


results_count_dfs = []

for i, df in enumerate(cleaned_result_visit_test_dfs):
    # Counting occurrences of each unique value in the 'result' column
    count_df = df.group_by('result').len()

    # Renaming the columns
    count_df = count_df.rename({'result': 'Result', 'len': 'Count'})

    # Optionally, add a column to identify the DataFrame
    count_df = count_df.with_columns(pl.lit(i).alias('DataFrameIndex'))

    results_count_dfs.append(count_df)

# Optional: Concatenating the results for all DataFrames, if needed
concatenated_results_count = pl.concat(results_count_dfs)


# In[ ]:


# DataFrame names in the desired order
desired_order = ['COVID', 'Flu', 'para', 'hCOV', 'RSV', 'ADV', 'RV', 'hMPV', 'pertussis', 'C. pna', 'M. pna']

# Convert the Polars DataFrame to Pandas DataFrame for plotting
concatenated_results_count_pd = concatenated_results_count.to_pandas()

# Pivot the DataFrame for plotting
pivot_df = concatenated_results_count_pd.pivot(index='DataFrameIndex', columns='Result', values='Count').fillna(0)

# Replace the index with DataFrame names
pivot_df.index = test_df_names

# Define a mapping for the column names (0, 1, 2) to more descriptive names
column_name_mapping = {
    0: 'indeterminate',
    1: 'positive',
    2: 'negative'
}

# Rename the columns using the mapping
pivot_df = pivot_df.rename(columns=column_name_mapping)

# Ensure the DataFrame index is in the desired order
pivot_df = pivot_df.reindex(desired_order)

# Create the figure for the plot
fig, ax = plt.subplots(figsize=(12, 8), dpi=300)  # Adjust the figure size to accommodate the table

# Plot the stacked bar chart
pivot_df.plot(kind='bar', stacked=True, ax=ax)

# Set labels for the y-axis
ax.set_ylabel('Count', color='k')
ax.tick_params(axis='y', labelcolor='k')
ax.set_title('Total Number of Results by Pathogen')
ax.set_xticklabels([])

# Transpose the DataFrame for the table
pivot_df_transposed = pivot_df.T

# Create a table showing the transposed numerical values below the chart
table_data = pivot_df_transposed.values
columns = pivot_df_transposed.columns
rows = pivot_df_transposed.index
table = ax.table(cellText=table_data, rowLabels=rows, colLabels=columns, loc='bottom', cellLoc='center')

# Set the color of table cells to match the tick color
tick_color = 'gray'  # Replace with the desired color
for key, cell in table.get_celld().items():
    cell.set_edgecolor(tick_color)
    cell.set_linewidth(0.5)

# Increase the row heights of the table
table.set_fontsize(10)  # Adjust fontsize as needed
table.scale(1, 2.8)  # Scale the rows (1 in x-direction, 1.5 in y-direction)

# Adjust layout to make room for the table
plt.subplots_adjust(left=0.2, bottom=0.5)  # Increase bottom parameter to make more room for the table
plt.legend(loc='upper right')
plt.show()


# In[ ]:


# Convert counts to percentages
pivot_df_percent = pivot_df.div(pivot_df.sum(axis=1), axis=0) * 100

# Define a mapping for the column names (0, 1, 2) to more descriptive names
column_name_mapping = {
    0: 'indeterminate',
    1: 'positive',
    2: 'negative'
}

# Rename the columns using the mapping
pivot_df_percent = pivot_df_percent.rename(columns=column_name_mapping)

# Replace the index with DataFrame names
pivot_df_percent.index = test_df_names

# DataFrame names in the desired order
desired_order = ['COVID', 'Flu', 'para', 'hCOV', 'RSV', 'ADV', 'RV', 'hMPV', 'pertussis', 'C. pna', 'M. pna']

# Ensure the DataFrame index is in the desired order
pivot_df_percent = pivot_df_percent.reindex(desired_order)

# Create a stacked bar chart of percentages
pivot_df_percent.plot(kind='bar', stacked=True, figsize=(12, 8))

plt.title('Stacked Bar Chart of Result Percentages by DataFrame')
plt.xlabel('DataFrame Name')
plt.ylabel('Percentage (%)')
plt.xticks(rotation=45)  # Adjust rotation for better label visibility
plt.legend(title='Result', bbox_to_anchor=(1, 1), loc='upper left')
plt.tight_layout()  # Adjust layout to fit the legend
plt.show()


# In[ ]:


# DataFrame names
test_df_names = ['para', 'hMPV', 'RV',
                 'ADV', 'RSV', 'hCOV', 'Flu', 
                 'COVID', 'pertussis', 'M. pna', 'C. pna']

# Aggregate 'result' counts by 'test_type' for each DataFrame
aggregated_data = []

for i, (df, name) in enumerate(zip(updated_result_visit_test_dfs, test_df_names)):
    # Filter for positive tests
    positive_tests = df.filter(pl.col('result') == 1)

    # Aggregate counts by 'test_type'
    aggregated_df = positive_tests.group_by('test_type').len().to_pandas()
    aggregated_df['DataFrameName'] = name

    aggregated_data.append(aggregated_df)
    
# Combine the aggregated data
combined_df = pd.concat(aggregated_data)

# Pivot for plotting
pivot_df = combined_df.pivot_table(index='DataFrameName', columns='test_type', values='len', aggfunc='sum', fill_value=0)

# DataFrame names in the desired order
desired_order = ['COVID', 'Flu', 'para', 'hCOV', 'RSV', 'ADV', 'RV', 'hMPV', 'pertussis', 'C. pna', 'M. pna']

# Ensure the DataFrame index is in the desired order
pivot_df = pivot_df.reindex(desired_order)

# Plotting
pivot_df.plot(kind='bar', stacked=True, figsize=(10, 6))

plt.title('Total Positive Tests by Pathogen and Test Type')
plt.ylabel('Count')
plt.xticks(rotation=45)
plt.legend(title='Test Type')
plt.show()


# In[ ]:


pertussis_test_df = updated_result_visit_test_dfs[8].to_pandas()

# Convert 'datetime' to a specific time zone (e.g., UTC) or make time zone-naive
pertussis_test_df['year'] = pd.to_datetime(pertussis_test_df['datetime']).dt.year

# Group by year and result, then count occurrences
grouped_data = pertussis_test_df.groupby(['year', 'result']).size().unstack(fill_value=0)

# Reset index to make 'year' a column
grouped_data.reset_index(inplace=True)

# Plotting
plt.figure(figsize=(12, 5), dpi=300)

# Width of a bar
width = 0.25

# Create an array with the position of each bar along the x-axis
positions = np.arange(len(grouped_data['year']))

for i, result_type in enumerate([0, 1, 2]):
    # Get the count for this result type
    result_count = grouped_data[result_type]

    # Plot the bars for this result type
    plt.bar(positions + i * width, result_count, width=width, label=f'Result {result_type}')

# Set the position of the x-ticks to be in the center of the group of bars
plt.xticks(positions + width, grouped_data['year'])

# Set labels and title
plt.ylabel('Count')
plt.title('Pertussis Test Results by Year')

# Add a legend
plt.legend(title='', loc='upper left', labels=['Indeterminate', 'Positive', 'Negative'])

# Show the plot
plt.show()


# # Create Cohorts

# ## Merge dfs and clean

# In[ ]:


# For cleaned_visit_code_dfs, add a column 'result' set to 3
for i, df in enumerate(cleaned_visit_code_dfs):
    cleaned_visit_code_dfs[i] = df.with_columns(pl.lit(3, pl.Int64).alias("result"))

# For cleaned_visit_drug_dfs, add a column 'result' set to 4
for i, df in enumerate(cleaned_visit_drug_dfs):
    cleaned_visit_drug_dfs[i] = df.with_columns(pl.lit(4, pl.Int64).alias("result"))


# In[ ]:


merged_dfs = []

for code_df, test_df in zip(cleaned_visit_code_dfs, updated_result_visit_test_dfs):
    # Merge code_df and test_df, keeping all rows
    merged_df = code_df.join(test_df, on=['person_id', 'concept_id', 'concept_name', 'code', 
                                          'vocab', 'datetime', 'visit_id', 'visit_start', 'visit_name',
                                          'visit_acuity', 'visit_group', 'result'], how='full', coalesce=True)

    # Sort by person_id and then datetime
    merged_df = merged_df.sort(['person_id', 'datetime'])

    # Append the sorted DataFrame to the merged_dfs list
    merged_dfs.append(merged_df)


# In[ ]:


# Merge in drug_dfs
merged_dfs[6] = merged_dfs[6].join(cleaned_visit_drug_dfs[0], on=['person_id', 'concept_id', 'concept_name', 
                                                                  'code', 'vocab', 'datetime', 'visit_id', 'visit_start', 
                                                                  'visit_name', 'visit_acuity', 'visit_group', 'result'], how='full', coalesce=True).sort(['person_id', 'datetime'])
merged_dfs[7] = merged_dfs[7].join(cleaned_visit_drug_dfs[1], on=['person_id', 'concept_id', 'concept_name', 
                                                                  'code', 'vocab', 'datetime', 'visit_id', 'visit_start', 
                                                                  'visit_name', 'visit_acuity', 'visit_group', 'result'], how='full', coalesce=True).sort(['person_id', 'datetime'])


# In[ ]:


merged_dfs[6].head()


# In[ ]:


row_counts = {df_name: df.height for df_name, df in zip(test_df_names, merged_dfs)}

# Convert the row_counts dictionary to a DataFrame
row_counts_df = pl.DataFrame({
    "DataFrame": list(row_counts.keys()),
    "Merged RowCount": list(row_counts.values())
})

row_counts_df


# In[ ]:


deduplicated_dfs = [df.unique() for df in merged_dfs]


# In[ ]:


# Calculate new row counts for deduplicated_dfs
new_row_counts = {df_name: df.height for df_name, df in zip(test_df_names, deduplicated_dfs)}

# Convert the new row counts to a Polars Series
new_row_counts_series = pl.Series("Deduplicated RowCount", list(new_row_counts.values()))

# Add this series as a new column to row_counts_df
row_counts_df = row_counts_df.with_columns(new_row_counts_series)

row_counts_df
# "para"	119104	118891
# "hMPV"	31804	31618
# "RV"	30684	30612
# "ADV"	37635	37562
# "RSV"	76919	76724
# "hCOV"	110471	110258
# "Flu"	405300	403946
# "COVID"	602101	577818
# "pertussis"	50548	50450
# "M. pna"	22263	22189
# "C. pna"	30520	30446


# In[ ]:


sorted_deduplicated_dfs = []

for df in deduplicated_dfs:
    # Create a custom order column using conditional statements
    df = df.with_columns(
        pl.when(pl.col("result") == 3).then(0)
         .when(pl.col("result") == 4).then(1)
         .when(pl.col("result") == 1).then(2)
         .when(pl.col("result") == 2).then(3)
         .when(pl.col("result") == 0).then(4)
        .otherwise(5)  # Catch-all for any unexpected values
         .alias("result_order")
    )

    # Sort by the custom order and then by other columns if needed
    df = df.sort(['person_id', 'datetime', 'result_order'])

    # Optionally drop the auxiliary column
    df = df.drop("result_order")

    sorted_deduplicated_dfs.append(df)


# ## Calculate time_zero (distinct episodes)
# Use 90 days to combine distinct episodes. The earliest positive result, ICD code, or antiviral will set time_zero. This will then be carried forward for all results 90 days after.<br><br>
# This takes a long time to process.

# In[ ]:


def calculate_time_zero_pd(grouped):
    # Initialize time_zero column with NaT values
    grouped['time_zero'] = pd.NaT

    # Convert datetime column to datetime type if it's not already
    grouped['datetime'] = pd.to_datetime(grouped['datetime']).dt.tz_localize(None)

    # Initialize the first date for the first row
    if grouped['result'].iloc[0] in [1, 3, 4]:
        grouped.loc[grouped.index[0], 'time_zero'] = grouped['datetime'].iloc[0]

    # Start from the second row
    for i in range(1, len(grouped)):
        # When time_zero has not been set, set time_zero if result is 1, 3, or 4
        if pd.isna(grouped['time_zero'].iloc[i-1]):
            if grouped['result'].iloc[i] in [1, 3, 4]:
                grouped.loc[grouped.index[i], 'time_zero'] = grouped['datetime'].iloc[i]
        else: 
            # When time_zero is set, calculate time difference (in hours)
            time_difference = (grouped['datetime'].iloc[i] - grouped['time_zero'].iloc[i-1]).total_seconds() / 3600

            # If more than 90 days apart, reset time_zero
            if time_difference > 2160:
                # Update time_zero if result is 1, 3, or 4
                if grouped['result'].iloc[i] in [1, 3, 4]:
                    grouped.loc[grouped.index[i], 'time_zero'] = grouped['datetime'].iloc[i]
                # Reset time_zero if result is 0 or 2
                else:
                    grouped.loc[grouped.index[i], 'time_zero'] = pd.NaT
            # If less than 90 days apart, carry forward time_zero
            else:
                grouped.loc[grouped.index[i], 'time_zero'] = grouped['time_zero'].iloc[i-1]
               
    return grouped        
        
# Assuming merged_dfs is a list of Polars DataFrames
transformed_dfs_pd = []

for df in sorted_deduplicated_dfs:
    # Convert the Polars DataFrame to a Pandas DataFrame
    df_pd = df.to_pandas()

    # Ensure datetime column is in datetime format
    df_pd['datetime'] = pd.to_datetime(df_pd['datetime']).dt.tz_localize(None)

    # Apply the calculate_time_zero_pd function
    grouped = df_pd.groupby('person_id')
    df_pd = pd.concat([calculate_time_zero_pd(group) for _, group in grouped])

    # Reset index
    df_pd = df_pd.reset_index(drop=True)
    
    # Append the transformed Pandas DataFrame to the list
    transformed_dfs_pd.append(df_pd)


# In[ ]:


for i, df in enumerate(transformed_dfs_pd):
    # Add a new column 'time_zero_1' that is a copy of 'time_zero'
    transformed_dfs_pd[i] = df.assign(time_zero_1=df['time_zero'])


# ## Extend time_zero to include earlier negative results
# If there are results (negative or indeterminate), include them in the episode if less than five days prior to positive result, billing code or treatment.<br><br>
# This also takes a long time to process.

# In[ ]:


def update_time_zero_based_on_future(grouped):
    
    result = grouped.copy()
    result = result.sort_values('datetime')

    for i in range(len(result)):
        if pd.isna(result['time_zero_1'].iloc[i]):
            current_datetime = result['datetime'].iloc[i]

            # Look ahead up to 5 days
            for j in range(i + 1, len(result)):
                if (result['datetime'].iloc[j] - current_datetime).days <= 5:
                    if not pd.isna(result['time_zero_1'].iloc[j]):
                        result.loc[result.index[i], 'time_zero_1'] = result['time_zero_1'].iloc[j]
                        break

    return result

# Apply the function to each DataFrame in your list of DataFrames
transformed_dfs_pd_updated = []

for df in transformed_dfs_pd:
    # Split the dataframe into groups
    groups = [group for _, group in df.groupby('person_id')]
    
    # Apply the function to each group
    updated_groups = [update_time_zero_based_on_future(group) for group in groups]
    
    # Concatenate the updated groups back into a single dataframe
    df_updated = pd.concat(updated_groups)
    
    # Sort by person_id and datetime to maintain original order
    df_updated = df_updated.sort_values(['person_id', 'datetime'])
    
    transformed_dfs_pd_updated.append(df_updated)


# In[ ]:


unique_person_ids = set()

for i, df in enumerate(transformed_dfs_pd_updated):
    # Create the boolean mask
    mask = (pd.notna(df['time_zero']) | pd.notna(df['time_zero_1'])) & (df['time_zero'] != df['time_zero_1'])
    
    # Apply the mask to the same DataFrame
    filtered_df = df[mask]

    unique_person_ids.update(filtered_df['person_id'].unique())

    print(f"Total number of entries where time_zero ≠ time_zero_1 for DataFrame {i}: {filtered_df['person_id'].nunique()}")


# In[ ]:


result_t0_not_t1 = transformed_dfs_pd_updated[6][transformed_dfs_pd_updated[6]['person_id'].isin(unique_person_ids)]


# In[ ]:


result_t0_not_t1.iloc[0:5]


# ## Calculate time_difference
# Determine difference in time between row and time_zero (90-day episode).

# In[ ]:


# DataFrame names
df_names = ['para', 'hMPV', 'RV', 'ADV', 'RSV', 'hCOV', 'Flu', 'COVID', 'pert', 'M_pna', 'C_pna']

# Assign names to DataFrames
for i, df in enumerate(transformed_dfs_pd_updated):
    df.name = df_names[i]

# Iterate over each DataFrame in the list
for i, df in enumerate(transformed_dfs_pd_updated):
    # Drop the 'time_zero' column
    df = df.drop(columns=['time_zero'])
    
    # Create 'time_difference' column as the difference in days between 'datetime' and 'time_zero_1'
    df['time_difference'] = (df['datetime'] - df['time_zero_1']).dt.days

    # Rename 'time_zero_1' column to 'time_zero'
    df = df.rename(columns={'time_zero_1': 'time_zero'})

    # Update the DataFrame in the list
    transformed_dfs_pd_updated[i] = dfmodified_dfs = []

# Loop through each DataFrame in the list, using their names
for i, df in enumerate(transformed_dfs_pd_updated):
    # Filter out rows based on conditions
    df_filtered = df.loc[(df['result'].isin([1, 3])) & (df['time_difference'].notna())].copy()
    df_filtered.loc[:, 'type'] = df_filtered['result'].apply(lambda x: 'lab' if x == 1 else 'icd')
    df_filtered.loc[:, 'df_id'] = i  # Add an identifier for each DataFrame
    modified_dfs.append(df_filtered)

concat_df = pd.concat(modified_dfs)

# Set the size of the plot
plt.figure(figsize=(12, 8), dpi=300)

# Plot the violin plot using the concatenated DataFrame
sns.violinplot(
    data=concat_df, x='time_difference', y='df_id', orient='h',
    split=True, gap=0.1, linewidth=1, inner='quart', hue='type'
)

# Set labels, title, and legend
plt.xlim(-5, 30)
plt.title('Time Difference for Result, ICD, or Treatment')

# Show the plot
plt.show()# Set the size of the plot
plt.figure(figsize=(12, 12), dpi=300)

# Plot the violin plot using the concatenated DataFrame
sns.stripplot(
    data=concat_df, x='time_difference', y='df_id', orient='h', hue='type', dodge=True, jitter=0.3, alpha=0.1
)

# Set labels, title, and legend
#plt.xlim(-5, 30)
plt.title('Time Difference for Result, ICD, or Treatment')

# Show the plot
plt.show()concat_df.head()
# An attempt to visualize changes in visit_type over time (i.e. more to less acuity).

# ## Save
# Save `transformed_dfs_pd_updated`

# In[ ]:


# for df in transformed_dfs_pd_updated:
#     # Convert to UTC if not already in UTC
#     df['datetime'] = df['datetime'].dt.tz_convert('UTC')
#     df['visit_start'] = df['visit_start'].dt.tz_convert('UTC')
#     df['time_zero'] = df['time_zero'].dt.tz_convert('UTC')


# In[ ]:


# get the bucket name
my_bucket = os.getenv('WORKSPACE_BUCKET')


# In[ ]:


print(subprocess.check_output(f"gsutil ls -r {my_bucket}/data/cohorts/", shell=True).decode('utf-8'))


# In[ ]:


# Removing all files within /data/cohorts/
subprocess.run(f"gsutil -m rm {my_bucket}/data/cohorts/icd_only/*", shell=True)


# In[ ]:


for df in transformed_dfs_pd_updated:
    # Construct the destination filename using the DataFrame's name
    destination_filename = df.name + ".csv"  # Assuming DataFrame has a 'name' attribute

    # save dataframe in a csv file in the same workspace as the notebook
    df.to_csv(f"{my_bucket}/data/cohorts/icd_only/{destination_filename}", index=False)

