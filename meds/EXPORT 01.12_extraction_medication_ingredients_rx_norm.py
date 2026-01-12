#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import polars as pl
from google.cloud import bigquery
import re


# In[ ]:


# Load environment variables
my_bucket = os.getenv('WORKSPACE_BUCKET')
CDR = os.getenv('WORKSPACE_CDR')
BILLING_PROJECT_ID = os.getenv('GOOGLE_PROJECT')


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


# In[ ]:


####this query is the best so far if dont need atc. works well to map drugs to rx norm
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
        cr.concept_id AS cb_criteria_root_ancestor_id,  -- Root ancestor concept ID
        cr_concept.concept_name AS cb_criteria_root_ancestor_name  -- Root ancestor concept name (ingredient name)
    FROM {CDR}.drug_exposure AS d_exposure
    LEFT JOIN {CDR}.concept AS d_standard_concept 
        ON d_exposure.drug_concept_id = d_standard_concept.concept_id
    LEFT JOIN {CDR}.concept AS d_route  -- Fix: Join to get route name
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
        WHERE cr.full_text LIKE '%_rank1]%'  -- Find only root-level (rank 1) ancestors
    ) cr 
        ON ca.ancestor_id = cr.concept_id
    LEFT JOIN {CDR}.concept AS cr_concept  -- Join to get the ingredient name
        ON cr.concept_id = cr_concept.concept_id;
"""
medication_df = polars_gbq(is_descendants_q)


# In[ ]:


# # Get all medications - query that combines getting both rxnorm and atc top level concepts
# is_descendants_q = f"""
#     WITH atc AS (
#         SELECT DISTINCT 
#             concept_id,
#             concept_code AS atc_concept_code,
#             concept_name AS atc_concept_name
#         FROM
#             {CDR}.concept
#         WHERE
#             vocabulary_id = 'ATC'
#             AND LENGTH(concept_code) = 1
#     ),
#     ranked_atc_concept_descendant AS (
#         SELECT DISTINCT 
#             ca.descendant_concept_id AS concept_id,
#             a.atc_concept_code AS ancestor_atc_code,
#             a.atc_concept_name AS ancestor_atc_name,
#             ROW_NUMBER() OVER (
#                 PARTITION BY ca.descendant_concept_id 
#                 ORDER BY ca.min_levels_of_separation ASC
#             ) AS rank
#         FROM 
#             {CDR}.concept_ancestor AS ca
#         JOIN 
#             atc AS a 
#         ON 
#             ca.ancestor_concept_id = a.concept_id
#     ),
#     filtered_data AS (
#         SELECT DISTINCT
#             de.person_id,
#             c1.concept_name AS concept_name,
#             c2.concept_name AS route,
#             a.ancestor_atc_code,
#             a.ancestor_atc_name,
#             de.drug_exposure_start_date,
#             de.drug_exposure_end_date,
#             de.quantity,
#             de.days_supply,
#             de.refills,
#             de.verbatim_end_date,
#             d_visit.concept_name AS visit_occurrence_concept_name,
#             cr.concept_id AS cb_criteria_root_ancestor_id,  -- Root ancestor concept ID
#             cr_concept.concept_name AS cb_criteria_root_ancestor_name  -- Root ancestor concept name (ingredient name)
#         FROM
#             {CDR}.drug_exposure AS de
#         LEFT JOIN
#             {CDR}.visit_occurrence AS v 
#             ON de.visit_occurrence_id = v.visit_occurrence_id
#         LEFT JOIN
#             {CDR}.concept AS d_visit 
#             ON v.visit_concept_id = d_visit.concept_id
#         LEFT JOIN
#             {CDR}.concept AS c1 
#             ON de.drug_concept_id = c1.concept_id
#         LEFT JOIN
#             {CDR}.concept AS c2 
#             ON de.route_concept_id = c2.concept_id
#         LEFT JOIN {CDR}.cb_criteria_ancestor AS cra 
#             ON de.drug_concept_id = cra.descendant_id   
#         LEFT JOIN (
#             SELECT DISTINCT cr.concept_id, CAST(cr.id AS STRING) AS id
#             FROM {CDR}.cb_criteria AS cr
#             WHERE cr.full_text LIKE '%_rank1]%'  -- Find only root-level (rank 1) ancestors
#         ) cr 
#             ON cra.ancestor_id = cr.concept_id
#         LEFT JOIN {CDR}.concept AS cr_concept  -- Join to get the ingredient name
#             ON cr.concept_id = cr_concept.concept_id
#         LEFT JOIN
#             ranked_atc_concept_descendant AS a 
#             ON de.drug_concept_id = a.concept_id
#         WHERE
#             a.rank = 1 -- Include only the top-ranked ATC ancestor for each drug concept
#     )
#     SELECT 
#         fd.person_id,
#         fd.concept_name,
#         fd.route,
#         fd.ancestor_atc_code,
#         fd.ancestor_atc_name,
#         fd.drug_exposure_start_datetime,
#         fd.drug_exposure_end_datetime,
#         fd.quantity,
#         fd.days_supply,
#         fd.refills,
#         fd.verbatim_end_date,
#         fd.visit_occurrence_concept_name,
#         fd.cb_criteria_root_ancestor_id,
#         fd.cb_criteria_root_ancestor_name
#     FROM filtered_data fd;
# """

# medication_df_atc = polars_gbq(is_descendants_q)


# In[ ]:


medication_df


# In[ ]:


# Remove these words altogether first (bc route is not systemic)
words_remove_rows_insensitive = [
    "topical", "spray", "lotion", "cream", "inhaler", "rectal", "suppository",
    "inhalation", "sponge", "patch", "metered", "ophthalmic", "ointment", 
    "shampoo", "gel", "pad", "soap", "mucosal", "vaginal", "insert", "vaccine", 
    "film", "meningococcal B", "virus", "hpv"
]

# Create regex pattern to match any of these words
remove_pattern = "|".join(map(re.escape, words_remove_rows_insensitive))  # Escape special characters

# Filter the DataFrame to exclude rows containing these words (case insensitive)
medication_df_filtered = medication_df.filter(
    ~pl.col("standard_concept_name").str.to_lowercase().str.contains(remove_pattern)
)


# In[ ]:


print(medication_df_filtered)


# In[ ]:


medication_df_filtered["route_concept_name"].value_counts().sort("count", descending=True).head(50)


# In[ ]:


###Replace route with "oral" etc (for cases where route is in name but is either miscoded or missing/NA)
oral_words = ["oral", "capsule", "softgel", "capsule", "tablet", "chewable"]

# regex pattern for case-insensitive matching
oral_pattern = "|".join(map(re.escape, oral_words))

# Update "route" column where "concept_name" contains any oral-related words
medication_df_updated = medication_df_filtered.with_columns(
    pl.when(pl.col("standard_concept_name").str.to_lowercase().str.contains(oral_pattern))
    .then(pl.lit("Oral route"))
    .otherwise(pl.col("route_concept_name"))
    .alias("route_concept_name")  # Ensure we modify the existing "route" column
)


# In[ ]:


medication_df_updated["route_concept_name"].value_counts().sort("count", descending=True).head(50)


# In[ ]:


# Keep rows w any of these words in the route name
route_keep_rows_insensitive = [
    "sublingual", "subcutaneous", "oral", "no matching concept", "NA", "Null", "intravenous", "intramuscular"
]

# Create regex pattern to match any of these words
keep_pattern = "|".join(map(re.escape, route_keep_rows_insensitive)) 

# route column contains one of the specified words (case insensitive) or null
medication_df_oral = medication_df_updated.filter(
    pl.col("route_concept_name").is_null() |
    pl.col("route_concept_name").str.to_lowercase().str.contains(keep_pattern)
)


# In[ ]:


medication_df_oral


# In[ ]:


medication_df_oral["route_concept_name"].value_counts().sort("count", descending=True)


# In[ ]:


medication_df_oral["visit_occurrence_concept_name"].value_counts().sort("count", descending=True)


# In[ ]:


# Remove rows w inpatient setting
words_remove_rows_insensitive = ["emergency", "inpatient", "intensive care"]
# words_remove_rows_insensitive_exact_match = "hospital"

# Create regex pattern for partial matches (contains)
remove_pattern = "|".join(map(re.escape, words_remove_rows_insensitive))

# remove rows where:
# - standard_concept_name contains word from words_remove_rows_insensitive (case insensitive)
# - OR standard_concept_name = exact match to words_remove_rows_insensitive_exact_match (case insensitive)
# keep rows where standard_concept_name is NULL
medication_df_oral_outpatient = medication_df_oral.filter(
    pl.col("visit_occurrence_concept_name").is_null() |  # Keep NULL values
    ~pl.col("visit_occurrence_concept_name").str.to_lowercase().str.contains(remove_pattern) #&  # Exclude partial matches
#     ~pl.col("visit_occurrence_concept_name").str.to_lowercase().eq(words_remove_rows_insensitive_exact_match)  # Exclude exact match
)


# In[ ]:


medication_df_oral_outpatient


# In[ ]:


medication_df_oral_outpatient["visit_occurrence_concept_name"].value_counts().sort("count", descending=True).slice(0, 10)


# In[ ]:


null_count = medication_df_oral_outpatient.filter(medication_df_oral_outpatient["cb_criteria_root_ancestor_name"].is_not_null()).shape[0]
print(f"Number of rows not NA for Rx norm mapping: {null_count}")


# In[ ]:


unique_concept_names = medication_df_oral_outpatient.group_by(
    ["standard_concept_name", "cb_criteria_root_ancestor_name", "route_concept_name"]
).agg(pl.len().alias("count")).sort("count", descending=True)

unique_concept_names.write_csv("meds_extract_grouped.csv")


# In[ ]:


unique_concept_names = medication_df_oral_outpatient.group_by(["cb_criteria_root_ancestor_name"]).agg(
    pl.col("person_id").n_unique().alias("count")
).sort("count", descending=True)

unique_concept_names.write_csv("meds_extract.csv")


# In[ ]:


# unique_ingredients_atc = medication_df_oral_outpatient.group_by(["ancestor_atc_code", "cb_criteria_root_ancestor_name"]).agg(
#     pl.col("person_id").n_unique().alias("count")
# ).sort("count", descending=True)

# unique_ingredients_atc.write_csv("ingredient_atc.csv")


# In[ ]:


# medication_df_atc["person_id"].n_unique()


# In[ ]:


medication_df_oral_outpatient


# In[ ]:


### Add dosage information to separate column for meds with strength in name
# Extract the numeric strength value
medication_df_doses = medication_df_oral_outpatient.with_columns(
    pl.col("standard_concept_name")
    .str.extract(r'(\d+(\.\d+)?)\s*MG', 1)  # Extracts numbers before "MG"
    .cast(pl.Float64)  # Convert to numeric
    .alias("med_strength_mg")
)


# In[ ]:


medication_df_doses


# In[ ]:


# medication_df_atc['ancestor_atc_code'].value_counts()


# In[ ]:


# medication_df_atc['ancestor_atc_name'].value_counts().sort("count", descending=True)


# ## Extract drug exposure end dates

# In[ ]:


### Apply start date as end date if NA

#in this section we calculate the end date to be the latest date for each medication record based on the 
#existing data. i.e. if end date that is coded is missing or is earlier than the verbatim end date, or that
#the quantity prescribed would indicate, then we extrapolate the end date


# In[ ]:


# Replace null values in 'drug_exposure_end_date' with corresponding values from 'drug_exposure_start_date'
medication_df_oral_outpatient = medication_df_oral_outpatient.with_columns(
    pl.when(pl.col("drug_exposure_end_date").is_null())
    .then(pl.col("drug_exposure_start_date"))
    .otherwise(pl.col("drug_exposure_end_date"))
    .alias("drug_exposure_end_date")
)


# In[ ]:


#Use verbatim end date if later than 'drug exposure end date'
# Update 'drug_exposure_end_date' where 'verbatim_end_date' is later and not null
medication_df = medication_df.with_columns(
    pl.when(
        (pl.col("drug_exposure_end_date") < pl.col("verbatim_end_date")) & 
        (pl.col("verbatim_end_date").is_not_null())
    )
    .then(pl.col("verbatim_end_date"))
    .otherwise(pl.col("drug_exposure_end_date"))
    .alias("drug_exposure_end_date")
)


# In[ ]:


# Use days supply if that + drug exposure start date is later than coded 'drug exposure end date' in the chart
# Compute the new end date by adding days_supply as an offset
new_end_date = pl.col("drug_exposure_start_date").dt.offset_by(pl.format("{}d", pl.col("days_supply")))

# Update 'drug_exposure_end_date' where the computed new end date is later
medication_df = medication_df.with_columns(
    pl.when(
        (pl.col("drug_exposure_end_date") < new_end_date) & 
        (pl.col("days_supply").is_not_null())
    )
    .then(new_end_date)
    .otherwise(pl.col("drug_exposure_end_date"))
    .alias("drug_exposure_end_date")
)


# In[ ]:


#if end date is earlier than start date + refills * quantity, then replace the end date w start date + 
#quantity*refills

# compute the number of days to add
days_to_add = (pl.col("refills").cast(pl.Float64) * pl.col("quantity").cast(pl.Float64)).cast(pl.Int64)

# Calculate new potential end date
new_end_date = pl.col("drug_exposure_start_date").dt.offset_by(pl.format("{}d", days_to_add))

# Update 'drug_exposure_end_date' where the new calculated end date is later
medication_df = medication_df.with_columns(
    pl.when(
        (pl.col("drug_exposure_end_date") < new_end_date) & 
        pl.col("refills").is_not_null() & 
        pl.col("quantity").is_not_null()
    )
    .then(new_end_date)
    .otherwise(pl.col("drug_exposure_end_date"))
    .alias("drug_exposure_end_date")
)


# In[ ]:


#### Using quantity if refills is not present

# Convert quantity to Int64 and add it as days to the start date
new_end_date_quantity = pl.col("drug_exposure_start_date").dt.offset_by(pl.col("quantity").cast(pl.Int64).alias("quantity_days"))

# Update 'drug_exposure_end_date' where refills is null and new end date based on quantity is later
medication_df = medication_df.with_columns(
    pl.when(
        (pl.col("drug_exposure_end_date") < new_end_date_quantity) & 
        pl.col("refills").is_null() & 
        pl.col("quantity").is_not_null()
    )
    .then(new_end_date_quantity)
    .otherwise(pl.col("drug_exposure_end_date"))
    .alias("drug_exposure_end_date")
)


# In[ ]:




