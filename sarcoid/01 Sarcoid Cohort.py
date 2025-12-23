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
import matplotlib, matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime, date
from typing import Dict, List, Optional, Union
from jinja2 import Template, Environment


# In[ ]:


from omop_unifier import Explorer, Mapper, Unifier


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


# # Class Definition

# In[ ]:


class AouQueries:
    def __init__(self, version=None, bucket=None):
        """
        Specialized query builder for common All of Us research patterns.
        
        """
        self.version = version or os.getenv('WORKSPACE_CDR')
        self.bucket = bucket or os.getenv('WORKSPACE_BUCKET')
        self.client = bigquery.Client()
        self.query_text = None

        # Validate required environment
        if not self.version:
            raise ValueError("CDR version not found. Check WORKSPACE_CDR environment variable if none provided.")
        if not self.bucket:
            raise ValueError("Workspace bucket not found. Check WORKSPACE_BUCKET environment variable if none provided.")

        print("version: " + self.version)
        print("bucket: " + self.bucket)
            
    def find_diagnosis_codes(self, 
                            vocabularies: List[str] = None,
                            search_terms: List[str] = None, 
                            exclude_terms: List[str] = None,
                            exact_codes: Dict[str, List[str]] = None,
                            pattern_codes: Dict[str, List[str]] = None,
                            exclude_codes: Dict[str, List[str]] = None,
                            person_ids: Union[List[int], str] = None):
        """
        Find diagnosis codes in condition_occurrence and observation tables.
        """

        if not vocabularies:
            # Default to standard vocabularies if nothing specified
            vocabularies = ["ICD9CM", "ICD10CM", "SNOMED"]

        # Add any vocabularies from exact_codes or pattern_codes that aren't already included
        vocab_additions = set()
        for vocab in list(exact_codes.keys() if exact_codes else []) + list(pattern_codes.keys() if pattern_codes else []):
            if vocab not in vocabularies:
                vocab_additions.add(vocab)
                vocabularies.append(vocab)

        # Notify if vocabularies were added
        if vocab_additions:
            print(f"Added vocabularies from code specifications: {', '.join(vocab_additions)}")

        # Store query parameters for summary
        self.last_query_params = {
            "type": "diagnosis_codes",
            "vocabularies": vocabularies,
            "search_terms": search_terms,
            "exclude_terms": exclude_terms,
            "exact_codes": exact_codes,
            "pattern_codes": pattern_codes,
            "exclude_codes": exclude_codes,
            "person_ids": person_ids
        }
        
        # Generate query summary text
        self.query_summary = self._generate_query_summary()
    
        # Create Jinja2 environment and add filter
        env = Environment()
        
        # Define and add filter to environment
        def join_quotes(items):
            return ', '.join(f"'{item}'" for item in items)

        env.filters['join_quotes'] = join_quotes

        template = env.from_string("""
        WITH filtered_concepts AS (
            SELECT DISTINCT concept_id, vocabulary_id, concept_code, concept_name
            FROM {{ version }}.concept
            WHERE 
                vocabulary_id IN ({{ vocabularies|join_quotes }})
                {%- if exclude_terms %}
                AND (
                    {%- for term in exclude_terms %}
                    {{ "AND " if not loop.first else "" }}LOWER(concept_name) NOT LIKE '%{{ term|lower }}%'
                    {%- endfor %}
                )
                {%- endif %}
                {%- for vocab, codes in exclude_codes.items() %}
                {%- if codes %}
                AND NOT (
                    vocabulary_id = '{{ vocab }}' AND concept_code IN ({{ codes|join_quotes }})
                )
                {%- endif %}
                {%- endfor %}
                {%- if search_terms or exact_codes or pattern_codes %}
                AND (
                    {%- if search_terms %}
                    (
                        {%- for term in search_terms %}
                        {{ "OR " if not loop.first else "" }}LOWER(concept_name) LIKE '%{{ term|lower }}%'
                        {%- endfor %}
                    )
                    {%- endif %}

                    {%- if search_terms and (exact_codes or pattern_codes) %}OR{% endif %}

                    {%- if exact_codes %}
                    (
                        {%- for vocab, codes in exact_codes.items() %}
                        {%- if codes %}
                        {{ "OR " if not loop.first else "" }}(vocabulary_id = '{{ vocab }}' AND concept_code IN ({{ codes|join_quotes }}))
                        {%- endif %}
                        {%- endfor %}
                    )
                    {%- endif %}

                    {%- if (search_terms or exact_codes) and pattern_codes %}OR{% endif %}

                    {%- if pattern_codes %}
                    (
                        {%- for vocab, patterns in pattern_codes.items() %}
                        {%- if patterns %}
                        {{ "OR " if not loop.first else "" }}(vocabulary_id = '{{ vocab }}' AND (
                            {%- for pattern in patterns %}
                            {{ "OR " if not loop.first else "" }}concept_code LIKE '{{ pattern }}'
                            {%- endfor %}
                        ))
                        {%- endif %}
                        {%- endfor %}
                    )
                    {%- endif %}
                )
                {%- endif %}
        ),
        -- Regular events (non-V codes)
        regular_codes AS (
            -- Get codes from condition_occurrence via source_value (non-V codes)
            SELECT
                co.person_id,
                co.condition_source_value AS concept_code,
                fc.vocabulary_id,
                fc.concept_name,
                'condition' AS domain
            FROM {{ version }}.condition_occurrence co
            JOIN filtered_concepts fc 
                ON co.condition_source_value = fc.concept_code
            WHERE NOT co.condition_source_value LIKE 'V%'
            {%- if person_ids %}
            AND co.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- Get codes from condition_occurrence via source_concept_id (non-V codes)
            SELECT
                co.person_id,
                fc.concept_code,
                fc.vocabulary_id,
                fc.concept_name,
                'condition' AS domain
            FROM {{ version }}.condition_occurrence co
            JOIN filtered_concepts fc 
                ON co.condition_source_concept_id = fc.concept_id
            WHERE NOT fc.concept_code LIKE 'V%'
            {%- if person_ids %}
            AND co.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- Get codes from observation via source_value (non-V codes)
            SELECT
                o.person_id,
                o.observation_source_value AS concept_code,
                fc.vocabulary_id,
                fc.concept_name,
                'observation' AS domain
            FROM {{ version }}.observation o
            JOIN filtered_concepts fc 
                ON o.observation_source_value = fc.concept_code
            WHERE NOT o.observation_source_value LIKE 'V%'
            {%- if person_ids %}
            AND o.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- Get codes from observation via source_concept_id (non-V codes)
            SELECT
                o.person_id,
                fc.concept_code,
                fc.vocabulary_id,
                fc.concept_name,
                'observation' AS domain
            FROM {{ version }}.observation o
            JOIN filtered_concepts fc 
                ON o.observation_source_concept_id = fc.concept_id
            WHERE NOT fc.concept_code LIKE 'V%'
            {%- if person_ids %}
            AND o.person_id IN ({{ person_ids }})
            {%- endif %}
        ),

        -- V code events from condition_occurrence with special handling
        v_codes AS (
            -- source_value path
            SELECT
                co.person_id,
                co.condition_source_value AS concept_code,
                co.condition_concept_id AS concept_id,
                'condition' AS domain
            FROM {{ version }}.condition_occurrence co
            JOIN filtered_concepts fc 
                ON co.condition_source_value = fc.concept_code
            WHERE co.condition_source_value LIKE 'V%'
            {%- if person_ids %}
            AND co.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- source_concept_id path
            SELECT
                co.person_id,
                fc.concept_code,
                co.condition_concept_id AS concept_id,
                'condition' AS domain
            FROM {{ version }}.condition_occurrence co
            JOIN filtered_concepts fc 
                ON co.condition_source_concept_id = fc.concept_id
            WHERE fc.concept_code LIKE 'V%'
            {%- if person_ids %}
            AND co.person_id IN ({{ person_ids }})
            {%- endif %}
            
            UNION ALL

            -- source_value path
            SELECT
                o.person_id,
                o.observation_source_value AS concept_code,
                o.observation_concept_id AS concept_id,
                'observation' AS domain
            FROM {{ version }}.observation o
            JOIN filtered_concepts fc 
                ON o.observation_source_value = fc.concept_code
            WHERE o.observation_source_value LIKE 'V%'
            {%- if person_ids %}
            AND o.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- source_concept_id path
            SELECT
                o.person_id,
                fc.concept_code,
                o.observation_concept_id AS concept_id,
                'observation' AS domain
            FROM {{ version }}.observation o
            JOIN filtered_concepts fc 
                ON o.observation_source_concept_id = fc.concept_id
            WHERE fc.concept_code LIKE 'V%'
            {%- if person_ids %}
            AND o.person_id IN ({{ person_ids }})
            {%- endif %}
        ),

        -- Apply correct vocabulary attribution for V codes
        v_codes_corrected AS (
            SELECT
                v.person_id,
                v.concept_code,
                c.vocabulary_id,
                c.concept_name,
                v.domain
            FROM v_codes v
            JOIN {{ version }}.concept_relationship cr
                ON v.concept_id = cr.concept_id_1
            JOIN {{ version }}.concept c
                ON cr.concept_id_2 = c.concept_id
            WHERE c.vocabulary_id IN ({{ vocabularies|join_quotes }})
            AND v.concept_code = c.concept_code
        ),

        -- Combine regular events with corrected V code events
        all_events AS (
            SELECT * FROM regular_codes
            UNION ALL
            SELECT * FROM v_codes_corrected
        )

        -- Final output with counts
        SELECT
            vocabulary_id,
            concept_code,
            concept_name,
            COUNT(DISTINCT person_id) AS unique_persons,
            COUNT(*) AS total_events,
        FROM all_events
        GROUP BY vocabulary_id, concept_code, concept_name
        ORDER BY vocabulary_id ASC, concept_code ASC
        """)       
        
        # Format person_ids if it's a list
        if isinstance(person_ids, list):
            person_ids = ', '.join(str(pid) for pid in person_ids)
            
        # Initialize dictionaries
        exact_codes = exact_codes or {}
        pattern_codes = pattern_codes or {}
        exclude_codes = exclude_codes or {}
        
        # Render the template
        self.query_text = template.render(
            version=self.version,
            vocabularies=vocabularies,
            search_terms=search_terms,
            exclude_terms=exclude_terms,
            exact_codes=exact_codes,
            pattern_codes=pattern_codes,
            exclude_codes=exclude_codes,
            person_ids=person_ids
        )
        
        return self
    
    def person_code_df(self,
                       name: str,
                       vocabularies: List[str] = None,
                       search_terms: List[str] = None, 
                       exclude_terms: List[str] = None,
                       exact_codes: Dict[str, List[str]] = None,
                       pattern_codes: Dict[str, List[str]] = None,
                       exclude_codes: Dict[str, List[str]] = None,
                       person_ids: Union[List[int], str] = None,
                       dates: bool = False):
        """
        Creates a person-level dataframe with diagnosis information.

        Args:
            name: String to use as column name prefix (e.g., 'htn')
            vocabularies: List of vocabulary IDs to search
            search_terms: List of terms to search in concept names
            exclude_terms: List of terms to exclude from concept names
            exact_codes: Dictionary mapping vocabularies to lists of exact codes
            pattern_codes: Dictionary mapping vocabularies to lists of code patterns
            exclude_codes: Dictionary mapping vocabularies to lists of codes to exclude
            person_ids: List of person IDs to filter to, or a string expression
            dates: Whether to include date information (first, second, last)

        Returns:
            Polars DataFrame with person_id, {name}_n, vocab, and date columns if requested
        """

        if not vocabularies:
            # Default to standard vocabularies if nothing specified
            vocabularies = ["ICD9CM", "ICD10CM", "SNOMED"]

        # Add any vocabularies from exact_codes or pattern_codes that aren't already included
        vocab_additions = set()
        for vocab in list(exact_codes.keys() if exact_codes else []) + list(pattern_codes.keys() if pattern_codes else []):
            if vocab not in vocabularies:
                vocab_additions.add(vocab)
                vocabularies.append(vocab)

        # Notify if vocabularies were added
        if vocab_additions:
            print(f"Added vocabularies from code specifications: {', '.join(vocab_additions)}")

        # Store query parameters for summary
        self.last_query_params = {
            "type": "person_diagnosis",
            "name": name,
            "vocabularies": vocabularies,
            "search_terms": search_terms,
            "exclude_terms": exclude_terms,
            "exact_codes": exact_codes,
            "pattern_codes": pattern_codes,
            "exclude_codes": exclude_codes,
            "person_ids": person_ids,
            "dates": dates
        }

        # Generate query summary text
        self.query_summary = self._generate_query_summary()

        # Create Jinja2 environment and add filter
        env = Environment()

        # Define and add filter to environment
        def join_quotes(items):
            return ', '.join(f"'{item}'" for item in items)

        env.filters['join_quotes'] = join_quotes

        template = env.from_string("""
        WITH filtered_concepts AS (
            SELECT DISTINCT concept_id, vocabulary_id, concept_code, concept_name
            FROM {{ version }}.concept
            WHERE 
                vocabulary_id IN ({{ vocabularies|join_quotes }})
                {%- if exclude_terms %}
                AND (
                    {%- for term in exclude_terms %}
                    {{ "AND " if not loop.first else "" }}LOWER(concept_name) NOT LIKE '%{{ term|lower }}%'
                    {%- endfor %}
                )
                {%- endif %}
                {%- for vocab, codes in exclude_codes.items() %}
                {%- if codes %}
                AND NOT (
                    vocabulary_id = '{{ vocab }}' AND concept_code IN ({{ codes|join_quotes }})
                )
                {%- endif %}
                {%- endfor %}
                {%- if search_terms or exact_codes or pattern_codes %}
                AND (
                    {%- if search_terms %}
                    (
                        {%- for term in search_terms %}
                        {{ "OR " if not loop.first else "" }}LOWER(concept_name) LIKE '%{{ term|lower }}%'
                        {%- endfor %}
                    )
                    {%- endif %}

                    {%- if search_terms and (exact_codes or pattern_codes) %}OR{% endif %}

                    {%- if exact_codes %}
                    (
                        {%- for vocab, codes in exact_codes.items() %}
                        {%- if codes %}
                        {{ "OR " if not loop.first else "" }}(vocabulary_id = '{{ vocab }}' AND concept_code IN ({{ codes|join_quotes }}))
                        {%- endif %}
                        {%- endfor %}
                    )
                    {%- endif %}

                    {%- if (search_terms or exact_codes) and pattern_codes %}OR{% endif %}

                    {%- if pattern_codes %}
                    (
                        {%- for vocab, patterns in pattern_codes.items() %}
                        {%- if patterns %}
                        {{ "OR " if not loop.first else "" }}(vocabulary_id = '{{ vocab }}' AND (
                            {%- for pattern in patterns %}
                            {{ "OR " if not loop.first else "" }}concept_code LIKE '{{ pattern }}'
                            {%- endfor %}
                        ))
                        {%- endif %}
                        {%- endfor %}
                    )
                    {%- endif %}
                )
                {%- endif %}
        ),
        -- Regular events (non-V codes)
        regular_events AS (
            -- Get codes from condition_occurrence via source_value (non-V codes)
            SELECT
                co.person_id,
                co.condition_source_value AS concept_code,
                fc.vocabulary_id,
                fc.concept_name,
                'condition' AS domain,
                co.condition_start_date AS event_date
            FROM {{ version }}.condition_occurrence co
            JOIN filtered_concepts fc 
                ON co.condition_source_value = fc.concept_code
            WHERE NOT co.condition_source_value LIKE 'V%'
            {%- if person_ids %}
            AND co.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- Get codes from condition_occurrence via source_concept_id (non-V codes)
            SELECT
                co.person_id,
                fc.concept_code,
                fc.vocabulary_id,
                fc.concept_name,
                'condition' AS domain,
                co.condition_start_date AS event_date
            FROM {{ version }}.condition_occurrence co
            JOIN filtered_concepts fc 
                ON co.condition_source_concept_id = fc.concept_id
            WHERE NOT fc.concept_code LIKE 'V%'
            {%- if person_ids %}
            AND co.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- Get codes from observation via source_value (non-V codes)
            SELECT
                o.person_id,
                o.observation_source_value AS concept_code,
                fc.vocabulary_id,
                fc.concept_name,
                'observation' AS domain,
                o.observation_date AS event_date
            FROM {{ version }}.observation o
            JOIN filtered_concepts fc 
                ON o.observation_source_value = fc.concept_code
            WHERE NOT o.observation_source_value LIKE 'V%'
            {%- if person_ids %}
            AND o.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- Get codes from observation via source_concept_id (non-V codes)
            SELECT
                o.person_id,
                fc.concept_code,
                fc.vocabulary_id,
                fc.concept_name,
                'observation' AS domain,
                o.observation_date AS event_date
            FROM {{ version }}.observation o
            JOIN filtered_concepts fc 
                ON o.observation_source_concept_id = fc.concept_id
            WHERE NOT fc.concept_code LIKE 'V%'
            {%- if person_ids %}
            AND o.person_id IN ({{ person_ids }})
            {%- endif %}
        ),

        -- V code events from condition_occurrence with special handling
        v_events AS (
            -- source_value path
            SELECT
                co.person_id,
                co.condition_source_value AS concept_code,
                co.condition_concept_id AS concept_id,
                'condition' AS domain,
                co.condition_start_date AS event_date
            FROM {{ version }}.condition_occurrence co
            JOIN filtered_concepts fc 
                ON co.condition_source_value = fc.concept_code
            WHERE co.condition_source_value LIKE 'V%'
            {%- if person_ids %}
            AND co.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- source_concept_id path
            SELECT
                co.person_id,
                fc.concept_code,
                co.condition_concept_id AS concept_id,
                'condition' AS domain,
                co.condition_start_date AS event_date
            FROM {{ version }}.condition_occurrence co
            JOIN filtered_concepts fc 
                ON co.condition_source_concept_id = fc.concept_id
            WHERE fc.concept_code LIKE 'V%'
            {%- if person_ids %}
            AND co.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- source_value path
            SELECT
                o.person_id,
                o.observation_source_value AS concept_code,
                o.observation_concept_id AS concept_id,
                'observation' AS domain,
                o.observation_date AS event_date
            FROM {{ version }}.observation o
            JOIN filtered_concepts fc 
                ON o.observation_source_value = fc.concept_code
            WHERE o.observation_source_value LIKE 'V%'
            {%- if person_ids %}
            AND o.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- source_concept_id path
            SELECT
                o.person_id,
                fc.concept_code,
                o.observation_concept_id AS concept_id,
                'observation' AS domain,
                o.observation_date AS event_date
            FROM {{ version }}.observation o
            JOIN filtered_concepts fc 
                ON o.observation_source_concept_id = fc.concept_id
            WHERE fc.concept_code LIKE 'V%'
            {%- if person_ids %}
            AND o.person_id IN ({{ person_ids }})
            {%- endif %}
        ),

        -- Apply correct vocabulary attribution for V codes
        v_events_corrected AS (
            SELECT
                v.person_id,
                v.concept_code,
                c.vocabulary_id,
                c.concept_name,
                v.domain,
                v.event_date
            FROM v_events v
            JOIN {{ version }}.concept_relationship cr
                ON v.concept_id = cr.concept_id_1
            JOIN {{ version }}.concept c
                ON cr.concept_id_2 = c.concept_id
            WHERE c.vocabulary_id IN ({{ vocabularies|join_quotes }})
            AND v.concept_code = c.concept_code
        ),

        -- Combine regular events with corrected V code events
        all_events AS (
            SELECT * FROM regular_events
            UNION ALL
            SELECT * FROM v_events_corrected
        ),

        -- Get distinct dates per person
        distinct_dates AS (
            SELECT 
                person_id,
                vocabulary_id,
                event_date
            FROM all_events
            GROUP BY person_id, vocabulary_id, event_date
        ),

        -- Calculate vocabulary summary per person
        vocab_summary AS (
            SELECT
                person_id,
                STRING_AGG(CASE 
                    WHEN vocabulary_id = 'ICD9CM' THEN '9'
                    WHEN vocabulary_id = 'ICD10CM' THEN '10'
                    WHEN vocabulary_id = 'SNOMED' THEN 'SNO'
                    ELSE SUBSTR(vocabulary_id, 1, 3)
                END, ' ' ORDER BY vocabulary_id) AS vocab
            FROM (
                SELECT DISTINCT person_id, vocabulary_id
                FROM all_events
            )
            GROUP BY person_id
        ),

        -- Count distinct dates per person
        date_counts AS (
            SELECT
                person_id,
                COUNT(DISTINCT event_date) AS {{ name }}_n
            FROM all_events
            GROUP BY person_id
        ){% if dates %},

        -- Get ordered dates per person
        ordered_dates AS (
            SELECT
                person_id,
                event_date,
                ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY event_date ASC) AS date_order,
                ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY event_date DESC) AS rev_date_order
            FROM (
                SELECT DISTINCT person_id, event_date
                FROM all_events
            )
        ),

        -- Get first, second, and last dates
        key_dates AS (
            SELECT
                person_id,
                MAX(CASE WHEN date_order = 1 THEN event_date END) AS {{ name }}_1,
                MAX(CASE WHEN date_order = 2 THEN event_date END) AS {{ name }}_2,
                MAX(CASE WHEN rev_date_order = 1 THEN event_date END) AS {{ name }}_last
            FROM ordered_dates
            GROUP BY person_id
        ){% endif %}

        -- Final person-level output
        SELECT
            dc.person_id,
            dc.{{ name }}_n,
            vs.vocab{% if dates %},
            kd.{{ name }}_1,
            kd.{{ name }}_2,
            kd.{{ name }}_last{% endif %}
        FROM date_counts dc
        JOIN vocab_summary vs ON dc.person_id = vs.person_id
        {% if dates %}JOIN key_dates kd ON dc.person_id = kd.person_id{% endif %}
        ORDER BY dc.person_id
        """)       

        # Format person_ids if it's a list
        if isinstance(person_ids, list):
            person_ids = ', '.join(str(pid) for pid in person_ids)

        # Initialize dictionaries
        exact_codes = exact_codes or {}
        pattern_codes = pattern_codes or {}
        exclude_codes = exclude_codes or {}

        # Render the template
        self.query_text = template.render(
            version=self.version,
            name=name,
            vocabularies=vocabularies,
            search_terms=search_terms,
            exclude_terms=exclude_terms,
            exact_codes=exact_codes,
            pattern_codes=pattern_codes,
            exclude_codes=exclude_codes,
            person_ids=person_ids,
            dates=dates
        )

        return self
    
    def find_procedure_codes(self, 
                            search_tables: List[str] = None,
                            vocabularies: List[str] = None,
                            concept_classes: List[str] = None,
                            search_terms: List[str] = None, 
                            exclude_terms: List[str] = None,
                            exact_codes: Dict[str, List[str]] = None,
                            pattern_codes: Dict[str, List[str]] = None,
                            exclude_codes: Dict[str, List[str]] = None,
                            person_ids: Union[List[int], str] = None):
        """
        Find procedure codes in procedure_occurrence and/or observation tables.

        Args:
            search_tables: Tables to search. Options: ['procedure_occurrence', 'observation']
                          Default: ['procedure_occurrence', 'observation'] 
            vocabularies: Procedure vocabularies to search 
                         Default: ["CPT4", "HCPCS", "ICD10PCS", "ICD9Proc", "LOINC", "SNOMED"]
            concept_classes: List of concept classes to filter to (optional)
                            Default: None (no filtering - returns all concept classes)
            search_terms: List of terms to search in concept names
            exclude_terms: List of terms to exclude from concept names  
            exact_codes: Dictionary mapping vocabularies to lists of exact codes
            pattern_codes: Dictionary mapping vocabularies to lists of code patterns (use % wildcards)
            exclude_codes: Dictionary mapping vocabularies to lists of codes to exclude
            person_ids: List of person IDs to filter to, or a string expression

        Returns:
            Self (for method chaining), call execute_gbq() to run the query
        """

        if not search_tables:
            search_tables = ['procedure_occurrence', 'observation']

        if not vocabularies:
            vocabularies = ["CPT4", "HCPCS", "ICD10PCS", "ICD9Proc", "Nebraska Lexicon", "LOINC", "SNOMED"]

        # concept_classes are not included by default, but these are the recommended classes for vocabularies
        # CPT4: "CPT4", "CPT4 Hierarchy"
        # HCPCS: "HCPCS"
        # ICD10: "ICD10PCS", "ICD10PCS Hierarchy",
        # ICD9: "2-dig nonbill code", "3-dig billing code", "3-dig nonbill code", "4-dig billing code", "Procedure"
        # LOINC, Nebraska Lexicon, SNOMED: "Procedure"
            
        # Add any vocabularies from exact_codes or pattern_codes that aren't already included
        vocab_additions = set()
        for vocab in list(exact_codes.keys() if exact_codes else []) + list(pattern_codes.keys() if pattern_codes else []):
            if vocab not in vocabularies:
                vocab_additions.add(vocab)
                vocabularies.append(vocab)

        # Notify if vocabularies were added
        if vocab_additions:
            print(f"Added vocabularies from code specifications: {', '.join(vocab_additions)}")

        # Store query parameters for summary
        self.last_query_params = {
            "type": "procedure_codes",
            "search_tables": search_tables,
            "vocabularies": vocabularies,
            "concept_classes": concept_classes,
            "search_terms": search_terms,
            "exclude_terms": exclude_terms,
            "exact_codes": exact_codes,
            "pattern_codes": pattern_codes,
            "exclude_codes": exclude_codes,
            "person_ids": person_ids
        }

        # Generate query summary text
        self.query_summary = self._generate_query_summary()

        # Create Jinja2 environment and add filter
        env = Environment()

        # Define and add filter to environment
        def join_quotes(items):
            return ', '.join(f"'{item}'" for item in items)

        env.filters['join_quotes'] = join_quotes

        template = env.from_string("""
        WITH filtered_concepts AS (
            SELECT DISTINCT concept_id, vocabulary_id, concept_class_id, concept_code, concept_name
            FROM {{ version }}.concept
            WHERE 
                vocabulary_id IN ({{ vocabularies|join_quotes }})
                {%- if concept_classes %}
                AND concept_class_id IN ({{ concept_classes|join_quotes}})
                {%- endif %}
                {%- if exclude_terms %}
                AND (
                    {%- for term in exclude_terms %}
                    {{ "AND " if not loop.first else "" }}LOWER(concept_name) NOT LIKE '%{{ term|lower }}%'
                    {%- endfor %}
                )
                {%- endif %}
                {%- for vocab, codes in exclude_codes.items() %}
                {%- if codes %}
                AND NOT (
                    vocabulary_id = '{{ vocab }}' AND concept_code IN ({{ codes|join_quotes }})
                )
                {%- endif %}
                {%- endfor %}
                {%- if search_terms or exact_codes or pattern_codes %}
                AND (
                    {%- if search_terms %}
                    (
                        {%- for term in search_terms %}
                        {{ "OR " if not loop.first else "" }}LOWER(concept_name) LIKE '%{{ term|lower }}%'
                        {%- endfor %}
                    )
                    {%- endif %}

                    {%- if search_terms and (exact_codes or pattern_codes) %}OR{% endif %}

                    {%- if exact_codes %}
                    (
                        {%- for vocab, codes in exact_codes.items() %}
                        {%- if codes %}
                        {{ "OR " if not loop.first else "" }}(vocabulary_id = '{{ vocab }}' AND concept_code IN ({{ codes|join_quotes }}))
                        {%- endif %}
                        {%- endfor %}
                    )
                    {%- endif %}

                    {%- if (search_terms or exact_codes) and pattern_codes %}OR{% endif %}

                    {%- if pattern_codes %}
                    (
                        {%- for vocab, patterns in pattern_codes.items() %}
                        {%- if patterns %}
                        {{ "OR " if not loop.first else "" }}(vocabulary_id = '{{ vocab }}' AND (
                            {%- for pattern in patterns %}
                            {{ "OR " if not loop.first else "" }}concept_code LIKE '{{ pattern }}'
                            {%- endfor %}
                        ))
                        {%- endif %}
                        {%- endfor %}
                    )
                    {%- endif %}
                )
                {%- endif %}
        ),

        -- Procedure events from procedure_occurrence
        {% if 'procedure_occurrence' in search_tables %}
        procedure_table_codes AS (
            -- Procedures via source_value
            SELECT
                po.person_id,
                po.procedure_source_value AS concept_code,
                fc.vocabulary_id,
                fc.concept_class_id,
                fc.concept_name,
                'procedure_occurrence' AS source_table
            FROM {{ version }}.procedure_occurrence po
            JOIN filtered_concepts fc 
                ON po.procedure_source_value = fc.concept_code
            {%- if person_ids %}
            WHERE po.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- Procedures via source_concept_id
            SELECT
                po.person_id,
                fc.concept_code,
                fc.vocabulary_id,
                fc.concept_class_id,
                fc.concept_name,
                'procedure_occurrence' AS source_table
            FROM {{ version }}.procedure_occurrence po
            JOIN filtered_concepts fc 
                ON po.procedure_source_concept_id = fc.concept_id
            {%- if person_ids %}
            WHERE po.person_id IN ({{ person_ids }})
            {%- endif %}
        ),
        {% endif %}

        -- Procedure events from observation
        {% if 'observation' in search_tables %}
        observation_procedure_codes AS (
            -- Procedures from observation via source_value
            SELECT
                o.person_id,
                o.observation_source_value AS concept_code,
                fc.vocabulary_id,
                fc.concept_class_id,
                fc.concept_name,
                'observation' AS source_table
            FROM {{ version }}.observation o
            JOIN filtered_concepts fc 
                ON o.observation_source_value = fc.concept_code
            {%- if person_ids %}
            WHERE o.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- Procedures from observation via source_concept_id
            SELECT
                o.person_id,
                fc.concept_code,
                fc.vocabulary_id,
                fc.concept_class_id,
                fc.concept_name,
                'observation' AS source_table
            FROM {{ version }}.observation o
            JOIN filtered_concepts fc 
                ON o.observation_source_concept_id = fc.concept_id
            {%- if person_ids %}
            WHERE o.person_id IN ({{ person_ids }})
            {%- endif %}
        ),
        {% endif %}

        -- Combine results from all searched tables
        all_procedure_events AS (
            {% if 'procedure_occurrence' in search_tables %}
            SELECT * FROM procedure_table_codes
            {% endif %}
            {% if 'observation' in search_tables %}
            {% if 'procedure_occurrence' in search_tables %}UNION ALL{% endif %}
            SELECT * FROM observation_procedure_codes
            {% endif %}
        )

        -- Final output with counts
        SELECT
            vocabulary_id,
            concept_class_id,
            concept_code,
            concept_name,
            source_table,
            COUNT(DISTINCT person_id) AS unique_persons,
            COUNT(*) AS total_events
        FROM all_procedure_events
        GROUP BY vocabulary_id, concept_class_id, concept_code, concept_name, source_table
        ORDER BY vocabulary_id ASC, concept_class_id ASC, concept_code ASC, source_table ASC
        """)       

        # Format person_ids if it's a list
        if isinstance(person_ids, list):
            person_ids = ', '.join(str(pid) for pid in person_ids)

        # Initialize dictionaries
        exact_codes = exact_codes or {}
        pattern_codes = pattern_codes or {}
        exclude_codes = exclude_codes or {}

        # Render the template
        self.query_text = template.render(
            version=self.version,
            search_tables=search_tables,
            vocabularies=vocabularies,
            concept_classes=concept_classes,
            search_terms=search_terms,
            exclude_terms=exclude_terms,
            exact_codes=exact_codes,
            pattern_codes=pattern_codes,
            exclude_codes=exclude_codes,
            person_ids=person_ids
        )

        return self
    
    def person_procedure_df(self,
                            name: str,
                            search_tables: List[str] = None,
                            vocabularies: List[str] = None,
                            concept_classes: List[str] = None,
                            search_terms: List[str] = None, 
                            exclude_terms: List[str] = None,
                            exact_codes: Dict[str, List[str]] = None,
                            pattern_codes: Dict[str, List[str]] = None,
                            exclude_codes: Dict[str, List[str]] = None,
                            person_ids: Union[List[int], str] = None,
                            dates: bool = False):
        """
        Creates a person-level dataframe with procedure information.

        Args:
            name: String to use as column name prefix (e.g., 'surgery')
            search_tables: Tables to search. Options: ['procedure_occurrence', 'observation']
                          Default: ['procedure_occurrence', 'observation'] 
            vocabularies: Procedure vocabularies to search 
                         Default: ["CPT4", "HCPCS", "ICD10PCS", "ICD9Proc", "Nebraska Lexicon", "LOINC", "SNOMED"]
            concept_classes: List of concept classes to filter to (optional)
                            Default: None (no filtering - returns all concept classes)
            search_terms: List of terms to search in concept names
            exclude_terms: List of terms to exclude from concept names  
            exact_codes: Dictionary mapping vocabularies to lists of exact codes
            pattern_codes: Dictionary mapping vocabularies to lists of code patterns (use % wildcards)
            exclude_codes: Dictionary mapping vocabularies to lists of codes to exclude
            person_ids: List of person IDs to filter to, or a string expression
            dates: Whether to include date information (first, second, last)

        Returns:
            Polars DataFrame with person_id, {name}_n, vocab, and date columns if requested
        """

        if not search_tables:
            search_tables = ['procedure_occurrence', 'observation']

        if not vocabularies:
            vocabularies = ["CPT4", "HCPCS", "ICD10PCS", "ICD9Proc", "Nebraska Lexicon", "LOINC", "SNOMED"]

        # Add any vocabularies from exact_codes or pattern_codes that aren't already included
        vocab_additions = set()
        for vocab in list(exact_codes.keys() if exact_codes else []) + list(pattern_codes.keys() if pattern_codes else []):
            if vocab not in vocabularies:
                vocab_additions.add(vocab)
                vocabularies.append(vocab)

        # Notify if vocabularies were added
        if vocab_additions:
            print(f"Added vocabularies from code specifications: {', '.join(vocab_additions)}")

        # Store query parameters for summary
        self.last_query_params = {
            "type": "person_procedure",
            "name": name,
            "search_tables": search_tables,
            "vocabularies": vocabularies,
            "concept_classes": concept_classes,
            "search_terms": search_terms,
            "exclude_terms": exclude_terms,
            "exact_codes": exact_codes,
            "pattern_codes": pattern_codes,
            "exclude_codes": exclude_codes,
            "person_ids": person_ids,
            "dates": dates
        }

        # Generate query summary text
        self.query_summary = self._generate_query_summary()

        # Create Jinja2 environment and add filter
        env = Environment()

        # Define and add filter to environment
        def join_quotes(items):
            return ', '.join(f"'{item}'" for item in items)

        env.filters['join_quotes'] = join_quotes

        template = env.from_string("""
        WITH filtered_concepts AS (
            SELECT DISTINCT concept_id, vocabulary_id, concept_class_id, concept_code, concept_name
            FROM {{ version }}.concept
            WHERE 
                vocabulary_id IN ({{ vocabularies|join_quotes }})
                {%- if concept_classes %}
                AND concept_class_id IN ({{ concept_classes|join_quotes}})
                {%- endif %}
                {%- if exclude_terms %}
                AND (
                    {%- for term in exclude_terms %}
                    {{ "AND " if not loop.first else "" }}LOWER(concept_name) NOT LIKE '%{{ term|lower }}%'
                    {%- endfor %}
                )
                {%- endif %}
                {%- for vocab, codes in exclude_codes.items() %}
                {%- if codes %}
                AND NOT (
                    vocabulary_id = '{{ vocab }}' AND concept_code IN ({{ codes|join_quotes }})
                )
                {%- endif %}
                {%- endfor %}
                {%- if search_terms or exact_codes or pattern_codes %}
                AND (
                    {%- if search_terms %}
                    (
                        {%- for term in search_terms %}
                        {{ "OR " if not loop.first else "" }}LOWER(concept_name) LIKE '%{{ term|lower }}%'
                        {%- endfor %}
                    )
                    {%- endif %}

                    {%- if search_terms and (exact_codes or pattern_codes) %}OR{% endif %}

                    {%- if exact_codes %}
                    (
                        {%- for vocab, codes in exact_codes.items() %}
                        {%- if codes %}
                        {{ "OR " if not loop.first else "" }}(vocabulary_id = '{{ vocab }}' AND concept_code IN ({{ codes|join_quotes }}))
                        {%- endif %}
                        {%- endfor %}
                    )
                    {%- endif %}

                    {%- if (search_terms or exact_codes) and pattern_codes %}OR{% endif %}

                    {%- if pattern_codes %}
                    (
                        {%- for vocab, patterns in pattern_codes.items() %}
                        {%- if patterns %}
                        {{ "OR " if not loop.first else "" }}(vocabulary_id = '{{ vocab }}' AND (
                            {%- for pattern in patterns %}
                            {{ "OR " if not loop.first else "" }}concept_code LIKE '{{ pattern }}'
                            {%- endfor %}
                        ))
                        {%- endif %}
                        {%- endfor %}
                    )
                    {%- endif %}
                )
                {%- endif %}
        ),

        -- Procedure events from procedure_occurrence
        {% if 'procedure_occurrence' in search_tables %}
        procedure_table_events AS (
            -- Procedures via source_value
            SELECT
                po.person_id,
                po.procedure_source_value AS concept_code,
                fc.vocabulary_id,
                fc.concept_class_id,
                fc.concept_name,
                'procedure_occurrence' AS source_table,
                po.procedure_date AS event_date
            FROM {{ version }}.procedure_occurrence po
            JOIN filtered_concepts fc 
                ON po.procedure_source_value = fc.concept_code
            {%- if person_ids %}
            WHERE po.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- Procedures via source_concept_id
            SELECT
                po.person_id,
                fc.concept_code,
                fc.vocabulary_id,
                fc.concept_class_id,
                fc.concept_name,
                'procedure_occurrence' AS source_table,
                po.procedure_date AS event_date
            FROM {{ version }}.procedure_occurrence po
            JOIN filtered_concepts fc 
                ON po.procedure_source_concept_id = fc.concept_id
            {%- if person_ids %}
            WHERE po.person_id IN ({{ person_ids }})
            {%- endif %}
        ),
        {% endif %}

        -- Procedure events from observation
        {% if 'observation' in search_tables %}
        observation_procedure_events AS (
            -- Procedures from observation via source_value
            SELECT
                o.person_id,
                o.observation_source_value AS concept_code,
                fc.vocabulary_id,
                fc.concept_class_id,
                fc.concept_name,
                'observation' AS source_table,
                o.observation_date AS event_date
            FROM {{ version }}.observation o
            JOIN filtered_concepts fc 
                ON o.observation_source_value = fc.concept_code
            {%- if person_ids %}
            WHERE o.person_id IN ({{ person_ids }})
            {%- endif %}

            UNION ALL

            -- Procedures from observation via source_concept_id
            SELECT
                o.person_id,
                fc.concept_code,
                fc.vocabulary_id,
                fc.concept_class_id,
                fc.concept_name,
                'observation' AS source_table,
                o.observation_date AS event_date
            FROM {{ version }}.observation o
            JOIN filtered_concepts fc 
                ON o.observation_source_concept_id = fc.concept_id
            {%- if person_ids %}
            WHERE o.person_id IN ({{ person_ids }})
            {%- endif %}
        ),
        {% endif %}

        -- Combine results from all searched tables
        all_procedure_events AS (
            {% if 'procedure_occurrence' in search_tables %}
            SELECT * FROM procedure_table_events
            {% endif %}
            {% if 'observation' in search_tables %}
            {% if 'procedure_occurrence' in search_tables %}UNION ALL{% endif %}
            SELECT * FROM observation_procedure_events
            {% endif %}
        ),

        -- Get distinct dates per person
        distinct_dates AS (
            SELECT 
                person_id,
                vocabulary_id,
                event_date
            FROM all_procedure_events
            GROUP BY person_id, vocabulary_id, event_date
        ),

        -- Calculate vocabulary summary per person
        vocab_summary AS (
            SELECT
                person_id,
                STRING_AGG(CASE 
                    WHEN vocabulary_id = 'CPT4' THEN 'CPT'
                    WHEN vocabulary_id = 'HCPCS' THEN 'HC'
                    WHEN vocabulary_id = 'ICD10PCS' THEN '10P'
                    WHEN vocabulary_id = 'ICD9Proc' THEN '9P'
                    WHEN vocabulary_id = 'LOINC' THEN 'LN'
                    WHEN vocabulary_id = 'SNOMED' THEN 'SNO'
                    ELSE SUBSTR(vocabulary_id, 1, 3)
                END, ' ' ORDER BY vocabulary_id) AS vocab
            FROM (
                SELECT DISTINCT person_id, vocabulary_id
                FROM all_procedure_events
            )
            GROUP BY person_id
        ),

        -- Count distinct dates per person
        date_counts AS (
            SELECT
                person_id,
                COUNT(DISTINCT event_date) AS {{ name }}_n
            FROM all_procedure_events
            GROUP BY person_id
        ){% if dates %},

        -- Get ordered dates per person
        ordered_dates AS (
            SELECT
                person_id,
                event_date,
                ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY event_date ASC) AS date_order,
                ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY event_date DESC) AS rev_date_order
            FROM (
                SELECT DISTINCT person_id, event_date
                FROM all_procedure_events
            )
        ),

        -- Get first, second, and last dates
        key_dates AS (
            SELECT
                person_id,
                MAX(CASE WHEN date_order = 1 THEN event_date END) AS {{ name }}_1,
                MAX(CASE WHEN date_order = 2 THEN event_date END) AS {{ name }}_2,
                MAX(CASE WHEN rev_date_order = 1 THEN event_date END) AS {{ name }}_last
            FROM ordered_dates
            GROUP BY person_id
        ){% endif %}

        -- Final person-level output
        SELECT
            dc.person_id,
            dc.{{ name }}_n,
            vs.vocab{% if dates %},
            kd.{{ name }}_1,
            kd.{{ name }}_2,
            kd.{{ name }}_last{% endif %}
        FROM date_counts dc
        JOIN vocab_summary vs ON dc.person_id = vs.person_id
        {% if dates %}JOIN key_dates kd ON dc.person_id = kd.person_id{% endif %}
        ORDER BY dc.person_id
        """)       

        # Format person_ids if it's a list
        if isinstance(person_ids, list):
            person_ids = ', '.join(str(pid) for pid in person_ids)

        # Initialize dictionaries
        exact_codes = exact_codes or {}
        pattern_codes = pattern_codes or {}
        exclude_codes = exclude_codes or {}

        # Render the template
        self.query_text = template.render(
            version=self.version,
            name=name,
            search_tables=search_tables,
            vocabularies=vocabularies,
            concept_classes=concept_classes,
            search_terms=search_terms,
            exclude_terms=exclude_terms,
            exact_codes=exact_codes,
            pattern_codes=pattern_codes,
            exclude_codes=exclude_codes,
            person_ids=person_ids,
            dates=dates
        )

        return self
 
    def count_participants_with_data(self,
                                     include_icd_codes: bool = True,
                                     include_snomed_codes: bool = True,
                                     include_loinc: bool = True,
                                     include_drugs: bool = True,
                                     include_procedures: bool = False,
                                     measurement_registration_exclusion: bool = True,
                                     custom_conditions: Dict[str, str] = None):
        """
        Count unique persons with data in specified domains.

        Args:
            include_icd_codes: Include ICD9/ICD10 codes from condition_occurrence and observation source fields
            include_snomed_codes: Include SNOMED codes from condition_occurrence and observation source fields
            include_loinc: Include LOINC codes from measurement
            include_drugs: Include drug exposures with domain_id = "Drug"
            include_procedures: Include procedures from procedure_occurrence and observation
                           Uses default vocabularies: CPT4, HCPCS, ICD10PCS, ICD9Proc, LOINC, SNOMED
            measurement_registration_exclusion: Exclude All of Us registration vitals from measurement count 
                (i.e. look for EHR measurements, not registration measurements)
            custom_conditions: Dictionary of {table_name: WHERE clause} for custom filtering

        Returns:
            The query object, call execute_gbq() to run the query and get the count
        """
        # Default vitals exclusions
        default_vitals_exclusion = [3022318, 3027018, 3031203, 40759207, 40765148, 3036277, 
                                    3025315, 3012888, 3004249, 3038553, 3022281]
        
        # Set string for measurement exclusion
        if measurement_registration_exclusion:
            exclude_concept_ids_str = ', '.join(map(str, default_vitals_exclusion))
            excluded_vitals = default_vitals_exclusion
        else:
            exclude_concept_ids_str = ""
            excluded_vitals = []

        custom_conditions = custom_conditions or {}

        # Determine which vocabularies to include
        condition_list = []
        if include_icd_codes:
            condition_list.extend(['ICD9CM', 'ICD10CM'])
        if include_snomed_codes:
            condition_list.append('SNOMED')

        # Default procedure vocabularies
        procedure_list = ["CPT4", "HCPCS", "ICD10PCS", "ICD9Proc", "LOINC", "SNOMED"] if include_procedures else []

        # Store query parameters for summary
        self.last_query_params = {
            "type": "person_count",
            "include_icd_codes": include_icd_codes,
            "include_snomed_codes": include_snomed_codes,
            "include_loinc": include_loinc,
            "include_drugs": include_drugs,
            "include_procedures": include_procedures,
            "procedure_vocabularies": procedure_list,
            "excluded_vitals": excluded_vitals
        }

        # Generate query summary text
        self.query_summary = self._generate_query_summary()

        # Create Jinja2 environment
        env = Environment()

        # Define and add filter to environment
        def join_quotes(items):
            return ', '.join(f"'{item}'" for item in items)

        env.filters['join_quotes'] = join_quotes
       
        template = env.from_string("""
        WITH combined AS (
            {% if condition_list %}
            -- Codes from observation (source_value)
            SELECT o.person_id
            FROM {{ version }}.observation AS o
            JOIN {{ version }}.concept AS c ON o.observation_source_value = c.concept_code
            WHERE c.vocabulary_id IN ({{ condition_list|join_quotes }})

            UNION ALL

            -- Codes from observation (source_concept_id)
            SELECT o.person_id
            FROM {{ version }}.observation AS o
            JOIN {{ version }}.concept AS c ON o.observation_source_concept_id = c.concept_id
            WHERE c.vocabulary_id IN ({{ condition_list|join_quotes }})

            UNION ALL

            -- Codes from condition_occurrence (source_value)
            SELECT co.person_id
            FROM {{ version }}.condition_occurrence AS co
            JOIN {{ version }}.concept AS c ON co.condition_source_value = c.concept_code
            WHERE c.vocabulary_id IN ({{ condition_list|join_quotes }})

            UNION ALL

            -- Codes from condition_occurrence (source_concept_id)
            SELECT co.person_id
            FROM {{ version }}.condition_occurrence AS co
            JOIN {{ version }}.concept AS c ON co.condition_source_concept_id = c.concept_id
            WHERE c.vocabulary_id IN ({{ condition_list|join_quotes }})
            {% endif %}

            {% if include_loinc %}
            {% if condition_list %}UNION ALL{% endif %}
            -- LOINC codes from measurement
            SELECT m.person_id
            FROM {{ version }}.measurement AS m
            JOIN {{ version }}.concept AS c ON m.measurement_concept_id = c.concept_id
            WHERE c.vocabulary_id = 'LOINC'
            {% if measurement_registration_exclusion %}
                AND c.concept_id NOT IN ({{ exclude_concept_ids_str }})
            {% endif %}
            {% endif %}
    
            {% if include_drugs %}
            {% if condition_list or include_loinc %}UNION ALL{% endif %}
            -- Drug exposures
            SELECT de.person_id
            FROM {{ version }}.drug_exposure AS de
            JOIN {{ version }}.concept AS c ON de.drug_concept_id = c.concept_id
            WHERE c.domain_id = 'Drug'
            {% endif %}
            
            {% if include_procedures %}
            {% if condition_list or include_loinc or include_drugs %}UNION ALL{% endif %}
            -- Procedures from procedure_occurrence (source_value)
            SELECT po.person_id
            FROM {{ version }}.procedure_occurrence AS po
            JOIN {{ version }}.concept AS c ON po.procedure_source_value = c.concept_code
            WHERE c.vocabulary_id IN ({{ procedure_list|join_quotes }})

            UNION ALL

            -- Procedures from procedure_occurrence (source_concept_id)
            SELECT po.person_id
            FROM {{ version }}.procedure_occurrence AS po
            JOIN {{ version }}.concept AS c ON po.procedure_source_concept_id = c.concept_id
            WHERE c.vocabulary_id IN ({{ procedure_list|join_quotes }})

            UNION ALL

            -- Procedures from observation (source_value)
            SELECT o.person_id
            FROM {{ version }}.observation AS o
            JOIN {{ version }}.concept AS c ON o.observation_source_value = c.concept_code
            WHERE c.vocabulary_id IN ({{ procedure_list|join_quotes }})

            UNION ALL

            -- Procedures from observation (source_concept_id)
            SELECT o.person_id
            FROM {{ version }}.observation AS o
            JOIN {{ version }}.concept AS c ON o.observation_source_concept_id = c.concept_id
            WHERE c.vocabulary_id IN ({{ procedure_list|join_quotes }})
            {% endif %}

        )

        SELECT COUNT(DISTINCT person_id) AS person_count
        FROM combined
        """)

        # Render the template
        self.query_text = template.render(
            version=self.version,
            condition_list=condition_list,
            procedure_list=procedure_list,
            include_icd_codes=include_icd_codes,
            include_snomed_codes=include_snomed_codes,
            include_loinc=include_loinc,
            include_drugs=include_drugs,
            include_procedures=include_procedures,
            measurement_registration_exclusion=measurement_registration_exclusion,
            exclude_concept_ids_str=exclude_concept_ids_str,
            custom_conditions=custom_conditions
        )
        
        return self

    def print_query(self):
        """Print the current query for inspection"""
        if not self.query_text:
            print("No query built yet")
            return self
            
        print(self.query_text)
        return self

    def query_to_variable(self, var_name):
        """Set a variable to the current query text"""
        if not self.query_text:
            raise ValueError("No query has been built yet")
        
        globals()[var_name] = self.query_text
        print(f"Query saved to variable '{var_name}'")
        return self

    def _generate_query_summary(self):
        """Generate a human-readable summary of the last query parameters"""
        if not hasattr(self, 'last_query_params'):
            return "No query parameters available"

        params = self.last_query_params
        query_type = params.get("type", "unknown")

        if query_type == "diagnosis_codes":
            summary = []

            # Describe data sources
            summary.append(f"Query Type: Diagnosis Codes")

            # Describe vocabularies
            vocabs = params.get("vocabularies")
            if vocabs:
                summary.append(f"Vocabularies: {', '.join(vocabs)}")

            # Describe tables
            summary.append("Tables: condition_occurrence, observation")
            summary.append("Fields: Using both source_value and source_concept_id paths")

            # Describe search criteria
            search_terms = params.get("search_terms")
            if search_terms:
                terms_with_quotes = [f'"{term}"' for term in search_terms]
                terms_joined = ', '.join(terms_with_quotes)
                summary.append(f"Search Terms: {terms_joined}")

            exclude_terms = params.get("exclude_terms")
            if exclude_terms:
                terms_with_quotes = [f'"{term}"' for term in exclude_terms]
                terms_joined = ', '.join(terms_with_quotes)
                summary.append(f"Exclude Terms: {terms_joined}")

            # Describe code criteria
            exact_codes = params.get("exact_codes")
            if exact_codes:
                code_lists = []
                for vocab, codes in exact_codes.items():
                    if codes:
                        sorted_codes = sorted(codes)  
                        code_lists.append(f"\n  {vocab}: {', '.join(sorted_codes)}")
                if code_lists:
                    summary.append(f"Exact Codes: {' '.join(code_lists)}")

            pattern_codes = params.get("pattern_codes")
            if pattern_codes:
                pattern_lists = []
                for vocab, patterns in pattern_codes.items():
                    if patterns:
                        sorted_patterns = sorted(patterns)
                        pattern_lists.append(f"\n  {vocab}: {', '.join(sorted_patterns)}")
                if pattern_lists:
                    summary.append(f"Pattern Codes: {' '.join(pattern_lists)}")

            exclude_codes = params.get("exclude_codes")
            if exclude_codes:
                exclude_lists = []
                for vocab, codes in exclude_codes.items():
                    if codes:
                        sorted_codes = sorted(codes) 
                        exclude_lists.append(f"\n  {vocab}: {', '.join(sorted_codes)}")
                if exclude_lists:
                    summary.append(f"Exclude Codes: {' '.join(exclude_lists)}")

            # Describe person filtering
            person_ids = params.get("person_ids")
            if person_ids:
                if isinstance(person_ids, list):
                    summary.append(f"Filtered to {len(person_ids)} specific person IDs")
                else:
                    summary.append(f"Filtered to specific person IDs")

            # V code handling
            summary.append("Special handling for V codes to ensure correct vocabulary attribution")

            return "\n".join(summary)
        
        elif query_type == "person_diagnosis":

            summary = []

            # Describe data sources
            summary.append(f"Query Type: Person-Level Diagnosis Data")
            summary.append(f"Column Name Prefix: {params.get('name', 'unknown')}")

            # Describe vocabularies
            vocabs = params.get("vocabularies")
            if vocabs:
                summary.append(f"Vocabularies: {', '.join(vocabs)}")

            # Describe tables
            summary.append("Tables: condition_occurrence, observation")
            summary.append("Fields: Using both source_value and source_concept_id paths")

            # Describe search criteria
            search_terms = params.get("search_terms")
            if search_terms:
                terms_with_quotes = [f'"{term}"' for term in search_terms]
                terms_joined = ', '.join(terms_with_quotes)
                summary.append(f"Search Terms: {terms_joined}")

            exclude_terms = params.get("exclude_terms")
            if exclude_terms:
                terms_with_quotes = [f'"{term}"' for term in exclude_terms]
                terms_joined = ', '.join(terms_with_quotes)
                summary.append(f"Exclude Terms: {terms_joined}")

            # Describe code criteria
            exact_codes = params.get("exact_codes")
            if exact_codes:
                code_lists = []
                for vocab, codes in exact_codes.items():
                    if codes:
                        sorted_codes = sorted(codes)  
                        code_lists.append(f"\n  {vocab}: {', '.join(sorted_codes)}")
                if code_lists:
                    summary.append(f"Exact Codes: {' '.join(code_lists)}")

            pattern_codes = params.get("pattern_codes")
            if pattern_codes:
                pattern_lists = []
                for vocab, patterns in pattern_codes.items():
                    if patterns:
                        sorted_patterns = sorted(patterns)
                        pattern_lists.append(f"\n  {vocab}: {', '.join(sorted_patterns)}")
                if pattern_lists:
                    summary.append(f"Pattern Codes: {' '.join(pattern_lists)}")

            exclude_codes = params.get("exclude_codes")
            if exclude_codes:
                exclude_lists = []
                for vocab, codes in exclude_codes.items():
                    if codes:
                        sorted_codes = sorted(codes) 
                        exclude_lists.append(f"\n  {vocab}: {', '.join(sorted_codes)}")
                if exclude_lists:
                    summary.append(f"Exclude Codes: {' '.join(exclude_lists)}")

            # Describe person filtering
            person_ids = params.get("person_ids")
            if person_ids:
                if isinstance(person_ids, list):
                    summary.append(f"Filtered to {len(person_ids)} specific person IDs")
                else:
                    summary.append(f"Filtered to specific person IDs")

            # Date columns
            dates = params.get("dates")
            if dates:
                summary.append(f"Including date columns: first, second, and last occurrences")

            # V code handling
            summary.append("Special handling for V codes to ensure correct vocabulary attribution")

            return "\n".join(summary)

        elif query_type == "procedure_codes":
            summary = []

            # Describe data sources
            summary.append(f"Query Type: Procedure Codes")

            # Describe search tables
            search_tables = params.get("search_tables")
            if search_tables:
                summary.append(f"Tables: {', '.join(search_tables)}")

            summary.append("Fields: Using both source_value and source_concept_id paths")

            # Describe vocabularies
            vocabs = params.get("vocabularies")
            if vocabs:
                summary.append(f"Vocabularies: {', '.join(vocabs)}")

            # Describe concept classes
            concept_classes = params.get("concept_classes")
            if concept_classes:
                summary.append(f"Concept Classes: {', '.join(concept_classes)}")

            # Describe search criteria
            search_terms = params.get("search_terms")
            if search_terms:
                terms_with_quotes = [f'"{term}"' for term in search_terms]
                terms_joined = ', '.join(terms_with_quotes)
                summary.append(f"Search Terms: {terms_joined}")

            exclude_terms = params.get("exclude_terms")
            if exclude_terms:
                terms_with_quotes = [f'"{term}"' for term in exclude_terms]
                terms_joined = ', '.join(terms_with_quotes)
                summary.append(f"Exclude Terms: {terms_joined}")

            # Describe code criteria
            exact_codes = params.get("exact_codes")
            if exact_codes:
                code_lists = []
                for vocab, codes in exact_codes.items():
                    if codes:
                        sorted_codes = sorted(codes)  
                        code_lists.append(f"\n  {vocab}: {', '.join(sorted_codes)}")
                if code_lists:
                    summary.append(f"Exact Codes: {' '.join(code_lists)}")

            pattern_codes = params.get("pattern_codes")
            if pattern_codes:
                pattern_lists = []
                for vocab, patterns in pattern_codes.items():
                    if patterns:
                        sorted_patterns = sorted(patterns)
                        pattern_lists.append(f"\n  {vocab}: {', '.join(sorted_patterns)}")
                if pattern_lists:
                    summary.append(f"Pattern Codes: {' '.join(pattern_lists)}")

            exclude_codes = params.get("exclude_codes")
            if exclude_codes:
                exclude_lists = []
                for vocab, codes in exclude_codes.items():
                    if codes:
                        sorted_codes = sorted(codes) 
                        exclude_lists.append(f"\n  {vocab}: {', '.join(sorted_codes)}")
                if exclude_lists:
                    summary.append(f"Exclude Codes: {' '.join(exclude_lists)}")

            # Describe person filtering
            person_ids = params.get("person_ids")
            if person_ids:
                if isinstance(person_ids, list):
                    summary.append(f"Filtered to {len(person_ids)} specific person IDs")
                else:
                    summary.append(f"Filtered to specific person IDs")

            return "\n".join(summary)

        elif query_type == "person_procedure":
            summary = []

            # Describe data sources
            summary.append(f"Query Type: Person-Level Procedure Data")
            summary.append(f"Column Name Prefix: {params.get('name', 'unknown')}")

            # Describe search tables
            search_tables = params.get("search_tables")
            if search_tables:
                summary.append(f"Tables: {', '.join(search_tables)}")

            summary.append("Fields: Using both source_value and source_concept_id paths")

            # Describe vocabularies
            vocabs = params.get("vocabularies")
            if vocabs:
                summary.append(f"Vocabularies: {', '.join(vocabs)}")

            # Describe concept classes
            concept_classes = params.get("concept_classes")
            if concept_classes:
                summary.append(f"Concept Classes: {', '.join(concept_classes)}")

            # Describe search criteria
            search_terms = params.get("search_terms")
            if search_terms:
                terms_with_quotes = [f'"{term}"' for term in search_terms]
                terms_joined = ', '.join(terms_with_quotes)
                summary.append(f"Search Terms: {terms_joined}")

            exclude_terms = params.get("exclude_terms")
            if exclude_terms:
                terms_with_quotes = [f'"{term}"' for term in exclude_terms]
                terms_joined = ', '.join(terms_with_quotes)
                summary.append(f"Exclude Terms: {terms_joined}")

            # Describe code criteria
            exact_codes = params.get("exact_codes")
            if exact_codes:
                code_lists = []
                for vocab, codes in exact_codes.items():
                    if codes:
                        sorted_codes = sorted(codes)  
                        code_lists.append(f"\n  {vocab}: {', '.join(sorted_codes)}")
                if code_lists:
                    summary.append(f"Exact Codes: {' '.join(code_lists)}")

            pattern_codes = params.get("pattern_codes")
            if pattern_codes:
                pattern_lists = []
                for vocab, patterns in pattern_codes.items():
                    if patterns:
                        sorted_patterns = sorted(patterns)
                        pattern_lists.append(f"\n  {vocab}: {', '.join(sorted_patterns)}")
                if pattern_lists:
                    summary.append(f"Pattern Codes: {' '.join(pattern_lists)}")

            exclude_codes = params.get("exclude_codes")
            if exclude_codes:
                exclude_lists = []
                for vocab, codes in exclude_codes.items():
                    if codes:
                        sorted_codes = sorted(codes) 
                        exclude_lists.append(f"\n  {vocab}: {', '.join(sorted_codes)}")
                if exclude_lists:
                    summary.append(f"Exclude Codes: {' '.join(exclude_lists)}")

            # Describe person filtering
            person_ids = params.get("person_ids")
            if person_ids:
                if isinstance(person_ids, list):
                    summary.append(f"Filtered to {len(person_ids)} specific person IDs")
                else:
                    summary.append(f"Filtered to specific person IDs")

            # Date columns
            dates = params.get("dates")
            if dates:
                summary.append(f"Including date columns: first, second, and last occurrences")

            return "\n".join(summary)

        elif query_type == "person_count":
            
            summary = []

            # Describe query type
            summary.append("Query Type: Participant Count")

            # Describe included data types
            data_types = []
            if params.get("include_icd_codes"):
                data_types.append("ICD codes (ICD9CM, ICD10CM)")
            if params.get("include_snomed_codes"):
                data_types.append("SNOMED codes")
            if params.get("include_loinc"):
                data_types.append("LOINC measurements")
            if params.get("include_drugs"):
                data_types.append("Drug exposures")
            if params.get("include_procedures"):
                procedure_vocabs = params.get("procedure_vocabularies", [])
                if procedure_vocabs:
                    summary.append(f"Procedure codes ({', '.join(procedure_vocabs)})")
                else:
                    data_types.append("Procedure codes")

            if data_types:
                summary.append(f"Included Data: {', '.join(data_types)}")

            # Describe procedure tables if procedures included
            if params.get("include_procedures"):
                summary.append("Procedure Tables: procedure_occurrence, observation")

            # Describe excluded concepts
            excluded_vitals = params.get("excluded_vitals", [])

            if params.get("include_loinc") and excluded_vitals:
                summary.append(f"{len(excluded_vitals)} common vitals concept IDs (from All of Us registration) excluded from measurement")

            return "\n".join(summary)
            
        else:
            return "Unknown query type"

    def print_summary(self):
        """Print a summary of the last query"""
        if hasattr(self, 'query_summary'):
            print(self.query_summary)
        else:
            print("No query summary available")
        return self

    def execute_gbq(self, quiet=False):
        """Execute the query and return a Polars DataFrame"""
        if not self.query_text:
            raise ValueError("No query has been built yet")
                
        if not quiet and hasattr(self, 'query_summary'):
            print(self.query_summary)

        # Execute query
        query_job = self.client.query(self.query_text)

        # Get results and convert to polars
        rows = query_job.result()
        df = pl.from_arrow(rows.to_arrow())

        return df

    def results_to_code_dict(self, df):
        """
        Convert a diagnosis codes dataframe to an exact_codes dictionary for use in find_diagnosis_codes.

        Args:
            df: Polars DataFrame with vocabulary_id and concept_code columns

        Returns:
            Dict mapping vocabularies to lists of concept codes
        """
        if not {'vocabulary_id', 'concept_code'}.issubset(df.columns):
            raise ValueError("DataFrame must contain 'vocabulary_id' and 'concept_code' columns")

        result = {}

        # Group by vocabulary_id and collect concept_codes
        for vocab in df['vocabulary_id'].unique():
            codes = df.filter(pl.col('vocabulary_id') == vocab)['concept_code'].to_list()
            result[vocab] = codes

        return result


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


# In[ ]:


def create_icd_dict_from_phecodes(phecodex_map, phecodes):
    # Convert single phecode to list if needed
    if isinstance(phecodes, str):
        phecodes = [phecodes]
    
    # Filter the DataFrame to only include rows with phecodes in the list
    filtered_df = phecodex_map.filter(pl.col("phecode").is_in(phecodes))
    
    # Group by vocabulary_id and aggregate unique ICD codes into lists
    result = (filtered_df
              .group_by("vocabulary_id")
              .agg(pl.col("ICD").unique().alias("icd_codes"))
              .to_dict(as_series=False))
    
    # Convert to the desired dictionary format
    icd_dict = {vocab: codes for vocab, codes in zip(result["vocabulary_id"], result["icd_codes"])}
    
    return icd_dict


# In[ ]:


phecodex_map = pl.read_csv('phecodeX_unrolled_ICD_CM.csv')


# # Initialize

# ## AouQueries

# In[ ]:


# Initialize
aou = AouQueries()


# In[ ]:


version = 'fc-aou-cdr-prod-ct.C2024Q3R5'
bucket = '{bucket or my_bucket}'


# ## Unifier

# In[ ]:


# Environment Variables
# bucket = os.getenv("WORKSPACE_BUCKET")
dataset = os.getenv('WORKSPACE_CDR')

# Global Variables
cohort = "allofus"
version = dataset.split('.')[1]
file_type = "parquet" # csv also possible
file_path = f"{bucket}/data/{cohort}/{version}"

prostate_dict = create_icd_dict_from_phecodes(phecodex_map, ["CA_107.2"])# Find sarcoid diagnoses
prostate_codes_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    exact_codes=prostate_dict, 
    exclude_terms=["family"]
).execute_gbq())prostate_codes_dfprostate_df = (aou.person_code_df(
    name="sarcoid",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    exact_codes=prostate_dict, 
    exclude_terms=["family"],
    dates=True
).execute_gbq())prostate_df = prostate_df.join(
    demographics_df.select('person_id', 'dob'),
    on='person_id',
    how='left'
)prostate_df = prostate_df.with_columns(
    ((pl.col("sarcoid_1") - pl.col("dob")).dt.total_days() / 365.25).floor().alias("age_at_sarcoid_dx")
)sns.histplot(prostate_df.filter(pl.col('sarcoid_n')>1)['age_at_sarcoid_dx'])prostate_df.filter((pl.col('sarcoid_n')>1) & (pl.col('age_at_sarcoid_dx')<50))
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
    {dataset}.person p
LEFT JOIN
    {dataset}.concept p_gender_concept 
        ON p.gender_concept_id = p_gender_concept.concept_id 
LEFT JOIN
    {dataset}.concept p_race_concept 
        ON p.race_concept_id = p_race_concept.concept_id 
LEFT JOIN
    {dataset}.concept p_ethnicity_concept 
        ON p.ethnicity_concept_id = p_ethnicity_concept.concept_id 
LEFT JOIN
    {dataset}.concept p_sex_at_birth_concept 
        ON p.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id  
LEFT JOIN
    {dataset}.death d
        ON p.person_id = d.person_id
WHERE
    p.person_id IN (SELECT person_id FROM distinct_participants)
"""

demographics_df = polars_gbq(demographics_q)


# # Identify Sarcoid Codes 

# In[ ]:


# Find sarcoid diagnoses
sarcoid_codes_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["sarcoid"], 
    exclude_terms=["history"]
).print_query())


# In[ ]:


# Find sarcoid diagnoses
sarcoid_codes_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["sarcoid"], 
    exclude_terms=["family"]
).execute_gbq())


# In[ ]:


sarcoid_codes_df


# In[ ]:


sarcoid_df = (aou.person_code_df(
    name="sarcoid",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["sarcoid"], 
    exclude_terms=["family"],
    dates=True
).execute_gbq())


# In[ ]:


print(f"Participants with 1 code: {sarcoid_df.filter(pl.col('sarcoid_n')==1).height}")
print(f"Participants with ≥ 2 codes: {sarcoid_df.filter(pl.col('sarcoid_n')>1).height}")


# In[ ]:


(aou.count_participants_with_data(
    include_icd_codes=True,
    include_snomed_codes=True,
    include_loinc=False,
    include_drugs=False,
).execute_gbq())


# In[ ]:


(aou.count_participants_with_data(
    include_icd_codes=False,
    include_snomed_codes=False,
    include_loinc=False,
    include_drugs=False,
    include_procedures=True
).execute_gbq())


# In[ ]:


sarcoid_df.head()


# # Identify Exclusions

# ## TB
# This excludes nonspecific reaction to TB skin test codes.

# In[ ]:


# Find TB diagnoses by string
tb_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["tuberculous", "tuberculosis", "tuberculoma"],
    exclude_terms=['encounter for screening', 'screening status', 'testing', 'contact', 
                   'without mention of effusion or current tuberculosis', 'other than tuberculosis',
                   'exposure to', 'bcg', 'vaccine', 'pleural effusion, except tuberculous',
                   'nonspecific reaction to'
                  ]
).execute_gbq())


# In[ ]:


tb_dict = create_icd_dict_from_phecodes(phecodex_map, ["ID_005.1", "RE_460.3"])


# In[ ]:


# Find sarcoid diagnoses by phecodes
tb_codes_phecode_df = (aou.find_diagnosis_codes(
    exact_codes=tb_dict
).execute_gbq())


# In[ ]:


# Create a combined dataframe with source labels
phecode_df = tb_codes_phecode_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("phecode").alias("source"))
string_df = tb_codes_string_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("string").alias("source"))

# Find rows in phecode_df not in string_df
unique_to_phecode = phecode_df.join(
    string_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Find rows in string_df not in phecode_df
unique_to_string = string_df.join(
    phecode_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Concatenate the unique rows
unique_rows_df = pl.concat([unique_to_phecode, unique_to_string])


# In[ ]:


unique_rows_df.head()


# In[ ]:


tb_code_df = (aou.person_code_df(
    name="tb",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["tuberculous", "tuberculosis", "tuberculoma"],
    exclude_terms=['encounter for screening', 'screening status', 'testing', 'contact', 
                   'without mention of effusion or current tuberculosis', 'other than tuberculosis',
                   'exposure to', 'bcg', 'vaccine', 'pleural effusion, except tuberculous',
                   'nonspecific reaction to'
                  ],
    dates=True
).execute_gbq())


# ## Abnormal TB Skin Test
# "ICD10CM"	"R76.11"	"Nonspecific reaction to tuberculin skin test without active tuberculosis"	3900	36006<br>
# "ICD10CM"	"R76.12"	"Nonspecific reaction to cell mediated immunity measurement of gamma interferon antigen response with…	1215	10858<br>
# "ICD9CM"	"795.5"	"Nonspecific reaction to tuberculin skin test without active tuberculosis"	2753	20920<br>
# "ICD9CM"	"795.51"	"Nonspecific reaction to tuberculin skin test without active tuberculosis"	2817	25598<br>
# "ICD9CM"	"795.52"	"Nonspecific reaction to cell mediated immunity measurement of gamma interferon antigen response with…	426	4204<br>

# In[ ]:


tb_skin_dict = {
    "ICD10CM": ["R76.1%"],
    "ICD9CM": ["795.5%"],
}

# Find TB skin test diagnoses
tb_skin_codes_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    pattern_codes=tb_skin_dict
).execute_gbq())


# In[ ]:


tb_skin_codes_df


# In[ ]:


tb_skin_code_df = (aou.person_code_df(
    name="tb_skin",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    pattern_codes=tb_skin_dict,
    dates=True
).execute_gbq())


# ## Other Mycobacterial Infections

# In[ ]:


ntm_dict = {
    "ICD10CM": ["A30%", "A31%"],
    "ICD9CM": ["V74.2", "030%", "031%"],
}

# Find NTM diagnoses by string
ntm_codes_icd_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    pattern_codes=ntm_dict,
).execute_gbq())


# In[ ]:


ntm_codes_icd_df


# In[ ]:


# Find NTM diagnoses by string
ntm_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["nontuberculous", "non-tuberculous", "avium-intracellulare", "abscessus", "chelonae", 
                  "fortuitum", "gordonae", "malmoense", "marinum", "ulcerans", "xenopi", "kansasii",
                 "nonchromogenicum", "avium complex", "genavense", "haemophilum", "leprosy", "mycobacteria"],
).execute_gbq())


# In[ ]:


ntm_code_df = (aou.person_code_df(
    name="ntm",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["nontuberculous", "non-tuberculous", "avium-intracellulare", "abscessus", "chelonae", 
                  "fortuitum", "gordonae", "malmoense", "marinum", "ulcerans", "xenopi", "kansasii",
                 "nonchromogenicum", "avium complex", "genavense", "haemophilum", "leprosy", "mycobacteria"],
    dates=True
).execute_gbq())


# ## Histo, Blasto, Cocci

# In[ ]:


dimorphic_dict = create_icd_dict_from_phecodes(phecodex_map, ["ID_071", "ID_072", "ID_073"])


# In[ ]:


# Find sarcoid diagnoses by phecodes
dimorphic_codes_phecode_df = (aou.find_diagnosis_codes(
    exact_codes=dimorphic_dict
).execute_gbq())


# In[ ]:


dimorphic_dict = {
    "ICD9CM": ["114%", "115%", "116%"],
    "ICD10CM": ["B38%", "B39%", "B40"],
}

# Find dimorphics diagnoses by code
dimorphic_codes_icd_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    pattern_codes=dimorphic_dict,
).execute_gbq())


# In[ ]:


# Find dimorphic diagnoses by string
dimorphic_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["histoplas", "blastomy", "coccidio"]
).execute_gbq())


# In[ ]:


# Create a combined dataframe with source labels
phecode_df = dimorphic_codes_phecode_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("phecode").alias("source"))
icd_df = dimorphic_codes_icd_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("phecode").alias("source"))
string_df = dimorphic_codes_string_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("string").alias("source"))

# Find rows in phecode_df not in string_df
unique_to_phecode = phecode_df.join(
    string_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Find rows in phecode_df not in string_df
unique_to_icd = icd_df.join(
    string_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

unique_rows_df = pl.concat([unique_to_phecode, unique_to_icd])


# In[ ]:


unique_rows_df


# In[ ]:


dimorphic_codes_string_df.slice(50,50)


# In[ ]:


histo_code_df = (aou.person_code_df(
    name="histo",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["histoplas"],
    dates=True
).execute_gbq())


# In[ ]:


blasto_code_df = (aou.person_code_df(
    name="blasto",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["blastomy"],
    dates=True
).execute_gbq())


# In[ ]:


cocci_code_df = (aou.person_code_df(
    name="cocci",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["coccidio"],
    dates=True
).execute_gbq())


# ## Aspergillus

# In[ ]:


asp_dict = create_icd_dict_from_phecodes(phecodex_map, ["ID_074"]) #RE_475.7 is redundant


# In[ ]:


asp_codes_phecode_df = (aou.find_diagnosis_codes(
    exact_codes=asp_dict
).execute_gbq())


# In[ ]:


aspergillus_dict = {
    "ICD10CM": ["B44%"],
    "ICD9CM": ["117.3%"],
}
# All included in above


# In[ ]:


asp_codes_phecode_df


# In[ ]:


# Find Aspergillus diagnoses by string
aspergillus_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["aspergill"],
).execute_gbq())


# In[ ]:


aspergillus_codes_string_df


# In[ ]:


aspergillus_code_df = (aou.person_code_df(
    name="asp",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["aspergill"],
    dates=True
).execute_gbq())


# ## ILD (incl. Eosinophilia Pneumoina)

# In[ ]:


ild_dict = create_icd_dict_from_phecodes(phecodex_map, ["RE_481"]) 


# In[ ]:


ild_codes_phecode_df = (aou.find_diagnosis_codes(
    exact_codes=ild_dict
).execute_gbq())


# In[ ]:


ild_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["interstitial pneumonia", "eosinophilic pneumonia", "eosinophilic asthma",
                 "pulmonary eosinophilia", "alveolar proteinosis", "alveolar microlithiasis",
                 "pulmonary hemosiderosis", "alveolar and parieto-alveolar", "interstitial pulmonary",
                 "interstitial lung", "pulmonary fibrosis", "interstitial pneumonitis", 
                 "organizing pneumonia", "lymphangioleiomyomatosis", "parietoalveolar pneumono", 
                  "neuroendocrine cell hyperplasia of infancy", "pulmonary interstitial", 
                  "alveolar and parieto-alveolar"
                 ],
).execute_gbq())


# In[ ]:


# Create a combined dataframe with source labels
phecode_df = ild_codes_phecode_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("phecode").alias("source"))
string_df = ild_codes_string_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("string").alias("source"))

# Find rows in phecode_df not in string_df
unique_to_phecode = phecode_df.join(
    string_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Find rows in string_df not in phecode_df
unique_to_string = string_df.join(
    phecode_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Concatenate the unique rows
unique_rows_df = pl.concat([unique_to_phecode, unique_to_string])


# In[ ]:


if unique_to_phecode.height == 0:
    print('No codes unique to phecode')
else:
    unique_to_phecode


# In[ ]:


ild_code_df = (aou.person_code_df(
    name="ild",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["interstitial pneumonia", "eosinophilic pneumonia", "eosinophilic asthma",
                 "pulmonary eosinophilia", "alveolar proteinosis", "alveolar microlithiasis",
                 "pulmonary hemosiderosis", "alveolar and parieto-alveolar", "interstitial pulmonary",
                 "interstitial lung", "pulmonary fibrosis", "interstitial pneumonitis", 
                 "organizing pneumonia", "lymphangioleiomyomatosis", "parietoalveolar pneumono", 
                  "neuroendocrine cell hyperplasia of infancy", "pulmonary interstitial", 
                  "alveolar and parieto-alveolar"
                 ],
    dates=True
).execute_gbq())


# ## Hypersensitivity Pneumonitis

# In[ ]:


hp_dict = create_icd_dict_from_phecodes(phecodex_map, ["RE_477.2"]) 


# In[ ]:


hp_codes_phecode_df = (aou.find_diagnosis_codes(
    exact_codes=hp_dict
).execute_gbq())


# In[ ]:


hp_codes_phecode_df


# In[ ]:


hp_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["s lung", "Farmers", "bagassosis", "Bird fancier", "Bird-fanciers", 
                 "Hypersensitivity pneumonitis", "Pneumonitis due to inhalation of oils", 
                 "Pneumonitis due to inhalation of other solids", "Suberosis", "Mushroom workers", 
                 "Maple bark-strippers", "Ventilation pneumonitis", "Other specified allergic alveolitis"
                 "Pneumonitis due to other solids", "Pneumonitis due to solids", "Pneumonitis due to other solids and liquids", 
                  "Extrinsic allergic alveolitis", "\"Ventilation\" pneumonitis"],
    exact_codes={'ICD9CM': ['495.8', '495.9']},
).execute_gbq())


# In[ ]:


# Create a combined dataframe with source labels
phecode_df = hp_codes_phecode_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("phecode").alias("source"))
string_df = hp_codes_string_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("string").alias("source"))

# Find rows in phecode_df not in string_df
unique_to_phecode = phecode_df.join(
    string_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Find rows in string_df not in phecode_df
unique_to_string = string_df.join(
    phecode_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Concatenate the unique rows
unique_rows_df = pl.concat([unique_to_phecode, unique_to_string])


# In[ ]:


if unique_to_phecode.height == 0:
    print('No codes unique to phecode')
else:
    display(unique_to_phecode)


# In[ ]:


hp_code_df = (aou.person_code_df(
    name="hp",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["s lung", "Farmers", "bagassosis", "Bird fancier", "Bird-fanciers", 
                 "Hypersensitivity pneumonitis", "Pneumonitis due to inhalation of oils", 
                 "Pneumonitis due to inhalation of other solids", "Suberosis", "Mushroom workers", 
                 "Maple bark-strippers", "Ventilation pneumonitis", "Other specified allergic alveolitis"
                 "Pneumonitis due to other solids", "Pneumonitis due to solids", "Pneumonitis due to other solids and liquids", 
                  "Extrinsic allergic alveolitis", "\"Ventilation\" pneumonitis"],
    exact_codes={'ICD9CM': ['495.8', '495.9']},
    dates=True
).execute_gbq())


# ## Pneumoconiosis (beryllium, titanium, aluminum, zirconium, cobalt, others)

# In[ ]:


pneumoconiosis_dict = create_icd_dict_from_phecodes(phecodex_map, ["RE_477.1"]) 


# In[ ]:


pneumoconiosis_codes_phecode_df = (aou.find_diagnosis_codes(
    exact_codes=pneumoconiosis_dict
).execute_gbq())


# In[ ]:


# Find diagnoses by string
pneumoconiosis_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["pneumoconiosis", "Asbestosis", "Stannosis", "Siderosis", "Graphite fibrosis",
                 "Bauxite fibrosis"],
).execute_gbq())


# In[ ]:


# Create a combined dataframe with source labels
phecode_df = pneumoconiosis_codes_phecode_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("phecode").alias("source"))
string_df = pneumoconiosis_codes_string_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("string").alias("source"))

# Find rows in phecode_df not in string_df
unique_to_phecode = phecode_df.join(
    string_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Find rows in string_df not in phecode_df
unique_to_string = string_df.join(
    phecode_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Concatenate the unique rows
unique_rows_df = pl.concat([unique_to_phecode, unique_to_string])


# In[ ]:


if unique_to_phecode.height == 0:
    print('No codes unique to phecode')
else:
    display(unique_to_phecode)


# In[ ]:


pneumoconiosis_codes_string_df


# In[ ]:


pneumoconiosis_code_df = (aou.person_code_df(
    name="pneumoconiosis",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["pneumoconiosis", "Asbestosis", "Stannosis", "Siderosis", "Graphite fibrosis",
                 "Bauxite fibrosis"],
    dates=True
).execute_gbq())


# ## Drugs

# In[ ]:


drug_exclusion_code_q = f"""
WITH filtered_concepts AS (
    SELECT DISTINCT concept_id, concept_code, concept_name
    FROM {dataset}.concept c
    WHERE vocabulary_id = 'RxNorm'
    AND (
    LOWER(concept_name) LIKE '%ipilimumab%'
    OR LOWER(concept_name) LIKE '%nivolumab%'
    OR LOWER(concept_name) LIKE '%pembrolizumab%'
    OR LOWER(concept_name) LIKE '%interferon-alpha%'
    OR LOWER(concept_name) LIKE '%interferon-beta%'
    OR LOWER(concept_name) LIKE '%etanercept%'
    OR LOWER(concept_name) LIKE '%adalimumab%'
    OR LOWER(concept_name) LIKE '%infliximab%'
    OR LOWER(concept_name) LIKE '%vemurafenib%'
    OR LOWER(concept_name) LIKE '%dabrafenib%'
    OR LOWER(concept_name) LIKE '%trametinib%'
    OR LOWER(concept_name) LIKE '%encorafenib%'
    )
),
combined AS (
    SELECT person_id, drug_exposure_start_date AS drug_date, 
        (SELECT concept_name FROM filtered_concepts WHERE concept_id = drug_concept_id) AS drug_name,
        (SELECT concept_code FROM filtered_concepts WHERE concept_id = drug_concept_id) AS drug_code,
    FROM {dataset}.drug_exposure
    WHERE drug_concept_id IN (SELECT concept_id FROM filtered_concepts)
),
distinct_dates AS (
    SELECT DISTINCT person_id, drug_date, drug_code, drug_name
    FROM combined
)
SELECT drug_code, drug_name, COUNT(*) AS count
FROM distinct_dates
GROUP BY drug_code, drug_name
ORDER BY count DESC
"""

exclusion_drug_code_df = polars_gbq(drug_exclusion_code_q)


# In[ ]:


drug_exclusion_q = f"""
WITH filtered_concepts AS (
    SELECT DISTINCT concept_id, concept_code, concept_name
    FROM {dataset}.concept c
    WHERE vocabulary_id = 'RxNorm'
    AND (
    LOWER(concept_name) LIKE '%ipilimumab%'
    OR LOWER(concept_name) LIKE '%nivolumab%'
    OR LOWER(concept_name) LIKE '%pembrolizumab%'
    OR LOWER(concept_name) LIKE '%interferon-alpha%'
    OR LOWER(concept_name) LIKE '%interferon-beta%'
    OR LOWER(concept_name) LIKE '%etanercept%'
    OR LOWER(concept_name) LIKE '%adalimumab%'
    OR LOWER(concept_name) LIKE '%infliximab%'
    OR LOWER(concept_name) LIKE '%vemurafenib%'
    OR LOWER(concept_name) LIKE '%dabrafenib%'
    OR LOWER(concept_name) LIKE '%trametinib%'
    OR LOWER(concept_name) LIKE '%encorafenib%'
    )
),
combined AS (
    SELECT de.person_id, drug_exposure_start_date AS drug_date
    FROM {dataset}.drug_exposure de
    WHERE drug_concept_id IN (SELECT concept_id FROM filtered_concepts)
),
distinct_dates AS (
    SELECT DISTINCT person_id, drug_date
    FROM combined
)
SELECT 
    person_id,
    MIN(drug_date) AS exclusion_drug
FROM distinct_dates
GROUP BY person_id
"""

exclusion_drug_df = polars_gbq(drug_exclusion_q)


# In[ ]:


exclusion_drug_df.head()


# ## ANCA-Associated Vasculitis (GPA, MPA, EGPA) 

# In[ ]:


vasculitis_phecode_dict = create_icd_dict_from_phecodes(phecodex_map, ["MS_704"]) 


# In[ ]:


vasculitis_codes_phecode_df = (aou.find_diagnosis_codes(
    exact_codes=vasculitis_phecode_dict
).execute_gbq())


# In[ ]:


vasculitis_codes_phecode_df


# In[ ]:


vasculitis_dict = {
    "ICD10CM": ["I77.6", "I77.82", "M30.1", "M31.3", "M31.30", "M31.31", "M31.7", "M31.8", "M31.9"],
    "ICD9CM": ["446.21", "446.4", "446.0"],
}

# Find diagnoses by string
vasculitis_codes_icd_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    exact_codes=vasculitis_dict,
    search_terms=["Arteritis, unspecified", "unspecified arteritis", 'ANCA', 
                  "Antineutrophilic cytoplasmic", "Churg-Strauss", "Polyarteritis", 
                  "Wegener", "granulomatosis", "Microscopic polyangiitis", "Necrotizing vasculopathy",
                  "Goodpasture"],
    exclude_terms=['blanca']
).execute_gbq())


# In[ ]:


# Create a combined dataframe with source labels
phecode_df = vasculitis_codes_phecode_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("phecode").alias("source"))
string_df = vasculitis_codes_icd_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("string").alias("source"))

# Find rows in phecode_df not in string_df
unique_to_phecode = phecode_df.join(
    string_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Find rows in string_df not in phecode_df
unique_to_string = string_df.join(
    phecode_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Concatenate the unique rows
unique_rows_df = pl.concat([unique_to_phecode, unique_to_string])


# In[ ]:


if unique_to_phecode.height == 0:
    print('No codes unique to phecode')
else:
    display(unique_to_phecode)


# In[ ]:


vasculitis_code_df = (aou.person_code_df(
    name="vasc",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    exact_codes=vasculitis_dict,
    search_terms=["Arteritis, unspecified", "unspecified arteritis", 'ANCA', 
                  "Antineutrophilic cytoplasmic", "Churg-Strauss", "Polyarteritis", 
                  "Wegener", "granulomatosis", "Microscopic polyangiitis", "Necrotizing vasculopathy",
                  "Goodpasture"],
    exclude_terms=['blanca'],
    dates=True
).execute_gbq())


# ## CGD

# In[ ]:


cgd_dict = {
    "ICD10CM": ["D71%"],
}

# Find CGD diagnoses dict
cgd_codes_icd_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    pattern_codes=cgd_dict,
).execute_gbq())


# In[ ]:


cgd_codes_icd_df


# In[ ]:


# Find CGD diagnoses by string
cgd_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["cgd", "chronic granulomatous", "functional disorders of polymorphonuclear neutrophils"],
).execute_gbq())


# In[ ]:


cgd_code_df = (aou.person_code_df(
    name="cgd",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["cgd", "chronic granulomatous", "functional disorders of polymorphonuclear neutrophils"],
    dates=True
).execute_gbq())


# ## CVID

# In[ ]:


cvid_dict = create_icd_dict_from_phecodes(phecodex_map, ["BI_179.7"]) 


# In[ ]:


cvid_codes_phecode_df = (aou.find_diagnosis_codes(
    exact_codes=cvid_dict
).execute_gbq())


# In[ ]:


# Find cvid diagnoses by phecode
cvid_codes_icd_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    exact_codes=cvid_dict,
).execute_gbq())


# In[ ]:


cvid_codes_icd_df


# In[ ]:


# Find cvid diagnoses by string
cvid_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["common variable immunodefici"],
).execute_gbq())


# In[ ]:


cvid_codes_string_df


# In[ ]:


cvid_code_df = (aou.person_code_df(
    name="cvid",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["common variable immunodefici"],
    dates=True
).execute_gbq())


# ## IgG4-Related Disease

# In[ ]:


igg4_dict = {
    "ICD10CM": ["D89.84"],
}

# Find igg4 diagnoses by string
igg4_codes_icd_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    exact_codes=igg4_dict,
    search_terms=["igg4"],
).execute_gbq())


# In[ ]:


igg4_codes_icd_df


# In[ ]:


igg4_code_df = (aou.person_code_df(
    name="igg4",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["igg4"],
    dates=True
).execute_gbq())


# In[ ]:


igg4_code_df


# ## Lymphoma

# In[ ]:


lymphoma_dict = create_icd_dict_from_phecodes(phecodex_map, ["CA_122"]) 


# In[ ]:


lymphoma_codes_phecode_df = (aou.find_diagnosis_codes(
    exact_codes=lymphoma_dict
).execute_gbq())


# In[ ]:


lymphoma_codes_phecode_df.filter(~pl.col('concept_name').str.contains('lymphoma')).slice(100,50)


# In[ ]:


# Find lymphoma diagnoses by string
lymphoma_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["lymphoma", "mycosis fungoides", "sezary", "t-cell proliferation", 
                 "waldenstrom macroglobulinemia", "reticulosarcoma", "lymphosarcoma",
                 'hodgkin', "malignant neoplasms of lymphoid and histiocytic tissue", 
                  "malignant histiocytosis", "leukemic reticuloendotheliosis", "letterer-siwe disease",
                 "malignant mast cell tumors", ],
).execute_gbq())


# In[ ]:


# Create a combined dataframe with source labels
phecode_df = lymphoma_codes_phecode_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("phecode").alias("source"))
string_df = lymphoma_codes_string_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("string").alias("source"))

# Find rows in phecode_df not in string_df
unique_to_phecode = phecode_df.join(
    string_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Find rows in string_df not in phecode_df
unique_to_string = string_df.join(
    phecode_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Concatenate the unique rows
unique_rows_df = pl.concat([unique_to_phecode, unique_to_string])


# In[ ]:


if unique_to_phecode.height == 0:
    print('No codes unique to phecode')
else:
    unique_to_phecode


# In[ ]:


if unique_to_string.height == 0:
    print('No codes unique to phecode')
else:
    display(unique_to_string.slice(150,50))


# In[ ]:


lymphoma_code_df = (aou.person_code_df(
    name="lymphoma",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["lymphoma", "mycosis fungoides", "sezary", "t-cell proliferation", 
                 "waldenstrom macroglobulinemia", "reticulosarcoma", "lymphosarcoma",
                 'hodgkin', "malignant neoplasms of lymphoid and histiocytic tissue", 
                  "malignant histiocytosis", "leukemic reticuloendotheliosis", "letterer-siwe disease",
                 "malignant mast cell tumors", ],
    exclude_terms=["family history"],
    dates=True
).execute_gbq())


# ## Lymphomatoid granulomatosis
# C83 already included above (lymphoma)

# In[ ]:


# Find NTM diagnoses by string
lg_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["lymphomatoid granulomatosis"],
).execute_gbq())


# In[ ]:


lg_code_df = (aou.person_code_df(
    name="lg",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["lymphomatoid granulomatosis"],
    dates=True
).execute_gbq())


# ## Lung Cancer

# In[ ]:


lung_ca_dict = create_icd_dict_from_phecodes(phecodex_map, ["CA_102.1"]) 


# In[ ]:


lung_ca_codes_phecode_df = (aou.find_diagnosis_codes(
    exact_codes=lung_ca_dict
).execute_gbq())


# In[ ]:


lung_ca_codes_phecode_df


# In[ ]:


# Find lung_ca diagnoses by string
lung_ca_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["malignant neoplasm of bronchus", "carcinoma in situ of bronchus", 
                  "malignant neoplasm of other parts of bronchus", "malignant neoplasm of lower lobe",
                  "malignant neoplasm of middle lobe", "malignant neoplasm of upper lobe", 
                  "malignant neoplasm of main bronchus", "malignant neoplasm of trachea, bronchus", 
                  "carcinoma in situ of left bronchus", "carcinoma in situ of right bronchus", 
                  "carcinoma in situ of unspecified bronchus", "malignant neoplasm of unspecified part of left bronchus",
                 "malignant neoplasm of unspecified part of right bronchus", 
                  "malignant neoplasm of unspecified part of unspecified bronchus", 
                  "malignant neoplasm of unspecified part of bronchus", 
                  "malignant neoplasm of overlapping sites of left bronchus", 
                  "malignant neoplasm of overlapping sites of right bronchus", 
                  "malignant neoplasm of overlapping sites of unspecified bronchus",
                  "malignant neoplasm of overlapping sites of bronchus",
                  "malignant neoplasm of left main bronchus",
                  "malignant neoplasm of right main bronchus",
                  "malignant neoplasm of unspecified main bronchus", "lung cancer"],
).execute_gbq())


# In[ ]:


# Create a combined dataframe with source labels
phecode_df = lung_ca_codes_phecode_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("phecode").alias("source"))
string_df = lung_ca_codes_string_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("string").alias("source"))

# Find rows in phecode_df not in string_df
unique_to_phecode = phecode_df.join(
    string_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Find rows in string_df not in phecode_df
unique_to_string = string_df.join(
    phecode_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Concatenate the unique rows
unique_rows_df = pl.concat([unique_to_phecode, unique_to_string])


# In[ ]:


if unique_to_phecode.height == 0:
    print('No codes unique to phecode')
else:
    unique_to_phecode


# In[ ]:


if unique_to_string.height == 0:
    print('No codes unique to phecode')
else:
    display(unique_to_string.slice(0,50))


# In[ ]:


lung_ca_code_df = (aou.person_code_df(
    name="lung_ca",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["malignant neoplasm of bronchus", "carcinoma in situ of bronchus", 
                  "malignant neoplasm of other parts of bronchus", "malignant neoplasm of lower lobe",
                  "malignant neoplasm of middle lobe", "malignant neoplasm of upper lobe", 
                  "malignant neoplasm of main bronchus", "malignant neoplasm of trachea, bronchus", 
                  "carcinoma in situ of left bronchus", "carcinoma in situ of right bronchus", 
                  "carcinoma in situ of unspecified bronchus", "malignant neoplasm of unspecified part of left bronchus",
                 "malignant neoplasm of unspecified part of right bronchus", 
                  "malignant neoplasm of unspecified part of unspecified bronchus", 
                  "malignant neoplasm of unspecified part of bronchus", 
                  "malignant neoplasm of overlapping sites of left bronchus", 
                  "malignant neoplasm of overlapping sites of right bronchus", 
                  "malignant neoplasm of overlapping sites of unspecified bronchus",
                  "malignant neoplasm of overlapping sites of bronchus",
                  "malignant neoplasm of left main bronchus",
                  "malignant neoplasm of right main bronchus",
                  "malignant neoplasm of unspecified main bronchus", "lung cancer"],
    exclude_terms=["family history", "screening declined", "suspected"],
    dates=True
).execute_gbq())


# ## Germ Cell Tumor
# No specific ICD, and somewhat covered in lung cancer

# ## Langerhans Cell Histiocytosis

# In[ ]:


lch_dict = {
    "ICD10CM": ["C96.0", "C96.5", "C96.6"],
    "ICD9CM": ["516.5"],
}

# Find lch diagnoses by string
lch_codes_icd_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    exact_codes=lch_dict,
).execute_gbq())


# In[ ]:


lch_codes_icd_df


# In[ ]:


# Find lch diagnoses by string
lch_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["langerhans-cell histiocytosis", "langerhans cell histiocytosis"],
).execute_gbq())


# In[ ]:


lch_codes_string_df


# In[ ]:


lch_code_df = (aou.person_code_df(
    name="lch",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["langerhans-cell histiocytosis", "langerhans cell histiocytosis"],
    dates=True
).execute_gbq())


# ## IBD

# In[ ]:


ibd_dict = create_icd_dict_from_phecodes(phecodex_map, ["GI_522.1"]) 


# In[ ]:


ibd_codes_phecode_df = (aou.find_diagnosis_codes(
    exact_codes=ibd_dict
).execute_gbq())


# In[ ]:


ibd_codes_phecode_df.slice(50,50)


# In[ ]:


# Find ibd diagnoses by string
ibd_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["crohn", "ulcerative colitis", "ulcerative (chronic) pancolitis", 
                 "ulcerative (chronic) proctitis", "ulcerative (chronic) rectosigmoiditis",
                 "inflammatory polyps of colon", "left sided colitis", "indeterminate colitis",
                 "collagenous colitis", "lymphocytic colitis", "microscopic colitis",
                  "regional enteritis", "ulcerative (chronic) enterocolitis", "pseudopolyposis of colon",
                  "ulcerative (chronic) ileocolitis", "ulcerative (chronic) colitis",
                  "ulcerative (chronic) proctosigmoiditis"],
).execute_gbq())


# In[ ]:


# Create a combined dataframe with source labels
phecode_df = ibd_codes_phecode_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("phecode").alias("source"))
string_df = ibd_codes_string_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("string").alias("source"))

# Find rows in phecode_df not in string_df
unique_to_phecode = phecode_df.join(
    string_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Find rows in string_df not in phecode_df
unique_to_string = string_df.join(
    phecode_df.select(["vocabulary_id", "concept_code"]),
    on=["vocabulary_id", "concept_code"],
    how="anti"
)

# Concatenate the unique rows
unique_rows_df = pl.concat([unique_to_phecode, unique_to_string])


# In[ ]:


if unique_to_phecode.height == 0:
    print('No codes unique to phecode')
else:
    display(unique_to_phecode)


# In[ ]:


if unique_to_string.height == 0:
    print('No codes unique to phecode')
else:
    display(unique_to_string.slice(0,50))


# In[ ]:


ibd_code_df = (aou.person_code_df(
    name="ibd",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["crohn", "ulcerative colitis", "ulcerative (chronic) pancolitis", 
                 "ulcerative (chronic) proctitis", "ulcerative (chronic) rectosigmoiditis",
                 "inflammatory polyps of colon", "left sided colitis", "indeterminate colitis",
                 "collagenous colitis", "lymphocytic colitis", "microscopic colitis",
                  "regional enteritis", "ulcerative (chronic) enterocolitis", "pseudopolyposis of colon",
                  "ulcerative (chronic) ileocolitis", "ulcerative (chronic) colitis",
                  "ulcerative (chronic) proctosigmoiditis"],
    exclude_terms=['fh: ', 'family history'],
    dates=True
).execute_gbq())


# ## Primary Biliary Cholangitis, Primary Sclerosing Cholangitis

# In[ ]:


pbc_dict = create_icd_dict_from_phecodes(phecodex_map, ["GI_542.81"]) 


# In[ ]:


pbc_codes_phecode_df = (aou.find_diagnosis_codes(
    exact_codes=pbc_dict
).execute_gbq())


# In[ ]:


pbc_codes_phecode_df


# In[ ]:


pbc_dict = {
    "ICD10CM": ["K74.3%", "K83.01%"],
    "ICD9CM": ["571.6%", "576.1%"],
}

# Find PBC/PSC diagnoses
pbc_codes_icd_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    pattern_codes=pbc_dict,
).execute_gbq())


# In[ ]:


pbc_codes_icd_df


# In[ ]:


# Find pbc diagnoses by string
pbc_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    exact_codes={"ICD9CM": ["576.1"], "SNOMED": ["82403002"]},
    search_terms=["primary biliary cirrhosis", "sclerosing cholangitis", "biliary cirrhosis"],
    exclude_terms=["secondary biliary cirrhosis", "cirrhosis, unspecified", "igg4"]
).execute_gbq())


# In[ ]:


pbc_code_df = (aou.person_code_df(
    name="pbc",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    exact_codes={"ICD9CM": ["576.1"], "SNOMED": ["82403002"]},
    search_terms=["primary biliary cirrhosis", "sclerosing cholangitis", "biliary cirrhosis"],
    exclude_terms=["secondary biliary cirrhosis", "cirrhosis, unspecified", "igg4"],
    dates=True
).execute_gbq())


# ## Autoimmune Hepatitis

# In[ ]:


aih_dict = {
    "ICD10CM": ["K75.4%"],
    "ICD9CM": ["571.42%"],
}

# Find aih diagnoses by string
aih_codes_icd_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    pattern_codes=aih_dict,
).execute_gbq())


# In[ ]:


aih_codes_icd_df


# In[ ]:


# Find aih diagnoses by string
aih_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["autoimmune hepatitis"],
).execute_gbq())


# In[ ]:


aih_code_df = (aou.person_code_df(
    name="aih",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["autoimmune hepatitis"],
    dates=True
).execute_gbq())


# ## HIV

# In[ ]:


hiv_dict = {
    'ICD9CM': ['042%', '043%', '044%', '079.53', '795.78', 'V08'],
    'ICD10CM': ['B20%', 'B21%', 'B22%', 'B23%', 'B24', 'B97.35', 'Z21', 'O98.7%']
}

# Find HIV diagnoses
hiv_codes_icd_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    pattern_codes=hiv_dict,
).execute_gbq())


# In[ ]:


hiv_codes_icd_df


# In[ ]:


# Find hiv diagnoses by string
hiv_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["human immunodeficiency virus"],
    exclude_terms=["family", "inconclusive", "encounter for screening", "human immunodeficiency virus antibody positive", 
                   "exposure to human immunodeficiency virus", "counseling", "nonspecific serologic evidence",]
).execute_gbq())


# In[ ]:


hiv_codes_string_df


# In[ ]:


hiv_code_df = (aou.person_code_df(
    name="hiv",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["human immunodeficiency virus"],
    exclude_terms=["family", "inconclusive", "encounter for screening", "human immunodeficiency virus antibody positive", 
                   "exposure to human immunodeficiency virus", "counseling", "nonspecific serologic evidence"],
    dates=True
).execute_gbq())


# In[ ]:


version = os.getenv('WORKSPACE_CDR')


# In[ ]:


hiv_drug_code_q = f"""
WITH filtered_concepts AS (
    SELECT DISTINCT concept_id, concept_code, concept_name
    FROM {dataset}.concept c
    WHERE vocabulary_id = 'RxNorm'
    AND (
    LOWER(concept_name) LIKE '%dideoxycytidine%'
    OR LOWER(concept_name) LIKE '%zalcitabine%'
    OR LOWER(concept_name) LIKE '%zidovudine%'
    OR LOWER(concept_name) LIKE '%azidothymidine%'
    OR LOWER(concept_name) LIKE '%didanosine%'
    OR LOWER(concept_name) LIKE '%dideoxyinosine%'
    OR LOWER(concept_name) LIKE '%stavudine%'
    OR LOWER(concept_name) LIKE '%abacavir%'
    OR LOWER(concept_name) LIKE '%rilpivirine%'
    OR LOWER(concept_name) LIKE '%etravirine%'
    OR LOWER(concept_name) LIKE '%delaviridine%'
    OR LOWER(concept_name) LIKE '%efavirenz%'
    OR LOWER(concept_name) LIKE '%nevirapine%'
    OR LOWER(concept_name) LIKE '%amprenavir%'
    OR LOWER(concept_name) LIKE '%tipranavir%'
    OR LOWER(concept_name) LIKE '%indinavir%'
    OR LOWER(concept_name) LIKE '%saquinavir%'
    OR LOWER(concept_name) LIKE '%fosamprenavir%'
    OR LOWER(concept_name) LIKE '%lopinavir%'
    OR LOWER(concept_name) LIKE '%darunavir%'
    OR LOWER(concept_name) LIKE '%atazanavir%'
    OR LOWER(concept_name) LIKE '%nelfinavir%'
    OR LOWER(concept_name) LIKE '%enfuvirtide%'
    OR LOWER(concept_name) LIKE '%maraviroc%'
    OR LOWER(concept_name) LIKE '%raltegravir%'
    OR LOWER(concept_name) LIKE '%dolutegravir%'
    OR LOWER(concept_name) LIKE '%elvitegravir%'
    OR LOWER(concept_name) LIKE '%cobicistat%'
    OR LOWER(concept_name) LIKE '%bictegravir%'
    OR LOWER(concept_name) LIKE '%ibalizumab%'
    OR LOWER(concept_name) LIKE '%fostemsavir%'
    OR LOWER(concept_name) LIKE '%lenacapavir%'
    OR LOWER(concept_name) LIKE '%hivid%'
    OR LOWER(concept_name) LIKE '%retrovir%'
    OR LOWER(concept_name) LIKE '%videx%'
    OR LOWER(concept_name) LIKE '%zerit%'
    OR LOWER(concept_name) LIKE '%ziagen%'
    OR LOWER(concept_name) LIKE '%edurant%'
    OR LOWER(concept_name) LIKE '%intelence%'
    OR LOWER(concept_name) LIKE '%rescriptor%'
    OR LOWER(concept_name) LIKE '%sustiva%'
    OR LOWER(concept_name) LIKE '%viramune%'
    OR LOWER(concept_name) LIKE '%agenerase%'
    OR LOWER(concept_name) LIKE '%aptivus%'
    OR LOWER(concept_name) LIKE '%crixivan%'
    OR LOWER(concept_name) LIKE '%fortovase%'
    OR LOWER(concept_name) LIKE '%invirase%'
    OR LOWER(concept_name) LIKE '%lexiva%'
    OR LOWER(concept_name) LIKE '%prezista%'
    OR LOWER(concept_name) LIKE '%reyataz%'
    OR LOWER(concept_name) LIKE '%viracept%'
    OR LOWER(concept_name) LIKE '%fuzeon%'
    OR LOWER(concept_name) LIKE '%selzentry%'
    OR LOWER(concept_name) LIKE '%isentress%'
    OR LOWER(concept_name) LIKE '%tivicay%'
    OR LOWER(concept_name) LIKE '%vitekta%'
    OR LOWER(concept_name) LIKE '%tybost%'
    OR LOWER(concept_name) LIKE '%atripla%'
    OR LOWER(concept_name) LIKE '%complera%'
    OR LOWER(concept_name) LIKE '%evotaz%'
    OR LOWER(concept_name) LIKE '%prezcobix%'
    OR LOWER(concept_name) LIKE '%stribild%'
    OR LOWER(concept_name) LIKE '%combivir%'
    OR LOWER(concept_name) LIKE '%epzicom%'
    OR LOWER(concept_name) LIKE '%trizivir%'
    OR LOWER(concept_name) LIKE '%kaletra%'
    OR LOWER(concept_name) LIKE '%odefsey%'
    OR LOWER(concept_name) LIKE '%genvoya%'
    OR LOWER(concept_name) LIKE '%triumeq%'
    OR LOWER(concept_name) LIKE '%delstrigo%'
    OR LOWER(concept_name) LIKE '%symfi%'
    OR LOWER(concept_name) LIKE '%dovato%'
    OR LOWER(concept_name) LIKE '%dutrebis%'
    OR LOWER(concept_name) LIKE '%symtuza%'
    OR LOWER(concept_name) LIKE '%biktarvy%'
    OR LOWER(concept_name) LIKE '%juluca%'
    OR LOWER(concept_name) LIKE '%cabenuva%'
    OR LOWER(concept_name) LIKE '%doravirine%'
    OR LOWER(concept_name) LIKE '%pifeltro%'
    ) -- not included: emtricitabine, lamivudine, tenofovir, disoproxil, alafenamide, ritonavir
      -- norvir, emtriva, epivir, viread, truvada, descovy, vemlidy, apretude, cabotegravir
),
combined AS (
    SELECT person_id, drug_exposure_start_date AS drug_date, 
        (SELECT concept_name FROM filtered_concepts WHERE concept_id = drug_concept_id) AS drug_name,
        (SELECT concept_code FROM filtered_concepts WHERE concept_id = drug_concept_id) AS drug_code,
    FROM {dataset}.drug_exposure
    WHERE drug_concept_id IN (SELECT concept_id FROM filtered_concepts)
),
distinct_dates AS (
    SELECT DISTINCT person_id, drug_date, drug_code, drug_name
    FROM combined
)
SELECT drug_code, drug_name, COUNT(*) AS count
FROM distinct_dates
GROUP BY drug_code, drug_name
ORDER BY count DESC
"""

hiv_drug_code_df = polars_gbq(hiv_drug_code_q)


# In[ ]:


hiv_drug_code_df.slice(0,15)


# In[ ]:


hiv_drug_q = f"""
WITH filtered_concepts AS (
    SELECT DISTINCT concept_id, concept_code, concept_name
    FROM {dataset}.concept c
    WHERE vocabulary_id = 'RxNorm'
    AND (
    LOWER(concept_name) LIKE '%dideoxycytidine%'
    OR LOWER(concept_name) LIKE '%zalcitabine%'
    OR LOWER(concept_name) LIKE '%zidovudine%'
    OR LOWER(concept_name) LIKE '%azidothymidine%'
    OR LOWER(concept_name) LIKE '%didanosine%'
    OR LOWER(concept_name) LIKE '%dideoxyinosine%'
    OR LOWER(concept_name) LIKE '%stavudine%'
    OR LOWER(concept_name) LIKE '%abacavir%'
    OR LOWER(concept_name) LIKE '%rilpivirine%'
    OR LOWER(concept_name) LIKE '%etravirine%'
    OR LOWER(concept_name) LIKE '%delaviridine%'
    OR LOWER(concept_name) LIKE '%efavirenz%'
    OR LOWER(concept_name) LIKE '%nevirapine%'
    OR LOWER(concept_name) LIKE '%amprenavir%'
    OR LOWER(concept_name) LIKE '%tipranavir%'
    OR LOWER(concept_name) LIKE '%indinavir%'
    OR LOWER(concept_name) LIKE '%saquinavir%'
    OR LOWER(concept_name) LIKE '%fosamprenavir%'
    OR LOWER(concept_name) LIKE '%lopinavir%'
    OR LOWER(concept_name) LIKE '%darunavir%'
    OR LOWER(concept_name) LIKE '%atazanavir%'
    OR LOWER(concept_name) LIKE '%nelfinavir%'
    OR LOWER(concept_name) LIKE '%enfuvirtide%'
    OR LOWER(concept_name) LIKE '%maraviroc%'
    OR LOWER(concept_name) LIKE '%raltegravir%'
    OR LOWER(concept_name) LIKE '%dolutegravir%'
    OR LOWER(concept_name) LIKE '%elvitegravir%'
    OR LOWER(concept_name) LIKE '%cobicistat%'
    OR LOWER(concept_name) LIKE '%bictegravir%'
    OR LOWER(concept_name) LIKE '%ibalizumab%'
    OR LOWER(concept_name) LIKE '%fostemsavir%'
    OR LOWER(concept_name) LIKE '%lenacapavir%'
    OR LOWER(concept_name) LIKE '%hivid%'
    OR LOWER(concept_name) LIKE '%retrovir%'
    OR LOWER(concept_name) LIKE '%videx%'
    OR LOWER(concept_name) LIKE '%zerit%'
    OR LOWER(concept_name) LIKE '%ziagen%'
    OR LOWER(concept_name) LIKE '%edurant%'
    OR LOWER(concept_name) LIKE '%intelence%'
    OR LOWER(concept_name) LIKE '%rescriptor%'
    OR LOWER(concept_name) LIKE '%sustiva%'
    OR LOWER(concept_name) LIKE '%viramune%'
    OR LOWER(concept_name) LIKE '%agenerase%'
    OR LOWER(concept_name) LIKE '%aptivus%'
    OR LOWER(concept_name) LIKE '%crixivan%'
    OR LOWER(concept_name) LIKE '%fortovase%'
    OR LOWER(concept_name) LIKE '%invirase%'
    OR LOWER(concept_name) LIKE '%lexiva%'
    OR LOWER(concept_name) LIKE '%prezista%'
    OR LOWER(concept_name) LIKE '%reyataz%'
    OR LOWER(concept_name) LIKE '%viracept%'
    OR LOWER(concept_name) LIKE '%fuzeon%'
    OR LOWER(concept_name) LIKE '%selzentry%'
    OR LOWER(concept_name) LIKE '%isentress%'
    OR LOWER(concept_name) LIKE '%tivicay%'
    OR LOWER(concept_name) LIKE '%vitekta%'
    OR LOWER(concept_name) LIKE '%tybost%'
    OR LOWER(concept_name) LIKE '%atripla%'
    OR LOWER(concept_name) LIKE '%complera%'
    OR LOWER(concept_name) LIKE '%evotaz%'
    OR LOWER(concept_name) LIKE '%prezcobix%'
    OR LOWER(concept_name) LIKE '%stribild%'
    OR LOWER(concept_name) LIKE '%combivir%'
    OR LOWER(concept_name) LIKE '%epzicom%'
    OR LOWER(concept_name) LIKE '%trizivir%'
    OR LOWER(concept_name) LIKE '%kaletra%'
    OR LOWER(concept_name) LIKE '%odefsey%'
    OR LOWER(concept_name) LIKE '%genvoya%'
    OR LOWER(concept_name) LIKE '%triumeq%'
    OR LOWER(concept_name) LIKE '%delstrigo%'
    OR LOWER(concept_name) LIKE '%symfi%'
    OR LOWER(concept_name) LIKE '%dovato%'
    OR LOWER(concept_name) LIKE '%dutrebis%'
    OR LOWER(concept_name) LIKE '%symtuza%'
    OR LOWER(concept_name) LIKE '%biktarvy%'
    OR LOWER(concept_name) LIKE '%juluca%'
    OR LOWER(concept_name) LIKE '%cabenuva%'
    OR LOWER(concept_name) LIKE '%doravirine%'
    OR LOWER(concept_name) LIKE '%pifeltro%'
    ) -- not included: emtricitabine, lamivudine, tenofovir, disoproxil, alafenamide, ritonavir
      -- norvir, emtriva, epivir, viread, truvada, descovy, vemlidy, apretude, cabotegravir
),
combined AS (
    SELECT de.person_id, drug_exposure_start_date AS drug_date
    FROM {dataset}.drug_exposure de
    WHERE drug_concept_id IN (SELECT concept_id FROM filtered_concepts)
),
distinct_dates AS (
    SELECT DISTINCT person_id, drug_date
    FROM combined
)
SELECT 
    person_id,
    MIN(drug_date) AS hiv_drug
FROM distinct_dates
GROUP BY person_id
"""

hiv_drug_df = polars_gbq(hiv_drug_q)


# In[ ]:


hiv_df = hiv_code_df.join(
    hiv_drug_df, 
    on='person_id',
    how='full', coalesce=True 
)


# In[ ]:


# Select positive Ab tests using value_as_unit_concept_id
hiv_ab_q = f"""
WITH filtered_concepts AS (
    SELECT DISTINCT concept_id, concept_code, concept_name
    FROM {dataset}.concept c
    WHERE (vocabulary_id = 'CPT4' OR vocabulary_id = 'LOINC')
    AND concept_code IN ('86701', '86703', '86702', '29327-4', '33866-5', 
    '41144-7', '35437-3', '7917-8', '14092-1', '29893-5', '5221-7', '86233-4', 
    '68961-2', '49905-3', '43599-0', '5220-9', '13499-9', '16977-1', '24012-7', 
    '44531-2', '5222-5', '16976-3', '44607-0', '69668-2', '80203-3', '42768-2', 
    '89365-1', '9661-0', '14126-7', '9660-2', '35452-2', '9662-8', '40438-4', 
    '12859-5', '9664-4', '42339-2', '9821-0', '18396-2', '9665-1', '9666-9', 
    '35564-4', '35565-1', '9667-7', '9668-5', '12856-1', '7918-6', '44873-8', 
    '44533-8', '31201-7', '80387-4', '43010-8', '42600-7', '49580-4', '22357-8', 
    '5223-3', '75666-8', '56888-1', '58900-2', '73906-0', '43009-0', '7919-4', 
    '5225-8', '30361-0', '81641-3', '5224-1', '10902-5', '21338-9', '21339-7', 
    '11081-7', '87390', '87389', '87391', '87806', '22356-0', '40439-2', '53601-1', 
    '22358-6', '28004-0', '35450-6', '73905-2', '85380-4', '33806-1', '40732-0')
),
combined AS (
    SELECT person_id, measurement_date AS test_date, 
    FROM {dataset}.measurement m
    WHERE measurement_concept_id IN (SELECT concept_id FROM filtered_concepts) 
    AND value_as_concept_id IN (36032716, 45878745, 45884084, 36310979, 45876384, 45879438, 45881802, 45877985)
        -- Presumptive positive, Abnormal, Positive, HIV-1 Positive, High, Present, Reactive, Detected
),
distinct_dates AS (
    SELECT DISTINCT person_id, test_date
    FROM combined
)
SELECT 
    person_id,
    MIN(test_date) AS hiv_pos_ab
FROM distinct_dates
GROUP BY person_id
"""

hiv_ab_df = polars_gbq(hiv_ab_q)


# In[ ]:


hiv_df = hiv_df.join(
    hiv_ab_df, 
    on='person_id',
    how='full', coalesce=True 
)


# In[ ]:


# Select positive VL tests
hiv_pos_vl_q = f"""
WITH filtered_concepts AS (
    SELECT DISTINCT concept_id, concept_code, concept_name
    FROM {dataset}.concept c
    WHERE (vocabulary_id = 'CPT4' OR vocabulary_id = 'LOINC')
    AND concept_code IN ('49890-7', '41514-1', '41516-6', '41498-7', '70241-5', 
    '59419-2', '20447-9', '21008-8', '48551-6', '48511-0', '41515-8', '25836-8', 
    '21333-0', '24013-5', '41497-9', '29539-4', '29541-0', '51780-5', '48510-2', 
    '5017-9', '25835-0', '47359-5', '23876-6', '62469-2', '10351-5', '69353-1', 
    '69354-9', '87535', '87534', '87536', '87538', '87537', '48552-4', '41513-3', 
    '81652-0', '86547-7', '86548-5', '86549-3', '50624-6', '81246-1', '5018-7', 
    '77369-7', '85361-4', '85368-9', '79379-4', '44871-2', '9837-6', '48023-6', 
    '34699-9', '9836-8', '25841-8', '25842-6', '73659-5', '73658-7')
),
combined AS (
    SELECT person_id, measurement_date AS test_date, 
    FROM {dataset}.measurement m
    WHERE measurement_concept_id IN (SELECT concept_id FROM filtered_concepts) 
    AND (
        (value_as_concept_id IN (45878745, 45876384, 45877985, 1620615, 1620405, 36303220, 1620483, 
        36308305, 45884084, 45884759, 45880601, 45878162))
            -- Abnormal, High, Detected, 50, 60, 100, 90, 101, Positive, 3000, 575, 450
        OR 
        (value_as_number > 20 AND unit_concept_id IN (8799, 0, 45756935, 8784, 4211671, 8510))
            --copies per milliliter, no matching concept, no value, cells per microliter, cpy/mL, unit
        OR 
        (value_as_number > log10(20) AND unit_concept_id IN (8873, 9348, 9084, 8492, 4302813, 45757552))
            --log copies per milliliter, log10 unit per milliliter, log international unit per milliliter
            --Log, log10, Log copies/mL
    )
),
distinct_dates AS (
    SELECT DISTINCT person_id, test_date
    FROM combined
)
SELECT 
    person_id,
    MIN(test_date) AS hiv_pos_vl
FROM distinct_dates
GROUP BY person_id
"""

hiv_pos_vl_df = polars_gbq(hiv_pos_vl_q)


# In[ ]:


hiv_df = hiv_df.join(
    hiv_pos_vl_df, 
    on='person_id',
    how='full', coalesce=True 
)


# In[ ]:


hiv_vl_perf_q = f"""
WITH filtered_concepts AS (
    SELECT DISTINCT concept_id, concept_code, concept_name
    FROM {dataset}.concept c
    WHERE (vocabulary_id = 'CPT4' OR vocabulary_id = 'LOINC')
    AND concept_code IN ('49890-7', '41514-1', '41516-6', '41498-7', '70241-5', 
    '59419-2', '20447-9', '21008-8', '48551-6', '48511-0', '41515-8', '25836-8', 
    '21333-0', '24013-5', '41497-9', '29539-4', '29541-0', '51780-5', '48510-2', 
    '5017-9', '25835-0', '47359-5', '23876-6', '62469-2', '10351-5', '69353-1', 
    '69354-9', '87535', '87534', '87536', '87538', '87537', '48552-4', '41513-3', 
    '81652-0', '86547-7', '86548-5', '86549-3', '50624-6', '81246-1', '5018-7', 
    '77369-7', '85361-4', '85368-9', '79379-4', '44871-2', '9837-6', '48023-6', 
    '34699-9', '9836-8', '25841-8', '25842-6', '73659-5', '73658-7')
),
combined AS (
    SELECT person_id, measurement_date AS test_date, 
    FROM {dataset}.measurement m
    WHERE measurement_concept_id IN (SELECT concept_id FROM filtered_concepts)
    AND value_as_concept_id NOT IN (45884091, 45877494, 45878602, 45878680, 45880107, 45884199, 45884087)
    -- Indeterminate, DNR, Not tested, Refused, N/A, Test not performed, Equivocal
    -- reduces repeats within 5 days from 221 to 68
    AND NOT (value_as_concept_id = 0 AND value_as_number IS NULL)
    -- reduces repeats within 5 days from 68 to 22
),
distinct_dates AS (
    SELECT DISTINCT person_id, test_date
    FROM combined
),
ranked AS (
    SELECT 
        person_id, 
        test_date,
        ROW_NUMBER() OVER(PARTITION BY person_id ORDER BY test_date) AS rn
    FROM distinct_dates
)
SELECT
    person_id, 
    MAX(CASE WHEN rn = 1 THEN test_date ELSE NULL END) AS hiv_1_vl,
    MAX(CASE WHEN rn = 2 THEN test_date ELSE NULL END) AS hiv_2_vl
FROM ranked
GROUP BY person_id
"""

hiv_vl_perf_df = polars_gbq(hiv_vl_perf_q)


# In[ ]:


hiv_df = hiv_df.join(
    hiv_vl_perf_df, 
    on='person_id',
    how='full', coalesce=True 
)


# In[ ]:


hiv_df.head()


# In[ ]:


hiv_df = hiv_df.with_columns([
    pl.when(pl.col("hiv_pos_ab").is_null())
      .then(pl.col("hiv_pos_vl"))
      .when(pl.col("hiv_pos_vl").is_null())
      .then(pl.col("hiv_pos_ab"))
      .otherwise(pl.when((pl.col("hiv_pos_ab") - pl.col("hiv_pos_vl")) > pl.duration(days=0))
                .then(pl.col("hiv_pos_vl"))
                .otherwise(pl.col("hiv_pos_ab")))
      .alias("cond_1"),

    pl.when(pl.col("hiv_1_vl").is_not_null() & pl.col("hiv_drug").is_not_null())
      .then(
          pl.when((pl.col("hiv_1_vl") - pl.col("hiv_drug")).dt.total_days() > 0)
            .then(pl.col("hiv_1_vl"))
            .otherwise(pl.col("hiv_drug"))
      )
      .otherwise(pl.lit(None))
      .alias("cond_2"),
    
    pl.when(pl.col("hiv_1").is_not_null() & pl.col("hiv_drug").is_not_null())
      .then(
          pl.when((pl.col("hiv_1") - pl.col("hiv_drug")).dt.total_days() > 0)
            .then(pl.col("hiv_1"))
            .otherwise(pl.col("hiv_drug"))
      )
      .otherwise(pl.lit(None))
      .alias("cond_3"),
    
    pl.when(pl.col("hiv_1").is_not_null() & pl.col("hiv_2_vl").is_not_null())
      .then(
          pl.when((pl.col("hiv_1") - pl.col("hiv_2_vl")).dt.total_days() > 0)
            .then(pl.col("hiv_1"))
            .otherwise(pl.col("hiv_2_vl"))
      )
      .otherwise(pl.lit(None))
      .alias("cond_4"),
    
    pl.when(pl.col("hiv_2").is_not_null())
      .then(pl.col("hiv_2"))
      .otherwise(pl.lit(None))
      .alias("cond_5")
])


# In[ ]:


hiv_pd_df = hiv_df.select('person_id', 'cond_1', 'cond_2', 'cond_3', 'cond_4', 'cond_5').to_pandas()


# In[ ]:


hiv_pd_df['hiv_2'] = hiv_pd_df[['cond_1', 'cond_2', 'cond_3', 'cond_4', 'cond_5']].min(axis=1)


# In[ ]:


hiv_final = pl.from_pandas(hiv_pd_df)
hiv_final = hiv_final.drop('cond_1', 'cond_2', 'cond_3', 'cond_4', 'cond_5')
hiv_final = hiv_final.with_columns(
    hiv_final['hiv_2'].cast(pl.Date)
)
hiv_final = hiv_final.with_columns(
    pl.lit("1990-01-01").str.strptime(pl.Date, "%Y-%m-%d").alias("hiv_1")
).select(['person_id', 'hiv_1', 'hiv_2'])


# In[ ]:


venn_df = pl.from_pandas(hiv_pd_df)


# In[ ]:


set_cond_1 = set(venn_df.filter(pl.col("cond_1").is_not_null())["person_id"].to_numpy())
set_cond_2 = set(venn_df.filter(pl.col("cond_2").is_not_null())["person_id"].to_numpy())
set_cond_3 = set(venn_df.filter(pl.col("cond_3").is_not_null())["person_id"].to_numpy())
set_cond_4 = set(venn_df.filter(pl.col("cond_4").is_not_null())["person_id"].to_numpy())
set_cond_5 = set(venn_df.filter(pl.col("cond_5").is_not_null())["person_id"].to_numpy())


# In[ ]:


from venny4py.venny4py import *

sets = {
    '+Ab or +VL': set_cond_1,
    'VL done and Tx': set_cond_2,
    'ICD and Tx': set_cond_3,
#     'ICD and two VL done': set_cond_4,
    'two ICD': set_cond_5,
}

venny4py(sets=sets)


# In[ ]:


filtered_df = hiv_df.filter(pl.col('person_id').is_in(hiv_final.filter(pl.col('hiv_2').is_null())['person_id']))


# In[ ]:


# counts for those who don't meed conditions
filtered_df = filtered_df.with_columns([
    (pl.when(pl.col('hiv_1').is_not_null()).then(pl.lit('hiv_1')).otherwise(pl.lit('')).cast(pl.Utf8) + '_' +
     pl.when(pl.col('hiv_2').is_not_null()).then(pl.lit('hiv_2')).otherwise(pl.lit('')).cast(pl.Utf8) + '_' +
     pl.when(pl.col('hiv_drug').is_not_null()).then(pl.lit('hiv_drug')).otherwise(pl.lit('')).cast(pl.Utf8) + '_' +
     pl.when(pl.col('hiv_pos_ab').is_not_null()).then(pl.lit('hiv_pos_ab')).otherwise(pl.lit('')).cast(pl.Utf8) + '_' +
     pl.when(pl.col('hiv_pos_vl').is_not_null()).then(pl.lit('hiv_pos_vl')).otherwise(pl.lit('')).cast(pl.Utf8) + '_' +
     pl.when(pl.col('hiv_1_vl').is_not_null()).then(pl.lit('hiv_1_vl')).otherwise(pl.lit('')).cast(pl.Utf8) + '_' +
     pl.when(pl.col('hiv_2_vl').is_not_null()).then(pl.lit('hiv_2_vl')).otherwise(pl.lit('')).cast(pl.Utf8)).alias('non_null_combination')
])

filtered_df.group_by('non_null_combination').agg(pl.len().alias('counts')).sort('counts', descending=True)


# ## Merge

# In[ ]:


# Function to select only person_id and columns ending with '_2'
def select_columns(df):
    # Get all columns that end with '_2'
    cols_ending_with_2 = [col for col in df.columns if col.endswith('_2')]
    # Return dataframe with only person_id and those columns
    if cols_ending_with_2:
        return df.select(['person_id'] + cols_ending_with_2)
    else:
        return df.select('person_id')

# Apply the function to each dataframe before merging
tb_code_filtered = select_columns(tb_code_df)

# List of all dataframes to merge (filtered)
dfs_to_merge_filtered = [
    select_columns(tb_skin_code_df),
    select_columns(ntm_code_df),
    select_columns(histo_code_df),
    select_columns(blasto_code_df),
    select_columns(cocci_code_df),
    select_columns(aspergillus_code_df),
    select_columns(ild_code_df),
    select_columns(hp_code_df),
    select_columns(pneumoconiosis_code_df),
    select_columns(vasculitis_code_df),
    select_columns(cvid_code_df),
    select_columns(igg4_code_df),
    select_columns(lymphoma_code_df),
    select_columns(lg_code_df),
    select_columns(lung_ca_code_df),
    select_columns(lch_code_df),
    select_columns(ibd_code_df),
    select_columns(pbc_code_df),
    select_columns(aih_code_df),
    select_columns(hiv_final)
]

# Start with the first filtered dataframe
exclude_df = tb_code_filtered

# Merge each filtered dataframe one by one
for df in dfs_to_merge_filtered:
    exclude_df = exclude_df.join(df, on='person_id', how='full', coalesce=True)
    
exclude_df = exclude_df.join(exclusion_drug_df, on='person_id', how='full', coalesce=True)


# In[ ]:


exclude_df.head()


# # Labs and Procedures

# ## CXR

# In[ ]:


# Find CXRs
cxr_procedures_df = (aou.find_procedure_codes(
    search_terms=["chest"]
).execute_gbq())


# In[ ]:


# Find CXRs
cxr_procedures_df = (aou.find_procedure_codes(
    search_terms=["chest x-ray", "radiologic examination, chest", "radiography of chest", "xr chest"],
    exclude_terms=['of the following', 'special views']
).execute_gbq())


# In[ ]:


cxr_df = (aou.person_procedure_df(
    name="cxr",
    search_terms=["chest x-ray", "radiologic examination, chest", "radiography of chest", "xr chest"],
    exclude_terms=['of the following', 'special views'],
    dates=True
).execute_gbq())


# In[ ]:


print(f"Participants with 1 procedure: {cxr_df.filter(pl.col('cxr_n')==1).height}")
print(f"Participants with ≥ 2 procedures: {cxr_df.filter(pl.col('cxr_n')>1).height}")


# In[ ]:


cxr_df.head()


# ## CT Chest

# In[ ]:


# Find sarcoid diagnoses
ct_chest_procedures_df = (aou.find_procedure_codes(
    search_terms=["ct chest", "computed tomography of chest", "(ct scan) of chest",
                "ct of chest", "computed tomography, thorax", "(ct scan) of thorax"], #computed tomographic angiography, chest
    exclude_terms=["low dose for lung cancer screening", "abdomen"]
).execute_gbq())


# In[ ]:


ct_chest_df = (aou.person_procedure_df(
    name="ct_chest",
    search_terms=["ct chest", "computed tomography of chest", "(ct scan) of chest",
                "ct of chest", "computed tomography, thorax", "(ct scan) of thorax"],
    # computed tomographic angiography, chest NOT included
    # no high-resolution CT found
    exclude_terms=["low dose for lung cancer screening", "abdomen"],
    dates=True
).execute_gbq())


# In[ ]:


print(f"Participants with 1 procedure: {ct_chest_df.filter(pl.col('ct_chest_n')==1).height}")
print(f"Participants with ≥ 2 procedures: {ct_chest_df.filter(pl.col('ct_chest_n')>1).height}")


# In[ ]:


ct_chest_df.head()


# ## Biopsy

# In[ ]:


# Find * procedures
biopsy_procedures_df = (aou.find_procedure_codes(
    search_terms=["biopsy"],
    exclude_terms = [
    # Reproductive/Gynecological
    "cervix", "cervical", "endocervical", "uterus", "uterine", "endometrium", "endometrial",
    "ovary", "ovarian", "fallopian", "vagina", "vaginal", "vulva", "vulval", 
    "prostate", "testis", "penis", "penile",
    
    # Breast
    "breast", "mammary", "mammogram",
        
    # GI
    "biliary", "cholangiography", "pancreas", "pancreatic", "exploratory laparotomy", 
    "exploration, retroperitoneal", "abdominal or retroperitoneal mass", 
    "laparoscopy, surgical", "perirectal", "perianal", "anus", "abdominal wall", 
    "peritoneum", "intra-abdominal mass", "cul-de-sac",
        
    # Genitourinary (except kidney)
    "bladder", "ureter", "urethra", "periurethral", "cystourethroscopy",
    
    # Head/Neck (except salivary glands)
    "tongue", "lip", "palate", "uvula", "pharynx", "larynx", "trachea",
    "thyroid", "parathyroid", "adrenal", "pituitary", "pineal",
    "nasal", "sinus", "ear", "auditory", "laryngoscopy", "of mouth",
    "salivary gland", "oral tissue", "nose", "gum", "salivary", "tonsils", "mouth",
    "pharyngeal", "vocal cord", 
    
    # Extremities/Joints
    "joint", "arthro", "metacarpo", "metatarsal", "phalang", "carpal",
    "shoulder", "elbow", "wrist", "hip", "knee", "ankle", "foot", "toe",
    "clavicle", "scapula", "sternum", "rib",
    
    # MSK (except bone marrow)
    "bone biopsy", "bone", "vertebral", "femur", "tibia", "fibula",
    "humerus", "radius", "ulna", "skull", "facial bone", "muscle",
    "soft tissue", "chest wall", "blood vessel", "lymphatic structure",
    "subcutaneous tissue",
    
    # Anesthesia procedures
    "anesthesia for",
    
    # Fetal/Obstetric
    "fetal", "embryo", "oocyte", "polar body", "blastomere",
    
    # Transplant
    "transplantation medicine",
        
    # Other
    "nail unit", "transcatheter biopsy", "temporal artery", "fluoroscopy, physician",
    "fluoroscopic guidance", "radiological supervision", "guidance for ", "tumor",
    "cancer", "carcinoma", "Biopsy results not reviewed", "previous biopsy",
    "sentinel", "biopsy results reviewed", "biopsy results were not reviewed",
    "specimen narrative", "tomography guidance", "hysteroscopy", 
    "fine needle aspiration biopsy - action", "ct guidance",
        
    # CNS
    "brain", "intracranial", "spinal cord", "intraspinal neoplasm", "nerve",
    "cornea", "extraocular muscle", "orbitotomy", "cerebral meninges", "iris"
]
).execute_gbq())


# In[ ]:


# Find * procedures
pulmonary_biopsy_procedures_df = (aou.find_procedure_codes(
    search_terms=['endobronchial biopsy', 'transbronchial lung biopsy', 
                  'transbronchial needle aspiration', 'bronchial brush biopsy', 
                  'biopsy of lung', 'biopsy(ies) of lung' 'biopsy(ies) of pleura',
                  'biopsy, pleura', 'biopsy, lung or mediastinum',
                  'pleural space, with biopsy', 'pericardial sac, with biopsy',
                  'mediastinal space, with biopsy', 'mediastinoscopy, includes biopsy(ies)',
                  'mediastinotomy with exploration, drainage, removal of foreign body, or biopsy',
                  'mediastinoscopy; includes biopsy(ies)', 'transbronchial biopsy', 
                  'mediastinoscopy; with lymph node biopsy(ies)', 'thoracoscopic lung biopsy',
                  'biopsy of bronchus', 'pleural biopsy', 'biopsy of mediastinum',
                  'mediastinal biopsy', 'mediastinoscopy with biopsy', 
                  'endoscopy of bronchus with biopsy', 'bronchoscopy with biopsy'],
    exclude_terms=['fluoroscopy, physician']
).execute_gbq())


# In[ ]:


pulmonary_biopsy_df = (aou.person_procedure_df(
    name="pulm_bx",
    search_terms=['endobronchial biopsy', 'transbronchial lung biopsy', 
                  'transbronchial needle aspiration', 'bronchial brush biopsy', 
                  'biopsy of lung', 'biopsy(ies) of lung' 'biopsy(ies) of pleura',
                  'biopsy, pleura', 'biopsy, lung or mediastinum',
                  'pleural space, with biopsy', 'pericardial sac, with biopsy',
                  'mediastinal space, with biopsy', 'mediastinoscopy, includes biopsy(ies)',
                  'mediastinotomy with exploration, drainage, removal of foreign body, or biopsy',
                  'mediastinoscopy; includes biopsy(ies)', 'transbronchial biopsy', 
                  'mediastinoscopy; with lymph node biopsy(ies)', 'thoracoscopic lung biopsy',
                  'biopsy of bronchus', 'pleural biopsy', 'biopsy of mediastinum',
                  'mediastinal biopsy', 'mediastinoscopy with biopsy', 
                  'endoscopy of bronchus with biopsy', 'bronchoscopy with biopsy'],
    exclude_terms=['fluoroscopy, physician'],
    dates=True
).execute_gbq())


# In[ ]:


print(f"Participants with 1 procedure: {pulmonary_biopsy_df.filter(pl.col('pulm_bx_n')==1).height}")
print(f"Participants with ≥ 2 procedures: {pulmonary_biopsy_df.filter(pl.col('pulm_bx_n')>1).height}")


# In[ ]:


pulmonary_biopsy_df.head()


# In[ ]:


# Find * procedures
non_pulmonary_biopsy_procedures_df = (aou.find_procedure_codes(
    search_terms=[
        # GI
        'esophageal biopsy', 'transoral; with biopsy', 'biopsy of stomach',
        'anoscopy, high resolution (hra) (with magnification and chemical agent enhancement); with biopsy(ies)',
        'esophagoscopy, flexible, transoral; with transendoscopic ultrasound-guided intramural or transmural fine needle aspiration/biopsy(s)',
        'esophagogastroduodenoscopy, flexible, transoral; with transendoscopic ultrasound-guided intramural or transmural fine needle aspiration/biopsy(s)',
        '(ErCP); with biopsy', 'duodenotomy, for exploration, biopsy', 'small intestine, other than duodenum; for exploration, biopsy',
        'colotomy, for exploration, biopsy', 'biopsy of intestine by capsule', 'ileum; with biopsy', 
        'through stoma; with biopsy', 'ileal reservoir [s or j]); with biopsy',
        'transendoscopic ultrasound guided intramural or transmural fine needle aspiration/biopsy',
        'proctosigmoidoscopy, rigid; with biopsy', 'sigmoidoscopy, flexible; with biopsy',
        'colonoscopy, flexible; with biopsy', 'anoscopy; with biopsy', 'anoscopy; with high-resolution magnification (hra) (eg, colposcope, operating microscope) and chemical agent enhancement, with biopsy',
        'anoscopy, high resolution (hra) (with magnification and chemical agent enhancement); with biopsy(ies)',
        'biopsy of spleen', 'biopsy of esophagus', 'biopsy of stomach', 'biopsy of small intestine', 
        '[egd] with closed biopsy', 'biopsy of large intestine', 'intestinal biopsy', 'biopsy of rectum', 
        'biopsy of liver', 'liver biopsy', 'endoscopy and biopsy', 'sigmoidoscopy with biopsy', 'biopsy of colon',
        'endoscopic biopsy',
        
        # General
        'fine needle aspiration biopsy',
        
        # Skin
        'biopsy of skin', 'biopsy of lesion of skin',
        
        # Bone Marrow
        'bone marrow; biopsy', 'bone marrow biopsy',
        
        # LN
        'biopsy or excision of lymph node', 'lymph node sampling', 'biopsy of lymph node',
        'biopsy of axillary lymph node', 'biopsy of abdominal lymph node', 'biopsy of inguinal lymph node',

        # renal
        'renal biopsy', 'biopsy of kidney', 'renal needle biopsy', 'kidney biopsy',

        # Ocular
        'eyelid skin', 'biopsy of conjunctiva', 'biopsy of eye', 'biopsy of lacrimal gland',
        'biopsy of lacrimal sac', 

        # Cardiac
        'endomyocardial biopsy', 'biopsy of pericardium',
    ],
    exclude_terms=['anesthesia for', 'hysterectomy', 'vaginectomy', 'trachelectomy'],
).execute_gbq())


# In[ ]:


non_pulmonary_biopsy_df = (aou.person_procedure_df(
    name="non_pulm_bx",
    search_terms=[
        # GI
        'esophageal biopsy', 'transoral; with biopsy', 'biopsy of stomach',
        'anoscopy, high resolution (hra) (with magnification and chemical agent enhancement); with biopsy(ies)',
        'esophagoscopy, flexible, transoral; with transendoscopic ultrasound-guided intramural or transmural fine needle aspiration/biopsy(s)',
        'esophagogastroduodenoscopy, flexible, transoral; with transendoscopic ultrasound-guided intramural or transmural fine needle aspiration/biopsy(s)',
        '(ErCP); with biopsy', 'duodenotomy, for exploration, biopsy', 'small intestine, other than duodenum; for exploration, biopsy',
        'colotomy, for exploration, biopsy', 'biopsy of intestine by capsule', 'ileum; with biopsy', 
        'through stoma; with biopsy', 'ileal reservoir [s or j]); with biopsy',
        'transendoscopic ultrasound guided intramural or transmural fine needle aspiration/biopsy',
        'proctosigmoidoscopy, rigid; with biopsy', 'sigmoidoscopy, flexible; with biopsy',
        'colonoscopy, flexible; with biopsy', 'anoscopy; with biopsy', 'anoscopy; with high-resolution magnification (hra) (eg, colposcope, operating microscope) and chemical agent enhancement, with biopsy',
        'anoscopy, high resolution (hra) (with magnification and chemical agent enhancement); with biopsy(ies)',
        'biopsy of spleen', 'biopsy of esophagus', 'biopsy of stomach', 'biopsy of small intestine', 
        '[egd] with closed biopsy', 'biopsy of large intestine', 'intestinal biopsy', 'biopsy of rectum', 
        'biopsy of liver', 'liver biopsy', 'endoscopy and biopsy', 'sigmoidoscopy with biopsy', 'biopsy of colon',
        'endoscopic biopsy',
        
        # General
        'fine needle aspiration biopsy',
        
        # Skin
        'biopsy of skin', 'biopsy of lesion of skin',
        
        # Bone Marrow
        'bone marrow; biopsy', 'bone marrow biopsy',
        
        # LN
        'biopsy or excision of lymph node', 'lymph node sampling', 'biopsy of lymph node',
        'biopsy of axillary lymph node', 'biopsy of abdominal lymph node', 'biopsy of inguinal lymph node',

        # renal
        'renal biopsy', 'biopsy of kidney', 'renal needle biopsy', 'kidney biopsy',

        # Ocular
        'eyelid skin', 'biopsy of conjunctiva', 'biopsy of eye', 'biopsy of lacrimal gland',
        'biopsy of lacrimal sac', 

        # Cardiac
        'endomyocardial biopsy', 'biopsy of pericardium',
    ],
    exclude_terms=['anesthesia for', 'hysterectomy', 'vaginectomy', 'trachelectomy'],
    dates=True
).execute_gbq())


# In[ ]:


print(f"Participants with 1 procedure: {non_pulmonary_biopsy_df.filter(pl.col('non_pulm_bx_n')==1).height}")
print(f"Participants with ≥ 2 procedures: {non_pulmonary_biopsy_df.filter(pl.col('non_pulm_bx_n')>1).height}")


# In[ ]:


non_pulmonary_biopsy_df.head()


# ## PFTs

# In[ ]:


# Find * procedures
pft_procedures_df = (aou.find_procedure_codes(
    search_terms=["pulmonary function", "spirometry"],
    exclude_terms=['(copd)', 'copd symptoms', 'prior to surgery', '(als)', 'exercise', '6-minute walk test']    
).execute_gbq())


# In[ ]:


pft_procedures_df['vocabulary_id', 'concept_name', 'unique_persons', 'total_events']


# In[ ]:


pft_df = (aou.person_procedure_df(
    name="pft",
    search_terms=["pulmonary function", "spirometry"],
    exclude_terms=['(copd)', 'copd symptoms', 'prior to surgery', '(als)', 'exercise', '6-minute walk test'],    
    dates=True
).execute_gbq())


# In[ ]:


print(f"Participants with 1 procedure: {pft_df.filter(pl.col('pft_n')==1).height}")
print(f"Participants with ≥ 2 procedures: {pft_df.filter(pl.col('pft_n')>1).height}")


# In[ ]:


pft_df.head()


# ## ECG Procedures

# In[ ]:


# Find * procedures
ecg_procedures_df = (aou.find_procedure_codes(
    search_terms=["electrocardiogram", 'electrokardiogram', 'ecg', 'ekg'],
    exclude_terms=["coronary artery obstruction", "telephonic transmission",
                  "mobile cardiovascular telemetry", "signal-averaged electrocardiography",
                   "tilt table evaluation", "loop recorder", "wearable cardioverter", 
                   "cardiac rehabilitation", "sleep", "clinical data stored in computers",
                   "blood pressure", "pump with oxygenator", "initial preventive", "not performed",
                   "portable ekg to facility", "fetal", "phonocardiogram", "carotid pulse tracing",
                   "benzoyl"]
).execute_gbq())


# In[ ]:


ecg_procedures_df[
    'vocabulary_id',
    'concept_name',
    'unique_persons',
    'total_events'
]


# In[ ]:


ecg_df = (aou.person_procedure_df(
    name="ecg",
    search_terms=["electrocardiogram", 'electrokardiogram', 'ecg', 'ekg'],
    exclude_terms=["coronary artery obstruction", "telephonic transmission",
                  "mobile cardiovascular telemetry", "signal-averaged electrocardiography",
                   "tilt table evaluation", "loop recorder", "wearable cardioverter", 
                   "cardiac rehabilitation", "sleep", "clinical data stored in computers",
                   "blood pressure", "pump with oxygenator", "initial preventive", "not performed",
                   "portable ekg to facility", "fetal", "phonocardiogram", "carotid pulse tracing",
                   "benzoyl"],
    dates=True
).execute_gbq())


# In[ ]:


print(f"Participants with 1 procedure: {ecg_df.filter(pl.col('ecg_n')==1).height}")
print(f"Participants with ≥ 2 procedures: {ecg_df.filter(pl.col('ecg_n')>1).height}")


# In[ ]:


ecg_df.head()


# ## ECG Abnormalities (AV block, PVCs, SVT [AF, flutter, AT])

# In[ ]:


# Find * procedures
*_procedures_df = (aou.find_procedure_codes(
    search_terms=[""]
).execute_gbq())


# In[ ]:


# Find * procedures
*_procedures_df = (aou.find_procedure_codes(
    search_terms=[],
).execute_gbq())


# In[ ]:


*_df = (aou.person_procedure_df(
    name="*",
    search_terms=[],
    dates=True
).execute_gbq())


# In[ ]:


*_df = (aou.person_procedure_df(
    name="*",
    search_terms=[],
    dates=True
).execute_gbq())


# In[ ]:


print(f"Participants with 1 procedure: {*_df.filter(pl.col('*_n')==1).height}")
print(f"Participants with ≥ 2 procedures: {*_df.filter(pl.col('*_n')>1).height}")


# In[ ]:


*_df.head()


# ## CBC
# Leukopenia (5 to 10 percent) [69], eosinophilia (3 percent) [70], and thrombocytopenia (rare) can be seen

# In[ ]:


# Initialize the explorer object 
explorer = Explorer(variable_type = "measurements", cohort = cohort, version = version)


# In[ ]:


# Display annotated variables
metadata = explorer.variables()


# ### WBC

# In[ ]:


metadata.filter(pl.col('measurement_name').str.contains('leukoc'))


# In[ ]:


# Initialize the unifier object
measurement = Unifier(measurement_cid = 3000905, 
                      file_path = file_path,
                      cohort = cohort,
                      version = version,
                      file_type = file_type,
                      save_annotation = True, 
                      drop_sites = True)

# Execute the unify() function to annotate and save measurements for the given measurement_cid
measurement.unify()


# In[ ]:


get_ipython().system('gsutil ls {bucket}/data/allofus/C2024Q3R5/measurement_data/annotations/')


# In[ ]:


wbc_time_series = pl.read_parquet(f'{bucket}/data/allofus/C2024Q3R5/measurement_data/timeseries/timeseries_3000905.parquet')


# In[ ]:


leukopenia_df = (
   wbc_time_series
   .filter(pl.col("value_as_number") < 4.4)
   .group_by("person_id")
   .agg(
       pl.col("measurement_datetime").sort().first().cast(pl.Date).alias("leukopenia_1"),
       pl.col("measurement_datetime").sort().slice(1, 1).first().cast(pl.Date).alias("leukopenia_2")
   )
)


# ### Plt

# In[ ]:


metadata.filter(pl.col('measurement_name').str.contains('platelet'))


# In[ ]:


# Initialize the unifier object
measurement = Unifier(measurement_cid = 3024929, 
                      file_path = file_path,
                      cohort = cohort,
                      version = version,
                      file_type = file_type,
                      save_annotation = True, 
                      drop_sites = True)

# Execute the unify() function to annotate and save measurements for the given measurement_cid
measurement.unify()


# In[ ]:


plt_time_series = pl.read_parquet(f'{bucket}/data/allofus/C2024Q3R5/measurement_data/timeseries/timeseries_3024929.parquet')


# In[ ]:


thrombocytopenia_df = (
   plt_time_series
   .filter(pl.col("value_as_number") < 150)
   .group_by("person_id")
   .agg(
       pl.col("measurement_datetime").sort().first().cast(pl.Date).alias("thrombocytopenia_1"),
       pl.col("measurement_datetime").sort().slice(1, 1).first().cast(pl.Date).alias("thrombocytopenia_2")
   )
)


# ### Eosinophils

# In[ ]:


metadata.filter(pl.col('measurement_name').str.contains('eosinophil'))


# In[ ]:


# Initialize the unifier object
measurement = Unifier(measurement_cid = 3028615, 
                      file_path = file_path,
                      cohort = cohort,
                      version = version,
                      file_type = file_type,
                      save_annotation = True, 
                      drop_sites = True)

# Execute the unify() function to annotate and save measurements for the given measurement_cid
measurement.unify()


# In[ ]:


aec_time_series = pl.read_parquet(f'{bucket}/data/allofus/C2024Q3R5/measurement_data/timeseries/timeseries_3028615.parquet')


# In[ ]:


aec_time_series.head()


# In[ ]:


pbe_df = (
   aec_time_series
   .filter(pl.col("value_as_number") > 0.5)
   .group_by("person_id")
   .agg(
       pl.col("measurement_datetime").sort().first().cast(pl.Date).alias("pbe_1"),
       pl.col("measurement_datetime").sort().slice(1, 1).first().cast(pl.Date).alias("pbe_2")
   )
)


# ## Ca
# hypercalcemia
# https://www.ncbi.nlm.nih.gov/books/NBK430714/

# In[ ]:


metadata.filter(pl.col('measurement_name').str.contains('calcium'))


# In[ ]:


# Initialize the unifier object
measurement = Unifier(measurement_cid = 3006906, 
                      file_path = file_path,
                      cohort = cohort,
                      version = version,
                      file_type = file_type,
                      save_annotation = True, 
                      drop_sites = True)

# Execute the unify() function to annotate and save measurements for the given measurement_cid
measurement.unify()


# In[ ]:


calcium_time_series = pl.read_parquet(f'{bucket}/data/allofus/C2024Q3R5/measurement_data/timeseries/timeseries_3006906.parquet')


# In[ ]:


hypercalcemia_df = (
   calcium_time_series
   .filter(pl.col("value_as_number") > 10.5)
   .group_by("person_id")
   .agg(
       pl.col("measurement_datetime").sort().first().cast(pl.Date).alias("hypercalcemia_1"),
       pl.col("measurement_datetime").sort().slice(1, 1).first().cast(pl.Date).alias("hypercalcemia_2")
   )
)


# ## ACE

# In[ ]:


ace_q = f"""
WITH filtered_concepts AS (
    SELECT DISTINCT concept_id, concept_code, concept_name
    FROM {dataset}.concept c
    WHERE vocabulary_id = 'LOINC'
        AND (LOWER(concept_name) LIKE '%angiotensin converting enzyme%')
)
SELECT concept_name, concept_id, COUNT(measurement_id) AS count
FROM {dataset}.measurement m
LEFT JOIN filtered_concepts fc ON fc.concept_id = m.measurement_concept_id
WHERE measurement_concept_id IN (SELECT concept_id FROM filtered_concepts) 
GROUP BY concept_name, concept_id
ORDER BY count DESC
LIMIT 5
"""

polars_gbq(ace_q)


# In[ ]:


measurement_cid = 3034780                    # OMOP measurement_concept_id 

test_concept_q = f"""
SELECT DISTINCT
    c.concept_name as lab_concept_name,
    c.concept_id as lab_concept_id, 
    c.concept_code as lab_concept_code,
    c.vocabulary_id as lab_vocab,
    c1.concept_name as standard_unit,
    count(measurement_id) AS count
FROM {dataset}.measurement m
JOIN {dataset}.concept c ON m.measurement_concept_id = c.concept_id
JOIN {dataset}.concept c1 ON m.unit_concept_id = c1.concept_id
WHERE c.concept_id = {measurement_cid}
GROUP BY lab_concept_name, lab_concept_id, lab_concept_code, lab_vocab, standard_unit
ORDER BY count DESC
LIMIT 5
"""

polars_gbq(test_concept_q)


# In[ ]:


measurement_cid = 3034780                    # OMOP measurement_concept_id 
measurement_name = "ace"                     # Self-assigned measurement_name
standard_unit = "unit per liter"             # Standard unit for the measurement
min_value = 1                                # Minimum physiologic value which is feasible for the measurement
max_value = 500                              # Maximum physiologic value which is feasible for the measurement
group = "sarcoid"                            # Grouping of the measurement


# In[ ]:


units = Mapper(measurement_cid = measurement_cid, 
               measurement_name = measurement_name,
               standard_unit = standard_unit,
               min_value = min_value,
               max_value = max_value,
               group = group,
               cohort = cohort,
               version=version,
               )


# In[ ]:


units.metadata


# In[ ]:


units.init_map()


# In[ ]:


units.predominant_unit_histogram()


# In[ ]:


current_unit = units.next_unit()    # .previous_unit() if user wants to return to the last unit
units.unit_histogram()              # Display the histogram for the current unit_concept_name


# In[ ]:


conversion_factor = 1      # conversion factor to convert to the standard unit
multimodal = 0             # 1 if distribution is multimodal or 0 if not

units.unit_update(current_unit, 
                  conversion_factor = conversion_factor, 
                  multimodal = multimodal)


# In[ ]:


units.save()


# In[ ]:


measurement_cid = 3036335                    # OMOP measurement_concept_id 

test_concept_q = f"""
SELECT DISTINCT
    c.concept_name as lab_concept_name,
    c.concept_id as lab_concept_id, 
    c.concept_code as lab_concept_code,
    c.vocabulary_id as lab_vocab,
    c1.concept_name as standard_unit,
    count(measurement_id) AS count
FROM {dataset}.measurement m
JOIN {dataset}.concept c ON m.measurement_concept_id = c.concept_id
JOIN {dataset}.concept c1 ON m.unit_concept_id = c1.concept_id
WHERE c.concept_id = {measurement_cid}
GROUP BY lab_concept_name, lab_concept_id, lab_concept_code, lab_vocab, standard_unit
ORDER BY count DESC
LIMIT 5
"""

polars_gbq(test_concept_q)


# In[ ]:


measurement_cid = 3036335                    # OMOP measurement_concept_id 
measurement_name = "ace"                     # Self-assigned measurement_name
standard_unit = "unit per liter"             # Standard unit for the measurement
min_value = 1                                # Minimum physiologic value which is feasible for the measurement
max_value = 500                              # Maximum physiologic value which is feasible for the measurement
group = "sarcoid"                            # Grouping of the measurement


# In[ ]:


units = Mapper(measurement_cid = measurement_cid, 
               measurement_name = measurement_name,
               standard_unit = standard_unit,
               min_value = min_value,
               max_value = max_value,
               group = group,
               cohort = cohort,
               version=version,
               )


# In[ ]:


units.metadata


# In[ ]:


units.init_map()


# In[ ]:


units.predominant_unit_histogram()


# In[ ]:


current_unit = units.next_unit()    # .previous_unit() if user wants to return to the last unit
units.unit_histogram()              # Display the histogram for the current unit_concept_name


# In[ ]:


conversion_factor = 1      # conversion factor to convert to the standard unit
multimodal = 0             # 1 if distribution is multimodal or 0 if not

units.unit_update(current_unit, 
                  conversion_factor = conversion_factor, 
                  multimodal = multimodal)


# In[ ]:


units.save()


# In[ ]:


# Display annotated variables
metadata = explorer.variables()
metadata.filter(pl.col('measurement_name').str.contains('ace'))


# In[ ]:


# Initialize the unifier object
measurement = Unifier(measurement_cid = 3036335, 
                      file_path = file_path,
                      cohort = cohort,
                      version = version,
                      file_type = file_type,
                      save_annotation = True, 
                      drop_sites = True)

# Execute the unify() function to annotate and save measurements for the given measurement_cid
measurement.unify()


# In[ ]:


ace_1_time_series = pl.read_parquet(f'{bucket}/data/allofus/C2024Q3R5/measurement_data/timeseries/timeseries_3036335.parquet')


# In[ ]:


# Initialize the unifier object
measurement = Unifier(measurement_cid = 3034780, 
                      file_path = file_path,
                      cohort = cohort,
                      version = version,
                      file_type = file_type,
                      save_annotation = True, 
                      drop_sites = True)

# Execute the unify() function to annotate and save measurements for the given measurement_cid
measurement.unify()


# In[ ]:


ace_2_time_series = pl.read_parquet(f'{bucket}/data/allofus/C2024Q3R5/measurement_data/timeseries/timeseries_3034780.parquet')


# In[ ]:


ace_time_series = pl.concat([ace_1_time_series, ace_2_time_series])


# In[ ]:


elevated_ace_df = (
   ace_time_series
   .filter(pl.col("value_as_number") > 65)
   .group_by("person_id")
   .agg(
       pl.col("measurement_datetime").sort().first().cast(pl.Date).alias("elevated_ace_1"),
       pl.col("measurement_datetime").sort().slice(1, 1).first().cast(pl.Date).alias("elevated_ace_2")
   )
)


# ## Alk Phos
# Elevated

# In[ ]:


metadata.filter(pl.col('measurement_name').str.contains('phosphatase'))


# In[ ]:


# Initialize the unifier object
measurement = Unifier(measurement_cid = 3035995, 
                      file_path = file_path,
                      cohort = cohort,
                      version = version,
                      file_type = file_type,
                      save_annotation = True, 
                      drop_sites = True)

# Execute the unify() function to annotate and save measurements for the given measurement_cid
measurement.unify()


# In[ ]:


ap_time_series = pl.read_parquet(f'{bucket}/data/allofus/C2024Q3R5/measurement_data/timeseries/timeseries_3035995.parquet')


# In[ ]:


elevated_alk_phos_df = (
   ap_time_series
   .filter(pl.col("value_as_number") > 147)
   .group_by("person_id")
   .agg(
       pl.col("measurement_datetime").sort().first().cast(pl.Date).alias("elevated_ap_1"),
       pl.col("measurement_datetime").sort().slice(1, 1).first().cast(pl.Date).alias("elevated_ap_2")
   )
)


# ## 1,25-OH vitamin D
# Deficiency of 25-hydroxyvitamin D is nearly universal among patients with sarcoidosis, although 1,25-dihydroxyvitamin D is sufficient in 70 percent.
# Granulomas make 1,25 OH D

# In[ ]:


ace_q = f"""
WITH filtered_concepts AS (
    SELECT DISTINCT concept_id, concept_code, concept_name
    FROM {dataset}.concept c
    WHERE vocabulary_id = 'LOINC'
        AND (LOWER(concept_name) LIKE '%vitamin d%')
)
SELECT concept_name, concept_id, COUNT(measurement_id) AS count
FROM {dataset}.measurement m
LEFT JOIN filtered_concepts fc ON fc.concept_id = m.measurement_concept_id
WHERE measurement_concept_id IN (SELECT concept_id FROM filtered_concepts) 
GROUP BY concept_name, concept_id
ORDER BY count DESC
LIMIT 5
"""

polars_gbq(ace_q)


# In[ ]:


measurement_cid = 40765040                    # OMOP measurement_concept_id 

test_concept_q = f"""
SELECT DISTINCT
    c.concept_name as lab_concept_name,
    c.concept_id as lab_concept_id, 
    c.concept_code as lab_concept_code,
    c.vocabulary_id as lab_vocab,
    c1.concept_name as standard_unit,
    count(measurement_id) AS count
FROM {dataset}.measurement m
JOIN {dataset}.concept c ON m.measurement_concept_id = c.concept_id
JOIN {dataset}.concept c1 ON m.unit_concept_id = c1.concept_id
WHERE c.concept_id = {measurement_cid}
GROUP BY lab_concept_name, lab_concept_id, lab_concept_code, lab_vocab, standard_unit
ORDER BY count DESC
LIMIT 5
"""

polars_gbq(test_concept_q)


# In[ ]:


measurement_cid = 40765040                    # OMOP measurement_concept_id 
measurement_name = "25-Hydroxyvitamin D3+25-Hydroxyvitamin D2"                     # Self-assigned measurement_name
standard_unit = "nanogram per milliliter"             # Standard unit for the measurement
min_value = 0                                # Minimum physiologic value which is feasible for the measurement
max_value = 250                              # Maximum physiologic value which is feasible for the measurement
group = "sarcoid"                            # Grouping of the measurement


# In[ ]:


units = Mapper(measurement_cid = measurement_cid, 
               measurement_name = measurement_name,
               standard_unit = standard_unit,
               min_value = min_value,
               max_value = max_value,
               group = group,
               cohort = cohort,
               version=version,
               )


# In[ ]:


units.metadata


# In[ ]:


units.init_map()


# In[ ]:


units.predominant_unit_histogram()


# In[ ]:


current_unit = units.next_unit()    # .previous_unit() if user wants to return to the last unit
units.unit_histogram()              # Display the histogram for the current unit_concept_name


# In[ ]:


conversion_factor = 1      # conversion factor to convert to the standard unit
multimodal = 0             # 1 if distribution is multimodal or 0 if not

units.unit_update(current_unit, 
                  conversion_factor = conversion_factor, 
                  multimodal = multimodal)


# In[ ]:


units.save(overwrite=True)


# In[ ]:


measurement_cid = 3049536                    # OMOP measurement_concept_id 

test_concept_q = f"""
SELECT DISTINCT
    c.concept_name as lab_concept_name,
    c.concept_id as lab_concept_id, 
    c.concept_code as lab_concept_code,
    c.vocabulary_id as lab_vocab,
    c1.concept_name as standard_unit,
    count(measurement_id) AS count
FROM {dataset}.measurement m
JOIN {dataset}.concept c ON m.measurement_concept_id = c.concept_id
JOIN {dataset}.concept c1 ON m.unit_concept_id = c1.concept_id
WHERE c.concept_id = {measurement_cid}
GROUP BY lab_concept_name, lab_concept_id, lab_concept_code, lab_vocab, standard_unit
ORDER BY count DESC
LIMIT 5
"""

polars_gbq(test_concept_q)


# In[ ]:


measurement_cid = 3049536                    # OMOP measurement_concept_id 
measurement_name = "25-hydroxyvitamin D2"                     # Self-assigned measurement_name
standard_unit = "nanogram per milliliter"             # Standard unit for the measurement
min_value = 0                                # Minimum physiologic value which is feasible for the measurement
max_value = 250                              # Maximum physiologic value which is feasible for the measurement
group = "sarcoid"                            # Grouping of the measurement


# In[ ]:


units = Mapper(measurement_cid = measurement_cid, 
               measurement_name = measurement_name,
               standard_unit = standard_unit,
               min_value = min_value,
               max_value = max_value,
               group = group,
               cohort = cohort,
               version=version,
               )


# In[ ]:


units.metadata


# In[ ]:


units.init_map()


# In[ ]:


units.predominant_unit_histogram()


# In[ ]:


current_unit = units.next_unit()    # .previous_unit() if user wants to return to the last unit
units.unit_histogram()              # Display the histogram for the current unit_concept_name


# In[ ]:


conversion_factor = 0      # conversion factor to convert to the standard unit
multimodal = 0             # 1 if distribution is multimodal or 0 if not

units.unit_update(current_unit, 
                  conversion_factor = conversion_factor, 
                  multimodal = multimodal)


# In[ ]:


units.save(overwrite=True)


# In[ ]:


measurement_cid = 40765038                    # OMOP measurement_concept_id 

test_concept_q = f"""
SELECT DISTINCT
    c.concept_name as lab_concept_name,
    c.concept_id as lab_concept_id, 
    c.concept_code as lab_concept_code,
    c.vocabulary_id as lab_vocab,
    c1.concept_name as standard_unit,
    count(measurement_id) AS count
FROM {dataset}.measurement m
JOIN {dataset}.concept c ON m.measurement_concept_id = c.concept_id
JOIN {dataset}.concept c1 ON m.unit_concept_id = c1.concept_id
WHERE c.concept_id = {measurement_cid}
GROUP BY lab_concept_name, lab_concept_id, lab_concept_code, lab_vocab, standard_unit
ORDER BY count DESC
LIMIT 5
"""

polars_gbq(test_concept_q)


# In[ ]:


measurement_cid = 40765038                    # OMOP measurement_concept_id 
measurement_name = "1,25-Dihydroxyvitamin D"                     # Self-assigned measurement_name
standard_unit = "picogram per milliliter"             # Standard unit for the measurement
min_value = 0                                # Minimum physiologic value which is feasible for the measurement
max_value = 250                              # Maximum physiologic value which is feasible for the measurement
group = "sarcoid"                            # Grouping of the measurement


# In[ ]:


units = Mapper(measurement_cid = measurement_cid, 
               measurement_name = measurement_name,
               standard_unit = standard_unit,
               min_value = min_value,
               max_value = max_value,
               group = group,
               cohort = cohort,
               version=version,
               )


# In[ ]:


units.metadata


# In[ ]:


units.init_map()


# In[ ]:


units.predominant_unit_histogram()


# In[ ]:


current_unit = units.next_unit()    # .previous_unit() if user wants to return to the last unit
units.unit_histogram()              # Display the histogram for the current unit_concept_name


# In[ ]:


conversion_factor = 1      # conversion factor to convert to the standard unit
multimodal = 0             # 1 if distribution is multimodal or 0 if not

units.unit_update(current_unit, 
                  conversion_factor = conversion_factor, 
                  multimodal = multimodal)


# In[ ]:


units.save(overwrite=True)


# In[ ]:


# Display annotated variables
metadata = explorer.variables()
metadata.filter(pl.col('measurement_name').str.contains('vitamin D'))


# In[ ]:


# Initialize the unifier object
measurement = Unifier(measurement_cid = 40765038, 
                      file_path = file_path,
                      cohort = cohort,
                      version = version,
                      file_type = file_type,
                      save_annotation = True, 
                      drop_sites = True)

# Execute the unify() function to annotate and save measurements for the given measurement_cid
measurement.unify()


# In[ ]:


vitamin_d_125_time_series = pl.read_parquet(f'{bucket}/data/allofus/C2024Q3R5/measurement_data/timeseries/timeseries_40765038.parquet')


# In[ ]:


high_1_25_vitamin_d_df = (
   vitamin_d_125_time_series
   .filter(pl.col("value_as_number") > 78)
   .group_by("person_id")
   .agg(
       pl.col("measurement_datetime").sort().first().cast(pl.Date).alias("high_1_25_d_1"),
       pl.col("measurement_datetime").sort().slice(1, 1).first().cast(pl.Date).alias("high_1_25_d_2")
   )
)


# https://emedicine.medscape.com/article/2088672-overview#:~:text=Increased%201%2C25%2Ddihydroxyvitamin%20D,increases%20in%201%CE%B1%2Dhydroxylase%20activity.
# '> 78'

# ## 25-OH vitamin D
# Vitamin D sufficiency is defined as a 25(OH)D concentration ≥20 ng/mL (50 nmol/L)<br>
# Vitamin D insufficiency is defined as a 25(OH)D concentration of 12 to <20 ng/mL (30 to 50 nmol/L)<br>
# Vitamin D deficiency is defined as a 25(OH)D level <12 ng/mL (30 nmol/L)<br>
# https://www.uptodate.com/contents/vitamin-d-deficiency-in-adults-definition-clinical-manifestations-and-treatment?search=vitamin%20d%20deficiency&source=search_result&selectedTitle=1~150&usage_type=default&display_rank=1#H2

# In[ ]:


# Display annotated variables
metadata.filter(pl.col('measurement_name').str.contains('vitamin D'))


# In[ ]:


# Initialize the unifier object
measurement = Unifier(measurement_cid = 3020149, 
                      file_path = file_path,
                      cohort = cohort,
                      version = version,
                      file_type = file_type,
                      save_annotation = True, 
                      drop_sites = True)

# Execute the unify() function to annotate and save measurements for the given measurement_cid
measurement.unify()


# In[ ]:


vid_d3_25_time_series = pl.read_parquet(f'{bucket}/data/allofus/C2024Q3R5/measurement_data/timeseries/timeseries_3020149.parquet')


# In[ ]:


# Initialize the unifier object
measurement = Unifier(measurement_cid = 3049536, 
                      file_path = file_path,
                      cohort = cohort,
                      version = version,
                      file_type = file_type,
                      save_annotation = True, 
                      drop_sites = True)

# Execute the unify() function to annotate and save measurements for the given measurement_cid
measurement.unify()


# In[ ]:


vid_d2_25_time_series = pl.read_parquet(f'{bucket}/data/allofus/C2024Q3R5/measurement_data/timeseries/timeseries_3049536.parquet')


# In[ ]:


# Initialize the unifier object
measurement = Unifier(measurement_cid = 40765040, 
                      file_path = file_path,
                      cohort = cohort,
                      version = version,
                      file_type = file_type,
                      save_annotation = True, 
                      drop_sites = True)

# Execute the unify() function to annotate and save measurements for the given measurement_cid
measurement.unify()


# In[ ]:


vid_d23_25_time_series = pl.read_parquet(f'{bucket}/data/allofus/C2024Q3R5/measurement_data/timeseries/timeseries_40765040.parquet')


# In[ ]:


# Add date column to both datasets
vid_d2_25_time_series = vid_d2_25_time_series.with_columns(
    pl.col("measurement_datetime").dt.date().alias("measurement_date")
)

vid_d3_25_time_series = vid_d3_25_time_series.with_columns(
    pl.col("measurement_datetime").dt.date().alias("measurement_date")
)


# In[ ]:


# Find same-day measurements for both D2 and D3
vid_d23_same_day = (
    vid_d2_25_time_series
    .join(
        vid_d3_25_time_series, 
        on=["person_id", "measurement_date"], 
        how="inner",
        suffix="_d3"
    )
    .with_columns([
        (pl.col("value_as_number") + pl.col("value_as_number_d3")).alias("value_as_number")
    ])
    .select([
        "person_id", 
        "measurement_date",
        "value_as_number",
        "measurement_datetime",
        "src_id",
        "standard_concept_name", 
        "measurement_name",
        "unit"
    ])
)


# In[ ]:


vid_d23_25_time_series.head()


# In[ ]:


vitamin_d_combined = pl.concat([vid_d23_same_day.select('person_id', 'measurement_datetime',
                                                       'value_as_number'),
                                vid_d23_25_time_series.select('person_id', 'measurement_datetime',
                                                       'value_as_number')])


# In[ ]:


low_25_vitamin_d_df = (
   vitamin_d_combined
   .filter(pl.col("value_as_number") < 12)
   .group_by("person_id")
   .agg(
       pl.col("measurement_datetime").sort().first().cast(pl.Date).alias("low_25_d_1"),
       pl.col("measurement_datetime").sort().slice(1, 1).first().cast(pl.Date).alias("low_25_d_2")
   )
)


# ## Vitamin D Deficiency

# In[ ]:


vit_d_deficiency_dict = create_icd_dict_from_phecodes(phecodex_map, ["EM_232.4"])


# In[ ]:


vit_d_deficiency_codes_phecode_df = (aou.find_diagnosis_codes(
    exact_codes=vit_d_deficiency_dict
).execute_gbq())


# In[ ]:


vit_d_deficiency_codes_phecode_df


# In[ ]:


# Find diagnoses by string
vit_d_deficiency_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["vitamin d deficiency", 'rickets'],
    exclude_terms=['hereditary', 'family history', 'familial']
).execute_gbq())


# In[ ]:


vit_d_deficiency_code_df = (aou.person_code_df(
    name="vit_d_deficiency",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["vitamin d deficiency", 'rickets'],
    exclude_terms=['hereditary', 'family history', 'familial'],
    dates=True
).execute_gbq())


# ## Ergocalciferol or cholecalciferol

# ## A1AT
# https://www.ncbi.nlm.nih.gov/books/NBK1519/#:~:text=Demonstration%20of%20Low%20Serum%20Concentration%20of%20the,with%20lung%20disease%20are%20usually%20%3C57%20mg/dL.
# Normal serum levels are 20-53 µmol/L or approximately 100-220 mg/dL by nephelometry.

# In[ ]:


a1at_q = f"""
WITH filtered_concepts AS (
    SELECT DISTINCT concept_id, concept_code, concept_name
    FROM {dataset}.concept c
    WHERE vocabulary_id = 'LOINC'
        AND (LOWER(concept_name) LIKE '%alpha 1 antitrypsin%')
)
SELECT concept_name, concept_id, COUNT(measurement_id) AS count
FROM {dataset}.measurement m
LEFT JOIN filtered_concepts fc ON fc.concept_id = m.measurement_concept_id
WHERE measurement_concept_id IN (SELECT concept_id FROM filtered_concepts) 
GROUP BY concept_name, concept_id
ORDER BY count DESC
LIMIT 5
"""

polars_gbq(a1at_q)


# In[ ]:


measurement_cid = 3026285                    # OMOP measurement_concept_id 

test_concept_q = f"""
SELECT DISTINCT
    c.concept_name as lab_concept_name,
    c.concept_id as lab_concept_id, 
    c.concept_code as lab_concept_code,
    c.vocabulary_id as lab_vocab,
    c1.concept_name as standard_unit,
    count(measurement_id) AS count
FROM {dataset}.measurement m
JOIN {dataset}.concept c ON m.measurement_concept_id = c.concept_id
JOIN {dataset}.concept c1 ON m.unit_concept_id = c1.concept_id
WHERE c.concept_id = {measurement_cid}
GROUP BY lab_concept_name, lab_concept_id, lab_concept_code, lab_vocab, standard_unit
ORDER BY count DESC
LIMIT 5
"""

polars_gbq(test_concept_q)


# In[ ]:


measurement_cid = 3026285                    # OMOP measurement_concept_id 
measurement_name = "a1at"                     # Self-assigned measurement_name
standard_unit = "milligram per deciliter"             # Standard unit for the measurement
min_value = 0                                # Minimum physiologic value which is feasible for the measurement
max_value = 500                              # Maximum physiologic value which is feasible for the measurement
group = "sarcoid"                            # Grouping of the measurement


# In[ ]:


units = Mapper(measurement_cid = measurement_cid, 
               measurement_name = measurement_name,
               standard_unit = standard_unit,
               min_value = min_value,
               max_value = max_value,
               group = group,
               cohort = cohort,
               version=version,
               )


# In[ ]:


units.metadata


# In[ ]:


units.init_map()


# In[ ]:


units.predominant_unit_histogram()


# In[ ]:


current_unit = units.next_unit()    # .previous_unit() if user wants to return to the last unit
units.unit_histogram()              # Display the histogram for the current unit_concept_name


# In[ ]:


conversion_factor = 0      # conversion factor to convert to the standard unit
multimodal = 0             # 1 if distribution is multimodal or 0 if not

units.unit_update(current_unit, 
                  conversion_factor = conversion_factor, 
                  multimodal = multimodal)


# In[ ]:


units.save(overwrite=True)


# In[ ]:


measurement_cid = 3001127                    # OMOP measurement_concept_id 

test_concept_q = f"""
SELECT DISTINCT
    c.concept_name as lab_concept_name,
    c.concept_id as lab_concept_id, 
    c.concept_code as lab_concept_code,
    c.vocabulary_id as lab_vocab,
    c1.concept_name as standard_unit,
    count(measurement_id) AS count
FROM {dataset}.measurement m
JOIN {dataset}.concept c ON m.measurement_concept_id = c.concept_id
JOIN {dataset}.concept c1 ON m.unit_concept_id = c1.concept_id
WHERE c.concept_id = {measurement_cid}
GROUP BY lab_concept_name, lab_concept_id, lab_concept_code, lab_vocab, standard_unit
ORDER BY count DESC
LIMIT 5
"""

polars_gbq(test_concept_q)


# In[ ]:


measurement_cid = 3001127                    # OMOP measurement_concept_id 
measurement_name = "a1at"                     # Self-assigned measurement_name
standard_unit = "gram per deciliter"             # Standard unit for the measurement
min_value = 0                                # Minimum physiologic value which is feasible for the measurement
max_value = 500                              # Maximum physiologic value which is feasible for the measurement
group = "sarcoid"                            # Grouping of the measurement


# In[ ]:


units = Mapper(measurement_cid = measurement_cid, 
               measurement_name = measurement_name,
               standard_unit = standard_unit,
               min_value = min_value,
               max_value = max_value,
               group = group,
               cohort = cohort,
               version=version,
               )


# In[ ]:


units.metadata


# In[ ]:


units.init_map()


# In[ ]:


units.predominant_unit_histogram()


# In[ ]:


current_unit = units.next_unit()    # .previous_unit() if user wants to return to the last unit
units.unit_histogram()              # Display the histogram for the current unit_concept_name


# In[ ]:


conversion_factor = 1      # conversion factor to convert to the standard unit
multimodal = 0             # 1 if distribution is multimodal or 0 if not

units.unit_update(current_unit, 
                  conversion_factor = conversion_factor, 
                  multimodal = multimodal)


# In[ ]:


units.save()


# In[ ]:


# Display annotated variables
metadata = explorer.variables()
metadata.filter(pl.col('measurement_name').str.contains('a1at'))


# In[ ]:


# Initialize the unifier object
measurement = Unifier(measurement_cid = 3001127, 
                      file_path = file_path,
                      cohort = cohort,
                      version = version,
                      file_type = file_type,
                      save_annotation = True, 
                      drop_sites = True)

# Execute the unify() function to annotate and save measurements for the given measurement_cid
measurement.unify()


# In[ ]:


a1at_1_time_series = pl.read_parquet(f'{bucket}/data/allofus/C2024Q3R5/measurement_data/timeseries/timeseries_3001127.parquet')


# In[ ]:


# Initialize the unifier object
measurement = Unifier(measurement_cid = 3026285, 
                      file_path = file_path,
                      cohort = cohort,
                      version = version,
                      file_type = file_type,
                      save_annotation = True, 
                      drop_sites = True)

# Execute the unify() function to annotate and save measurements for the given measurement_cid
measurement.unify()


# In[ ]:


a1at_2_time_series = pl.read_parquet(f'{bucket}/data/allofus/C2024Q3R5/measurement_data/timeseries/timeseries_3026285.parquet')


# In[ ]:


a1at_time_series = pl.concat([a1at_1_time_series, a1at_2_time_series])


# In[ ]:


a1at_time_series.head()


# In[ ]:


sns.histplot(data=a1at_time_series, x='value_as_number', hue='standard_concept_name')


# In[ ]:


high_a1at_df = (
   a1at_time_series
   .filter(pl.col("value_as_number") > 220)
   .group_by("person_id")
   .agg(
       pl.col("measurement_datetime").sort().first().cast(pl.Date).alias("high_a1at_1"),
       pl.col("measurement_datetime").sort().slice(1, 1).first().cast(pl.Date).alias("high_a1at_2")
   )
)


# ## soluble IL-2

# In[ ]:


sol_il2_q = f"""
WITH filtered_concepts AS (
    SELECT DISTINCT concept_id, concept_code, concept_name
    FROM {dataset}.concept c
    WHERE vocabulary_id = 'LOINC'
        AND (LOWER(concept_name) LIKE '%interleukin 2%')
)
SELECT concept_name, concept_id, COUNT(measurement_id) AS count
FROM {dataset}.measurement m
LEFT JOIN filtered_concepts fc ON fc.concept_id = m.measurement_concept_id
WHERE measurement_concept_id IN (SELECT concept_id FROM filtered_concepts) 
GROUP BY concept_name, concept_id
ORDER BY count DESC
LIMIT 5
"""

polars_gbq(sol_il2_q)


# In[ ]:


measurement_cid = 46235736                    # OMOP measurement_concept_id 

test_concept_q = f"""
SELECT DISTINCT
    c.concept_name as lab_concept_name,
    c.concept_id as lab_concept_id, 
    c.concept_code as lab_concept_code,
    c.vocabulary_id as lab_vocab,
    c1.concept_name as standard_unit,
    count(measurement_id) AS count
FROM {dataset}.measurement m
JOIN {dataset}.concept c ON m.measurement_concept_id = c.concept_id
JOIN {dataset}.concept c1 ON m.unit_concept_id = c1.concept_id
WHERE c.concept_id = {measurement_cid}
GROUP BY lab_concept_name, lab_concept_id, lab_concept_code, lab_vocab, standard_unit
ORDER BY count DESC
LIMIT 5
"""

polars_gbq(test_concept_q)


# In[ ]:


measurement_cid = 46235736                    # OMOP measurement_concept_id 
measurement_name = "soluble_il2"                     # Self-assigned measurement_name
standard_unit = "picogram per milliliter"             # Standard unit for the measurement
min_value = 0                               # Minimum physiologic value which is feasible for the measurement
max_value = 50000                              # Maximum physiologic value which is feasible for the measurement
group = "sarcoid"                            # Grouping of the measurement


# In[ ]:


units = Mapper(measurement_cid = measurement_cid, 
               measurement_name = measurement_name,
               standard_unit = standard_unit,
               min_value = min_value,
               max_value = max_value,
               group = group,
               cohort = cohort,
               version=version,
               )


# In[ ]:


units.metadata


# In[ ]:


units.init_map()


# In[ ]:


units.predominant_unit_histogram()


# In[ ]:


current_unit = units.next_unit()    # .previous_unit() if user wants to return to the last unit
units.unit_histogram()              # Display the histogram for the current unit_concept_name


# In[ ]:


conversion_factor = 1      # conversion factor to convert to the standard unit
multimodal = 0             # 1 if distribution is multimodal or 0 if not

units.unit_update(current_unit, 
                  conversion_factor = conversion_factor, 
                  multimodal = multimodal)


# In[ ]:


units.save()


# In[ ]:


# Display annotated variables
metadata = explorer.variables()
metadata.filter(pl.col('measurement_name').str.contains('soluble'))


# In[ ]:


# Initialize the unifier object
measurement = Unifier(measurement_cid = 46235736, 
                      file_path = file_path,
                      cohort = cohort,
                      version = version,
                      file_type = file_type,
                      save_annotation = True, 
                      drop_sites = True)

# Execute the unify() function to annotate and save measurements for the given measurement_cid
measurement.unify()


# In[ ]:


soluble_il2_time_series = pl.read_parquet(f'{bucket}/data/allofus/C2024Q3R5/measurement_data/timeseries/timeseries_46235736.parquet')


# In[ ]:


sns.histplot(soluble_il2_time_series['value_as_number'])


# In[ ]:


elevated_sil2_df = (
   soluble_il2_time_series
   .filter(pl.col("value_as_number") > 1000)
   .group_by("person_id")
   .agg(
       pl.col("measurement_datetime").sort().first().cast(pl.Date).alias("elevated_sil2_1"),
       pl.col("measurement_datetime").sort().slice(1, 1).first().cast(pl.Date).alias("elevated_sil2_2")
   )
)


# ## Merge Procedures and Labs

# In[ ]:


# Function to select only person_id and columns ending with '_1'
def select_columns_1(df):
    # Get all columns that end with '_1'
    cols_ending_with_1 = [col for col in df.columns if col.endswith('_1')]
    # Return dataframe with only person_id and those columns
    if cols_ending_with_1:
        return df.select(['person_id'] + cols_ending_with_1)
    else:
        return df.select('person_id')
    
def select_columns_2(df):
    # Get all columns that end with '_1'
    cols_ending_with_1 = [col for col in df.columns if col.endswith('_1')]
    # Return dataframe with only person_id and those columns
    if cols_ending_with_1:
        return df.select(['person_id'] + cols_ending_with_1)
    else:
        return df.select('person_id')

# Apply the function to each dataframe before merging
cxr_filtered = select_columns_1(cxr_df)

# List of all dataframes to merge (filtered)
dfs_to_merge_filtered = [
    select_columns_1(ct_chest_df),
    select_columns_1(pulmonary_biopsy_df),
    select_columns_1(non_pulmonary_biopsy_df),
    select_columns_1(pft_df),
    select_columns_1(ecg_df),
    select_columns_1(leukopenia_df),
    select_columns_1(thrombocytopenia_df),
    select_columns_1(pbe_df),
    select_columns_1(hypercalcemia_df),
    select_columns_1(elevated_ace_df),
    select_columns_1(elevated_alk_phos_df),
    select_columns_1(low_25_vitamin_d_df),
    select_columns_1(high_1_25_vitamin_d_df),
    select_columns_2(vit_d_deficiency_code_df),
    select_columns_1(high_a1at_df),
    select_columns_1(elevated_sil2_df)
]

# Start with the first filtered dataframe
procedure_df = cxr_filtered

# Merge each filtered dataframe one by one
for df in dfs_to_merge_filtered:
    procedure_df = procedure_df.join(df, on='person_id', how='full', coalesce=True)


# In[ ]:


procedure_df.head()


# In[ ]:


procedure_df.write_csv(f'{bucket}/data/cohorts/all_procedures.csv')


# # Exclude Sarcoid Patients

# In[ ]:


sarcoid_final_df = sarcoid_df.filter(pl.col('sarcoid_n')>0).select('person_id', 'sarcoid_n', 'vocab', 'sarcoid_1')


# In[ ]:


# Left join to merge exclusion conditions with sarcoid patients
sarcoid_with_exclusions = sarcoid_final_df.join(
    exclude_df,
    on='person_id',
    how='left'
)

# Create a new column drug_2 based on the exclusion_drug date and our logic
sarcoid_with_exclusions = sarcoid_with_exclusions.with_columns(
    drug_2=pl.when(
        (pl.col('exclusion_drug').is_not_null()) & 
        (pl.col('exclusion_drug') < (pl.col('sarcoid_1') + pl.duration(days=365)))
    )
    .then(pl.col('exclusion_drug'))
    .otherwise(None)
)

# Drop the original exclusion_drug column
sarcoid_with_exclusions = sarcoid_with_exclusions.drop('exclusion_drug')

# Create the has_exclusion_condition flag using all exclusion columns including the new drug_2
exclusion_cols = [c for c in sarcoid_with_exclusions.columns 
                  if c != 'person_id' and c.endswith('_2')]

# Set the has_exclusion_condition flag
sarcoid_with_exclusions = sarcoid_with_exclusions.with_columns(
    has_exclusion_condition=pl.any_horizontal([
        pl.col(col).is_not_null() for col in exclusion_cols
    ])
)


# In[ ]:


from upsetplot import plot
from upsetplot import from_memberships
from upsetplot import from_contents


# In[ ]:


# Create binary indicators for each exclusion condition
exclusion_cols = [col for col in sarcoid_with_exclusions.columns 
                  if col.endswith('_2')]

# Convert to pandas for UpSetPlot
pandas_df = sarcoid_with_exclusions.to_pandas()

# Create a list of sets where each set contains person_ids that have a specific condition
sets = {}
for col in exclusion_cols:
    condition = col.replace('_2', '')
    # Get person_ids for patients with this condition
    person_ids = pandas_df.loc[pandas_df[col].notnull(), 'person_id'].tolist()
    sets[condition] = set(person_ids)  # Using set for efficiency

# Create a list of memberships (each entry represents a patient's conditions)
memberships = []
for _, row in pandas_df.iterrows():
    patient_conditions = []
    for col in exclusion_cols:
        if pd.notnull(row[col]):
            condition = col.replace('_2', '')
            patient_conditions.append(condition)
    memberships.append(patient_conditions)

# Create data for UpSetPlot using the correct syntax
membership_data = from_memberships(memberships)

plt.rcParams.update({'font.size': 12})  # General font size

# Create the UpSet plot
plot(membership_data, subset_size="count", sort_by='cardinality', 
     sort_categories_by='-cardinality', 
     max_subset_rank=20,
     show_counts=True,
    )
plt.title('')
plt.show()


# In[ ]:


# Create a list to store results
count_data = []

# Calculate counts for each exclusion condition
for col in exclusion_cols:
    condition_name = col.replace('_2', '')
    # Count non-null values
    count = sarcoid_with_exclusions.filter(pl.col(col).is_not_null()).height
    count_data.append({"condition": condition_name, "count": count})

# Convert to dataframe
exclusion_counts = pl.DataFrame(count_data)

# Sort by count in descending order
exclusion_counts = exclusion_counts.sort("count", descending=True)


# In[ ]:


exclusion_counts


# # Describe Sarcoid Exclusion Cohort

# In[ ]:


sarcoid_without_exclusions = sarcoid_with_exclusions.filter(~pl.col('has_exclusion_condition')).select(['person_id', 'sarcoid_n', 'vocab', 'sarcoid_1'])


# In[ ]:


print(sns.color_palette("colorblind").as_hex())
sns.color_palette("colorblind")


# In[ ]:


# Assuming you're working with both dataframes
df = sarcoid_df.to_pandas()
df_without_exclusions = sarcoid_without_exclusions.to_pandas()

# Create a figure with subplots
fig = plt.figure(figsize=(10, 5))

# 1. Histogram of sarcoid_n with percentage annotation
plt.subplot(1, 2, 1)
sns.histplot(data=df, x='sarcoid_n', binwidth=1, color='#0173b2', kde=True)
plt.xlim(0, 75)
plt.title('Distribution of Sarcoid Counts per Person', fontsize=12)
plt.xlabel('Sarcoid Count')
plt.ylabel('Frequency')

# Calculate percentage of data within the x-limits
total_rows = len(df)
rows_in_range = len(df[df['sarcoid_n'] <= 75])
percentage = (rows_in_range / total_rows) * 100

# Add percentage annotation to top right
plt.annotate(f'{percentage:.1f}% of data', 
             xy=(0.95, 0.95), 
             xycoords='axes fraction',
             ha='right', va='top',
             bbox=dict(boxstyle='round', fc='white', alpha=0.8))

# Remove top and right spines
sns.despine(ax=plt.gca())

plt.gca().tick_params(which='major', width=1, left=True, bottom=True, color='#999', grid_alpha=0.2)

# 2. Bar plot of counts per year from sarcoid_1 with improved x-axis
plt.subplot(1, 2, 2)
# Extract year from sarcoid_1 dates
df_without_exclusions['year'] = pd.to_datetime(df_without_exclusions['sarcoid_1']).dt.year

# Create empty dataframe with all years from 1982 to 2023
all_years = pd.DataFrame({'year': range(1982, 2024)})
year_counts = df_without_exclusions['year'].value_counts().reset_index()
year_counts.columns = ['year', 'count']

# Merge to ensure all years are represented
complete_years = all_years.merge(year_counts, on='year', how='left').fillna(0)
complete_years = complete_years.sort_values('year')

# Plot the bar chart
plt.bar(complete_years['year'], complete_years['count'], color='#0173b2', alpha=0.7)
plt.title('Sarcoidosis Cases by Year', fontsize=12)
plt.xlabel('Year')
plt.ylabel('Number of Cases')

# Set up x-axis ticks
plt.xticks(np.arange(1985, 2024, 5), rotation=45)
plt.xlim(1981.5, 2023.5)

# Remove top and right spines
sns.despine(ax=plt.gca())
plt.gca().tick_params(which='major', width=1, left=True, bottom=True, color='#999', grid_alpha=0.2)

plt.minorticks_on()
plt.gca().tick_params(which='minor', width=0.5, bottom=True, color='#999', grid_alpha=0.2)

plt.tight_layout()
plt.show()


# In[ ]:


df = sarcoid_df.to_pandas()

# UpSet plot of vocab
# First, create binary indicators for each vocab category
vocab_categories = ["9", "10", "SNO"]

# Create dictionary of sets for each category
vocab_dict = {}
for category in vocab_categories:
    # Find patients with this category in their vocab
    vocab_dict[category] = set(df[df['vocab'].str.contains(category, na=False)]['person_id'])

# Create memberships list for UpSet plot
memberships = []
for _, row in df.iterrows():
    if pd.isna(row['vocab']):
        memberships.append([])
        continue
        
    patient_vocabs = []
    for category in vocab_categories:
        if category in str(row['vocab']):
            patient_vocabs.append(category)
    memberships.append(patient_vocabs)

membership_data = from_memberships(memberships)

# Create the UpSet plot

fig = plt.figure(figsize=(5, 5))
plot(membership_data, subset_size="count", sort_by='cardinality', 
     sort_categories_by='-cardinality', 
     max_subset_rank=20,
     show_counts=True,
    )
plt.title('')
plt.show()


# # Add Covariates

# In[ ]:


all_sarcoid_cases = sarcoid_with_exclusions.select(['person_id', 'sarcoid_n', 'vocab', 'sarcoid_1'])


# In[ ]:


sarcoid_cases_without_exclusions = sarcoid_without_exclusions


# In[ ]:


all_person_ids_str = ', '.join(all_sarcoid_cases['person_id'].cast(pl.Utf8).to_list())


# In[ ]:


test_q = f"""
SELECT MAX(observation_date)
FROM {dataset}.observation
"""

polars_gbq(test_q)


# In[ ]:


# Heart rate rhythm, Heart rate, Blood pressure panel, Adult Waist Circumference Protocol, PhenX - hip circumference protocol 020801
# Body height, Body weight, Diastolic blood pressure, Systolic blood pressure, Body mass index (BMI) [Ratio], Body weight Measured --pre pregnancy
vitals_exclusion = [3022318, 3027018, 3031203, 40759207, 40765148, 3036277, 3025315, 3012888, 3004249, 3038553, 3022281]
vitals_exclusion_str = ', '.join(map(str, vitals_exclusion))

# unique person_id w/ ICD or SNOMED in source_ values in observation/condition_occurrence
distinct_participants_cte = f"""
WITH distinct_participants AS (
    SELECT DISTINCT o.person_id
    FROM {dataset}.observation AS o
    JOIN {dataset}.concept AS c 
        ON o.observation_source_value = c.concept_code
        OR o.observation_source_concept_id = c.concept_id
    WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM', 'SNOMED')

    UNION DISTINCT
    
    SELECT DISTINCT co.person_id
    FROM {dataset}.condition_occurrence AS co
    JOIN {dataset}.concept AS c 
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
# ICD only: 350,427
# ICD + SNOMED: 358,180


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
    {dataset}.person p
LEFT JOIN
    {dataset}.concept p_gender_concept 
        ON p.gender_concept_id = p_gender_concept.concept_id 
LEFT JOIN
    {dataset}.concept p_race_concept 
        ON p.race_concept_id = p_race_concept.concept_id 
LEFT JOIN
    {dataset}.concept p_ethnicity_concept 
        ON p.ethnicity_concept_id = p_ethnicity_concept.concept_id 
LEFT JOIN
    {dataset}.concept p_sex_at_birth_concept 
        ON p.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id  
LEFT JOIN
    {dataset}.death d
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
        {dataset}.observation AS obs
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
        {dataset}.observation AS obs
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
        {dataset}.observation AS obs
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
        {dataset}.observation AS obs
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
        {dataset}.observation AS obs
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


# total (ICD): 350,427
# total (ICD/SNOMED): 358,180
demographics_df.null_count()


# In[ ]:


demographics_df = demographics_df.to_pandas()


# In[ ]:


demographics_df[['race', 'ethnicity']].value_counts()


# In[ ]:


demographics_df = pl.from_pandas(demographics_df)


# In[ ]:


all_sarcoid_cases = all_sarcoid_cases.with_columns(
    pl.lit(1).alias('case')
)
sarcoid_cases_without_exclusions = sarcoid_cases_without_exclusions.with_columns(
    pl.lit(1).alias('case')
)


# In[ ]:


all_sarcoid_cases = demographics_df.join(
    all_sarcoid_cases,
    on='person_id',
    how='left'
)
sarcoid_cases_without_exclusions = demographics_df.join(
    sarcoid_cases_without_exclusions,
    on='person_id',
    how='left'
)


# In[ ]:


all_sarcoid_cases = all_sarcoid_cases.with_columns(
    (
        (pl.col("sarcoid_1") - pl.col("dob")).dt.total_days() / 365.25
    ).alias("age_at_diagnosis")
)
sarcoid_cases_without_exclusions = sarcoid_cases_without_exclusions.with_columns(
    (
        (pl.col("sarcoid_1") - pl.col("dob")).dt.total_days() / 365.25
    ).alias("age_at_diagnosis")
)


# In[ ]:


all_sarcoid_cases = all_sarcoid_cases.with_columns(
    pl.col('case').fill_null(0)
)
sarcoid_cases_without_exclusions = sarcoid_cases_without_exclusions.with_columns(
    pl.col('case').fill_null(0)
)


# In[ ]:


demographics_df = demographics_df.join(procedure_df, on='person_id', how='left')


# In[ ]:


all_sarcoid_cases = all_sarcoid_cases.join(procedure_df, on='person_id', how='left')
sarcoid_cases_without_exclusions = sarcoid_cases_without_exclusions.join(procedure_df, on='person_id', how='left')


# In[ ]:


all_sarcoid_cases.write_csv(f'{bucket}/data/cohorts/all_sarcoid_cases.csv')
sarcoid_cases_without_exclusions.write_csv(f'{bucket}/data/cohorts/sarcoid_cases_without_exclusions.csv')


# In[ ]:


all_sarcoid_cases_pd = all_sarcoid_cases.to_pandas()
sarcoid_cases_without_exclusions_pd = sarcoid_cases_without_exclusions.to_pandas()


# In[ ]:


# all_sarcoid_cases_pd = pd.read_csv(f'{bucket}/data/cohorts/all_sarcoid_cases.csv')
# sarcoid_cases_without_exclusions_pd = pd.read_csv(f'{bucket}/data/cohorts/sarcoid_cases_without_exclusions.csv')


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

    # Sarcoid finginds
    conditions = ['cxr', 'ct_chest', 'pulm_bx', 'non_pulm_bx', 'pft', 'ecg', 'leukopenia', 
                  'thrombocytopenia', 'pbe', 'hypercalcemia', 'elevated_ace', 'elevated_ap', 
                  'low_25_d', 'high_1_25_d', 'vit_d_deficiency', 'high_a1at', 'elevated_sil2']

    for cond in conditions:
        # Check if condition case date exists (not null) - 1 if ever had condition, 0 if never
        df[cond] = df[f'{cond}_1'].notna().astype(int)

    return df

# Apply the function
all_sarcoid_cases_pd = update_covariates(all_sarcoid_cases_pd, sarcoid=True)
sarcoid_cases_without_exclusions_pd = update_covariates(sarcoid_cases_without_exclusions_pd, sarcoid=True)


# In[ ]:


demographics_df = demographics_df.to_pandas()
demographics_df = update_covariates(demographics_df, sarcoid=False)


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
    # Sarcoid conditions
    conditions_percentages = {col: calc_percentage(col) for col in ['cxr', 'ct_chest', 'pulm_bx', 
                                                                    'non_pulm_bx', 'pft', 'ecg', 'leukopenia', 
                                                                    'thrombocytopenia', 'pbe', 'hypercalcemia', 
                                                                    'elevated_ace', 'elevated_ap',
                                                                    'low_25_d', 'high_1_25_d', 
                                                                    'vit_d_deficiency', 'high_a1at', 
                                                                    'elevated_sil2']}

    # Building the results dictionary
    if sarcoid == True:
        results = {
            'Total': len(df),
            'End of Study Age, median (IQR)': f"{age_desc['median']} ({age_desc.iloc[1]}, {age_desc.iloc[2]})",
            'Sarcoid Dx Age, median (IQR)': f"{dx_age_desc['median']} ({dx_age_desc.iloc[1]}, {dx_age_desc.iloc[2]})",
            **{f"{key.capitalize()} (n (%))": value for key, value in demographics_percentages.items()},
            **{f"{key.capitalize()} (n (%))": value for key, value in ses_percentages.items()},
            **{f"{key.capitalize()} (n (%))": value for key, value in conditions_percentages.items()}
        }    
    else:
        results = {
            'Total': len(df),
            'End of Study Age, median (IQR)': f"{age_desc['median']} ({age_desc.iloc[1]}, {age_desc.iloc[2]})",
            'Sarcoid Dx Age, median (IQR)': f"NA",
            **{f"{key.capitalize()} (n (%))": value for key, value in demographics_percentages.items()},
            **{f"{key.capitalize()} (n (%))": value for key, value in ses_percentages.items()},
            **{f"{key.capitalize()} (n (%))": value for key, value in conditions_percentages.items()}
        }           

    return pd.DataFrame(results, index=[0])

total_ehr_cohort = describe_group(demographics_df, sarcoid=False)
        
all_summaries = pd.concat([
    describe_group(all_sarcoid_cases_pd[(all_sarcoid_cases_pd["case"]==1) & (all_sarcoid_cases_pd["sarcoid_n"]==1)]), 
    describe_group(all_sarcoid_cases_pd[(all_sarcoid_cases_pd["case"]==1) & (all_sarcoid_cases_pd["sarcoid_n"]>=1)]), 
    describe_group(all_sarcoid_cases_pd[(all_sarcoid_cases_pd["case"]==1) & (all_sarcoid_cases_pd["sarcoid_n"]>=2)]), 
    describe_group(sarcoid_cases_without_exclusions_pd[(sarcoid_cases_without_exclusions_pd["case"]==1) & (all_sarcoid_cases_pd["sarcoid_n"]==1)]), 
    describe_group(sarcoid_cases_without_exclusions_pd[(sarcoid_cases_without_exclusions_pd["case"]==1) & (all_sarcoid_cases_pd["sarcoid_n"]>=1)]), 
    describe_group(sarcoid_cases_without_exclusions_pd[(sarcoid_cases_without_exclusions_pd["case"]==1) & (all_sarcoid_cases_pd["sarcoid_n"]>=2)]), 
    ],
    keys=["All Sarcoid Cases (1 Phecode Only)", 
        "All Sarcoid Cases (≥1 Phecode)", 
        "All Sarcoid Cases (≥2 Phecodes)",
        "Sarcoid Cases Without Exclusions (Phecode Only)",
        "Sarcoid Cases Without Exclusions (≥1 Phecode)",
        "Sarcoid Cases Without Exclusions (≥2 Phecodes)"]
)

# Then concatenate with total_ehr_cohort using a single-level index
table_1 = pd.concat([
    pd.DataFrame(total_ehr_cohort).assign(category='Total Cohort'),
    pd.DataFrame(all_summaries).assign(category='')
]).set_index('category', append=True)


# In[ ]:


transposed_table = table_1.transpose()
transposed_table


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


total_ehr_df.head()


# In[ ]:


all_sarcoid_cases_pd["case"] = 1
sarcoid_cases_without_exclusions_pd["case"] = 1


# In[ ]:


# Merge and fillna for all_sarcoid_cases
all_sarcoid_cases_cohort = demographics_df.merge(
    all_sarcoid_cases_pd[['person_id', 'case']],
    on='person_id',
    how='left'
)
all_sarcoid_cases_cohort['case'] = all_sarcoid_cases_cohort['case'].fillna(0)
all_sarcoid_cases_cohort = all_sarcoid_cases_cohort[all_sarcoid_cases_cohort['female'] != -9]
all_sarcoid_cases_cohort['female'] = all_sarcoid_cases_cohort['female'].map({1: 0, 0: 1})

# Save to CSV
all_sarcoid_cases_cohort.to_csv('all_sarcoid_cases_cohort.csv', index=False)

# Merge and fillna for sarcoid_cases_without_exclusions
sarcoid_cases_without_exclusions_df = demographics_df.merge(
    sarcoid_cases_without_exclusions_pd[['person_id', 'case']],
    on='person_id',
    how='left'
)
sarcoid_cases_without_exclusions_df['case'] = sarcoid_cases_without_exclusions_df['case'].fillna(0)
sarcoid_cases_without_exclusions_df = sarcoid_cases_without_exclusions_df[sarcoid_cases_without_exclusions_df['female'] != -9]
sarcoid_cases_without_exclusions_df['female'] = sarcoid_cases_without_exclusions_df['female'].map({1: 0, 0: 1})

# Save to CSV
sarcoid_cases_without_exclusions_df.to_csv('sarcoid_cases_without_exclusions_cohort.csv', index=False)


# In[ ]:


sarcoid_cases_without_exclusions_df.value_counts('female')


# In[ ]:


all_sarcoid_phewas = PheWAS(
    phecode_version="X",
    phecode_count_csv_path="phecode_X_counts.csv",
    cohort_csv_path="all_sarcoid_cases_cohort.csv",
    sex_at_birth_col="female",
    male_as_one=True,
    covariate_cols=["end_of_study_age", "female", "white", "black", "asian", "hispanic"],
    independent_variable_of_interest="case",
    min_cases=50,
    min_phecode_count=2,
    output_file_name="all_sarcoid_phewas_results.csv"
)
all_sarcoid_phewas.run()


# In[ ]:


sarcoid_without_exclusions_phewas = PheWAS(
    phecode_version="X",
    phecode_count_csv_path="phecode_X_counts.csv",
    cohort_csv_path="sarcoid_cases_without_exclusions_cohort.csv",
    sex_at_birth_col="female",
    male_as_one=True,
    covariate_cols=["end_of_study_age", "female", "white", "black", "asian", "hispanic"],
    independent_variable_of_interest="case",
    min_cases=50,
    min_phecode_count=2,
    output_file_name="sarcoid_without_exclusions_phewas_results.csv"
)
sarcoid_without_exclusions_phewas.run()


# In[ ]:


p = Plot("all_sarcoid_phewas_results.csv")
p.manhattan(label_values="p_value", label_count=25, save_plot=False)


# In[ ]:


p = Plot("sarcoid_without_exclusions_phewas_results.csv")
p.manhattan(label_values="p_value", label_count=25, save_plot=False)


# In[ ]:


pgrm_bucket = '{bucket or my_bucket}'
acaf_person_id_df = pl.read_csv(f'{pgrm_bucket}/data/hail/pgrm_filtered_acaf_plink.fam', 
                                separator="\t", 
                                new_columns=["A","B","C","D","E","F"])

acaf_person_id_df = acaf_person_id_df.select('B')


# In[ ]:


all_sarcoid_cases = pl.read_csv('all_sarcoid_cases.csv', try_parse_dates=True)
sarcoid_cases_without_exclusions = pl.read_csv('sarcoid_cases_without_exclusions.csv', try_parse_dates=True)


# In[ ]:


# Count people who are both in sarcoid_cases_without_exclusions with case=1 and in acaf_person_id_df
overlap_count = (
    sarcoid_cases_without_exclusions
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


# # Zip3
exclusion_cols = [
    'tb_2', 'tb_skin_2', 'ntm_2', 'histo_2', 'blasto_2', 'cocci_2', 
    'asp_2', 'ild_2', 'hp_2', 'pneumoconiosis_2', 'vasc_2', 'cvid_2', 
    'igg4_2', 'lymphoma_2', 'lg_2', 'lung_ca_2', 'lch_2', 'ibd_2', 
    'pbc_2', 'aih_2', 'hiv_2', 'drug_2'
]# Map of column names to more readable labels (remove '_2' suffix)
exclusion_labels = {col: col.replace('_2', '') for col in exclusion_cols}exclusion_labels = {
    'tb_2': 'tb',
    'tb_skin_2': 'abnormal_tb_skin_test',
    'ntm_2': 'ntm',
    'histo_2': 'histo',
    'blasto_2': 'blasto',
    'cocci_2': 'cocci',
    'asp_2': 'aspergillus',
    'ild_2': 'ild',
    'hp_2': 'hp',
    'pneumoconiosis_2': 'pneumoconiosis',
    'vasc_2': 'vasculitis',
    'cvid_2': 'cvid',
    'igg4_2': 'igg4',
    'lymphoma_2': 'lymphoma',
    'lg_2': 'lymphomatoid granulomatosis',
    'lung_ca_2': 'lung_ca',
    'lch_2': 'langerhans_cell_histiocytosis',
    'ibd_2': 'ibd',
    'pbc_2': 'pbc',
    'aih_2': 'autoimmune_hepatitis',
    'hiv_2': 'hiv',
    'drug_2': 'drug'
}# Create the 'exclusion' column with a list of labels where the date is not null
sarcoid_with_exclusions = sarcoid_with_exclusions.with_columns([
    pl.concat_list([
        pl.when(pl.col(col).is_not_null()).then(pl.lit(exclusion_labels[col])).otherwise(pl.lit(None))
        for col in exclusion_cols
    ]).alias('exclusion')
])sarcoid_with_exclusions = sarcoid_with_exclusions.with_columns([
    pl.col('exclusion').list.eval(pl.element().filter(pl.element().is_not_null())).alias('exclusion')
])# sarcoid_with_exclusions = sarcoid_with_exclusions.select(['person_id', 'sarcoid_1', 'exclusion'])total_cohort_q = f"""
        WITH combined AS (
            
            -- Codes from observation (source_value)
            SELECT o.person_id
            FROM {dataset}.observation AS o
            JOIN {dataset}.concept AS c ON o.observation_source_value = c.concept_code
            WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM', 'SNOMED')

            UNION ALL

            -- Codes from observation (source_concept_id)
            SELECT o.person_id
            FROM {dataset}.observation AS o
            JOIN {dataset}.concept AS c ON o.observation_source_concept_id = c.concept_id
            WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM', 'SNOMED')

            UNION ALL

            -- Codes from condition_occurrence (source_value)
            SELECT co.person_id
            FROM {dataset}.condition_occurrence AS co
            JOIN {dataset}.concept AS c ON co.condition_source_value = c.concept_code
            WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM', 'SNOMED')

            UNION ALL

            -- Codes from condition_occurrence (source_concept_id)
            SELECT co.person_id
            FROM {dataset}.condition_occurrence AS co
            JOIN {dataset}.concept AS c ON co.condition_source_concept_id = c.concept_id
            WHERE c.vocabulary_id IN ('ICD9CM', 'ICD10CM', 'SNOMED')            
        ),
        zip AS (
            SELECT
                observation.person_id,
                value_as_string as zip_code
            FROM
                {dataset}.observation observation
            WHERE
                observation_source_concept_id = 1585250
        )
  
        SELECT
            DISTINCT c.person_id,
            zip_code
        FROM combined AS c
        JOIN zip AS z ON c.person_id = z.person_id
"""

zip_df = polars_gbq(total_cohort_q)zip_df = zip_df.join(
    sarcoid_with_exclusions,
    on='person_id',
    how='left'
)# Trim zip_code to first three digits and create a new column zip_3
zip_df = zip_df.with_columns([
    pl.col("zip_code").str.slice(0, 3).alias("zip_3")
]).drop('zip_code')

zip_df = zip_df.with_columns([
    pl.when(pl.col("sarcoid_1").is_not_null())
    .then(pl.lit(1))
    .otherwise(pl.lit(0))
    .alias("case")
])

# Fill NA values in exclusion column with empty list
zip_df = zip_df.with_columns([
    pl.when(pl.col("exclusion").is_null())
    .then(pl.lit([]))
    .otherwise(pl.col("exclusion"))
    .alias("exclusion")
])
# Select the desired columns
result_df = zip_df.select(['case', 'sarcoid_1', 'exclusion', 'zip_3'])

result_df = result_df.with_columns([
    pl.col("exclusion").list.join(", ").alias("exclusion_str")
])

# Drop the original list column
result_df = result_df.drop("exclusion")

# Now save to TSV
result_df.write_csv("sarcoid_zip3.tsv", separator="\t")# Create sarcoid_zip3 DataFrame (cases by zip code)
sarcoid_zip3 = (
    result_df
    .filter(pl.col("case") == 1)
    .group_by("zip_3")
    .agg(pl.count().alias("case_count"))
    .sort("zip_3")
)

# Create control_zip3 DataFrame (controls by zip code)
control_zip3 = (
    result_df
    .filter(pl.col("case") == 0)
    .group_by("zip_3")
    .agg(pl.count().alias("control_count"))
    .sort("zip_3")
)sarcoid_zip3.write_csv('sarcoid_zip3_trim.tsv', separator='\t')
control_zip3.write_csv('control_zip3_trim.tsv', separator='\t')