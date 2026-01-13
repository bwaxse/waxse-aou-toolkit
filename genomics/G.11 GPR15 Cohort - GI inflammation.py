#!/usr/bin/env python
# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><ul class="toc-item"><li><span><a href="#README" data-toc-modified-id="README-0.1"><span class="toc-item-num">0.1&nbsp;&nbsp;</span>README</a></span></li></ul></li><li><span><a href="#Initial-setup" data-toc-modified-id="Initial-setup-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Initial setup</a></span></li><li><span><a href="#Gather-Homozygote/Heterozygote-Counts-(Pathfinder)" data-toc-modified-id="Gather-Homozygote/Heterozygote-Counts-(Pathfinder)-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Gather Homozygote/Heterozygote Counts (Pathfinder)</a></span><ul class="toc-item"><li><span><a href="#Read-in-Files-for-Analysis" data-toc-modified-id="Read-in-Files-for-Analysis-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>Read in Files for Analysis</a></span></li><li><span><a href="#Define-hyperparameters" data-toc-modified-id="Define-hyperparameters-2.2"><span class="toc-item-num">2.2&nbsp;&nbsp;</span>Define hyperparameters</a></span></li><li><span><a href="#Run" data-toc-modified-id="Run-2.3"><span class="toc-item-num">2.3&nbsp;&nbsp;</span>Run</a></span></li><li><span><a href="#Generate-participant-lists" data-toc-modified-id="Generate-participant-lists-2.4"><span class="toc-item-num">2.4&nbsp;&nbsp;</span>Generate participant lists</a></span></li><li><span><a href="#Variant-Counts" data-toc-modified-id="Variant-Counts-2.5"><span class="toc-item-num">2.5&nbsp;&nbsp;</span>Variant Counts</a></span></li></ul></li><li><span><a href="#Generate-Cohort-(PheTK)---check-related" data-toc-modified-id="Generate-Cohort-(PheTK)---check-related-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Generate Cohort (PheTK) - check related</a></span><ul class="toc-item"><li><span><a href="#Generate-Cohorts" data-toc-modified-id="Generate-Cohorts-3.1"><span class="toc-item-num">3.1&nbsp;&nbsp;</span>Generate Cohorts</a></span></li><li><span><a href="#Compare-Pathfinder-and-PheTK" data-toc-modified-id="Compare-Pathfinder-and-PheTK-3.2"><span class="toc-item-num">3.2&nbsp;&nbsp;</span>Compare Pathfinder and PheTK</a></span></li><li><span><a href="#Add-covariate-data" data-toc-modified-id="Add-covariate-data-3.3"><span class="toc-item-num">3.3&nbsp;&nbsp;</span>Add covariate data</a></span></li><li><span><a href="#Add-IBD-data" data-toc-modified-id="Add-IBD-data-3.4"><span class="toc-item-num">3.4&nbsp;&nbsp;</span>Add IBD data</a></span><ul class="toc-item"><li><span><a href="#Initialize" data-toc-modified-id="Initialize-3.4.1"><span class="toc-item-num">3.4.1&nbsp;&nbsp;</span>Initialize</a></span></li><li><span><a href="#IBD-by-AoUQueries-(ICD-and-SNOMED)" data-toc-modified-id="IBD-by-AoUQueries-(ICD-and-SNOMED)-3.4.2"><span class="toc-item-num">3.4.2&nbsp;&nbsp;</span>IBD by AoUQueries (ICD and SNOMED)</a></span></li><li><span><a href="#IBD-by-PheTK-Method-(Phecode-Only)" data-toc-modified-id="IBD-by-PheTK-Method-(Phecode-Only)-3.4.3"><span class="toc-item-num">3.4.3&nbsp;&nbsp;</span>IBD by PheTK Method (Phecode Only)</a></span></li><li><span><a href="#Add-Age-at-IBD" data-toc-modified-id="Add-Age-at-IBD-3.4.4"><span class="toc-item-num">3.4.4&nbsp;&nbsp;</span>Add Age at IBD</a></span></li></ul></li></ul></li><li><span><a href="#Assess-Genetic-Ancestry" data-toc-modified-id="Assess-Genetic-Ancestry-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Assess Genetic Ancestry</a></span></li></ul></div>

# ## README ##
# 
# This notebook creates a cohort using hail starting from a variant list.
# 
# __For VM configurations__, __if you are only run PheWAS and Manhattan plot__, a __standard VM__ with __96 CPUs, 86.4 GB RAM (3.51 USD/hour)__ is sufficient.__dataproc VM__ with __main instance__ having __64 CPUs, 240 GB RAM__, and __10 workers (2 standard, 8 preemptible)__ having __4 CPUs, 15GB RAM (5.76 USD/hour)__ is suggested. Storage just needs to be enough to meet requirement, e.g., ~ 150GB.

# # Initial setup

# In[ ]:


# !pip install polars && pip install connectorx && pip install adjustText


# In[ ]:


# !pip install PheTK


# In[ ]:


pip show PheTK | grep Version


# In[ ]:


from google.cloud import bigquery
import pandas as pd
import polars as pl
import numpy as np
import scipy
import matplotlib as mpl 
import matplotlib.pyplot as plt
import seaborn as sns
import os
import subprocess
import math
# import hail as hl
import logging
from datetime import datetime
from typing import Dict, List, Optional, Union
from jinja2 import Template, Environment


# In[ ]:


bucket = os.getenv('WORKSPACE_BUCKET')
print(bucket)
CDR = os.getenv('WORKSPACE_CDR')
print(CDR)


# In[ ]:


# This line allows for the plots to be displayed inline in the Jupyter notebook
get_ipython().run_line_magic('matplotlib', 'inline')

sns.set(style="ticks",font_scale=1)


# In[ ]:


pl.Config.set_fmt_str_lengths(128)

# Set the row limit to a higher value
pl.Config.set_tbl_rows(50)

# show all columns in pandas
pd.set_option("display.max_columns", None)

# show full column width
pd.set_option('display.max_colwidth', 100)


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


def initialize():
    """
    Initialize Hail
    :param query: BigQuery SQL query
    :return: polars dataframe
    """
    start = datetime.now()
    start

    hl.init(idempotent = True)
    hl.default_reference('GRCh38')
    logging.getLogger('hail').setLevel(logging.ERROR)


# In[ ]:


def filter_exome(mt, interval):
    """
    Filter exome matrix table to gene interval of interest
    :param mt: Hail matrix table
    :param interval: gene interval
    :return: Hail matrix table of filtered exome, with n_alt column
    """
    
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(x) for x in interval])
    mt = mt.annotate_entries(n_alt = mt.GT.n_alt_alleles())
    mt = mt.drop(mt.filters, mt.variant_qc, mt.info)

    return mt 


# In[ ]:


def process_variants(variants):
    """
    convert Pandas dataframe uploaded_var_pd to Hail mt
    :param variants: pandas df of variants; key by Name if column with transcripts/c-dots 
                     already in ENST00000332972.9:c.49C>T format (or key by vid)
    :return: Hail matrix table of filtered exome, with n_alt column
    """

    variants_ht = hl.Table.from_pandas(variants, key = 'vid')
    
    return variants_ht


# In[ ]:


def annotate_exome(vat_ht, interval, variants_ht, filtered_exome):
    """
    Annotates Hail mt using VAT
    :param vat_ht: Hail mt of VAT
    :param interval: gene interval
    :param variants: Hail mt of variants
    :param filtered_exome: filtered exome
    :return: Hail matrix table of annotated exome
    """
    
    vat_cols = ['vid', 'transcript', 'gene_symbol', 'gnomad_all_af', 
                'aa_change', 'consequence', 'dna_change_in_transcript', 
                'exon_number', 'intron_number', 'genomic_location', 'dbsnp_rsid', 
                'is_canonical_transcript', 'revel', 'splice_ai_acceptor_gain_score',
               'splice_ai_acceptor_loss_score', 'splice_ai_donor_gain_score',
               'splice_ai_donor_loss_score']
    
    # filter variant annotation table (VAT) to gene interval of interest
    vat_ht = hl.filter_intervals(vat_ht, [hl.parse_locus_interval(x) for x in interval])

    # filter VAT by canonical transcript
    vat_ht = vat_ht.filter(vat_ht.is_canonical_transcript)
    
    # subset VAT to list of columns
    vat_ht = vat_ht.select(*vat_cols)
 
    # annotate VAT with column "uploaded_var", which contains information from "src" column
    if variants_ht is not None:
        # Create a lookup dictionary
        variant_lookup = variants_ht.to_pandas().set_index('vid')['src'].to_dict()

        # Annotate using the lookup
        vat_ht = vat_ht.annotate(
            uploaded_var = hl.literal(variant_lookup).get(vat_ht.vid)
        )
#     # annotate VAT with column "uploaded_var", which contains information from "src" column
#     if variants_ht is not None:
#         variants_for_join = variants_ht.select('src').key_by('vid')
#         vat_ht_keyed = vat_ht.key_by('vid')
        
#         # annotate
#         vat_ht_annotated = vat_ht_keyed.annotate(uploaded_var = variants_for_join[vat_ht_keyed.vid].src)
                
#         # Re-key back to original
#         vat_ht = vat_ht_annotated.key_by('locus', 'alleles')
#     else:
#         vat_ht = vat_ht.annotate(uploaded_var = hl.missing(hl.tstr))


    # annotate (filtered) exome matrix table with VAT
    ann_exome_mt = filtered_exome.annotate_rows(annotations = vat_ht[filtered_exome.row_key])

    # add allele_freq column from gnomad
    ann_exome_mt = ann_exome_mt.annotate_rows(
        allele_freq = hl.if_else(
            hl.is_defined(ann_exome_mt.annotations.gnomad_all_af),
            ann_exome_mt.annotations.gnomad_all_af,
            0
        )
    )
    
    # add pLOF column - possible values are 0 and 1
    lof_consequence_set = hl.set([
                            'start_lost', 'stop_lost', 
                            'stop_gained','frameshift_variant'
                        ])
    
    ann_exome_mt = ann_exome_mt.annotate_rows(
        pLOF = hl.if_else(
            (ann_exome_mt.allele_freq <= MAF_freq)
             &(lof_consequence_set.contains(ann_exome_mt.annotations.consequence)),
            1,
            0
        )
    )
    
    # add clinvar column - possible values are 0 and 1
    ann_exome_mt = ann_exome_mt.annotate_rows(
        clinvar = hl.if_else(
                    (hl.is_defined(ann_exome_mt.annotations.uploaded_var))
                     &(ann_exome_mt.annotations.uploaded_var == "clinvar"),
                    1,
                    0
        )
    )

    # add lab column - possible values are 0 and 1
    ann_exome_mt = ann_exome_mt.annotate_rows(
        lab = hl.if_else(
                    (hl.is_defined(ann_exome_mt.annotations.uploaded_var))
                     &(ann_exome_mt.annotations.uploaded_var == "lab"),
                    1,
                    0
        )
    )
    
    # add gnomad_plof column - possible values are 0 and 1
    ann_exome_mt = ann_exome_mt.annotate_rows(
        gnomad_plof = hl.if_else(
                    (hl.is_defined(ann_exome_mt.annotations.uploaded_var))
                     &(ann_exome_mt.annotations.uploaded_var == "gnomad_plof"),
                    1,
                    0
        )
    )
    
    # add gnomad_other column - possible values are 0 and 1
    ann_exome_mt = ann_exome_mt.annotate_rows(
        gnomad_other = hl.if_else(
                    (hl.is_defined(ann_exome_mt.annotations.uploaded_var))
                     &(ann_exome_mt.annotations.uploaded_var == "gnomad_other"),
                    1,
                    0
        )
    )

    # add gnomad_missense column - possible values are 0 and 1
    ann_exome_mt = ann_exome_mt.annotate_rows(
        gnomad_missense = hl.if_else(
                    (hl.is_defined(ann_exome_mt.annotations.uploaded_var))
                     &(ann_exome_mt.annotations.uploaded_var == "gnomad_missense"),
                    1,
                    0
        )
    )
    
    return ann_exome_mt


# In[ ]:


def pLOF_variants (ann_exome_mt):
    """
    Filter the interested exome matrix table to only include pLOF variants
    :param ann_exome_mt: Hail matrix table of annotated exome
    :return: pLOF mt
    """
    plof_mt = ann_exome_mt.filter_rows(ann_exome_mt.pLOF == 1)
            
    return plof_mt


# In[ ]:


def clinvar_variants (ann_exome_mt):
    """
    Filter the interested exome matrix table to only include clinvar variants
    :param ann_exome_mt: Hail matrix table of annotated exome
    :return: clinvar mt
    """
    clinvar_mt = ann_exome_mt.filter_rows(ann_exome_mt.clinvar == 1)
            
    return clinvar_mt


# In[ ]:


def lab_variants (ann_exome_mt):
    """
    Filter the interested exome matrix table to only include lab variants
    :param ann_exome_mt: Hail matrix table of annotated exome
    :return: lab mt
    """
    lab_mt = ann_exome_mt.filter_rows(ann_exome_mt.lab == 1)
            
    return lab_mt


# In[ ]:


def gnomad_plof_variants (ann_exome_mt):
    """
    Filter the interested exome matrix table to only include gnomad pLOF variants
    :param ann_exome_mt: Hail matrix table of annotated exome
    :return: gnomad_plof mt
    """
    gnomad_plof_mt = ann_exome_mt.filter_rows(ann_exome_mt.gnomad_plof == 1)
            
    return gnomad_plof_mt


# In[ ]:


def gnomad_other_variants (ann_exome_mt):
    """
    Filter the interested exome matrix table to only include gnomad other variants
    :param ann_exome_mt: Hail matrix table of annotated exome
    :return: gnomad_other mt
    """
    gnomad_other_mt = ann_exome_mt.filter_rows(ann_exome_mt.gnomad_other == 1)
            
    return gnomad_other_mt


# In[ ]:


def gnomad_missense_variants (ann_exome_mt):
    """
    Filter the interested exome matrix table to only include gnomad missense variants
    :param ann_exome_mt: Hail matrix table of annotated exome
    :return: gnomad_missense mt
    """
    gnomad_missense_mt = ann_exome_mt.filter_rows(ann_exome_mt.gnomad_missense == 1)
            
    return gnomad_missense_mt


# In[ ]:


def combine_components(components):
    """
    Takes in a list of component matrix tables and filters to remove duplicate variants
    :param components: component matrix tables
    :return: combined exome mt
    """
    
    int_exome_mt = components[0]
    for component in components[1:]:
        int_exome_mt = int_exome_mt.union_rows(component)
        
    int_exome_mt = int_exome_mt.distinct_by_row() # filter to remove duplicate variants
    
    return int_exome_mt


# In[ ]:


def find_biallelic_participants(int_exome_mt):
    """
    Searches AoU for participants with 2 alternate alleles among variants
    :return: sorted_participant_info, all participant + variant information for all sorted participants
    :return: signal_set, set of participant IDs ('s' values) that the program finds to be biallelic for a variant of interest
    """
    print("searching biallelic patients now!")

    signal_table = int_exome_mt.entries() # entries step
    print("converted matrix table to table")
    
    # aggregated participant list
    geno_agg_ht = signal_table.group_by(signal_table.s).aggregate(n_alt_all = hl.agg.sum(signal_table.n_alt))
    signal_set = geno_agg_ht.filter(geno_agg_ht.n_alt_all >= 2)
    
    print("Biallelic participants:")
    signal_set.show()
    
    cnt = signal_set.count()
    print(f"Number of biallelic participants: {cnt}")
    
#     # estimated prevalence
#     div_cnt = 245394/cnt
#     print(f"Estimated prevalence based on this feature is 1 in {div_cnt}")
    
    # entries table (variant + participant information)
    participant_info = signal_table.filter(signal_table.n_alt > 0)
    participant_info = participant_info.filter(hl.is_defined(signal_set.index(participant_info['s'])))
    
    # @Kate DELETE THIS TO BE FASTER... 
    sorted_participant_info = participant_info.order_by(participant_info.s)

    sorted_participant_info.show(100)
    
    return sorted_participant_info, signal_set


# In[ ]:


def find_monallelic_participants(int_exome_mt):
    """
    Searches AoU for participants with 1 (or more) alternate alleles among likely pathogenic variants
    :return: sorted_participant_info, all participant + variant information for all sorted participants
    :return: signal_set, set of participant IDs ('s' values) that the program finds to be monallelic for a variant of interest. 
    """
    
    print("searching monallelic patients now!")

    signal_table = int_exome_mt.entries()
    print("turned matrix table to table")
    
    # filter entries table to only include entries with n_alt > 0
    signal_set = signal_table.filter(signal_table.n_alt >= 1)
    
    # export hail table to pandas (~5 minutes) -- alternative is to make a Hail checkpoint...
    signal_set_pd = signal_set.to_pandas()
    
    # group by 's' and sum 'n_alt' (aggregation step)
    geno_agg = signal_set_pd.groupby('s')['n_alt'].sum().reset_index()

    # merge the aggregated sum back to the original dataframe
    geno_agg_pd = pd.merge(signal_set_pd, geno_agg, on='s', how='left', suffixes=('', '_sum'))

    # rename columns
    geno_agg_pd.rename(columns={'n_alt_sum': 'n_alt_agg', 'annotations.vid':'vid', 
                           'annotations.gene_symbol':'gene', 'annotations.gnomad_all_af':'gnomad_all_af',
                           'annotations.aa_change':'aa_change', 'annotations.dna_change_in_transcript':'dna_change_in_transcript',
                           'annotations.genomic_location':'genomic_location', 'annotations.dbsnp_rsid':'rsid',
                           'annotations.is_canonical_transcript':'is_canonical_transcript',
                           'annotations.revel':'revel',
                            'annotations.splice_ai_acceptor_gain_score':'splice_ai_acceptor_gain_score',
                           'annotations.splice_ai_acceptor_loss_score':'splice_ai_acceptor_loss_score',
                            'annotations.splice_ai_donor_gain_score':'splice_ai_donor_gain_score',
                            'annotations.splice_ai_donor_loss_score':'splice_ai_donor_loss_score'}, inplace=True)
    
    # reduce number of columns (can modify this)
    geno_agg_short_pd = geno_agg_pd[['locus', 'alleles', 'vid', 'gene', 'gnomad_all_af', 'aa_change',
                                     'dna_change_in_transcript', 'rsid', 'revel', 'splice_ai_acceptor_gain_score',
                                     'splice_ai_acceptor_loss_score', 'splice_ai_donor_gain_score',
                                     'splice_ai_donor_loss_score', 'pLOF', 'clinvar', 'lab', 'gnomad_plof', 'gnomad_other',
                                     'gnomad_missense', 's', 'n_alt_agg']]

    geno_agg_short_pd.head()

    return geno_agg_short_pd


# In[ ]:


def variant_counts(geno_agg_short_pd):
    pLOF_pd = geno_agg_short_pd[geno_agg_short_pd['pLOF'] == 1]
    pLOF_count = pLOF_pd['vid'].nunique()
    print(f"The pLOF stream contributes: {str(pLOF_count)} variants")

    lab_pd = geno_agg_short_pd[geno_agg_short_pd['lab'] == 1]
    lab_count = lab_pd['vid'].nunique()
    print(f"The lab stream contributes: {str(lab_count)} variants")

    clinvar_pd = geno_agg_short_pd[geno_agg_short_pd['clinvar'] == 1]
    clinvar_count = clinvar_pd['vid'].nunique()
    print(f"The ClinVar stream contributes: {str(clinvar_count)} variants")
 
    gnomad_plof_pd = geno_agg_short_pd[geno_agg_short_pd['gnomad_plof'] == 1]
    gnomad_plof_count = gnomad_plof_pd['vid'].nunique()
    print(f"The gnomAD pLOF stream contributes: {str(gnomad_plof_count)} variants")
 
    gnomad_other_pd = geno_agg_short_pd[geno_agg_short_pd['gnomad_other'] == 1]
    gnomad_other_count = gnomad_other_pd['vid'].nunique()
    print(f"The gnomAD other stream contributes: {str(gnomad_other_count)} variants")
 
    gnomad_missense_pd = geno_agg_short_pd[geno_agg_short_pd['gnomad_missense'] == 1]
    gnomad_missense_count = gnomad_missense_pd['vid'].nunique()
    print(f"The gnomAD missense stream contributes: {str(gnomad_missense_count)} variants")
    
    geno_agg_unique_count = geno_agg_short_pd['vid'].nunique()
    print(f"The total unique variant count is: {str(geno_agg_unique_count)} variants")


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
 
    def count_participants_with_data(self,
                                     include_icd_codes: bool = True,
                                     include_snomed_codes: bool = True,
                                     include_loinc: bool = True,
                                     include_drugs: bool = True,
                                     measurement_registration_exclusion: bool = True,
                                     custom_conditions: Dict[str, str] = None):
        """
        Count unique persons with data in specified domains.

        Args:
            include_icd_codes: Include ICD9/ICD10 codes from condition_occurrence and observation source fields
            include_snomed_codes: Include SNOMED codes from condition_occurrence and observation source fields
            include_loinc: Include LOINC codes from measurement
            include_drugs: Include drug exposures with domain_id = "Drug"
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

        # Store query parameters for summary
        self.last_query_params = {
            "type": "person_count",
            "include_icd_codes": include_icd_codes,
            "include_snomed_codes": include_snomed_codes,
            "include_loinc": include_loinc,
            "include_drugs": include_drugs,
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
        )

        SELECT COUNT(DISTINCT person_id) AS person_count
        FROM combined
        """)

        # Render the template
        self.query_text = template.render(
            version=self.version,
            condition_list=condition_list,
            include_icd_codes=include_icd_codes,
            include_snomed_codes=include_snomed_codes,
            include_loinc=include_loinc,
            include_drugs=include_drugs,
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

            if data_types:
                summary.append(f"Included Data: {', '.join(data_types)}")

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


# # Gather Homozygote/Heterozygote Counts (Pathfinder)

# ## Read in Files for Analysis

# In[ ]:


# import subprocess

# # copy csv file to the bucket
# args = ["gsutil", "cp", f"./GPR15_variants_v0.csv", f"{bucket}/gpr15/variants/"]
# output = subprocess.run(args, capture_output=True)

# # print output from gsutil
# output.stderr


# In[ ]:


# # List objects in the bucket
# print(subprocess.check_output(f"gsutil -u $GOOGLE_PROJECT ls -r gs://fc-aou-datasets-controlled/v8", shell=True).decode('utf-8'))


# In[ ]:


#Setup/initialize the program
initialize()

#exome matrix table: contains all participants, all variants in exonic regions (+ 15bp into the introns)
exome_mt = hl.read_matrix_table("gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/splitMT/hail.mt")

#variant annotation table
vat_8_ht = hl.read_table('{bucket or my_bucket}/wgs_v8/vat_v8_all.ht')

# df of variants with vid
uploaded_var_pd = pd.read_csv(f'{bucket}/gpr15/variants/GPR15_variants_v0.csv', sep = ',')

# add hail_id
uploaded_var_pd['hail_id'] = 'chr' + uploaded_var_pd['vid'].str.replace('-', ':')


# ## Define hyperparameters

# In[ ]:


#define gene interval for GPR15
gene_interval = ['chr3:98531978-98534681']

#define cut-off for minor allele frequency (MAF) in gnomad (all populations) 
MAF_freq = 0.001


# ## Run

# In[ ]:


filtered_exome = filter_exome(exome_mt, gene_interval)

# mt with processed variants
processed_variants = process_variants(uploaded_var_pd)

# mt with filtered variants and all participants
annotated_exome = annotate_exome(vat_8_ht, gene_interval, processed_variants, filtered_exome)

# pLOF_mt is a matrix table containing only variants from the pLOF stream
pLOF_mt = pLOF_variants(annotated_exome)

# matrix tables containing specific variants
clinvar_mt = clinvar_variants(annotated_exome)
lab_mt = lab_variants(annotated_exome)
gnomad_plof_variants_mt = gnomad_plof_variants(annotated_exome)
gnomad_other_variants_mt = gnomad_other_variants(annotated_exome)
gnomad_missense_variants_mt = gnomad_missense_variants(annotated_exome)

# components you want to include in the filtering process
components = [clinvar_mt, lab_mt, gnomad_plof_variants_mt, gnomad_other_variants_mt, gnomad_missense_variants_mt, pLOF_mt]

# matrix table containing only variants of interest (as specified in "components")
interested_exome = combine_components(components)


# In[ ]:


annotated_exome.filter_rows(hl.is_defined(annotated_exome.annotations.uploaded_var)).rows().show()


# ## Generate participant lists

# In[ ]:


# this is the slow step, takes about 5 min for PKD1
# filters to only include participants with 1 or more variant of interest
# converts hail mt to ht; then, converts ht to pandas df

entries_pd = find_monallelic_participants(interested_exome)


# ## Variant Counts

# In[ ]:


# Step 2: Find people who appear exactly twice in the original data
duplicate_people = entries_pd['s'].value_counts()
people_appearing_twice = duplicate_people[duplicate_people == 2].index


# In[ ]:


entries_pd[entries_pd['s'].isin(people_appearing_twice)].sort_values('s')


# In[ ]:


entries_pd.to_csv(f'{bucket}/data/gpr15/cohorts/pathfinder_results.csv')


# In[ ]:


# Step 1: Create summary by vid
summary_df = entries_pd.groupby('vid').agg({
    'locus': 'first',
    'aa_change': 'first', 
    'dna_change_in_transcript': 'first',
    'rsid': 'first',
    'pLOF': 'first',
    'lab': 'first',
    's': lambda x: (entries_pd.loc[x.index, 'n_alt_agg'] == 1).sum(),  # het count
}).rename(columns={'s': 'het_count'})

# Add homozygote counts
summary_df['hom_count'] = entries_pd.groupby('vid').apply(
    lambda x: (x['n_alt_agg'] == 2).sum()
).values

summary_df = summary_df.reset_index()


# In[ ]:


summary_df.to_csv(f'{bucket}/data/gpr15/cohorts/pathfinder_summary.csv')


# In[ ]:


summary_df


# In[ ]:


from upsetplot import plot
from upsetplot import from_memberships


# In[ ]:


# Create boolean columns for each annotation 
annotation_cols = ['pLOF', 'lab', 'clinvar', 'gnomad_plof', 'gnomad_other', 'gnomad_missense']

# Plot 1: All variants
data_all = entries_pd[annotation_cols].astype(bool).reset_index(drop=True)
result_list = [
    list(data_all.columns[row].to_list()) 
    for row in data_all.values.astype(bool)
]

# Create the UpSet plot data
data = from_memberships(result_list)

# Create the UpSet plot
plot(data, subset_size="count", sort_by='cardinality', sort_categories_by='-cardinality', show_counts=True)
plt.title(f'UpSet Plot of Variant Annotations: All Participant Variants ({len(entries_pd)})')
plt.show()

# Plot 2: Unique variants
# Create composite key for unique variants
# Get unique variants
entries_pd_unique = entries_pd.drop_duplicates(subset=['dna_change_in_transcript'])
data_all_unique = entries_pd_unique[annotation_cols].astype(bool).reset_index(drop=True)
result_list_unique = [
    list(data_all_unique.columns[row].to_list()) 
    for row in data_all_unique.values.astype(bool)
]

# Create the UpSet plot data
data_unique = from_memberships(result_list_unique)

# Create the UpSet plot
plot(data_unique, subset_size="count", sort_by='cardinality', sort_categories_by='-cardinality', show_counts=True)
plt.title(f'UpSet Plot of Variant Annotations: Unique Variants ({len(entries_pd_unique)})')
plt.show()


# In[ ]:


tyr132ser_pathfinder_df = entries_pd[entries_pd['vid']=='3-98532428-A-C']
tyr215ter_pathfinder_df = entries_pd[entries_pd['vid']=='3-98532678-C-G']
asp306asn_pathfinder_df = entries_pd[entries_pd['vid']=='3-98532949-G-A']


# # Generate Cohort (PheTK) - check related

# In[ ]:


# !pip install PheTK


# In[ ]:


from PheTK.Cohort import Cohort
from PheTK.Phecode import Phecode
from PheTK.PheWAS import PheWAS
from PheTK.Plot import Plot
import PheTK._queries as _queries
from PheTK import _utils


# In[ ]:


# instantiate class Cohort object for _All of Us_
cohort = Cohort(platform="aou", aou_db_version=8)


# ## Generate Cohorts

# In[ ]:


# Select one variant above to create the entire AoU cohort, with case == 1 set for that variant
cohort.by_genotype(
    chromosome_number=3,
    genomic_position=98532428,
    ref_allele="A",
    alt_allele="C",
    case_gt="1/1",
    control_gt="0/0",
    reference_genome="GRCh38",
    mt_path="gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/multiMT/hail.mt"
)


# In[ ]:


# Select one variant above to create the entire AoU cohort, with case == 1 set for that variant
cohort.by_genotype(
    chromosome_number=3,
    genomic_position=98532678,
    ref_allele="C",
    alt_allele="G",
    case_gt="0/1",
    control_gt="0/0",
    reference_genome="GRCh38",
    mt_path="gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/multiMT/hail.mt"
)


# In[ ]:


# Select one variant above to create the entire AoU cohort, with case == 1 set for that variant
cohort.by_genotype(
    chromosome_number=3,
    genomic_position=98532949,
    ref_allele="G",
    alt_allele="A",
    case_gt="0/1",
    control_gt="0/0",
    reference_genome="GRCh38",
    mt_path="gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/multiMT/hail.mt"
)


# In[ ]:


# # copy csv file to the bucket
# args = ["gsutil", "cp", f"./aou_chr3_98532428_A_C.csv", f"{bucket}/data/gpr15/cohorts/"]
# output = subprocess.run(args, capture_output=True)

# # print output from gsutil
# output.stderr


# In[ ]:


# # copy csv file to the bucket
# args = ["gsutil", "cp", f"./aou_chr3_98532678_C_G.csv", f"{bucket}/data/gpr15/cohorts/"]
# output = subprocess.run(args, capture_output=True)

# # print output from gsutil
# output.stderr


# In[ ]:


# # copy csv file to the bucket
# args = ["gsutil", "cp", f"./aou_chr3_98532949_G_A.csv", f"{bucket}/data/gpr15/cohorts/"]
# output = subprocess.run(args, capture_output=True)

# # print output from gsutil
# output.stderr


# In[ ]:


tyr132ser_df = pd.read_csv(f'{bucket}/data/gpr15/cohorts/aou_chr3_98532428_A_C.csv')
tyr215ter_df = pd.read_csv(f'{bucket}/data/gpr15/cohorts/aou_chr3_98532678_C_G.csv')
asp306asn_df = pd.read_csv(f'{bucket}/data/gpr15/cohorts/aou_chr3_98532949_G_A.csv')


# ## Compare Pathfinder and PheTK

# In[ ]:


tyr132ser_pathfinder_df['s'] = tyr132ser_pathfinder_df['s'].astype(str)
tyr132ser_df['person_id'] = tyr132ser_df['person_id'].astype(str)

tyr132ser_pathfinder_df = tyr132ser_pathfinder_df.merge(
    tyr132ser_df,
    left_on='s',
    right_on='person_id',
    how='left'
)


# In[ ]:


# tyr132ser_pathfinder_df[tyr132ser_pathfinder_df['n_alt_agg']==2]


# In[ ]:


# Direct filtering and selection
genotype_info = annotated_exome.filter_entries(
    (annotated_exome.s == "1591765") & 
    (annotated_exome.locus == hl.parse_locus("chr3:98532428", reference_genome='GRCh38'))
).select_entries('GT', 'n_alt').entries().collect()

print("Genotype info:")
for entry in genotype_info:
    print(f"Locus: {entry.locus}")
    print(f"Alleles: {entry.alleles}")
    print(f"Genotype: {entry.GT}")
    print(f"n_alt: {entry.n_alt}")


# In[ ]:


# Direct filtering and selection
genotype_info = annotated_exome.filter_entries(
    (annotated_exome.s == "1639589") & 
    (annotated_exome.locus == hl.parse_locus("chr3:98532428", reference_genome='GRCh38'))
).select_entries('GT', 'n_alt').entries().collect()

print("Genotype info:")
for entry in genotype_info:
    print(f"Locus: {entry.locus}")
    print(f"Alleles: {entry.alleles}")
    print(f"Genotype: {entry.GT}")
    print(f"n_alt: {entry.n_alt}")


# In[ ]:


tyr215ter_pathfinder_df['s'] = tyr215ter_pathfinder_df['s'].astype(str)
tyr215ter_df['person_id'] = tyr215ter_df['person_id'].astype(str)

tyr215ter_pathfinder_df = tyr215ter_pathfinder_df.merge(
    tyr215ter_df,
    left_on='s',
    right_on='person_id',
    how='left'
)


# In[ ]:


# tyr215ter_pathfinder_df[tyr215ter_pathfinder_df['n_alt_agg']==1]


# In[ ]:


asp306asn_pathfinder_df['s'] = asp306asn_pathfinder_df['s'].astype(str)
asp306asn_df['person_id'] = asp306asn_df['person_id'].astype(str)

asp306asn_df = asp306asn_df.merge(
    asp306asn_pathfinder_df,
    right_on='s',
    left_on='person_id',
    how='left'
)


# In[ ]:


asp306asn_df[(asp306asn_df['case']==1) & (asp306asn_df['n_alt_agg']!=1)]


# In[ ]:


# Direct filtering and selection
genotype_info = annotated_exome.filter_entries(
    (annotated_exome.s == "1639589") & 
    (annotated_exome.locus == hl.parse_locus("chr3:98532949", reference_genome='GRCh38'))
).select_entries('GT', 'n_alt').entries().collect()

print("Genotype info:")
for entry in genotype_info:
    print(f"Locus: {entry.locus}")
    print(f"Alleles: {entry.alleles}")
    print(f"Genotype: {entry.GT}")
    print(f"n_alt: {entry.n_alt}")


# ## Add covariate data
# Covariate descriptions:<br>
# 
# natural_age: current age or age at death<br>
# age_at_last_event: age at last diagnosis event (ICD or SNOMED) in EHR.<br>
# sex_at_birth: sex at birth<br>
# ehr_length: EHR duration, in year, from first to last diagnosis code<br>
# dx_code_occurrence_count: counts the occurrences of diagnosis codes throughout EHR of each participant. For example: person 1 having R50 (fever) code on 5 different dates, R05 (cough) code on 3 different dates, and R05.1 (acute cough) code on 2 different dates, will have a dx_code_occurrence_count = 10.<br>
# dx_condition_count: counts the number of unique conditions occurred throughout EHR of each participant. For example, for the same person 1 above, the dx_condition_count = 3 (R05 - cough, R05.1 - acute cough, R50 - fever).<br>
# genetic_ancestry: returns string values of predicted ancestries, e.g., "eur", "afr", etc. These are only useful if user would like to filter data by genetic ancestries.<br>
# first_n_pcs: retrieves first n genetic PC components from genetic PCA data generated by All of Us.<br>
# drop_nulls: remove rows containing null values in any column.

# In[ ]:


get_ipython().run_cell_magic('time', '', 'cohort.add_covariates(\n    cohort_csv_path=f"{bucket}/data/gpr15/cohorts/aou_chr3_98532428_A_C.csv",\n    natural_age=False,\n    age_at_last_event=True,\n    sex_at_birth=True,\n    ehr_length=True,\n    dx_code_occurrence_count=False,\n    dx_condition_count=False,\n    genetic_ancestry=True,\n    first_n_pcs=10,\n    drop_nulls=True,\n    output_file_name=f"{bucket}/data/gpr15/cohorts/aou_chr3_98532428_A_C_with_covariates.csv"\n)\n')


# In[ ]:


get_ipython().run_cell_magic('time', '', 'cohort.add_covariates(\n    cohort_csv_path=f"{bucket}/data/gpr15/cohorts/aou_chr3_98532678_C_G.csv",\n    natural_age=False,\n    age_at_last_event=True,\n    sex_at_birth=True,\n    ehr_length=True,\n    dx_code_occurrence_count=False,\n    dx_condition_count=False,\n    genetic_ancestry=True,\n    first_n_pcs=10,\n    drop_nulls=True,\n    output_file_name=f"{bucket}/data/gpr15/cohorts/aou_chr3_98532678_C_G_with_covariates.csv"\n)\n')


# In[ ]:


get_ipython().run_cell_magic('time', '', 'cohort.add_covariates(\n    cohort_csv_path=f"{bucket}/data/gpr15/cohorts/aou_chr3_98532949_G_A.csv",\n    natural_age=False,\n    age_at_last_event=True,\n    sex_at_birth=True,\n    ehr_length=True,\n    dx_code_occurrence_count=False,\n    dx_condition_count=False,\n    genetic_ancestry=True,\n    first_n_pcs=10,\n    drop_nulls=True,\n    output_file_name=f"{bucket}/data/gpr15/cohorts/aou_chr3_98532949_G_A_with_covariates.csv"\n)\n')


# ## Add IBD data

# ### Initialize

# In[ ]:


phecodex_map = pl.read_csv('phecodeX_unrolled_ICD_CM.csv')


# In[ ]:


# Initialize
aou = AouQueries()


# ### IBD by AoUQueries (ICD and SNOMED)

# In[ ]:


# Find IBD diagnoses by string
gi_inflammation_codes_string_df = (aou.find_diagnosis_codes(
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["ulcerative colitis", "crohn", 'ulcerative (chronic) pancolitis',
                  'ulcerative (chronic) rectosigmoiditis', 'ulcerative (chronic) proctitis', 
                  'left sided colitis', 'regional enteritis', 'ulcerative (chronic) enterocolitis',
                  'ulcerative (chronic) colitis', 'ulcerative (chronic) proctosigmoiditis', 
                  'ulcerative (chronic) ileocolitis', 'pseudopolyposis', 'inflammatory polyps of colon',
                  'lymphocytic colitis', 'collagenous colitis', 'microscopic colitis', 'enteritis and colitis',
                  'protein-induced enteropathy', 'protein-induced enterocolitis', 'eosinophilic esophagitis',
                  'eosinophilic gastr', 'eosinophilic colitis', 'ulceration of intestine', 'ulcer of anus', 
                  'ulcer of intestine', 'ulcer of anus', 'duodenitis', 'gastric mucosal hypertrophy', 
                  'specified gastritis', 'alcoholic gastritis', 'gastritis, unspecified', 'other gastritis',
                  'chronic gastritis', 'superficial gastritis', 'acute gastritis', 'atrophic gastritis', 
                  'noninfective gastroenteritis', 'indeterminate colitis'
                 ],
    exclude_terms=['fh:', 'h/o:', 'radiation', 'infectious gastroenteritis and colitis, unspecified', 
                  'helicobacter', 'history of']
).execute_gbq())


# In[ ]:


gi_inflammation_dict = create_icd_dict_from_phecodes(phecodex_map, ['GI_522'])


# In[ ]:


# Find IBD diagnoses by phecodes
gi_inflammation_codes_phecode_df = (aou.find_diagnosis_codes(
    exact_codes=gi_inflammation_dict
).execute_gbq())


# In[ ]:


# Create a combined dataframe with source labels
phecode_df = gi_inflammation_codes_phecode_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("phecode").alias("source"))
string_df = gi_inflammation_codes_string_df.unique(subset=["vocabulary_id", "concept_code"]).with_columns(pl.lit("string").alias("source"))

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


gi_inflammation_code_df = (aou.person_code_df(
    name="gi_inflammation",
    vocabularies=["ICD9CM", "ICD10CM", "SNOMED"],
    search_terms=["ulcerative colitis", "crohn", 'ulcerative (chronic) pancolitis',
                  'ulcerative (chronic) rectosigmoiditis', 'ulcerative (chronic) proctitis', 
                  'left sided colitis', 'regional enteritis', 'ulcerative (chronic) enterocolitis',
                  'ulcerative (chronic) colitis', 'ulcerative (chronic) proctosigmoiditis', 
                  'ulcerative (chronic) ileocolitis', 'pseudopolyposis', 'inflammatory polyps of colon',
                  'lymphocytic colitis', 'collagenous colitis', 'microscopic colitis', 'enteritis and colitis',
                  'protein-induced enteropathy', 'protein-induced enterocolitis', 'eosinophilic esophagitis',
                  'eosinophilic gastr', 'eosinophilic colitis', 'ulceration of intestine', 'ulcer of anus', 
                  'ulcer of intestine', 'ulcer of anus', 'duodenitis', 'gastric mucosal hypertrophy', 
                  'specified gastritis', 'alcoholic gastritis', 'gastritis, unspecified', 'other gastritis',
                  'chronic gastritis', 'superficial gastritis', 'acute gastritis', 'atrophic gastritis', 
                  'noninfective gastroenteritis', 'indeterminate colitis'
                 ],
    exclude_terms=['fh:', 'h/o:', 'radiation', 'infectious gastroenteritis and colitis, unspecified', 
                  'helicobacter', 'history of'],
    dates=True
).execute_gbq())


# In[ ]:


gi_inflammation_code_df.head()


# In[ ]:


print(f'{gi_inflammation_code_df.height} person_id with at least 1 GI inflammation code')
print(f'{gi_inflammation_code_df.filter(pl.col("gi_inflammation_n")>1).height} person_id with at least 2 GI inflammation codes')

# Keep those with at least 1 code
result_df = gi_inflammation_code_df.filter(pl.col('gi_inflammation_n')>0).select(['person_id', 'gi_inflammation_1', 'gi_inflammation_2'])


# ### IBD by PheTK Method (Phecode Only)

# In[ ]:


# phecode = Phecode(platform="aou")


# In[ ]:


# phecode.icd_events.head()


# In[ ]:


# phecode_df = _utils.get_phecode_mapping_table(
#     phecode_version="X",
#     icd_version="US",
#     phecode_map_file_path=None,
#     keep_all_columns=False
# )


# In[ ]:


# phecode.icd_events = phecode.icd_events.join(
#     phecode_df,
#     how='inner',
#     on=['ICD', 'flag']
# )


# In[ ]:


# # Filter for relevant phecodes
# gi_inflammation_codes = phecode.icd_events.filter(
#    pl.col('phecode').is_in(['GI_522.11', 'GI_522.12'])
# )

# # Sort by person_id and date to get chronological order
# gi_inflammation_codes = gi_inflammation_codes.sort(['person_id', 'date'])

# # Get first two IBD events (regardless of CD/UC)
# gi_inflammation_events = (
#     gi_inflammation_codes
#     .group_by('person_id')
#     .agg([
#         pl.col('date').head(2).alias('dates'),
#         pl.col('phecode').head(2).alias('phecodes')
#     ])
#     .with_columns([
#         pl.col('dates').list.get(0, null_on_oob=True).alias('gi_inflammation_1'),
#         pl.col('dates').list.get(1, null_on_oob=True).alias('gi_inflammation_2'),
#         pl.col('phecodes').list.get(0, null_on_oob=True).alias('gi_inflammation_1_code'),
#         pl.col('phecodes').list.get(1, null_on_oob=True).alias('gi_inflammation_2_code')
#     ])
#     .drop(['dates', 'phecodes'])
# )


# In[ ]:


# # Count occurrences for percentage
# code_counts = (
#    gi_inflammation_codes
#    .group_by('person_id')
#    .agg([
#        pl.len().alias('total_events'),
#        pl.col('phecode').filter(pl.col('phecode') == 'GI_522.11').len().alias('cd_count'),
#        pl.col('phecode').filter(pl.col('phecode') == 'GI_522.12').len().alias('uc_count')
#    ])
#    .with_columns([
#        pl.when(pl.col('total_events') >= 1)
#        .then((pl.col('cd_count') / pl.col('total_events') * 100).round(1))
#        .otherwise(None)
#        .alias('crohns_pct')
#    ])
#    .with_columns([
#        pl.when(pl.col('total_events') >= 1)
#        .then((pl.col('uc_count') / pl.col('total_events') * 100).round(1))
#        .otherwise(None)
#        .alias('uc_pct')
#    ])
# )

# # Join the results
# result_df = gi_inflammation_events.join(code_counts, on='person_id', how='left')

# print(f'{result_df.height} person_id with at least 1 IBD code')
# print(f'{result_df.filter(pl.col("gi_inflammation_2").is_not_null()).height} person_id with at least 2 IBD codes')

# # Create breakdown column
# result_df = result_df.with_columns([
#     pl.when(pl.col('crohns_pct').is_not_null())
#     .then(pl.format("{}% Crohn's, {}% UC", pl.col('crohns_pct'), pl.col('uc_pct')))
#     .otherwise(pl.lit("Insufficient events"))
#     .alias('breakdown')
# ])

# result_df = result_df.filter(pl.col('gi_inflammation_1').is_not_null()).select(['person_id', 'gi_inflammation_1', 'total_events', 
#                                                                     'cd_count', 'uc_count', 'breakdown'])


# ### Add Age at IBD

# In[ ]:


version = get_ipython().run_line_magic('env', 'WORKSPACE_CDR')

dob_q = f"""
SELECT DISTINCT
    p.person_id,
    CAST(p.birth_datetime AS DATE) AS dob,
FROM
    {version}.person p
"""

dob_df = polars_gbq(dob_q)


# In[ ]:


result_df = result_df.join(
    dob_df,
    on='person_id',
    how='left'
).with_columns([
   pl.when(pl.col('gi_inflammation_1').is_not_null())
   .then((pl.col('gi_inflammation_1') - pl.col('dob')).dt.total_days() / 365.25)
   .otherwise(None)
   .alias('age_at_gi_inflammation_1'),
   pl.when(pl.col('gi_inflammation_2').is_not_null())
   .then((pl.col('gi_inflammation_2') - pl.col('dob')).dt.total_days() / 365.25)
   .otherwise(None)
   .alias('age_at_gi_inflammation_2')

])


# In[ ]:


result_df.head()


# In[ ]:


### Write Files
# files = {
#    'Tyr132Ser 1/1': f"{bucket}/data/gpr15/cohorts/aou_chr3_98532428_A_C_with_covariates",
#    'Tyr215Ter 0/1': f"{bucket}/data/gpr15/cohorts/aou_chr3_98532678_C_G_with_covariates", 
#    'Asp306Asn 0/1': f"{bucket}/data/gpr15/cohorts/aou_chr3_98532949_G_A_with_covariates"
# }

# for _, file in files.items():
#     df = pl.read_csv(f'{file}.csv', separator=',')
#     df = df.join(
#         result_df.select(['person_id', 'gi_inflammation_1', 'age_at_gi_inflammation_1', 'total_events', 'cd_count',
#                           'uc_count', 'breakdown']),
#         on='person_id',
#         how='left'
#     ).with_columns([
#         pl.col('gi_inflammation_1').fill_null(pl.lit("1901-01-01").str.to_date()),
#         pl.col('age_at_gi_inflammation_1').fill_null(-9),
#         pl.col('total_events').fill_null(-9),
#         pl.col('cd_count').fill_null(-9),
#         pl.col('uc_count').fill_null(-9),
#         pl.col('breakdown').fill_null("NA")
#     ])
#     df.write_csv(f'{file}_and_gi_inflammation_phetk.csv')


# In[ ]:


### Load Files
files = {
   'Tyr132Ser 1/1': f"{bucket}/data/gpr15/cohorts/aou_chr3_98532428_A_C_with_covariates",
   'Tyr215Ter 0/1': f"{bucket}/data/gpr15/cohorts/aou_chr3_98532678_C_G_with_covariates", 
   'Asp306Asn 0/1': f"{bucket}/data/gpr15/cohorts/aou_chr3_98532949_G_A_with_covariates"
}

for _, file in files.items():
    df = pl.read_csv(f'{file}.csv', separator=',')
    df = df.join(
        result_df.select(['person_id', 'gi_inflammation_1', 'age_at_gi_inflammation_1', 
                          'gi_inflammation_2', 'age_at_gi_inflammation_2']),
        on='person_id',
        how='left'
    ).with_columns([
        pl.col('gi_inflammation_1').fill_null(pl.lit("1901-01-01").str.to_date()),
        pl.col('age_at_gi_inflammation_1').fill_null(-9),
        pl.col('gi_inflammation_2').fill_null(pl.lit("1901-01-01").str.to_date()),
        pl.col('age_at_gi_inflammation_2').fill_null(-9),
    ])
    df.write_csv(f'{file}_and_gi_inflammation_aouq.csv')


# # Assess Genetic Ancestry

# In[ ]:


# Hard code consistent color palette before the loop
ancestry_colors = {
   'eur': '#1f77b4',  # blue
   'afr': '#ff7f0e',  # orange
   'amr': '#2ca02c',  # green
   'eas': '#d62728',  # red
   'mid': '#9467bd',  # purple
   'sas': '#8c564b'   # brown
}

for variant_name, file_path in files.items():
    df = pl.read_csv(f'{file_path}.csv', separator=',')

    # 1. Display breakdown for genetic_ancestry where case == 1
    cases_df = df.filter(pl.col('case') == 1)
    ancestry_breakdown = cases_df.group_by('genetic_ancestry').agg(pl.len().alias('count'))
    print(f"\nVariant: {variant_name}")
    print("Cases by genetic ancestry:")
    display(ancestry_breakdown)

    # Convert to pandas for plotting
    df_pd = df.to_pandas()
    cases_pd = df_pd[df_pd['case'] == 1]

    # 2. Create scatterplots
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), dpi=300)
    fig.suptitle(f'{variant_name} - PC Analysis', fontsize=16)

    # PC1 vs PC2
    sns.scatterplot(data=df_pd, x='pc1', y='pc2', hue='genetic_ancestry', 
                  alpha=0.1, ax=axes[0], legend=False, palette=color_map)
    sns.scatterplot(data=cases_pd, x='pc1', y='pc2', hue='genetic_ancestry', 
                  s=50, ax=axes[0], legend=False, palette=color_map)
    axes[0].set_title('PC1 vs PC2')

    # PC2 vs PC3  
    sns.scatterplot(data=df_pd, x='pc2', y='pc3', hue='genetic_ancestry', 
                  alpha=0.1, ax=axes[1], legend=False, palette=color_map)
    sns.scatterplot(data=cases_pd, x='pc2', y='pc3', hue='genetic_ancestry', 
                  s=50, ax=axes[1], legend=False, palette=color_map)
    axes[1].set_title('PC2 vs PC3')

    # PC1 vs PC3
    sns.scatterplot(data=df_pd, x='pc1', y='pc3', hue='genetic_ancestry', 
                  alpha=0.1, ax=axes[2], legend=False, palette=color_map)
    sns.scatterplot(data=cases_pd, x='pc1', y='pc3', hue='genetic_ancestry', 
                  s=50, ax=axes[2], legend=False, palette=color_map)
    axes[2].set_title('PC1 vs PC3')

    # Add legend to bottom left of third facet
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map[ancestry], 
                        markersize=8, label=ancestry) for ancestry in unique_ancestries]
    axes[2].legend(handles=handles, loc='lower left', title='Genetic Ancestry')

    plt.tight_layout()
    plt.show()

