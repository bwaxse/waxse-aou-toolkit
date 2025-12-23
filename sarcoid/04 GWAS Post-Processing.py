#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# # this must be installed before gwaslab
# !pip install pypandoc==1.5


# In[ ]:


# # install gwaslab
# !pip install gwaslab


# In[ ]:


# shouldn't need
# !pip install --force-reinstall "numpy<2"


# In[ ]:


# shouldn't need
# !pip install --force-reinstall "matplotlib<3.9" "seaborn<0.13"


# In[ ]:


# Load libraries
from google.cloud import bigquery
from google.cloud import storage
from IPython.display import display, HTML
import pandas as pd
import polars as pl
import scipy
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
from pathlib import Path
import gzip
import re
import io

import os
import subprocess

version = get_ipython().run_line_magic('env', 'WORKSPACE_CDR')
my_bucket = os.getenv('WORKSPACE_BUCKET')


# In[ ]:


# Polars string length to 100
pl.Config.set_fmt_str_lengths(100)

# Set the row limit to a higher value
pl.Config.set_tbl_rows(50)


# # Functions

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


def add_chr_pos_from_snpid(df: pd.DataFrame, snp_col: str = "SNPID") -> pd.DataFrame:
    """Extract CHR and POS from SNPID column (format CHR-POS-...)."""
    parts = df[snp_col].str.split("-", n=2, expand=True)
    df = df.copy()
    df["CHR"] = pd.to_numeric(parts[0], errors="coerce")
    df["POS"] = pd.to_numeric(parts[1], errors="coerce")
    return df


# In[ ]:


# Function to convert SNP format
def convert_snp_format(row):
    snp_id = row['SNPID']
    
    if not isinstance(snp_id, str):
        return snp_id
    
    if snp_id == '.':
        return f"chr{row['CHR']}:{row['POS']}:{row['NEA']}:{row['EA']}"
    
    pattern = r'(\d+|X|Y)_(\d+)_([ACGT])_([ACGT])'
    match = re.match(pattern, snp_id)
    
    if match:
        chrom, pos, ref, alt = match.groups()
        return f"chr{chrom}:{pos}:{ref}:{alt}"
    else:
        return snp_id


# In[ ]:


def parse_loci(loci_list):
    result = []
    for locus in loci_list:
        chr_str, pos_alleles = locus.split(':')
        pos = int(pos_alleles.split('_')[0])
        
        # Handle X/Y chromosomes
        if chr_str == 'X':
            chr_num = 23
        elif chr_str == 'Y':
            chr_num = 24
        else:
            chr_num = int(chr_str)
            
        result.append((chr_num, pos))
    return result


# # Get GWAS Catalog Associations

# In[ ]:


# Define the mapping of conditions to their file(s)
gwas_catalog_files = {
    'condition__sarcoid': ['MONDO_0019338_associations_export.tsv']
}

# Create DataFrames
gwas_catalog_dfs = {}

for condition, files in gwas_catalog_files.items():
    dfs_to_concat = []
    
    for file in files:
        # Use 'condition - file' format
        file_path = f"{my_bucket}/data/saige_gwas/{file}"
        df = pl.read_csv(file_path, separator='\t', columns=['riskAllele', 'efoTraits', 'locations'])
        dfs_to_concat.append(df)
        print(f"Loaded {condition} - {file}: {len(df)} rows")
    
    # Concatenate if multiple files, otherwise use single df
    if len(dfs_to_concat) > 1:
        gwas_catalog_dfs[condition] = pl.concat(dfs_to_concat)
        print(f"Concatenated {condition}: {len(gwas_catalog_dfs[condition])} total rows")
    else:
        gwas_catalog_dfs[condition] = dfs_to_concat[0]
        print(f"Single file for {condition}: {len(gwas_catalog_dfs[condition])} rows")

# Access individual DataFrames
sarcoid_catalog_df = gwas_catalog_dfs['condition__sarcoid']


# In[ ]:


sarcoid_catalog_df.group_by('efoTraits').len()


# In[ ]:


sarcoid_traits = ['sarcoidosis']


# In[ ]:


catalog_data = {
    'condition__sarcoid': (sarcoid_catalog_df, sarcoid_traits),
}

filtered_dfs = {}
for condition, (df, traits) in catalog_data.items():
    filtered_df = df.filter(pl.col('efoTraits').is_in(traits))

    # Remove SNP interactions
    filtered_df = filtered_df.filter(~(pl.col('locations').str.contains(',')))
    
    # Fix missing locations using riskAllele formats:
    # "chr2:211694960-...", "10:94433455-A", "chrX:21534509-ACC"
    # This will still leave null locations for riskAllele == rsID with no location information
    filtered_df = filtered_df.with_columns(
        pl.when(pl.col('locations') == "-")
        .then(
            pl.when(pl.col('riskAllele').str.contains(r"^chr([XY]|\d+):"))
            .then(
                # Handle "chrX:21534509-ACC" or "chr2:211694960-..."
                pl.col('riskAllele').str.extract(r"^chr([XY]|\d+):(\d+)", 1) + ":" + 
                pl.col('riskAllele').str.extract(r"^chr([XY]|\d+):(\d+)", 2)
            )
            .when(pl.col('riskAllele').str.contains(r"^([XY]|\d+):"))
            .then(
                # Handle "10:94433455-A" 
                pl.col('riskAllele').str.extract(r"^([XY]|\d+):(\d+)", 1) + ":" + 
                pl.col('riskAllele').str.extract(r"^([XY]|\d+):(\d+)", 2)
            )
            .otherwise(pl.col('locations'))  # Keep original if no pattern matches
        )
        .otherwise(pl.col('locations'))
        .alias('locations')
    )

    # Add risk nucleotide to locations (only where locations != "-")
    filtered_df = filtered_df.with_columns(
        pl.when((pl.col('locations') != "-") & (~pl.col('riskAllele').str.ends_with("-?")))
        .then(
            pl.col('locations') + "-" + 
            pl.col('riskAllele').str.extract(r"-([ATCG]+)$", 1)
        )
        .otherwise(pl.col('locations'))
        .alias('locations')
    )

    # Sort by chromosome: 1-22, then X, Y
    filtered_df = filtered_df.with_columns(
        pl.col('locations')
        .str.extract(r"^([^:]+):")  # Extract chromosome (including X, Y)
        .map_elements(
            lambda x: int(x) if x and x.isdigit() else (23 if x == 'X' else (24 if x == 'Y' else 99)), 
            return_dtype=pl.Int32
        )
        .alias('chr_sort_order')
    ).sort(['chr_sort_order', 'locations']).drop('chr_sort_order')

    filtered_dfs[condition] = filtered_df
    print(f"{condition}: Filtered from {len(df)} to {len(filtered_df)} rows")


# In[ ]:


# Extract known loci for all conditions
known_loci = {}
for condition in ['condition__sarcoid']:
    valid_loci_df = filtered_dfs[condition].filter(pl.col('locations') != '-')
    total_count = len(valid_loci_df)
    unique_loci = valid_loci_df['locations'].unique().to_list()
    unique_count = len(unique_loci)
    
    known_loci[condition] = unique_loci
    print(f"{condition}: {total_count} loci ≠ '-', {unique_count} unique loci")


# # gwaslab

# In[ ]:


import gwaslab as gl


# In[ ]:


# # Only need to do once
# gl.download_ref('1kg_eur_hg38')
# gl.download_ref('1kg_afr_hg38')
# gl.download_ref('1kg_amr_hg38')
# gl.download_ref('1kg_pan_hg38')


# In[ ]:


ref_populations = {
    'eur': '1kg_eur_hg38',
    'afr': '1kg_afr_hg38', 
    'amr': '1kg_amr_hg38',
    'metal': '1kg_pan_hg38'
}


# In[ ]:


def get_tsv_file(query_dir: str) -> str | None:
    """Return the first .tsv file (not .tsv.info) in a GCS directory."""
    result = subprocess.run(f'gsutil ls {query_dir}', capture_output=True, shell=True)
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip())

    files = result.stdout.decode('utf-8').split('\n')
    tsvs = [f for f in files if f.endswith("1.tsv") and not f.endswith("1.tsv.info")]

    if not tsvs:
        return None  # no .tsv found
    if len(tsvs) > 1:
        print(f"Warning: multiple 1.tsv files in {query_dir}, returning first")

    return tsvs[0]


# In[ ]:


def file_exists_in_gcs(bucket_path: str) -> bool:
    try:
        result = subprocess.run(['gsutil', 'ls', bucket_path], 
                              capture_output=True, text=True, check=False)
        return result.returncode == 0
    except FileNotFoundError:
        raise RuntimeError("gsutil not found - is Google Cloud SDK installed?")


# In[ ]:


def annotate_with_known_loci(lead_snps_df, known_loci_list, window_kb=500):
    """
    Annotate lead SNPs with closest known loci within specified window.
    
    Args:
        lead_snps_df: DataFrame with columns CHR, POS, SNPID
        known_loci_list: List of 'CHR:POS_REF_ALT' strings
        window_kb: Window size in kb (default 500)
    """
    lead_snps_df = lead_snps_df.copy()

    # Parse known loci
    known_loci = []
    excluded_count = 0

    for locus in known_loci_list:
        try:
            # Split on '-' but handle cases where there might be no allele
            if '-' in locus:
                chr_pos, alleles = locus.split('-', 1)  
            else:
                chr_pos = locus
                alleles = None

            chr_str, pos_str = chr_pos.split(':')
            
            # Try to convert to int, skip if it's X/Y/MT
            try:
                chr_num = int(chr_str)
                locus_dict = {
                    'chr': chr_num,
                    'pos': int(pos_str),
                    'original': locus
                }
                # Optionally store allele info if present
                if alleles:
                    locus_dict['alleles'] = alleles

                known_loci.append(locus_dict)
            except ValueError:
                # X, Y, MT chromosomes
                excluded_count += 1
                continue

        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse known locus '{locus}': {e}")
            continue
    
    if excluded_count > 0:
        print(f"Info: Excluded {excluded_count} known loci on non-autosomal chromosomes")
    
    known_df = pd.DataFrame(known_loci)
    window_bp = window_kb * 1000
    
    def find_closest_known(row):
        # Filter to same chromosome
        chr_matches = known_df[known_df['chr'] == row['CHR']]
        if chr_matches.empty:
            return 'None'
        
        # Calculate distances and filter by window
        distances = np.abs(chr_matches['pos'] - row['POS'])
        within_window = distances <= window_bp
        
        if not within_window.any():
            return 'None'
        
        # Return closest match
        closest_idx = distances[within_window].idxmin()
        return chr_matches.loc[closest_idx, 'original']
    
    lead_snps_df['KNOWN'] = lead_snps_df.apply(find_closest_known, axis=1)
    return lead_snps_df


# In[ ]:


def process_gwas(base_dir, ancestry, trait, known_loci_list, plot_params):
    """Process a single GWAS condition: load, parse, analyze, and plot."""
    
    print(f"\nProcessing {trait} {ancestry} GWAS...")
    if ancestry == 'metal':
        # Load summary statistics
        sumstats = gl.Sumstats(
            get_tsv_file(f'{base_dir}/metal/{trait}/'),
            fmt="metal",
            snpid='MarkerName',
            ea='Allele1', 
            nea='Allele2',
            eaf='Freq1',
            se='StdErr',
            p='P-value',
            beta='Effect',
            direction='Direction',
            build='38'
        )
        # Get CHR POS from SNPID
        sumstats.data = add_chr_pos_from_snpid(sumstats.data, snp_col='SNPID')
    else: 
        sumstats = gl.Sumstats(
            f'{base_dir}/{ancestry}/condition__sarcoid/gwas_results.tsv.gz',
            fmt="saige",
            snpid='vid',
            rsid=None,
            chrom='chromosome',
            pos='base_pair_location',
            ea='effect_allele',
            nea='other_allele',
            eaf='effect_allele_frequency',
            se='standard_error',
            p='p_value',
            beta='beta',
            build='38'
        )
        
        # Assuming eur_sumstats.data is a pandas Dateurame
        print("Converting SNPID")
        sumstats.data['SNPID'] = sumstats.data.apply(convert_snp_format, axis=1)
        # Verify the conversion with a sample
        print("Conversion done.")
        
        sumstats.basic_check()
        
    # Get or load lead SNPs (500 kb window)    
    lead_snps = sumstats.get_lead(sig_level=1e-5) #5e-8
    lead_snps = annotate_with_known_loci(lead_snps, known_loci_list)
    lead_snps.to_csv(f'{base_dir}/{ancestry}/{trait}/{trait}_lead_snps.tsv', 
                     sep='\t', index=False)
    
    # Convert known loci matches to highlight format for plotting
    known_matches = lead_snps[lead_snps['KNOWN'] != 'None']['SNPID'].tolist()
    unknown_matches = lead_snps[lead_snps['KNOWN'] == 'None']['SNPID'].tolist()

    # Build base plotting kwargs
    base_plot_kwargs = {
        'skip': 2,
        'check': False,
        'font_family': 'DejaVu Sans',
        'marker_size': (5,5),
        'ylabels': plot_params['ylabels'],
        'fontsize': 11,
        'verbose': False,
        'fig_args': {"figsize":(20,6), "dpi":300},
        'save_args': {"dpi":300, "facecolor":"white"}
    }
    
    # Add cut and cutfactor if they exist and are not None
    if plot_params.get('cut') is not None:
        base_plot_kwargs['cut'] = plot_params['cut']
    if plot_params.get('cutfactor') is not None:
        base_plot_kwargs['cutfactor'] = plot_params['cutfactor']
    
    print("Plotting mqq plots (verbose mode off)")
    
    # Plot highlighted MQQ
    highlighted_kwargs = base_plot_kwargs.copy()
    highlighted_kwargs.update({
        'highlight_windowkb': 500,
        'highlight': [known_matches, unknown_matches],
        'highlight_color': ['#DC143C', '#FFA500'],
        'save': f'manhattan_plots/{trait}_{ancestry}_min_2_highlighted_mqq.png'
    })
    sumstats.plot_mqq(**highlighted_kwargs)
    
    # Plot regular MQQ
    regular_kwargs = base_plot_kwargs.copy()
    regular_kwargs['save'] = f'manhattan_plots/{trait}_{ancestry}_min_2_mqq.png'
    sumstats.plot_mqq(**regular_kwargs)
    
    return sumstats, lead_snps


# In[ ]:


def consolidate_lead_snps_across_ancestries(results, trait, ancestries):
    """
    Consolidate lead SNPs across all ancestries for a trait, keeping the SNP 
    with the lowest p-value for each locus.
    
    Returns:
        dict: {'known': [(snp, known_locus, min_p)], 'novel': [(snp, min_p)]}
    """
    # Collect all lead SNPs across ancestries
    all_leads = []

    for ancestry, data in results[trait].items():
        if 'lead_snps' not in data:
            continue

        lead_snps = data['lead_snps']
        for _, row in lead_snps.iterrows():
            all_leads.append({
                'snpid': row['SNPID'],
                'chr': row['CHR'], 
                'pos': row['POS'],
                'p': row['P'],
                'known': row['KNOWN'],
                'ancestry': ancestry
            })
                
    if not all_leads:
        return {'known': [], 'novel': []}
    
    leads_df = pd.DataFrame(all_leads)
    
    # Group nearby SNPs (within 500kb) and keep best p-value
    consolidated = {'known': [], 'novel': []}
    processed_positions = set()
    
    # Sort by p-value to process best hits first
    leads_df = leads_df.sort_values('p')
    
    for _, row in leads_df.iterrows():
        # Check if we've already processed a nearby SNP
        pos_key = f"{row['chr']}:{row['pos']}"
        if any(abs(row['pos'] - int(pos.split(':')[1])) < 500000 
               for pos in processed_positions 
               if pos.split(':')[0] == str(row['chr'])):
            continue
            
        processed_positions.add(pos_key)
        
        if row['known'] != 'None':
            consolidated['known'].append((row['snpid'], row['known'], row['p']))
        else:
            consolidated['novel'].append((row['snpid'], row['p']))
    
    return consolidated


# In[ ]:


def prepare_regional_plot_snps(consolidated_leads, results, trait, ancestry, known_loci_list):
    """
    Prepare SNPs for regional plotting with proper annotation handling.
    
    Returns:
        list: [(snp, ancestry, anno_alias)]
    """
    plot_snps = []
    
    # Get all lead SNPs from all conditions for reference
    all_trait_snps = set()
    for ancestry, data in results[trait].items():
        if 'lead_snps' in data:
            all_trait_snps.update(data['lead_snps']['SNPID'].tolist())
    
    # Process known loci
    for snp, known_locus, p_val in consolidated_leads['known']:
        # Extract chr:pos from known_locus (format: 'CHR:POS-')
        try:
            known_chr_pos = known_locus.split('-')[0]  # Gets 'CHR:POS'
        except (ValueError, IndexError):
            print(f"Warning: Could not parse known_locus format: {known_locus}")
            continue
        
        # Extract chr:pos from snp (format: 'CHR-POS-REF-ALT')
        try:
            snp_parts = snp.split('-')
            snp_chr_pos = f"{snp_parts[0]}:{snp_parts[1]}"  # Convert to 'CHR:POS'
        except (ValueError, IndexError):
            print(f"Warning: Could not parse SNP format: {snp}")
            continue
        
        if snp_chr_pos == known_chr_pos:
            # SNP is at same chr:pos as known locus
            anno_alias = {snp: f'  {snp}\n  Reported: {known_locus}'}
        else:
            # Find closest SNP in dataset to the known locus
            closest_snp = find_closest_snp_in_dataset(known_locus, all_trait_snps)
            if closest_snp:
                anno_alias = {closest_snp: f'  {closest_snp}\n  Reported: {known_locus}'}
            else:
                print(f"Warning: No close SNP found for known locus {known_locus}")
                continue
                
        plot_snps.append((snp, 'reported', anno_alias, p_val))
    
    # Process novel loci  
    for snp, p_val in consolidated_leads['novel']:
        anno_alias = {snp: f'  Unreported: {snp}'}
        plot_snps.append((snp, 'unreported', anno_alias, p_val))
    
    return plot_snps


# In[ ]:


def find_closest_snp_in_dataset(known_locus, dataset_snps, max_distance_kb=500):
    """
    Find the closest SNP in dataset to a known locus.
    
    Args:
        known_locus: str in format 'CHR:POS-ALT'
        dataset_snps: set of SNPs in format 'CHR-POS-REF-ALT'
        max_distance_kb: maximum distance in kb
    """
    try:
        # Parse known locus
        if '-' in known_locus:
            chr_pos = known_locus.split('-')[0]
        else:
            chr_pos = known_locus
        target_chr, target_pos = chr_pos.split(':')
        target_chr = int(target_chr)
        target_pos = int(target_pos)
    except (ValueError, IndexError):
        return None
    
    closest_snp = None
    min_distance = float('inf')
    max_distance_bp = max_distance_kb * 1000
    
    # Parse dataset SNPs and find closest
    for snp in dataset_snps:
        match = re.match(r'^(\d+|X|Y)-(\d+)-([ATCG])-([ATCG])$', snp)
        if not match:
            continue
            
        try:
            snp_chr = int(match.group(1))
            snp_pos = int(match.group(2))
        except ValueError:
            continue  # Skip X, Y chromosomes
            
        if snp_chr != target_chr:
            continue
            
        distance = abs(snp_pos - target_pos)
        if distance < min_distance and distance <= max_distance_bp:
            min_distance = distance
            closest_snp = snp
    
    return closest_snp


# In[ ]:


def create_stacked_regional_plots(results, trait, ancestries, trait_sumstats, plot_titles, ref_populations):
    """
    Create stacked regional plots for consolidated lead SNPs.
    """
    print(f"\nCreating consolidated regional plots")
    
    # Consolidate lead SNPs across conditions
    consolidated = consolidate_lead_snps_across_ancestries(results, trait, ancestries)
        
    if not consolidated['known'] and not consolidated['novel']:
        print(f"No lead SNPs found for {trait}")
        return
        
    # Get known loci for this trait
    known_loci_list = known_loci[trait]
    
    # Prepare SNPs for plotting
    plot_snps = prepare_regional_plot_snps(
        consolidated, results, trait, ancestry, known_loci_list
    )
    
    if not plot_snps:
        print(f"No valid SNPs for regional plotting in {trait}")
        return
    
    # Create regional plots
    for plot_snp, plot_type, anno_alias, p_val in plot_snps:
        try:
            # Parse coordinates for region
            match = re.match(r'^(\d+|X|Y)-(\d+)-([ATCG])-([ATCG])$', plot_snp)
            if not match:
                print(f"Warning: Could not parse SNP format: {plot_snp}")
                continue

            chrom = int(match.group(1))
            pos = int(match.group(2))

            print(f"Plotting {plot_type} locus: {plot_snp}")

            all_paths = [gl.get_path(ref_populations.get(ancestry)) for ancestry in ancestries]
            
            # Calculate the maximum -log10(p) value
            neg_log10_p = -np.log10(p_val)

            # Round up to next multiple of 5, then add 5
            max_y = (np.ceil(neg_log10_p / 5) * 5) + 5

            # Create labels in groups of 5
            ylabels = list(range(0, int(max_y) + 1, 5))

            gl.plot_stacked_mqq(
                objects=trait_sumstats,
                vcfs=all_paths,
                region=(chrom, pos - 250000, pos + 250000),
                region_ref=[plot_snp],
                build="38",
                mode="r",
                anno=True,
                anno_set=[plot_snp],
                anno_alias=anno_alias,
                anno_style="tight",
                fontfamily='DejaVu Sans',
                marker_size=(50, 50),
                titles=plot_titles,
                track_n=3,
                track_fontsize_ratio=1.2,
                ylabels=ylabels,
                title_args={"size": 20},
                anno_args={"rotation": 0},
                fig_args={"figsize": (14, 14), "dpi": 300},
                save_args={"dpi": 300, "facecolor": "white"},
                save=f'regional_plots/{trait}_{plot_type}_{chrom}_{pos}_stacked_regional.png',
                verbose=True,
                check=False
            )            
        except Exception as e:
            print(f"Error plotting regional for {plot_snp}: {e}")
            continue


# # Run All

# In[ ]:


base_dir = f'{my_bucket}/saige_gwas/min_2' # no trailing '/'


# In[ ]:


traits = ['condition__sarcoid']

# Define disease-specific plotting parameters
trait_plot_params = {
    'condition__sarcoid': {'cut': None, 'cutfactor': None, 'ylabels': [2,4,6,8,10]}
}


# In[ ]:


# Process all conditions for all traits
results = {}
for trait in traits:
    print(f"\n=== Processing trait: {trait} ===")
    results[trait] = {}
    
    # Get trait-specific parameters
    plot_params = trait_plot_params.get(trait, {'ylabels': [2,8,14,20,50,80,100]})

    trait_sumstats = []
    plot_titles = []
    all_trait_known = []
    all_trait_new = []
    
    # Single Ancestries
    for ancestry in ['eur', 'afr']:
        try:
            sumstats, lead_snps = process_gwas(
                base_dir, ancestry, trait, known_loci[trait], plot_params
            )
            results[trait][ancestry] = {
                'sumstats': sumstats,
                'lead_snps': lead_snps
            }
            trait_sumstats.append(sumstats)
            plot_titles.append(f"{ancestry.upper()}")

            # Collect loci for stacked plot highlighting
            ancestry_known_loci = lead_snps[lead_snps['KNOWN'] != 'None']['SNPID'].tolist()
            all_trait_known.extend(ancestry_known_loci)
            new_loci = lead_snps[lead_snps['KNOWN'] == 'None']['SNPID'].tolist()
            all_trait_new.extend(new_loci)

        except Exception as e:
            print(f"ERROR processing {ancestry}: {e}")
            continue

#     # METAL
#     try:
#         sumstats, lead_snps = process_gwas(
#             base_dir, 'metal', trait, known_loci[trait], plot_params
#         )
#         results[trait]['metal'] = {
#             'sumstats': sumstats,
#             'lead_snps': lead_snps
#         }
#         trait_sumstats.append(sumstats)
#         plot_titles.append(f"Meta-analysis")

#         # Collect matches for stacked plot highlighting
#         ancestry_known_loci = lead_snps[lead_snps['KNOWN'] != 'None']['SNPID'].tolist()
#         all_trait_known.extend(ancestry_known_loci)
#         new_loci = lead_snps[lead_snps['KNOWN'] == 'None']['SNPID'].tolist()
#         all_trait_new.extend(new_loci)
        
#     except Exception as e:
#         print(f"ERROR processing Meta-analysis: {e}")
#         continue

    # Create stacked Manhattan plot for this trait (only if we have data)
    if trait_sumstats:
        try:
            unique_all_trait_known = list(set(all_trait_known))
            unique_all_trait_new = list(set(all_trait_new))
            # Build kwargs dictionary for stacked plot
            plot_kwargs = {
                'objects': trait_sumstats,
                'mode': 'm',
                'build': "38", 
                'check': False,
                'skip': 2,
                'ylabels': plot_params['ylabels'],
                'titles': plot_titles,
                'title_args': {"size": 16},
                'fontsize': 11,
                'font_family': 'DejaVu Sans',
                'marker_size': (10,10),
                'verbose': False,
                'fig_args': {"figsize":(16,12), "dpi":300},
                'common_ylabel': True,
                'save_args': {"dpi":300, "facecolor":"white"}
            }
            
            # Only add cut and cutfactor if they're not None
            if plot_params.get('cut') is not None:
                plot_kwargs['cut'] = plot_params['cut']
            if plot_params.get('cutfactor') is not None:
                plot_kwargs['cutfactor'] = plot_params['cutfactor']

            # Plot regular stacked plot
            regular_kwargs = plot_kwargs.copy()
            regular_kwargs['save'] = f'manhattan_plots/{trait}_min_2_stacked_m.png'

            print("Plotting stacked mqq plots (verbose mode off)")
            gl.plot_stacked_mqq(**regular_kwargs)
            
            # Plot highlighted stacked plot
            highlighted_kwargs = plot_kwargs.copy()
            highlighted_kwargs.update({
                'highlight_windowkb': 500,
                'highlight': [unique_all_trait_known, unique_all_trait_new],
                'highlight_color': ['#DC143C', '#FFA500'],
                'save': f'manhattan_plots/{trait}_min_2_highlighted_stacked_m.png'
            })
            gl.plot_stacked_mqq(**highlighted_kwargs)
            plt.close('all')  # Clear all figures
            
            print(f"✓ Stacked plots saved for {trait}")
            
#             # Create consolidated regional plots
#             create_stacked_regional_plots(
#                 results, trait, ['eur', 'afr', 'metal'], trait_sumstats, plot_titles, ref_populations
#             )

        except Exception as e:
            print(f"ERROR creating stacked plot for {trait}: {e}")

print("\n=== Processing complete! ===")


# In[ ]:


get_ipython().system('gsutil ls {base_dir}/manhattan_plots/')


# In[ ]:


# # move all plots from manhattan_plots to bucket metal folder
# args = ["gsutil", "cp", "-r", "./manhattan_plots/", f"{my_bucket}/saige_gwas/min_1/"]
# output = subprocess.run(args, capture_output=True)

# # print output from gsutil
# output.stderr


# # METAL N by SE

# In[ ]:


eur_df = pl.read_csv(f'{base_dir}/eur/condition__sarcoid/gwas_results.tsv.gz',
                     separator='\t')
afr_df = pl.read_csv(f'{base_dir}/afr/condition__sarcoid/gwas_results.tsv.gz',
                     separator='\t')


# In[ ]:


afr_df


# In[ ]:


afr_df.columns


# In[ ]:


eur_df.filter(pl.col('vid')=='3-76477573-T-C')


# In[ ]:


metal_df = pl.read_csv('{bucket or my_bucket}/saige_gwas/min_1/metal/condition__sarcoid/condition__sarcoid_afr_eur1.tsv',
                      separator='\t')


# In[ ]:


import seaborn as sns
import matplotlib.pyplot as plt

# Calculate n_effective for both datasets
eur_df = eur_df.with_columns([
    (pl.col("n_cases") + pl.col("n_controls")).alias("n_effective")
])

afr_df = afr_df.with_columns([
    (pl.col("n_cases") + pl.col("n_controls")).alias("n_effective")
])

# Determine common axis limits
x_min = min(eur_df["standard_error"].min(), afr_df["standard_error"].min())
x_max = max(eur_df["standard_error"].max(), afr_df["standard_error"].max())
y_min = min(eur_df["n_effective"].min(), afr_df["n_effective"].min())
y_max = max(eur_df["n_effective"].max(), afr_df["n_effective"].max())

# Create subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# EUR plot
sns.scatterplot(
    data=eur_df.to_pandas(), 
    x="standard_error", 
    y="n_effective", 
    alpha=0.6, 
    s=20,
    ax=ax1
)
ax1.set_xlim(x_min, x_max)
ax1.set_ylim(y_min, y_max)
ax1.set_title("EUR")
ax1.set_xlabel("Standard Error")
ax1.set_ylabel("Effective Sample Size")

# AFR plot
sns.scatterplot(
    data=afr_df.to_pandas(), 
    x="standard_error", 
    y="n_effective", 
    alpha=0.6, 
    s=20,
    ax=ax2
)
ax2.set_xlim(x_min, x_max)
ax2.set_ylim(y_min, y_max)
ax2.set_title("AFR")
ax2.set_xlabel("Standard Error")
ax2.set_ylabel("Effective Sample Size")

plt.tight_layout()
plt.show()

