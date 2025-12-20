#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Load libraries
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


def combine_gwas_results(in_dir, trait_type, out_file):
    df = None
    for chrom in range(1, 23):
        chrom_df = pd.read_csv(
            f'{in_dir}/gwas_results_chr{chrom}.txt',
            sep='\t',
            header=0
        )
        
        if df is None:
            df = chrom_df
        else:
            df = pd.concat([df, chrom_df])
    
    col_dict = {
        'CHR': 'chromosome',
        'POS': 'base_pair_location',
        'MarkerID': 'variant_id',
        'Allele1': 'other_allele',
        'Allele2': 'effect_allele',
        'AF_Allele2': 'effect_allele_frequency',
        'imputationInfo': 'r2',
        'BETA': 'beta',
        'SE': 'standard_error',
        'Tstat': 'test_statistic_value',
        'p.value': 'p_value'
    }
    
    if trait_type == 'quantitative':
        col_dict['N'] = 'n'
    else:
        col_dict['p.value.NA'] = 'p_value__nospa'
        col_dict['Is.SPA'] = 'spa_converged'
        col_dict['AF_case'] = 'effect_allele_frequency__cases'
        col_dict['AF_ctrl'] = 'effect_allele_frequency__control'
        col_dict['N_case'] = 'n_cases'
        col_dict['N_ctrl'] = 'n_controls'
        col_dict['N_case_hom'] = 'n_case__alt_homs'
        col_dict['N_case_het'] = 'n_case__hets'
        col_dict['N_ctrl_hom'] = 'n_controls__alt_homs'
        col_dict['N_ctrl_het'] = 'n_controls__hets'
    
    df = df.rename(columns=col_dict)
    
    # Re-order based on below
    col_order = [
        'chromosome', 'base_pair_location', 'variant_id', 'other_allele', 'effect_allele',
        'effect_allele_frequency', 'r2', 'beta', 'standard_error', 'test_statistic_value',
        'p_value', 'p_value__nospa', 'spa_converged',
        'effect_allele_frequency__cases', 'effect_allele_frequency__control',
        'n', 'n_cases', 'n_controls', 'n_case__alt_homs', 'n_case__hets', 'n_controls__alt_homs', 'n_controls__hets'
    ]
    df = df[[ x for x in col_order if x in df.columns.tolist() ]]
    
    # write and save
    df.to_csv(
        f'{out_file}',
        sep='\t',
        index=False,
        header=True,
        compression='gzip'
    )
    return(df)


# # Merge and plot GWAS results

# In[ ]:


ancestries_considered = ['eur', 'afr', 'amr', 'transancestry'] #, 'eas', 'sas', 'mid']
traits = {
    'condition__hpv': 'binary',
}

init_combine = False


# In[ ]:


if init_combine:
    for anc in ancestries_considered:
        print(f'combining {anc} files')
        for trait in traits.keys():
            rez = combine_gwas_results(
                f'{my_bucket}/saige_gwas/v1_redo/{anc}/{trait}/gwas',
                traits[trait],
                f'{my_bucket}/saige_gwas/v1_redo/{anc}/{trait}/gwas_results.tsv.gz'
            )


# # gwaslab

# In[ ]:


# # !pip install pypandoc==1.5
# !pip install gwaslab


# In[ ]:


import gwaslab as gl


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


# gl.download_ref('1kg_eur_hg38')
# gl.download_ref('1kg_afr_hg38')
# gl.download_ref('1kg_amr_hg38')
# gl.download_ref('1kg_pan_hg38')


# # EUR

# In[ ]:


eur_sumstats = gl.Sumstats(f'{my_bucket}/saige_gwas/v1_redo/eur/condition__hpv/gwas_results.tsv.gz',
                           fmt="saige",
                           snpid='variant_id',
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


# In[ ]:


# Assuming eur_sumstats.data is a pandas Dateurame
# Make a copy to avoid modifying the original data unexpectedly
df = eur_sumstats.data.copy()

df['SNPID'] = df.apply(convert_snp_format, axis=1)

# Update eur_sumstats.data with the modified Dateurame
eur_sumstats.data = df

# Verify the conversion with a sample
print("\nSample of converted data:")
print(df['SNPID'].head())


# In[ ]:


eur_sumstats.get_lead(sig_level=5e-8)


# In[ ]:


eur_sumstats.basic_check()


# In[ ]:


eur_sumstats.get_lead(sig_level=5e-8)


# In[ ]:


eur_sumstats.plot_mqq(skip=2,
                     stratified=True,
                     check=False,
                     fontfamily='DejaVu Sans')


# In[ ]:


eur_lead_snps = eur_sumstats.get_lead(sig_level=1e-5)


# In[ ]:


eur_sig_snps = eur_sumstats.get_lead(sig_level=1e-5, windowsizekb=.001)


# In[ ]:


eur_sig_snps.to_csv(f'{my_bucket}/saige_gwas/v1_redo/eur/condition__hpv/eur_sig_snps.tsv', index=False, sep='\t')


# In[ ]:


# Save regional plot for each lead SNP
ancestry = 'eur'

for index, row in eur_lead_snps.iterrows():
    snp = row['SNPID']
    
    pattern = r'chr(\d+|X|Y):(\d+):'
    match = re.match(pattern, snp)

    if match:
        chrom = int(match.group(1))
        pos = int(match.group(2))

        # Define region (± 500kb)
        start_pos = max(1, pos - 500000)  # Ensure start position is at least 1
        end_pos = pos + 500000
        
        # Create regional plot
        eur_sumstats.plot_mqq(
            mode='r',
            build="99", #99 skips gene annotations
            region=(chrom, start_pos, end_pos),
            vcf_path=gl.get_path("1kg_eur_hg38"),
            # anno_set=[snp],
            check=False,
            fontfamily='DejaVu Sans'
        )

        plt.title(f"Regional plot for {ancestry} - {snp}")
        plt.savefig(f"regional_plots/eur/regional_plot_{ancestry}_{chrom}_{pos}.png")
        plt.close()


# In[ ]:


eur_lead_snps['SNPID']


# In[ ]:


# eur_plot_snps = [
#     '',
# ]


# In[ ]:


# eur_sumstats.plot_mqq(skip=2,
#                       anno_set=eur_plot_snps,
#                       anno=True,
#                       highlight=eur_plot_snps,
#                       stratified=True,
#                       check=False,
#                       fontfamily='DejaVu Sans')


# # AFR

# In[ ]:


afr_sumstats = gl.Sumstats(f'{my_bucket}/saige_gwas/v1_redo/afr/condition__hpv/gwas_results.tsv.gz',
                           fmt="saige",
                           snpid='variant_id',
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


# In[ ]:


# Assuming afr_sumstats.data is a pandas Dataframe
# Make a copy to avoid modifying the original data unexpectedly
df = afr_sumstats.data.copy()

df['SNPID'] = df.apply(convert_snp_format, axis=1)

# Update afr_sumstats.data with the modified Dataframe
afr_sumstats.data = df

# Verify the conversion with a sample
print("\nSample of converted data:")
print(df['SNPID'].head())


# In[ ]:


afr_sumstats.get_lead(sig_level=5e-8)


# In[ ]:


afr_sumstats.basic_check()


# In[ ]:


afr_sumstats.get_lead(sig_level=5e-8)


# In[ ]:


afr_sumstats.plot_mqq(skip=2,
                     stratified=True,
                     check=False,
                     fontfamily='DejaVu Sans')


# In[ ]:


afr_lead_snps = eur_sumstats.get_lead(sig_level=1e-5)


# In[ ]:


afr_sig_snps = afr_sumstats.get_lead(sig_level=1e-5, windowsizekb=.001)


# In[ ]:


afr_sig_snps.to_csv(f'{my_bucket}/saige_gwas/v1_redo/afr/condition__hpv/afr_sig_snps.tsv', index=False, sep='\t')


# In[ ]:


# Save regional plot for each lead SNP
ancestry = 'afr'

for index, row in afr_lead_snps.iterrows():
    snp = row['SNPID']
    
    pattern = r'chr(\d+|X|Y):(\d+):'
    match = re.match(pattern, snp)

    if match:
        chrom = int(match.group(1))
        pos = int(match.group(2))

        # Define region (± 500kb)
        start_pos = max(1, pos - 500000)  # Ensure start position is at least 1
        end_pos = pos + 500000
        
        # Create regional plot
        afr_sumstats.plot_mqq(
            mode='r',
            build="99", #99 skips gene annotations
            region=(chrom, start_pos, end_pos),
            vcf_path=gl.get_path("1kg_afr_hg38"),
            # anno_set=[snp],
            check=False,
            fontfamily='DejaVu Sans'
        )

        plt.title(f"Regional plot for {ancestry} - {snp}")
        plt.savefig(f"regional_plots/afr/regional_plot_{ancestry}_{chrom}_{pos}.png")
        plt.close()


# In[ ]:


afr_lead_snps['SNPID']


# In[ ]:


afr_plot_snps = [
    '',
]


# In[ ]:


# afr_sumstats.plot_mqq(skip=2,
#                       anno_set=afr_plot_snps,
#                       anno=True,
#                       highlight=afr_plot_snps,
#                       stratified=True,
#                       check=False,
#                       fontfamily='DejaVu Sans')


# # AMR

# In[ ]:


amr_sumstats = gl.Sumstats(f'{my_bucket}/saige_gwas/v1_redo/amr/condition__hpv/gwas_results.tsv.gz',
                           fmt="saige",
                           snpid='variant_id',
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


# In[ ]:


# Assuming amr_sumstats.data is a pandas Datamrame
# Make a copy to avoid modifying the original data unexpectedly
df = amr_sumstats.data.copy()

df['SNPID'] = df.apply(convert_snp_format, axis=1)

# Update amr_sumstats.data with the modified Datamrame
amr_sumstats.data = df

# Verify the conversion with a sample
print("\nSample of converted data:")
print(df['SNPID'].head())


# In[ ]:


amr_sumstats.get_lead(sig_level=5e-8)


# In[ ]:


amr_sumstats.basic_check()


# In[ ]:


amr_sumstats.get_lead(sig_level=5e-8)


# In[ ]:


amr_sumstats.plot_mqq(skip=2,
                     stratified=True,
                     check=False,
                     fontfamily='DejaVu Sans')


# In[ ]:


amr_lead_snps = amr_sumstats.get_lead(sig_level=1e-5)


# In[ ]:


amr_sig_snps = amr_sumstats.get_lead(sig_level=1e-5, windowsizekb=.001)


# In[ ]:


amr_sig_snps.to_csv(f'{my_bucket}/saige_gwas/v1_redo/amr/condition__hpv/amr_sig_snps.tsv', index=False, sep='\t')


# In[ ]:


# Save regional plot for each lead SNP
ancestry = 'amr'

for index, row in amr_lead_snps.iterrows():
    snp = row['SNPID']
    
    pattern = r'chr(\d+|X|Y):(\d+):'
    match = re.match(pattern, snp)

    if match:
        chrom = int(match.group(1))
        pos = int(match.group(2))

        # Define region (± 500kb)
        start_pos = max(1, pos - 500000)  # Ensure start position is at least 1
        end_pos = pos + 500000
        
        # Create regional plot
        amr_sumstats.plot_mqq(
            mode='r',
            build="99", #99 skips gene annotations
            region=(chrom, start_pos, end_pos),
            vcf_path=gl.get_path("1kg_amr_hg38"),
            # anno_set=[snp],
            check=False,
            fontfamily='DejaVu Sans'
        )

        plt.title(f"Regional plot for {ancestry} - {snp}")
        plt.savefig(f"regional_plots/amr/regional_plot_{ancestry}_{chrom}_{pos}.png")
        plt.close()


# In[ ]:


amr_lead_snps['SNPID']


# In[ ]:


amr_plot_snps = [
    '',
]


# In[ ]:


# amr_sumstats.plot_mqq(skip=2,
#                       anno_set=amr_plot_snps,
#                       anno=True,
#                       highlight=amr_plot_snps,
#                       stratified=True,
#                       check=False,
#                       fontfamily='DejaVu Sans')


# # Transancestry

# In[ ]:


transancestry_sumstats = gl.Sumstats(f'{my_bucket}/saige_gwas/v1_redo/transancestry/condition__hpv/gwas_results.tsv.gz',
                           fmt="saige",
                           snpid='variant_id',
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


# In[ ]:


# Assuming transancestry_sumstats.data is a pandas DataFrame
# Make a copy to avoid modifying the original data unexpectedly
df = transancestry_sumstats.data.copy()

df['SNPID'] = df.apply(convert_snp_format, axis=1)

# Update transancestry_sumstats.data with the modified DataFrame
transancestry_sumstats.data = df

# Verify the conversion with a sample
print("\nSample of converted data:")
print(df['SNPID'].head())


# In[ ]:


transancestry_sumstats.get_lead(sig_level=5e-8)


# In[ ]:


transancestry_sumstats.basic_check()


# In[ ]:


transancestry_sumstats.get_lead(sig_level=5e-8)


# In[ ]:


transancestry_sumstats.plot_mqq(skip=2,
                     stratified=True,
                     check=False,
                     fontfamily='DejaVu Sans')


# In[ ]:


transancestry_lead_snps = transancestry_sumstats.get_lead(sig_level=1e-5)


# In[ ]:


transancestry_sig_snps = transancestry_sumstats.get_lead(sig_level=1e-5, windowsizekb=.001)


# In[ ]:


# Save regional plot for each lead SNP
ancestry = 'transancestry'

for index, row in transancestry_lead_snps.iterrows():
    snp = row['SNPID']
    
    pattern = r'chr(\d+|X|Y):(\d+):'
    match = re.match(pattern, snp)

    if match:
        chrom = int(match.group(1))
        pos = int(match.group(2))

        # Define region (± 500kb)
        start_pos = max(1, pos - 500000)  # Ensure start position is at least 1
        end_pos = pos + 500000
        
        # Create regional plot
        transancestry_sumstats.plot_mqq(
            mode='r',
            build="99", #99 skips gene annotations
            region=(chrom, start_pos, end_pos),
            vcf_path=gl.get_path("1kg_pan_hg38"),
            # anno_set=[snp],
            check=False,
            fontfamily='DejaVu Sans'
        )

        plt.title(f"Regional plot for {ancestry} - {snp}")
        plt.savefig(f"regional_plots/transancestry/regional_plot_{ancestry}_{chrom}_{pos}.png")
        plt.close()


# In[ ]:


transancestry_lead_snps['SNPID']


# In[ ]:


transancestry_plot_snps = [
    '',
]


# In[ ]:


# transancestry_sumstats.plot_mqq(skip=2,
#                       anno_set=transancestry_plot_snps,
#                       anno=True,
#                       highlight=transancestry_plot_snps,
#                       stratified=True,
#                       check=False,
#                       fontfamily='DejaVu Sans')


# # Grouped Regional Plots

# In[ ]:


snp_dict = {
    "eur": eur_plot_snps,
    "afr": afr_plot_snps,
    "amr": amr_plot_snps,
    "transancestry": transancestry_plot_snps
}


# In[ ]:


# For each SNP in each population
for pop, snp_list in snp_dict.items():
    for snp in snp_list:
        # Extract chromosome and position using regex
        pattern = r'chr(\d+|X|Y):(\d+):'
        match = re.match(pattern, snp)
        if match:
            chrom = int(match.group(1)) if match.group(1).isdigit() else match.group(1)
            pos = int(match.group(2))
            
            # Set LD legends based on which population's SNP we're plotting
            if pop == "eur":
                region_ld_legends = [0]  # Show legend in EUR panel
            elif pop == "afr":
                region_ld_legends = [1]  # Show legend in AFR panel
            elif pop == "amr":
                region_ld_legends = [2]  # Show legend in AMR panel
            elif pop == "transancestry":
                region_ld_legends = [3]  # Show legend in transancestry panel

            # Create stacked regional plot with +/- 500kb
            gl.plot_stacked_mqq(
                objects=[eur_sumstats, afr_sumstats, transancestry_sumstats],
                vcfs=[
                    gl.get_path("1kg_eur_hg38"), 
                    gl.get_path("1kg_afr_hg38"), 
                    gl.get_path("1kg_amr_hg38"), 
                    gl.get_path("1kg_pan_hg38"), 
                ],
                region=(chrom, max(1, pos - 500000), pos + 500000),
                mode='r',
                build="38",
                anno=True,
                region_lead_grids=[],
                region_ld_legends=region_ld_legends,
                anno_source='ensembl',
                anno_style="tight",
                anno_set=[snp],
                check=False,
                ylim=(0,8),
                titles=["EUR", "AFR", "AMR", 'Transancestry'],
                title_args={"size": 20},
                anno_args={"fontsize": 12, "rotation": 0},  # Increased font size
                fontfamily='DejaVu Sans',
            )
            plt.savefig(f"stacked_plots/stacked_plot{chrom}_{pos}.png", dpi=300, bbox_inches='tight')
            plt.close()
        else:
            print(f"Could not parse position from SNP: {snp}")

