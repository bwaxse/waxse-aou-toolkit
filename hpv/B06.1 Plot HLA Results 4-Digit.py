#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Load libraries
from IPython.display import display, HTML
import pandas as pd
import polars as pl
import scipy
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from datetime import date
from pathlib import Path
import gzip
import re

import os
import subprocess

import warnings
warnings.filterwarnings('ignore')

version = get_ipython().run_line_magic('env', 'WORKSPACE_CDR')
my_bucket = os.getenv('WORKSPACE_BUCKET')


# In[ ]:


# Polars string length to 100
pl.Config.set_fmt_str_lengths(100)

# Set the row limit to a higher value
pl.Config.set_tbl_rows(50)


# # Load Data

# In[ ]:


import gwaslab as gl


# ## Gene Annotations

# In[ ]:


hla_start = 29602238
hla_end = 33409896


# In[ ]:


def extract_gene_annotations(gtf_path, chromosome, start_pos, end_pos, 
                           build="hg38", protein_coding_only=False):
    """Extract gene annotations for Manhattan plot"""
    from gwaslab.viz_plot_regionalplot import process_gtf
    
    # Set up parameters like gwaslab does
    region = [chromosome, start_pos, end_pos]
    gtf_chr_dict = {chromosome: str(chromosome)}  # or "chr" + str(chromosome)
    
    uniq_genes, exons = process_gtf(
        gtf_path=gtf_path,
        region=region,
        region_flank_factor=0.1,  # 10% flank
        build=build,
        region_protein_coding=protein_coding_only,
        gtf_chr_dict=gtf_chr_dict,
        gtf_gene_name='gene_name'
    )
    
    # Return simplified gene annotation DataFrame
    gene_annotations = uniq_genes[['name', 'start', 'end']].copy()
    gene_annotations['mid_point'] = (gene_annotations['start'] + gene_annotations['end']) // 2
    
    return gene_annotations


# In[ ]:


hla_region_genes = pl.from_pandas(extract_gene_annotations(
    gtf_path=gl.get_path('ensembl_hg38_gtf'),
    chromosome=6,
    start_pos=hla_start,  # HLA region start
    end_pos=hla_end,    # HLA region end
    build="hg38",
    protein_coding_only=True
))


# In[ ]:


hla_region_genes = pl.concat([
    hla_region_genes,
    pl.DataFrame({
        'name': ['HLA-H', 'HLA-J', 'HLA-L', 'HLA-DRB3'],
        'start': [29887752, 30006606, 30259625, 32345227], 
        'end': [29890482, 30009539, 30261703, 32358367],
        'mid_point': [29889117, 30008072, 30260664, 32351797]
    })
]).sort('start')


# In[ ]:


hla_region_genes.filter(pl.col('name').str.contains('HLA'))


# # Anal Disease

# ## HLA and GWAS Data

# In[ ]:


afr_anal_hla_df = pl.read_csv(f'{my_bucket}/HLA/four_digit_HPV_AFR_HLA_results_freq.csv')
eur_anal_hla_df = pl.read_csv(f'{my_bucket}/HLA/four_digit_HPV_EUR_HLA_results_freq.csv')
amr_anal_hla_df = pl.read_csv(f'{my_bucket}/HLA/four_digit_HPV_AMR_HLA_results_freq.csv')
metal_anal_hla_df = pl.read_csv(f'{my_bucket}/HLA/meta_HPV_4D_1.txt', separator='\t')


# In[ ]:


eur_anal_sumstats = gl.Sumstats(f'{my_bucket}/saige_gwas/v1_redo/eur/condition__hpv/gwas_results.tsv.gz',
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
eur_anal_sumstats.basic_check()
eur_hla_anal_sumstats = pl.from_pandas(eur_anal_sumstats.data).filter(
    (pl.col('CHR')==6) &
    ((pl.col('POS')>hla_start) | (pl.col('CHR')<hla_end))
)


# In[ ]:


afr_anal_sumstats = gl.Sumstats(f'{my_bucket}/saige_gwas/v1_redo/afr/condition__hpv/gwas_results.tsv.gz',
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
afr_anal_sumstats.basic_check()
afr_hla_anal_sumstats = pl.from_pandas(afr_anal_sumstats.data).filter(
    (pl.col('CHR')==6) &
    ((pl.col('POS')>hla_start) | (pl.col('CHR')<hla_end))
)


# In[ ]:


amr_anal_sumstats = gl.Sumstats(f'{my_bucket}/saige_gwas/v1_redo/amr/condition__hpv/gwas_results.tsv.gz',
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
amr_anal_sumstats.basic_check()
amr_hla_anal_sumstats = pl.from_pandas(amr_anal_sumstats.data).filter(
    (pl.col('CHR')==6) &
    ((pl.col('POS')>hla_start) | (pl.col('CHR')<hla_end))
)


# In[ ]:


def add_chr_pos_from_snpid(df: pd.DataFrame, snp_col: str = "SNPID") -> pd.DataFrame:
    """Extract CHR and POS from SNPID column (format CHR-POS-...)."""
    parts = df[snp_col].str.split("-", n=2, expand=True)
    df = df.copy()
    df["CHR"] = pd.to_numeric(parts[0], errors="coerce")
    df["POS"] = pd.to_numeric(parts[1], errors="coerce")
    return df

metal_anal_sumstats = gl.Sumstats(f'{my_bucket}/saige_gwas/v1_redo/metal/condition__hpv/condition__hpv_afr_amr_eur1.tsv',
                           fmt="metal",
                           snpid='MarkerName',
                           rsid=None,
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
metal_anal_sumstats.data = add_chr_pos_from_snpid(metal_anal_sumstats.data, snp_col='SNPID')
metal_anal_sumstats.basic_check()
metal_hla_anal_sumstats = pl.from_pandas(metal_anal_sumstats.data).filter(
    (pl.col('CHR')==6) &
    ((pl.col('POS')>hla_start) | (pl.col('CHR')<hla_end))
)


# ## Plot

# In[ ]:


metal_anal_hla_df = metal_anal_hla_df.rename({'MarkerName':'Marker', 'P-value': 'P'})


# In[ ]:


hla_genes_positions = hla_region_genes.filter(pl.col('name').str.contains('HLA')).sort('name')
    
hla_genes_positions = hla_genes_positions.with_columns(
    pl.col("name").str.split("-").list.last().alias("hla_name"),
)

locus_positions = dict(zip(hla_genes_positions['hla_name'], hla_genes_positions['mid_point']))


# In[ ]:


# Remove DMA and DRA from analysis because AC == 1 for both
eur_anal_hla_df = eur_anal_hla_df.filter(~pl.col('Marker').str.contains('DMA|DRA'))
afr_anal_hla_df = afr_anal_hla_df.filter(~pl.col('Marker').str.contains('DMA|DRA'))
amr_anal_hla_df = amr_anal_hla_df.filter(~pl.col('Marker').str.contains('DMA|DRA'))
metal_anal_hla_df = metal_anal_hla_df.filter(~pl.col('Marker').str.contains('DMA|DRA'))


# In[ ]:


# Prepare your data dictionary
data_dict = {
    'eur': (eur_anal_hla_df, eur_hla_anal_sumstats),
    'afr': (afr_anal_hla_df, afr_hla_anal_sumstats), 
    'amr': (amr_anal_hla_df, amr_hla_anal_sumstats),
    'metal': (metal_anal_hla_df, metal_hla_anal_sumstats)
}

# Calculate Bonferroni thresholds
bonf_thresholds = {
    'eur': 0.05 / len(eur_anal_hla_df),
    'afr': 0.05 / len(afr_anal_hla_df),
    'amr': 0.05 / len(amr_anal_hla_df),
    'metal': 0.05 / len(metal_anal_hla_df)
}


# ### GWAS Overlay

# In[ ]:


import matplotlib.colors as mcolors

def create_hla_manhattan_plot(data_dict, hla_region_genes, locus_positions, bonf_thresholds, figsize=(12, 12)):
    """
    Create a multi-panel HLA Manhattan plot with gene annotations
    """
    
    # Define HLA class I and II genes for color mapping
    hla_class_i = ['A', 'B', 'C', 'E', 'F', 'G', 'H', 'J', 'L']
    hla_class_ii = [ 'DMB', 'DOA', 'DOB', 'DPA1',                         # 'DMA',
                    'DPB1', 'DQA1', 'DQB1', 'DRB1', 'DRB3', 'DRB5']       #'DRA', 
 
    # Color palettes for HLA classes - using darker colors
    class_i_colors = sns.color_palette("Reds", n_colors=len(hla_class_i)+3)[3:]  # Skip 3 lightest
    class_ii_colors = sns.color_palette("Blues", n_colors=len(hla_class_ii)+3)[3:]  # Skip 3 lightest
    
    # Gene → color mapping
    hla_gene_colors = {g: c for g, c in zip(hla_class_i, class_i_colors)}
    hla_gene_colors.update({g: c for g, c in zip(hla_class_ii, class_ii_colors)})
    hla_gene_colors = {g: mcolors.to_hex(c) for g, c in hla_gene_colors.items()}

    # Hard-coded label positions to avoid overlap
    label_positions = {
        'HLA-F': 29.728,      
        'HLA-G': 29.825,      
        'HLA-H': 29.887,      
        'HLA-A': 29.95,      
        'HLA-J': 30.005,      
        'HLA-B': 31.36,        
        'HLA-C': 31.27,      
#         'HLA-E': 30.495,      
        'HLA-L': 30.26,      
        'HLA-DRB3': 32.325,   
#         'HLA-DRA': 32.413,    
        'HLA-DRB5': 32.498,   
        'HLA-DRB1': 32.565,   
        'HLA-DQA1': 32.638,   
        'HLA-DQB1': 32.688,   
        'HLA-DOB': 32.815,    
        'HLA-DMB': 32.905,    
#         'HLA-DMA': 32.962,    
        'HLA-DOA': 33.008,    
        'HLA-DPA1': 33.125,   
        'HLA-DPB1': 33.118,   
    }
    
    # Create figure with subplots
    fig, axes = plt.subplots(4, 1, figsize=figsize, sharex=True)
    fig.suptitle('HLA Region Association Results by Ancestry', fontsize=16, y=0.98)
    
    ancestries = ['eur', 'afr', 'amr', 'metal']
    ancestry_labels = ['EUR', 'AFR', 'AMR', 'METAL']
    
    for idx, (ancestry, label) in enumerate(zip(ancestries, ancestry_labels)):
        ax = axes[idx]
        hla_results, hla_anal_sumstats = data_dict[ancestry]
        bonferroni = bonf_thresholds[ancestry]
        
        # Plot background GWAS results
        ax.scatter(hla_anal_sumstats['POS'], -np.log10(hla_anal_sumstats['P']),
                   c='lightgray', alpha=0.4, s=8, zorder=1)

        df = (
            hla_results
            .with_columns([
                pl.col("Marker").str.split(".").list.first().alias("HLA_gene"),
                pl.when(pl.col("P") > 0)
                  .then(-pl.col("P").log10())
                  .otherwise(None)
                  .alias("neg_log10P"),
            ])
            .with_columns(
                pl.col("HLA_gene").replace(locus_positions, default=None).alias("pos")
            )
        )

        # Map HLA gene → color/alpha/size
        df = df.with_columns([
            pl.col("HLA_gene").replace(hla_gene_colors, default="darkgray").alias("color"),
            pl.when(pl.col("HLA_gene").is_in(list(hla_gene_colors)))
              .then(1).otherwise(0.6).alias("alpha"),
            pl.when(pl.col("HLA_gene").is_in(list(hla_gene_colors)))
              .then(20).otherwise(8).alias("size"),
        ])

        # Regular scatter plot 
        ax.scatter(
            df["pos"].to_numpy(),
            df["neg_log10P"].to_numpy(),
            c=df["color"].to_list(),
            alpha=df["alpha"].to_numpy(),
            s=df["size"].to_numpy(),
            zorder=3,
            edgecolors="white",
            linewidth=0.5
        )

        # Gene track below x-axis
        gene_track_y = -2.5
        
        for gene in hla_region_genes.iter_rows(named=True):
            gene_width = gene['end'] - gene['start']
            if gene['name'].startswith('HLA-'):
                hla_gene = gene['name']
                if hla_gene in hla_gene_colors:
                    color, alpha = hla_gene_colors[hla_gene], 0.8
                else:
                    color, alpha = 'gray', 0.6
                    
                rect = Rectangle((gene['start'], gene_track_y), gene_width, 0.6,
                                 facecolor=color, alpha=alpha, zorder=2)
                ax.add_patch(rect)
                
                # Draw label at fixed position if we have one
                if hla_gene in label_positions:
                    label_x = label_positions[hla_gene] * 1e6  # Convert Mb to bp
                    original_x = gene['mid_point']
                    
                    # Check if it's an HLA-D gene for alternation logic
                    if hla_gene.startswith('HLA-D'):
                        # Alternate HLA-D genes: even indices on top, odd on bottom
                        gene_idx = list(label_positions.keys()).index(hla_gene)
                        if gene_idx % 2 == 0:  # Top
                            label_y = gene_track_y + .9
                            line_start_y = gene_track_y + 0.6
                            line_end_y = gene_track_y + 0.8
                            va = 'bottom'
                        else:  # Bottom
                            label_y = gene_track_y - .3
                            line_start_y = gene_track_y
                            line_end_y = gene_track_y - 0.4
                            va = 'top'
                        
                        # Draw connecting line for HLA-D genes only
                        ax.plot([original_x, label_x], 
                               [line_start_y, line_end_y], 
                               color='gray', linestyle='-', linewidth=0.5, alpha=0.7)
                    else:
                        # All non-HLA-D genes go to bottom (no connecting lines)
                        label_y = gene_track_y - .3
                        va = 'top'
                    
                    # Draw label
                    ax.text(label_x, label_y, hla_gene[4:],
                            ha='center', va=va, fontsize=8, 
                            color='black')
                    
            else:
                rect = Rectangle((gene['start'], gene_track_y), gene_width, 0.6,
                                 facecolor='lightgray', alpha=0.3, zorder=1)
                ax.add_patch(rect)

        # Styling
        ax.set_ylabel(f'-log₁₀(P)\n{label}', fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(hla_region_genes['start'].min() - 100000,
                    hla_region_genes['end'].max() + 100000)
        ax.axhline(y=-np.log10(5e-8), color='gray', linestyle='--', alpha=0.75, zorder=0)
        ax.axhline(y=-np.log10(bonferroni), color='purple', linestyle='--', alpha=0.75, zorder=2)
        
        # Set y-axis with gene track below x-axis
        ax.set_ylim(-4.5, 23)
        ax.set_yticks(np.arange(0, 20.1, 2.5))
        
        # Despine and remove ticks from all but bottom plot
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        if idx < len(ancestries) - 1:  # Not the last plot
            ax.spines['bottom'].set_visible(False)
            ax.tick_params(bottom=False, labelbottom=False)  # Remove x-axis ticks and labels
    
    # axis formatting (only on last plot)
    axes[-1].set_xlabel('Chromosome 6 Position (Mb)', fontsize=12)
    axes[-1].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.1f}'))

    # Separate legends for each HLA class
    class_i_legend = plt.Line2D([0], [0], marker='o', color='w',
                               markerfacecolor=class_i_colors[1], markersize=8, alpha=0.8)
    class_ii_legend = plt.Line2D([0], [0], marker='o', color='w', 
                                markerfacecolor=class_ii_colors[1], markersize=8, alpha=0.8)

    # Create both legend objects
    legend1 = axes[0].legend([class_i_legend], ['HLA Class I (A,B,C)'],
                            loc='upper left', frameon=True, fancybox=True)
    legend2 = axes[0].legend([class_ii_legend], ['HLA Class II (DR,DQ,DP)'],
                            loc='upper right', frameon=True, fancybox=True)

    # Add the first legend back as an artist
    axes[0].add_artist(legend1)

    plt.tight_layout()
    return fig, axes


# In[ ]:


# Create the plot
fig, axes = create_hla_manhattan_plot(data_dict, hla_region_genes, locus_positions, bonf_thresholds)
plt.show()


# ### Forest Plot

# In[ ]:


from typing import Tuple

def create_hla_forest_plot(
    data_dict: dict,
    locus_positions: dict,
    bonf_thresholds: dict,
    max_markers: int = 40,
    show_marker_ids: bool = True,
    x_min: float = 0.1,
    x_max: float = 10,
) -> None:
    """
    Create forest plot comparing HLA association results between ancestries.
    """
    
    combined_data = []
    
    for ancestry, (results_df, _) in data_dict.items():
        if ancestry == 'metal':
            # Handle METAL meta-analysis format
            df_with_or = (
                results_df
                .select([
                    "Marker",
                    pl.col("Effect").alias("Beta"),
                    pl.col("StdErr").alias("SE"),
                    "P"
                ])
                .with_columns([
                    (pl.col("P") < bonf_thresholds[ancestry]).alias("bonf_significant")
                ])
                .with_columns([
                    pl.col("Beta").exp().alias("odds_ratio"),
                    (pl.col("Beta") - 1.96 * pl.col("SE")).exp().alias("odds_ratio_ci_low"),
                    (pl.col("Beta") + 1.96 * pl.col("SE")).exp().alias("odds_ratio_ci_high"),
                    pl.col("Marker").str.split(".").list.first().alias("HLA_locus"),
                    pl.lit(ancestry.upper()).alias("ancestry")
                ])
            )
        else:
            # Handle ancestry-specific format
            df_with_or = (
                results_df
                .select(["Marker", "Beta", "SE", "P", "bonf_significant"])
                .with_columns([
                    pl.col("Beta").exp().alias("odds_ratio"),
                    (pl.col("Beta") - 1.96 * pl.col("SE")).exp().alias("odds_ratio_ci_low"),
                    (pl.col("Beta") + 1.96 * pl.col("SE")).exp().alias("odds_ratio_ci_high"),
                    pl.col("Marker").str.split(".").list.first().alias("HLA_locus"),
                    pl.lit(ancestry.upper()).alias("ancestry")
                ])
            )
        combined_data.append(df_with_or)
        
    # Combine and process data
    all_data = pl.concat(combined_data)
    df_pd = all_data.to_pandas()
    
    # Color scheme by ancestry
    ancestry_colors = {
        'EUR': '#1f77b4', 'AFR': '#d62728', 'AMR': '#ff7f0e', 'METAL': '#666'
    }
    
    df_pd['color'] = df_pd['ancestry'].map(ancestry_colors).fillna('gray')
    df_pd['locus_position'] = df_pd['HLA_locus'].map(locus_positions)
    
    # Select significant markers and sort
    significant_markers = (df_pd.groupby('Marker')['bonf_significant']
                          .any().loc[lambda x: x].index.tolist())
    df_plot = df_pd[df_pd['Marker'].isin(significant_markers)].copy()
    
    marker_positions = df_plot.groupby('Marker')['locus_position'].first()
    marker_order = (marker_positions.reset_index()
                    .sort_values(['locus_position', 'Marker'])['Marker'].tolist())
    
    # Create figure with dynamic sizing
    n_markers = len(marker_order)
    figsize = (10, max(6, min(20, n_markers * 0.25 + 3)))
    fig, ax = plt.subplots(figsize=figsize)
    
    # Add backgrounds and group separators
    for i in range(len(marker_order)):
        if i % 2 == 1:
            ax.axhspan(i - 0.5, i + 0.5, color='lightgray', alpha=0.15, zorder=0)

    # Add group separators and labels
    locus_y_positions = {}
    prev_locus = None
    
    for i, marker in enumerate(marker_order):
        locus = marker.split('.')[0]
        y_pos = len(marker_order) - i - 1
        
        # Track locus positions for labels
        if locus not in locus_y_positions:
            locus_y_positions[locus] = []
        locus_y_positions[locus].append(y_pos)
        
        # Add separator line between groups (extend into label area)
        if prev_locus is not None and locus != prev_locus:
            y_line = len(marker_order) - i - 0.5
            ax.plot([-0.3, 1], [y_line, y_line], color='gray', linestyle='-',
                   alpha=0.6, linewidth=0.5, zorder=2, 
                   transform=ax.get_yaxis_transform(), clip_on=False)
        prev_locus = locus
        
    # Add locus labels
    for locus, positions in locus_y_positions.items():
        center_y = (max(positions) + min(positions)) / 2
        ax.text(-0.3, center_y, f'HLA-{locus}', 
               transform=ax.get_yaxis_transform(),
               ha='right', va='center', fontsize=9, 
               color='black')

    # Create y-positions for plotting
    ancestries = df_plot['ancestry'].unique()
    offsets = {4: [0.25, 0.08, -0.08, -0.25], 3: [-0.2, 0, 0.2], 
               2: [-0.15, 0.15]}.get(len(ancestries), [0])
    
    y_positions = {}
    for i, marker in enumerate(marker_order):
        base_y = len(marker_order) - i - 1
        for j, ancestry in enumerate(sorted(ancestries)):
            y_positions[(marker, ancestry)] = base_y + offsets[j % len(offsets)]
    
    # Plot data points
    ancestry_styles = {
        'EUR': {'marker': 'o'}, 'AFR': {'marker': 's'}, 
        'AMR': {'marker': '^'}, 'METAL': {'marker': 'D'}
    }
    
    for _, row in df_plot.iterrows():
        y_pos = y_positions.get((row['Marker'], row['ancestry']), 0)
        alpha = 1.0 if row['bonf_significant'] else 0.2
        style = ancestry_styles.get(row['ancestry'], {'marker': 'o'})
        
        # Plot CI and point
        ax.plot([row['odds_ratio_ci_low'], row['odds_ratio_ci_high']], 
               [y_pos, y_pos], color=row['color'], linewidth=2, alpha=alpha)
        ax.scatter(row['odds_ratio'], y_pos, color=row['color'], 
                  marker=style['marker'], s=20, alpha=alpha, zorder=5, 
                  edgecolors='white', linewidth=0.5)
    
    # Create legend
    legend_elements = []
    for ancestry in sorted(ancestries):
        style = ancestry_styles.get(ancestry, {'marker': 'o'})
        color = ancestry_colors.get(ancestry, 'gray')
        legend_elements.append(plt.Line2D([0], [0], color=color, linewidth=2, 
                                        marker=style['marker'], markersize=6, 
                                        alpha=1.0, label=f'{ancestry}'))
    
    legend_elements.extend([
        plt.Line2D([0], [0], color='gray', linewidth=2, alpha=1.0, 
                   label='Bonferroni Significant'),
        plt.Line2D([0], [0], color='gray', linewidth=2, alpha=0.2, 
                   label='Non-significant')
    ])

    # Customize plot
    ax.set_xscale('log')
    ax.axvline(x=1, color='black', linestyle='-', alpha=0.8, linewidth=2, zorder=1)
    
    # Add reference lines
    for ref_or in [0.5, 1.5, 2, 3, 5]:
        if x_min <= ref_or <= x_max:
            ax.axvline(x=ref_or, color='lightgray', linestyle=':', 
                      alpha=0.5, linewidth=1)
    
    # Format axes
    ax.set_ylim(-0.5, len(marker_order) - 0.5)
    ax.set_yticks(range(len(marker_order)))
    ax.set_yticklabels(list(reversed(marker_order)), fontsize=9)
    
    ax.set_xlim(x_min, x_max)
    x_ticks = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50]
    x_ticks_filtered = [tick for tick in x_ticks if x_min <= tick <= x_max]
    ax.set_xticks(x_ticks_filtered)
    ax.set_xticklabels([str(tick) for tick in x_ticks_filtered])
    
    # Labels and styling
    ax.set_xlabel('Odds Ratio', fontsize=10)
    ax.set_ylabel('', fontsize=12)
    ax.set_title(f'HLA Association Forest Plot\n({len(marker_order)} Alleles Significant in at Least One Group)', 
                fontsize=10)
    
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), 
             loc='upper left', fontsize=8)
    ax.grid(True, alpha=0.3, axis='x', which='both')
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    return fig, ax


# In[ ]:


fig, ax = create_hla_forest_plot(data_dict, locus_positions, bonf_thresholds, x_min=0.1, x_max=5)


# # Vulvovaginal Disease

# ## HLA and GWAS Data

# In[ ]:


afr_vv_hla_df = pl.read_csv(f'{my_bucket}/HLA/four_digit_HPV_AFR_VV_HLA_results_freq.csv')
eur_vv_hla_df = pl.read_csv(f'{my_bucket}/HLA/four_digit_HPV_EUR_VV_HLA_results_freq.csv')
amr_vv_hla_df = pl.read_csv(f'{my_bucket}/HLA/four_digit_HPV_AMR_VV_HLA_results_freq.csv')
metal_vv_hla_df = pl.read_csv(f'{my_bucket}/HLA/meta_VV_4D_1.txt', separator='\t')


# In[ ]:


eur_vv_sumstats = gl.Sumstats(f'{my_bucket}/saige_gwas/v2/eur/condition__hpv/gwas_results.tsv.gz',
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
eur_vv_sumstats.basic_check()
eur_hla_vv_sumstats = pl.from_pandas(eur_vv_sumstats.data).filter(
    (pl.col('CHR')==6) &
    ((pl.col('POS')>hla_start) | (pl.col('CHR')<hla_end))
)


# In[ ]:


afr_vv_sumstats = gl.Sumstats(f'{my_bucket}/saige_gwas/v2/afr/condition__hpv/gwas_results.tsv.gz',
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
afr_vv_sumstats.basic_check()
afr_hla_vv_sumstats = pl.from_pandas(afr_vv_sumstats.data).filter(
    (pl.col('CHR')==6) &
    ((pl.col('POS')>hla_start) | (pl.col('CHR')<hla_end))
)


# In[ ]:


amr_vv_sumstats = gl.Sumstats(f'{my_bucket}/saige_gwas/v2/amr/condition__hpv/gwas_results.tsv.gz',
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
amr_vv_sumstats.basic_check()
amr_hla_vv_sumstats = pl.from_pandas(amr_vv_sumstats.data).filter(
    (pl.col('CHR')==6) &
    ((pl.col('POS')>hla_start) | (pl.col('CHR')<hla_end))
)


# In[ ]:


def add_chr_pos_from_snpid(df: pd.DataFrame, snp_col: str = "SNPID") -> pd.DataFrame:
    """Extract CHR and POS from SNPID column (format CHR-POS-...)."""
    parts = df[snp_col].str.split("-", n=2, expand=True)
    df = df.copy()
    df["CHR"] = pd.to_numeric(parts[0], errors="coerce")
    df["POS"] = pd.to_numeric(parts[1], errors="coerce")
    return df

metal_vv_sumstats = gl.Sumstats(f'{my_bucket}/saige_gwas/v2/metal/condition__hpv/condition__hpv_afr_amr_eur1.tsv',
                           fmt="metal",
                           snpid='MarkerName',
                           rsid=None,
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
metal_vv_sumstats.data = add_chr_pos_from_snpid(metal_vv_sumstats.data, snp_col='SNPID')
metal_vv_sumstats.basic_check()
metal_hla_vv_sumstats = pl.from_pandas(metal_vv_sumstats.data).filter(
    (pl.col('CHR')==6) &
    ((pl.col('POS')>hla_start) | (pl.col('CHR')<hla_end))
)


# ## Plot

# In[ ]:


metal_vv_hla_df = metal_vv_hla_df.rename({'MarkerName':'Marker', 'P-value': 'P'})


# In[ ]:


hla_genes_positions = hla_region_genes.filter(pl.col('name').str.contains('HLA')).sort('name')
    
hla_genes_positions = hla_genes_positions.with_columns(
    pl.col("name").str.split("-").list.last().alias("hla_name"),
)

locus_positions = dict(zip(hla_genes_positions['hla_name'], hla_genes_positions['mid_point']))


# In[ ]:


# Remove DMA and DRA from analysis because AC == 1 for both
eur_vv_hla_df = eur_vv_hla_df.filter(~pl.col('Marker').str.contains('DMA|DRA'))
afr_vv_hla_df = afr_vv_hla_df.filter(~pl.col('Marker').str.contains('DMA|DRA'))
amr_vv_hla_df = amr_vv_hla_df.filter(~pl.col('Marker').str.contains('DMA|DRA'))
metal_vv_hla_df = metal_vv_hla_df.filter(~pl.col('Marker').str.contains('DMA|DRA'))


# In[ ]:


# Prepare your data dictionary
data_dict = {
    'eur': (eur_vv_hla_df, eur_hla_vv_sumstats),
    'afr': (afr_vv_hla_df, afr_hla_vv_sumstats), 
    'amr': (amr_vv_hla_df, amr_hla_vv_sumstats),
    'metal': (metal_vv_hla_df, metal_hla_vv_sumstats)
}

# Calculate Bonferroni thresholds
bonf_thresholds = {
    'eur': 0.05 / len(eur_vv_hla_df),
    'afr': 0.05 / len(afr_vv_hla_df),
    'amr': 0.05 / len(amr_vv_hla_df),
    'metal': 0.05 / len(metal_vv_hla_df)
}


# ### GWAS Overlay

# In[ ]:


# Create the plot
fig, axes = create_hla_manhattan_plot(data_dict, hla_region_genes, locus_positions, bonf_thresholds)
plt.show()


# ### Forest Plot

# In[ ]:


fig, ax = create_hla_forest_plot(data_dict, locus_positions, bonf_thresholds, x_min=0.1, x_max=5)

