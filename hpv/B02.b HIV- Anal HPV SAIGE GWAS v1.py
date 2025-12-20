#!/usr/bin/env python
# coding: utf-8

# # Load and Setup

# In[ ]:


# Load libraries
from IPython.display import display, HTML
import pandas as pd
import fastparquet
import polars as pl
import scipy
import numpy as np
import copy

import os
import subprocess

version = get_ipython().run_line_magic('env', 'WORKSPACE_CDR')
my_bucket = os.getenv('WORKSPACE_BUCKET')
src_bucket = '{bucket}' # ancestry-specific PCs

# huan_bucket = '{bucket}'  # Reference files (Huan Mo)
# henry_bucket = '{bucket}' # Original SAIGE implementation (Henry) 


# In[ ]:


my_bucket


# In[ ]:


src_bucket


# In[ ]:


# !gsutil ls {src_bucket}/data/stg003/pruned_genotypes/afr/


# In[ ]:


ancestries_considered = ['eur', 'afr', 'amr', 'eas', 'sas', 'mid']

# GWAS info
traits = {
    'condition__hpv': 'binary'
}
covariates = ['sex_binary', 'end_of_study_age'] + ['ancPC{}'.format(str(x)) for x in range(1, 21)]
all_covariates = ['sex_binary', 'end_of_study_age'] + ['PC{}'.format(str(x)) for x in range(1, 17)]
covariates_discrete = []


# In[ ]:


## output folders
output_folder = f'{my_bucket}/saige_gwas/v1'


# # Load Data

# In[ ]:


# !gsutil -u $GOOGLE_PROJECT -m cp gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/qc/genomic_metrics.tsv .


# In[ ]:


sex_df = pd.read_csv(f'genomic_metrics.tsv', 
                     sep='\t',
                     dtype={'research_id': str},
                     usecols=['research_id', 'sex_at_birth', 'dragen_sex_ploidy'])


# In[ ]:


hpv_df = pd.read_csv(f'{my_bucket}/data/cohorts/hiv_negative_hpv_gwas_cohort_v1_demographics.tsv', 
                     sep='\t',
                     dtype={'person_id': str})


# In[ ]:


pc_df = pd.read_parquet(f'{src_bucket}/data/ancestry_specific_pcs.parquet').astype({'research_id': str})


# # Function

# In[ ]:


def get_file_list(query_dir):
    tmp = subprocess.run(
        f'gsutil ls {query_dir}',
        shell=True,
        capture_output=True
    )
    files = tmp.stdout.decode('utf-8').split('\n')
    return(files)

def dsub_script(
    label,
    machine_type,
    envs,
    in_params,
    out_params,
    boot_disk = 100,
    disk_size = 150,
    image = 'us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.2.14',
    script = 'run_saige_null_model.sh',
    preemptible = True
):
    
    # get useful info
    dsub_user_name = os.getenv("OWNER_EMAIL").split('@')[0]
    user_name = os.getenv("OWNER_EMAIL").split('@')[0].replace('.','-')

    job_name = f'{label}_{script.replace(".sh", "")}'
    
    dsub_cmd = 'dsub '
    dsub_cmd += '--provider google-cls-v2 '
    dsub_cmd += '--machine-type "{}" '.format(machine_type)
    
    if preemptible:
        dsub_cmd += '--preemptible '
        
    if 'c4' in machine_type:
        # c4 doesn't use pd-ssd
        dsub_cmd += '--disk-type "hyperdisk-balanced" '
    else:
        dsub_cmd += '--disk-type "pd-ssd" '
        
    dsub_cmd += '--boot-disk-size {} '.format(boot_disk)
    dsub_cmd += '--disk-size {} '.format(disk_size)
    dsub_cmd += '--user-project "${GOOGLE_PROJECT}" '
    dsub_cmd += '--project "${GOOGLE_PROJECT}" '
    dsub_cmd += '--image "{}" '.format(image)
    dsub_cmd += '--network "network" '
    dsub_cmd += '--subnetwork "subnetwork" '
    dsub_cmd += '--service-account "$(gcloud config get-value account)" '
    dsub_cmd += '--user "{}" '.format(dsub_user_name)
    dsub_cmd += '--logging "${WORKSPACE_BUCKET}/dsub/logs/{job-name}/{user-id}/$(date +\'%Y%m%d\')/{job-id}-{task-id}-{task-attempt}.log" '
    dsub_cmd += ' "$@" '
    dsub_cmd += '--name "{}" '.format(job_name)
    dsub_cmd += '--env GOOGLE_PROJECT="${GOOGLE_PROJECT}" '
    dsub_cmd += '--script "{}" '.format(script)
    
    # Assign any environmental conditions
    for env_key in envs.keys():
        dsub_cmd += '--env {}="{}" '.format(env_key, envs[env_key])
        
    # Assign any inputs
    for in_key in in_params.keys():
        dsub_cmd += '--input {}="{}" '.format(in_key, in_params[in_key])
        
    # Assign any outputs
    for out_key in out_params.keys():
        dsub_cmd += '--output {}="{}" '.format(out_key, out_params[out_key])
        
    os.system(dsub_cmd)
    print()


# In[ ]:


def check_dsub_status(user=None, full=False):
    """Check status of dsub jobs for the specified user"""
    if user is None:
        # Get current user if not specified
        user = os.getenv("OWNER_EMAIL").split('@')[0]
    
    project = os.getenv("GOOGLE_PROJECT")

    if full:
        make_full = ' --full'
    else:
        make_full = ''
    
    cmd = f"dstat --provider google-cls-v2 --user {user} --status '*' --project {project}{make_full}"
    # cmd = f"ddel --provider google-cls-v2 --project terra-vpc-sc-840afe1e --location us-central1 --jobs 'transances--bwaxse--250319-022343-75' --users 'bwaxse'"
    print(f"Running: {cmd}")
    return subprocess.run(cmd, shell=True, capture_output=False)


# In[ ]:


def cancel_running_jobs():
    """Cancel only running/pending jobs (safer)"""
    project = os.getenv("GOOGLE_PROJECT")
    
    # Cancel only running jobs
    cancel_cmd = f"ddel --provider google-cls-v2 --project {project} --users 'bwaxse' --jobs '*'"
    print(f"Canceling running jobs: {cancel_cmd}")
    
    return subprocess.run(cancel_cmd, shell=True, capture_output=False)


# In[ ]:


def job_details(user=None, job=None):
    """List all jobs for the user, including failed ones"""
    project = os.getenv("GOOGLE_PROJECT")
    
    if user is None:
        user = os.getenv("OWNER_EMAIL").split('@')[0]
        
    if job is None:
        job = "'*' "
    else:
        job = f'--jobs {job} '
    
    cmd = f"dstat --provider google-cls-v2 --project {project} --user {user} --status {job}--full"
    print(f"Running: {cmd}")
    return subprocess.run(cmd, shell=True, capture_output=False)


# In[ ]:


def binarize_columns(df, col_refs):
    # Iterate through each column, binarize it. Col_ref = trait1::ref
    for col_ref in col_refs:
        col = col_ref.split('::')[0]
        ref = col_ref.split('::')[1]

        col_vals = df[col].dropna().unique()
        if ((len(col_vals) > 2) and (ref in col_vals)):
            print('The term {} contains >2 levels. Skipping binarize...'.format(
                col
            ))
            continue
        elif (ref not in col_vals):
            raise Exception('The reference {} is not a level for term {}. Please fix!'.format(
                ref, col
            ))

        col_alt = np.setdiff1d(col_vals, [ref]).tolist()[0]
        val_map = {ref: 0, col_alt: 1}
        df[col] = [ val_map[x] if not pd.isna(x) else x for x in df[col] ]
    return df

def binarize_sex_chromosomes(df, sex_col='dragen_sex_ploidy'):
    """
    Binarize sex chromosome data: XX=0 (female), XY=1 (male), others=NaN
    
    Parameters:
    df: DataFrame with sex chromosome data
    sex_col: column name containing sex chromosome calls
    
    Returns:
    DataFrame with new binary sex column
    """
    df = df.copy()
      
    # Create binary sex variable
    df['sex_binary'] = df[sex_col].map({
        'XX': 0,  # Female reference
        'XY': 1   # Male
        # All others (XO, XXY, XYY, XXX, etc.) become NaN automatically
    })
    
    print(f"\nAfter binarization:")
    print(f"XX (female): {(df['sex_binary'] == 0).sum():,}")
    print(f"XY (male): {(df['sex_binary'] == 1).sum():,}")
    print(f"Missing/Other: {df['sex_binary'].isna().sum():,}\n")
    
    return df

def _rank_inverse_normalize(data_ser, c=3.0/8):
    '''
    Inverse rank-normalize phenotype. Takes Series.
    '''
    # Shuffle by index
    orig_idx = data_ser.index

    data_ser = data_ser.loc[~pd.isnull(data_ser)]
    alg_input = data_ser.loc[np.random.permutation(data_ser.index.tolist())].copy()

    # Get rank, ties are determined by their position in the series (hence
    # why we randomised the series)
    rank = ss.rankdata(alg_input, method='ordinal')
    rank = pd.Series(rank, index=alg_input.index)

    # Convert rank to normal distribution
    norm_ser = rank.apply(
        lambda x, c, n: ss.norm.ppf((x - c) / (n - 2*c + 1)),
        n=len(rank),
        c=c
    )
    final = pd.Series(
        [ norm_ser[x] if x in norm_ser.index else pd.NA for x in orig_idx ],
        index=orig_idx
    )
    return(final)

def _zscore_standardize(data_ser):
    '''
    Standardize the phenotype using Z-score standardization. Takes Series.
    '''
    data_ser -= data_ser.mean(skipna=True)
    data_ser /= data_ser.std(skipna=True)
    return data_ser


def normalize_columns(df, col_method_dict):
    # Iterate through each column, normalize
    for col in col_method_dict.keys():
        df['{}__untransformed'.format(col)] = df[col]

        col_vals = df[col]
        method = col_method_dict[col]
        
        if method == 'inv_rank_norm':
            col_vals = _rank_inverse_normalize(col_vals, 3.0/8)
        elif method == 'standardize':
            col_vals = _zscore_standardize(col_vals)
        else:
            print('Incorrect normalization method `{}`. Skipping...'.format(method))

        df[col] = col_vals
    return df


# In[ ]:


def print_ancestry_summary(anc_metadata, anc):
    """Print clean summary of cases/controls by sex for an ancestry group"""
    
    # Get counts by case status and sex
    summary = anc_metadata.groupby(['condition__hpv', 'sex_binary']).size().reset_index(name='count')
    
    # Calculate totals and percentages
    case_total = anc_metadata[anc_metadata['condition__hpv'] == 1].shape[0]
    ctrl_total = anc_metadata[anc_metadata['condition__hpv'] == 0].shape[0]
    
    # Get male counts (sex_binary = 1)
    case_female = summary[(summary['condition__hpv'] == 1) & (summary['sex_binary'] == 0)]['count'].sum()
    ctrl_female = summary[(summary['condition__hpv'] == 0) & (summary['sex_binary'] == 0)]['count'].sum()
    
    # Calculate percentages (handle division by zero)
    case_female_pct = (case_female / case_total * 100) if case_total > 0 else 0
    ctrl_female_pct = (ctrl_female / ctrl_total * 100) if ctrl_total > 0 else 0
    
    print(f'{anc} ancestry: {case_total:,} cases ({case_female:,}, {case_female_pct:.1f}% female), {ctrl_total:,} controls ({ctrl_female:,}, {ctrl_female_pct:.1f}% female)')


# In[ ]:


def parse_pcs(df, pc_col, pc_prefix, num_pcs=20):
    """
    Parse comma-separated PC string into individual PC columns
    
    Parameters:
    df: DataFrame with PC data
    pc_col: column name containing comma-separated PC values
    pc_prefix: prefix for new column names (e.g., 'PC', 'ancPC')
    num_pcs: number of PCs to extract (default first 20)
    
    Returns:
    DataFrame with new PC columns
    """
    df = df.copy()
    
    # Convert list column to DataFrame and select first num_pcs columns
    pc_data = pd.DataFrame(df[pc_col].tolist(), index=df.index).iloc[:, :num_pcs]
    
    # Create column names: PC1, PC2, etc. or ancPC1, ancPC2, etc.
    pc_columns = [f'{pc_prefix}{i+1}' for i in range(num_pcs)]
    pc_data.columns = pc_columns
    
    # Concatenate with original dataframe
    for col in pc_columns:
        df[col] = pc_data[col]
    
    print(f"Parsed {num_pcs} PCs from {pc_col} as {pc_prefix}1-{pc_prefix}{num_pcs}\n")
    
    return df


# # Prepare metadata

# In[ ]:


# # ## This metadata file contains both demographic covariates 
# # ## and examples of phenotypes of interest provided by Henry
# # ## You can join this dataframe to your phenotype tables

# h_metadata_file = f'{src_bucket}/data/cohort_metadata__pcs__phenotypes.tsv.gz'
# h_metadata = pd.read_csv(metadata_file, sep='\t')
# h_metadata['imputed_sex'].value_counts()


# In[ ]:


hpv_df = hpv_df[['person_id', 'end_of_study_age', 'case']]


# In[ ]:


hpv_df = hpv_df.rename(columns={'case': 'condition__hpv'})


# In[ ]:


metadata = hpv_df.merge(
    sex_df[['research_id', 'dragen_sex_ploidy']],
    left_on='person_id',
    right_on='research_id',
    how='left'
).drop(columns=['research_id'])


# In[ ]:


metadata = metadata.merge(
    pc_df[['research_id', 'ancestry_pred', 'pca_features', 'ancestry_pred_other', 'anc_pca_features']],
    left_on='person_id',
    right_on='research_id',
    how='left'
).drop(columns=['research_id'])


# In[ ]:


metadata.isnull().sum()


# In[ ]:


# SAIGE can handle sample IDs not found in genomic dataset, but here we'll drop to minimize later processing
non_null_metadata = metadata.dropna(subset=['dragen_sex_ploidy', 'ancestry_pred_other', 'anc_pca_features'])


# In[ ]:


non_null_metadata.info()


# In[ ]:


for anc in ancestries_considered:
    print(f'=== Processing {anc} ancestry ===')
    anc_metadata = non_null_metadata[non_null_metadata['ancestry_pred_other'] == anc]
    
    # This sets XX=0 (female), XY=1 (male), others (e.g., XXY, X0) = NaN
    # For SAIGE Step 1 (fitNULLGLMM), missing` covariate data is excluded during null model fitting 
    # and therefore, Step 2
    anc_metadata = binarize_sex_chromosomes(anc_metadata, 'dragen_sex_ploidy')

    # If normalizing is required:
#     anc_metadata = normalize_columns(anc_metadata, covariates_normalize)
    
    # Parse PCs
    anc_metadata = parse_pcs(anc_metadata, 'anc_pca_features', 'ancPC', 20)

    print_ancestry_summary(anc_metadata, anc)
    
    # write
    out_file = f'{output_folder}/{anc}/gwas_metadata.tsv'
    anc_metadata.to_csv(out_file, sep='\t', index=False, header=True)
    
    print(f'\nSaved to {out_file}\n')


# In[ ]:


print(f'=== Processing Trans-ancestry ===')
transancestry_metadata = non_null_metadata

# This sets XX=0 (female), XY=1 (male), others (e.g., XXY, X0) = NaN
# For SAIGE Step 1 (fitNULLGLMM), missing` covariate data is excluded during null model fitting 
# and therefore, Step 2
transancestry_metadata = binarize_sex_chromosomes(transancestry_metadata, 'dragen_sex_ploidy')

# If normalizing is required:
#     transancestry_metadata = normalize_columns(transancestry_metadata, covariates_normalize)

# Parse PCs
transancestry_metadata = parse_pcs(transancestry_metadata, 'pca_features', 'PC', 16)

print_ancestry_summary(transancestry_metadata, 'all')

# write
out_file = f'{output_folder}/transancestry/gwas_metadata.tsv'
transancestry_metadata.to_csv(out_file, sep='\t', index=False, header=True)

print(f'\nSaved to {out_file}\n')


# # Fit null model

# In[ ]:


get_ipython().run_cell_magic('writefile', 'debug_saige_chrom_multi.sh', '#!/bin/bash\n\n# Set debug mode if DEBUG=true\nif [ "$DEBUG" = "true" ]; then\n    set -x  # Print all commands\n    echo "=== DEBUG MODE: Testing variable parsing ==="\nfi\n\nplink_base=$(echo $INPUT_BED | sed \'s/.bed//g\');\necho "Plink base: $plink_base"\n\n# Parse comma-separated ancestries\nIFS=\',\' read -ra ANCS <<< "$ANCESTRIES"\necho "Parsed ancestries: ${ANCS[@]}"\n\nfor anc in "${ANCS[@]}"; do\n    echo "=== Processing ancestry: $anc ==="\n    \n    # Get ancestry-specific inputs using variable indirection\n    rda_var="INPUT_NULL_RDA_${anc}"\n    varrat_var="INPUT_NULL_VARRAT_${anc}"\n    out_var="OUTPUT_FILE_${anc}"\n    \n    echo "RDA file: ${!rda_var}"\n    echo "VarRatio file: ${!varrat_var}"\n    echo "Output file: ${!out_var}"\n    \n    if [ "$DEBUG" = "true" ]; then\n        echo "=== SAIGE command that would be run ==="\n        echo "step2_SPAtests.R \\\\"\n        echo "    --bedFile=\\"${plink_base}.bed\\" \\\\"\n        echo "    --bimFile=\\"${plink_base}.bim\\" \\\\"\n        echo "    --famFile=\\"${plink_base}.fam\\" \\\\"\n        echo "    --chrom=\\"${CHR}\\" \\\\"\n        echo "    --GMMATmodelFile=\\"${!rda_var}\\" \\\\"\n        echo "    --varianceRatioFile=\\"${!varrat_var}\\" \\\\"\n        echo "    --SAIGEOutputFile=\\"${!out_var}\\""\n        echo ""\n    else\n        # Normal SAIGE execution\n        step2_SPAtests.R \\\n            --bedFile="${plink_base}.bed" \\\n            --bimFile="${plink_base}.bim" \\\n            --famFile="${plink_base}.fam" \\\n            --chrom="${CHR}" \\\n            --is_imputed_data="FALSE" \\\n            --AlleleOrder="alt-first" \\\n            --GMMATmodelFile="${!rda_var}" \\\n            --varianceRatioFile="${!varrat_var}" \\\n            --is_Firth_beta="TRUE" \\\n            --pCutoffforFirth="0.05" \\\n            --minMAC=20 \\\n            --is_output_moreDetails="TRUE" \\\n            --SAIGEOutputFile="${!out_var}" \\\n            --LOCO="TRUE";\n    fi\ndone\n')


# In[ ]:


get_ipython().run_cell_magic('writefile', 'run_saige_null_model.sh', "\n#!/bin/bash\n\nin_base=$(echo $INPUT_BED | sed 's/.bed//g');\nout_base=$(echo $OUTPUT_NULL_RDA | sed 's/.rda//g');\n\nstep1_fitNULLGLMM.R \\\n    --plinkFile=${in_base} \\\n    --phenoFile=${INPUT_METADATA} \\\n    --phenoCol=${TRAIT} \\\n    --covarColList=${COVARIATES} \\\n    --qCovarColList=${COVARIATES_DISCRETE} \\\n    --sampleIDColinphenoFile=person_id \\\n    --traitType=${TRAIT_TYPE} \\\n    --invNormalize=FALSE \\\n    --nThreads=${THREADS} \\\n    --IsOverwriteVarianceRatioFile=TRUE \\\n    --skipVarianceRatioEstimation=FALSE \\\n    --outputPrefix=${out_base} \\\n    --useSparseGRMtoFitNULL=FALSE;\n")


# In[ ]:


def run_saige_null(
    anc,
    trait,
    trait_type,
    covs,
    covs_discrete,
    script
):
    artifact_registry = os.getenv('ARTIFACT_REGISTRY_DOCKER_REPO', '')
    
    # get base files
    ref_base = f'{src_bucket}/data/stg005/merged_genotypes/{anc}_merged_genotypes'
    if anc == 'all':
        anc_folder = 'transancestry'
    else:
        anc_folder = anc
        
    out_dir = f'{output_folder}/{anc_folder}/{trait}'
    
    env_dict = {
        'TRAIT': trait,
        'TRAIT_TYPE': trait_type,
        'COVARIATES': ','.join(covs),
        'COVARIATES_DISCRETE': ','.join(covs_discrete),
        'THREADS': 8
    }
    
    in_dict = {
        'INPUT_BED': f'{ref_base}.bed',
        'INPUT_BIM': f'{ref_base}.bim',
        'INPUT_FAM': f'{ref_base}.fam',
        'INPUT_METADATA': f'{output_folder}/{anc_folder}/gwas_metadata.tsv'
    }
    
    out_dict = {
        'OUTPUT_NULL_RDA': f'{out_dir}/saige_null_model.rda',
        'OUTPUT_NULL_VARRAT': f'{out_dir}/saige_null_model.varianceRatio.txt'
    }
            
    dsub_script(
        label=f'step1_{anc}',
        machine_type = 'c4-standard-8',
        envs = env_dict,
        in_params = in_dict,
        out_params = out_dict,
        boot_disk = 100,
        disk_size = 150,
        image=f'{artifact_registry}/wzhou88/saige:1.3.6',
        script = script,
        preemptible = False
    )


# In[ ]:


# for anc in ['eur', 'afr', 'amr', 'eas', 'sas', 'mid']:
#     for trait in traits.keys():
#         _ = run_saige_null(
#             anc,
#             trait,
#             traits[trait],
#             covariates,
#             covariates_discrete,
#             'run_saige_null_model.sh'
#         )


# In[ ]:


# Trans-ancestry null model
_ = run_saige_null(
    'all',
    'condition__hpv',
    'binary',
    all_covariates,
    covariates_discrete,
    'run_saige_null_model.sh'
)


# In[ ]:


check_dsub_status()


# In[ ]:


job_details(job='all-run-sa--bwaxse--250623-175151-01')


# In[ ]:


get_ipython().system('gsutil cat {bucket}/dsub/logs/all-run-saige-null-model/bwaxse/20250623/all-run-sa--bwaxse--250623-175151-01-task-None.log')


# # Run SAIGE hypothesis test

# In[ ]:


get_ipython().run_cell_magic('writefile', 'run_saige_chrom_multi.sh', '#!/bin/bash\nplink_base=$(echo $INPUT_BED | sed \'s/\\.bed$//\');\n\n# Parse comma-separated ancestries\nIFS=\',\' read -ra ANCS <<< "$ANCESTRIES"\n\nfor anc in "${ANCS[@]}"; do\n    echo "Processing ancestry: $anc"\n    \n    # Get ancestry-specific inputs using variable indirection\n    rda_var="INPUT_NULL_RDA_${anc}"\n    varrat_var="INPUT_NULL_VARRAT_${anc}"\n    out_var="OUTPUT_FILE_${anc}"\n\n    step2_SPAtests.R \\\n        --bedFile="${plink_base}.bed" \\\n        --bimFile="${plink_base}.bim" \\\n        --famFile="${plink_base}.fam" \\\n        --chrom="chr${CHR}" \\\n        --is_imputed_data="FALSE" \\\n        --AlleleOrder="alt-first" \\\n        --GMMATmodelFile="${!rda_var}" \\\n        --varianceRatioFile="${!varrat_var}" \\\n        --is_Firth_beta="TRUE" \\\n        --pCutoffforFirth="0.05" \\\n        --minMAC=20 \\\n        --is_output_moreDetails="TRUE" \\\n        --SAIGEOutputFile="${!out_var}" \\\n        --LOCO="TRUE";\ndone\n')


# In[ ]:


def run_saige_test_multiancestry(
    ancestries,
    trait,
    script,
    chroms = range(1, 23)
):
    
    artifact_registry = os.getenv('ARTIFACT_REGISTRY_DOCKER_REPO', '')

    # All of Us ACAF plink files
    plink_base = 'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/plink_bed/chr{}'

    for chrom in chroms:
        # Process all ancestries for this chromosome
        env_dict = {
            'CHR': chrom,
            'ANCESTRIES': ','.join(ancestries),
        }
        
        in_dict = {
            'INPUT_BED': plink_base.format(chrom) + '.bed',
            'INPUT_BIM': plink_base.format(chrom) + '.bim', 
            'INPUT_FAM': plink_base.format(chrom) + '.fam'
        }
        
        # Add null models and metadata for each ancestry
        for anc in ancestries:
            in_dict[f'INPUT_NULL_RDA_{anc}'] = f'{output_folder}/{anc}/{trait}/saige_null_model.rda'
            in_dict[f'INPUT_NULL_VARRAT_{anc}'] = f'{output_folder}/{anc}/{trait}/saige_null_model.varianceRatio.txt'
            in_dict[f'INPUT_METADATA_{anc}'] = f'{output_folder}/{anc}/gwas_metadata.tsv'
        
        # Output files for each ancestry
        out_dict = {}
        for anc in ancestries:
            out_dict[f'OUTPUT_FILE_{anc}'] = f'{output_folder}/{anc}/{trait}/gwas/gwas_results_chr{chrom}.txt'
            out_dict[f'OUTPUT_FILE_{anc}_INDEX'] = f'{output_folder}/{anc}/{trait}/gwas/gwas_results_chr{chrom}.txt.index'

        dsub_script(
            label=f'step2_chr{chrom}',
            machine_type='c4-standard-16',  # Larger for multiple ancestries
            envs=env_dict,
            in_params=in_dict,
            out_params=out_dict,
            boot_disk=100,
            disk_size=1200,
            image=f'{artifact_registry}/wzhou88/saige:1.3.6',
            script=script,
            preemptible=True
        )


# In[ ]:


ancestries_considered


# In[ ]:


run_saige_test_multiancestry(ancestries_considered, 'condition__hpv', 'run_saige_chrom_multi.sh', [1, 3])


# In[ ]:


# Check All Statuses
check_dsub_status(full=False)


# In[ ]:


job_details(job='step2-chr1--bwaxse--250624-004953-11')


# In[ ]:


get_ipython().system('gsutil cat {bucket}/dsub/logs/step2-chr1-run-saige-chrom-multi/bwaxse/20250624/step2-chr1--bwaxse--250624-004953-11-task-None.log')


# In[ ]:


get_ipython().system(' gsutil ls {my_bucket}/saige_gwas/v1/afr/condition__hpv/gwas/')

