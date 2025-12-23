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
import re

import os
import subprocess

version = get_ipython().run_line_magic('env', 'WORKSPACE_CDR')
my_bucket = os.getenv('WORKSPACE_BUCKET')
src_bucket = '{bucket or my_bucket}' # ancestry-specific PCs and biallelic GTs


# In[ ]:


my_bucket


# In[ ]:


src_bucket


# In[ ]:


ancestries_considered = ['eur', 'afr', 'amr', 'eas', 'sas']

# GWAS info
traits = {
    'condition__sarcoid': 'binary'
}
# Ancestry-specific
covarColList = ['sex_binary', 'end_of_study_age'] + ['ancPC{}'.format(str(x)) for x in range(1, 21)]
qCovarColList = ['sex_binary'] # Categorical covariates to be used in the null model. 
    # All categorical covariates listed in qCovarCol must be also in covarColList,  e,g c("Sex"). 
sex_qCovarColList = [] # when running sex-specific GWAS, sex_binary must be removed from this list 
    
# Transancestry
all_covarColList = ['sex_binary', 'end_of_study_age'] + ['PC{}'.format(str(x)) for x in range(1, 17)]
all_qCovarColList = ['sex_binary'] # Categorical covariates to be used in the null model.    
    # All categorical covariates listed in qCovarCol must be also in covarColList,  e,g c("Sex"). 
sex_all_qCovarColList = [] # when running sex-specific GWAS, sex_binary must be removed from this list

# Columns to manipulate
covariates_binarize = ['imputed_sex::F']


# In[ ]:


## output folders
output_folder = f'{my_bucket}/saige_gwas/min_1'


# # Load Data

# In[ ]:


# !gsutil -u $GOOGLE_PROJECT -m cp gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/qc/genomic_metrics.tsv .


# In[ ]:


sex_df = pd.read_csv(f'genomic_metrics.tsv', 
                     sep='\t',
                     dtype={'research_id': str},
                     usecols=['research_id', 'sex_at_birth', 'dragen_sex_ploidy'])


# In[ ]:


pc_df = pd.read_parquet(f'{src_bucket}/data/ancestry_specific_pcs.parquet').astype({'research_id': str})


# In[ ]:


sarcoid_df = pd.read_csv(f'{my_bucket}/data/cohorts/all_sarcoid_cases.csv',
                     dtype={'person_id': str})
sarcoid_df.groupby('case', as_index = False).agg({'person_id' : 'nunique'})


# In[ ]:


sarcoid_without_exclusions_df = pd.read_csv(f'{my_bucket}/data/cohorts/sarcoid_cases_without_exclusions.csv',
                     dtype={'person_id': str})
sarcoid_without_exclusions_df.groupby('case', as_index = False).agg({'person_id' : 'nunique'})


# In[ ]:


sarcoid_1_rows = sarcoid_without_exclusions_df[sarcoid_without_exclusions_df['sarcoid_n']==1]
if not sarcoid_1_rows.empty:
    print(f"FYI, found {len(sarcoid_1_rows)} people with sarcoid_n==1 (i.e. this is sarcoid_n ≥ 1)")
else:
    print("No people found with sarcoid_n==1")


# ## Remove Sarcoid Cases with Exclusions From Controls 

# In[ ]:


# remove person_ids from sarcoid_without_exclusions_df that have case == 1 in sarcoid_df (i.e. with exclusions)
# these sarcoid cases had other diagnoses similar to sarcoid and should be removed

# Get person_ids that are case == 1 in sarcoid_df
case_1_in_sarcoid = set(sarcoid_df[sarcoid_df['case'] == 1]['person_id'])

# Get person_ids that are case == 1 in sarcoid_without_exclusions_df  
case_1_in_without_exclusions = set(sarcoid_without_exclusions_df[sarcoid_without_exclusions_df['case'] == 1]['person_id'])

# Find person_ids that are case == 1 in sarcoid_df but NOT case == 1 in sarcoid_without_exclusions_df
ids_to_remove = case_1_in_sarcoid - case_1_in_without_exclusions

# Filter out these person_ids from sarcoid_without_exclusions_df
filtered_sarcoid_df = sarcoid_without_exclusions_df[~sarcoid_without_exclusions_df['person_id'].isin(ids_to_remove)]


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

def gcs_file_exists(gs_path):
    """Check if specific GCS file exists using gsutil ls"""
    try:
        result = subprocess.run(f'gsutil ls {gs_path}', 
                              shell=True, capture_output=True)
        return result.returncode == 0
    except:
        return False
    
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
    dsub_cmd += '--provider google-batch '
    dsub_cmd += '--user-project "${GOOGLE_PROJECT}" '
    dsub_cmd += '--project "${GOOGLE_PROJECT}" '
    dsub_cmd += '--image "{}" '.format(image)
    dsub_cmd += '--network "global/networks/network" '
    dsub_cmd += '--subnetwork "regions/us-central1/subnetworks/subnetwork" '
    dsub_cmd += '--service-account "$(gcloud config get-value account)" '
    dsub_cmd += '--use-private-address '
    dsub_cmd += '--user "{}" '.format(dsub_user_name)
    dsub_cmd += '--regions us-central1 '
    dsub_cmd += '--logging "${WORKSPACE_BUCKET}/dsub/logs/{job-name}/{user-id}/$(date +\'%Y%m%d\')/{job-id}-{task-id}-{task-attempt}.log" '
    dsub_cmd += ' "$@" '
    dsub_cmd += '--name "{}" '.format(job_name)
    dsub_cmd += '--machine-type "{}" '.format(machine_type)
    
    if preemptible:
        dsub_cmd += '--preemptible '
        
    if 'c4' in machine_type:
        raise ValueError(
            f"c4 machine types ('{machine_type}') are not supported with dsub. "
            f"c4 requires hyperdisk-balanced boot disks, but dsub doesn't allow "
            f"setting boot disks. Use c2 or n2 instead."
        )

#        # c4 doesn't use pd-ssd
#         dsub_cmd += '--disk-type "hyperdisk-balanced" '
#     else:
#         dsub_cmd += '--disk-type "pd-ssd" '
        
    dsub_cmd += '--boot-disk-size {} '.format(boot_disk)
    dsub_cmd += '--disk-size {} '.format(disk_size)
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
    print('')


# In[ ]:


def validate_age_format(age: str) -> bool:
    """
    Validate age format for dsub dstat command
    """
    # Pattern: one or more digits followed by exactly one valid unit
    pattern = r'^\d+[smhdw]$'
    return bool(re.match(pattern, age.lower()))


# In[ ]:


def check_dsub_status(user: str = None, full: bool = False, age: str = '1d') -> subprocess.CompletedProcess:
    """
    Check status of dsub jobs for the specified user
    
    Parameters:
    -----------
    user : str, optional
        Username to check jobs for. Defaults to current user from OWNER_EMAIL
    full : bool, default False
        Include full job details in output
    age : str, default '1d'
        Maximum age of jobs to display. Format: <integer><unit>
        Units: s (seconds), m (minutes), h (hours), d (days), w (weeks)
        Examples: '3d', '12h', '30m', '7w'
        
    Returns:
    --------
    subprocess.CompletedProcess
        Result of the dstat command
        
    Examples:
    ---------
    >>> check_dsub_status(age='3d', full=True)  # Last 3 days, full details
    >>> check_dsub_status()  # Default: last day, summary view
    """
    
    if user is None:
        # Get current user if not specified
        user = os.getenv("OWNER_EMAIL").split('@')[0]
    
    project = os.getenv("GOOGLE_PROJECT")
    
    # Validate age parameter
    if age is not None:
        if not validate_age_format(age):
            raise ValueError(
                f"Invalid age format: '{age}'. "
                "Expected format: <integer><unit> where unit is one of: s, m, h, d, w. "
                "Examples: '3d', '12h', '30m', '7w'"
            )
    
    # Build command
    cmd_parts = [
        "dstat",
        "--provider google-batch",
        f"--user {user}",
        "--status '*'",
        f"--project {project}"
    ]
    
    if full:
        cmd_parts.append("--full")
    
    if age:
        cmd_parts.append(f"--age {age}")
    
    cmd = " ".join(cmd_parts)
    print(f"Running: {cmd}")
    return subprocess.run(cmd, shell=True, capture_output=False)


# In[ ]:


def cancel_running_jobs():
    """Cancel only running/pending jobs (safer)"""
    project = os.getenv("GOOGLE_PROJECT")
    
    # Cancel only running jobs
    cancel_cmd = f"ddel --provider google-batch --project {project} --users 'margaret.rencher' --jobs '*'"
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
    
    cmd = f"dstat --provider google-batch --project {project} --user {user} --status {job}--full"
    print(f"Running: {cmd}")
    return subprocess.run(cmd, shell=True, capture_output=False)


# In[ ]:


def print_ancestry_summary(anc_metadata, anc, condition):
    """Print clean summary of cases/controls by sex for an ancestry group"""
           
    # Get counts by case status and sex
    summary = anc_metadata.groupby([condition, 'sex_binary']).size().reset_index(name='count')
    
    # Calculate totals and percentages
    case_total = anc_metadata[anc_metadata[condition] == 1].shape[0]
    ctrl_total = anc_metadata[anc_metadata[condition] == 0].shape[0]
    
    # Get male counts (sex_binary = 1)
    case_female = summary[(summary[condition] == 1) & (summary['sex_binary'] == 0)]['count'].sum()
    ctrl_female = summary[(summary[condition] == 0) & (summary['sex_binary'] == 0)]['count'].sum()
    
    # Calculate percentages (handle division by zero)
    case_female_pct = (case_female / case_total * 100) if case_total > 0 else 0
    ctrl_female_pct = (ctrl_female / ctrl_total * 100) if ctrl_total > 0 else 0
    
    print(f'{anc} ancestry: {case_total:,} cases ({case_female:,}, {case_female_pct:.1f}% female), {ctrl_total:,} controls ({ctrl_female:,}, {ctrl_female_pct:.1f}% female)')


# In[ ]:


def view_dsub_logs(log_path):
    base_path = log_path.replace('.log', '')
    
    print("=== STDOUT ===")
    subprocess.run(['gsutil', 'cat', f'{base_path}-stdout.log'])
    
    print("\n=== STDERR ===") 
    subprocess.run(['gsutil', 'cat', f'{base_path}-stderr.log'])


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


# In[ ]:


# def binarize_columns(df, col_refs):
#     # Iterate through each column, binarize it. Col_ref = trait1::ref
#     for col_ref in col_refs:
#         col = col_ref.split('::')[0]
#         ref = col_ref.split('::')[1]

#         col_vals = df[col].dropna().unique()
#         if ((len(col_vals) > 2) and (ref in col_vals)):
#             print('The term {} contains >2 levels. Skipping binarize...'.format(
#                 col
#             ))
#             continue
#         elif (ref not in col_vals):
#             raise Exception('The reference {} is not a level for term {}. Please fix!'.format(
#                 ref, col
#             ))

#         col_alt = np.setdiff1d(col_vals, [ref]).tolist()[0]
#         val_map = {ref: 0, col_alt: 1}
#         df[col] = [ val_map[x] if not pd.isna(x) else x for x in df[col] ]
#     return df

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

# def _rank_inverse_normalize(data_ser, c=3.0/8):
#     '''
#     Inverse rank-normalize phenotype. Takes Series.
#     '''
#     # Shuffle by index
#     orig_idx = data_ser.index

#     data_ser = data_ser.loc[~pd.isnull(data_ser)]
#     alg_input = data_ser.loc[np.random.permutation(data_ser.index.tolist())].copy()

#     # Get rank, ties are determined by their position in the series (hence
#     # why we randomised the series)
#     rank = ss.rankdata(alg_input, method='ordinal')
#     rank = pd.Series(rank, index=alg_input.index)

#     # Convert rank to normal distribution
#     norm_ser = rank.apply(
#         lambda x, c, n: ss.norm.ppf((x - c) / (n - 2*c + 1)),
#         n=len(rank),
#         c=c
#     )
#     final = pd.Series(
#         [ norm_ser[x] if x in norm_ser.index else pd.NA for x in orig_idx ],
#         index=orig_idx
#     )
#     return(final)

# def _zscore_standardize(data_ser):
#     '''
#     Standardize the phenotype using Z-score standardization. Takes Series.
#     '''
#     data_ser -= data_ser.mean(skipna=True)
#     data_ser /= data_ser.std(skipna=True)
#     return data_ser


# def normalize_columns(df, col_method_dict):
#     # Iterate through each column, normalize
#     for col in col_method_dict.keys():
#         df['{}__untransformed'.format(col)] = df[col]

#         col_vals = df[col]
#         method = col_method_dict[col]
        
#         if method == 'inv_rank_norm':
#             col_vals = _rank_inverse_normalize(col_vals, 3.0/8)
#         elif method == 'standardize':
#             col_vals = _zscore_standardize(col_vals)
#         else:
#             print('Incorrect normalization method `{}`. Skipping...'.format(method))

#         df[col] = col_vals
#     return df


# # Prepare metadata

# In[ ]:


filtered_sarcoid_df = filtered_sarcoid_df[['person_id', 'end_of_study_age', 'case']]


# In[ ]:


filtered_sarcoid_df = filtered_sarcoid_df.rename(columns={'case': 'condition__sarcoid'})


# In[ ]:


metadata = filtered_sarcoid_df.merge(
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

    print_ancestry_summary(anc_metadata, anc, 'condition__sarcoid')
    
    # write
    out_file = f'{output_folder}/{anc}/gwas_metadata.tsv'
    anc_metadata.to_csv(out_file, sep='\t', index=False, header=True)
    
    print(f'\nSaved to {out_file}\n')


# In[ ]:


# print(f'=== Processing Trans-ancestry ===')
# transancestry_metadata = non_null_metadata

# # This sets XX=0 (female), XY=1 (male), others (e.g., XXY, X0) = NaN
# # For SAIGE Step 1 (fitNULLGLMM), missing` covariate data is excluded during null model fitting 
# # and therefore, Step 2
# transancestry_metadata = binarize_sex_chromosomes(transancestry_metadata, 'dragen_sex_ploidy')

# # If normalizing is required:
# #     transancestry_metadata = normalize_columns(transancestry_metadata, covariates_normalize)

# # Parse PCs
# transancestry_metadata = parse_pcs(transancestry_metadata, 'pca_features', 'PC', 16)

# print_ancestry_summary(transancestry_metadata, 'all', 'condition__sarcoid')

# # write
# out_file = f'{output_folder}/transancestry/gwas_metadata.tsv'
# transancestry_metadata.to_csv(out_file, sep='\t', index=False, header=True)

# print(f'\nSaved to {out_file}\n')


# # Fit null model

# In[ ]:


get_ipython().run_cell_magic('writefile', 'run_saige_null_model.sh', "\n#!/bin/bash\n\nin_base=$(echo $INPUT_BED | sed 's/.bed//g');\nout_base=$(echo $OUTPUT_NULL_RDA | sed 's/.rda//g');\n\nstep1_fitNULLGLMM.R \\\n    --plinkFile=${in_base} \\\n    --phenoFile=${INPUT_METADATA} \\\n    --phenoCol=${TRAIT} \\\n    --covarColList=${COVARCOLLIST} \\\n    --qCovarColList=${QCOVARCOLLIST} \\\n    --sampleIDColinphenoFile=person_id \\\n    --traitType=${TRAIT_TYPE} \\\n    --invNormalize=FALSE \\\n    --nThreads=${THREADS} \\\n    --IsOverwriteVarianceRatioFile=TRUE \\\n    --skipVarianceRatioEstimation=FALSE \\\n    --outputPrefix=${out_base} \\\n    --useSparseGRMtoFitNULL=FALSE;\n")


# In[ ]:


def run_saige_null(
    anc,
    trait,
    trait_type,
    covarColList,
    qCovarColList,
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

    # check if both required files already exist
    existing_files = [x.split('/')[-1] for x in get_file_list(out_dir)]
    required_files = ['saige_null_model.rda', 'saige_null_model.varianceRatio.txt']
    all_exist = all(file in existing_files for file in required_files)

    if not all_exist:
        env_dict = {
            'TRAIT': trait,
            'TRAIT_TYPE': trait_type,
            'COVARCOLLIST': ','.join(covarColList),
            'QCOVARCOLLIST': ','.join(qCovarColList),
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
            machine_type = 'n2d-standard-8', # recommended for 300K samples × 135K SNPs ≈ 40 billion genotypes,
                # ~ 15-20GB RAM
            envs = env_dict,
            in_params = in_dict,
            out_params = out_dict,
            boot_disk = 100,
            disk_size = 150,
            image=f'{artifact_registry}/wzhou88/saige:1.3.6',
            script = script,
            preemptible = True
        )
    else:
        print(f"Skipping {anc}/{trait} - null model files already exist")


# In[ ]:


for anc in ['eur', 'afr']: 
    for trait in traits.keys():
        _ = run_saige_null(
            anc,
            trait,
            traits[trait],
            covarColList,
            sex_qCovarColList if trait.startswith('fibroids') else qCovarColList,
            'run_saige_null_model_f_only.sh' if trait.startswith('fibroids') else 'run_saige_null_model.sh'
        )


# In[ ]:


check_dsub_status()


# In[ ]:


# job_details(job='step1-eur---bwaxse--250630-173928-69')


# In[ ]:


log=""
view_dsub_logs(log)


# In[ ]:


get_ipython().system('gsutil ls {output_folder}/*')


# # Run SAIGE hypothesis test

# In[ ]:


get_ipython().run_cell_magic('writefile', 'run_saige_chrom_multi.sh', '\n# Single ancestry, single chromosome SAIGE step2\necho "Processing ancestry: $ANCESTRY, chromosome: $CHR"\n\n# Extract plink base path from bed file\nplink_base=$(echo $INPUT_BED | sed \'s/\\.bed$//\')\n\necho "Using plink base: $plink_base"\n\nstep2_SPAtests.R \\\n    --bedFile="${plink_base}.bed" \\\n    --bimFile="${plink_base}.bim" \\\n    --famFile="${plink_base}.fam" \\\n    --chrom="${CHR}" \\\n    --is_imputed_data="FALSE" \\\n    --AlleleOrder="alt-first" \\\n    --GMMATmodelFile="${INPUT_NULL_RDA}" \\\n    --varianceRatioFile="${INPUT_NULL_VARRAT}" \\\n    --is_Firth_beta="TRUE" \\\n    --pCutoffforFirth="0.05" \\\n    --minMAC=20 \\\n    --is_output_moreDetails="TRUE" \\\n    --SAIGEOutputFile="${OUTPUT_FILE}" \\\n    --LOCO="TRUE"\n')


# In[ ]:


def run_saige_step2(
    ancestries,
    trait,
    script,
    chroms=range(1, 23),
    machine_type='n2d-standard-8',
    disk_size=200
):
    """
    Run SAIGE step2 parallelized by ancestry and chromosome.
    
    Args:
        ancestries: List of ancestry codes (e.g., ['eur', 'afr', 'amr'])
        trait: Trait name
        script: Path to the step2 script
        chroms: Chromosomes to process
        machine_type: VM machine type
        disk_size: Disk size in GB
    """
    
    artifact_registry = os.getenv('ARTIFACT_REGISTRY_DOCKER_REPO', '')
    jobs_submitted = 0
    jobs_skipped = 0

    # Pre-fetch all file lists and null model existence for all ancestries
    print(f"Fetching file lists for {trait}...")
    ancestry_files = {}
    null_models_exist = {}

    for anc in ancestries:
        out_dir = f'{output_folder}//{anc}/{trait}/gwas'
        null_model_file = f'{output_folder}//{anc}/{trait}/saige_null_model.rda'
        
        try:
            ancestry_files[anc] = set(x.split('/')[-1] for x in get_file_list(out_dir))
        except:
            ancestry_files[anc] = set()
            
        try:
            null_models_exist[anc] = gcs_file_exists(null_model_file)
        except:
            null_models_exist[anc] = False

        print(f"Submitting jobs for {anc}")
        existing_files = ancestry_files[anc]
        null_model_exists = null_models_exist[anc]
       
        for chrom in chroms:
            # Check if output files already exist
            required_files = {f'gwas_results_chr{chrom}.txt', f'gwas_results_chr{chrom}.txt.index'}
            results_exist = required_files.issubset(existing_files)

            if results_exist and null_model_exists:
                print(f"  Skipping {anc} chr{chrom} - results already exist")
                jobs_skipped += 1
            elif results_exist and not null_model_exists:
                print(f"  Skipping {anc} chr{chrom} - results exist but null model missing")
                jobs_skipped += 1
            elif not results_exist and not null_model_exists:
                print(f"  Skipping {anc} chr{chrom} - null model missing")
                jobs_skipped += 1
            else:  # not results_exist and null_model_exists
                # Submit the job
                plink_base = f'{src_bucket}/data/stg009/{anc}/genotypes_chr{chrom}'

                env_dict = {
                    'CHR': chrom,
                    'ANCESTRY': anc,
                }

                in_dict = {
                    'INPUT_BED': f'{plink_base}.bed',
                    'INPUT_BIM': f'{plink_base}.bim',
                    'INPUT_FAM': f'{plink_base}.fam',
                    'INPUT_NULL_RDA': f'{output_folder}//{anc}/{trait}/saige_null_model.rda',
                    'INPUT_NULL_VARRAT': f'{output_folder}//{anc}/{trait}/saige_null_model.varianceRatio.txt',
                }

                out_dict = {
                    'OUTPUT_FILE': f'{output_folder}//{anc}/{trait}/gwas/gwas_results_chr{chrom}.txt',
                    'OUTPUT_FILE_INDEX': f'{output_folder}//{anc}/{trait}/gwas/gwas_results_chr{chrom}.txt.index'
                }

                dsub_script(
                    label=f's2_{anc}_{chrom}_{trait}',
                    machine_type=machine_type,
                    envs=env_dict,
                    in_params=in_dict,
                    out_params=out_dict,
                    boot_disk=100,
                    disk_size=disk_size,
                    image=f'{artifact_registry}/wzhou88/saige:1.3.6',
                    script=script,
                    preemptible=True
                )
                jobs_submitted += 1
    
    print(f"SAIGE step2 summary: {jobs_submitted} jobs submitted, {jobs_skipped} jobs skipped")


# In[ ]:


## Run ancestry-specific Step 2 for all traits
for trait in traits.keys():
    run_saige_step2(
        ancestries=['eur', 'afr'],
        trait=trait,
        script='run_saige_chrom_multi.sh',
#         chroms=[22] # [22] would run only chr 22; not setting 'chroms' would run all
    )


# In[ ]:


# Check All Statuses
check_dsub_status(full=False)


# In[ ]:


job_details(job='s2-afr-2-c--bwaxse--250827-210514-71')


# In[ ]:


log="{bucket or my_bucket}/dsub/logs/s2-afr-2-condition--sarcoid-run-saige-chrom-multi/bwaxse/20250827/s2-afr-2-c--bwaxse--250827-210514-71-task-None.log"
view_dsub_logs(log)


# In[ ]:


get_ipython().system('gsutil cat {bucket or my_bucket}/dsub/logs/step2-eur-chr13-condition--sarcoid-run-saige-chrom-multi/bwaxse/20250701/step2-eur---bwaxse--250701-140545-85-task-None.log')


# In[ ]:


get_ipython().system(' gsutil ls {my_bucket}/saige_gwas/afr/condition__sarcoid/gwas/')

