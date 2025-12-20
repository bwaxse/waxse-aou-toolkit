#!/usr/bin/env python
# coding: utf-8

# # Load and Setup

# In[ ]:


# Load libraries
from IPython.display import display, HTML
import pandas as pd
import polars as pl
import scipy
import numpy as np
import copy

import os
import subprocess

version = get_ipython().run_line_magic('env', 'WORKSPACE_CDR')
my_bucket = os.getenv('WORKSPACE_BUCKET')
src_bucket = '{bucket}'  # Reference files (Huan Mo)
henry_bucket = '{bucket}' # Original SAIGE implementation (Henry) 


# In[ ]:


my_bucket


# In[ ]:


src_bucket


# In[ ]:


henry_bucket


# In[ ]:


ancestries_considered = ['eur', 'afr', 'amr', 'eas', 'sas']

# GWAS info
traits = {
    'condition__hpv': 'binary'
}
covariates = ['imputed_sex', 'age_at_last_ehr'] + ['ancPC{}'.format(str(x)) for x in range(1, 21)]
covariates_discrete = []

# Columns to manipulate
covariates_binarize = ['imputed_sex::F']


# In[ ]:


## output folders
output_folder = f'{my_bucket}/saige_gwas'


# # Copy Henry's Data

# In[ ]:


# ! gsutil cp {henry_bucket}/data/cohort_metadata__pcs__phenotypes.tsv.gz {src_bucket}/data/


# In[ ]:


## For stg303 data, Henry's folders contain many files/results that we don't need
# ! gsutil ls {henry_bucket}/


# In[ ]:


## my other bucket with useful data from Henry
# ! gsutil -m cp -r {bucket}/data/stg303 {src_bucket}/data/


# In[ ]:


get_ipython().system(' gsutil ls {src_bucket}/data/stg303/eur/')


# In[ ]:


get_ipython().system(' gsutil ls {henry_bucket}/data/stg105/')


# In[ ]:


# ! gsutil ls {henry_bucket}/data/stg105/eur/


# In[ ]:


# ! gsutil -m cp -r {henry_bucket}/data/stg105/ {src_bucket}/data/


# In[ ]:


## This folders are the reference folders that we need to run
get_ipython().system(' gsutil ls {src_bucket}/data/')


# In[ ]:


hpv_df = pd.read_csv('hpv_gwas_cohort.csv')


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
    machine_type,
    envs,
    in_params,
    out_params,
    boot_disk = 100,
    disk_size = 150,
    image = 'us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.2.14',
    script = 'run_impute.sh',
    preemptible = True
):
    
    # get useful info
    dsub_user_name = os.getenv("OWNER_EMAIL").split('@')[0]
    user_name = os.getenv("OWNER_EMAIL").split('@')[0].replace('.','-')

    
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
    dsub_cmd += '--name "{}" '.format(machine_type)
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


# # Merge pruned variants for step1 input

# In[ ]:


# %%writefile merge_pruned_genotypes.sh

# #!/bin/bash

# in_base=$(echo $INPUT_CHR1_BED | sed 's/chr1.bed//g');
# out_base=$(echo $OUTPUT_BED | sed 's/.bed//g');

# merge_lst='merge_input_beds.txt';
# touch $merge_lst;
# for i in ${in_base}*.bed; do
#     echo $i | sed 's/.bed//g' >> $merge_lst; 
# done;

# plink2 \
#     --pmerge-list merge_input_beds.txt bfile \
#     --indiv-sort none \
#     --delete-pmerge-result \
#     --remove $REF_SAMPLES \
#     --make-bed \
#     --out $out_base;


# In[ ]:


# def merge_bed_files(
#     anc,
#     script,
#     chroms = list(range(1,23))
# ):
    
#     # get base files
#     mb = os.getenv('WORKSPACE_BUCKET')
#     in_base = '{}/data/stg201/pruned_genotypes/{}/genotypes_chr{{}}'.format(mb, anc)
    
#     # base out
#     out_dir = '{}/data/stg303/{}'.format(mb, anc)
    
#     in_dict = {'REF_SAMPLES': '{}/data/stg105/{}/reference_samples.txt'.format(mb, anc)}
#     for chrom in chroms:
#         in_dict['INPUT_CHR{}_BED'.format(chrom)] = in_base.format(chrom) + '.bed'
#         in_dict['INPUT_CHR{}_BIM'.format(chrom)] = in_base.format(chrom) + '.bim'
#         in_dict['INPUT_CHR{}_FAM'.format(chrom)] = in_base.format(chrom) + '.fam'
    
#     env_dict = {}
#     out_dict = {
#         'OUTPUT_BED': '{}/pruned_genotypes.bed'.format(out_dir),
#         'OUTPUT_BIM': '{}/pruned_genotypes.bim'.format(out_dir),
#         'OUTPUT_PED': '{}/pruned_genotypes.fam'.format(out_dir)
#     }
            
#     dsub_script(
#         machine_type = 'c4-standard-8',
#         envs = env_dict,
#         in_params = in_dict,
#         out_params = out_dict,
#         boot_disk = 100,
#         disk_size = 150,
#         image = 'us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.2.14',
#         script = script,
#         preemptible = True
#     )


# In[ ]:


# for anc in ancestries_considered:
#     _ = merge_bed_files(
#         anc,
#         'merge_pruned_genotypes.sh',
#         chroms = list(range(1,23))
#     )


# # Prepare metadata

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


## This metadata file contains both demographic covariates 
## and examples of phenotypes of interest provided by Henry
## You can join this dataframe to your phenotype tables

metadata_file = f'{src_bucket}/data/cohort_metadata__pcs__phenotypes.tsv.gz'

metadata = pd.read_csv(
    metadata_file,
    sep='\t',
    header=0,
    dtype={'person_id' : str }
)

metadata.head()


# In[ ]:


hpv_df[['person_id', 'case']]


# In[ ]:


# Create an explicit copy of the filtered dataframe
hpv_df_copy = hpv_df.copy()

# Convert person_id to string
hpv_df_copy['person_id'] = hpv_df_copy['person_id'].astype(str)


# In[ ]:


metadata = metadata.merge(
    hpv_df_copy[['person_id', 'case']].rename(columns={'case': 'condition__hpv'}),
    on='person_id',
    how='left'
)


# In[ ]:


metadata.groupby('ancestry_pred', as_index = False).agg({'person_id' : 'nunique'})


# In[ ]:


metadata.groupby('condition__hpv', as_index = False).agg({'person_id' : 'nunique'})


# In[ ]:


len(metadata[metadata['condition__hpv'].isna()])


# In[ ]:


metadata.info()


# In[ ]:


for anc in ancestries_considered:
    anc_metadata = metadata[metadata['ancestry_pred_other'] == anc]
    
    # Manipulate any columns
    # some abnormal imputed sex -- limit to M or F
    # this will make n of metadata < genotypes but saige can handle just fine
    anc_metadata = anc_metadata[anc_metadata.imputed_sex.isin(['M', 'F'])]
    
    anc_metadata = binarize_columns(anc_metadata, covariates_binarize)
#     anc_metadata = normalize_columns(anc_metadata, covariates_normalize)
    
    # write
    anc_metadata.to_csv(
        f'{output_folder}/{anc}/gwas_metadata.tsv',
        sep='\t',
        index=False,
        header=True
    )


# In[ ]:


get_ipython().system(' gsutil ls {output_folder}')


# In[ ]:


amr_df = pd.read_csv(f'{output_folder}/amr/gwas_metadata.tsv', sep='\t')


# In[ ]:


amr_df.groupby('condition__hpv', as_index=False).agg({'person_id' : 'nunique'})


# In[ ]:


eur_df = pd.read_csv(f'{output_folder}/eur/gwas_metadata.tsv', sep='\t')


# In[ ]:


eur_df.groupby('condition__hpv', as_index=False).agg({'person_id' : 'nunique'})


# In[ ]:


afr_df = pd.read_csv(f'{output_folder}/afr/gwas_metadata.tsv', sep='\t')


# In[ ]:


afr_df.groupby('condition__hpv', as_index=False).agg({'person_id' : 'nunique'})


# In[ ]:


sas_df = pd.read_csv(f'{output_folder}/sas/gwas_metadata.tsv', sep='\t')


# In[ ]:


sas_df.groupby('condition__hpv', as_index=False).agg({'person_id' : 'nunique'})


# In[ ]:


eas_df = pd.read_csv(f'{output_folder}/eas/gwas_metadata.tsv', sep='\t')


# In[ ]:


eas_df.groupby('condition__hpv', as_index=False).agg({'person_id' : 'nunique'})


# # Fit null model

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
    # get base files

    in_base = f'{output_folder}/{anc}'
    ref_base = f'{src_bucket}/data/stg303/{anc}'
    out_dir = f'{output_folder}/{anc}/{trait}'
    
    env_dict = {
        'TRAIT': trait,
        'TRAIT_TYPE': trait_type,
        'COVARIATES': ','.join(covs),
        'COVARIATES_DISCRETE': ','.join(covs_discrete),
        'THREADS': 8
    }
    
    in_dict = {
        'INPUT_BED': f'{ref_base}/pruned_genotypes.bed',
        'INPUT_BIM': f'{ref_base}/pruned_genotypes.bim',
        'INPUT_FAM': f'{ref_base}/pruned_genotypes.fam',
        'INPUT_METADATA': f'{in_base}/gwas_metadata.tsv'
    }
    
    out_dict = {
        'OUTPUT_NULL_RDA': f'{out_dir}/saige_null_model.rda',
        'OUTPUT_NULL_VARRAT': f'{out_dir}/saige_null_model.varianceRatio.txt'
    }
            
    dsub_script(
        machine_type = 'c4-standard-8',
        envs = env_dict,
        in_params = in_dict,
        out_params = out_dict,
        boot_disk = 100,
        disk_size = 150,
        image = 'wzhou88/saige:1.3.6',
        script = script,
        preemptible = False
    )


# In[ ]:


for anc in ['eur', 'afr', 'amr', 'eas', 'sas']:
    for trait in traits.keys():
        _ = run_saige_null(
            anc,
            trait,
            traits[trait],
            covariates,
            covariates_discrete,
            'run_saige_null_model.sh'
        )


# In[ ]:


check_dsub_status()


# In[ ]:


get_ipython().system(' gsutil ls {my_bucket}/dsub/logs/c4-standard-8/bwaxse/20250514/')


# # Run SAIGE hypothesis test

# In[ ]:


get_ipython().run_cell_magic('writefile', 'run_saige_chrom.sh', '\n#!/bin/bash\nstep2_SPAtests.R \\\n    --vcfFile="${INPUT_VCF}" \\\n    --vcfFileIndex="${INPUT_VCF_IX}" \\\n    --vcfField="DS" \\\n    --chrom="${CHR}" \\\n    --is_imputed_data="TRUE" \\\n    --AlleleOrder="alt-first" \\\n    --GMMATmodelFile="${INPUT_NULL_RDA}" \\\n    --varianceRatioFile="${INPUT_NULL_VARRAT}" \\\n    --is_Firth_beta="TRUE" \\\n    --pCutoffforFirth="0.05" \\\n    --minMAC=20 \\\n    --is_output_moreDetails="TRUE" \\\n    --SAIGEOutputFile="${OUTPUT_FILE}" \\\n    --LOCO="TRUE";\n')


# In[ ]:


def run_saige_test(
    anc,
    trait,
    script,
    chroms = range(1, 23)
):
    # get base files
    mb = os.getenv('WORKSPACE_BUCKET')
    vcf_in_base = f'{src_bucket}/data/stg105/{anc}/targimp_genotypes_chr{{}}.vcf.gz'
    nm_in_base = f'{output_folder}/{anc}/{trait}'
    out_dir = f'{output_folder}/{anc}/{trait}/swarm_gwas'
    
    # get existing files to avoid repeat work
    efs = [ x.split('/')[-1].replace('.txt', '') for x in get_file_list(out_dir) if x.endswith('.txt') ]
    
    for chrom in chroms:
        if 'gwas_results_chr{}'.format(chrom) not in efs:
            env_dict = {'CHR': chrom}

            in_dict = {
                'INPUT_VCF': vcf_in_base.format(str(chrom)),
                'INPUT_VCF_IX': vcf_in_base.format(str(chrom)) + '.csi',
                'INPUT_NULL_RDA': f'{nm_in_base}/saige_null_model.rda',
                'INPUT_NULL_VARRAT': f'{nm_in_base}/saige_null_model.varianceRatio.txt'
            }

            out_dict = {
                'OUTPUT_FILE': f'{out_dir}/gwas_results_chr{chrom}.txt',
                'OUTPUT_FILE_INDEX': f'{out_dir}/gwas_results_chr{chrom}.txt.index'
            }

            dsub_script(
                machine_type = 'c4-standard-8',
                envs = env_dict,
                in_params = in_dict,
                out_params = out_dict,
                boot_disk = 100,
                disk_size = 150,
                image = 'wzhou88/saige:1.3.6',
                script = script,
                preemptible = True
            )


# In[ ]:


for anc in ['eur', 'afr', 'amr', 'sas', 'eas']:
    for trait in traits.keys():
        _ = run_saige_test(
            anc,
            trait,
            'run_saige_chrom.sh',
            chroms = range(1, 23)
        )


# In[ ]:


# Check All Statuses
check_dsub_status(full=False)


# In[ ]:


get_ipython().system(" dstat --provider google-cls-v2 --project terra-vpc-sc-e05b4e1b --location us-central1  --jobs 'c4-standar--bwaxse--250506-190742-82' --users 'bwaxse' --status '*'")


# In[ ]:


get_ipython().system(' gsutil ls {my_bucket}/saige_gwas/afr/condition__hpv/swarm_gwas/')

