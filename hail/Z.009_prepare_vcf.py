#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Load libraries
from IPython.display import display, HTML
import pandas as pd
import polars as pl

import os
import subprocess

version = get_ipython().run_line_magic('env', 'WORKSPACE_CDR')
my_bucket = os.getenv('WORKSPACE_BUCKET')


# # Function

# In[ ]:


# !gsutil ls {my_bucket}/data/stg009/eur/


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


def cancel_running_jobs():
    """Cancel only running/pending jobs (safer)"""
    project = os.getenv("GOOGLE_PROJECT")
    
    # Cancel only running jobs
    cancel_cmd = f"ddel --provider google-cls-v2 --project {project} --users 'bwaxse' --jobs '*'"
    print(f"Canceling running jobs: {cancel_cmd}")
    
    return subprocess.run(cancel_cmd, shell=True, capture_output=False)


# In[ ]:


def cancel_job(job_id):
    """Cancel a specific job"""
    project = os.getenv("GOOGLE_PROJECT")
    
    cmd = f"ddel --provider google-cls-v2 --project {project} --jobs {job_id}"
    print(f"Running: {cmd}")
    return subprocess.run(cmd, shell=True, capture_output=False)


# # Variant filter script and function

# In[ ]:


get_ipython().run_cell_magic('writefile', 'run_vcf_export.sh', '\n#!/bin/bash\n\n# VCF export for single chromosome\n# Input: Single chromosome pgen files by ancestry\n# Output: Single chromosome vcf file by ancestry\n\nINPUT_PGEN_BASE="${INPUT_PGEN_PGEN%.*}"  # Remove .pgen extension\necho "Derived INPUT_PGEN_BASE: $INPUT_PGEN_BASE"\n\nnthread=$(python -c "import os; print(len(os.sched_getaffinity(0)))");\necho "Running with $nthread threads";\n\nOUTPUT_PREFIX="${OUTPUT_RESULTS%\\*}"\n# Set ancestry-specific output\n\necho "Processing ancestry: $ANCESTRY, chromosome: $CHROM"\necho "Output prefix: $OUTPUT_PREFIX"\n\necho "Converting pgen to vcf"\nplink2 \\\n    --pfile $INPUT_PGEN_BASE \\\n    --export vcf bgz \\\n    --threads $nthread \\\n    --out $OUTPUT_PREFIX\n\n# Index the VCF (creates .tbi and .csi)\necho "Creating tabix index"\ntabix -p vcf ${OUTPUT_PREFIX}.vcf.gz\nbcftools index -c ${OUTPUT_PREFIX}.vcf.gz\n\necho "vcf output complete for $ANCESTRY"\necho "Files created with prefix: $OUTPUT_PREFIX"\nls -la ${OUTPUT_PREFIX}*\n')


# In[ ]:


def get_file_list(query_dir):
    tmp = subprocess.run(
        f'gsutil ls {query_dir}',
        shell=True,
        capture_output=True
    )
    files = tmp.stdout.decode('utf-8').split('\n')
    return(files)


# In[ ]:


def dsub_script(
    machine_type,
    out_base,
    anc,
    chrom=1,
    boot_disk=100,
    disk_size=100,
    memory=12000,
    script='run_vcf_export.sh'
):
    
    # get useful info
    dsub_user_name = os.getenv("OWNER_EMAIL").split('@')[0]
    user_name = os.getenv("OWNER_EMAIL").split('@')[0].replace('.', '-')

    job_name = f'{anc}_{chrom}_vcf_export'

    # Template for input files (will be substituted in script)
    my_bucket = os.getenv('WORKSPACE_BUCKET') 
    input_pgen_base = f'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/pgen/acaf_threshold.chr{chrom}'
    
    # Build dsub command
    cmd = [
        'dsub',
        '--provider', 'google-cls-v2',
        '--machine-type', machine_type,
        '--disk-type', 'pd-ssd',
        '--boot-disk-size', str(boot_disk),
        '--disk-size', str(disk_size),
        '--user-project', os.environ['GOOGLE_PROJECT'],
        '--project', os.environ['GOOGLE_PROJECT'],
        '--image', 'us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.2.14',
        '--network', 'network',
        '--subnetwork', 'subnetwork',
        '--service-account', subprocess.check_output(['gcloud', 'config', 'get-value', 'account']).decode().strip(),
        '--user', dsub_user_name,
        '--logging', f"{os.environ['WORKSPACE_BUCKET']}/dsub/logs/{{job-name}}/{{user-id}}/{{job-id}}-{{task-id}}-{{task-attempt}}.log",
        '--name', job_name,
        '--env', f'GOOGLE_PROJECT={os.environ["GOOGLE_PROJECT"]}',
        '--env', f'ANCESTRY={anc}',
        '--env', f'CHROM={chrom}',
        # Input files
        '--input', f'INPUT_PGEN_PGEN={input_pgen_base}.pgen',
        '--input', f'INPUT_PGEN_PSAM={input_pgen_base}.psam', 
        '--input', f'INPUT_PGEN_PVAR={input_pgen_base}.pvar',
        # Output files
        '--output', f'OUTPUT_RESULTS={out_base}*',
        '--script', script
    ]
            
    subprocess.run(cmd)


# In[ ]:


def run_vcf_export(
    my_bucket,
    ancestries,
    script='run_vcf_export.sh',
):
    """
    Run VCF export for each ancestry and chromosome
    Output: one vcf.gz file with index per ancestry-chromosome
    """

    # Process each chromosome for each ancestry
    for chrom in [str(i) for i in range(1, 23)] + ['X', 'Y']:
        for anc in ancestries:
            # Output directory
            out_dir = f'{my_bucket}/data/stg009/{anc}'
    
            # Check if already exists (check for a few chromosome files)
            existing_files = get_file_list(out_dir)
            check_chroms = [1, 10, 22]  # Check beginning, middle, end
            if any(f'{anc}_genotypes_chr{chrom}.vcf.gz' in f for f in existing_files):
                print(f"VCF files already exist for {anc} chr{chrom}")
                continue  # Skip this combination
        
            print(f"Starting serial vcf export for {anc} ancestry...")

            dsub_script(
                machine_type='c3-standard-8',
                out_base=f'{out_dir}/{anc}_genotypes_chr{chrom}',
                anc=anc,
                chrom=chrom,
                boot_disk=100,
                disk_size=100,
                script=script
            )


# # Run

# In[ ]:


ancestries_considered = ['eur', 'afr', 'amr', 'eas', 'sas', 'mid']


# In[ ]:


# # # Test
run_vcf_export(my_bucket, ['mid'])


# In[ ]:


run_vcf_export(my_bucket, ancestries_considered)


# # Check dsub

# In[ ]:


check_dsub_status(full=False)


# In[ ]:


job_id = 'mid-22-vcf--bwaxse--250618-175310-96'


# In[ ]:


job_details(job=job_id)


# In[ ]:


get_ipython().system('gsutil cat {bucket or my_bucket}/dsub/logs/mid-22-vcf-export/bwaxse/mid-22-vcf--bwaxse--250618-150520-32-task-None.log')


# In[ ]:


get_ipython().system('gsutil du -h {bucket or my_bucket}/data/stg009/mid/')

