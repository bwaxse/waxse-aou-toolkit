#!/usr/bin/env python
# coding: utf-8

# # METAL
# https://github.com/bwaxse/metal

# In[ ]:


# Load libraries
from google.cloud import storage
from IPython.display import display, HTML
import polars as pl
import pandas as pd
import gzip
import re
import io
from pathlib import Path
import os
import subprocess

version = get_ipython().run_line_magic('env', 'WORKSPACE_CDR')
my_bucket = os.getenv('WORKSPACE_BUCKET')


# # Function

# In[ ]:


def combine_gwas_results(my_bucket, in_dir, trait_type, output_path):
    """
    Combine chromosome-wise GWAS results from GCS.
    
    Parameters:
    my_bucket: GCS bucket (with or without gs:// prefix)
    in_dir: '{output_folder}//{anc}/{trait}/gwas'
    trait_type: 'quantitative' or 'binary'
    output_path: '{output_folder}/{anc}/{trait}/gwas_results.tsv.gz'
    """
    
    client = storage.Client()
    bucket = client.bucket(my_bucket.removeprefix('gs://'))
    
    # Column mapping
    col_dict = {
        'CHR': 'chromosome',
        'POS': 'base_pair_location', 
        'MarkerID': 'vid',
        'Allele1': 'other_allele',
        'Allele2': 'effect_allele',
        'AF_Allele2': 'effect_allele_frequency',
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
    
    # Column order
    col_order = [
        'chromosome', 'base_pair_location', 'vid', 'rsID', 'other_allele', 'effect_allele',
        'effect_allele_frequency', 'beta', 'standard_error', 'test_statistic_value',
        'p_value', 'p_value__nospa', 'spa_converged',
        'effect_allele_frequency__cases', 'effect_allele_frequency__control',
        'n', 'n_cases', 'n_controls', 'n_case__alt_homs', 'n_case__hets', 'n_controls__alt_homs', 'n_controls__hets'
    ]
        
    # Combine chromosome files
    combined_dfs = []
    
    for chrom in range(1, 23): 
        blob_path = f'{in_dir}/gwas_results_chr{chrom}.txt'
        blob = bucket.blob(blob_path)
        
        if not blob.exists():
            print(f"Skipping chromosome {chrom}: file not found")
            continue
            
        print(f"Processing chromosome {chrom}...")
        
        # Download and read
        data = blob.download_as_text()
        chrom_df = pl.read_csv(io.StringIO(data), separator='\t')
        combined_dfs.append(chrom_df)
    
    if not combined_dfs:
        raise FileNotFoundError(f"No chromosome files found in {my_bucket}/{in_dir}")
    
    # Combine all chromosomes
    print("Combining chromosomes...")
    df = pl.concat(combined_dfs)
    
    # Rename columns
    df = df.rename(col_dict)
    
    # Construct variant_id if it's missing/empty
    vid_missing_count = df.filter(pl.col('vid')=='.').height
    if vid_missing_count > 0:
        print(f"Found {vid_missing_count} missing vid values, constructing from components when missing...")
        df = df.with_columns(
            pl.when(pl.col('vid')=='.')
            .then(
                pl.concat_str([
                    pl.col('chromosome').cast(pl.Utf8),
                    pl.col('base_pair_location').cast(pl.Utf8),
                    pl.col('other_allele'),
                    pl.col('effect_allele')
                ], separator='-')
            )
            .otherwise(pl.col('vid'))
            .alias('vid')
        )
    
    # Select and reorder columns
    available_cols = [col for col in col_order if col in df.columns]
    df = df.select(available_cols)
    
    # Save to GCS
    output_blob = bucket.blob(output_path)
    
    output_buffer = io.BytesIO()
    df.write_csv(output_buffer, separator='\t')
    output_buffer.seek(0)
    compressed_output = gzip.compress(output_buffer.getvalue())
    output_blob.upload_from_string(compressed_output, content_type='application/gzip')
    
    print(f"Combined file saved to: {my_bucket}/{output_path}")
    return df


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


# In[ ]:


def validate_saige_inputs(trait, ancestries, base_output_folder):
    """
    Validate SAIGE output files for METAL meta-analysis
    """
    required_file = "gwas_results.tsv.gz"
    validated_inputs = {}
    
    print(f"Validating inputs for {trait}")
    
    for anc in ancestries:
        in_dir = f"{base_output_folder}/{anc}/{trait}"
        file_path = f"{in_dir}/{required_file}"
    
        # Check file exists
        if not gcs_file_exists(file_path):
            print(f"Missing {required_file} for ancestry {anc}")
            continue
            
        # Extract sample sizes from the file
        sample_info = extract_sample_sizes(file_path)
        if not sample_info:
            print(f"Could not extract sample sizes for ancestry {anc}")
            continue

        validated_inputs[anc] = {
            'path': file_path,
            'n_cases': sample_info['n_cases'],
            'n_controls': sample_info['n_controls'], 
            'n_total': sample_info['n_total']
        }
        
        print(f"Ancestry {anc}: {sample_info['n_cases']} cases, {sample_info['n_controls']} controls (total: {sample_info['n_total']})")
    
    if len(validated_inputs) < 2:
        print(f"Error: Need at least 2 ancestries for meta-analysis. Found: {len(validated_inputs)}")
        return None
    
    print(f"{len(validated_inputs)} validated studies available for {trait}")    
    return validated_inputs

def extract_sample_sizes(file_path):
    """
    Extract n_cases and n_controls from first data row of SAIGE file
    Returns dict with sample size info or None if failed
    """
    try:
        # Read just the first few lines to get sample sizes
        cmd = f"gsutil cat '{file_path}' | gunzip | head -2 2>/dev/null"
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        
        if not result.stdout:
            print(f"No output from command for {file_path}")
            return None
        lines = result.stdout.strip().split('\n')

        if len(lines) < 2:
            print(f"File {file_path} has insufficient data")
            return None
            
        header = lines[0].split('\t')
        data = lines[1].split('\t')
        
        # Find column indices
        try:
            n_cases_idx = header.index('n_cases')
            n_controls_idx = header.index('n_controls')
        except ValueError as e:
            print(f"Missing required columns in {file_path}: {e}")
            return None
        
        # Extract values
        n_cases = int(float(data[n_cases_idx]))
        n_controls = int(float(data[n_controls_idx]))
        n_total = n_cases + n_controls
        
        return {
            'n_cases': n_cases,
            'n_controls': n_controls,
            'n_total': n_total
        }
        
    except Exception as e:
        print(f"Error extracting sample sizes from {file_path}: {e}")
        return None


# In[ ]:


def dsub_script(
    label,
    machine_type,
    envs,
    in_params,
    out_params,
    out_dirs,
    boot_disk = 100,
    disk_size = 150,
    image = 'us.gcr.io/broad-dsp-gcr-public/terra-jupyter-aou:2.2.14',
    script = 'run_metal.sh',
    preemptible = True
):
    
    # get useful info
    dsub_user_name = os.getenv("OWNER_EMAIL").split('@')[0]
    
    user_name = os.getenv("OWNER_EMAIL").split('@')[0].replace('.','-')

    job_name = f'{label}'
    
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
    if out_params != None:
        for out_key in out_params.keys():
            dsub_cmd += '--output {}="{}" '.format(out_key, out_params[out_key])
        
    for out_key in out_dirs.keys():
        dsub_cmd += '--output-recursive {}="{}" '.format(out_key, out_dirs[out_key])

    os.system(dsub_cmd)
#     print_dsub_readable(dsub_cmd)


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
    cancel_cmd = f"ddel --provider google-batch --project {project} --users 'bwaxse' --jobs '*'"
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


def view_dsub_logs(log_path):
    base_path = log_path.replace('.log', '')
    
    print("=== STDOUT ===")
    subprocess.run(['gsutil', 'cat', f'{base_path}-stdout.log'])
    
    print("\n=== STDERR ===") 
    subprocess.run(['gsutil', 'cat', f'{base_path}-stderr.log'])


# In[ ]:


def print_dsub_readable(cmd):
    """
    Simple readable format - newline before each --
    """
    readable_cmd = cmd.replace(' --', ' \\\n    --')
    print(readable_cmd)
    print()  # Extra line for separation


# # Merge GWAS Results

# In[ ]:


get_ipython().system('gsutil ls {my_bucket}/saige_gwas/v2/amr/condition__hpv/')


# In[ ]:


ancestries_considered = ['afr', 'amr', 'eur'] 
traits = {
    'condition__hpv': 'binary',
}

init_combine = True


# In[ ]:


## output folders
output_folder = f'{my_bucket}/saige_gwas/v2' # no trailing /


# In[ ]:


if init_combine:
    for anc in ancestries_considered:
        print(f'combining {anc} files')
        for trait in traits.keys():
            rez = combine_gwas_results(
                my_bucket=my_bucket,
                in_dir=f'saige_gwas/v2/{anc}/{trait}/gwas',
                trait_type=traits[trait],
                output_path=f'saige_gwas/v2/{anc}/{trait}/gwas_results.tsv.gz'
            )


# # Setup

# In[ ]:


get_ipython().system('gsutil ls {output_folder}/amr/condition__hpv/gwas_results.tsv.gz')


# # Run METAL

# In[ ]:


get_ipython().run_cell_magic('writefile', 'run_metal.sh', '\n#!/bin/bash\nset -euo pipefail\n\necho "Starting METAL meta-analysis for ${TRAIT}"\necho "Number of input files: ${N_FILES}"\n\n# Parse ancestries\nIFS=\',\' read -ra INPUTS <<< "${INPUTS}"\n\n# Parse sample sizes\nIFS=\',\' read -ra SIZES <<< "${SAMPLE_SIZES}"\n\necho "Preprocessing files to add ancestry-specific n_total columns..."\nPROCESSED_FILES=()\n\n# Pre-processing adds sample size via n_cases and n_controls extracted from SAIGE GWAS results \nfor i in "${!INPUTS[@]}"; do\n    var_name="${INPUTS[$i]}"\n    input_file="${!var_name}"\n    n_total="${SIZES[$i]}"\n    processed_file="/tmp/processed_file_${i}.tsv"\n    \n    echo "Processing file $((i+1)): ${input_file} (n_total=${n_total})"\n    \n    # Check if file exists\n    if [[ ! -f "${input_file}" ]]; then\n        echo "ERROR: Input file not found: ${input_file}"\n        echo "Available files in /mnt/data/input/:"\n        ls -la /mnt/data/input/ || echo "Input directory not accessible"\n        exit 1\n    fi\n    \n    # Add constant n_total column for this ancestry\n    if [[ "${input_file}" == *.gz ]]; then\n        gzip -dc "${input_file}"\n    else\n        cat "${input_file}"\n    fi | awk -F\'\\t\' -v OFS=\'\\t\' -v n_total="${n_total}" \'\n        NR==1 { print $0, "n_total"; next }\n               { print $0, n_total }\' > "${processed_file}"\n    \n    PROCESSED_FILES+=("/tmp/processed_file_${i}.tsv")\n    \n    # Verify output\n    n_lines=$(wc -l < "${processed_file}")\n    echo "  Created ${processed_file} with ${n_lines} lines"\ndone\n\n# Convert processed files array to comma-separated string\nPROCESSED_FILES_STR=$(IFS=\',\'; echo "${PROCESSED_FILES[*]}")\necho "Running METAL with processed files: ${PROCESSED_FILES_STR}"\n\n# Create a METAL script dynamically\nMETAL_SCRIPT="${OUTPUT_PATH}/metal_script.txt"\ncat > "${METAL_SCRIPT}" <<EOF\nMARKER vid\nWEIGHT n_total\nALLELE effect_allele other_allele\nFREQ effect_allele_frequency\nPVAL p_value\nEFFECT beta\nSTDERR standard_error\nSEPARATOR TAB\nSCHEME STDERR\n\n# Auto-flip alleles based on frequency\nAVERAGEFREQ ON\nMINMAXFREQ ON\n\n$(for file in "${PROCESSED_FILES[@]}"; do echo "PROCESS $file"; done)\n\nOUTFILE ${OUT_PREF} .tsv\nANALYZE HETEROGENEITY\nQUIT\nEOF\n\necho "Metal script created"\n\n# Run METAL\nmetal "${METAL_SCRIPT}"\n\nmv ${OUT_PREF}1.tsv "${OUTPUT_PATH}/"\nmv ${OUT_PREF}1.tsv.info "${OUTPUT_PATH}/"\n\n# List all output files\necho "Final contents of OUTPUT_PATH (${OUTPUT_PATH}):"\nls -lh "${OUTPUT_PATH}" || echo "OUTPUT_PATH not accessible"\n\necho "METAL analysis completed successfully for ${TRAIT}"\n')


# In[ ]:


def run_metal(trait, ancestries, base_output_folder):
    """
    Run METAL meta-analysis with SAIGE-formatted inputs
    """
    artifact_registry = os.getenv('ARTIFACT_REGISTRY_DOCKER_REPO', '')
            
    # Validate SAIGE inputs
    validated_inputs = validate_saige_inputs(trait, ancestries, base_output_folder)
    if not validated_inputs:
        return None

    out_dir = f'{output_folder}/metal/{trait}'
    existing_files = set(x.split('/')[-1] for x in get_file_list(out_dir))

    # Create ancestry suffix from validated inputs
    ancestry_suffix = '_'.join(sorted(validated_inputs.keys()))

    result_file = f'{trait}_{ancestry_suffix}1.tsv'
    results_exist = {result_file}.issubset(existing_files) # True if result_file is in existing_files

    if not results_exist:
        # METAL configuration
        env_dict = {
            'TRAIT': trait,
            'OUT_PREF': f'{trait}_{ancestry_suffix}',
        }

        in_dict = {}
        file_list = []
        sample_sizes = []

        for i, (anc, file_info) in enumerate(validated_inputs.items(), 1):
            in_dict[f'INPUT_{anc.upper()}'] = file_info['path']
            file_list.append(file_info['path'])
            sample_sizes.append(str(file_info['n_total']))

        # Pass file list and sample sizes as comma-separated strings
        env_dict['INPUTS'] = ",".join(in_dict.keys())
        env_dict['SAMPLE_SIZES'] = ','.join(sample_sizes)  # n_total per ancestry
        env_dict['N_FILES'] = str(len(file_list))

        out_dirs = {
            'OUTPUT_PATH': out_dir
        }

        dsub_script(
            label=f'metal_{trait}',
            machine_type='n2d-standard-8',
            envs=env_dict,
            in_params=in_dict,
            out_params=None,
            out_dirs=out_dirs,
            boot_disk=100,
            disk_size=150,
            image=f'{artifact_registry}/bwaxse/metal',
            script='run_metal.sh',
            preemptible=True
        )
    else:
        print(f'skipped METAL meta-analysis for {trait}_{ancestry_suffix}\n Results exist in {out_dir}\n')


# ## Run METAL

# In[ ]:


# Run METAL
_ = run_metal(
    trait='condition__hpv',
    ancestries=['amr', 'afr', 'eur'],
    base_output_folder=output_folder
)


# In[ ]:


check_dsub_status(full=False)


# In[ ]:


job_details(job='metal-cond--bwaxse--250908-184321-05')


# In[ ]:


log="{bucket}/dsub/logs/metal-condition--hpv/bwaxse/20250908/metal-cond--bwaxse--250908-184321-05-task-None.log"
view_dsub_logs(log)


# In[ ]:


get_ipython().system('gsutil ls {output_folder}/metal/condition__hpv/')


# In[ ]:


# cancel_running_jobs()


# In[ ]:


# """Cancel only running/pending jobs (safer)"""
# project = os.getenv("GOOGLE_PROJECT")

# # Cancel only running jobs
# cancel_cmd = f"ddel --provider google-batch --project {project} --users 'bwaxse' --jobs 'metal-alle--bwaxse--250825-193010-56'"
# print(f"Canceling running jobs: {cancel_cmd}")

# subprocess.run(cancel_cmd, shell=True, capture_output=False)

