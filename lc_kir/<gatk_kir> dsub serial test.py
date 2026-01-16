#!/usr/bin/env python
# coding: utf-8

# # KIR-Mapper Serial Execution Test Script
# Test multiple machine types for serial (non-parallel) execution to determine cost-optimal configuration

# ## Initial setup

# In[ ]:

get_ipython().system('pip install tctk --upgrade')

# In[ ]:

from tctk import AoUTools as at
import os
import polars as pl
from datetime import datetime
from google.cloud import storage
import json
from pathlib import Path

# In[ ]:

# Environment setup
WORKSPACE_BUCKET = os.getenv('WORKSPACE_BUCKET')
WORKSPACE_CDR = os.getenv('WORKSPACE_CDR')
GOOGLE_PROJECT = os.getenv('GOOGLE_PROJECT')

# In[ ]:

# ============================================================================
# CONFIGURATION: Machine Type Comparison for Serial Execution
# ============================================================================

TEST_CONFIGURATIONS = [
    # Format: (num_samples, machine_type)

    # c2d family (AMD EPYC, compute-optimized)
    (20, 'c2d-highcpu-4'),    # Colleague's baseline machine
    (20, 'c2d-highcpu-8'),    # More cores

    # c2 family (Intel, compute-optimized)
    (20, 'c2-standard-4'),    # More memory per core

    # n2 family (Intel Cascade Lake, balanced)
    (20, 'n2-standard-4'),    # Standard balanced
]

DOCKER_IMAGE = "phetk/gatk-kirmapper:0.1"
OUTPUT_BASE = f"{WORKSPACE_BUCKET}/dsub/results/kir_serial"

# In[ ]:

# ============================================================================
# HELPER FUNCTIONS: Data Loading and Processing
# ============================================================================

def load_ancestry_metadata(bucket):
    """Load ancestry metadata and filter for EUR samples."""
    ancestry_path = f'{bucket}/data/ancestry_metadata.tsv'
    print(f"Loading ancestry metadata from {ancestry_path}...")
    ancestry_df = pl.read_csv(ancestry_path, separator='\t')
    print(f"Total samples: {len(ancestry_df)}")
    print(f"Ancestry distribution:\n{ancestry_df.group_by('ancestry_pred_other').len()}")
    return ancestry_df

def get_eur_person_ids(ancestry_df):
    """Extract EUR person IDs from ancestry dataframe."""
    eur_ids = ancestry_df.filter(pl.col('ancestry_pred_other') == 'eur')['research_id'].to_list()
    print(f"EUR samples available: {len(eur_ids)}")
    return eur_ids

def load_and_filter_manifest(manifest_local_path, eur_person_ids, n_samples=20):
    """Load manifest CSV and filter to EUR samples."""
    print(f"Loading manifest from {manifest_local_path}...")
    manifest_df = pl.read_csv(manifest_local_path)

    print(f"Manifest total rows: {len(manifest_df)}")
    print(f"Manifest columns: {manifest_df.columns}")

    # Filter to EUR samples
    eur_manifest = manifest_df.filter(pl.col('person_id').is_in(eur_person_ids))
    print(f"EUR samples in manifest: {len(eur_manifest)}")

    # Select first n_samples
    selected = eur_manifest.head(n_samples)
    print(f"Selected {len(selected)} samples for testing")

    return selected

def extract_cram_files(manifest_df):
    """Extract CRAM/CRAI file paths from manifest dataframe."""
    cram_files = []
    for row in manifest_df.iter_rows(named=True):
        cram_uri = row['cram_uri']
        crai_uri = row['cram_index_uri']
        cram_files.append((cram_uri, crai_uri))
    return cram_files

# In[ ]:

# ============================================================================
# BASH SCRIPT GENERATION: Serial Execution with Timing
# ============================================================================

def generate_instrumented_script(num_samples):
    """Generate bash script for serial execution with timing instrumentation."""

    script = f'''#!/bin/bash
set -e

NUM_SAMPLES=${{NUM_SAMPLES:-{num_samples}}}

# Initialize timing CSV
TIMING_CSV="/tmp/timing_results.csv"
echo "step,start_epoch,end_epoch,duration_sec" > $TIMING_CSV

echo "=== Starting GATK conversions (sequential) ==="
GATK_START=$(date +%s)

for idx in $(seq 1 $NUM_SAMPLES); do
    cram_var="INPUT_CRAM_$idx"
    crai_var="INPUT_CRAI_$idx"

    start=$(date +%s)
    echo "  [GATK] Processing sample $idx..."
    gatk PrintReads \\
        -R gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \\
        -I "${{!cram_var}}" \\
        --read-index "${{!crai_var}}" \\
        -L chr19:54000000-55100000 \\
        --gcs-project-for-requester-pays $GOOGLE_PROJECT \\
        -O "sample${{idx}}_chr19.bam" 2>&1

    end=$(date +%s)
    echo "  [GATK] sample${{idx}}: $((end - start))s"
done

GATK_END=$(date +%s)
GATK_DURATION=$((GATK_END - GATK_START))
echo "gatk_total,$GATK_START,$GATK_END,$GATK_DURATION" >> $TIMING_CSV
echo "GATK conversions completed in $GATK_DURATION seconds"

# Step 2: Sequential kir-mapper map
echo "=== Starting kir-mapper map (sequential, using all cores) ==="
KIRMAP_START=$(date +%s)

for idx in $(seq 1 $NUM_SAMPLES); do
    start=$(date +%s)
    echo "  [KIR-MAP] Processing sample $idx..."
    kir-mapper map \\
        -bam "sample${{idx}}_chr19.bam" \\
        -sample "sample${{idx}}" \\
        -output ./kir_output \\
        -threads 8 2>&1

    end=$(date +%s)
    echo "  [KIR-MAP] sample${{idx}}: $((end - start))s"
done

KIRMAP_END=$(date +%s)
KIRMAP_DURATION=$((KIRMAP_END - KIRMAP_START))
echo "kirmap_total,$KIRMAP_START,$KIRMAP_END,$KIRMAP_DURATION" >> $TIMING_CSV
echo "kir-mapper map completed in $KIRMAP_DURATION seconds"

# Step 3: ncopy (batch operation)
echo "=== Running ncopy ==="
NCOPY_START=$(date +%s)
kir-mapper ncopy -output ./kir_output -threads 16 2>&1
NCOPY_END=$(date +%s)
NCOPY_DURATION=$((NCOPY_END - NCOPY_START))
echo "ncopy,$NCOPY_START,$NCOPY_END,$NCOPY_DURATION" >> $TIMING_CSV
echo "ncopy completed in $NCOPY_DURATION seconds"

# Step 4: genotype (batch operation)
echo "=== Running genotype ==="
GENO_START=$(date +%s)
kir-mapper genotype -output ./kir_output -threads 16 2>&1
GENO_END=$(date +%s)
GENO_DURATION=$((GENO_END - GENO_START))
echo "genotype,$GENO_START,$GENO_END,$GENO_DURATION" >> $TIMING_CSV
echo "genotype completed in $GENO_DURATION seconds"

# Copy results
echo "=== Copying results to output directory ==="
mkdir -p $OUTPUT_DIR

# Copy timing CSV first (most important for analysis)
if [ -f "$TIMING_CSV" ]; then
    cp $TIMING_CSV $OUTPUT_DIR/timing_results.csv
    echo "✓ Copied timing results CSV"
fi

cp -r kir_output/* $OUTPUT_DIR/ 2>/dev/null || echo "kir_output directory empty or not found"

# Calculate totals (now includes genotype)
TOTAL_TIME=$((GENO_END - GATK_START))
PER_SAMPLE_AVG=$((TOTAL_TIME / NUM_SAMPLES))

# Write summary row
echo "total,$GATK_START,$GENO_END,$TOTAL_TIME" >> $TIMING_CSV

echo "=== Test Complete ==="
echo "Total time: $TOTAL_TIME seconds"
echo "Per-sample average: $PER_SAMPLE_AVG seconds/sample"
'''

    return script

# In[ ]:

# ============================================================================
# DSUB JOB CONFIGURATION
# ============================================================================

def build_env_dict(config, cram_files):
    """Build environment variables dictionary for dsub job."""
    num_samples, machine_type = config

    env_dict = {
        'GOOGLE_PROJECT': GOOGLE_PROJECT,
        'NUM_SAMPLES': str(num_samples),
    }

    # Add CRAM files (cycle through available CRAMs if num_samples > len(cram_files))
    for i in range(1, num_samples + 1):
        cram_idx = (i - 1) % len(cram_files)
        cram_uri, _ = cram_files[cram_idx]
        env_dict[f'INPUT_CRAM_{i}'] = cram_uri

    return env_dict

def build_input_dict(config, cram_files):
    """Build input files dictionary for dsub job (CRAI index files)."""
    num_samples, machine_type = config

    input_dict = {}

    # Add CRAI index files
    for i in range(1, num_samples + 1):
        crai_idx = (i - 1) % len(cram_files)
        _, crai_uri = cram_files[crai_idx]
        input_dict[f'INPUT_CRAI_{i}'] = crai_uri

    return input_dict

# In[ ]:

# ============================================================================
# MAIN EXECUTION: Job Submission
# ============================================================================

def submit_test_job(config, cram_files, output_base, dry_run=False):
    """Submit a dsub job for given configuration."""
    num_samples, machine_type = config

    # Generate script
    script_content = generate_instrumented_script(num_samples)
    script_name = f"kir_serial_{num_samples}s_{machine_type}.sh"

    # Write script to file
    with open(script_name, 'w') as f:
        f.write(script_content)

    # Build env and input dicts
    env_dict = build_env_dict(config, cram_files)
    input_dict = build_input_dict(config, cram_files)

    job_name = f"kir_serial_{machine_type}"
    output_dir = f"{output_base}/{num_samples}s_{machine_type}"

    print(f"\n{'='*70}")
    print(f"Configuration: {num_samples} samples (serial) on {machine_type}")
    print(f"Job name: {job_name}")
    print(f"Output: {output_dir}")
    print(f"{'='*70}")

    # Build custom args with optional --dry-run flag
    custom_args = f"--output-recursive OUTPUT_DIR={output_dir}/"
    if dry_run:
        custom_args = f"--dry-run {custom_args}"

    # Create dsub job
    job = at.Dsub(
        provider="google-batch",
        machine_type=machine_type,
        docker_image=DOCKER_IMAGE,
        job_script_name=script_name,
        job_name=job_name,
        env_dict=env_dict,
        input_dict=input_dict,
        output_dict={},
        custom_args=custom_args
    )

    if not dry_run:
        print("Submitting job...")
        job.run(show_command=False)
        print(f"✓ Job submitted successfully")
    else:
        print("DRY RUN: Previewing job configuration (not submitting)")
        job.run(show_command=True)

    return job

# In[ ]:

# Download manifest from GCS
print("Downloading manifest from GCS...")
get_ipython().system('gsutil -u $GOOGLE_PROJECT cp gs://fc-aou-datasets-controlled/v8/wgs/cram/manifest.csv .')
get_ipython().system('head -n 3 manifest.csv')

# In[ ]:

# Load ancestry metadata and filter for EUR samples
ancestry_df = load_ancestry_metadata(WORKSPACE_BUCKET)
eur_person_ids = get_eur_person_ids(ancestry_df)

# In[ ]:

# Load and filter manifest to EUR samples
eur_manifest = load_and_filter_manifest('manifest.csv', eur_person_ids, n_samples=20)

print("\nSelected EUR samples for testing:")
print(eur_manifest[['person_id', 'cram_uri']].head(10))

# In[ ]:

# Extract CRAM file list
cram_files = extract_cram_files(eur_manifest)
print(f"\nTotal CRAM files selected: {len(cram_files)}")
print("\nFirst 3 CRAM file paths:")
for i, (cram, crai) in enumerate(cram_files[:3], 1):
    print(f"  {i}. {cram}")

# In[ ]:

# ## Submit Test Jobs

# In[ ]:

# Submit jobs for each configuration
print("\n" + "="*70)
print("SUBMITTING SERIAL EXECUTION TEST JOBS")
print("="*70)

submitted_jobs = []

for config in TEST_CONFIGURATIONS:
    num_samples, machine_type = config

    try:
        job = submit_test_job(config, cram_files, OUTPUT_BASE, dry_run=False)
        submitted_jobs.append({
            'config': config,
            'job': job,
            'timestamp': datetime.now().isoformat()
        })
    except Exception as e:
        print(f"ERROR submitting job for {config}: {e}")

print(f"\n✓ Submitted {len(submitted_jobs)} jobs")

# In[ ]:

# ## Monitor Job Status

# In[ ]:

# Check status of all submitted jobs
print("\n" + "="*70)
print("JOB STATUS")
print("="*70)

for job_info in submitted_jobs:
    config = job_info['config']
    job = job_info['job']
    num_samples, machine_type = config

    print(f"\nConfiguration: {machine_type}")
    try:
        job.check_status(streaming=False)
    except Exception as e:
        print(f"Could not retrieve status: {e}")

# In[ ]:

# ## Retrieve and Parse Timing Results

print("\n" + "="*70)
print("TIMING RESULTS")
print("="*70)

def retrieve_timing_csv(bucket, num_samples, machine_type):
    """Download timing CSV from GCS."""
    storage_client = storage.Client(project=GOOGLE_PROJECT)
    bucket_obj = storage_client.bucket(bucket.replace('gs://', ''))

    csv_path = f"dsub/results/kir_serial/{num_samples}s_{machine_type}/timing_results.csv"
    blob = bucket_obj.blob(csv_path)

    if not blob.exists():
        return None

    try:
        csv_content = blob.download_as_string().decode('utf-8')
        from io import StringIO
        df = pl.read_csv(StringIO(csv_content))
        return df
    except Exception as e:
        print(f"Error parsing CSV from {csv_path}: {e}")
        return None

# Collect timing results from all completed jobs
all_results = []

for config in TEST_CONFIGURATIONS:
    num_samples, machine_type = config

    print(f"\nRetrieving timing for {machine_type}...")
    timing_df = retrieve_timing_csv(WORKSPACE_BUCKET, num_samples, machine_type)

    if timing_df is not None:
        # Extract total time
        total_row = timing_df.filter(pl.col('step') == 'total')
        if len(total_row) > 0:
            total_time = total_row['duration_sec'][0]
            per_sample_time = total_time / num_samples

            all_results.append({
                'machine_type': machine_type,
                'total_time_sec': total_time,
                'per_sample_sec': per_sample_time,
                'per_sample_min': per_sample_time / 60
            })
            print(f"  ✓ Total: {total_time}s ({per_sample_time:.1f}s/sample, {per_sample_time/60:.2f}m/sample)")
        else:
            print(f"  ✗ No total row found in CSV")
    else:
        print(f"  ✗ Could not retrieve timing CSV (job may still be running)")

# Display comparison table
if all_results:
    print("\n" + "="*70)
    print("MACHINE TYPE COMPARISON")
    print("="*70)

    results_df = pl.DataFrame(all_results)

    # Sort by per-sample time
    results_df = results_df.sort('per_sample_sec')

    # Display table
    print(f"\n{'Machine Type':<20} {'Total Time':<15} {'Per-Sample (s)':<15} {'Per-Sample (m)':<15}")
    print("-" * 65)
    for row in results_df.iter_rows(named=True):
        print(f"{row['machine_type']:<20} {row['total_time_sec']:>6.0f}s       {row['per_sample_sec']:>6.1f}s        {row['per_sample_min']:>6.2f}m")

    # Identify best option
    best = results_df.row(0, named=True)
    baseline = 6.4 * 60  # 6.4 min in seconds (colleague's baseline)
    speedup = baseline / best['per_sample_sec']

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\n✓ Best configuration: {best['machine_type']}")
    print(f"  - Per-sample time: {best['per_sample_sec']:.1f}s ({best['per_sample_min']:.2f}m)")
    print(f"  - Speedup vs colleague baseline (6.4m/sample): {speedup:.2f}x")
    print(f"  - Total time for 20 samples: {best['total_time_sec']/60:.1f}m")
else:
    print("\n✗ No timing results retrieved yet")
    print("Jobs may still be running. Check back in a few hours.")

print(f"\n✓ Test complete!")
print(f"Full results available at: {OUTPUT_BASE}")
