#!/usr/bin/env python
# coding: utf-8

# # KIR-Mapper Parallelization Test Script
# Test multiple dsub configurations to determine optimal parallelization strategy

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

print(f"Bucket: {WORKSPACE_BUCKET}")
print(f"CDR: {WORKSPACE_CDR}")
print(f"Project: {GOOGLE_PROJECT}")

# In[ ]:

# ============================================================================
# CONFIGURATION: Test Matrix for Parallelization Sweep
# ============================================================================

TEST_CONFIGURATIONS = [
    # Format: (num_samples, parallel_samples, threads_per_sample, machine_type)

    # Baseline comparisons
    (4, 2, 2, 'n2-standard-4'),       # Small: 4 samples, 2 parallel, baseline
    (4, 4, 2, 'n2-standard-8'),       # Small: 4 samples, full parallel

    # Medium scale with escalating parallelization
    (20, 5, 2, 'n2-standard-8'),      # Medium: 20 samples, conservative parallel
    (20, 10, 2, 'n2-standard-16'),    # Medium: 20 samples, moderate parallel (current winner)
    (20, 20, 1, 'n2-standard-16'),    # Medium: 20 samples, high parallel, low threads (cost optimization test)
    (20, 20, 2, 'n2-standard-32'),    # Medium: 20 samples, full parallel

    # Production scale with high parallelization
    (64, 32, 2, 'n2-standard-64'),    # Large: 64 samples, half parallel
    (64, 64, 2, 'n2-standard-64'),    # Large: 64 samples, full parallel (optimal throughput)
]

DOCKER_IMAGE = "phetk/gatk-kirmapper:0.1"
OUTPUT_BASE = f"{WORKSPACE_BUCKET}/dsub/results/kir_parallel"

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
# BASH SCRIPT GENERATION: Instrumented Test Script with Timing
# ============================================================================

def generate_instrumented_script(num_samples, parallel_samples, threads_per_sample):
    """Generate bash script with parallelization and timing instrumentation."""

    script = f'''#!/bin/bash
set -e

# Configuration (passed as env vars)
THREADS_PER_SAMPLE=${{THREADS_PER_SAMPLE:-{threads_per_sample}}}
PARALLEL_SAMPLES=${{PARALLEL_SAMPLES:-{parallel_samples}}}
NUM_SAMPLES=${{NUM_SAMPLES:-{num_samples}}}

# Helper function to run GATK with environment variable lookup
run_gatk_with_env() {{
    local idx=$1
    local cram_var="INPUT_CRAM_$idx"
    local crai_var="INPUT_CRAI_$idx"
    run_gatk "$idx" "${{!cram_var}}" "${{!crai_var}}"
}}

echo "=== Configuration ==="
echo "NUM_SAMPLES: $NUM_SAMPLES"
echo "PARALLEL_SAMPLES: $PARALLEL_SAMPLES"
echo "THREADS_PER_SAMPLE: $THREADS_PER_SAMPLE"

# Function: GATK conversion for one sample
run_gatk() {{
    local idx=$1
    local cram=$2
    local crai=$3
    local start=$(date +%s)

    echo "  [GATK] Processing sample $idx..."
    gatk PrintReads \\
        -R gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \\
        -I "$cram" \\
        --read-index "$crai" \\
        -L chr19:54000000-55100000 \\
        --gcs-project-for-requester-pays $GOOGLE_PROJECT \\
        -O "sample${{idx}}_chr19.bam" 2>&1

    local end=$(date +%s)
    echo "  [GATK] sample${{idx}}: $((end - start))s"
}}

# Function: kir-mapper map for one sample
run_kirmap() {{
    local idx=$1
    local start=$(date +%s)

    echo "  [KIR-MAP] Processing sample $idx..."
    kir-mapper map \\
        -bam "sample${{idx}}_chr19.bam" \\
        -sample "sample${{idx}}" \\
        -output ./kir_output \\
        -threads $THREADS_PER_SAMPLE 2>&1

    local end=$(date +%s)
    echo "  [KIR-MAP] sample${{idx}}: $((end - start))s"
}}

export -f run_gatk run_gatk_with_env run_kirmap
export THREADS_PER_SAMPLE GOOGLE_PROJECT

# Step 1: Parallel GATK conversions
echo "=== Starting GATK conversions (parallel=$PARALLEL_SAMPLES) ==="
GATK_START=$(date +%s)

# Check if GNU parallel is available
if command -v parallel &> /dev/null; then
    echo "Using GNU parallel for GATK conversions..."
    seq 1 $NUM_SAMPLES | \\
        parallel -j $PARALLEL_SAMPLES 'run_gatk_with_env {{}}'
else
    echo "GNU parallel not found, using sequential processing..."
    for idx in $(seq 1 $NUM_SAMPLES); do
        run_gatk_with_env "$idx"
    done
fi

GATK_END=$(date +%s)
echo "GATK conversions completed in $((GATK_END - GATK_START)) seconds"

# Step 2: Parallel kir-mapper map
echo "=== Starting kir-mapper map (parallel=$PARALLEL_SAMPLES) ==="
KIRMAP_START=$(date +%s)

if command -v parallel &> /dev/null; then
    echo "Using GNU parallel for kir-mapper map..."
    seq 1 $NUM_SAMPLES | \\
        parallel -j $PARALLEL_SAMPLES 'run_kirmap {{}}'
else
    echo "GNU parallel not found, using sequential processing..."
    for idx in $(seq 1 $NUM_SAMPLES); do
        run_kirmap "$idx"
    done
fi

KIRMAP_END=$(date +%s)
echo "kir-mapper map completed in $((KIRMAP_END - KIRMAP_START)) seconds"

# Step 3: ncopy (batch operation)
echo "=== Running ncopy ==="
NCOPY_START=$(date +%s)
kir-mapper ncopy -output ./kir_output -threads $((PARALLEL_SAMPLES * THREADS_PER_SAMPLE))
NCOPY_END=$(date +%s)
echo "ncopy completed in $((NCOPY_END - NCOPY_START)) seconds"

# Step 4: genotype (batch operation)
echo "=== Running genotype ==="
GENO_START=$(date +%s)
kir-mapper genotype -output ./kir_output -threads $((PARALLEL_SAMPLES * THREADS_PER_SAMPLE))
GENO_END=$(date +%s)
echo "genotype completed in $((GENO_END - GENO_START)) seconds"

# Step 5: Haplotype (only if NUM_SAMPLES >= 50)
if [ $NUM_SAMPLES -ge 50 ]; then
    echo "=== Running haplotype ==="
    HAPLO_START=$(date +%s)
    kir-mapper haplotype -output ./kir_output -threads $((PARALLEL_SAMPLES * THREADS_PER_SAMPLE))
    HAPLO_END=$(date +%s)
    echo "haplotype completed in $((HAPLO_END - HAPLO_START)) seconds"
else
    echo "=== Skipping haplotype (need >=50 samples, have $NUM_SAMPLES) ==="
fi

echo "=== Output statistics ==="
du -sh ./kir_output/ 2>/dev/null || echo "kir_output not yet created"
find ./kir_output -type f -name "*.txt" -o -name "*.tsv" -o -name "*.vcf*" 2>/dev/null | wc -l || true

# Copy results
echo "=== Copying results to output directory ==="
mkdir -p $OUTPUT_DIR
cp -r kir_output/* $OUTPUT_DIR/ 2>/dev/null || echo "kir_output directory empty or not found"

echo "=== Test Complete ==="
'''

    return script

# In[ ]:

# ============================================================================
# DSUB JOB CONFIGURATION
# ============================================================================

def build_env_dict(config, cram_files):
    """Build environment variables dictionary for dsub job."""
    num_samples, parallel_samples, threads_per_sample, machine_type = config

    # Dynamically build INPUT_CRAM_* and INPUT_CRAI_* env vars
    env_dict = {
        'GOOGLE_PROJECT': GOOGLE_PROJECT,
        'NUM_SAMPLES': str(num_samples),
        'PARALLEL_SAMPLES': str(parallel_samples),
        'THREADS_PER_SAMPLE': str(threads_per_sample),
    }

    # Add CRAM files (cycle through available CRAMs if num_samples > len(cram_files))
    for i in range(1, num_samples + 1):
        cram_idx = (i - 1) % len(cram_files)
        cram_uri, _ = cram_files[cram_idx]
        env_dict[f'INPUT_CRAM_{i}'] = cram_uri

    return env_dict

def build_input_dict(config, cram_files):
    """Build input files dictionary for dsub job (CRAI index files)."""
    num_samples, _, _, _ = config

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
    num_samples, parallel_samples, threads_per_sample, machine_type = config

    # Generate script
    script_content = generate_instrumented_script(num_samples, parallel_samples, threads_per_sample)
    script_name = f"kir_parallel_{num_samples}s_{parallel_samples}p_{threads_per_sample}t.sh"

    # Write script to file
    with open(script_name, 'w') as f:
        f.write(script_content)

    # Build env and input dicts
    env_dict = build_env_dict(config, cram_files)
    input_dict = build_input_dict(config, cram_files)

    job_name = f"kir_parallel_n{num_samples}_p{parallel_samples}"
    output_dir = f"{output_base}/n{num_samples}_p{parallel_samples}_{machine_type}"

    print(f"\n{'='*70}")
    print(f"Configuration: {num_samples} samples, {parallel_samples} parallel, {machine_type}")
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

# ## Verify Docker Image (Optional - Quick Check)
# This submits a quick test job to verify the Docker image has required tools

# In[ ]:

verify_docker = False  # Set to True to run verification

if verify_docker:
    test_script = """#!/bin/bash
echo "=== Checking Docker Image Contents ==="
which parallel && echo "✓ GNU parallel found" || echo "✗ GNU parallel NOT found"
parallel --version 2>/dev/null | head -1 || true
echo ""
echo "=== Tool versions ==="
gatk --version 2>/dev/null | head -1 || echo "GATK version check failed"
kir-mapper --help 2>/dev/null | head -1 || echo "kir-mapper check failed"
which samtools bcftools || echo "samtools/bcftools missing"
echo ""
echo "=== Available commands ==="
ls /usr/local/bin/ | head -20
"""

    # Write script to file
    script_name = "verify_docker.sh"
    with open(script_name, 'w') as f:
        f.write(test_script)

    # Submit quick verification job
    verify_job = at.Dsub(
        provider="google-batch",
        machine_type="n2-standard-2",
        docker_image=DOCKER_IMAGE,
        job_script_name=script_name,
        job_name="verify_docker_tools",
        env_dict={"GOOGLE_PROJECT": GOOGLE_PROJECT},
        input_dict={},
        output_dict={},
        custom_args=f"--output-recursive OUTPUT_DIR={WORKSPACE_BUCKET}/dsub/verify_docker/"
    )

    print("Submitting quick Docker verification job...")
    verify_job.run(show_command=False)
    print("✓ Verification job submitted")
    print("\nCheck results in ~5-10 minutes at:")
    print(f"  gsutil cat {WORKSPACE_BUCKET}/dsub/verify_docker/stdout")
    print("")
    print("What to look for:")
    print("  ✓ GNU parallel found → Ready for parallelization")
    print("  ✗ GNU parallel NOT found → Will use sequential fallback (still works)")

# In[ ]:

# ## Download and Filter Manifest

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
# Note: Need 64 samples for largest test configuration
eur_manifest = load_and_filter_manifest('manifest.csv', eur_person_ids, n_samples=64)

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
print("SUBMITTING PARALLELIZATION TEST JOBS")
print("="*70)

submitted_jobs = []

for config in TEST_CONFIGURATIONS:
    num_samples, parallel_samples, threads_per_sample, machine_type = config

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
    num_samples, parallel_samples, threads_per_sample, machine_type = config

    print(f"\nConfiguration: {num_samples}s_{parallel_samples}p")
    try:
        job.check_status(streaming=False)
    except Exception as e:
        print(f"Could not retrieve status: {e}")

# In[ ]:

# ## Retrieve and Parse Results

# In[ ]:

def parse_timing_log(gcs_path):
    """Download and parse timing log from GCS."""
    local_file = 'timing_log.csv'
    try:
        get_ipython().system(f'gsutil cp {gcs_path} {local_file}')
        df = pl.read_csv(local_file)
        return df
    except Exception as e:
        print(f"Could not retrieve timing log from {gcs_path}: {e}")
        return None

# In[ ]:

# Timing information is logged to stdout by the bash scripts
print("\n" + "="*70)
print("TIMING RESULTS")
print("="*70)

print("\n✓ Timing information captured in job stdout logs")
print("\nTo retrieve timing for a completed job:")
print("  gsutil cat {WORKSPACE_BUCKET}/dsub/logs/{JOB_ID}.log | grep 'completed in'")
print("\nExample timings will show:")
print("  - GATK conversions completed in X seconds")
print("  - kir-mapper map completed in Y seconds")
print("  - ncopy completed in Z seconds")
print("  - genotype completed in A seconds")
print("\nPer-job timing calculation:")
print("  Total time ≈ GATK time + Map time + ncopy time + genotype time")

# In[ ]:

print("\n✓ Test job submission and monitoring complete!")
print(f"Results available at: {OUTPUT_BASE}")
