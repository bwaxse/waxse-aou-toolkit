#!/usr/bin/env python
# coding: utf-8

"""
KIR Copy Number Calibration: Pipeline Submission
================================================

Submits dsub jobs for KIR copy number calibration on ~200 samples per ancestry.
- Machine: n2-standard-16 (4GB RAM per vCPU)
- Parallelization: 10 samples
- Preemptible: Yes (70% cost savings)

Runs independently for each ancestry (EUR, AFR, AMR, EAS, SAS).
"""

get_ipython().system('pip install tctk polars google-cloud-storage --upgrade -q')

from tctk import AoUTools as at
import os
import polars as pl
from google.cloud import storage
from pathlib import Path

# Environment setup
WORKSPACE_BUCKET = os.getenv('WORKSPACE_BUCKET')
WORKSPACE_CDR = os.getenv('WORKSPACE_CDR')
GOOGLE_PROJECT = os.getenv('GOOGLE_PROJECT')

print(f"Bucket: {WORKSPACE_BUCKET}")
print(f"Project: {GOOGLE_PROJECT}")

# ============================================================================
# Configuration
# ============================================================================

ANCESTRIES = ['EUR', 'AFR', 'AMR', 'EAS', 'SAS']
MACHINE_TYPE = 'n2-standard-16'  # Chosen based on test results (20s_10p winner)
PARALLEL_SAMPLES = 10
THREADS_PER_SAMPLE = 2
DOCKER_IMAGE = "phetk/gatk-kirmapper:0.1"
BASH_SCRIPT = "bash_scripts/kir_ncopy_calibration.sh"
OUTPUT_BASE = f"{WORKSPACE_BUCKET}/kir/calibration"
MANIFEST_BASE = f"{WORKSPACE_BUCKET}/kir/calibration/sample_manifests"

print(f"\n{'='*80}")
print(f"KIR Calibration: Pipeline Submission")
print(f"{'='*80}")
print(f"Machine type: {MACHINE_TYPE}")
print(f"Parallelization: {PARALLEL_SAMPLES} samples")
print(f"Docker image: {DOCKER_IMAGE}")
print(f"Bash script: {BASH_SCRIPT}")

# ============================================================================
# Helper Functions
# ============================================================================

def check_ancestry_complete(bucket, ancestry):
    """Check if ancestry calibration is already complete in GCS."""
    output_path = f"kir/calibration/{ancestry}"
    required_files = [
        f"{output_path}/ncopy/copy_numbers.table.txt",
        f"{output_path}/ncopy/presence.table.txt",
        f"{output_path}/timing_log.csv"
    ]

    storage_client = storage.Client(project=GOOGLE_PROJECT)
    bucket_obj = storage_client.bucket(bucket.replace('gs://', ''))

    for file_path in required_files:
        blob = bucket_obj.blob(file_path)
        if not blob.exists():
            return False
    return True

def load_ancestry_manifest(bucket, ancestry):
    """Load sample manifest for ancestry from GCS."""
    manifest_path = f"kir/calibration/sample_manifests/{ancestry}_samples.csv"
    print(f"  Loading manifest from gs://{bucket}/{manifest_path}...")

    storage_client = storage.Client(project=GOOGLE_PROJECT)
    bucket_obj = storage_client.bucket(bucket.replace('gs://', ''))
    blob = bucket_obj.blob(manifest_path)

    if not blob.exists():
        raise FileNotFoundError(f"Manifest not found: gs://{bucket}/{manifest_path}")

    csv_content = blob.download_as_string().decode('utf-8')
    manifest_df = pl.read_csv(csv_content.encode())
    return manifest_df

def build_env_dict(manifest_df, ancestry):
    """Build environment variables from sample manifest."""
    env_dict = {
        'NUM_SAMPLES': str(len(manifest_df)),
        'PARALLEL_SAMPLES': str(PARALLEL_SAMPLES),
        'THREADS_PER_SAMPLE': str(THREADS_PER_SAMPLE),
        'ANCESTRY': ancestry,
        'GOOGLE_PROJECT': GOOGLE_PROJECT,
    }

    # Add CRAM/CRAI paths as numbered environment variables
    for i, row in enumerate(manifest_df.iter_rows(named=True), start=1):
        env_dict[f'INPUT_CRAM_{i}'] = row['cram_uri']
        env_dict[f'INPUT_CRAI_{i}'] = row['cram_index_uri']

    return env_dict

def submit_calibration_job(manifest_df, ancestry):
    """Submit dsub calibration job for ancestry."""
    print(f"\n{ancestry}:")
    print(f"  Samples: {len(manifest_df)}")

    # Check if already complete
    if check_ancestry_complete(WORKSPACE_BUCKET, ancestry):
        print(f"  ✓ Already complete, skipping")
        return None

    # Build environment variables
    env_dict = build_env_dict(manifest_df, ancestry)

    # Output directory
    output_path = f"{OUTPUT_BASE}/{ancestry}"

    print(f"  Submitting dsub job...")
    print(f"    Machine: {MACHINE_TYPE}")
    print(f"    Parallelization: {PARALLEL_SAMPLES}")
    print(f"    Output: {output_path}")

    # Submit job
    job = at.Dsub(
        provider="google-batch",
        machine_type=MACHINE_TYPE,
        docker_image=DOCKER_IMAGE,
        job_script_name=BASH_SCRIPT,
        job_name=f"kir_calib_{ancestry}",
        env_dict=env_dict,
        input_dict={},
        output_dict={},
        custom_args=f"--output-recursive OUTPUT_DIR={output_path}",
        preemptible=True  # 70% cost savings
    )

    job.run(show_command=False)
    return job

# ============================================================================
# Main Pipeline
# ============================================================================

print(f"\n{'='*80}")
print("Processing Ancestries")
print(f"{'='*80}")

submitted_jobs = {}

for ancestry in ANCESTRIES:
    print(f"\n--- {ancestry} ---")

    try:
        # Load manifest
        manifest_df = load_ancestry_manifest(WORKSPACE_BUCKET, ancestry)

        # Submit job
        job = submit_calibration_job(manifest_df, ancestry)
        if job:
            submitted_jobs[ancestry] = job

    except FileNotFoundError as e:
        print(f"  ✗ {e}")
    except Exception as e:
        print(f"  ✗ Error: {e}")

# ============================================================================
# Summary and Monitoring
# ============================================================================

print(f"\n{'='*80}")
print("Summary")
print(f"{'='*80}")
print(f"\nSubmitted {len(submitted_jobs)} job(s):")
for ancestry, job in submitted_jobs.items():
    print(f"  - {ancestry}: {job.job_id if hasattr(job, 'job_id') else 'submitted'}")

if submitted_jobs:
    print(f"\nEstimated time per job: ~3.8 hours (200 samples × 69.5s/sample)")
    print(f"Total cost: ~$4.50 (5 jobs × 0.386hr × $0.234/hr preemptible)")

    print(f"\nMonitor jobs:")
    print(f"  1. In Verily Workbench: Monitor > Batch Submission Jobs")
    print(f"  2. Or run: gsutil ls gs://{WORKSPACE_BUCKET}/dsub/logs/")
else:
    print(f"\n✓ All ancestries already complete!")

print(f"\nNext step: Run 03_extract_thresholds.py to extract ancestry-specific thresholds")
