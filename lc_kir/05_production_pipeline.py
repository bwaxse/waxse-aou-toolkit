#!/usr/bin/env python
# coding: utf-8

"""
KIR Copy Number Production: Pipeline Submission
==============================================

Submits dsub jobs for production KIR copy number calling on full dataset.
- Batches: 200 samples per job
- Machine: n2-standard-16 (chosen from test results)
- Parallelization: 10 samples per job
- Thresholds: Pre-calibrated ancestry-specific files
- Preemptible: Yes (70% cost savings)

Runs independently for each batch and ancestry.
Idempotent: Check GCS before submitting; skip completed batches.
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
MACHINE_TYPE = 'n2-standard-16'  # Chosen based on test results
PARALLEL_SAMPLES = 10
THREADS_PER_SAMPLE = 2
DOCKER_IMAGE = "phetk/gatk-kirmapper:0.1"
BASH_SCRIPT = "bash_scripts/kir_ncopy_production.sh"
BATCH_MANIFEST_BASE = f"{WORKSPACE_BUCKET}/kir/production/batch_manifests"
THRESHOLDS_BASE = f"{WORKSPACE_BUCKET}/kir/thresholds"
OUTPUT_BASE = f"{WORKSPACE_BUCKET}/kir/production"

print(f"\n{'='*80}")
print(f"KIR Production: Pipeline Submission")
print(f"{'='*80}")
print(f"Machine type: {MACHINE_TYPE}")
print(f"Parallelization: {PARALLEL_SAMPLES} samples")
print(f"Docker image: {DOCKER_IMAGE}")
print(f"Batch manifests: {BATCH_MANIFEST_BASE}/")
print(f"Thresholds: {THRESHOLDS_BASE}/")

# ============================================================================
# Helper Functions
# ============================================================================

def check_batch_complete(bucket, ancestry, batch_id):
    """Check if batch is already complete in GCS."""
    output_path = f"kir/production/{ancestry}/{batch_id}"
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

def load_batch_manifest(bucket, ancestry, batch_id):
    """Load batch manifest from GCS."""
    manifest_path = f"kir/production/batch_manifests/{ancestry}/{batch_id}.csv"

    storage_client = storage.Client(project=GOOGLE_PROJECT)
    bucket_obj = storage_client.bucket(bucket.replace('gs://', ''))
    blob = bucket_obj.blob(manifest_path)

    if not blob.exists():
        return None

    csv_content = blob.download_as_string().decode('utf-8')
    manifest_df = pl.read_csv(csv_content.encode())
    return manifest_df

def get_thresholds_path(bucket, ancestry):
    """Get path to pre-calibrated thresholds for ancestry."""
    # Check if thresholds exist
    storage_client = storage.Client(project=GOOGLE_PROJECT)
    bucket_obj = storage_client.bucket(bucket.replace('gs://', ''))
    blob_path = f"kir/thresholds/thresholds_{ancestry}.txt"
    blob = bucket_obj.blob(blob_path)

    if not blob.exists():
        return None

    return f"gs://{bucket}/{blob_path}"

def build_env_dict(manifest_df, ancestry, batch_id, thresholds_path):
    """Build environment variables from batch manifest."""
    env_dict = {
        'NUM_SAMPLES': str(len(manifest_df)),
        'PARALLEL_SAMPLES': str(PARALLEL_SAMPLES),
        'THREADS_PER_SAMPLE': str(THREADS_PER_SAMPLE),
        'ANCESTRY': ancestry,
        'BATCH_ID': batch_id,
        'INPUT_THRESHOLDS': thresholds_path if thresholds_path else '',
        'GOOGLE_PROJECT': GOOGLE_PROJECT,
    }

    # Add CRAM/CRAI paths as numbered environment variables
    for i, row in enumerate(manifest_df.iter_rows(named=True), start=1):
        env_dict[f'INPUT_CRAM_{i}'] = row['cram_uri']
        env_dict[f'INPUT_CRAI_{i}'] = row['cram_index_uri']

    return env_dict

def submit_batch_job(manifest_df, ancestry, batch_id, thresholds_path):
    """Submit dsub production job for batch."""
    # Check if already complete
    if check_batch_complete(WORKSPACE_BUCKET, ancestry, batch_id):
        return None

    # Build environment variables
    env_dict = build_env_dict(manifest_df, ancestry, batch_id, thresholds_path)

    # Output directory
    output_path = f"{OUTPUT_BASE}/{ancestry}/{batch_id}"

    # Submit job
    job = at.Dsub(
        provider="google-batch",
        machine_type=MACHINE_TYPE,
        docker_image=DOCKER_IMAGE,
        job_script_name=BASH_SCRIPT,
        job_name=f"kir_prod_{ancestry}_{batch_id}",
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
print("Verifying Pre-Calibrated Thresholds")
print(f"{'='*80}")

thresholds_available = {}
for ancestry in ANCESTRIES:
    thresholds_path = get_thresholds_path(WORKSPACE_BUCKET, ancestry)
    if thresholds_path:
        print(f"✓ {ancestry}: {thresholds_path}")
        thresholds_available[ancestry] = thresholds_path
    else:
        print(f"✗ {ancestry}: Thresholds not found - run calibration first!")

if len(thresholds_available) < len(ANCESTRIES):
    print(f"\n⚠ Not all thresholds available - cannot proceed")
    print(f"Run 03_extract_thresholds.py first to extract calibration results")
    import sys
    sys.exit(1)

print(f"\n✓ All thresholds available")

# ============================================================================
# List and Submit Batches
# ============================================================================

print(f"\n{'='*80}")
print("Submitting Production Batches")
print(f"{'='*80}")

storage_client = storage.Client(project=GOOGLE_PROJECT)
bucket_obj = storage_client.bucket(WORKSPACE_BUCKET.replace('gs://', ''))

submitted_jobs = []
skipped_batches = []
failed_batches = []

for ancestry in ANCESTRIES:
    print(f"\n{ancestry}:")
    thresholds_path = thresholds_available[ancestry]

    # List all batch manifests for this ancestry
    prefix = f"kir/production/batch_manifests/{ancestry}/"
    blobs = storage_client.list_blobs(bucket_obj, prefix=prefix)
    batch_files = [b.name.split('/')[-1].replace('.csv', '') for b in blobs if b.name.endswith('.csv')]

    if not batch_files:
        print(f"  ⚠ No batch manifests found")
        continue

    batch_files.sort()
    print(f"  Found {len(batch_files)} batches")

    for batch_id in batch_files:
        # Load manifest
        manifest_df = load_batch_manifest(WORKSPACE_BUCKET, ancestry, batch_id)
        if manifest_df is None:
            failed_batches.append(f"{ancestry}/{batch_id}")
            continue

        # Submit job
        try:
            job = submit_batch_job(manifest_df, ancestry, batch_id, thresholds_path)
            if job:
                submitted_jobs.append({
                    'ancestry': ancestry,
                    'batch_id': batch_id,
                    'num_samples': len(manifest_df)
                })
                print(f"  → {batch_id}: {len(manifest_df)} samples submitted")
            else:
                skipped_batches.append(f"{ancestry}/{batch_id}")
                print(f"  ✓ {batch_id}: Already complete (skipped)")
        except Exception as e:
            failed_batches.append(f"{ancestry}/{batch_id}")
            print(f"  ✗ {batch_id}: Error - {e}")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*80}")
print("Summary")
print(f"{'='*80}")

total_samples = sum(j['num_samples'] for j in submitted_jobs)

print(f"\nSubmitted: {len(submitted_jobs)} jobs ({total_samples:,} samples)")
if submitted_jobs:
    for job in submitted_jobs[:5]:
        print(f"  - {job['ancestry']}/{job['batch_id']}: {job['num_samples']} samples")
    if len(submitted_jobs) > 5:
        print(f"  ... and {len(submitted_jobs) - 5} more")

if skipped_batches:
    print(f"\nSkipped (already complete): {len(skipped_batches)} batches")
    for batch in skipped_batches[:3]:
        print(f"  - {batch}")
    if len(skipped_batches) > 3:
        print(f"  ... and {len(skipped_batches) - 3} more")

if failed_batches:
    print(f"\nFailed: {len(failed_batches)} batches")
    for batch in failed_batches:
        print(f"  - {batch}")

# ============================================================================
# Cost Estimation
# ============================================================================

print(f"\n{'='*80}")
print("Cost Estimation (Preemptible n2-standard-16)")
print(f"{'='*80}")

cost_per_sample = 0.0045  # From test results
estimated_cost = total_samples * cost_per_sample
estimated_hours = (total_samples / 10) * (69.5 / 3600)  # 10 parallel, 69.5s/sample

print(f"  Submitted samples: {total_samples:,}")
print(f"  Cost per sample: ${cost_per_sample:.4f}")
print(f"  Estimated cost: ${estimated_cost:,.2f}")
print(f"  Estimated time: {estimated_hours:,.1f} hours (if all jobs run in parallel)")

# ============================================================================
# Next Steps
# ============================================================================

print(f"\n{'='*80}")
print("Next Steps")
print(f"{'='*80}")

print(f"\n1. Monitor job progress:")
print(f"   - Verily Workbench: Monitor > Batch Submission Jobs")
print(f"   - Or: gsutil ls {WORKSPACE_BUCKET}/dsub/logs/ | head -20")

print(f"\n2. Once jobs complete (typically 3-5 days):")
print(f"   - Run 06_collect_ncopy_results.py to aggregate results")
print(f"   - Combine copy_numbers.table.txt and presence.table.txt per ancestry")
print(f"\n3. If jobs fail:")
print(f"   - Check logs: gsutil cat {WORKSPACE_BUCKET}/dsub/logs/{{JOB_ID}}.log")
print(f"   - Re-run this script (idempotent - only resubmits missing batches)")
