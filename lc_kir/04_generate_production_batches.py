#!/usr/bin/env python
# coding: utf-8

"""
KIR Copy Number Production: Generate Batch Manifests
===================================================

Generates batch manifests (200 samples per batch) for all WGS participants,
stratified by ancestry. Uses pre-calibrated thresholds.txt from calibration phase.

Output: CSV manifests with person_id, ancestry, cram_uri, crai_uri for each batch
"""

get_ipython().system('pip install tctk polars google-cloud-storage --upgrade -q')

import os
import polars as pl
from google.cloud import bigquery, storage
import math

# Environment setup
WORKSPACE_BUCKET = os.getenv('WORKSPACE_BUCKET')
WORKSPACE_CDR = os.getenv('WORKSPACE_CDR')
GOOGLE_PROJECT = os.getenv('GOOGLE_PROJECT')

print(f"Bucket: {WORKSPACE_BUCKET}")
print(f"CDR: {WORKSPACE_CDR}")
print(f"Project: {GOOGLE_PROJECT}")

# ============================================================================
# Configuration
# ============================================================================

ANCESTRIES = ['EUR', 'AFR', 'AMR', 'EAS', 'SAS']
BATCH_SIZE = 200
MANIFEST_URL = "gs://fc-aou-datasets-controlled/v8/wgs/cram/manifest.csv"
OUTPUT_BASE = f"{WORKSPACE_BUCKET}/kir/production/batch_manifests"

print(f"\n{'='*80}")
print(f"KIR Production: Generate Batch Manifests")
print(f"{'='*80}")
print(f"Ancestries: {', '.join(ANCESTRIES)}")
print(f"Batch size: {BATCH_SIZE} samples")
print(f"Output base: {OUTPUT_BASE}")

# ============================================================================
# Load Data
# ============================================================================

print(f"\n{'='*80}")
print("Loading Data")
print(f"{'='*80}")

# Load ancestry metadata
ancestry_path = f'{WORKSPACE_BUCKET}/data/ancestry_metadata.tsv'
print(f"Loading ancestry metadata from {ancestry_path}...")
ancestry_df = pl.read_csv(ancestry_path, separator='\t')
print(f"✓ Loaded {len(ancestry_df):,} samples")

# Load manifest
print(f"\nDownloading manifest from {MANIFEST_URL}...")
import subprocess
result = subprocess.run(['gsutil', '-u', GOOGLE_PROJECT, 'cp', MANIFEST_URL, 'manifest_production.csv'],
                       capture_output=True, text=True)
if result.returncode != 0:
    print(f"Error downloading manifest: {result.stderr}")
    raise RuntimeError(f"Failed to download manifest")

manifest_df = pl.read_csv('manifest_production.csv')
print(f"✓ Loaded {len(manifest_df):,} samples from manifest")

# ============================================================================
# Generate Batches per Ancestry
# ============================================================================

print(f"\n{'='*80}")
print("Generating Batch Manifests")
print(f"{'='*80}")

storage_client = storage.Client(project=GOOGLE_PROJECT)
bucket_obj = storage_client.bucket(WORKSPACE_BUCKET.replace('gs://', ''))

batch_summary = {}

for ancestry in ANCESTRIES:
    print(f"\n{ancestry}:")

    # Get person IDs for this ancestry
    ancestry_col = 'ancestry_pred_other'
    ancestry_ids = ancestry_df.filter(
        pl.col(ancestry_col).str.to_uppercase() == ancestry
    )['research_id'].to_list()

    if not ancestry_ids:
        print(f"  ⚠ No samples found")
        continue

    print(f"  Available: {len(ancestry_ids):,}")

    # Filter manifest to this ancestry
    ancestry_manifest = manifest_df.filter(pl.col('person_id').is_in(ancestry_ids))
    print(f"  In manifest: {len(ancestry_manifest):,}")

    if len(ancestry_manifest) == 0:
        print(f"  ⚠ No samples in manifest, skipping")
        continue

    # Split into batches
    num_batches = math.ceil(len(ancestry_manifest) / BATCH_SIZE)
    print(f"  Batches: {num_batches}")

    ancestry_batches = []

    for batch_num in range(num_batches):
        start_idx = batch_num * BATCH_SIZE
        end_idx = min((batch_num + 1) * BATCH_SIZE, len(ancestry_manifest))
        batch_df = ancestry_manifest.slice(start_idx, end_idx - start_idx)

        batch_id = f"batch_{batch_num+1:04d}"
        ancestry_batches.append({
            'batch_id': batch_id,
            'num_samples': len(batch_df),
            'manifest': batch_df
        })

        # Save batch manifest to GCS
        csv_content = batch_df.write_csv()
        blob_path = f"kir/production/batch_manifests/{ancestry}/{batch_id}.csv"
        blob = bucket_obj.blob(blob_path)
        blob.upload_from_string(csv_content)

    batch_summary[ancestry] = {
        'num_batches': num_batches,
        'total_samples': len(ancestry_manifest),
        'batches': ancestry_batches
    }

    print(f"  ✓ Generated {num_batches} batches")
    total_samples = sum(b['num_samples'] for b in ancestry_batches)
    print(f"    Total samples: {total_samples:,}")
    if total_samples != len(ancestry_manifest):
        print(f"    ⚠ Sample count mismatch!")

# ============================================================================
# Summary Report
# ============================================================================

print(f"\n{'='*80}")
print("Summary")
print(f"{'='*80}")

total_samples = 0
total_batches = 0

for ancestry, info in batch_summary.items():
    total_samples += info['total_samples']
    total_batches += info['num_batches']
    print(f"\n{ancestry}:")
    print(f"  Samples: {info['total_samples']:,}")
    print(f"  Batches: {info['num_batches']}")
    print(f"  Output: gs://{WORKSPACE_BUCKET}/kir/production/batch_manifests/{ancestry}/")

print(f"\n{'='*80}")
print(f"Total: {total_samples:,} samples in {total_batches} batches")
print(f"{'='*80}")

# ============================================================================
# Cost Estimation
# ============================================================================

print(f"\nCost Estimation (Preemptible n2-standard-16):")
cost_per_sample = 0.0045  # From test results
estimated_cost = total_samples * cost_per_sample
print(f"  Per sample: ${cost_per_sample:.4f}")
print(f"  Total: ${estimated_cost:,.2f}")

# ============================================================================
# Next Steps
# ============================================================================

print(f"\n{'='*80}")
print("Next Steps")
print(f"{'='*80}")
print(f"\n1. Wait for calibration phase to complete (04_generate_production_batches.py)")
print(f"2. Verify thresholds.txt exist in: gs://{WORKSPACE_BUCKET}/kir/thresholds/")
print(f"3. Run 05_production_pipeline.py to submit production jobs")
print(f"\nBatch manifests are ready at:")
print(f"  gs://{WORKSPACE_BUCKET}/kir/production/batch_manifests/")
