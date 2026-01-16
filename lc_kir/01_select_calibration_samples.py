#!/usr/bin/env python
# coding: utf-8

"""
KIR Copy Number Calibration: Sample Selection
==============================================

Selects 200 random samples per ancestry (EUR, AFR, AMR, EAS, SAS) for threshold calibration.

Outputs sample manifests with person_id, ancestry, cram_uri, crai_uri for each ancestry.
"""

get_ipython().system('pip install tctk polars google-cloud-storage --upgrade -q')

import os
import polars as pl
from google.cloud import storage
import random

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
SAMPLES_PER_ANCESTRY = 200
MANIFEST_URL = "gs://fc-aou-datasets-controlled/v8/wgs/cram/manifest.csv"
OUTPUT_BASE = f"{WORKSPACE_BUCKET}/kir/calibration/sample_manifests"

print(f"\n{'='*80}")
print(f"KIR Calibration: Sample Selection")
print(f"{'='*80}")
print(f"Ancestries: {', '.join(ANCESTRIES)}")
print(f"Samples per ancestry: {SAMPLES_PER_ANCESTRY}")
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
print(f"\nAncestry distribution:")
ancestry_summary = ancestry_df.group_by('ancestry_pred_other').agg(pl.len().alias('count')).sort('count', descending=True)
print(ancestry_summary)

# Load manifest
print(f"\nDownloading manifest from {MANIFEST_URL}...")
import subprocess
subprocess.run(['gsutil', '-u', GOOGLE_PROJECT, 'cp', MANIFEST_URL, 'manifest.csv'], check=True)
manifest_df = pl.read_csv('manifest.csv')
print(f"✓ Loaded {len(manifest_df):,} samples from manifest")
print(f"Manifest columns: {manifest_df.columns}")

# ============================================================================
# Sample Selection per Ancestry
# ============================================================================

print(f"\n{'='*80}")
print("Selecting Calibration Samples")
print(f"{'='*80}")

ancestry_manifests = {}

for ancestry in ANCESTRIES:
    print(f"\n{ancestry}:")

    # Get person IDs for this ancestry
    ancestry_col = 'ancestry_pred_other'
    # Handle case-insensitive matching
    eur_ids = ancestry_df.filter(
        pl.col(ancestry_col).str.to_uppercase() == ancestry
    )['research_id'].to_list()

    if not eur_ids:
        print(f"  ⚠ No samples found with ancestry_pred_other == '{ancestry}'")
        print(f"  Available values: {ancestry_df[ancestry_col].unique().to_list()}")
        continue

    print(f"  Available samples: {len(eur_ids):,}")

    # Filter manifest to this ancestry
    ancestry_manifest = manifest_df.filter(pl.col('person_id').is_in(eur_ids))
    print(f"  In manifest: {len(ancestry_manifest)}")

    # Random sample
    if len(ancestry_manifest) < SAMPLES_PER_ANCESTRY:
        print(f"  ⚠ Warning: Only {len(ancestry_manifest)} samples available (need {SAMPLES_PER_ANCESTRY})")
        selected = ancestry_manifest
    else:
        selected = ancestry_manifest.sample(n=SAMPLES_PER_ANCESTRY, seed=42)

    print(f"  Selected: {len(selected)}")

    # Save manifest
    ancestry_manifests[ancestry] = selected

# ============================================================================
# Validate and Save Manifests
# ============================================================================

print(f"\n{'='*80}")
print("Saving Manifests to GCS")
print(f"{'='*80}")

storage_client = storage.Client(project=GOOGLE_PROJECT)
bucket = storage_client.bucket(WORKSPACE_BUCKET.replace('gs://', ''))

for ancestry, manifest in ancestry_manifests.items():
    if len(manifest) == 0:
        print(f"{ancestry}: Skipping (no samples)")
        continue

    # Validate required columns
    required_cols = ['person_id', 'cram_uri', 'cram_index_uri']
    if not all(col in manifest.columns for col in required_cols):
        print(f"✗ {ancestry}: Missing required columns. Have: {manifest.columns}")
        continue

    # Save as CSV
    csv_content = manifest.write_csv()
    blob_path = f"kir/calibration/sample_manifests/{ancestry}_samples.csv"
    blob = bucket.blob(blob_path)
    blob.upload_from_string(csv_content)

    print(f"✓ {ancestry}: Saved {len(manifest)} samples to gs://{WORKSPACE_BUCKET}/{blob_path}")

print(f"\n{'='*80}")
print("Complete")
print(f"{'='*80}")
print(f"\nNext steps:")
print(f"1. Run 02_calibration_pipeline.py to submit calibration jobs")
print(f"2. Wait for jobs to complete (~4 hours)")
print(f"3. Run 03_extract_thresholds.py to extract ancestry-specific thresholds")
