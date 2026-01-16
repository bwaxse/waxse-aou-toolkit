#!/usr/bin/env python
# coding: utf-8

"""
KIR Copy Number Production: Collect and Aggregate Results
========================================================

After all production jobs complete, collects and combines ncopy results.
- Downloads copy_numbers.table.txt from all batches
- Downloads presence.table.txt from all batches
- Combines into single ancestry-wide tables
- Generates QC summary (missing samples, row counts, etc.)

Output: Final results in gs://{WORKSPACE_BUCKET}/kir/final_results/
"""

get_ipython().system('pip install tctk polars google-cloud-storage --upgrade -q')

import os
import polars as pl
from google.cloud import storage
from io import StringIO

# Environment setup
WORKSPACE_BUCKET = os.getenv('WORKSPACE_BUCKET')
GOOGLE_PROJECT = os.getenv('GOOGLE_PROJECT')

print(f"Bucket: {WORKSPACE_BUCKET}")
print(f"Project: {GOOGLE_PROJECT}")

# ============================================================================
# Configuration
# ============================================================================

ANCESTRIES = ['EUR', 'AFR', 'AMR', 'EAS', 'SAS']
PRODUCTION_BASE = f"{WORKSPACE_BUCKET}/kir/production"
OUTPUT_BASE = f"{WORKSPACE_BUCKET}/kir/final_results"

print(f"\n{'='*80}")
print(f"KIR Production: Collect and Aggregate Results")
print(f"{'='*80}")
print(f"Input base: {PRODUCTION_BASE}/")
print(f"Output base: {OUTPUT_BASE}/")

# ============================================================================
# Helper Functions
# ============================================================================

def list_batches(storage_client, bucket_name, ancestry):
    """List all completed batches for ancestry."""
    prefix = f"kir/production/{ancestry}/"
    bucket_obj = storage_client.bucket(bucket_name)
    blobs = storage_client.list_blobs(bucket_obj, prefix=prefix, delimiter='/')

    batches = []
    for page in blobs.pages:
        for prefix_obj in page.prefixes:
            batch_id = prefix_obj.rstrip('/').split('/')[-1]
            batches.append(batch_id)

    return sorted(batches)

def download_file(storage_client, bucket_name, gcs_path):
    """Download file from GCS as string."""
    bucket_obj = storage_client.bucket(bucket_name)
    blob = bucket_obj.blob(gcs_path)

    if not blob.exists():
        return None

    return blob.download_as_string().decode('utf-8')

def collect_batch_results(storage_client, bucket_name, ancestry, batch_id, file_type):
    """Download ncopy results from a single batch."""
    # file_type: 'copy_numbers' or 'presence'
    file_name = f"{file_type}.table.txt"
    gcs_path = f"kir/production/{ancestry}/{batch_id}/ncopy/{file_name}"

    content = download_file(storage_client, bucket_name, gcs_path)
    if content is None:
        return None

    try:
        df = pl.read_csv(StringIO(content), separator='\t')
        return df
    except Exception as e:
        print(f"    ⚠ Error parsing {gcs_path}: {e}")
        return None

# ============================================================================
# Main Collection Process
# ============================================================================

print(f"\n{'='*80}")
print("Collecting Results by Ancestry")
print(f"{'='*80}")

storage_client = storage.Client(project=GOOGLE_PROJECT)
bucket_name = WORKSPACE_BUCKET.replace('gs://', '')

final_results = {}
qc_summary = {}

for ancestry in ANCESTRIES:
    print(f"\n{ancestry}:")

    # List batches
    batches = list_batches(storage_client, bucket_name, ancestry)
    if not batches:
        print(f"  ⚠ No batches found")
        continue

    print(f"  Batches: {len(batches)}")

    # Collect copy_numbers.table.txt
    print(f"  Collecting copy_numbers.table.txt...")
    copy_numbers_dfs = []
    missing_batches = []

    for batch_id in batches:
        df = collect_batch_results(storage_client, bucket_name, ancestry, batch_id, 'copy_numbers')
        if df is not None:
            copy_numbers_dfs.append(df)
        else:
            missing_batches.append(batch_id)

    if copy_numbers_dfs:
        # Combine
        combined_copy_numbers = pl.concat(copy_numbers_dfs)
        print(f"    ✓ Combined {len(copy_numbers_dfs)} batches ({len(combined_copy_numbers)} rows)")
        final_results[f"{ancestry}_copy_numbers"] = combined_copy_numbers
    else:
        print(f"    ✗ No copy_numbers data found")

    # Collect presence.table.txt
    print(f"  Collecting presence.table.txt...")
    presence_dfs = []

    for batch_id in batches:
        df = collect_batch_results(storage_client, bucket_name, ancestry, batch_id, 'presence')
        if df is not None:
            presence_dfs.append(df)

    if presence_dfs:
        # Combine
        combined_presence = pl.concat(presence_dfs)
        print(f"    ✓ Combined {len(presence_dfs)} batches ({len(combined_presence)} rows)")
        final_results[f"{ancestry}_presence"] = combined_presence
    else:
        print(f"    ✗ No presence data found")

    # QC Summary
    qc_summary[ancestry] = {
        'total_batches': len(batches),
        'batches_with_data': len(copy_numbers_dfs),
        'missing_batches': len(missing_batches),
        'total_rows': len(combined_copy_numbers) if copy_numbers_dfs else 0,
    }

# ============================================================================
# Save Results to GCS
# ============================================================================

print(f"\n{'='*80}")
print("Saving Combined Results")
print(f"{'='*80}")

bucket_obj = storage_client.bucket(bucket_name)

saved_files = {}

for ancestry in ANCESTRIES:
    print(f"\n{ancestry}:")

    # Save copy_numbers
    df = final_results.get(f"{ancestry}_copy_numbers")
    if df is not None and len(df) > 0:
        csv_content = df.write_csv(separator='\t')
        blob_path = f"kir/final_results/{ancestry}/copy_numbers.table.txt"
        blob = bucket_obj.blob(blob_path)
        blob.upload_from_string(csv_content)
        saved_files[f"{ancestry}_copy_numbers"] = f"gs://{WORKSPACE_BUCKET}/{blob_path}"
        print(f"  ✓ copy_numbers.table.txt: {len(df)} rows")

    # Save presence
    df = final_results.get(f"{ancestry}_presence")
    if df is not None and len(df) > 0:
        csv_content = df.write_csv(separator='\t')
        blob_path = f"kir/final_results/{ancestry}/presence.table.txt"
        blob = bucket_obj.blob(blob_path)
        blob.upload_from_string(csv_content)
        saved_files[f"{ancestry}_presence"] = f"gs://{WORKSPACE_BUCKET}/{blob_path}"
        print(f"  ✓ presence.table.txt: {len(df)} rows")

# ============================================================================
# QC Summary Report
# ============================================================================

print(f"\n{'='*80}")
print("QC Summary")
print(f"{'='*80}")

for ancestry in ANCESTRIES:
    if ancestry not in qc_summary:
        continue

    summary = qc_summary[ancestry]
    print(f"\n{ancestry}:")
    print(f"  Total batches: {summary['total_batches']}")
    print(f"  Batches with data: {summary['batches_with_data']}")
    print(f"  Missing batches: {summary['missing_batches']}")
    print(f"  Total samples: {summary['total_rows']:,}")

    if summary['missing_batches'] > 0:
        print(f"  ⚠ WARNING: {summary['missing_batches']} batches missing - check job logs")

# ============================================================================
# Summary
# ============================================================================

print(f"\n{'='*80}")
print("Collection Complete")
print(f"{'='*80}")

print(f"\nSaved files:")
for file_id, path in saved_files.items():
    print(f"  {file_id}: {path}")

total_rows = sum(qc_summary[a]['total_rows'] for a in ANCESTRIES if a in qc_summary)
print(f"\nTotal samples processed: {total_rows:,}")

# ============================================================================
# Next Steps
# ============================================================================

print(f"\n{'='*80}")
print("Next Steps")
print(f"{'='*80}")

print(f"\n1. Download results for analysis:")
print(f"   gsutil -m cp -r gs://{WORKSPACE_BUCKET}/kir/final_results/ ./")

print(f"\n2. Final checks:")
print(f"   - Verify row counts match expected sample counts")
print(f"   - Check for duplicate samples across batches")
print(f"   - Review KIR copy number distributions by ancestry")

print(f"\n3. Archive intermediate results (optional):")
print(f"   - Delete batch directories: gsutil -m rm -r {PRODUCTION_BASE}/**/batch_*/")
print(f"   - Keep calibration results for future reference")

print(f"\n✓ KIR copy number calling complete!")
print(f"  Results: gs://{WORKSPACE_BUCKET}/kir/final_results/")
