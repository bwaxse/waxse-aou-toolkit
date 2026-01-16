#!/usr/bin/env python
# coding: utf-8

"""
KIR Copy Number Calibration: Extract and Validate Thresholds
============================================================

Extracts thresholds.txt from completed calibration jobs.
Validates thresholds and saves ancestry-specific versions for production use.
"""

get_ipython().system('pip install tctk polars google-cloud-storage --upgrade -q')

import os
import polars as pl
from google.cloud import storage
import json

# Environment setup
WORKSPACE_BUCKET = os.getenv('WORKSPACE_BUCKET')
GOOGLE_PROJECT = os.getenv('GOOGLE_PROJECT')

print(f"Bucket: {WORKSPACE_BUCKET}")
print(f"Project: {GOOGLE_PROJECT}")

# ============================================================================
# Configuration
# ============================================================================

ANCESTRIES = ['EUR', 'AFR', 'AMR', 'EAS', 'SAS']
CALIBRATION_BASE = f"{WORKSPACE_BUCKET}/kir/calibration"
THRESHOLDS_OUTPUT_BASE = f"{WORKSPACE_BUCKET}/kir/thresholds"

print(f"\n{'='*80}")
print(f"KIR Calibration: Extract and Validate Thresholds")
print(f"{'='*80}")

# ============================================================================
# Helper Functions
# ============================================================================

def download_thresholds(bucket, ancestry):
    """Download thresholds.txt from calibration output."""
    source_path = f"kir/calibration/{ancestry}/thresholds.txt"
    print(f"  Downloading thresholds from gs://{bucket}/{source_path}...")

    storage_client = storage.Client(project=GOOGLE_PROJECT)
    bucket_obj = storage_client.bucket(bucket.replace('gs://', ''))
    blob = bucket_obj.blob(source_path)

    if not blob.exists():
        raise FileNotFoundError(f"Thresholds not found: gs://{bucket}/{source_path}")

    content = blob.download_as_string().decode('utf-8')
    return content

def validate_thresholds(content, ancestry):
    """Validate thresholds for reasonable values."""
    lines = [l for l in content.strip().split('\n') if l]
    print(f"  Lines: {len(lines)}")

    # Basic validation
    errors = []

    if 'nan' in content.lower():
        errors.append("Contains NaN values")

    if 'inf' in content.lower():
        errors.append("Contains Inf values")

    # Check for empty file
    if len(lines) < 2:
        errors.append("File appears empty or too short")

    if errors:
        for error in errors:
            print(f"    ⚠ {error}")
        return False
    else:
        print(f"  ✓ Validation passed")
        return True

def save_thresholds(bucket, ancestry, content):
    """Save thresholds to standard location."""
    dest_path = f"kir/thresholds/thresholds_{ancestry}.txt"
    print(f"  Saving to gs://{bucket}/{dest_path}...")

    storage_client = storage.Client(project=GOOGLE_PROJECT)
    bucket_obj = storage_client.bucket(bucket.replace('gs://', ''))
    blob = bucket_obj.blob(dest_path)
    blob.upload_from_string(content)

    print(f"  ✓ Saved")
    return dest_path

# ============================================================================
# Main Process
# ============================================================================

print(f"\n{'='*80}")
print("Processing Ancestries")
print(f"{'='*80}")

summary = {}

for ancestry in ANCESTRIES:
    print(f"\n{ancestry}:")

    try:
        # Download
        content = download_thresholds(WORKSPACE_BUCKET, ancestry)

        # Validate
        valid = validate_thresholds(content, ancestry)

        # Save
        if valid:
            dest_path = save_thresholds(WORKSPACE_BUCKET, ancestry, content)
            summary[ancestry] = {
                'status': 'complete',
                'path': f"gs://{WORKSPACE_BUCKET}/{dest_path}",
                'valid': True
            }
        else:
            summary[ancestry] = {
                'status': 'complete_with_warnings',
                'valid': False
            }

    except FileNotFoundError as e:
        print(f"  ✗ {e}")
        summary[ancestry] = {
            'status': 'missing',
            'error': str(e)
        }
    except Exception as e:
        print(f"  ✗ Error: {e}")
        summary[ancestry] = {
            'status': 'error',
            'error': str(e)
        }

# ============================================================================
# Summary Report
# ============================================================================

print(f"\n{'='*80}")
print("Summary")
print(f"{'='*80}")

complete_count = sum(1 for s in summary.values() if s['status'].startswith('complete'))
print(f"\nExtracted thresholds: {complete_count}/{len(ANCESTRIES)}")

print(f"\nStatus by ancestry:")
for ancestry, status_dict in summary.items():
    status = status_dict['status']
    if status == 'complete':
        print(f"  ✓ {ancestry}: {status_dict['path']}")
    elif status == 'complete_with_warnings':
        print(f"  ⚠ {ancestry}: Complete with warnings (not valid)")
    else:
        print(f"  ✗ {ancestry}: {status_dict.get('error', status)}")

# ============================================================================
# Next Steps
# ============================================================================

print(f"\n{'='*80}")
print("Next Steps")
print(f"{'='*80}")

if complete_count == len(ANCESTRIES):
    print(f"\n✓ All ancestry thresholds extracted and validated!")
    print(f"\nThresholds are ready for production use.")
    print(f"Location: {THRESHOLDS_OUTPUT_BASE}/")
    print(f"\nNext: Run 04_generate_production_batches.py to prepare production jobs")
else:
    print(f"\n⚠ Not all thresholds available ({complete_count}/{len(ANCESTRIES)})")
    print(f"Check calibration job logs for failures:")
    print(f"  gsutil ls {WORKSPACE_BUCKET}/dsub/logs/ | head -20")
    print(f"\nAfter debugging, re-run this script to extract completed thresholds")

# Save summary as JSON
summary_path = f"{WORKSPACE_BUCKET}/kir/calibration/thresholds_extraction_summary.json"
print(f"\nSummary saved to: gs://{summary_path}")
