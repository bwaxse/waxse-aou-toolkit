# KIR-Mapper Parallelization Test Suite

## Overview

This test suite orchestrates multiple dsub jobs to determine the optimal parallelization strategy for KIR genotyping on All of Us data. It tests different combinations of:
- **Sample counts**: 4, 10, 20 samples
- **Parallelization factors**: 2, 4, 5, 10 parallel jobs per VM
- **Machine types**: n2-standard-4, n2-standard-8, n2-standard-16

Each configuration runs GATK PrintReads and kir-mapper with detailed timing instrumentation to measure:
- Per-step runtime (GATK, kir-mapper map, ncopy, genotype)
- Time per sample
- Total execution time
- Resource utilization (via machine type)

## Files

- **`<gatk_kir> dsub parallelization test.py`**: Main test orchestration script
- **`kir_instrumented_test.sh`**: Bash script with GNU parallel orchestration and timing (referenced by main script)
- **`<gatk_kir> dsub test.py`**: Original working single-job submission (reference/comparison)

## Prerequisites

1. **All of Us Research Workbench access** with genomic data
2. **Environment variables** configured (run `00_setup_workspace.ipynb` first):
   - `WORKSPACE_BUCKET`
   - `WORKSPACE_CDR`
   - `GOOGLE_PROJECT`

3. **Python packages** (installed automatically):
   - `tctk` (AoUTools for dsub)
   - `polars`
   - `google-cloud-storage`

4. **GCS manifest file** (automatically downloaded):
   - `gs://fc-aou-datasets-controlled/v8/wgs/cram/manifest.csv`

## Quick Start

### Step 1: Run Setup Notebook
Before running any analysis, initialize your workspace:
```python
# In Verily Workbench, run this first:
# Open and execute: _reference/verily/00_setup_workspace.ipynb
```

### Step 2: Run Test Script in Notebook
Convert the Python script to a notebook and upload to Verily:
```bash
# Claude: Convert `lc_kir/<gatk_kir> dsub parallelization test.py` to a notebook
```

Then in Verily:
1. Upload the generated notebook to your workspace
2. Open in JupyterLab
3. Run all cells sequentially

### Step 3: Monitor Jobs
The script will submit 6 dsub jobs with different configurations. Monitor using:
```python
# Check status
job.check_status(streaming=True)

# Or manually in Verily interface: Monitor > Batch Submission Jobs
```

### Step 4: Collect and Analyze Results
After jobs complete (~1-5 hours depending on machine type and queue):
1. Script automatically retrieves timing logs
2. Generates summary table with per-sample metrics
3. Saves results to: `{WORKSPACE_BUCKET}/dsub/results/kir_parallel/`

## Configuration Details

### Test Configurations

```
Size  Samples  Parallel  Machine       vCPU  RAM    Focus
─────────────────────────────────────────────────────────────────
Small    4        2      n2-std-4       4    16GB   Sequential vs Low Parallel
         4        4      n2-std-8       8    32GB   Max Parallel (small)

Medium  10        2      n2-std-4       4    16GB   Conservative parallel
        10        5      n2-std-8       8    32GB   Mid-range parallel

Large   20        4      n2-std-8       8    32GB   High sample count, moderate parallel
        20       10      n2-std-16     16    64GB   High sample count, high parallel
```

### Why n2-standard?

- **Balanced CPU/memory ratio**: 4GB RAM per vCPU (vs ~1GB for c2d-highcpu)
- **Better for GATK**: CRAM decompression is memory-intensive
- **Safer**: Reduces OOM (Out of Memory) risk during testing
- **Performance comparable**: CPU performance is similar to c2d for this workload
- **Future comparison**: Can test c2d-highcpu series later if memory usage proves low

## Expected Outputs

Each configuration creates a timestamped results directory:
```
{WORKSPACE_BUCKET}/dsub/results/kir_parallel/
├── n4_p2_n2-standard-4/
│   ├── timing_log.csv              # Per-step timing data
│   ├── map/                        # Per-sample BAM mappings
│   ├── ncopy/                      # Copy number calls
│   ├── genotype/                   # SNP/allele calls
│   └── haplotype/                  # (skipped for <50 samples)
├── n4_p4_n2-standard-8/
├── n10_p2_n2-standard-4/
├── n10_p5_n2-standard-8/
├── n20_p4_n2-standard-8/
└── n20_p10_n2-standard-16/
```

### Timing Log Format

`timing_log.csv` contains:
```csv
step,sample,start_epoch,end_epoch,duration_sec
gatk,sample01,1705335600,1705335720,120
gatk,sample02,1705335600,1705335840,240
...
gatk_total,all,1705335600,1705336200,600
kirmap,sample01,1705336200,1705336380,180
...
kirmap_total,all,1705336200,1705338000,1800
ncopy,all,1705338000,1705338300,300
genotype,all,1705338300,1705338900,600
```

## Analysis and Interpretation

### Key Metrics to Compare

1. **Total Pipeline Time**: `gatk_total + kirmap_total + ncopy + genotype`
   - Lower is better
   - Dominated by per-sample operations (GATK + kir-mapper map)

2. **Time Per Sample**: `Total Time / Number of Samples`
   - Shows scaling efficiency
   - Helps project time for larger cohorts

3. **Parallelization Efficiency**: `(Serial Time) / (Parallel Time × Parallel_Factor)`
   - ~1.0 = perfect scaling
   - <1.0 = overhead from parallelization

4. **Resource Utilization**: `(CPU Time Used) / (Available CPU Time)`
   - Shows if machine is truly utilizing all vCPUs

5. **Cost Per Sample**: `(Machine Hourly Rate × Job Duration) / Number of Samples`
   - Lowest cost option may not have best wall-clock time
   - Choose based on your priority: speed vs. cost

### Example Analysis

```
Configuration           Total Time  Per-Sample  Cost/Sample
────────────────────────────────────────────────────────────
4 samples, 2p, n2-std4  1500 sec    375 sec     $0.15
4 samples, 4p, n2-std8  900 sec     225 sec     $0.18
10 samples, 2p, std4    3000 sec    300 sec     $0.12  ← Best cost/sample
10 samples, 5p, std8    1500 sec    150 sec     $0.19
20 samples, 4p, std8    5000 sec    250 sec     $0.25
20 samples, 10p, std16  2500 sec    125 sec     $0.40
```

**Recommendations:**
- For speed-critical: Use configuration n20_p10 (2.5 min per sample)
- For cost-efficient: Use configuration n10_p2 (5 min per sample, 60% cheaper)
- For balanced: Use configuration n10_p5 (2.5 min per sample, moderate cost)

## Pre-Test Verification (Optional)

Before running the full test suite, verify the Docker image has required tools:

```bash
# From your local machine (if Docker installed):
docker pull phetk/gatk-kirmapper:0.1

# Check for GNU parallel
docker run --rm phetk/gatk-kirmapper:0.1 which parallel
docker run --rm phetk/gatk-kirmapper:0.1 parallel --version

# List all available tools
docker run --rm phetk/gatk-kirmapper:0.1 bash -c "which gatk kir-mapper samtools bcftools"
```

**If GNU parallel is missing:**
- The script has a fallback: it will use sequential processing if `parallel` is not found
- Timing instrumentation still works
- Results will be valid but sequential configuration

## Ancestry Filtering

The script **automatically filters to EUR ancestry only** by:
1. Loading `{WORKSPACE_BUCKET}/data/ancestry_metadata.tsv`
2. Extracting EUR person IDs: `ancestry_df.filter(col('ancestry_pred_other') == 'EUR')`
3. Selecting first 20 EUR samples from manifest
4. Using these for all test configurations

This ensures consistent, ancestry-specific results aligned with your research focus.

## Troubleshooting

### Job Submission Fails
- Check `WORKSPACE_BUCKET`, `GOOGLE_PROJECT` are set correctly
- Verify manifest.csv was downloaded: `gsutil ls gs://fc-aou-datasets-controlled/v8/wgs/cram/manifest.csv`
- Check GCS bucket access permissions

### Jobs Timeout or OOM (Out of Memory)
- GATK CRAM decompression can be memory-intensive
- Increase machine type (e.g., n2-standard-8 → n2-standard-16)
- Reduce parallel factor per VM
- Check if CRAM files are corrupted or unusually large

### Timing Log Not Generated
- Check job's stderr for errors (available in dsub interface)
- Verify `kir-mapper` ran successfully (check `map/`, `ncopy/`, `genotype/` directories)
- If pipeline failed, timing_log.csv may be incomplete

### GNU Parallel Not Found
- Script automatically falls back to sequential processing
- Results are still valid (just not parallelized)
- Consider updating Docker image or submitting parallel via separate tool

## Extending the Test

### Add More Sample Sizes
Edit `TEST_CONFIGURATIONS` in the script:
```python
TEST_CONFIGURATIONS = [
    # Add to existing...
    (50, 10, 2, 'n2-standard-16'),   # Large: enables haplotype step
    (100, 20, 2, 'n2-highcpu-32'),   # Production scale
]
```

### Test Different Machine Types
```python
(20, 4, 2, 'c2d-highcpu-8'),   # Compute-optimized
(20, 4, 2, 'n2-highmem-16'),   # Memory-optimized
```

### Test Different Reference Regions
Edit the bash script generator to change the genomic region:
```bash
# Current: chr19:54000000-55100000
# Change to: -L chr19  (entire chromosome)
# Or: -L chr1,chr2,chr3  (multiple chromosomes)
```

## Output Files Summary

| File | Description | Access |
|------|-------------|--------|
| `timing_log.csv` | Per-step timing data | Download from GCS |
| `map/*.bam` | GATK-converted BAM files | GCS (intermediate, can delete) |
| `ncopy/copy_numbers.table.txt` | KIR copy number calls | Download for analysis |
| `ncopy/presence.table.txt` | KIR presence/absence calls | Download for analysis |
| `genotype/*.vcf*` | SNP/allele calls in VCF format | Download for analysis |
| `genotype/alleles.txt` | Assigned KIR alleles | Download for analysis |

## Notes

- **Haplotype step**: Only runs if `NUM_SAMPLES >= 50` (current tests stop at 20, so this step is skipped)
- **Preemptible instances**: Test suite uses non-preemptible for reliability; production can use preemptible for 60-80% cost savings
- **Cost estimate**: Typical test run costs $2-5 USD (varies by queue time and machine type)
- **Timeline**: 4-6 hours from submission to completion (dominated by GATK and individual CRAM file processing)

## References

- **KIR-Mapper Documentation**: `kir-mapper/MANUAL.md`
- **All of Us Genomic Data**: [All of Us Documentation](https://docs.researchallofus.org/)
- **GATK PrintReads**: [GATK Tool Documentation](https://gatk.broadinstitute.org/)
- **dsub**: [Docker Submit Tool](https://github.com/DataBiosphere/dsub)

## Questions?

If you encounter issues or have questions:
1. Check the `timing_log.csv` for per-step performance
2. Review dsub job logs: `gsutil cat gs://{BUCKET}/dsub/logs/{JOB_ID}.log`
3. Examine kir-mapper outputs for data quality issues
4. Reach out to Bennett Waxse (bennett.waxse@gmail.com) or open a GitHub issue
