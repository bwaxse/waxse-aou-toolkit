# KIR-Mapper Parallelization Test Suite - Implementation Summary

## Deliverables

### 1. Main Test Script
**File**: `<gatk_kir> dsub parallelization test.py`

A comprehensive Python script for orchestrating parallelization testing across multiple dsub configurations. The script:

- **Loads data**:
  - Ancestry metadata from `{WORKSPACE_BUCKET}/data/ancestry_metadata.tsv`
  - Manifest from `gs://fc-aou-datasets-controlled/v8/wgs/cram/manifest.csv`
  - Filters to EUR ancestry samples only

- **Submits 6 test configurations**:
  ```
  (4 samples, 2 parallel, n2-standard-4)
  (4 samples, 4 parallel, n2-standard-8)
  (10 samples, 2 parallel, n2-standard-4)
  (10 samples, 5 parallel, n2-standard-8)
  (20 samples, 4 parallel, n2-standard-8)
  (20 samples, 10 parallel, n2-standard-16)
  ```

- **Generates instrumented bash scripts** with:
  - GATK PrintReads for CRAM→BAM conversion
  - kir-mapper map, ncopy, genotype steps
  - Per-step timing instrumentation (start_epoch, end_epoch, duration_sec)
  - Fallback to sequential processing if GNU parallel unavailable

- **Submits via dsub** using TCTK (AoUTools) library:
  - Google Batch provider
  - Docker image: `phetk/gatk-kirmapper:0.1`
  - Dynamic environment variable construction
  - Output collection to GCS

- **Analyzes results**:
  - Retrieves timing logs from completed jobs
  - Calculates per-sample metrics
  - Generates comparison table

### 2. Documentation
**File**: `PARALLELIZATION_TEST_README.md`

Comprehensive user guide including:
- Quick start workflow
- Configuration details and rationale
- Expected output structure
- Analysis and interpretation guidance
- Troubleshooting section
- Extension examples
- References

### 3. Implementation Summary
**File**: `IMPLEMENTATION_SUMMARY.md` (this file)

Overview of deliverables and technical approach.

## Key Features

### Architecture
- **Modular design**: Separate functions for data loading, script generation, job submission, results parsing
- **EUR ancestry filtering**: Automatically selects EUR samples from manifest
- **Dynamic configuration**: Test matrix easily extensible for more configurations
- **Robust error handling**: Graceful degradation if GNU parallel missing
- **Comprehensive timing**: Logs each step separately for detailed analysis

### Test Coverage
- **Small samples** (4): Tests overhead of parallelization setup
- **Medium samples** (10): Tests efficiency at modest scale
- **Large samples** (20): Tests scaling to production-relevant size
- **Parallelization factors** (2-10): Tests efficiency curve
- **Machine types** (n2-standard-4/8/16): Tests cost-performance tradeoff

### Automation
- Manifest download and caching
- Dynamic CRAM file cycling (if more samples than available CRAMs)
- Automatic timing log retrieval
- Summary table generation
- Per-sample metric calculation

## Technical Decisions

### Machine Type Selection: n2-standard over c2d-highcpu

**Rationale**:
- GATK PrintReads: Memory-intensive CRAM decompression (typically 5-10GB for chr19 region)
- kir-mapper: Alignment processing with temporary index files
- n2-standard: 4GB RAM per vCPU (vs ~1GB for c2d)
- Result: Reduced OOM risk, comparable CPU performance

**Testing approach**:
- n2-standard for reliability during testing
- Can compare c2d-highcpu later if memory usage proves low

### Parallelization Strategy

**Parallelizable steps**:
- GATK PrintReads: Per-sample (parallelizable)
- kir-mapper map: Per-sample (parallelizable)

**Batch operations** (limit parallelization benefit):
- ncopy: Batch operation across all samples
- genotype: Batch operation across all samples
- haplotype: Batch operation (requires ≥50 samples)

**Result**: Focus testing on GATK + map (per-sample) steps where parallelization provides most benefit

### Ancestry Filtering

**Approach**: EUR-only from start
- Reduces variability from ancestry effects
- Aligns with common GWAS pattern (ancestry-stratified)
- Simplifies manifest handling (known EUR count)
- Can extend to other ancestries easily

**Implementation**:
1. Load `ancestry_pred_other` column from metadata
2. Filter manifest to matching person_ids
3. Select first N EUR samples for each test

### Instrumentation

**Timing approach**:
- Unix epoch timestamps (resistant to time zone issues)
- Per-step granularity (identify bottlenecks)
- CSV format (easy parsing and analysis)
- Includes both individual and aggregate timings

**Fallback for missing GNU parallel**:
- Sequential processing with same timing instrumentation
- Results still valid (just not parallelized)
- Easy to identify if parallelization is missing

## Implementation Complexity

### What's Included
- ✅ Ancestry filtering (EUR-only)
- ✅ Manifest processing with cycling
- ✅ 6-configuration test matrix
- ✅ Dynamic bash script generation
- ✅ dsub job submission loop
- ✅ Timing instrumentation
- ✅ Results retrieval and parsing
- ✅ Summary table generation
- ✅ Comprehensive documentation

### What's Not Included (by design)
- ❌ Database integration (relies on GCS)
- ❌ Automatic cost optimization (user analyzes and decides)
- ❌ Machine learning for prediction (empirical measurement instead)
- ❌ Multi-ancestry testing (can extend easily)
- ❌ Preemptible instance testing (focus on reliability first)

## Usage Workflow

### Phase 1: Setup
1. Run `00_setup_workspace.ipynb` to initialize environment
2. Clone/pull this repository

### Phase 2: Execute
1. Convert `.py` to notebook: `jupyter nbconvert --to notebook`
2. Upload to Verily Workbench
3. Run notebook cells sequentially
4. Monitor job submissions via dsub interface

### Phase 3: Analyze
1. Wait for jobs to complete (1-5 hours)
2. Script automatically retrieves timing logs
3. Review generated comparison table
4. Download detailed results from GCS

### Phase 4: Decision
1. Compare configurations by time/cost
2. Choose optimal setting for full cohort
3. Scale up with selected configuration

## File Structure

```
lc_kir/
├── <gatk_kir> dsub parallelization test.py    # Main script
├── kir_instrumented_test.sh                    # Referenced bash script
├── <gatk_kir> dsub test.py                     # Original single-job reference
├── PARALLELIZATION_TEST_README.md              # User guide
└── IMPLEMENTATION_SUMMARY.md                   # This file
```

## Verification Checklist

- [x] Syntax validation: Python script parseable
- [x] Import validation: All required packages listed
- [x] Function validation: All helper functions present
- [x] Configuration validation: Test matrix well-defined
- [x] Script generation: Bash script template complete
- [x] Job submission: dsub integration correct
- [x] Results parsing: Timing log parsing implemented
- [x] Documentation: Comprehensive README provided
- [ ] End-to-end test: Awaiting user execution

## Next Steps

### For User
1. Read `PARALLELIZATION_TEST_README.md`
2. Convert script to notebook
3. Upload to Verily Workbench
4. Run all cells
5. Monitor and analyze results

### Optional Enhancements (if needed)
- Add other ancestry groups to test matrix
- Compare c2d-highcpu vs n2-standard
- Test preemptible instances for cost optimization
- Test different genomic regions
- Automate cost calculation with GCP pricing API
- Create visualization dashboard for results

## Known Limitations

1. **GNU parallel dependency**: If not in Docker image, falls back to sequential
2. **Haplotype step skipped**: Testing up to 20 samples (requires ≥50)
3. **Single region testing**: chr19:54000000-55100000 (can extend)
4. **No real-time monitoring**: Results available after job completion
5. **Manual cost analysis**: User must compare configurations

## Support

- Refer to `PARALLELIZATION_TEST_README.md` for troubleshooting
- Check dsub logs in GCS for job-specific errors
- Review kir-mapper outputs for data quality issues
- Contact Bennett Waxse for further assistance

---

**Created**: January 2026
**Script Purpose**: Determine optimal KIR genotyping parallelization strategy for All of Us data
**Status**: Ready for testing in Verily Workbench
