#!/bin/bash
set -e

# ============================================================================
# KIR Copy Number Calibration Script
# ============================================================================
# Runs GATK PrintReads → kir-mapper map → kir-mapper ncopy on calibration samples
# Generates ancestry-specific thresholds.txt for use in production runs
#
# Environment variables (set by dsub):
#   - NUM_SAMPLES: Number of samples to process (200)
#   - PARALLEL_SAMPLES: Samples to process in parallel (10)
#   - THREADS_PER_SAMPLE: Threads per kir-mapper map (2)
#   - ANCESTRY: EUR, AFR, AMR, EAS, or SAS
#   - OUTPUT_DIR: GCS path for results
#   - INPUT_CRAM_1, INPUT_CRAM_2, ...: CRAM file GCS paths
#   - INPUT_CRAI_1, INPUT_CRAI_2, ...: CRAI file GCS paths
# ============================================================================

# Configuration (from environment variables)
NUM_SAMPLES=${NUM_SAMPLES:-200}
PARALLEL_SAMPLES=${PARALLEL_SAMPLES:-10}
THREADS_PER_SAMPLE=${THREADS_PER_SAMPLE:-2}
ANCESTRY=${ANCESTRY:-EUR}

# Verify critical environment variables
if [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: OUTPUT_DIR not set"
    exit 1
fi

if [ -z "$GOOGLE_PROJECT" ]; then
    echo "ERROR: GOOGLE_PROJECT not set"
    exit 1
fi

# Timing log - per-sample to avoid race conditions during parallel execution
TIMING_DIR="/tmp/kir_timing_$$"  # Use process ID for uniqueness
mkdir -p $TIMING_DIR

log_time() {
    local step=$1
    local sample=$2
    local start=$3
    local end=$4
    local duration=$((end - start))
    local sample_log="$TIMING_DIR/timing_${sample}.csv"

    # Create header if file doesn't exist
    if [ ! -f "$sample_log" ]; then
        echo "step,sample,start_epoch,end_epoch,duration_sec" > "$sample_log"
    fi

    echo "${step},${sample},${start},${end},${duration}" >> "$sample_log"
}

# ============================================================================
# Helper Functions
# ============================================================================

run_gatk() {
    local idx=$1
    local cram=$2
    local crai=$3

    local start=$(date +%s)
    echo "  [GATK] Processing sample ${idx}..."

    gatk PrintReads \
        -R gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta \
        -I "$cram" \
        --read-index "$crai" \
        -L chr19:54000000-55100000 \
        --gcs-project-for-requester-pays $GOOGLE_PROJECT \
        -O "sample${idx}_chr19.bam" 2>&1

    local end=$(date +%s)
    log_time "gatk" "sample${idx}" "$start" "$end"
}

run_kirmap() {
    local idx=$1

    local start=$(date +%s)
    echo "  [KIR-MAP] Processing sample ${idx}..."

    kir-mapper map \
        -bam "sample${idx}_chr19.bam" \
        -sample "sample${idx}" \
        -output ./kir_output \
        -threads $THREADS_PER_SAMPLE 2>&1

    local end=$(date +%s)
    log_time "kirmap" "sample${idx}" "$start" "$end"
}

run_gatk_with_env() {
    local idx=$1
    local cram_var="INPUT_CRAM_$idx"
    local crai_var="INPUT_CRAI_$idx"
    run_gatk "$idx" "${!cram_var}" "${!crai_var}"
}

run_kirmap_with_env() {
    local idx=$1
    run_kirmap "$idx"
}

export -f run_gatk run_kirmap run_gatk_with_env run_kirmap_with_env log_time
export TIMING_LOG THREADS_PER_SAMPLE GOOGLE_PROJECT

# ============================================================================
# Main Pipeline
# ============================================================================

echo "============================================================================"
echo "KIR Copy Number Calibration"
echo "============================================================================"
echo "Ancestry: $ANCESTRY"
echo "Samples: $NUM_SAMPLES"
echo "Parallel: $PARALLEL_SAMPLES"
echo "Threads/sample: $THREADS_PER_SAMPLE"
echo "============================================================================"

# Step 1: Parallel GATK conversions
echo ""
echo "=== Step 1: GATK PrintReads (parallel=$PARALLEL_SAMPLES) ==="
GATK_START=$(date +%s)

if command -v parallel &> /dev/null; then
    seq 1 $NUM_SAMPLES | parallel -j $PARALLEL_SAMPLES 'run_gatk_with_env {}'
else
    echo "GNU parallel not found, falling back to sequential processing"
    for idx in $(seq 1 $NUM_SAMPLES); do
        run_gatk_with_env "$idx"
    done
fi

GATK_END=$(date +%s)
log_time "gatk_total" "all" "$GATK_START" "$GATK_END"

# Step 2: Parallel kir-mapper map
echo ""
echo "=== Step 2: kir-mapper map (parallel=$PARALLEL_SAMPLES) ==="
KIRMAP_START=$(date +%s)

if command -v parallel &> /dev/null; then
    seq 1 $NUM_SAMPLES | parallel -j $PARALLEL_SAMPLES 'run_kirmap_with_env {}'
else
    echo "GNU parallel not found, falling back to sequential processing"
    for idx in $(seq 1 $NUM_SAMPLES); do
        run_kirmap_with_env "$idx"
    done
fi

KIRMAP_END=$(date +%s)
log_time "kirmap_total" "all" "$KIRMAP_START" "$KIRMAP_END"

# ============================================================================
# Consolidate Timing Logs
# ============================================================================

echo ""
echo "=== Consolidating Timing Logs ==="
TIMING_LOG="/tmp/timing_log_combined.csv"
echo "step,sample,start_epoch,end_epoch,duration_sec" > $TIMING_LOG

# Combine all per-sample logs (skip headers)
for log_file in $TIMING_DIR/timing_*.csv; do
    if [ -f "$log_file" ]; then
        tail -n +2 "$log_file" >> $TIMING_LOG
    fi
done

echo "✓ Consolidated $(wc -l < $TIMING_LOG) timing entries"

# Cleanup per-sample logs
rm -rf $TIMING_DIR

# Step 3: kir-mapper ncopy (batch operation, generates thresholds.txt)
echo ""
echo "=== Step 3: kir-mapper ncopy (generates thresholds.txt) ==="
NCOPY_START=$(date +%s)

kir-mapper ncopy -output ./kir_output -threads 20 2>&1

NCOPY_END=$(date +%s)
echo "ncopy completed in $((NCOPY_END - NCOPY_START)) seconds"

# ============================================================================
# Output Collection
# ============================================================================

echo ""
echo "=== Copying Results to Output Directory ==="
mkdir -p $OUTPUT_DIR

# Copy kir-mapper outputs
if [ -d "./kir_output" ]; then
    cp -r kir_output/* $OUTPUT_DIR/
    echo "✓ Copied kir-mapper outputs"
fi

# Copy timing log
if [ -f "$TIMING_LOG" ]; then
    cp $TIMING_LOG $OUTPUT_DIR/timing_log.csv
    echo "✓ Copied timing log"
fi

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "=== Calibration Complete ==="
echo ""
echo "Timing Summary:"
cat $TIMING_LOG | column -t -s ','

echo ""
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Key outputs:"
if [ -f "$OUTPUT_DIR/thresholds.txt" ]; then
    echo "✓ thresholds.txt: $OUTPUT_DIR/thresholds.txt"
    echo "  (Use this in production runs)"
fi
if [ -f "$OUTPUT_DIR/timing_log.csv" ]; then
    echo "✓ timing_log.csv: $OUTPUT_DIR/timing_log.csv"
fi

echo ""
echo "Next: Extract thresholds.txt from all ancestries using 03_extract_thresholds.py"
