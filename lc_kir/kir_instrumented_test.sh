%%writefile kir_instrumented_test.sh
#!/bin/bash
set -e

# Configuration (passed as env vars)
THREADS_PER_SAMPLE=${THREADS_PER_SAMPLE:-2}
PARALLEL_SAMPLES=${PARALLEL_SAMPLES:-2}
NUM_SAMPLES=${NUM_SAMPLES:-4}

# Timing log
TIMING_LOG="/tmp/timing_log.csv"
echo "step,sample,start_epoch,end_epoch,duration_sec" > $TIMING_LOG

log_time() {
    local step=$1
    local sample=$2
    local start=$3
    local end=$4
    local duration=$((end - start))
    echo "${step},${sample},${start},${end},${duration}" >> $TIMING_LOG
}

# Create sample list from env vars (INPUT_CRAM_1, INPUT_CRAM_2, etc.)
SAMPLE_LIST=()
for i in $(seq 1 $NUM_SAMPLES); do
    cram_var="INPUT_CRAM_$i"
    crai_var="INPUT_CRAI_$i"
    if [ -n "${!cram_var}" ]; then
        SAMPLE_LIST+=("$i:${!cram_var}:${!crai_var}")
    fi
done

# Function: GATK conversion for one sample
run_gatk() {
    local idx=$1
    local cram=$2
    local crai=$3
    local start=$(date +%s)
    
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

# Function: kir-mapper map for one sample
run_kirmap() {
    local idx=$1
    local start=$(date +%s)
    
    kir-mapper map \
        -bam "sample${idx}_chr19.bam" \
        -sample "sample${idx}" \
        -output ./kir_output \
        -threads $THREADS_PER_SAMPLE 2>&1
    
    local end=$(date +%s)
    log_time "kirmap" "sample${idx}" "$start" "$end"
}

export -f run_gatk run_kirmap log_time
export TIMING_LOG THREADS_PER_SAMPLE GOOGLE_PROJECT

# Step 1: Parallel GATK conversions
echo "=== Starting GATK conversions (parallel=$PARALLEL_SAMPLES) ==="
GATK_START=$(date +%s)

printf '%s\n' "${SAMPLE_LIST[@]}" | \
    parallel -j $PARALLEL_SAMPLES --colsep ':' \
    'run_gatk {1} {2} {3}'

GATK_END=$(date +%s)
log_time "gatk_total" "all" "$GATK_START" "$GATK_END"

# Step 2: Parallel kir-mapper map
echo "=== Starting kir-mapper map (parallel=$PARALLEL_SAMPLES) ==="
KIRMAP_START=$(date +%s)

seq 1 $NUM_SAMPLES | \
    parallel -j $PARALLEL_SAMPLES 'run_kirmap {}'

KIRMAP_END=$(date +%s)
log_time "kirmap_total" "all" "$KIRMAP_START" "$KIRMAP_END"

# Step 3: ncopy (batch operation)
echo "=== Running ncopy ==="
NCOPY_START=$(date +%s)
kir-mapper ncopy -output ./kir_output -threads $((PARALLEL_SAMPLES * THREADS_PER_SAMPLE))
NCOPY_END=$(date +%s)
log_time "ncopy" "all" "$NCOPY_START" "$NCOPY_END"

# Step 4: genotype (batch operation)
echo "=== Running genotype ==="
GENO_START=$(date +%s)
kir-mapper genotype -output ./kir_output -threads $((PARALLEL_SAMPLES * THREADS_PER_SAMPLE))
GENO_END=$(date +%s)
log_time "genotype" "all" "$GENO_START" "$GENO_END"

# Step 5: Haplotype (only if NUM_SAMPLES >= 50)
if [ $NUM_SAMPLES -ge 50 ]; then
    echo "=== Running haplotype ==="
    HAPLO_START=$(date +%s)
    kir-mapper haplotype -output ./kir_output -threads $((PARALLEL_SAMPLES * THREADS_PER_SAMPLE))
    HAPLO_END=$(date +%s)
    log_time "haplotype" "all" "$HAPLO_START" "$HAPLO_END"
else
    echo "=== Skipping haplotype (need >=50 samples, have $NUM_SAMPLES) ==="
fi

echo "=== Output sizes ==="
du -sh ./kir_output/
du -sh ./kir_output/*/ 2>/dev/null || true
find ./kir_output -type f -name "*.txt" -o -name "*.tsv" -o -name "*.vcf*" | head -20 | xargs ls -lh

# Copy results
cp -r kir_output/* $OUTPUT_DIR/
cp $TIMING_LOG $OUTPUT_DIR/timing_log.csv

echo "=== Complete ==="
cat $TIMING_LOG