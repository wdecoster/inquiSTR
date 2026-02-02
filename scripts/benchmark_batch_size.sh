#!/usr/bin/env bash
#
# Benchmark inquiSTR performance across different thread counts and batch sizes
# to determine optimal parameters for different hardware configurations.
#
# Usage: ./benchmark_batch_size.sh <cram_file> <bed_file> <reference> [output_dir]

set -euo pipefail

# Configuration
CRAM_FILE="${1:?Usage: $0 <cram_file> <bed_file> <reference> [output_dir]}"
BED_FILE="${2:?Usage: $0 <cram_file> <bed_file> <reference> [output_dir]}"
REFERENCE="${3:?Usage: $0 <cram_file> <bed_file> <reference> [output_dir]}"
OUTPUT_DIR="${4:-benchmark_results}"
BINARY="./target/x86_64-unknown-linux-musl/release/inquiSTR"

# Test parameters - adjust these based on your system
THREAD_COUNTS=(1 2 4 8 16)
BATCH_SIZES=(5 10 20 30 50 75 100)

# Number of repeats per combination (for statistical stability)
REPEATS=3

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Create output directory
mkdir -p "$OUTPUT_DIR"
RESULTS_FILE="$OUTPUT_DIR/benchmark_results.tsv"
LOG_FILE="$OUTPUT_DIR/benchmark.log"

# Initialize results file
echo -e "threads\tbatch_size\trepeat\twall_time_s\tuser_time_s\tsys_time_s\tcpu_percent\tmem_max_kb\tvol_ctx_switches\tinvol_ctx_switches\tfs_inputs" > "$RESULTS_FILE"

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}inquiSTR Batch Size Benchmark${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Configuration:"
echo "  CRAM:      $CRAM_FILE"
echo "  BED:       $BED_FILE"
echo "  Reference: $REFERENCE"
echo "  Binary:    $BINARY"
echo "  Output:    $OUTPUT_DIR"
echo ""
echo "Test matrix:"
echo "  Threads:     ${THREAD_COUNTS[*]}"
echo "  Batch sizes: ${BATCH_SIZES[*]} KB"
echo "  Repeats:     $REPEATS"
echo ""

# Calculate total tests
TOTAL_TESTS=$((${#THREAD_COUNTS[@]} * ${#BATCH_SIZES[@]} * REPEATS))
CURRENT_TEST=0

START_TIME=$(date +%s)

for threads in "${THREAD_COUNTS[@]}"; do
    for batch_size in "${BATCH_SIZES[@]}"; do
        for repeat in $(seq 1 $REPEATS); do
            CURRENT_TEST=$((CURRENT_TEST + 1))
            
            echo -e "${BLUE}[${CURRENT_TEST}/${TOTAL_TESTS}]${NC} Testing: threads=${threads}, batch_size=${batch_size}KB, repeat=${repeat}/${REPEATS}"
            
            # Run with /usr/bin/time to capture detailed metrics
            OUTPUT_FILE="$OUTPUT_DIR/time_t${threads}_b${batch_size}_r${repeat}.txt"
            
            # Clear page cache (requires sudo) - comment out if not available
            # sudo sh -c 'echo 3 > /proc/sys/vm/drop_caches' 2>/dev/null || true
            
            if /usr/bin/time -v -o "$OUTPUT_FILE" \
                "$BINARY" call "$CRAM_FILE" \
                -R "$BED_FILE" \
                --reference "$REFERENCE" \
                --threads "$threads" \
                --batch-size "$batch_size" \
                > "$OUTPUT_DIR/output_t${threads}_b${batch_size}_r${repeat}.tsv" 2>&1; then
                
                # Parse the time output
                wall_time=$(grep "Elapsed (wall clock)" "$OUTPUT_FILE" | awk '{print $NF}' | awk -F: '{if (NF==3) print ($1*3600)+($2*60)+$3; else if (NF==2) print ($1*60)+$2; else print $1}')
                user_time=$(grep "User time" "$OUTPUT_FILE" | awk '{print $4}')
                sys_time=$(grep "System time" "$OUTPUT_FILE" | awk '{print $4}')
                cpu_percent=$(grep "Percent of CPU" "$OUTPUT_FILE" | awk '{print $7}' | tr -d '%')
                mem_max=$(grep "Maximum resident" "$OUTPUT_FILE" | awk '{print $6}')
                vol_ctx=$(grep "Voluntary context switches" "$OUTPUT_FILE" | awk '{print $5}')
                invol_ctx=$(grep "Involuntary context switches" "$OUTPUT_FILE" | awk '{print $5}')
                fs_inputs=$(grep "File system inputs" "$OUTPUT_FILE" | awk '{print $4}')
                
                # Write to results file
                echo -e "${threads}\t${batch_size}\t${repeat}\t${wall_time}\t${user_time}\t${sys_time}\t${cpu_percent}\t${mem_max}\t${vol_ctx}\t${invol_ctx}\t${fs_inputs}" >> "$RESULTS_FILE"
                
                echo -e "  ${GREEN}✓${NC} Wall time: ${wall_time}s, CPU: ${cpu_percent}%, Memory: ${mem_max}KB"
            else
                echo -e "  ${RED}✗ FAILED${NC}" | tee -a "$LOG_FILE"
                echo -e "${threads}\t${batch_size}\t${repeat}\tFAILED\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> "$RESULTS_FILE"
            fi
            
            # Brief pause between tests
            sleep 1
        done
    done
done

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Benchmark Complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo "Total time: ${ELAPSED}s ($((ELAPSED / 60)) minutes)"
echo "Results saved to: $RESULTS_FILE"
echo ""
echo "Next steps:"
echo "  1. Analyze results: python3 scripts/analyze_benchmark_results.py $RESULTS_FILE"
echo "  2. View detailed logs in: $OUTPUT_DIR"
