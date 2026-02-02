#!/usr/bin/env bash
#
# Quick auto-tuning script that users can run to find optimal batch_size
# for their specific hardware and dataset.
#
# This runs a reduced test matrix (faster than full benchmark) and provides
# immediate recommendations.
#
# Usage: ./auto_tune_batch_size.sh <cram_file> <bed_file> <reference> [threads]

set -euo pipefail

CRAM_FILE="${1:?Usage: $0 <cram_file> <bed_file> <reference> [threads]}"
BED_FILE="${2:?Usage: $0 <cram_file> <bed_file> <reference> [threads]}"
REFERENCE="${3:?Usage: $0 <cram_file> <bed_file> <reference> [threads]}"
THREADS="${4:-$(nproc)}"
BINARY="./target/x86_64-unknown-linux-musl/release/inquiSTR"

# Test fewer batch sizes for quick tuning
BATCH_SIZES=(10 20 30 50)

echo "========================================"
echo "inquiSTR Auto-Tuning"
echo "========================================"
echo ""
echo "System: $(nproc) cores available"
echo "Testing with: $THREADS threads"
echo "Batch sizes to test: ${BATCH_SIZES[*]} KB"
echo ""

BEST_TIME=999999
BEST_BATCH=50

for batch_size in "${BATCH_SIZES[@]}"; do
    echo "Testing batch_size=${batch_size}KB..."
    
    OUTPUT=$(/usr/bin/time -f "Wall: %e, CPU: %P, Sys: %S" \
        "$BINARY" call "$CRAM_FILE" \
        -R "$BED_FILE" \
        --reference "$REFERENCE" \
        --threads "$THREADS" \
        --batch-size "$batch_size" \
        > /dev/null 2>&1)
    
    WALL_TIME=$(echo "$OUTPUT" | grep -oP 'Wall: \K[0-9.]+')
    CPU_PCT=$(echo "$OUTPUT" | grep -oP 'CPU: \K[0-9]+')
    SYS_TIME=$(echo "$OUTPUT" | grep -oP 'Sys: \K[0-9.]+')
    
    echo "  → ${WALL_TIME}s (CPU: ${CPU_PCT}%, System: ${SYS_TIME}s)"
    
    if (( $(echo "$WALL_TIME < $BEST_TIME" | bc -l) )); then
        BEST_TIME=$WALL_TIME
        BEST_BATCH=$batch_size
    fi
done

echo ""
echo "========================================"
echo "RECOMMENDATION"
echo "========================================"
echo ""
echo "  Optimal batch_size: ${BEST_BATCH}KB"
echo "  Runtime: ${BEST_TIME}s"
echo ""
echo "To use this setting, add to your command:"
echo "  --batch-size ${BEST_BATCH}"
echo ""
