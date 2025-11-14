#!/usr/bin/env bash
# Run inquiSTR call sweeping --minlen from 1 to 10
# Fixed parameters:
#   --max-locus 1000
#   --threads 4
#   -R /home/wdecoster/bubblebath-backup/wdecoster/optimize_inquiSTR/adotto/HG002_GRCh38_TandemRepeats_v1.0.bed.gz
#   --reference ~/GRCh38.fa
#   <BAM/CRAM> /home/wdecoster/bubblebath-backup/wdecoster/optimize_inquiSTR/benchmark/ont.cram

set -euo pipefail

# Configure paths (override via environment if desired)
INQUISTR_BIN=${INQUISTR_BIN:-inquiSTR}
REGION_FILE=${REGION_FILE:-/home/wdecoster/bubblebath-backup/wdecoster/optimize_inquiSTR/adotto/HG002_GRCh38_TandemRepeats_v1.0.bed.gz}
REFERENCE=${REFERENCE:-$HOME/GRCh38.fa}
INPUT_CRAM=${INPUT_CRAM:-/home/wdecoster/bubblebath-backup/wdecoster/optimize_inquiSTR/benchmark/ont.cram}
THREADS=${THREADS:-4}
MAX_LOCUS=${MAX_LOCUS:-1000}
OUTDIR=${OUTDIR:-results_minlen}
START_MINLEN=${START_MINLEN:-1}
END_MINLEN=${END_MINLEN:-10}
REPLICATES=${REPLICATES:-3}

# Basic checks
if ! command -v "$INQUISTR_BIN" >/dev/null 2>&1; then
  echo "Error: '$INQUISTR_BIN' not found in PATH. Set INQUISTR_BIN to your inquiSTR binary." >&2
  exit 1
fi

if [ ! -f "$REFERENCE" ]; then
  echo "Error: reference fasta not found at '$REFERENCE'" >&2
  exit 1
fi

if [ ! -f "$INPUT_CRAM" ]; then
  echo "Error: input CRAM/BAM not found at '$INPUT_CRAM'" >&2
  exit 1
fi

mkdir -p "$OUTDIR"

echo "Running inquiSTR call sweep: minlen ${START_MINLEN}..${END_MINLEN}"
echo "Replicates per minlen: ${REPLICATES}"
echo "Output directory: $OUTDIR"

# Create results file for tabular output
RESULTS_TABLE="$OUTDIR/results_summary.tsv"
echo -e "minlen\tloci\tzero_excluded\tr\tr2\texact_matches\texact_pct\twithin_tol\twithin_tol_pct\tcall_time_mean\tcall_time_stdev" > "$RESULTS_TABLE"

# Create array of all minlen-replicate combinations and shuffle them
COMBINATIONS=()
for MINLEN in $(seq "$START_MINLEN" "$END_MINLEN"); do
  for REP in $(seq 1 "$REPLICATES"); do
    COMBINATIONS+=("${MINLEN}:${REP}")
  done
done

# Shuffle the combinations using shuf if available, otherwise use a simple sort -R
if command -v shuf >/dev/null 2>&1; then
  mapfile -t SHUFFLED < <(printf '%s\n' "${COMBINATIONS[@]}" | shuf)
else
  mapfile -t SHUFFLED < <(printf '%s\n' "${COMBINATIONS[@]}" | sort -R)
fi

# Create associative arrays to store timing results
declare -A CALL_TIMES_BY_MINLEN

echo "Running ${#SHUFFLED[@]} timing measurements in randomized order..."
echo ""

# Run all timing measurements in randomized order
for i in "${!SHUFFLED[@]}"; do
  COMBO="${SHUFFLED[$i]}"
  MINLEN="${COMBO%%:*}"
  REP="${COMBO##*:}"
  
  echo "[${i}/${#SHUFFLED[@]}] minlen=${MINLEN}, replicate=${REP}"
  
  # Time the call
  CALL_START=$(date +%s)
  "$INQUISTR_BIN" call "$INPUT_CRAM" \
    -R "$REGION_FILE" \
    --reference "$REFERENCE" \
    --max-locus "$MAX_LOCUS" \
    --threads "$THREADS" \
    --minlen "$MINLEN" \
    > /dev/null
  CALL_END=$(date +%s)
  CALL_TIME=$((CALL_END - CALL_START))
  
  # Store timing (append to existing values for this minlen)
  if [ -z "${CALL_TIMES_BY_MINLEN[$MINLEN]+x}" ]; then
    CALL_TIMES_BY_MINLEN[$MINLEN]="$CALL_TIME"
  else
    CALL_TIMES_BY_MINLEN[$MINLEN]="${CALL_TIMES_BY_MINLEN[$MINLEN]} $CALL_TIME"
  fi
  
  echo "  -> ${CALL_TIME}s"
done

echo ""
echo "Timing measurements complete. Now generating output files and running benchmarks..."
echo ""

# Now generate the actual output files and run benchmarks (once per minlen)
for MINLEN in $(seq "$START_MINLEN" "$END_MINLEN"); do
  OUT_FILE="$OUTDIR/calls_minlen_${MINLEN}.inq"
  BENCHMARK_FILE="$OUTDIR/benchmark_minlen_${MINLEN}.txt"
  PLOT_FILE="$OUTDIR/plot_minlen_${MINLEN}.html"
  DIFF_FILE="$OUTDIR/diff_minlen_${MINLEN}.tsv"
  
  echo "[minlen=$MINLEN] Generating output file..."
  
  # Generate the output file (just once, no timing)
  if [ ! -s "$OUT_FILE" ]; then
    "$INQUISTR_BIN" call "$INPUT_CRAM" \
      -R "$REGION_FILE" \
      --reference "$REFERENCE" \
      --max-locus "$MAX_LOCUS" \
      --threads "$THREADS" \
      --minlen "$MINLEN" \
      > "$OUT_FILE"
    echo "  -> Wrote: $OUT_FILE"
  else
    echo "  -> Output exists: $OUT_FILE"
  fi
  
  # Run benchmark (just once, no timing)
  if [ ! -s "$BENCHMARK_FILE" ]; then
    echo "[minlen=$MINLEN] Running benchmark..."
    "$INQUISTR_BIN" benchmark "$OUT_FILE" \
      --bed "$REGION_FILE" \
      --max-locus "$MAX_LOCUS" \
      --plot "$PLOT_FILE" \
      --diff-out "$DIFF_FILE" \
      --tier1 \
      --tolerance 3 \
      --nonzero \
      > "$BENCHMARK_FILE"
    echo "  -> Wrote: $BENCHMARK_FILE, $PLOT_FILE, and $DIFF_FILE"
  else
    echo "  -> Benchmark exists: $BENCHMARK_FILE"
  fi
  
  # Calculate mean and stdev from the timing data collected earlier
  CALL_TIMES_STR="${CALL_TIMES_BY_MINLEN[$MINLEN]:-}"
  if [ -n "$CALL_TIMES_STR" ]; then
    # Convert space-separated string to array
    read -ra CALL_TIMES_ARRAY <<< "$CALL_TIMES_STR"
    
    # Calculate mean
    CALL_MEAN=$(awk 'BEGIN{s=0; for(i=1; i<ARGC; i++) s+=ARGV[i]; print s/(ARGC-1)}' "${CALL_TIMES_ARRAY[@]}")
    
    # Calculate stdev
    CALL_STDEV=$(awk -v mean="$CALL_MEAN" 'BEGIN{s=0; for(i=1; i<ARGC; i++) s+=(ARGV[i]-mean)^2; print sqrt(s/(ARGC-1))}' "${CALL_TIMES_ARRAY[@]}")
    
    CALL_MEAN=$(printf "%.1f" "$CALL_MEAN")
    CALL_STDEV=$(printf "%.1f" "$CALL_STDEV")
  else
    CALL_MEAN="N/A"
    CALL_STDEV="N/A"
  fi
  
  # Extract metrics from benchmark output (nonzero mode format)
  if [ -f "$BENCHMARK_FILE" ]; then
    LOCI=$(grep -oP 'LOCI_ASSESSED:\s*\K[0-9]+' "$BENCHMARK_FILE" || echo "N/A")
    ZERO_EXCLUDED=$(grep -oP 'ZERO_ZERO_PAIRS_EXCLUDED:\s*\K[0-9]+' "$BENCHMARK_FILE" || echo "N/A")
    R=$(grep -oP 'PEARSON_R:\s*\K[0-9.-]+' "$BENCHMARK_FILE" || echo "N/A")
    R2=$(grep -oP 'R_SQUARED:\s*\K[0-9.]+' "$BENCHMARK_FILE" || echo "N/A")
    EXACT=$(grep -oP 'EXACT_MATCHES:\s*\K[0-9]+' "$BENCHMARK_FILE" || echo "N/A")
    EXACT_PCT=$(grep -oP 'EXACT_MATCHES_PERCENT:\s*\K[0-9.]+' "$BENCHMARK_FILE" || echo "N/A")
    WITHIN_TOL=$(grep -oP 'WITHIN_TOLERANCE:\s*\K[0-9]+' "$BENCHMARK_FILE" || echo "N/A")
    WITHIN_TOL_PCT=$(grep -oP 'WITHIN_TOLERANCE_PERCENT:\s*\K[0-9.]+' "$BENCHMARK_FILE" || echo "N/A")
    echo -e "$MINLEN\t$LOCI\t$ZERO_EXCLUDED\t$R\t$R2\t$EXACT\t$EXACT_PCT\t$WITHIN_TOL\t$WITHIN_TOL_PCT\t$CALL_MEAN\t$CALL_STDEV" >> "$RESULTS_TABLE"
  fi
done

echo ""
echo "All runs completed. Results in: $OUTDIR"
echo ""
echo "Benchmark Summary:"
column -t -s $'\t' "$RESULTS_TABLE"