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
echo "Output directory: $OUTDIR"

autoskip() {
  local out="$1"
  if [ -s "$out" ]; then
    echo "  -> Skipping, output exists: $out"
    return 0
  fi
  return 1
}

# Create results file for tabular output
RESULTS_TABLE="$OUTDIR/benchmark_results.tsv"
echo -e "minlen\tloci_assessed\tpearson_r\tr_squared" > "$RESULTS_TABLE"

for MINLEN in $(seq "$START_MINLEN" "$END_MINLEN"); do
  OUT_FILE="$OUTDIR/calls_minlen_${MINLEN}.inq"
  BENCHMARK_FILE="$OUTDIR/benchmark_minlen_${MINLEN}.json"
  
  echo "[minlen=$MINLEN] Running inquiSTR call..."
  if autoskip "$OUT_FILE"; then
    # If call output exists, check if benchmark exists too
    if [ ! -f "$BENCHMARK_FILE" ]; then
      echo "  -> Call output exists, but running benchmark..."
    else
      echo "  -> Both outputs exist, skipping"
      continue
    fi
  fi
  
  if [ "${VERBOSE:-0}" != "0" ]; then
    echo "+ $INQUISTR_BIN call \"$INPUT_CRAM\" \\
      -R \"$REGION_FILE\" \\
      --reference \"$REFERENCE\" \\
      --max-locus \"$MAX_LOCUS\" \\
      --threads \"$THREADS\" \\
      --minlen \"$MINLEN\" > \"$OUT_FILE\""
  fi
  
  # Run inquiSTR call
  "$INQUISTR_BIN" call "$INPUT_CRAM" \
    -R "$REGION_FILE" \
    --reference "$REFERENCE" \
    --max-locus "$MAX_LOCUS" \
    --threads "$THREADS" \
    --minlen "$MINLEN" \
    > "$OUT_FILE"
  echo "  -> Wrote: $OUT_FILE"
  
  # Run inquiSTR benchmark
  echo "[minlen=$MINLEN] Running inquiSTR benchmark..."
  PLOT_FILE="$OUTDIR/plot_minlen_${MINLEN}.html"
  if [ "${VERBOSE:-0}" != "0" ]; then
    echo "+ $INQUISTR_BIN benchmark \"$OUT_FILE\" \\
      --bed \"$REGION_FILE\" \\
      --max-locus \"$MAX_LOCUS\" \\
      --plot \"$PLOT_FILE\" > \"$BENCHMARK_FILE\""
  fi
  
  "$INQUISTR_BIN" benchmark "$OUT_FILE" \
    --bed "$REGION_FILE" \
    --max-locus "$MAX_LOCUS" \
    --plot "$PLOT_FILE" \
    > "$BENCHMARK_FILE"
  echo "  -> Wrote: $BENCHMARK_FILE and $PLOT_FILE"
  
  # Extract metrics from benchmark output and append to results table
  if [ -f "$BENCHMARK_FILE" ]; then
    LOCI=$(grep -oP 'LOCI_ASSESSED:\s*\K[0-9]+' "$BENCHMARK_FILE" || echo "N/A")
    PEARSON_R=$(grep -oP 'PEARSON_R:\s*\K[0-9.-]+' "$BENCHMARK_FILE" || echo "N/A")
    R_SQUARED=$(grep -oP 'R_SQUARED:\s*\K[0-9.]+' "$BENCHMARK_FILE" || echo "N/A")
    echo -e "$MINLEN\t$LOCI\t$PEARSON_R\t$R_SQUARED" >> "$RESULTS_TABLE"
  fi
done

echo ""
echo "All runs completed. Results in: $OUTDIR"
echo ""
echo "Benchmark Summary:"
column -t -s $'\t' "$RESULTS_TABLE"