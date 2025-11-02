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

# Create results file for tabular output
RESULTS_TABLE="$OUTDIR/benchmark_results.tsv"
echo -e "minlen\tloci_assessed\tzero_pairs_excluded\tr\tr2" > "$RESULTS_TABLE"

for MINLEN in $(seq "$START_MINLEN" "$END_MINLEN"); do
  OUT_FILE="$OUTDIR/calls_minlen_${MINLEN}.inq"
  BENCHMARK_FILE="$OUTDIR/benchmark_minlen_${MINLEN}.txt"
  
  # Check if both outputs already exist
  if [ -s "$OUT_FILE" ] && [ -s "$BENCHMARK_FILE" ]; then
    echo "[minlen=$MINLEN] Both call and benchmark outputs exist, skipping"
    continue
  fi
  
  # Run inquiSTR call if needed
  if [ ! -s "$OUT_FILE" ]; then
    echo "[minlen=$MINLEN] Running inquiSTR call..."
    if [ "${VERBOSE:-0}" != "0" ]; then
      echo "+ $INQUISTR_BIN call \"$INPUT_CRAM\" \\
        -R \"$REGION_FILE\" \\
        --reference \"$REFERENCE\" \\
        --max-locus \"$MAX_LOCUS\" \\
        --threads \"$THREADS\" \\
        --minlen \"$MINLEN\" > \"$OUT_FILE\""
    fi
    
    "$INQUISTR_BIN" call "$INPUT_CRAM" \
      -R "$REGION_FILE" \
      --reference "$REFERENCE" \
      --max-locus "$MAX_LOCUS" \
      --threads "$THREADS" \
      --minlen "$MINLEN" \
      > "$OUT_FILE"
    echo "  -> Wrote: $OUT_FILE"
  else
    echo "[minlen=$MINLEN] Call output exists: $OUT_FILE"
  fi
  
  # Run inquiSTR benchmark if needed
  PLOT_FILE="$OUTDIR/plot_minlen_${MINLEN}.html"
  DIFF_FILE="$OUTDIR/diff_minlen_${MINLEN}.tsv"
  if [ ! -s "$BENCHMARK_FILE" ]; then
    echo "[minlen=$MINLEN] Running inquiSTR benchmark..."
    if [ "${VERBOSE:-0}" != "0" ]; then
      echo "+ $INQUISTR_BIN benchmark \"$OUT_FILE\" \\
        --bed \"$REGION_FILE\" \\
        --max-locus \"$MAX_LOCUS\" \\
        --plot \"$PLOT_FILE\" \\
        --diff-out \"$DIFF_FILE\" \\
        --tier1 \\
        --nonzero > \"$BENCHMARK_FILE\""
    fi
    
    "$INQUISTR_BIN" benchmark "$OUT_FILE" \
      --bed "$REGION_FILE" \
      --max-locus "$MAX_LOCUS" \
      --plot "$PLOT_FILE" \
      --diff-out "$DIFF_FILE" \
      --tier1 \
      --nonzero \
      > "$BENCHMARK_FILE"
    echo "  -> Wrote: $BENCHMARK_FILE, $PLOT_FILE, and $DIFF_FILE"
  else
    echo "[minlen=$MINLEN] Benchmark output exists: $BENCHMARK_FILE"
  fi
  
  # Extract metrics from benchmark output (nonzero mode format)
  if [ -f "$BENCHMARK_FILE" ]; then
    LOCI=$(grep -oP 'LOCI_ASSESSED:\s*\K[0-9]+' "$BENCHMARK_FILE" || echo "N/A")
    ZERO_EXCLUDED=$(grep -oP 'ZERO_ZERO_PAIRS_EXCLUDED:\s*\K[0-9]+' "$BENCHMARK_FILE" || echo "N/A")
    R=$(grep -oP 'PEARSON_R:\s*\K[0-9.-]+' "$BENCHMARK_FILE" || echo "N/A")
    R2=$(grep -oP 'R_SQUARED:\s*\K[0-9.]+' "$BENCHMARK_FILE" || echo "N/A")
    echo -e "$MINLEN\t$LOCI\t$ZERO_EXCLUDED\t$R\t$R2" >> "$RESULTS_TABLE"
  fi
done

echo ""
echo "All runs completed. Results in: $OUTDIR"
echo ""
echo "Benchmark Summary:"
column -t -s $'\t' "$RESULTS_TABLE"