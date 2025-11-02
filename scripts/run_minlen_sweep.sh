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
echo -e "minlen\tloci_assessed\tzero_pairs\tnonzero_loci\tr_all\tr2_all\tr_nonzero\tr2_nonzero" > "$RESULTS_TABLE"

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
  if [ ! -s "$BENCHMARK_FILE" ]; then
    echo "[minlen=$MINLEN] Running inquiSTR benchmark..."
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
  else
    echo "[minlen=$MINLEN] Benchmark output exists: $BENCHMARK_FILE"
  fi
  
  # Extract metrics from benchmark output and append to results table
  if [ -f "$BENCHMARK_FILE" ]; then
    LOCI=$(grep -oP 'LOCI_ASSESSED:\s*\K[0-9]+' "$BENCHMARK_FILE" || echo "N/A")
    ZERO_PAIRS=$(grep -oP 'ZERO_ZERO_PAIRS:\s*\K[0-9]+' "$BENCHMARK_FILE" || echo "N/A")
    NONZERO=$(grep -oP 'NONZERO_LOCI:\s*\K[0-9]+' "$BENCHMARK_FILE" || echo "N/A")
    R_ALL=$(grep -oP 'PEARSON_R_ALL:\s*\K[0-9.-]+' "$BENCHMARK_FILE" || echo "N/A")
    R2_ALL=$(grep -oP 'R_SQUARED_ALL:\s*\K[0-9.]+' "$BENCHMARK_FILE" || echo "N/A")
    R_NONZERO=$(grep -oP 'PEARSON_R_NONZERO:\s*\K[0-9.-]+' "$BENCHMARK_FILE" || echo "N/A")
    R2_NONZERO=$(grep -oP 'R_SQUARED_NONZERO:\s*\K[0-9.]+' "$BENCHMARK_FILE" || echo "N/A")
    echo -e "$MINLEN\t$LOCI\t$ZERO_PAIRS\t$NONZERO\t$R_ALL\t$R2_ALL\t$R_NONZERO\t$R2_NONZERO" >> "$RESULTS_TABLE"
  fi
done

echo ""
echo "All runs completed. Results in: $OUTDIR"
echo ""
echo "Benchmark Summary:"
column -t -s $'\t' "$RESULTS_TABLE"