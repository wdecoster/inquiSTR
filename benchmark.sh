#!/bin/bash
# Complete inquiSTR benchmarking script with optional Valgrind memory profiling
# Usage: ./benchmark.sh [--skip-valgrind] [--help]
# Environment: SKIP_VALGRIND=1 to disable Valgrind profiling

set -e  # Exit on any error

# Parse command line arguments
SKIP_VALGRIND_ARG=false
while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-valgrind)
            SKIP_VALGRIND_ARG=true
            shift
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "OPTIONS:"
            echo "  --skip-valgrind    Disable Valgrind memory profiling (faster execution)"
            echo "  --help, -h         Show this help message"
            echo ""
            echo "ENVIRONMENT VARIABLES:"
            echo "  SKIP_VALGRIND=1    Alternative way to disable Valgrind profiling"
            echo ""
            echo "Examples:"
            echo "  $0                           # Full benchmark with Valgrind"
            echo "  $0 --skip-valgrind          # Fast benchmark, no Valgrind"
            echo "  SKIP_VALGRIND=1 $0          # Fast benchmark using env var"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Determine if Valgrind should be skipped (command line arg takes precedence)
if [[ "$SKIP_VALGRIND_ARG" == true || "$SKIP_VALGRIND" == "1" ]]; then
    SKIP_VALGRIND_ENABLED=true
    echo "=== inquiSTR Fast Benchmarking Script ==="
    echo "Comparing v0.14.0 vs v0.15.0 (Valgrind disabled for speed)"
else
    SKIP_VALGRIND_ENABLED=false
    echo "=== inquiSTR Performance Benchmarking Script ==="
    echo "Comparing v0.14.0 vs v0.15.0 with Valgrind memory profiling"
fi

echo "Started at: $(date)"
echo ""

# Function to convert time to seconds for calculations
time_to_seconds() {
    local time_str="$1"
    if [[ $time_str == *:*:* ]]; then
        # Format: h:mm:ss
        echo "$time_str" | awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3 }'
    elif [[ $time_str == *:* ]]; then
        # Format: m:ss
        echo "$time_str" | awk -F: '{ print ($1 * 60) + $2 }'
    else
        # Already in seconds
        echo "$time_str"
    fi
}

# Function to calculate speedup ratio
calculate_speedup() {
    local old_time="$1"
    local new_time="$2"
    
    if [[ "$BC_AVAILABLE" != true ]]; then
        echo "N/A (bc not available)"
        return
    fi
    
    local old_sec
    local new_sec
    old_sec=$(time_to_seconds "$old_time")
    new_sec=$(time_to_seconds "$new_time")
    echo "scale=2; $old_sec / $new_sec" | bc -l
}

# Set up binary paths
INQUISTR_014_PATH="./inquiSTR_014"
INQUISTR_015_PATH="./inquiSTR_015"

# Check if required files exist
if [[ ! -f "$INQUISTR_014_PATH" ]]; then
    echo "ERROR: inquiSTR_014 binary not found!"
    exit 1
fi

if [[ ! -f "$INQUISTR_015_PATH" ]]; then
    echo "ERROR: inquiSTR_015 binary not found!"
    exit 1
fi

if [[ ! -f "alignment.cram" ]]; then
    echo "ERROR: alignment.cram not found!"
    exit 1
fi

if [[ ! -f "benchmark_repeats.bed" ]]; then
    echo "ERROR: benchmark_repeats.bed not found!"
    exit 1
fi

# Check if valgrind should be used (skip if disabled or not available)
if [[ "$SKIP_VALGRIND_ENABLED" == true ]]; then
    echo "⚡ Valgrind disabled for faster execution"
    VALGRIND_AVAILABLE=false
elif ! command -v valgrind &> /dev/null; then
    echo "WARNING: Valgrind not found! Memory profiling will be skipped."
    echo "To install: sudo apt-get install valgrind (Ubuntu/Debian) or brew install valgrind (macOS)"
    VALGRIND_AVAILABLE=false
else
    echo "✓ Valgrind found - memory profiling enabled"
    VALGRIND_AVAILABLE=true
fi

# Check if bc is available for calculations
if ! command -v bc &> /dev/null; then
    echo "WARNING: bc calculator not found! Performance analysis will be limited."
    echo "To install: sudo apt-get install bc (Ubuntu/Debian) or brew install bc (macOS)"
    BC_AVAILABLE=false
else
    BC_AVAILABLE=true
fi

echo "✓ All required files found"
echo "Bed file regions: $(wc -l < benchmark_repeats.bed)"
echo ""

# Log binary versions for tracking
echo "=== BINARY VERSION INFORMATION ==="
echo "inquiSTR v0.14.0 version details:"
"$INQUISTR_014_PATH" --version 2>/dev/null || echo "  Version information not available"

echo ""
echo "inquiSTR v0.15.x version details:"  
"$INQUISTR_015_PATH" --version 2>/dev/null || echo "  Version information not available"

echo ""
echo "Build information:"
echo "  inquiSTR_014 size: $(ls -lhL "$INQUISTR_014_PATH" | awk '{print $5}')"
echo "  inquiSTR_015 size: $(ls -lhL "$INQUISTR_015_PATH" | awk '{print $5}')"
echo "  inquiSTR_014 modified: $(stat -L -c %y "$INQUISTR_014_PATH" 2>/dev/null | cut -d. -f1 || echo "unknown")"
echo "  inquiSTR_015 modified: $(stat -L -c %y "$INQUISTR_015_PATH" 2>/dev/null | cut -d. -f1 || echo "unknown")"
echo ""

# Clean previous results
echo "Cleaning previous benchmark results..."
rm -f output_014*.tsv output_015*.tsv benchmark_*.log memory_*.txt massif_*.out valgrind_*.log
echo ""

# =============================================================================
# TIMING BENCHMARKS
# =============================================================================
echo "=== TIMING BENCHMARKS ==="

echo "📊 Testing v0.14.0 (timing)..."
/usr/bin/time -v "$INQUISTR_014_PATH" call alignment.cram --region-file benchmark_repeats.bed --threads 1 > output_014_time.tsv 2> benchmark_014_time.log
echo "✓ v0.14.0 timing complete"

echo "📊 Testing v0.15.0 (timing)..."
/usr/bin/time -v "$INQUISTR_015_PATH" call alignment.cram --region-file benchmark_repeats.bed --threads 1 > output_015_time.tsv 2> benchmark_015_time.log
echo "✓ v0.15.0 timing complete"

echo ""

# =============================================================================
# VALGRIND MEMORY PROFILING (if available and enabled)
# =============================================================================
if [[ "$VALGRIND_AVAILABLE" == true ]]; then
    echo "=== VALGRIND MEMORY PROFILING ==="
    
    echo "🧠 Running v0.14.0 with Valgrind Massif (memory profiling)..."
    echo "   This may take significantly longer than normal execution..."
    valgrind --tool=massif \
             --massif-out-file=massif_014.out \
             --time-unit=ms \
             --detailed-freq=1 \
             --max-snapshots=200 \
             "$INQUISTR_014_PATH" call alignment.cram --region-file benchmark_repeats.bed --threads 1 \
             > output_014_valgrind.tsv 2> valgrind_014.log
    
    # Generate memory report for v0.14.0
    ms_print massif_014.out > memory_report_014.txt
    echo "✓ v0.14.0 Valgrind profiling complete"
    
    echo "🧠 Running v0.15.0 with Valgrind Massif (memory profiling)..."
    echo "   This may take significantly longer than normal execution..."
    valgrind --tool=massif \
             --massif-out-file=massif_015.out \
             --time-unit=ms \
             --detailed-freq=1 \
             --max-snapshots=200 \
             "$INQUISTR_015_PATH" call alignment.cram --region-file benchmark_repeats.bed --threads 1 \
             > output_015_valgrind.tsv 2> valgrind_015.log
    
    # Generate memory report for v0.15.0
    ms_print massif_015.out > memory_report_015.txt
    echo "✓ v0.15.0 Valgrind profiling complete"
    
    echo ""
    
    # Additional Valgrind checks for memory leaks
    echo "🔍 Running Valgrind Memcheck for memory leaks..."
    
    echo "   Checking v0.14.0 for memory leaks..."
    valgrind --tool=memcheck \
             --leak-check=full \
             --show-leak-kinds=all \
             --track-origins=yes \
             --log-file=memcheck_014.log \
             "$INQUISTR_014_PATH" call alignment.cram --region-file benchmark_repeats.bed --threads 1 \
             > /dev/null 2>&1
    
    echo "   Checking v0.15.0 for memory leaks..."
    valgrind --tool=memcheck \
             --leak-check=full \
             --show-leak-kinds=all \
             --track-origins=yes \
             --log-file=memcheck_015.log \
             "$INQUISTR_015_PATH" call alignment.cram --region-file benchmark_repeats.bed --threads 1 \
             > /dev/null 2>&1
    
    echo "✓ Memory leak analysis complete"
    echo ""
else
    echo "=== FAST BENCHMARKING MODE ==="
    echo "⚡ Valgrind disabled - running faster timing-only benchmarks"
    echo ""
fi

# =============================================================================
# MULTITHREADED PERFORMANCE TEST
# =============================================================================
echo "=== MULTITHREADED PERFORMANCE TEST ==="

for threads in 2 4 8; do
    if [[ $(nproc) -ge $threads ]]; then
        echo "🔄 Testing with $threads threads..."
        
        echo "   v0.14.0 ($threads threads)..."
        /usr/bin/time -v "$INQUISTR_014_PATH" call alignment.cram --region-file benchmark_repeats.bed --threads $threads > output_014_t${threads}.tsv 2> benchmark_014_t${threads}.log
        
        echo "   v0.15.0 ($threads threads)..."
        /usr/bin/time -v "$INQUISTR_015_PATH" call alignment.cram --region-file benchmark_repeats.bed --threads $threads > output_015_t${threads}.tsv 2> benchmark_015_t${threads}.log
        
        echo "✓ $threads thread test complete"
    else
        echo "⏭️  Skipping $threads threads (only $(nproc) CPUs available)"
    fi
done

echo ""

# =============================================================================
# RESULTS ANALYSIS
# =============================================================================
echo "=== BENCHMARK RESULTS ANALYSIS ==="

# Create benchmark summary with version information
echo "🔬 BENCHMARK ANALYSIS REPORT" > benchmark_summary.txt
echo "============================" >> benchmark_summary.txt
echo "Generated on: $(date)" >> benchmark_summary.txt
echo "" >> benchmark_summary.txt

# Version information for inquiSTR v0.14.0
echo "📦 BINARY VERSION INFORMATION:" >> benchmark_summary.txt
echo "==============================" >> benchmark_summary.txt
echo "" >> benchmark_summary.txt
echo "inquiSTR v0.14.0:" >> benchmark_summary.txt
echo "  Location: $INQUISTR_014_PATH" >> benchmark_summary.txt
echo "  File size: $(stat -L -f%z "$INQUISTR_014_PATH" 2>/dev/null || stat -L -c%s "$INQUISTR_014_PATH")B" >> benchmark_summary.txt
echo "  Modified: $(stat -L -f%Sm -t%Y-%m-%d\ %H:%M:%S "$INQUISTR_014_PATH" 2>/dev/null || stat -L -c%y "$INQUISTR_014_PATH")" >> benchmark_summary.txt
echo "  Version output:" >> benchmark_summary.txt
"$INQUISTR_014_PATH" --version 2>&1 | sed 's/^/    /' >> benchmark_summary.txt
echo "" >> benchmark_summary.txt

echo "inquiSTR v0.15.0:" >> benchmark_summary.txt
echo "  Location: $INQUISTR_015_PATH" >> benchmark_summary.txt
echo "  File size: $(stat -L -f%z "$INQUISTR_015_PATH" 2>/dev/null || stat -L -c%s "$INQUISTR_015_PATH")B" >> benchmark_summary.txt
echo "  Modified: $(stat -L -f%Sm -t%Y-%m-%d\ %H:%M:%S "$INQUISTR_015_PATH" 2>/dev/null || stat -L -c%y "$INQUISTR_015_PATH")" >> benchmark_summary.txt
echo "  Version output:" >> benchmark_summary.txt
"$INQUISTR_015_PATH" --version 2>&1 | sed 's/^/    /' >> benchmark_summary.txt
echo "" >> benchmark_summary.txt

# Extract timing information
echo "📈 TIMING RESULTS:" >> benchmark_summary.txt
echo "==================" >> benchmark_summary.txt
echo "" >> benchmark_summary.txt

echo "v0.14.0 (Single Thread):" >> benchmark_summary.txt
grep -E "(Elapsed.*wall clock|Maximum resident set size)" benchmark_014_time.log >> benchmark_summary.txt
echo "" >> benchmark_summary.txt

echo "v0.15.0 (Single Thread):" >> benchmark_summary.txt
grep -E "(Elapsed.*wall clock|Maximum resident set size)" benchmark_015_time.log >> benchmark_summary.txt
echo "" >> benchmark_summary.txt

# Multithreaded results
for threads in 2 4 8; do
    if [[ -f "benchmark_014_t${threads}.log" ]]; then
        echo "v0.14.0 ($threads threads):" >> benchmark_summary.txt
        grep -E "(Elapsed.*wall clock|Maximum resident set size)" benchmark_014_t${threads}.log >> benchmark_summary.txt
        echo "" >> benchmark_summary.txt
        
        echo "v0.15.0 ($threads threads):" >> benchmark_summary.txt
        grep -E "(Elapsed.*wall clock|Maximum resident set size)" benchmark_015_t${threads}.log >> benchmark_summary.txt
        echo "" >> benchmark_summary.txt
    fi
done

# Enhanced performance analysis
echo "🚀 PERFORMANCE ANALYSIS:" >> benchmark_summary.txt
echo "=========================" >> benchmark_summary.txt
echo "" >> benchmark_summary.txt

# Dataset characteristics
echo "Dataset Characteristics:" >> benchmark_summary.txt
echo "Total chromosomes: $(cut -f1 benchmark_repeats.bed | sort | uniq | wc -l)" >> benchmark_summary.txt
echo "Total STR loci: $(wc -l < benchmark_repeats.bed)" >> benchmark_summary.txt
echo "STRs per chromosome:" >> benchmark_summary.txt
cut -f1 benchmark_repeats.bed | sort | uniq -c | sort -nr | head -10 | sed 's/^/  /' >> benchmark_summary.txt
echo "" >> benchmark_summary.txt

# Extract timing data for analysis
v014_1t=$(grep -E "Elapsed.*wall clock" benchmark_014_time.log | sed 's/.*time.*: //' || echo "N/A")
v015_1t=$(grep -E "Elapsed.*wall clock" benchmark_015_time.log | sed 's/.*time.*: //' || echo "N/A")

if [[ "$v014_1t" != "N/A" && "$v015_1t" != "N/A" ]]; then
    speedup=$(calculate_speedup "$v014_1t" "$v015_1t")
    echo "Single Thread Speedup: ${speedup}x faster (${v014_1t} → ${v015_1t})" >> benchmark_summary.txt
fi

echo "" >> benchmark_summary.txt
echo "Thread Scaling Analysis:" >> benchmark_summary.txt

# Calculate v0.15.0 internal thread scaling
if [[ "$BC_AVAILABLE" == true ]]; then
    v015_1t_sec=$(time_to_seconds "$v015_1t" 2>/dev/null || echo "0")
    echo "v0.15.0 Thread Scaling vs Single Thread:" >> benchmark_summary.txt
    
    for threads in 2 4 8; do
        if [[ -f "benchmark_015_t${threads}.log" ]]; then
            v015_time=$(grep -E "Elapsed.*wall clock" benchmark_015_t${threads}.log | sed 's/.*time.*: //')
            if [[ -n "$v015_time" && "$v015_1t_sec" != "0" ]]; then
                v015_time_sec=$(time_to_seconds "$v015_time")
                internal_speedup=$(echo "scale=2; $v015_1t_sec / $v015_time_sec" | bc -l)
                efficiency=$(echo "scale=1; ($internal_speedup / $threads) * 100" | bc -l)
                echo "  $threads threads: ${internal_speedup}x speedup (${efficiency}% efficiency)" >> benchmark_summary.txt
            fi
        fi
    done
    echo "" >> benchmark_summary.txt
fi

echo "Cross-Version Thread Comparison:" >> benchmark_summary.txt
for threads in 2 4 8; do
    if [[ -f "benchmark_014_t${threads}.log" && -f "benchmark_015_t${threads}.log" ]]; then
        v014_time=$(grep -E "Elapsed.*wall clock" benchmark_014_t${threads}.log | sed 's/.*time.*: //')
        v015_time=$(grep -E "Elapsed.*wall clock" benchmark_015_t${threads}.log | sed 's/.*time.*: //')
        
        if [[ -n "$v014_time" && -n "$v015_time" ]]; then
            ratio=$(calculate_speedup "$v014_time" "$v015_time")
            if (( $(echo "$ratio > 1" | bc -l) )); then
                echo "$threads threads: v0.15.0 is ${ratio}x faster" >> benchmark_summary.txt
            else
                inverse_ratio=$(calculate_speedup "$v015_time" "$v014_time")
                echo "$threads threads: v0.14.0 is ${inverse_ratio}x faster" >> benchmark_summary.txt
            fi
        fi
    fi
done

echo "" >> benchmark_summary.txt

# Memory usage comparison
echo "Memory Usage Comparison:" >> benchmark_summary.txt
v014_mem=$(grep -E "Maximum resident set size" benchmark_014_time.log | sed 's/.*: //' | sed 's/ .*//')
v015_mem=$(grep -E "Maximum resident set size" benchmark_015_time.log | sed 's/.*: //' | sed 's/ .*//')

if [[ -n "$v014_mem" && -n "$v015_mem" ]]; then
    mem_diff=$(echo "scale=1; ($v015_mem - $v014_mem) / 1024" | bc -l)
    if (( $(echo "$mem_diff > 0" | bc -l) )); then
        echo "v0.15.0 uses ${mem_diff}MB more memory (single thread)" >> benchmark_summary.txt
    else
        mem_diff_abs=$(echo "$mem_diff * -1" | bc -l)
        echo "v0.15.0 uses ${mem_diff_abs}MB less memory (single thread)" >> benchmark_summary.txt
    fi
fi

echo "" >> benchmark_summary.txt

# Memory profiling results
if [[ "$VALGRIND_AVAILABLE" == true ]]; then
    echo "📊 MEMORY PROFILING RESULTS:" >> benchmark_summary.txt
    echo "=============================" >> benchmark_summary.txt
    echo "" >> benchmark_summary.txt
    
    # Extract peak memory usage from Massif reports
    echo "Peak Memory Usage:" >> benchmark_summary.txt
    
    # Debug: Check if memory report files exist
    echo "   [Debug] Memory report files:" >> benchmark_summary.txt
    echo "   memory_report_014.txt: $([ -f "memory_report_014.txt" ] && echo "exists ($(wc -l < memory_report_014.txt) lines)" || echo "missing")" >> benchmark_summary.txt
    echo "   memory_report_015.txt: $([ -f "memory_report_015.txt" ] && echo "exists ($(wc -l < memory_report_015.txt) lines)" || echo "missing")" >> benchmark_summary.txt
    
    # Show sample lines for debugging
    if [[ -f "memory_report_014.txt" ]]; then
        echo "   [Debug] Sample from memory_report_014.txt:" >> benchmark_summary.txt
        echo "   $(head -5 memory_report_014.txt | tail -1 | sed 's/^/   /')" >> benchmark_summary.txt
    fi
    echo "" >> benchmark_summary.txt
    
    # Try multiple patterns to extract peak memory from ms_print output
    v014_peak="N/A"
    v015_peak="N/A"
    
    # Function to extract peak memory from ms_print output
    extract_peak_memory() {
        local file="$1"
        local peak="N/A"
        
        if [[ -f "$file" ]]; then
            # Try different parsing strategies
            
            # Strategy 1: Find the maximum heap size in typical ms_print format
            peak=$(awk '
                /^ *[0-9]/ && NF >= 2 {
                    gsub(/,/, "", $2)
                    if ($2 ~ /^[0-9]+$/ && $2 > max) max = $2
                }
                END { print max > 0 ? max : "N/A" }
            ' "$file" 2>/dev/null)
            
            # Strategy 2: If that failed, try looking for explicit mem_heap_B lines
            if [[ "$peak" == "N/A" ]]; then
                peak=$(grep -o "mem_heap_B=[0-9,]*" "$file" 2>/dev/null | \
                       sed 's/mem_heap_B=//' | sed 's/,//g' | \
                       sort -n | tail -1 || echo "N/A")
            fi
            
            # Strategy 3: Try parsing snapshot lines with MB/KB units
            if [[ "$peak" == "N/A" ]]; then
                peak=$(grep -E "[0-9,]+\s*(MB|KB|bytes)" "$file" 2>/dev/null | \
                       sed 's/[^0-9,]//g' | sed 's/,//g' | \
                       sort -n | tail -1 || echo "N/A")
            fi
        fi
        
        echo "$peak"
    }
    
    # Extract peak memory for both versions
    v014_peak=$(extract_peak_memory "memory_report_014.txt")
    v015_peak=$(extract_peak_memory "memory_report_015.txt")
    
    # Fallback: Extract directly from massif files if ms_print parsing failed
    if [[ "$v014_peak" == "N/A" && -f "massif_014.out" ]]; then
        v014_peak=$(grep "^mem_heap_B=" massif_014.out 2>/dev/null | sed 's/mem_heap_B=//' | sort -n | tail -1 || echo "N/A")
    fi
    
    if [[ "$v015_peak" == "N/A" && -f "massif_015.out" ]]; then
        v015_peak=$(grep "^mem_heap_B=" massif_015.out 2>/dev/null | sed 's/mem_heap_B=//' | sort -n | tail -1 || echo "N/A")
    fi
    
    echo "v0.14.0: ${v014_peak} bytes" >> benchmark_summary.txt
    echo "v0.15.0: ${v015_peak} bytes" >> benchmark_summary.txt
    echo "" >> benchmark_summary.txt
    
    # Check for memory leaks
    echo "Memory Leak Analysis:" >> benchmark_summary.txt
    echo "v0.14.0 leaks: $(grep -c "definitely lost\|possibly lost" memcheck_014.log || echo "0") issues" >> benchmark_summary.txt
    echo "v0.15.0 leaks: $(grep -c "definitely lost\|possibly lost" memcheck_015.log || echo "0") issues" >> benchmark_summary.txt
    echo "" >> benchmark_summary.txt
fi

# Output verification
echo "🔍 OUTPUT VERIFICATION:" >> benchmark_summary.txt
echo "=======================" >> benchmark_summary.txt
echo "" >> benchmark_summary.txt
echo "v0.14.0 output lines: $(wc -l < output_014_time.tsv)" >> benchmark_summary.txt
echo "v0.15.0 output lines: $(wc -l < output_015_time.tsv)" >> benchmark_summary.txt

# Compare outputs for correctness
sort output_014_time.tsv > output_014_sorted.tsv
sort output_015_time.tsv > output_015_sorted.tsv
diff --side-by-side --suppress-common-lines output_014_sorted.tsv output_015_sorted.tsv > output_diff.txt || true
echo "Output differences: $(wc -l < output_diff.txt) lines" >> benchmark_summary.txt

if [[ $(wc -l < output_diff.txt) -eq 0 ]]; then
    echo "✅ Outputs are identical (after sorting)" >> benchmark_summary.txt
else
    echo "⚠️  Outputs differ - check output_diff.txt for details" >> benchmark_summary.txt
    echo "   Use: diff --side-by-side --suppress-common-lines output_014_sorted.tsv output_015_sorted.tsv" >> benchmark_summary.txt
fi

echo "" >> benchmark_summary.txt

# =============================================================================
# DISPLAY FINAL RESULTS
# =============================================================================
echo ""
echo "🎉 BENCHMARKING COMPLETE!"
echo "=========================="
cat benchmark_summary.txt

echo ""
echo "📁 Generated Files:"
echo "   • benchmark_summary.txt    - Complete results summary"
echo "   • output_014_time.tsv      - v0.14.0 output (timing run)"
echo "   • output_015_time.tsv      - v0.15.0 output (timing run)"
echo "   • benchmark_014_time.log   - v0.14.0 detailed timing"
echo "   • benchmark_015_time.log   - v0.15.0 detailed timing"

if [[ "$VALGRIND_AVAILABLE" == true ]]; then
    echo "   • memory_report_014.txt    - v0.14.0 memory profile"
    echo "   • memory_report_015.txt    - v0.15.0 memory profile"
    echo "   • massif_014.out          - v0.14.0 Massif data"
    echo "   • massif_015.out          - v0.15.0 Massif data"
    echo "   • memcheck_014.log        - v0.14.0 memory leak check"
    echo "   • memcheck_015.log        - v0.15.0 memory leak check"
    echo "   • valgrind_014.log        - v0.14.0 Valgrind output"
    echo "   • valgrind_015.log        - v0.15.0 Valgrind output"
fi

echo "   • output_diff.txt          - Output comparison"
echo "   • benchmark_*_t*.log       - Multithreaded results"
echo ""

echo "Completed at: $(date)"
echo ""

echo "💡 TIPS:"
if [[ "$VALGRIND_AVAILABLE" == true ]]; then
    echo "   📊 Visualize memory usage: massif-visualizer massif_014.out"
    echo "   🚀 For faster runs on large datasets: $0 --skip-valgrind"
else
    echo "   🧠 For memory profiling: Install valgrind and re-run without --skip-valgrind"
    echo "   ⚡ Current mode optimized for speed on large datasets"
fi