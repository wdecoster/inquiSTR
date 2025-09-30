#!/bin/bash
# Complete inquiSTR benchmarking script with optional Valgrind memory profiling
# Usage: ./benchmark.sh --current <binary> [--prev <binary>] [--skip-valgrind] [--help]

set -e  # Exit on any error

# Initialize variables
SKIP_VALGRIND_ARG=false
CURRENT_BINARY=""
PREV_BINARY=""
COMPARISON_MODE=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --current)
            CURRENT_BINARY="$2"
            shift 2
            ;;
        --prev)
            PREV_BINARY="$2"
            shift 2
            ;;
        --skip-valgrind)
            SKIP_VALGRIND_ARG=true
            shift
            ;;
        --help|-h)
            echo "Usage: $0 --current <binary> [--prev <binary>] [OPTIONS]"
            echo ""
            echo "REQUIRED:"
            echo "  --current <path>   Path to current binary to benchmark"
            echo ""
            echo "OPTIONAL:"
            echo "  --prev <path>      Path to previous binary for comparison"
            echo "                     If provided, enables comparison mode"
            echo ""
            echo "OPTIONS:"
            echo "  --skip-valgrind    Disable Valgrind memory profiling (faster execution)"
            echo "  --help, -h         Show this help message"
            echo ""
            echo "ENVIRONMENT VARIABLES:"
            echo "  SKIP_VALGRIND=1    Alternative way to disable Valgrind profiling"
            echo ""
            echo "Examples:"
            echo "  $0 --current ./target/release/inquiSTR"
            echo "  $0 --current ./target/release/inquiSTR --prev ./inquiSTR_old"
            echo "  $0 --current ./target/release/inquiSTR --skip-valgrind"
            echo "  $0 --current ./target/release/inquiSTR --prev ./inquiSTR_old --skip-valgrind"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$CURRENT_BINARY" ]]; then
    echo "ERROR: --current <binary> is required"
    echo "Use --help for usage information"
    exit 1
fi

# Check if comparison mode is enabled
if [[ -n "$PREV_BINARY" ]]; then
    COMPARISON_MODE=true
fi

# Determine if Valgrind should be skipped (command line arg takes precedence)
if [[ "$SKIP_VALGRIND_ARG" == true || "$SKIP_VALGRIND" == "1" ]]; then
    SKIP_VALGRIND_ENABLED=true
    if [[ "$COMPARISON_MODE" == true ]]; then
        echo "=== inquiSTR Fast Benchmarking Script (Comparison Mode) ==="
        echo "Comparing $(basename "$PREV_BINARY") vs $(basename "$CURRENT_BINARY") (Valgrind disabled for speed)"
    else
        echo "=== inquiSTR Fast Benchmarking Script (Single Binary) ==="
        echo "Testing $(basename "$CURRENT_BINARY") (Valgrind disabled for speed)"
    fi
else
    SKIP_VALGRIND_ENABLED=false
    if [[ "$COMPARISON_MODE" == true ]]; then
        echo "=== inquiSTR Performance Benchmarking Script (Comparison Mode) ==="
        echo "Comparing $(basename "$PREV_BINARY") vs $(basename "$CURRENT_BINARY") with Valgrind memory profiling"
    else
        echo "=== inquiSTR Performance Benchmarking Script (Single Binary) ==="
        echo "Testing $(basename "$CURRENT_BINARY") with Valgrind memory profiling"
    fi
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

# Check if required binaries exist
if [[ ! -f "$CURRENT_BINARY" ]]; then
    echo "ERROR: Current binary not found: $CURRENT_BINARY"
    exit 1
fi

if [[ "$COMPARISON_MODE" == true && ! -f "$PREV_BINARY" ]]; then
    echo "ERROR: Previous binary not found: $PREV_BINARY"
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

if [[ "$COMPARISON_MODE" == true ]]; then
    echo "Previous binary ($(basename "$PREV_BINARY")):"
    "$PREV_BINARY" --version 2>/dev/null || echo "  Version information not available"
    echo ""
fi

echo "Current binary ($(basename "$CURRENT_BINARY")):"
"$CURRENT_BINARY" --version 2>/dev/null || echo "  Version information not available"

echo ""
echo "Build information:"
if [[ "$COMPARISON_MODE" == true ]]; then
    echo "  Previous binary size: $(ls -lhL "$PREV_BINARY" | awk '{print $5}')"
    echo "  Previous binary modified: $(stat -L -c %y "$PREV_BINARY" 2>/dev/null | cut -d. -f1 || echo "unknown")"
fi
echo "  Current binary size: $(ls -lhL "$CURRENT_BINARY" | awk '{print $5}')"
echo "  Current binary modified: $(stat -L -c %y "$CURRENT_BINARY" 2>/dev/null | cut -d. -f1 || echo "unknown")"
echo ""

# Clean previous results
echo "Cleaning previous benchmark results..."
rm -f output_prev*.tsv output_current*.tsv benchmark_*.log memory_*.txt massif_*.out valgrind_*.log
echo ""

# =============================================================================
# TIMING BENCHMARKS
# =============================================================================
echo "=== TIMING BENCHMARKS ==="

if [[ "$COMPARISON_MODE" == true ]]; then
    echo "📊 Testing previous binary ($(basename "$PREV_BINARY")) (timing)..."
    /usr/bin/time -v "$PREV_BINARY" call alignment.cram --region-file benchmark_repeats.bed --threads 1 > output_prev_time.tsv 2> benchmark_prev_time.log
    echo "✓ Previous binary timing complete"
fi

echo "📊 Testing current binary ($(basename "$CURRENT_BINARY")) (timing)..."
/usr/bin/time -v "$CURRENT_BINARY" call alignment.cram --region-file benchmark_repeats.bed --threads 1 > output_current_time.tsv 2> benchmark_current_time.log
echo "✓ Current binary timing complete"

echo ""

# =============================================================================
# VALGRIND MEMORY PROFILING (if available and enabled)
# =============================================================================
if [[ "$VALGRIND_AVAILABLE" == true ]]; then
    echo "=== VALGRIND MEMORY PROFILING ==="
    
    if [[ "$COMPARISON_MODE" == true ]]; then
        echo "🧠 Running previous binary ($(basename "$PREV_BINARY")) with Valgrind Massif (memory profiling)..."
        echo "   This may take significantly longer than normal execution..."
        valgrind --tool=massif \
                 --massif-out-file=massif_prev.out \
                 --time-unit=ms \
                 --detailed-freq=1 \
                 --max-snapshots=200 \
                 "$PREV_BINARY" call alignment.cram --region-file benchmark_repeats.bed --threads 1 \
                 > output_prev_valgrind.tsv 2> valgrind_prev.log
        
        # Generate memory report for previous binary
        ms_print massif_prev.out > memory_report_prev.txt
        echo "✓ Previous binary Valgrind profiling complete"
    fi
    
    echo "🧠 Running current binary ($(basename "$CURRENT_BINARY")) with Valgrind Massif (memory profiling)..."
    echo "   This may take significantly longer than normal execution..."
    valgrind --tool=massif \
             --massif-out-file=massif_current.out \
             --time-unit=ms \
             --detailed-freq=1 \
             --max-snapshots=200 \
             "$CURRENT_BINARY" call alignment.cram --region-file benchmark_repeats.bed --threads 1 \
             > output_current_valgrind.tsv 2> valgrind_current.log
    
    # Generate memory report for current binary
    ms_print massif_current.out > memory_report_current.txt
    echo "✓ Current binary Valgrind profiling complete"
    
    echo ""
    
    # Additional Valgrind checks for memory leaks
    echo "🔍 Running Valgrind Memcheck for memory leaks..."
    
    if [[ "$COMPARISON_MODE" == true ]]; then
        echo "   Checking previous binary ($(basename "$PREV_BINARY")) for memory leaks..."
        valgrind --tool=memcheck \
                 --leak-check=full \
                 --show-leak-kinds=all \
                 --track-origins=yes \
                 --log-file=memcheck_prev.log \
                 "$PREV_BINARY" call alignment.cram --region-file benchmark_repeats.bed --threads 1 \
                 > /dev/null 2>&1
    fi
    
    echo "   Checking current binary ($(basename "$CURRENT_BINARY")) for memory leaks..."
    valgrind --tool=memcheck \
             --leak-check=full \
             --show-leak-kinds=all \
             --track-origins=yes \
             --log-file=memcheck_current.log \
             "$CURRENT_BINARY" call alignment.cram --region-file benchmark_repeats.bed --threads 1 \
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
        
        if [[ "$COMPARISON_MODE" == true ]]; then
            echo "   Previous binary ($threads threads)..."
            /usr/bin/time -v "$PREV_BINARY" call alignment.cram --region-file benchmark_repeats.bed --threads $threads > output_prev_t${threads}.tsv 2> benchmark_prev_t${threads}.log
        fi
        
        echo "   Current binary ($threads threads)..."
        /usr/bin/time -v "$CURRENT_BINARY" call alignment.cram --region-file benchmark_repeats.bed --threads $threads > output_current_t${threads}.tsv 2> benchmark_current_t${threads}.log
        
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

# Binary version information
echo "📦 BINARY VERSION INFORMATION:" >> benchmark_summary.txt
echo "==============================" >> benchmark_summary.txt
echo "" >> benchmark_summary.txt

if [[ "$COMPARISON_MODE" == true ]]; then
    echo "Previous Binary ($(basename "$PREV_BINARY")):" >> benchmark_summary.txt
    echo "  Location: $PREV_BINARY" >> benchmark_summary.txt
    echo "  File size: $(stat -L -f%z "$PREV_BINARY" 2>/dev/null || stat -L -c%s "$PREV_BINARY")B" >> benchmark_summary.txt
    echo "  Modified: $(stat -L -f%Sm -t%Y-%m-%d\ %H:%M:%S "$PREV_BINARY" 2>/dev/null || stat -L -c%y "$PREV_BINARY")" >> benchmark_summary.txt
    echo "  Version output:" >> benchmark_summary.txt
    "$PREV_BINARY" --version 2>&1 | sed 's/^/    /' >> benchmark_summary.txt
    echo "" >> benchmark_summary.txt
fi

echo "Current Binary ($(basename "$CURRENT_BINARY")):" >> benchmark_summary.txt
echo "  Location: $CURRENT_BINARY" >> benchmark_summary.txt
echo "  File size: $(stat -L -f%z "$CURRENT_BINARY" 2>/dev/null || stat -L -c%s "$CURRENT_BINARY")B" >> benchmark_summary.txt
echo "  Modified: $(stat -L -f%Sm -t%Y-%m-%d\ %H:%M:%S "$CURRENT_BINARY" 2>/dev/null || stat -L -c%y "$CURRENT_BINARY")" >> benchmark_summary.txt
echo "  Version output:" >> benchmark_summary.txt
"$CURRENT_BINARY" --version 2>&1 | sed 's/^/    /' >> benchmark_summary.txt
echo "" >> benchmark_summary.txt

# Extract timing information
echo "📈 TIMING RESULTS:" >> benchmark_summary.txt
echo "==================" >> benchmark_summary.txt
echo "" >> benchmark_summary.txt

if [[ "$COMPARISON_MODE" == true ]]; then
    echo "Previous Binary (Single Thread):" >> benchmark_summary.txt
    grep -E "(Elapsed.*wall clock|Maximum resident set size)" benchmark_prev_time.log >> benchmark_summary.txt
    echo "" >> benchmark_summary.txt
fi

echo "Current Binary (Single Thread):" >> benchmark_summary.txt
grep -E "(Elapsed.*wall clock|Maximum resident set size)" benchmark_current_time.log >> benchmark_summary.txt
echo "" >> benchmark_summary.txt

# Multithreaded results
for threads in 2 4 8; do
    if [[ "$COMPARISON_MODE" == true && -f "benchmark_prev_t${threads}.log" ]]; then
        echo "Previous Binary ($threads threads):" >> benchmark_summary.txt
        grep -E "(Elapsed.*wall clock|Maximum resident set size)" benchmark_prev_t${threads}.log >> benchmark_summary.txt
        echo "" >> benchmark_summary.txt
    fi
    
    if [[ -f "benchmark_current_t${threads}.log" ]]; then
        echo "Current Binary ($threads threads):" >> benchmark_summary.txt
        grep -E "(Elapsed.*wall clock|Maximum resident set size)" benchmark_current_t${threads}.log >> benchmark_summary.txt
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
if [[ "$COMPARISON_MODE" == true ]]; then
    prev_1t=$(grep -E "Elapsed.*wall clock" benchmark_prev_time.log | sed 's/.*time.*: //' || echo "N/A")
fi
current_1t=$(grep -E "Elapsed.*wall clock" benchmark_current_time.log | sed 's/.*time.*: //' || echo "N/A")

if [[ "$COMPARISON_MODE" == true && "$prev_1t" != "N/A" && "$current_1t" != "N/A" ]]; then
    speedup=$(calculate_speedup "$prev_1t" "$current_1t")
    echo "Single Thread Speedup: ${speedup}x faster (${prev_1t} → ${current_1t})" >> benchmark_summary.txt
fi

echo "" >> benchmark_summary.txt
echo "Thread Scaling Analysis:" >> benchmark_summary.txt

# Calculate current binary internal thread scaling
if [[ "$BC_AVAILABLE" == true && "$current_1t" != "N/A" ]]; then
    current_1t_sec=$(time_to_seconds "$current_1t" 2>/dev/null || echo "0")
    echo "Current Binary Thread Scaling vs Single Thread:" >> benchmark_summary.txt
    
    for threads in 2 4 8; do
        if [[ -f "benchmark_current_t${threads}.log" ]]; then
            current_time=$(grep -E "Elapsed.*wall clock" benchmark_current_t${threads}.log | sed 's/.*time.*: //')
            if [[ -n "$current_time" && "$current_1t_sec" != "0" ]]; then
                current_time_sec=$(time_to_seconds "$current_time")
                internal_speedup=$(echo "scale=2; $current_1t_sec / $current_time_sec" | bc -l)
                efficiency=$(echo "scale=1; ($internal_speedup / $threads) * 100" | bc -l)
                echo "  $threads threads: ${internal_speedup}x speedup (${efficiency}% efficiency)" >> benchmark_summary.txt
            fi
        fi
    done
    echo "" >> benchmark_summary.txt
fi

if [[ "$COMPARISON_MODE" == true ]]; then
    echo "Cross-Binary Thread Comparison:" >> benchmark_summary.txt
    for threads in 2 4 8; do
        if [[ -f "benchmark_prev_t${threads}.log" && -f "benchmark_current_t${threads}.log" ]]; then
            prev_time=$(grep -E "Elapsed.*wall clock" benchmark_prev_t${threads}.log | sed 's/.*time.*: //')
            current_time=$(grep -E "Elapsed.*wall clock" benchmark_current_t${threads}.log | sed 's/.*time.*: //')
            
            if [[ -n "$prev_time" && -n "$current_time" ]]; then
                ratio=$(calculate_speedup "$prev_time" "$current_time")
                if (( $(echo "$ratio > 1" | bc -l) )); then
                    echo "$threads threads: Current binary is ${ratio}x faster" >> benchmark_summary.txt
                else
                    inverse_ratio=$(calculate_speedup "$current_time" "$prev_time")
                    echo "$threads threads: Previous binary is ${inverse_ratio}x faster" >> benchmark_summary.txt
                fi
            fi
        fi
    done
fi

echo "" >> benchmark_summary.txt

# Memory usage comparison
echo "Memory Usage Comparison:" >> benchmark_summary.txt
if [[ "$COMPARISON_MODE" == true ]]; then
    prev_mem=$(grep -E "Maximum resident set size" benchmark_prev_time.log | sed 's/.*: //' | sed 's/ .*//')
fi
current_mem=$(grep -E "Maximum resident set size" benchmark_current_time.log | sed 's/.*: //' | sed 's/ .*//')

if [[ "$COMPARISON_MODE" == true && -n "$prev_mem" && -n "$current_mem" ]]; then
    mem_diff=$(echo "scale=1; ($current_mem - $prev_mem) / 1024" | bc -l)
    if (( $(echo "$mem_diff > 0" | bc -l) )); then
        echo "Current binary uses ${mem_diff}MB more memory (single thread)" >> benchmark_summary.txt
    else
        mem_diff_abs=$(echo "$mem_diff * -1" | bc -l)
        echo "Current binary uses ${mem_diff_abs}MB less memory (single thread)" >> benchmark_summary.txt
    fi
elif [[ "$COMPARISON_MODE" == false && -n "$current_mem" ]]; then
    echo "Current binary memory usage: $(echo "scale=1; $current_mem / 1024" | bc -l)MB (single thread)" >> benchmark_summary.txt
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
    if [[ "$COMPARISON_MODE" == true ]]; then
        echo "   memory_report_prev.txt: $([ -f "memory_report_prev.txt" ] && echo "exists ($(wc -l < memory_report_prev.txt) lines)" || echo "missing")" >> benchmark_summary.txt
    fi
    echo "   memory_report_current.txt: $([ -f "memory_report_current.txt" ] && echo "exists ($(wc -l < memory_report_current.txt) lines)" || echo "missing")" >> benchmark_summary.txt
    
    # Show sample lines for debugging
    if [[ "$COMPARISON_MODE" == true && -f "memory_report_prev.txt" ]]; then
        echo "   [Debug] Sample from memory_report_prev.txt:" >> benchmark_summary.txt
        echo "   $(head -5 memory_report_prev.txt | tail -1 | sed 's/^/   /')" >> benchmark_summary.txt
    fi
    echo "" >> benchmark_summary.txt
    
    # Try multiple patterns to extract peak memory from ms_print output
    
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
    
    # Extract peak memory for binaries
    if [[ "$COMPARISON_MODE" == true ]]; then
        prev_peak=$(extract_peak_memory "memory_report_prev.txt")
        # Fallback: Extract directly from massif files if ms_print parsing failed
        if [[ "$prev_peak" == "N/A" && -f "massif_prev.out" ]]; then
            prev_peak=$(grep "^mem_heap_B=" massif_prev.out 2>/dev/null | sed 's/mem_heap_B=//' | sort -n | tail -1 || echo "N/A")
        fi
        echo "Previous binary: ${prev_peak} bytes" >> benchmark_summary.txt
    fi
    
    current_peak=$(extract_peak_memory "memory_report_current.txt")
    # Fallback: Extract directly from massif files if ms_print parsing failed
    if [[ "$current_peak" == "N/A" && -f "massif_current.out" ]]; then
        current_peak=$(grep "^mem_heap_B=" massif_current.out 2>/dev/null | sed 's/mem_heap_B=//' | sort -n | tail -1 || echo "N/A")
    fi
    echo "Current binary: ${current_peak} bytes" >> benchmark_summary.txt
    echo "" >> benchmark_summary.txt
    
    # Check for memory leaks
    echo "Memory Leak Analysis:" >> benchmark_summary.txt
    if [[ "$COMPARISON_MODE" == true ]]; then
        echo "Previous binary leaks: $(grep -c "definitely lost\|possibly lost" memcheck_prev.log || echo "0") issues" >> benchmark_summary.txt
    fi
    echo "Current binary leaks: $(grep -c "definitely lost\|possibly lost" memcheck_current.log || echo "0") issues" >> benchmark_summary.txt
    echo "" >> benchmark_summary.txt
fi

# Output verification
echo "🔍 OUTPUT VERIFICATION:" >> benchmark_summary.txt
echo "=======================" >> benchmark_summary.txt
echo "" >> benchmark_summary.txt

if [[ "$COMPARISON_MODE" == true ]]; then
    echo "Previous binary output lines: $(wc -l < output_prev_time.tsv)" >> benchmark_summary.txt
fi
echo "Current binary output lines: $(wc -l < output_current_time.tsv)" >> benchmark_summary.txt

# Compare outputs for correctness if in comparison mode
if [[ "$COMPARISON_MODE" == true ]]; then
    sort output_prev_time.tsv > output_prev_sorted.tsv
    sort output_current_time.tsv > output_current_sorted.tsv
    diff --side-by-side --suppress-common-lines output_prev_sorted.tsv output_current_sorted.tsv > output_diff.txt || true
    echo "Output differences: $(wc -l < output_diff.txt) lines" >> benchmark_summary.txt

    if [[ $(wc -l < output_diff.txt) -eq 0 ]]; then
        echo "✅ Outputs are identical (after sorting)" >> benchmark_summary.txt
    else
        echo "⚠️  Outputs differ - check output_diff.txt for details" >> benchmark_summary.txt
        echo "   Use: diff --side-by-side --suppress-common-lines output_prev_sorted.tsv output_current_sorted.tsv" >> benchmark_summary.txt
    fi
else
    echo "✅ Single binary mode - no output comparison performed" >> benchmark_summary.txt
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
echo "   • benchmark_summary.txt       - Complete results summary"

if [[ "$COMPARISON_MODE" == true ]]; then
    echo "   • output_prev_time.tsv        - Previous binary output (timing run)"
    echo "   • benchmark_prev_time.log     - Previous binary detailed timing"
fi
echo "   • output_current_time.tsv     - Current binary output (timing run)"
echo "   • benchmark_current_time.log  - Current binary detailed timing"

if [[ "$VALGRIND_AVAILABLE" == true ]]; then
    if [[ "$COMPARISON_MODE" == true ]]; then
        echo "   • memory_report_prev.txt      - Previous binary memory profile"
        echo "   • massif_prev.out            - Previous binary Massif data"
        echo "   • memcheck_prev.log          - Previous binary memory leak check"
        echo "   • valgrind_prev.log          - Previous binary Valgrind output"
    fi
    echo "   • memory_report_current.txt   - Current binary memory profile"
    echo "   • massif_current.out         - Current binary Massif data"
    echo "   • memcheck_current.log       - Current binary memory leak check"
    echo "   • valgrind_current.log       - Current binary Valgrind output"
fi

if [[ "$COMPARISON_MODE" == true ]]; then
    echo "   • output_diff.txt             - Output comparison"
fi
echo "   • benchmark_*_t*.log          - Multithreaded results"
echo ""

echo "Completed at: $(date)"
echo ""

echo "💡 TIPS:"
if [[ "$VALGRIND_AVAILABLE" == true ]]; then
    if [[ "$COMPARISON_MODE" == true ]]; then
        echo "   📊 Visualize memory usage: massif-visualizer massif_prev.out massif_current.out"
    else
        echo "   📊 Visualize memory usage: massif-visualizer massif_current.out"
    fi
    echo "   🚀 For faster runs on large datasets: add --skip-valgrind"
else
    echo "   🧠 For memory profiling: Install valgrind and re-run without --skip-valgrind"
    echo "   ⚡ Current mode optimized for speed on large datasets"
fi
if [[ "$COMPARISON_MODE" == false ]]; then
    echo "   📊 For comparison benchmarks: use --prev <binary_path> --current <binary_path>"
fi