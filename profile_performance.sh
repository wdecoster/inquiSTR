#!/bin/bash
# Advanced performance profiling for inquiSTR with flamegraphs
# Usage: ./profile_performance.sh [--threads N] [--help]

set -e

# Default configuration
THREADS=1
DURATION=30
SAMPLE_RATE=997  # Prime number for better sampling distribution

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --duration)
            DURATION="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "OPTIONS:"
            echo "  --threads N        Number of threads to use (default: 1)"
            echo "  --duration N       Profiling duration in seconds (default: 30)"
            echo "  --help, -h         Show this help message"
            echo ""
            echo "REQUIREMENTS:"
            echo "  - perf (Linux performance tools)"
            echo "  - FlameGraph (https://github.com/brendangregg/FlameGraph)"
            echo ""
            echo "INSTALLATION:"
            echo "  Ubuntu/Debian: sudo apt install linux-tools-common linux-tools-generic"
            echo "  FlameGraph: git clone https://github.com/brendangregg/FlameGraph.git"
            echo ""
            echo "Examples:"
            echo "  $0                           # Profile single-thread"
            echo "  $0 --threads 4              # Profile 4-thread execution"
            echo "  $0 --threads 1 --duration 60  # Profile for 60 seconds"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo "=== inquiSTR Advanced Performance Profiling ==="
echo "Threads: $THREADS"
echo "Duration: ${DURATION}s"
echo "Started at: $(date)"
echo ""

# Check dependencies
MISSING_DEPS=()

if ! command -v perf &> /dev/null; then
    MISSING_DEPS+=("perf")
fi

if [[ ! -d "FlameGraph" ]]; then
    if ! command -v flamegraph.pl &> /dev/null; then
        MISSING_DEPS+=("FlameGraph")
    fi
fi

if [[ ! -f "./inquiSTR_015" ]]; then
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

# Install missing dependencies
if [[ ${#MISSING_DEPS[@]} -gt 0 ]]; then
    echo "⚠️  Missing dependencies: ${MISSING_DEPS[*]}"
    echo ""
    
    for dep in "${MISSING_DEPS[@]}"; do
        case $dep in
            perf)
                echo "Installing perf..."
                echo "Run: sudo apt install linux-tools-common linux-tools-generic linux-tools-\$(uname -r)"
                echo "Or: sudo yum install perf  # RHEL/CentOS"
                ;;
            FlameGraph)
                echo "Installing FlameGraph..."
                if command -v git &> /dev/null; then
                    git clone https://github.com/brendangregg/FlameGraph.git || echo "FlameGraph clone failed"
                else
                    echo "Git not found. Please install FlameGraph manually:"
                    echo "git clone https://github.com/brendangregg/FlameGraph.git"
                fi
                ;;
        esac
    done
    
    echo ""
    echo "📋 REQUIRED DEPENDENCIES:"
    echo "========================"
    echo ""
    echo "Missing: ${MISSING_DEPS[*]}"
    echo ""
    echo "To install all dependencies, run:"
    echo ""
    echo "1. Install perf (Linux performance tools):"
    echo "   sudo apt install linux-tools-common linux-tools-generic linux-tools-\$(uname -r)"
    echo "   # Or on RHEL/CentOS: sudo yum install perf"
    echo ""
    echo "2. Install FlameGraph:"
    echo "   git clone https://github.com/brendangregg/FlameGraph.git"
    echo ""
    echo "3. Fix perf permissions (if needed):"
    echo "   echo -1 | sudo tee /proc/sys/kernel/perf_event_paranoid"
    echo ""
    echo "Then re-run this script: ./profile_performance.sh"
    echo ""
    exit 1
fi

# Set up FlameGraph path
if [[ -d "FlameGraph" ]]; then
    FLAMEGRAPH_PATH="./FlameGraph"
else
    FLAMEGRAPH_PATH=""  # Assume it's in PATH
fi

echo "✓ All dependencies found"
echo ""

# Create profiling output directory
PROFILE_DIR="profile_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$PROFILE_DIR"
echo "📁 Results will be saved to: $PROFILE_DIR"
echo ""

# Function to run perf with better error handling
run_perf_profile() {
    local output_name="$1"
    local extra_events="$2"
    
    echo "🔍 Running CPU profiling: $output_name"
    
    # Check if we can run perf record
    if ! perf list | grep -q "cpu-cycles"; then
        echo "⚠️  CPU profiling not available (try: echo -1 | sudo tee /proc/sys/kernel/perf_event_paranoid)"
        return 1
    fi
    
    local perf_cmd="perf record"
    perf_cmd="$perf_cmd -F $SAMPLE_RATE"  # Sample rate
    perf_cmd="$perf_cmd -g"               # Call graph
    perf_cmd="$perf_cmd --call-graph dwarf"  # DWARF call graph for better Rust support
    
    if [[ -n "$extra_events" ]]; then
        perf_cmd="$perf_cmd -e $extra_events"
    else
        perf_cmd="$perf_cmd -e cpu-cycles"
    fi
    
    perf_cmd="$perf_cmd -o $PROFILE_DIR/perf_${output_name}.data"
    perf_cmd="$perf_cmd -- ./inquiSTR_015 call alignment.cram --region-file benchmark_repeats.bed --threads $THREADS"
    
    echo "Command: $perf_cmd"
    eval "$perf_cmd" 2>&1 | tee "$PROFILE_DIR/perf_${output_name}.log"
    
    if [[ -f "$PROFILE_DIR/perf_${output_name}.data" ]]; then
        echo "✓ Profiling data collected: perf_${output_name}.data"
        return 0
    else
        echo "❌ Profiling failed for $output_name"
        return 1
    fi
}

# Function to generate flamegraph
generate_flamegraph() {
    local perf_file="$1"
    local output_name="$2"
    
    if [[ ! -f "$perf_file" ]]; then
        echo "❌ Perf data file not found: $perf_file"
        return 1
    fi
    
    echo "🔥 Generating flamegraph: $output_name"
    
    # Generate perf script output
    perf script -i "$perf_file" > "$PROFILE_DIR/${output_name}_perf.out"
    
    # Generate flamegraph
    if [[ -n "$FLAMEGRAPH_PATH" ]]; then
        cat "$PROFILE_DIR/${output_name}_perf.out" | \
        "$FLAMEGRAPH_PATH/stackcollapse-perf.pl" | \
        "$FLAMEGRAPH_PATH/flamegraph.pl" --title "inquiSTR $output_name Profile" \
        > "$PROFILE_DIR/${output_name}_flamegraph.svg"
    else
        cat "$PROFILE_DIR/${output_name}_perf.out" | \
        stackcollapse-perf.pl | \
        flamegraph.pl --title "inquiSTR $output_name Profile" \
        > "$PROFILE_DIR/${output_name}_flamegraph.svg"
    fi
    
    if [[ -f "$PROFILE_DIR/${output_name}_flamegraph.svg" ]]; then
        echo "✓ Flamegraph generated: ${output_name}_flamegraph.svg"
        return 0
    else
        echo "❌ Flamegraph generation failed for $output_name"
        return 1
    fi
}

# Function to analyze perf data
analyze_perf_data() {
    local perf_file="$1"
    local output_name="$2"
    
    echo "📊 Analyzing performance data: $output_name"
    
    # Generate detailed reports
    perf report -i "$perf_file" --stdio > "$PROFILE_DIR/${output_name}_report.txt"
    perf annotate -i "$perf_file" --stdio > "$PROFILE_DIR/${output_name}_annotate.txt"
    
    # Top functions
    echo "Top 20 functions by CPU time:" > "$PROFILE_DIR/${output_name}_top_functions.txt"
    perf report -i "$perf_file" --stdio | grep -A 30 "# Overhead.*Command.*Shared Object.*Symbol" >> "$PROFILE_DIR/${output_name}_top_functions.txt"
    
    echo "✓ Performance analysis complete: ${output_name}_report.txt"
}

# =============================================================================
# MAIN PROFILING EXECUTION
# =============================================================================

echo "=== CPU PROFILING ==="

# Profile with different event types
declare -A PROFILE_CONFIGS=(
    ["cpu_cycles"]="cpu-cycles"
    ["instructions"]="instructions"
    ["cache_misses"]="cache-misses"
    ["branch_misses"]="branch-misses"
)

SUCCESSFUL_PROFILES=()

for profile_name in "${!PROFILE_CONFIGS[@]}"; do
    event="${PROFILE_CONFIGS[$profile_name]}"
    
    if run_perf_profile "$profile_name" "$event"; then
        SUCCESSFUL_PROFILES+=("$profile_name")
        analyze_perf_data "$PROFILE_DIR/perf_${profile_name}.data" "$profile_name"
        generate_flamegraph "$PROFILE_DIR/perf_${profile_name}.data" "$profile_name"
    fi
    echo ""
done

# =============================================================================
# MEMORY PROFILING (if available)
# =============================================================================

if command -v valgrind &> /dev/null; then
    echo "=== MEMORY PROFILING ==="
    
    echo "🧠 Running Valgrind Callgrind profiling..."
    valgrind --tool=callgrind \
             --callgrind-out-file="$PROFILE_DIR/callgrind.out" \
             --separate-threads=yes \
             --cache-sim=yes \
             --branch-sim=yes \
             ./inquiSTR_015 call alignment.cram --region-file benchmark_repeats.bed --threads "$THREADS" \
             > "$PROFILE_DIR/valgrind_output.txt" 2>&1
    
    if [[ -f "$PROFILE_DIR/callgrind.out" ]]; then
        echo "✓ Callgrind profiling complete"
        
        # Generate callgrind report
        if command -v callgrind_annotate &> /dev/null; then
            callgrind_annotate "$PROFILE_DIR/callgrind.out" > "$PROFILE_DIR/callgrind_report.txt"
            echo "✓ Callgrind report generated"
        fi
    fi
    echo ""
fi

# =============================================================================
# I/O PROFILING
# =============================================================================

echo "=== I/O PROFILING ==="

if run_perf_profile "io_profile" "syscalls:sys_enter_read,syscalls:sys_enter_write,syscalls:sys_enter_openat"; then
    SUCCESSFUL_PROFILES+=("io_profile")
    analyze_perf_data "$PROFILE_DIR/perf_io_profile.data" "io_profile"
    
    # Detailed I/O analysis
    echo "📊 Analyzing I/O patterns..."
    perf script -i "$PROFILE_DIR/perf_io_profile.data" > "$PROFILE_DIR/io_trace.txt"
    echo "✓ I/O trace generated: io_trace.txt"
fi

echo ""

# =============================================================================
# RESULTS SUMMARY
# =============================================================================

echo "=== PROFILING RESULTS SUMMARY ==="
echo ""

echo "📁 Generated Files in $PROFILE_DIR:"
echo "=================================="

if [[ ${#SUCCESSFUL_PROFILES[@]} -gt 0 ]]; then
    echo "🔥 Flamegraphs:"
    for profile in "${SUCCESSFUL_PROFILES[@]}"; do
        if [[ -f "$PROFILE_DIR/${profile}_flamegraph.svg" ]]; then
            echo "   • ${profile}_flamegraph.svg - Interactive CPU flamegraph"
        fi
    done
    echo ""
    
    echo "📊 Performance Reports:"
    for profile in "${SUCCESSFUL_PROFILES[@]}"; do
        echo "   • ${profile}_report.txt - Detailed function-level analysis"
        echo "   • ${profile}_top_functions.txt - Top CPU-consuming functions"
    done
    echo ""
    
    echo "📈 Raw Data:"
    for profile in "${SUCCESSFUL_PROFILES[@]}"; do
        echo "   • perf_${profile}.data - Raw perf data (use: perf report -i file)"
    done
fi

if [[ -f "$PROFILE_DIR/callgrind.out" ]]; then
    echo ""
    echo "🧠 Memory Profiling:"
    echo "   • callgrind.out - Valgrind callgrind data"
    echo "   • callgrind_report.txt - Function call analysis"
    echo "   • valgrind_output.txt - Valgrind execution log"
fi

echo ""
echo "💡 NEXT STEPS:"
echo "=============="
echo ""
echo "1. 🔍 View flamegraphs in browser:"
for profile in "${SUCCESSFUL_PROFILES[@]}"; do
    if [[ -f "$PROFILE_DIR/${profile}_flamegraph.svg" ]]; then
        echo "   firefox $PROFILE_DIR/${profile}_flamegraph.svg"
        break
    fi
done

echo ""
echo "2. 📊 Analyze top functions:"
echo "   cat $PROFILE_DIR/cpu_cycles_top_functions.txt"

echo ""
echo "3. 🔬 Identify optimization targets:"
echo "   - Look for wide flamegraph sections (high CPU time)"
echo "   - Check for excessive I/O syscalls in io_trace.txt"
echo "   - Examine cache miss patterns in cache_misses report"

echo ""
echo "4. 🎯 Common optimization opportunities:"
echo "   - Memory allocations (look for malloc/free patterns)"
echo "   - String operations (parsing, formatting)"
echo "   - File I/O patterns (seek operations, buffer sizes)"
echo "   - Algorithm efficiency (nested loops, redundant work)"

echo ""
echo "Profiling completed at: $(date)"
echo ""

echo "🚀 Ready to identify optimization opportunities!"