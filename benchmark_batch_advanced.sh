#!/bin/bash

# Advanced Batch Size Optimization Benchmark for inquiSTR
# This script provides comprehensive benchmarking with real-world considerations

set -euo pipefail

# Trap to handle script errors gracefully
trap 'echo "❌ Script failed at line $LINENO. Check your inputs and file paths." >&2; exit 1' ERR

# Script configuration
SCRIPT_DIR="$(pwd)"
BINARY="$HOME/repositories/inquiSTR/target/release/inquiSTR"
OUTPUT_DIR="${SCRIPT_DIR}/batch_benchmarks"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Default parameters (can be overridden)
BAM_FILE="${1:-test-data/small-test.bam}"
BED_FILE="${2:-test-data/test.bed}"
THREADS="${3:-2}"
RUNS_PER_SIZE="${4:-3}"  # Multiple runs for statistical significance

# Batch sizes to test (in KB) - expanded range
BATCH_SIZES=(20 35 50 65 80 100)

print_header() {
    echo "🔬 ADVANCED INQUISTR BATCH SIZE BENCHMARK"
    echo "========================================="
    echo "Generated: $(date)"
    echo "BAM file: $BAM_FILE"
    echo "BED file: $BED_FILE"
    echo "Threads: $THREADS"
    echo "Runs per size: $RUNS_PER_SIZE"
    echo "Testing batch sizes: ${BATCH_SIZES[*]} KB"
    echo "Execution: Randomized/interleaved to reduce system state bias"
    echo ""
}

setup_environment() {
    # Create output directory and temp directory
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$HOME/tmp"
    
    # Build release binary if needed
    if [[ ! -f "$BINARY" ]] || [[ src/ -nt "$BINARY" ]]; then
        echo "🔨 Building optimized release binary..."
        cargo build --release
        echo ""
    fi
    
    # Validate input files
    if [[ ! -f "$BAM_FILE" ]]; then
        echo "❌ ERROR: BAM file '$BAM_FILE' not found!"
        echo ""
        echo "Usage: $0 [BAM_FILE] [BED_FILE] [THREADS] [RUNS_PER_SIZE]"
        echo ""
        echo "Examples:"
        echo "  $0                                    # Use defaults"
        echo "  $0 my_sample.bam targets.bed 4 5     # Custom parameters"
        echo "  $0 large_file.bam big_targets.bed 8  # Production test"
        exit 1
    fi
    
    if [[ ! -f "$BED_FILE" ]]; then
        echo "❌ ERROR: BED file '$BED_FILE' not found!"
        exit 1
    fi
}

run_benchmark() {
    local batch_size=$1
    local run_number=$2
    local temp_output="$HOME/tmp/inquistr_${batch_size}kb_run${run_number}_${TIMESTAMP}.out"
    local timing_file="$HOME/tmp/timing_${batch_size}kb_run${run_number}_${TIMESTAMP}.txt"
    
    # Debug: Check if files exist
    if [[ ! -f "$BINARY" ]]; then
        echo "ERROR: Binary not found at $BINARY" >&2
        echo "FAILED,0,0,0,0"
        return
    fi
    
    if [[ ! -f "$BAM_FILE" ]]; then
        echo "ERROR: BAM file not found at $BAM_FILE" >&2
        echo "FAILED,0,0,0,0"
        return
    fi
    
    if [[ ! -f "$BED_FILE" ]]; then
        echo "ERROR: BED file not found at $BED_FILE" >&2
        echo "FAILED,0,0,0,0"
        return
    fi
    
    # Determine timing command
    local time_cmd
    if command -v /usr/bin/time >/dev/null 2>&1; then
        time_cmd="/usr/bin/time"
    elif command -v gtime >/dev/null 2>&1; then  # macOS with GNU time via homebrew
        time_cmd="gtime"
    else
        time_cmd="time"
    fi
    # Run benchmark with proper error handling
    set +e
    # Use a separate file for stderr to capture both timing and error messages
    local stderr_file="$HOME/tmp/stderr_${batch_size}kb_run${run_number}_${TIMESTAMP}.txt"
    
    "$time_cmd" -f "%e,%M,%P" "$BINARY" call "$BAM_FILE" \
        --region-file "$BED_FILE" \
        --threads "$THREADS" \
        --batch-size "$batch_size" \
        --unphased \
        > "$temp_output" 2> "$stderr_file"
    
    local exit_code=$?
    set -e
    
    # Check for errors and display them prominently
    if [[ $exit_code -ne 0 ]]; then
        echo ""
        echo "🚨 ERROR: inquiSTR failed with exit code $exit_code 🚨" >&2
        echo "📄 Command that failed:" >&2
        echo "   $BINARY call $BAM_FILE --region-file $BED_FILE --threads $THREADS --batch-size $batch_size --unphased" >&2
        echo "📋 Error output:" >&2
        echo "----------------------------------------" >&2
        # Show the error, filtering out timing info
        grep -v "^[0-9]*\.[0-9]*,[0-9]*,[0-9]*%" "$stderr_file" 2>/dev/null | head -20 >&2 || cat "$stderr_file" >&2
        echo "----------------------------------------" >&2
        echo "FAILED,0,0,0,0"
        rm -f "$temp_output" "$timing_file" "$stderr_file" 2>/dev/null
        return
    fi
    
    # Parse timing results (timing info should be the last line)
    if [[ -f "$stderr_file" ]]; then
        local timing_data
        timing_data=$(tail -1 "$stderr_file" 2>/dev/null || echo "0,0,0%")
        
        # Check if this line looks like timing data
        if [[ $timing_data =~ ^[0-9]*\.[0-9]*,[0-9]*,[0-9]*%$ ]]; then
            local wall_time max_memory cpu_percent
            wall_time=$(echo "$timing_data" | cut -d',' -f1 | tr -d '\n' | xargs)
            max_memory=$(echo "$timing_data" | cut -d',' -f2 | tr -d '\n' | xargs)
            cpu_percent=$(echo "$timing_data" | cut -d',' -f3 | tr -d '\n' | tr -d '%' | xargs)
            
            local output_lines
            output_lines=$(wc -l < "$temp_output" 2>/dev/null || echo "0")
            
            echo "SUCCESS,$wall_time,$max_memory,$cpu_percent,$output_lines"
        else
            echo "🚨 WARNING: Could not parse timing data from: $timing_data" >&2
            echo "FAILED,0,0,0,0"
        fi
    else
        echo "🚨 WARNING: No timing file found" >&2
        echo "FAILED,0,0,0,0"
    fi
    
    # Cleanup
    rm -f "$temp_output" "$timing_file" "$stderr_file"
}

run_interleaved_benchmarks() {
    local results_cache="/tmp/benchmark_results_${TIMESTAMP}.txt"
    true > "$results_cache"  # Clear the cache file
    
    # Create an interleaved execution plan
    local -a execution_plan=()
    
    # Create execution plan: randomize the order of (batch_size, run_number) pairs
    for ((run=1; run<=RUNS_PER_SIZE; run++)); do
        for batch_size in "${BATCH_SIZES[@]}"; do
            execution_plan+=("${batch_size},${run}")
        done
    done
    
    # Shuffle the execution plan using shuf if available, otherwise use a simple randomization
    local -a shuffled_plan=()
    if command -v shuf >/dev/null 2>&1; then
        readarray -t shuffled_plan < <(printf '%s\n' "${execution_plan[@]}" | shuf)
    else
        # Simple randomization without shuf
        shuffled_plan=("${execution_plan[@]}")
        for ((i=${#shuffled_plan[@]}-1; i>0; i--)); do
            local j=$((RANDOM % (i+1)))
            local temp="${shuffled_plan[i]}"
            shuffled_plan[i]="${shuffled_plan[j]}"
            shuffled_plan[j]="$temp"
        done
    fi
    
    local total_runs=${#shuffled_plan[@]}
    local current_run=0
    
    echo "🔀 Running ${total_runs} benchmarks in randomized order..."
    echo "📊 Progress:"
    
    # Execute benchmarks in randomized order
    for plan_item in "${shuffled_plan[@]}"; do
        current_run=$((current_run + 1))
        local batch_size run_number
        batch_size=$(echo "$plan_item" | cut -d',' -f1)
        run_number=$(echo "$plan_item" | cut -d',' -f2)
        printf "[%2d/%2d] Testing %3dKB (run %d)... " "$current_run" "$total_runs" "$batch_size" "$run_number"
        
        # Run benchmark with explicit error handling
        local result
        if result=$(run_benchmark "$batch_size" "$run_number"); then
            if [[ $result == SUCCESS,* ]]; then
                local time_val
                time_val=$(echo "$result" | cut -d',' -f2)
                printf "✅ %.3fs\n" "$time_val"
            else
                printf "❌ FAILED (invalid result: %s)\n" "$result"
                result="FAILED,0,0,0,0"
            fi
        else
            printf "❌ FAILED (benchmark execution error)\n"
            result="FAILED,0,0,0,0"
        fi
        
        # Cache the result
        echo "${batch_size},${result}" >> "$results_cache"
        
        # Add a small random delay to further reduce system state correlation
        local delay=$((RANDOM % 2 + 1))  # 1 to 2 seconds
        sleep "$delay"
    done
    
    echo ""
    echo "✅ Completed all randomized benchmark runs!"
    echo ""
}

calculate_statistics() {
    local batch_size=$1
    local -a times=()
    local -a memories=()
    local -a cpu_usages=()
    local success_count=0
    
    # Get results from the results cache that was populated by run_interleaved_benchmarks
    while IFS=',' read -r size result_data; do
        if [[ "$size" == "$batch_size" && "$result_data" == SUCCESS,* ]]; then
            ((success_count++))
            local time_val memory_val cpu_val
            time_val=$(echo "$result_data" | cut -d',' -f2)
            memory_val=$(echo "$result_data" | cut -d',' -f3)
            cpu_val=$(echo "$result_data" | cut -d',' -f4)
            
            times+=("$time_val")
            memories+=("$memory_val")
            cpu_usages+=("$cpu_val")
        fi
    done < "/tmp/benchmark_results_${TIMESTAMP}.txt"
    
    if [[ $success_count -gt 0 ]]; then
        # Calculate statistics using awk
        local avg_time min_time max_time
        local avg_memory min_memory max_memory
        local avg_cpu
        
        avg_time=$(printf '%s\n' "${times[@]}" | awk '{s+=$1} END {printf "%.3f", s/NR}')
        min_time=$(printf '%s\n' "${times[@]}" | awk 'NR==1{min=$1} {if($1<min) min=$1} END {printf "%.3f", min}')
        max_time=$(printf '%s\n' "${times[@]}" | awk 'NR==1{max=$1} {if($1>max) max=$1} END {printf "%.3f", max}')
        
        avg_memory=$(printf '%s\n' "${memories[@]}" | awk '{s+=$1} END {printf "%.0f", s/NR}')
        min_memory=$(printf '%s\n' "${memories[@]}" | awk 'NR==1{min=$1} {if($1<min) min=$1} END {printf "%.0f", min}')
        max_memory=$(printf '%s\n' "${memories[@]}" | awk 'NR==1{max=$1} {if($1>max) max=$1} END {printf "%.0f", max}')
        
        avg_cpu=$(printf '%s\n' "${cpu_usages[@]}" | awk '{s+=$1} END {printf "%.1f", s/NR}')
        
        echo "$batch_size,$avg_time,$min_time,$max_time,$avg_memory,$min_memory,$max_memory,$avg_cpu,$success_count"
    else
        echo "$batch_size,FAILED,FAILED,FAILED,FAILED,FAILED,FAILED,FAILED,0"
    fi
}

main() {
    print_header
    setup_environment
    
    local results_file="${OUTPUT_DIR}/batch_optimization_${TIMESTAMP}.csv"
    local summary_file="${OUTPUT_DIR}/batch_summary_${TIMESTAMP}.txt"
    
    # Initialize results file
    {
        echo "# inquiSTR Batch Size Optimization Results"
        echo "# Generated: $(date)"
        echo "# BAM: $BAM_FILE"
        echo "# BED: $BED_FILE"
        echo "# Threads: $THREADS"
        echo "# Runs per size: $RUNS_PER_SIZE"
        echo "BatchSize_KB,AvgTime_s,MinTime_s,MaxTime_s,AvgMemory_KB,MinMemory_KB,MaxMemory_KB,AvgCPU_percent,SuccessfulRuns"
    } > "$results_file"
    
    # Run all benchmarks in randomized/interleaved order
    run_interleaved_benchmarks
    
    echo "📊 ANALYZING RESULTS:"
    echo "====================="
    printf "%-12s %-12s %-15s %-12s %-10s\n" "Batch Size" "Avg Time" "Memory (KB)" "CPU %" "Status"
    echo "---------------------------------------------------------------------"
    
    # Calculate and display statistics for each batch size
    for batch_size in "${BATCH_SIZES[@]}"; do
        local stats
        stats=$(calculate_statistics "$batch_size")
        
        # Parse and display results
        if [[ $stats == *,FAILED,* ]]; then
            printf "%-12s %-12s %-15s %-12s %-10s\n" "${batch_size}KB" "FAILED" "FAILED" "FAILED" "❌"
        else
            local avg_time avg_memory avg_cpu success_count
            avg_time=$(echo "$stats" | cut -d',' -f2)
            avg_memory=$(echo "$stats" | cut -d',' -f5)
            avg_cpu=$(echo "$stats" | cut -d',' -f8)
            success_count=$(echo "$stats" | cut -d',' -f9)
            
            printf "%-12s %-12ss %-15s %-12s%% %-10s\n" "${batch_size}KB" "$avg_time" "$avg_memory" "$avg_cpu" "✅ ${success_count}/${RUNS_PER_SIZE}"
        fi
        
        # Save to results file
        echo "$stats" >> "$results_file"
    done
    
    # Generate summary and recommendations
    generate_summary "$results_file" "$summary_file"
    
    # Cleanup temporary files
    rm -f "/tmp/benchmark_results_${TIMESTAMP}.txt"
}

generate_summary() {
    local results_file=$1
    local summary_file=$2
    
    echo ""
    echo "📈 ANALYSIS & RECOMMENDATIONS:"
    echo "==============================="
    
    {
        echo "inquiSTR Batch Size Optimization Summary"
        echo "========================================"
        echo "Generated: $(date)"
        echo "Test configuration:"
        echo "  - BAM file: $BAM_FILE"
        echo "  - BED file: $BED_FILE"
        echo "  - Threads: $THREADS"
        echo "  - Runs per configuration: $RUNS_PER_SIZE"
        echo ""
    } > "$summary_file"
    
    if command -v awk >/dev/null 2>&1; then
        # Find optimal configurations
        echo "🏆 OPTIMAL CONFIGURATIONS:" | tee -a "$summary_file"
        echo "" | tee -a "$summary_file"
        
        # Fastest configuration
        local fastest
        fastest=$(awk -F',' 'NR>8 && $2!="FAILED" {print $2 "," $1}' "$results_file" | sort -n | head -1 | cut -d',' -f2)
        if [[ -n "$fastest" ]]; then
            local fastest_stats
            fastest_stats=$(grep "^$fastest," "$results_file")
            local f_time f_memory f_cpu
            f_time=$(echo "$fastest_stats" | cut -d',' -f2)
            f_memory=$(echo "$fastest_stats" | cut -d',' -f5)
            f_cpu=$(echo "$fastest_stats" | cut -d',' -f8)
            echo "🚀 Fastest: ${fastest}KB (${f_time}s avg, ${f_memory}KB memory, ${f_cpu}% CPU)" | tee -a "$summary_file"
        fi
        
        # Most memory efficient
        local most_efficient
        most_efficient=$(awk -F',' 'NR>8 && $5!="FAILED" {print $5 "," $1}' "$results_file" | sort -n | head -1 | cut -d',' -f2)
        if [[ -n "$most_efficient" ]]; then
            local eff_stats
            eff_stats=$(grep "^$most_efficient," "$results_file")
            local e_time e_memory e_cpu
            e_time=$(echo "$eff_stats" | cut -d',' -f2)
            e_memory=$(echo "$eff_stats" | cut -d',' -f5)
            e_cpu=$(echo "$eff_stats" | cut -d',' -f8)
            echo "💾 Most Memory Efficient: ${most_efficient}KB (${e_time}s avg, ${e_memory}KB memory, ${e_cpu}% CPU)" | tee -a "$summary_file"
        fi
        
        # Best balanced (considering both time and memory)
        echo "" | tee -a "$summary_file"
        echo "📊 Top 5 performing configurations:" | tee -a "$summary_file"
        awk -F',' 'NR>8 && $2!="FAILED" {
            score = $2 + ($5/1000)  # Simple scoring: time + memory/1000
            print score "," $1 "KB: " $2 "s, " $5 "KB memory"
        }' "$results_file" | sort -n | head -5 | cut -d',' -f2 | tee -a "$summary_file"
    fi
    
    echo "" | tee -a "$summary_file"
    echo "📊 Files generated:" | tee -a "$summary_file"
    echo "  - Detailed results: $results_file" | tee -a "$summary_file"
    echo "  - Summary: $summary_file" | tee -a "$summary_file"
    echo "" | tee -a "$summary_file"
    echo "🔧 Usage with optimized settings:" | tee -a "$summary_file"
    echo "   inquiSTR call your_file.bam --region-file targets.bed --batch-size OPTIMAL_SIZE --threads $THREADS" | tee -a "$summary_file"
    echo "" | tee -a "$summary_file"
    echo "💡 Performance tips:" | tee -a "$summary_file"
    echo "   - Test with your actual production data for final optimization" | tee -a "$summary_file"
    echo "   - Consider your available memory when choosing batch size" | tee -a "$summary_file"
    echo "   - Larger files may benefit from different optimal batch sizes" | tee -a "$summary_file"
    echo "   - SSD storage typically allows for larger batch sizes" | tee -a "$summary_file"
}

# Run the benchmark
main "$@"