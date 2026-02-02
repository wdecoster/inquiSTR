#!/usr/bin/env python3
"""
Analyze and visualize inquiSTR batch size benchmarking results.

Generates heatmaps and plots showing optimal parameter combinations
for different performance metrics.
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

def load_results(results_file):
    """Load benchmark results from TSV file."""
    df = pd.read_csv(results_file, sep='\t')
    
    # Filter out failed runs
    df = df[df['wall_time_s'] != 'FAILED'].copy()
    
    # Convert numeric columns
    numeric_cols = ['threads', 'batch_size', 'repeat', 'wall_time_s', 'user_time_s', 
                    'sys_time_s', 'cpu_percent', 'mem_max_kb', 'vol_ctx_switches', 
                    'invol_ctx_switches', 'fs_inputs']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    return df


def aggregate_results(df):
    """Compute mean and std dev for each thread/batch_size combination."""
    agg_df = df.groupby(['threads', 'batch_size']).agg({
        'wall_time_s': ['mean', 'std', 'min'],
        'user_time_s': ['mean'],
        'sys_time_s': ['mean'],
        'cpu_percent': ['mean'],
        'mem_max_kb': ['mean'],
        'vol_ctx_switches': ['mean'],
        'invol_ctx_switches': ['mean'],
        'fs_inputs': ['mean']
    }).reset_index()
    
    # Flatten column names
    agg_df.columns = ['_'.join(col).strip('_') for col in agg_df.columns.values]
    
    return agg_df


def create_heatmap(df, metric, title, output_file, vmin=None, vmax=None, cmap='RdYlGn_r', 
                    annotate=True, fmt='.1f'):
    """Create a heatmap for a specific metric."""
    # Pivot for heatmap
    pivot = df.pivot(index='batch_size', columns='threads', values=metric)
    
    # Sort for better visualization
    pivot = pivot.sort_index(ascending=False)
    
    plt.figure(figsize=(10, 8))
    
    # Create heatmap
    if annotate:
        sns.heatmap(pivot, annot=True, fmt=fmt, cmap=cmap, 
                    vmin=vmin, vmax=vmax, cbar_kws={'label': title})
    else:
        sns.heatmap(pivot, cmap=cmap, vmin=vmin, vmax=vmax, 
                    cbar_kws={'label': title})
    
    plt.title(f'{title} vs Threads and Batch Size', fontsize=14, fontweight='bold')
    plt.xlabel('Number of Threads', fontsize=12)
    plt.ylabel('Batch Size (KB)', fontsize=12)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")


def create_optimization_plot(df, output_file):
    """Create plot showing optimal batch_size for each thread count."""
    # Find best batch_size for each thread count (minimum wall time)
    best = df.loc[df.groupby('threads')['wall_time_s_mean'].idxmin()]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Optimal Configuration Analysis', fontsize=16, fontweight='bold')
    
    # 1. Optimal batch size by thread count
    ax = axes[0, 0]
    ax.plot(best['threads'], best['batch_size'], marker='o', linewidth=2, markersize=8)
    ax.set_xlabel('Number of Threads')
    ax.set_ylabel('Optimal Batch Size (KB)')
    ax.set_title('Optimal Batch Size vs Thread Count')
    ax.grid(True, alpha=0.3)
    for _, row in best.iterrows():
        ax.annotate(f"{row['batch_size']:.0f}KB", 
                    (row['threads'], row['batch_size']),
                    textcoords="offset points", xytext=(0,10), ha='center')
    
    # 2. Wall time improvement over baseline (1 thread, 50KB)
    ax = axes[0, 1]
    baseline_time = df[(df['threads'] == 1) & (df['batch_size'] == 50)]['wall_time_s_mean'].values[0]
    speedup = baseline_time / best['wall_time_s_mean']
    ax.plot(best['threads'], speedup, marker='o', linewidth=2, markersize=8, color='green')
    ax.axhline(y=1, color='red', linestyle='--', alpha=0.5, label='Baseline')
    ax.set_xlabel('Number of Threads')
    ax.set_ylabel('Speedup vs Baseline')
    ax.set_title('Speedup vs Baseline (1 thread, 50KB)')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # 3. CPU utilization efficiency (CPU% / threads)
    ax = axes[1, 0]
    efficiency = best['cpu_percent_mean'] / (best['threads'] * 100)
    ax.plot(best['threads'], efficiency, marker='o', linewidth=2, markersize=8, color='purple')
    ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, label='Perfect efficiency')
    ax.set_xlabel('Number of Threads')
    ax.set_ylabel('CPU Efficiency (CPU% / (threads × 100))')
    ax.set_title('Parallel Efficiency')
    ax.set_ylim([0, 1.1])
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    # 4. System overhead (sys_time / wall_time)
    ax = axes[1, 1]
    overhead = best['sys_time_s_mean'] / best['wall_time_s_mean']
    ax.plot(best['threads'], overhead, marker='o', linewidth=2, markersize=8, color='orange')
    ax.set_xlabel('Number of Threads')
    ax.set_ylabel('System Time Ratio (sys / wall)')
    ax.set_title('System Overhead')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"  Saved: {output_file}")
    
    return best


def print_recommendations(df, best):
    """Print performance recommendations based on analysis."""
    print("\n" + "="*70)
    print("PERFORMANCE RECOMMENDATIONS")
    print("="*70)
    
    print("\nOptimal batch_size by thread count:")
    print("-" * 50)
    for _, row in best.iterrows():
        threads = int(row['threads'])
        batch_size = int(row['batch_size'])
        wall_time = row['wall_time_s_mean']
        cpu_pct = row['cpu_percent_mean']
        sys_time = row['sys_time_s_mean']
        
        efficiency = cpu_pct / (threads * 100)
        
        print(f"\n  {threads:2d} threads: {batch_size:3d}KB batch_size")
        print(f"    Wall time:   {wall_time:6.1f}s")
        print(f"    CPU usage:   {cpu_pct:5.1f}%")
        print(f"    Efficiency:  {efficiency:5.1%}")
        print(f"    System time: {sys_time:6.1f}s")
    
    # General patterns
    print("\n" + "="*70)
    print("PATTERN ANALYSIS")
    print("="*70)
    
    # Does batch size decrease with more threads?
    if best['batch_size'].is_monotonic_decreasing:
        print("\n✓ Pattern: Smaller batch sizes are better with more threads")
        print("  → Confirm hypothesis that parallelism needs finer granularity")
    
    # Check if there's a sweet spot
    optimal_row = best.loc[best['wall_time_s_mean'].idxmin()]
    print(f"\n✓ Overall fastest: {int(optimal_row['threads'])} threads with {int(optimal_row['batch_size'])}KB batch_size")
    print(f"  → Wall time: {optimal_row['wall_time_s_mean']:.1f}s")
    print(f"  → CPU: {optimal_row['cpu_percent_mean']:.1f}%")
    
    # Diminishing returns check
    speedups = []
    baseline = df[(df['threads'] == 1)]['wall_time_s_mean'].min()
    for threads in sorted(df['threads'].unique()):
        fastest = df[df['threads'] == threads]['wall_time_s_mean'].min()
        speedups.append((threads, baseline / fastest))
    
    print("\n✓ Speedup by thread count:")
    for threads, speedup in speedups:
        print(f"  {threads:2d} threads: {speedup:5.2f}x")
    
    # Check for diminishing returns
    if len(speedups) >= 3:
        last_gain = speedups[-1][1] - speedups[-2][1]
        prev_gain = speedups[-2][1] - speedups[-3][1]
        if last_gain < prev_gain * 0.5:
            print(f"\n⚠ Warning: Diminishing returns observed at {speedups[-1][0]} threads")
            print(f"  → Consider using {speedups[-2][0]} threads for better efficiency")


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <results.tsv> [output_dir]")
        sys.exit(1)
    
    results_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else str(Path(results_file).parent / "analysis")
    
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    print("Loading benchmark results...")
    df = load_results(results_file)
    
    if df.empty:
        print("Error: No valid results found!")
        sys.exit(1)
    
    print(f"Loaded {len(df)} benchmark runs")
    print(f"  Thread counts: {sorted(df['threads'].unique())}")
    print(f"  Batch sizes: {sorted(df['batch_size'].unique())} KB")
    print(f"  Repeats per config: {df['repeat'].max()}")
    
    print("\nAggregating results...")
    agg_df = aggregate_results(df)
    
    # Save aggregated results
    agg_file = Path(output_dir) / "aggregated_results.tsv"
    agg_df.to_csv(agg_file, sep='\t', index=False)
    print(f"Saved aggregated results to: {agg_file}")
    
    print("\nGenerating visualizations...")
    
    # Create heatmaps for key metrics
    create_heatmap(agg_df, 'wall_time_s_mean', 'Wall Time (seconds)',
                   Path(output_dir) / 'heatmap_wall_time.png', cmap='RdYlGn_r')
    
    create_heatmap(agg_df, 'cpu_percent_mean', 'CPU Utilization (%)',
                   Path(output_dir) / 'heatmap_cpu_percent.png', cmap='RdYlGn')
    
    create_heatmap(agg_df, 'sys_time_s_mean', 'System Time (seconds)',
                   Path(output_dir) / 'heatmap_sys_time.png', cmap='YlOrRd')
    
    create_heatmap(agg_df, 'vol_ctx_switches_mean', 'Voluntary Context Switches',
                   Path(output_dir) / 'heatmap_ctx_switches.png', cmap='YlOrRd',
                   annotate=False)
    
    # Create optimization analysis
    best = create_optimization_plot(agg_df, Path(output_dir) / 'optimization_analysis.png')
    
    # Print recommendations
    print_recommendations(agg_df, best)
    
    print("\n" + "="*70)
    print("Analysis complete! Check the following files:")
    print(f"  - {output_dir}/heatmap_*.png")
    print(f"  - {output_dir}/optimization_analysis.png")
    print(f"  - {output_dir}/aggregated_results.tsv")
    print("="*70)


if __name__ == "__main__":
    main()
