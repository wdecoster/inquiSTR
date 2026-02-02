#!/usr/bin/env python3
"""
Analyze STR locus density from BED file to inform batch_size parameter tuning.

This script calculates inter-locus distances and visualizes the distribution
to help determine optimal batching strategies for inquiSTR.
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def read_bed(bed_file):
    """Read BED file and return DataFrame with chromosome, start, end."""
    df = pd.read_csv(
        bed_file,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
    )
    return df.sort_values(["chrom", "start"])


def calculate_distances(df):
    """Calculate inter-locus distances within each chromosome."""
    distances = []
    
    for chrom, group in df.groupby("chrom", sort=False):
        positions = group["start"].values
        if len(positions) > 1:
            # Distance from end of one locus to start of next
            ends = group["end"].values[:-1]
            starts = positions[1:]
            chrom_distances = starts - ends
            distances.extend(chrom_distances)
    
    return np.array(distances)


def print_statistics(distances):
    """Print summary statistics of inter-locus distances."""
    print("\n" + "="*70)
    print("INTER-LOCUS DISTANCE STATISTICS")
    print("="*70)
    print(f"Total loci: {len(distances) + 1:,}")
    print(f"Total inter-locus gaps: {len(distances):,}")
    print(f"\nDistance statistics (bp):")
    print(f"  Mean:       {np.mean(distances):>12,.1f}")
    print(f"  Median:     {np.median(distances):>12,.1f}")
    print(f"  Std Dev:    {np.std(distances):>12,.1f}")
    print(f"  Min:        {np.min(distances):>12,.0f}")
    print(f"  Max:        {np.max(distances):>12,.0f}")
    
    print(f"\nPercentiles:")
    for p in [10, 25, 50, 75, 90, 95, 99]:
        print(f"  {p:>2}th:      {np.percentile(distances, p):>12,.1f}")
    
    # Analysis for batch size recommendations
    print(f"\n" + "="*70)
    print("BATCH SIZE ANALYSIS")
    print("="*70)
    
    for threshold_kb in [10, 20, 50, 100]:
        threshold_bp = threshold_kb * 1000
        pct_below = (distances < threshold_bp).sum() / len(distances) * 100
        pct_above = 100 - pct_below
        print(f"\n{threshold_kb}KB batch_size ({threshold_bp:,}bp):")
        print(f"  {pct_below:5.1f}% of gaps are SMALLER (loci will be batched together)")
        print(f"  {pct_above:5.1f}% of gaps are LARGER  (loci will be in separate batches)")


def plot_distribution(distances, output_file):
    """Create visualization of distance distribution."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("STR Locus Density Analysis", fontsize=16, fontweight='bold')
    
    # 1. Histogram with log scale
    ax = axes[0, 0]
    bins = np.logspace(0, np.log10(distances.max()), 50)
    ax.hist(distances, bins=bins, edgecolor='black', alpha=0.7)
    ax.set_xscale('log')
    ax.set_xlabel('Inter-locus distance (bp, log scale)')
    ax.set_ylabel('Count')
    ax.set_title('Distance Distribution (log scale)')
    ax.grid(True, alpha=0.3)
    
    # Add batch size markers
    for threshold_kb, color in [(10, 'red'), (20, 'orange'), (50, 'green')]:
        ax.axvline(threshold_kb * 1000, color=color, linestyle='--', 
                   linewidth=2, label=f'{threshold_kb}KB', alpha=0.7)
    ax.legend()
    
    # 2. Cumulative distribution
    ax = axes[0, 1]
    sorted_dist = np.sort(distances)
    cumulative = np.arange(1, len(sorted_dist) + 1) / len(sorted_dist) * 100
    ax.plot(sorted_dist, cumulative, linewidth=2)
    ax.set_xscale('log')
    ax.set_xlabel('Inter-locus distance (bp, log scale)')
    ax.set_ylabel('Cumulative percentage')
    ax.set_title('Cumulative Distribution')
    ax.grid(True, alpha=0.3)
    
    # Add batch size markers
    for threshold_kb, color in [(10, 'red'), (20, 'orange'), (50, 'green')]:
        threshold_bp = threshold_kb * 1000
        pct = (distances < threshold_bp).sum() / len(distances) * 100
        ax.axvline(threshold_bp, color=color, linestyle='--', 
                   linewidth=2, label=f'{threshold_kb}KB ({pct:.1f}%)', alpha=0.7)
    ax.legend()
    
    # 3. Box plot showing distribution at different scales
    ax = axes[1, 0]
    data_for_boxplot = [
        distances[distances < 10000],
        distances[(distances >= 10000) & (distances < 50000)],
        distances[(distances >= 50000) & (distances < 100000)],
        distances[distances >= 100000],
    ]
    labels = ['<10KB', '10-50KB', '50-100KB', '>100KB']
    positions = range(1, 5)
    bp = ax.boxplot(data_for_boxplot, positions=positions, labels=labels,
                    patch_artist=True, showfliers=False)
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
    ax.set_ylabel('Inter-locus distance (bp)')
    ax.set_title('Distance Distribution by Range')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add counts
    for i, data in enumerate(data_for_boxplot, 1):
        ax.text(i, ax.get_ylim()[1] * 0.95, f'n={len(data):,}',
                ha='center', va='top', fontsize=9)
    
    # 4. Density plot showing batching impact
    ax = axes[1, 1]
    
    # Calculate expected batch counts for different thresholds
    batch_thresholds = range(5, 101, 5)  # 5KB to 100KB in 5KB steps
    batch_counts = []
    
    for threshold_kb in batch_thresholds:
        threshold_bp = threshold_kb * 1000
        # Approximate: gaps larger than threshold create new batches
        n_batches = (distances >= threshold_bp).sum() + 1
        batch_counts.append(n_batches)
    
    ax.plot(batch_thresholds, batch_counts, linewidth=2, marker='o', markersize=4)
    ax.set_xlabel('batch_size (KB)')
    ax.set_ylabel('Approximate number of batches')
    ax.set_title('Expected Parallelism vs Batch Size')
    ax.grid(True, alpha=0.3)
    
    # Highlight tested values
    for threshold_kb, color in [(10, 'red'), (20, 'orange'), (50, 'green')]:
        idx = (threshold_kb - 5) // 5
        ax.axvline(threshold_kb, color=color, linestyle='--', alpha=0.5)
        ax.plot(threshold_kb, batch_counts[idx], 'o', color=color, 
                markersize=10, label=f'{threshold_kb}KB')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n📊 Plot saved to: {output_file}")


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <bed_file> [output_plot.png]")
        sys.exit(1)
    
    bed_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "locus_density_analysis.png"
    
    print(f"Reading BED file: {bed_file}")
    df = read_bed(bed_file)
    
    print(f"Calculating inter-locus distances...")
    distances = calculate_distances(df)
    
    print_statistics(distances)
    
    print(f"\nGenerating plots...")
    plot_distribution(distances, output_file)
    
    print("\n" + "="*70)
    print("RECOMMENDATIONS")
    print("="*70)
    
    median_dist = np.median(distances)
    mean_dist = np.mean(distances)
    
    if median_dist < 20000:
        print("✓ DENSE catalog (median gap < 20KB)")
        print("  → Recommended batch_size: 20-50KB")
        print("    Rationale: Most loci are close together, batching saves I/O")
    elif median_dist < 50000:
        print("✓ MEDIUM density catalog (median gap 20-50KB)")
        print("  → Recommended batch_size: 10-20KB")
        print("    Rationale: Balance batching benefits with parallelism")
    else:
        print("✓ SPARSE catalog (median gap > 50KB)")
        print("  → Recommended batch_size: 10KB or less")
        print("    Rationale: Most loci will be fetched separately anyway")
    
    print(f"\n  Your test results:")
    print(f"    50KB: 172% CPU, 248s runtime")
    print(f"    20KB: 207% CPU, 207s runtime  ← 17% faster")
    print(f"    10KB: 519% CPU, 133s runtime  ← 46% faster, high overhead")


if __name__ == "__main__":
    main()
