/*!
 * Performance optimization module for inquiSTR
 *
 * Empirically determines optimal batch_size and thread count for a given system,
 * input file, and repeat catalog configuration.
 *
 * Key insights recognized by the optimizer:
 * - Optimal settings are system-dependent (CPU count, memory bandwidth)
 * - Optimal settings are dataset-dependent (catalog density, file format)
 * - Trade-offs exist between parallelism (small batches) and I/O efficiency (large batches)
 * - System overhead (context switches, coordination) increases with thread count
 * - Diminishing returns occur beyond a certain thread count
 */

use crate::repeats::{RepeatInterval, TargetConfig};
use plotly::{
    HeatMap, Plot, Scatter,
    common::Title,
    layout::{Axis, GridPattern, Layout, LayoutGrid},
};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

/// Results from a single benchmark run
#[derive(Debug, Clone)]
pub struct BenchmarkResult {
    pub threads: usize,
    pub batch_size: u32,
    pub repeat: usize,
    pub wall_time: Duration,
    pub success: bool,
}

/// Aggregated statistics for a configuration
#[derive(Debug, Clone)]
pub struct ConfigStats {
    pub threads: usize,
    pub batch_size: u32,
    pub mean_wall_time: f64,
    pub std_wall_time: f64,
    pub min_wall_time: f64,
    pub success_rate: f64,
}

/// Recommendation for optimal configuration
#[derive(Debug, Clone)]
pub struct Recommendation {
    pub threads: usize,
    pub batch_size: u32,
    pub mean_wall_time: f64,
    pub rationale: String,
    pub confidence: Confidence,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Confidence {
    High,   // Clear winner, statistically significant
    Medium, // Best option but close competitors
    Low,    // Multiple similar options, marginal differences
}

/// Main optimization function
#[allow(clippy::too_many_arguments)]
pub fn optimize_parameters(
    bam: PathBuf,
    region_config: TargetConfig,
    reference: Option<String>,
    min_threads: usize,
    max_threads: usize,
    batch_sizes: Vec<u32>,
    repeats: usize,
    output_dir: &Path,
) -> Result<Recommendation, String> {
    eprintln!("===========================================");
    eprintln!("inquiSTR Performance Optimization");
    eprintln!("===========================================\n");

    // Detect system capabilities
    let system_info = detect_system_info();
    eprintln!("System: {} CPUs, {:.1} GB RAM", system_info.cpu_count, system_info.memory_gb);

    // Load target regions
    let (target_regions, chrom_mapper) = load_regions(&region_config, &bam, &reference)?;
    eprintln!("Testing with {} target regions", target_regions.len());
    eprintln!("Tip: To test faster, use a BED file with only one chromosome\n");

    // Validate thread range
    let max_threads = max_threads.min(system_info.cpu_count);
    if min_threads > max_threads {
        return Err(format!("min_threads ({}) > max_threads ({})", min_threads, max_threads));
    }

    // Generate thread counts to test (exponential spacing)
    let thread_counts = generate_thread_sequence(min_threads, max_threads);
    eprintln!("Thread counts: {:?}", thread_counts);
    eprintln!("Batch sizes: {:?} KB", batch_sizes);
    eprintln!("Repeats per config: {}\n", repeats);

    let total_tests = thread_counts.len() * batch_sizes.len() * repeats;
    eprintln!("Total tests: {}\n", total_tests);

    // Run benchmarks
    let results = run_benchmarks(
        &bam,
        &target_regions,
        &chrom_mapper,
        &reference,
        &thread_counts,
        &batch_sizes,
        repeats,
    )?;

    // Aggregate results
    let stats = aggregate_results(&results);

    // Generate visualizations
    std::fs::create_dir_all(output_dir)
        .map_err(|e| format!("Failed to create output dir: {}", e))?;
    generate_visualizations(&stats, output_dir)?;

    // Determine optimal configuration
    let recommendation = determine_recommendation(&stats, &system_info)?;

    // Print summary
    print_summary(&stats, &recommendation);

    // Save detailed results
    save_results(&results, &stats, &recommendation, output_dir)?;

    Ok(recommendation)
}

/// Detect system information
#[derive(Debug)]
struct SystemInfo {
    cpu_count: usize,
    memory_gb: f64,
}

fn detect_system_info() -> SystemInfo {
    let cpu_count = num_cpus::get();

    // Attempt to get memory info (Linux-specific)
    let memory_gb = if cfg!(target_os = "linux") {
        std::fs::read_to_string("/proc/meminfo")
            .ok()
            .and_then(|s| {
                s.lines()
                    .find(|line| line.starts_with("MemTotal:"))
                    .and_then(|line| {
                        line.split_whitespace()
                            .nth(1)
                            .and_then(|kb| kb.parse::<f64>().ok())
                            .map(|kb| kb / 1024.0 / 1024.0)
                    })
            })
            .unwrap_or(0.0)
    } else {
        0.0
    };

    SystemInfo { cpu_count, memory_gb }
}

/// Generate exponential sequence of thread counts
fn generate_thread_sequence(min: usize, max: usize) -> Vec<usize> {
    let mut threads = vec![min];
    let mut current = min;

    while current < max {
        current = if current == 1 {
            2
        } else {
            (current * 2).min(max)
        };
        threads.push(current);
    }

    // Ensure max is included if not already
    if threads.last() != Some(&max) {
        threads.push(max);
    }

    threads
}

/// Load target regions from configuration
fn load_regions(
    config: &TargetConfig,
    bam: &Path,
    reference: &Option<String>,
) -> Result<(Vec<RepeatInterval>, crate::repeats::ChromosomeMapper), String> {
    config
        .get_targets(bam.to_str().unwrap(), reference)
        .map_err(|e| format!("Failed to load regions: {:?}", e))
}

/// Run all benchmark combinations
fn run_benchmarks(
    bam: &Path,
    regions: &[RepeatInterval],
    chrom_mapper: &crate::repeats::ChromosomeMapper,
    reference: &Option<String>,
    thread_counts: &[usize],
    batch_sizes: &[u32],
    repeats: usize,
) -> Result<Vec<BenchmarkResult>, String> {
    let mut results = Vec::new();
    let total = thread_counts.len() * batch_sizes.len() * repeats;
    let mut current = 0;

    for &threads in thread_counts {
        for &batch_size in batch_sizes {
            for repeat in 1..=repeats {
                current += 1;
                eprint!(
                    "[{}/{}] Testing: {} threads, {}KB batch, repeat {}/{}... ",
                    current, total, threads, batch_size, repeat, repeats
                );

                let result = run_single_benchmark(
                    bam,
                    regions,
                    chrom_mapper,
                    reference,
                    threads,
                    batch_size,
                )?;

                eprintln!("{:.1}s", result.wall_time.as_secs_f64());
                results.push(result);
            }
        }
    }

    Ok(results)
}

/// Run a single benchmark configuration
fn run_single_benchmark(
    bam: &Path,
    regions: &[RepeatInterval],
    chrom_mapper: &crate::repeats::ChromosomeMapper,
    reference: &Option<String>,
    threads: usize,
    batch_size: u32,
) -> Result<BenchmarkResult, String> {
    use std::io::Write;
    use tempfile::NamedTempFile;

    // Write regions to a temporary BED file
    // Need to write in BED format: chrom, start, end, name (optional), score, strand, motif (optional)
    let mut temp_bed =
        NamedTempFile::new().map_err(|e| format!("Failed to create temp BED file: {}", e))?;

    for region in regions {
        let chrom_name = chrom_mapper.get_name(region.chrom_id);
        writeln!(temp_bed, "{}\t{}\t{}\t{}", chrom_name, region.start, region.end, region.info)
            .map_err(|e| format!("Failed to write to temp BED: {}", e))?;
    }
    temp_bed
        .flush()
        .map_err(|e| format!("Failed to flush temp BED: {}", e))?;

    let start = Instant::now();

    // Suppress stdout during benchmarking (TSV output goes to stdout by default)
    let _stdout_redirect =
        gag::BufferRedirect::stdout().map_err(|e| format!("Failed to redirect stdout: {}", e))?;

    // Run genotyping with current configuration
    let result = crate::call::genotype_repeats(
        bam.to_str().unwrap().to_string(),
        crate::repeats::TargetConfig {
            region: None,
            region_file: Some(temp_bed.path().to_path_buf()),
            preset: None,
            max_locus: None,
        },
        crate::call::GenotypeConfig { minlen: 1, support: 3, unphased: false },
        crate::call::ProcessingConfig { threads, batch_size_kb: batch_size, output_vcf: false },
        None, // sample_name
        reference.clone(),
        false, // show_progress
    );

    // Drop stdout redirect to restore normal output
    drop(_stdout_redirect);

    let wall_time = start.elapsed();

    // If benchmark failed, return the error
    if let Err(e) = result {
        return Err(format!(
            "Benchmark failed for threads={}, batch_size={}KB: {:?}",
            threads, batch_size, e
        ));
    }

    Ok(BenchmarkResult {
        threads,
        batch_size,
        repeat: 0, // Will be set by caller
        wall_time,
        success: true,
    })
}

/// Aggregate benchmark results by configuration
fn aggregate_results(results: &[BenchmarkResult]) -> Vec<ConfigStats> {
    let mut grouped: HashMap<(usize, u32), Vec<&BenchmarkResult>> = HashMap::new();

    for result in results {
        grouped
            .entry((result.threads, result.batch_size))
            .or_default()
            .push(result);
    }

    let mut stats = Vec::new();
    for ((threads, batch_size), group) in grouped {
        let times: Vec<f64> = group.iter().map(|r| r.wall_time.as_secs_f64()).collect();
        let successes = group.iter().filter(|r| r.success).count();

        let mean = times.iter().sum::<f64>() / times.len() as f64;
        let variance = times.iter().map(|t| (t - mean).powi(2)).sum::<f64>() / times.len() as f64;
        let std = variance.sqrt();
        let min = times.iter().cloned().fold(f64::INFINITY, f64::min);

        stats.push(ConfigStats {
            threads,
            batch_size,
            mean_wall_time: mean,
            std_wall_time: std,
            min_wall_time: min,
            success_rate: successes as f64 / group.len() as f64,
        });
    }

    stats.sort_by(|a, b| {
        a.threads
            .cmp(&b.threads)
            .then(a.batch_size.cmp(&b.batch_size))
    });

    stats
}

/// Generate plotly visualizations
fn generate_visualizations(stats: &[ConfigStats], output_dir: &Path) -> Result<(), String> {
    // Create heatmap of wall times
    create_wall_time_heatmap(stats, output_dir)?;

    // Create optimization analysis plots
    create_optimization_plots(stats, output_dir)?;

    Ok(())
}

/// Create wall time heatmap
fn create_wall_time_heatmap(stats: &[ConfigStats], output_dir: &Path) -> Result<(), String> {
    let unique_threads: Vec<usize> = {
        let mut threads: Vec<_> = stats.iter().map(|s| s.threads).collect();
        threads.sort_unstable();
        threads.dedup();
        threads
    };

    let unique_batch_sizes: Vec<u32> = {
        let mut sizes: Vec<_> = stats.iter().map(|s| s.batch_size).collect();
        sizes.sort_unstable();
        sizes.dedup();
        sizes
    };

    // Build z-matrix for heatmap
    let mut z_matrix = vec![vec![f64::NAN; unique_threads.len()]; unique_batch_sizes.len()];

    for stat in stats {
        if let Some(row) = unique_batch_sizes
            .iter()
            .position(|&bs| bs == stat.batch_size)
            && let Some(col) = unique_threads.iter().position(|&t| t == stat.threads)
        {
            z_matrix[row][col] = stat.mean_wall_time;
        }
    }

    let heatmap = HeatMap::new_z(z_matrix)
        .x(unique_threads.iter().map(|t| *t as f64).collect())
        .y(unique_batch_sizes.iter().map(|bs| *bs as f64).collect())
        .color_scale(plotly::common::ColorScale::Palette(
            plotly::common::ColorScalePalette::Viridis,
        ));

    let layout = Layout::new()
        .title(Title::from("Wall Time (seconds) vs Threads and Batch Size"))
        .x_axis(Axis::new().title(Title::from("Number of Threads")))
        .y_axis(Axis::new().title(Title::from("Batch Size")));

    let mut plot = Plot::new();
    plot.add_trace(heatmap);
    plot.set_layout(layout);

    let output_file = output_dir.join("wall_time_heatmap.html");
    plot.write_html(&output_file);

    eprintln!("Saved: {}", output_file.display());
    Ok(())
}

/// Create optimization analysis plots (4-panel figure)
fn create_optimization_plots(stats: &[ConfigStats], output_dir: &Path) -> Result<(), String> {
    // Find optimal batch_size for each thread count
    let mut best_by_threads: HashMap<usize, &ConfigStats> = HashMap::new();

    for stat in stats {
        best_by_threads
            .entry(stat.threads)
            .and_modify(|current| {
                if stat.mean_wall_time < current.mean_wall_time {
                    *current = stat;
                }
            })
            .or_insert(stat);
    }

    let mut sorted_best: Vec<_> = best_by_threads.into_iter().collect();
    sorted_best.sort_by_key(|(threads, _)| *threads);

    // Extract data for plots
    let threads: Vec<usize> = sorted_best.iter().map(|(t, _)| *t).collect();
    let batch_sizes: Vec<u32> = sorted_best.iter().map(|(_, s)| s.batch_size).collect();
    let wall_times: Vec<f64> = sorted_best.iter().map(|(_, s)| s.mean_wall_time).collect();

    // Calculate speedup (vs 1 thread baseline)
    let baseline = wall_times[0];
    let speedups: Vec<f64> = wall_times.iter().map(|t| baseline / t).collect();

    // Create 4-panel plot using subplot grid
    let mut plot = Plot::new();

    // Plot 1: Optimal batch size vs threads
    let trace1 = Scatter::new(threads.clone(), batch_sizes)
        .mode(plotly::common::Mode::LinesMarkers)
        .name("Optimal Batch Size")
        .x_axis("x1")
        .y_axis("y1");
    plot.add_trace(trace1);

    // Plot 2: Speedup vs threads
    let trace2 = Scatter::new(threads.clone(), speedups.clone())
        .mode(plotly::common::Mode::LinesMarkers)
        .name("Speedup")
        .x_axis("x2")
        .y_axis("y2");
    plot.add_trace(trace2);

    // Calculate efficiency (speedup / threads)
    let efficiency: Vec<f64> = threads
        .iter()
        .zip(speedups.iter())
        .map(|(t, s)| s / (*t as f64))
        .collect();

    // Plot 3: Parallel efficiency
    let trace3 = Scatter::new(threads.clone(), efficiency)
        .mode(plotly::common::Mode::LinesMarkers)
        .name("Efficiency")
        .x_axis("x3")
        .y_axis("y3");
    plot.add_trace(trace3);

    // Plot 4: Wall time vs threads
    let trace4 = Scatter::new(threads, wall_times)
        .mode(plotly::common::Mode::LinesMarkers)
        .name("Wall Time")
        .x_axis("x4")
        .y_axis("y4");
    plot.add_trace(trace4);

    let layout = Layout::new()
        .title(Title::from("Optimization Analysis"))
        .grid(
            LayoutGrid::new()
                .rows(2)
                .columns(2)
                .pattern(GridPattern::Independent),
        )
        .x_axis(
            Axis::new()
                .title(Title::from("Threads"))
                .domain(&[0.0, 0.45]),
        )
        .y_axis(Axis::new().title(Title::from("Batch Size (KB)")))
        .x_axis2(
            Axis::new()
                .title(Title::from("Threads"))
                .domain(&[0.55, 1.0]),
        )
        .y_axis2(Axis::new().title(Title::from("Speedup")).anchor("x2"))
        .x_axis3(
            Axis::new()
                .title(Title::from("Threads"))
                .domain(&[0.0, 0.45]),
        )
        .y_axis3(Axis::new().title(Title::from("Efficiency")).anchor("x3"))
        .x_axis4(
            Axis::new()
                .title(Title::from("Threads"))
                .domain(&[0.55, 1.0]),
        )
        .y_axis4(Axis::new().title(Title::from("Wall Time (s)")).anchor("x4"));

    plot.set_layout(layout);

    let output_file = output_dir.join("optimization_analysis.html");
    plot.write_html(&output_file);

    eprintln!("Saved: {}", output_file.display());
    Ok(())
}

/// Determine optimal configuration with rationale
fn determine_recommendation(
    stats: &[ConfigStats],
    system_info: &SystemInfo,
) -> Result<Recommendation, String> {
    if stats.is_empty() {
        return Err("No successful benchmark results".to_string());
    }

    // Filter successful runs only
    let successful: Vec<_> = stats.iter().filter(|s| s.success_rate > 0.9).collect();

    if successful.is_empty() {
        return Err("No successful configurations found".to_string());
    }

    // Find absolute minimum wall time
    let best = successful
        .iter()
        .min_by(|a, b| a.mean_wall_time.partial_cmp(&b.mean_wall_time).unwrap())
        .unwrap();

    // Check if winner is clear (>5% faster than second best)
    let second_best = successful
        .iter()
        .filter(|s| s.threads != best.threads || s.batch_size != best.batch_size)
        .min_by(|a, b| a.mean_wall_time.partial_cmp(&b.mean_wall_time).unwrap());

    let confidence = if let Some(second) = second_best {
        let improvement = (second.mean_wall_time - best.mean_wall_time) / second.mean_wall_time;
        if improvement > 0.05 {
            Confidence::High
        } else if improvement > 0.02 {
            Confidence::Medium
        } else {
            Confidence::Low
        }
    } else {
        Confidence::High
    };

    // Build rationale
    let mut rationale = format!("Fastest configuration: {:.1}s wall time. ", best.mean_wall_time);

    // Check if using all CPUs
    if best.threads < system_info.cpu_count {
        rationale.push_str(&format!(
            "Note: Uses {}/{} CPUs. More threads showed diminishing returns. ",
            best.threads, system_info.cpu_count
        ));
    }

    // Comment on batch size
    if best.batch_size <= 10 {
        rationale.push_str("Small batch size provides fine-grained parallelism. ");
    } else if best.batch_size >= 50 {
        rationale.push_str("Large batch size minimizes I/O overhead. ");
    } else {
        rationale.push_str("Medium batch size balances parallelism and I/O efficiency. ");
    }

    Ok(Recommendation {
        threads: best.threads,
        batch_size: best.batch_size,
        mean_wall_time: best.mean_wall_time,
        rationale,
        confidence,
    })
}

/// Print summary to stderr
fn print_summary(stats: &[ConfigStats], recommendation: &Recommendation) {
    eprintln!("\n===========================================");
    eprintln!("OPTIMIZATION RESULTS");
    eprintln!("===========================================\n");

    eprintln!("🎯 RECOMMENDED CONFIGURATION:");
    eprintln!(
        "   --threads {}  --batch-size {}",
        recommendation.threads, recommendation.batch_size
    );
    eprintln!("   Expected runtime: {:.1}s", recommendation.mean_wall_time);
    eprintln!("   Confidence: {:?}", recommendation.confidence);
    eprintln!("\n   {}", recommendation.rationale);

    // Show top 5 configurations
    let mut sorted_stats = stats.to_vec();
    sorted_stats.sort_by(|a, b| a.mean_wall_time.partial_cmp(&b.mean_wall_time).unwrap());

    eprintln!("\n📊 TOP 5 CONFIGURATIONS:");
    eprintln!("   Rank  Threads  Batch   Time(s)");
    eprintln!("   ----  -------  -----   -------");
    for (i, stat) in sorted_stats.iter().take(5).enumerate() {
        eprintln!(
            "   {:4}  {:7}  {:3}KB   {:7.1}",
            i + 1,
            stat.threads,
            stat.batch_size,
            stat.mean_wall_time
        );
    }

    eprintln!("\n===========================================");
}

/// Save detailed results to files
fn save_results(
    results: &[BenchmarkResult],
    stats: &[ConfigStats],
    recommendation: &Recommendation,
    output_dir: &Path,
) -> Result<(), String> {
    // Save raw results
    let raw_file = output_dir.join("benchmark_results.tsv");
    let mut raw_content = String::from("threads\tbatch_size\trepeat\twall_time_s\tsuccess\n");
    for r in results {
        raw_content.push_str(&format!(
            "{}\t{}\t{}\t{:.3}\t{}\n",
            r.threads,
            r.batch_size,
            r.repeat,
            r.wall_time.as_secs_f64(),
            r.success
        ));
    }
    std::fs::write(&raw_file, raw_content)
        .map_err(|e| format!("Failed to write results: {}", e))?;

    // Save aggregated stats
    let stats_file = output_dir.join("aggregated_stats.tsv");
    let mut stats_content =
        String::from("threads\tbatch_size\tmean_time_s\tstd_time_s\tmin_time_s\tsuccess_rate\n");
    for s in stats {
        stats_content.push_str(&format!(
            "{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\n",
            s.threads,
            s.batch_size,
            s.mean_wall_time,
            s.std_wall_time,
            s.min_wall_time,
            s.success_rate
        ));
    }
    std::fs::write(&stats_file, stats_content)
        .map_err(|e| format!("Failed to write stats: {}", e))?;

    // Save recommendation
    let rec_file = output_dir.join("recommendation.txt");
    let rec_content = format!(
        "RECOMMENDED CONFIGURATION\n\
         ========================\n\n\
         Threads:     {}\n\
         Batch Size:  {}KB\n\
         Wall Time:   {:.1}s\n\
         Confidence:  {:?}\n\n\
         Rationale:\n{}\n\n\
         To use this configuration:\n\
         inquiSTR call <input> -R <regions> --threads {} --batch-size {}\n",
        recommendation.threads,
        recommendation.batch_size,
        recommendation.mean_wall_time,
        recommendation.confidence,
        recommendation.rationale,
        recommendation.threads,
        recommendation.batch_size,
    );
    std::fs::write(&rec_file, rec_content)
        .map_err(|e| format!("Failed to write recommendation: {}", e))?;

    eprintln!("\nResults saved to: {}", output_dir.display());

    Ok(())
}
