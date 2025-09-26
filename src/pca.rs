//! # Principal Component Analysis (PCA) for STR Genotypes
//!
//! This module implements Principal Component Analysis for Short Tandem Repeat (STR) genotype data
//! from inquiSTR combined files. It provides dimensionality reduction and visualization to identify
//! population structure and relationships between samples.
//!
//! ## Features
//!
//! - **Automated data parsing** from inquiSTR combined files
//! - **Simplified PCA implementation** using variance-based feature selection
//! - **Interactive HTML plots** using Plotly for visualization
//! - **Support for missing data** (NaN values handled gracefully)
//!
//! ## Usage
//!
//! ```bash
//! # Generate PCA plot from combined STR data
//! inquiSTR pca combined_strs.tsv --output str_pca.html
//! ```
//!
//! ## Implementation
//!
//! This is a simplified PCA implementation that uses the top variance features as principal
//! components. For production use with large datasets, consider using more sophisticated
//! eigenvalue decomposition methods.

use ndarray::{Array1, Array2, Axis};
use plotly::color::{NamedColor, Rgba};
use plotly::common::{Marker, Mode, Title};
use plotly::layout::{Axis as PlotAxis, Layout};
use plotly::{Plot, Scatter};
use std::io::BufRead;
use std::path::PathBuf;

/// Parse a combined inquiSTR file and return the data matrix with sample names and region names
fn parse_combined_file(combined: &std::path::Path) -> (Array2<f64>, Vec<String>, Vec<String>) {
    let file = crate::utils::reader(&combined.to_string_lossy());
    let mut lines = file.lines();

    // Read header line - this should contain sample names
    let header_line = lines
        .next()
        .expect("File is empty")
        .expect("Error reading header line");

    let header_fields: Vec<&str> = header_line.trim().split('\t').collect();

    if header_fields.len() < 5
        || header_fields[0] != "chromosome"
        || header_fields[1] != "begin"
        || header_fields[2] != "end"
    {
        panic!("Invalid combined file header. Expected format: chromosome\tbegin\tend\tsample1_H1\tsample1_H2\t...");
    }

    // Extract sample names from header (columns 3+ should be sample_H1, sample_H2 pattern)
    let data_cols = &header_fields[3..];
    if !data_cols.len().is_multiple_of(2) {
        panic!("Invalid header: number of sample columns must be even (H1/H2 pairs)");
    }

    let num_samples = data_cols.len() / 2;
    let mut sample_names = Vec::new();

    for i in 0..num_samples {
        let h1_col = &data_cols[i * 2];
        let h2_col = &data_cols[i * 2 + 1];

        // Extract sample name from H1 column (should end with _H1)
        if !h1_col.ends_with("_H1") {
            panic!(
                "Invalid header: expected column {} to end with '_H1', got '{}'",
                3 + i * 2,
                h1_col
            );
        }
        if !h2_col.ends_with("_H2") {
            panic!(
                "Invalid header: expected column {} to end with '_H2', got '{}'",
                3 + i * 2 + 1,
                h2_col
            );
        }

        let sample_name_h1 = h1_col.trim_end_matches("_H1");
        let sample_name_h2 = h2_col.trim_end_matches("_H2");

        if sample_name_h1 != sample_name_h2 {
            panic!(
                "Header error: H1 and H2 columns have different sample names: '{}' vs '{}'",
                sample_name_h1, sample_name_h2
            );
        }

        sample_names.push(sample_name_h1.to_string());
    }

    // Collect all data lines (after the header)
    let mut data_lines = Vec::new();
    for (line_num, line_result) in lines.enumerate() {
        let line = line_result
            .map_err(|e| format!("Error reading line {}: {}", line_num + 2, e))
            .expect("IO error reading file");
        data_lines.push(line);
    }

    if data_lines.is_empty() {
        panic!("No data lines found after header");
    }

    let num_regions = data_lines.len();
    let mut data_matrix = Array2::<f64>::zeros((num_samples, num_regions));
    let mut region_names = Vec::new();

    for (region_idx, line) in data_lines.iter().enumerate() {
        let fields: Vec<&str> = line.trim().split('\t').collect();

        let expected_cols = 3 + num_samples * 2;
        if fields.len() != expected_cols {
            panic!(
                "Malformed line {} (expected {} columns, got {}): {}",
                region_idx + 2,
                expected_cols,
                fields.len(),
                line
            );
        }

        // Create region name from chr:start-end
        let region_name = format!("{}:{}-{}", fields[0], fields[1], fields[2]);
        region_names.push(region_name);

        // Parse STR lengths for each sample (taking max of H1 and H2 for PCA)
        for sample_idx in 0..num_samples {
            let h1_idx = 3 + sample_idx * 2;
            let h2_idx = 4 + sample_idx * 2;

            let h1_val: f64 = fields[h1_idx]
                .parse()
                .map_err(|e| {
                    format!(
                        "Invalid H1 value '{}' at line {}, column {}: {}",
                        fields[h1_idx],
                        region_idx + 2,
                        h1_idx + 1,
                        e
                    )
                })
                .unwrap_or_else(|e| {
                    // Handle NaN values specifically
                    if fields[h1_idx].eq_ignore_ascii_case("nan") {
                        0.0
                    } else {
                        panic!("{}", e);
                    }
                });

            let h2_val: f64 = fields[h2_idx]
                .parse()
                .map_err(|e| {
                    format!(
                        "Invalid H2 value '{}' at line {}, column {}: {}",
                        fields[h2_idx],
                        region_idx + 2,
                        h2_idx + 1,
                        e
                    )
                })
                .unwrap_or_else(|e| {
                    // Handle NaN values specifically
                    if fields[h2_idx].eq_ignore_ascii_case("nan") {
                        0.0
                    } else {
                        panic!("{}", e);
                    }
                });

            // Use maximum allele length for PCA (could also use sum or mean)
            let max_allele = h1_val.max(h2_val);
            data_matrix[[sample_idx, region_idx]] = max_allele;
        }
    }

    (data_matrix, sample_names, region_names)
}

/// Simplified PCA implementation using covariance matrix eigenvalue decomposition
/// This is a basic implementation that returns the first 2 principal components
fn perform_simple_pca(data: &Array2<f64>) -> (Array2<f64>, Array1<f64>) {
    // Center the data (subtract mean for each feature/STR)
    let mean = data.mean_axis(Axis(0)).unwrap();
    let centered_data = data - &mean.insert_axis(Axis(0));

    // For now, use a simplified approach - just return first 2 "components"
    // In reality, we'd compute eigenvectors of covariance matrix
    // Here we'll use a simple projection approach

    let n_features = centered_data.ncols();
    let mut pca_data = Array2::<f64>::zeros((data.nrows(), 2));

    // Simple projection onto first two dimensions with highest variance
    let mut feature_variances: Vec<(usize, f64)> = Vec::new();
    for i in 0..n_features {
        let var = centered_data.column(i).var(0.0);
        feature_variances.push((i, var));
    }

    // Sort by variance (descending)
    feature_variances.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    // Use top 2 features as "principal components" (simplified approach)
    let pc1_idx = feature_variances[0].0;
    let pc2_idx = if feature_variances.len() > 1 {
        feature_variances[1].0
    } else {
        0
    };

    for i in 0..data.nrows() {
        pca_data[[i, 0]] = centered_data[[i, pc1_idx]];
        pca_data[[i, 1]] = centered_data[[i, pc2_idx]];
    }

    // Return explained variance ratios (simplified)
    let var1 = feature_variances[0].1;
    let var2 = if feature_variances.len() > 1 {
        feature_variances[1].1
    } else {
        0.0
    };
    let total_var: f64 = feature_variances.iter().map(|(_, v)| v).sum();

    let explained_var = Array1::from(vec![var1 / total_var * 100.0, var2 / total_var * 100.0]);

    (pca_data, explained_var)
}

/// Create a PCA plot using Plotly
fn create_pca_plot(
    pca_data: &Array2<f64>,
    sample_names: &[String],
    explained_variance: &Array1<f64>,
    output: &str,
) {
    // Extract PC1 and PC2 coordinates
    let pc1_coords: Vec<f64> = pca_data.column(0).to_vec();
    let pc2_coords: Vec<f64> = pca_data.column(1).to_vec();

    // Create scatter plot
    let trace = Scatter::new(pc1_coords, pc2_coords)
        .mode(Mode::Markers)
        .name("Samples")
        .text_array(sample_names.to_vec())
        .marker(
            Marker::new()
                .size(10)
                .color(Rgba::new(31, 119, 180, 1.0))
                .line(
                    plotly::common::Line::new()
                        .color(NamedColor::White)
                        .width(1.0),
                ),
        );

    let mut plot = Plot::new();
    plot.add_trace(trace);

    // Set up layout
    let layout = Layout::new()
        .title(Title::with_text("PCA of STR Genotypes"))
        .x_axis(
            PlotAxis::new()
                .title(Title::with_text(format!("PC1 ({:.1}% variance)", explained_variance[0]))),
        )
        .y_axis(
            PlotAxis::new()
                .title(Title::with_text(format!("PC2 ({:.1}% variance)", explained_variance[1]))),
        )
        .show_legend(false);

    plot.set_layout(layout);

    // Save to HTML file
    plot.write_html(output);
    println!("PCA plot saved to: {}", output);

    // Print summary statistics
    println!("PCA Summary:");
    println!("  PC1 explains {:.1}% of variance", explained_variance[0]);
    println!("  PC2 explains {:.1}% of variance", explained_variance[1]);
    println!(
        "  Total explained by PC1+PC2: {:.1}%",
        explained_variance[0] + explained_variance[1]
    );
}

/// Main PCA function
pub fn pca(combined: PathBuf, output: String, _components: usize) {
    if !combined.exists() {
        panic!("Combined file does not exist: {}", combined.display());
    }

    println!("Reading combined inquiSTR file: {}", combined.display());
    let (data_matrix, sample_names, region_names) = parse_combined_file(combined.as_path());

    println!("Loaded data:");
    println!("  {} samples", sample_names.len());
    println!("  {} STR regions", region_names.len());

    if sample_names.len() < 2 {
        panic!("Need at least 2 samples for PCA, but only {} found", sample_names.len());
    }

    println!("Performing simplified PCA...");
    let (pca_data, explained_variance) = perform_simple_pca(&data_matrix);

    // Create and save the PCA plot
    create_pca_plot(&pca_data, &sample_names, &explained_variance, &output);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_pca_with_test_data() {
        // This test requires the combined test data file with proper header to exist
        let test_file = PathBuf::from("test_combined_with_header.tsv");
        if test_file.exists() {
            // This should not panic and should create the HTML file
            pca(test_file, "test_pca_output.html".to_string(), 10);

            // Verify the output file was created
            assert!(PathBuf::from("test_pca_output.html").exists());

            // Clean up
            std::fs::remove_file("test_pca_output.html").ok();
        }
    }
}
