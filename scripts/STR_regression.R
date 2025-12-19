#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(argparser)))
suppressWarnings(suppressMessages(library(parallel)))

# Chunk-based parallel version of STR_regression.R that processes variants in parallel chunks
# for improved memory efficiency and faster processing

parse_arguments <- function() {
    p <- argparser::arg_parser("Run association testing for STRs with parallel chunk-based approach. Version 2.1, Dec 2025")
    p <- argparser::add_argument(p, "--input", help = "inquiSTR input STR file with header", type = "character", nargs = 1)
    p <- argparser::add_argument(p, "--phenocovar", help = "Phenotype and covariate file with header, first column is individual", type = "character", nargs = 1)
    p <- argparser::add_argument(p, "--covnames", help = "Covariate names, comma separated (optional)", type = "character", nargs = "*")
    p <- argparser::add_argument(p, "--phenotype", help = "Column name of phenotype in --phenocovar file", type = "character", nargs = 1)
    p <- argparser::add_argument(p, "--out", help = "Output file name", type = "character", nargs = 1)
    p <- argparser::add_argument(p, "--STRmode", help = "MEAN, MAX, or MIN for H1/H2 combination", type = "character", nargs = 1)
    p <- argparser::add_argument(p, "--missing_cutoff", help = "Call rate cutoff for variants (default 0.80)", type = "numeric", default = "0.80")
    p <- argparser::add_argument(p, "--minimal_length", help = "Minimum maximum STR length across samples for variant to be included", type = "numeric", nargs = "?")
    p <- argparser::add_argument(p, "--threads", help = "Number of threads for parallel processing (default 1)", type = "integer", default = 1)
    p <- argparser::add_argument(p, "--chunk_size", help = "Number of variants to process in each chunk (default 1000)", type = "integer", default = 1000)
    p <- argparser::add_argument(p, "--outcometype", help = "binary or continuous", type = "character", nargs = 1)
    p <- argparser::add_argument(p, "--binaryOrder", help = "Binary phenotype order, comma separated (e.g., Control,Patient)", type = "character", nargs = "*")
    p <- argparser::add_argument(p, "--quiet", help = "Do not print progress messages", flag = TRUE)
    p <- argparser::add_argument(p, "--plot", help = "Prefix for QQ plot and Manhattan plot filenames", type = "character", nargs = "?")
    p <- argparser::add_argument(p, "--sort", help = "Sort results by Bonferroni corrected p-value (requires loading full results into memory)", flag = TRUE)
    
    arg <- argparser::parse_args(p)
    
    # Detect input file type (STR vs kmer)
    arg$input_type <- detect_input_type(arg$input)
    
    # Validation
    if (is.na(arg$input) || is.na(arg$phenocovar) || is.na(arg$phenotype) || 
        is.na(arg$out) || is.na(arg$outcometype)) {
        stop("Error: Missing required arguments")
    }
    
    # STRmode is only required for STR data
    if (arg$input_type == "str" && is.na(arg$STRmode)) {
        stop("Error: --STRmode required for STR data")
    }
    
    if (arg$outcometype == "binary" && is.na(arg$binaryOrder)) {
        stop("Error: --binaryOrder required for binary outcomes")
    }
    
    if (arg$input_type == "str" && !arg$STRmode %in% c("MEAN", "MAX", "MIN")) {
        stop("Error: STRmode must be MEAN, MAX, or MIN for STR data")
    }
    
    if (!arg$outcometype %in% c("binary", "continuous")) {
        stop("Error: outcometype must be binary or continuous")
    }
    
    # Validate file existence
    if (!file.exists(arg$input)) {
        stop("Error: Input file does not exist: ", arg$input)
    }
    
    if (!file.exists(arg$phenocovar)) {
        stop("Error: Phenotype file does not exist: ", arg$phenocovar)
    }
    
    # Validate numeric arguments
    if (arg$missing_cutoff < 0 || arg$missing_cutoff > 1) {
        stop("Error: missing_cutoff must be between 0 and 1")
    }
    
    if (!is.na(arg$minimal_length) && arg$minimal_length < 0) {
        stop("Error: minimal_length must be non-negative")
    }
    
    if (arg$threads < 1) {
        stop("Error: threads must be >= 1")
    }
    
    if (arg$chunk_size < 1) {
        stop("Error: chunk_size must be >= 1")
    }
    
    # Warn if plot flag is used with kmer data (no chromosome/position info)
    if (arg$input_type == "kmer" && !is.na(arg$plot)) {
        if(!arg$quiet) {
            cat("Warning: --plot flag not supported for kmer frequency data (no chromosome/position information)\n", file = stderr())
            cat("Plots will be skipped. Continuing with association analysis...\n", file = stderr())
        }
        arg$plot <- NA  # Disable plotting
    }
    
    return(arg)
}

# Detect input file type (STR vs kmer)
detect_input_type <- function(input_file) {
    # Read first non-metadata line to determine file type
    con <- file(input_file, "r")
    
    while(TRUE) {
        line <- readLines(con, n = 1, warn = FALSE)
        if(length(line) == 0) {
            close(con)
            stop("Error: Empty input file")
        }
        
        # Skip metadata lines
        if(startsWith(line, "#")) {
            # Check file_type metadata if present
            if(grepl("^# file_type=", line)) {
                if(grepl("kmer", line)) {
                    close(con)
                    return("kmer")
                } else if(grepl("call", line)) {
                    close(con)
                    return("str")
                }
            }
            next
        }
        
        # This is the header line - check first column
        fields <- strsplit(line, "\t")[[1]]
        close(con)
        
        if(fields[1] == "kmer") {
            return("kmer")
        } else if(fields[1] == "chromosome") {
            return("str")
        } else {
            stop("Error: Unrecognized file format. Expected 'kmer' or 'chromosome' as first column, got: ", fields[1])
        }
    }
}

# Prepare phenotype data for analysis
prepare_phenotype_data <- function(phenocovar_file, phenotype_col, covariates_str, binaryOrder = NULL, quiet = FALSE) {
    # Read with explicit header=TRUE to ensure first row is used as column names
    pheno_data <- fread(phenocovar_file, header = TRUE)
    
    # Clean column names (remove any hidden characters)
    names(pheno_data) <- trimws(gsub("[\r\n\t]", "", names(pheno_data)))
    
    # Parse covariates
    covariates <- if(!is.na(covariates_str) && covariates_str != "") {
        trimws(unlist(strsplit(covariates_str, ",")))
    } else {
        character(0)
    }
    
    # Check required columns exist
    required_cols <- c(names(pheno_data)[1], phenotype_col, covariates)  # First column is sample ID
    missing_cols <- setdiff(required_cols, names(pheno_data))
    if(length(missing_cols) > 0) {
        stop("Missing columns in phenotype file: ", paste(missing_cols, collapse = ", "))
    }
    
    # Set sample ID column name
    sample_id_col <- names(pheno_data)[1]
    
    # Filter for binary outcomes if needed
    if(!is.null(binaryOrder)) {
        binary_levels <- trimws(unlist(strsplit(binaryOrder, ",")))
        pheno_data <- pheno_data[get(phenotype_col) %in% binary_levels]
        pheno_data[[phenotype_col]] <- factor(pheno_data[[phenotype_col]], levels = binary_levels)
        if(!quiet) {
            cat("Filtered to", nrow(pheno_data), "samples with binary phenotype levels:", 
                paste(binary_levels, collapse = ", "), "\n")
        }
        
        # Validate that we have samples in each binary group
        if(nrow(pheno_data) == 0) {
            stop("Error: No samples found with the specified binary phenotype levels: ", 
                 paste(binary_levels, collapse = ", "), 
                 ". Check your --binaryOrder specification and phenotype file.")
        }
        
        # Check that we have at least one sample in each binary group
        level_counts <- table(pheno_data[[phenotype_col]])
        missing_levels <- binary_levels[!binary_levels %in% names(level_counts)]
        if(length(missing_levels) > 0) {
            stop("Error: No samples found for binary phenotype level(s): ", 
                 paste(missing_levels, collapse = ", "), 
                 ". Available levels in data: ", 
                 paste(unique(pheno_data[[phenotype_col]]), collapse = ", "))
        }
        
        # Check minimum sample count per group
        min_samples_per_group <- 2  # Minimum for statistical analysis
        low_count_levels <- names(level_counts)[level_counts < min_samples_per_group]
        if(length(low_count_levels) > 0) {
            stop("Error: Insufficient samples in binary phenotype level(s): ", 
                 paste(paste(low_count_levels, " (n=", level_counts[low_count_levels], ")", sep=""), collapse = ", "),
                 ". Need at least ", min_samples_per_group, " samples per group for analysis.")
        }
        
        if(!quiet) {
            cat("Sample counts per group: ", paste(paste(names(level_counts), "=", level_counts), collapse = ", "), "\n")
        }
    }
    
    # Remove samples with missing phenotype
    pheno_data <- pheno_data[!is.na(get(phenotype_col))]
    
    # Final validation: ensure we have enough samples for analysis
    min_total_samples <- 10  # Minimum for meaningful statistical analysis
    if(nrow(pheno_data) < min_total_samples) {
        stop("Error: Insufficient samples for analysis. Found ", nrow(pheno_data), 
             " samples with valid phenotype data, but need at least ", min_total_samples, 
             " for meaningful statistical analysis.")
    }
    
    return(list(
        data = pheno_data,
        sample_id_col = sample_id_col,
        phenotype_col = phenotype_col,
        covariates = covariates,
        binary_levels = if(!is.null(binaryOrder)) trimws(unlist(strsplit(binaryOrder, ","))) else NULL
    ))
}

# Parse a single kmer line and extract frequency data
parse_kmer_line <- function(line_data, sample_names) {
    # Extract kmer info (first column)
    kmer_seq <- line_data[1]
    
    # Calculate number of samples
    n_samples <- length(sample_names)
    
    # Extract frequency data (columns 2 onwards)
    if(length(line_data) < (1 + n_samples)) {
        stop("Error: Kmer line has fewer columns than expected. Expected ", 1 + n_samples, ", got ", length(line_data))
    }
    
    frequency_data <- as.numeric(line_data[2:(1 + n_samples)])
    names(frequency_data) <- sample_names
    
    return(list(
        kmer_info = list(kmer = kmer_seq),
        kmer_values = frequency_data
    ))
}

# Parse a single variant line and extract H1/H2 data
parse_variant_line <- function(line_data, sample_names, str_mode) {
    # Extract variant info (first 3 columns: chr, start, end)
    variant_info <- list(
        chrom = line_data[1],
        start = as.numeric(line_data[2]),
        end = as.numeric(line_data[3])
    )
    
    # Calculate number of samples
    n_samples <- length(sample_names)
    
    # Extract H1 and H2 data (columns 4 onwards, alternating H1/H2)
    h1_indices <- seq(4, 3 + 2*n_samples, by = 2)
    h2_indices <- seq(5, 4 + 2*n_samples, by = 2)
    
    h1_data <- as.numeric(line_data[h1_indices])
    h2_data <- as.numeric(line_data[h2_indices])
    
    names(h1_data) <- sample_names
    names(h2_data) <- sample_names
    
    # Apply STR mode
    if(str_mode == "MEAN") {
        str_values <- (h1_data + h2_data) / 2
    } else if(str_mode == "MAX") {
        str_values <- pmax(h1_data, h2_data, na.rm = FALSE)
    } else if(str_mode == "MIN") {
        str_values <- pmin(h1_data, h2_data, na.rm = FALSE)
    }
    
    return(list(
        variant_info = variant_info,
        str_values = str_values,
        h1_data = h1_data,
        h2_data = h2_data
    ))
}

# Process a single variant through GLM analysis
process_single_variant <- function(variant_data, pheno_info, missing_cutoff, minimal_length, outcometype, quiet = FALSE) {
    str_values <- variant_data$str_values
    variant_info <- variant_data$variant_info
    
    # Create data frame for this variant
    variant_df <- data.frame(
        sample_id = names(str_values),
        variant_value = as.numeric(str_values),
        stringsAsFactors = FALSE
    )
    names(variant_df)[1] <- pheno_info$sample_id_col
    
    # Join with phenotype data
    analysis_data <- merge(variant_df, pheno_info$data, by = pheno_info$sample_id_col)
    
    # Check missingness for this variant AFTER merging with phenotype data
    # This ensures we only calculate call rate for samples that have phenotype data
    non_missing_count <- sum(!is.na(analysis_data$variant_value))
    non_missing_rate <- non_missing_count / nrow(analysis_data)
    if(non_missing_rate < missing_cutoff) {
        return(NULL)  # Skip this variant due to high missingness
    }
    
    # Remove samples with missing variant data
    analysis_data <- analysis_data[!is.na(analysis_data$variant_value), ]
    
    if(nrow(analysis_data) < 10) {
        return(NULL)  # Skip if too few samples
    }
    
    # Check minimal length filter if specified (after filtering to phenotyped samples)
    if(!is.na(minimal_length)) {
        max_length <- max(analysis_data$variant_value, na.rm = TRUE)
        if(is.na(max_length) || max_length <= minimal_length) {
            return(NULL)  # Skip this variant due to insufficient length
        }
    }
    
    # Build formula
    formula_str <- paste(pheno_info$phenotype_col, "~ variant_value")
    if(length(pheno_info$covariates) > 0) {
        formula_str <- paste(formula_str, "+", paste(pheno_info$covariates, collapse = " + "))
    }
    
    # Run GLM
    tryCatch({
        if(outcometype == "binary") {
            model <- glm(as.formula(formula_str), data = analysis_data, family = binomial(link = "logit"))
            coef_summary <- summary(model)$coefficients
            variant_coef <- coef_summary["variant_value", ]
            
            # Calculate OR and confidence intervals
            or <- exp(variant_coef["Estimate"])
            ci <- exp(confint.default(model)["variant_value", ])
            
            # Calculate group statistics if binary
            group_stats <- NULL
            if(!is.null(pheno_info$binary_levels)) {
                for(i in 1:length(pheno_info$binary_levels)) {
                    level <- pheno_info$binary_levels[i]
                    group_data <- analysis_data[analysis_data[[pheno_info$phenotype_col]] == level, ]
                    # Use actual phenotype level name instead of generic GroupA/GroupB
                    group_name <- level
                    group_stats[[paste0(group_name, "_N")]] <- nrow(group_data)
                    group_stats[[paste0(group_name, "_AvgSize")]] <- round(mean(group_data$variant_value, na.rm = TRUE), 3)
                }
            }
            
            # Create VariantID based on data type
            variant_id <- if(!is.null(variant_info$kmer)) {
                variant_info$kmer
            } else {
                paste(variant_info$chrom, variant_info$start, variant_info$end, sep = "_")
            }
            
            result <- data.frame(
                VariantID = variant_id,
                OR = round(or, 3),
                OR_L95 = round(ci[1], 3),
                OR_U95 = round(ci[2], 3),
                OR_stdErr = round(variant_coef["Std. Error"], 3),
                Pvalue = variant_coef["Pr(>|z|)"],
                N = nrow(analysis_data),
                AvgSize = round(mean(analysis_data$variant_value, na.rm = TRUE), 3),
                stringsAsFactors = FALSE
            )
            
            # Add group statistics if available
            if(!is.null(group_stats)) {
                for(name in names(group_stats)) {
                    result[[name]] <- group_stats[[name]]
                }
            }
            
        } else { # continuous
            model <- glm(as.formula(formula_str), data = analysis_data, family = gaussian())
            coef_summary <- summary(model)$coefficients
            variant_coef <- coef_summary["variant_value", ]
            
            # Calculate confidence intervals
            ci <- confint.default(model)["variant_value", ]
            
            # Create VariantID based on data type
            variant_id <- if(!is.null(variant_info$kmer)) {
                variant_info$kmer
            } else {
                paste(variant_info$chrom, variant_info$start, variant_info$end, sep = "_")
            }
            
            result <- data.frame(
                VariantID = variant_id,
                Beta = round(variant_coef["Estimate"], 3),
                Beta_L95 = round(ci[1], 3),
                Beta_U95 = round(ci[2], 3),
                Beta_stdErr = round(variant_coef["Std. Error"], 3),
                Pvalue = variant_coef["Pr(>|t|)"],
                N = nrow(analysis_data),
                AvgSize = round(mean(analysis_data$variant_value, na.rm = TRUE), 3),
                MinSize = round(min(analysis_data$variant_value, na.rm = TRUE), 3),
                MaxSize = round(max(analysis_data$variant_value, na.rm = TRUE), 3),
                stringsAsFactors = FALSE
            )
        }
        
        return(result)
        
    }, error = function(e) {
        if(!quiet) {
            cat("Error processing variant", variant_info$chrom, variant_info$start, ":", e$message, "\n", file = stderr())
        }
        return(NULL)
    })
}

# Process a chunk of variants in parallel
process_variant_chunk <- function(variant_lines, sample_names, str_mode, pheno_info, missing_cutoff, minimal_length, outcometype, quiet = FALSE, input_type = "str") {
    results <- list()
    
    for(i in 1:length(variant_lines)) {
        line_data <- strsplit(variant_lines[i], "\t")[[1]]
        
        # Parse variant data based on input type
        if(input_type == "kmer") {
            variant_data <- parse_kmer_line(line_data, sample_names)
            # Use kmer values directly (no H1/H2 combination needed)
            variant_data$str_values <- variant_data$kmer_values
            variant_data$variant_info <- variant_data$kmer_info
        } else {
            variant_data <- parse_variant_line(line_data, sample_names, str_mode)
        }
        
        # Process the variant
        result <- process_single_variant(
            variant_data, 
            pheno_info, 
            missing_cutoff,
            minimal_length,
            outcometype,
            quiet
        )
        
        if(!is.null(result)) {
            results[[length(results) + 1]] <- result
        }
    }
    
    return(results)
}

# Main streaming analysis function
run_streaming_analysis <- function(arg) {
    if(!arg$quiet) {
        cat("Starting streaming STR association analysis\n")
        cat("STR mode:", arg$STRmode, "\n")
        cat("Outcome type:", arg$outcometype, "\n")
        cat("Missing cutoff:", arg$missing_cutoff, "\n")
        if(!is.na(arg$minimal_length)) {
            cat("Minimal length filter:", arg$minimal_length, "\n")
        }
    }
    
    # Prepare phenotype data
    pheno_info <- prepare_phenotype_data(
        arg$phenocovar, 
        arg$phenotype, 
        arg$covnames,
        if(arg$outcometype == "binary") arg$binaryOrder else NULL,
        arg$quiet
    )
    
    if(!arg$quiet) {
        cat("Loaded", nrow(pheno_info$data), "samples with phenotype data\n")
    }
    
    # Open input file for streaming
    con <- file(arg$input, "r")
    
    # Skip metadata lines (lines starting with #) and read header
    header_line <- readLines(con, n = 1)
    while(length(header_line) > 0 && grepl("^#", header_line)) {
        header_line <- readLines(con, n = 1)
    }
    
    if(length(header_line) == 0) {
        close(con)
        stop("Error: Input file is empty or contains only metadata lines")
    }
    
    header_cols <- strsplit(header_line, "\t")[[1]]
    
    # Extract sample names based on input type
    if(arg$input_type == "kmer") {
        # Kmer format: first column is 'kmer', rest are sample names
        sample_names <- header_cols[-1]  # Remove first column (kmer)
    } else {
        # STR format: extract sample names from header (every other column starting from 5th, removing _H1 suffix)
        h1_cols <- header_cols[seq(5, length(header_cols), by = 2)]
        sample_names <- gsub("_H1$", "", h1_cols)
    }
    
    if(!arg$quiet) {
        cat("Found", length(sample_names), "samples in input file\n")
        if(arg$input_type == "kmer") {
            cat("Starting kmer frequency processing...\n")
        } else {
            cat("Starting variant processing...\n")
        }
    }
    
    # Write output header
    if(arg$outcometype == "binary") {
        if(!is.null(pheno_info$binary_levels)) {
            # Use actual phenotype labels from binary_levels
            # Match the order in which group_stats are added: label1_N, label1_AvgSize, label2_N, label2_AvgSize
            label1 <- pheno_info$binary_levels[1]
            label2 <- pheno_info$binary_levels[2]
            output_header <- paste0("VariantID\tOR\tOR_L95\tOR_U95\tOR_stdErr\tPvalue\tN\tAvgSize\t",
                                   label1, "_N\t", label1, "_AvgSize\t", label2, "_N\t", label2, "_AvgSize")
        } else {
            output_header <- "VariantID\tOR\tOR_L95\tOR_U95\tOR_stdErr\tPvalue\tN\tAvgSize"
        }
    } else {
        output_header <- "VariantID\tBeta\tBeta_L95\tBeta_U95\tBeta_stdErr\tPvalue\tN\tAvgSize\tMinSize\tMaxSize"
    }
    
    writeLines(output_header, arg$out)
    
    # Set up parallel processing
    num_cores <- min(arg$threads, detectCores())
    if(!arg$quiet) {
        cat("Using", num_cores, "cores for parallel processing\n")
    }
    
    # Process variants in chunks
    variants_processed <- 0
    variants_written <- 0
    chunk_lines <- character(0)
    
    if(!arg$quiet) {
        cat("Reading and processing variants in chunks of", arg$chunk_size, "...\n")
    }
    
    while(length(line <- readLines(con, n = 1)) > 0) {
        chunk_lines <- c(chunk_lines, line)
        variants_processed <- variants_processed + 1
        
        # Process chunk when it reaches chunk_size or end of file
        if(length(chunk_lines) >= arg$chunk_size) {
            # Split chunk into sub-chunks for parallel processing
            if(num_cores > 1) {
                # Split chunk_lines into smaller sub-chunks for each core
                sub_chunk_size <- ceiling(length(chunk_lines) / num_cores)
                sub_chunks <- split(chunk_lines, ceiling(seq_along(chunk_lines) / sub_chunk_size))
                
                # Process sub-chunks in parallel
                chunk_results <- mclapply(sub_chunks, function(sub_chunk) {
                    process_variant_chunk(
                        sub_chunk, 
                        sample_names, 
                        arg$STRmode, 
                        pheno_info, 
                        arg$missing_cutoff, 
                        arg$minimal_length, 
                        arg$outcometype, 
                        arg$quiet,
                        arg$input_type
                    )
                }, mc.cores = num_cores)
                
                # Flatten results
                all_results <- unlist(chunk_results, recursive = FALSE)
            } else {
                # Single-threaded processing
                all_results <- process_variant_chunk(
                    chunk_lines, 
                    sample_names, 
                    arg$STRmode, 
                    pheno_info, 
                    arg$missing_cutoff, 
                    arg$minimal_length, 
                    arg$outcometype, 
                    arg$quiet,
                    arg$input_type
                )
            }
            
            # Write final results
            if(length(all_results) > 0) {
                result_lines <- character(length(all_results))
                for(i in 1:length(all_results)) {
                    result_values <- as.character(unlist(all_results[[i]]))
                    result_lines[i] <- paste(result_values, collapse = "\t")
                }
                write(paste(result_lines, collapse = "\n"), file = arg$out, append = TRUE)
                variants_written <- variants_written + length(all_results)
            }
            
            # Progress reporting
            if(!arg$quiet) {
                pass_rate <- if(variants_processed > 0) round(variants_written / variants_processed * 100, 1) else 0
                cat("Processed", variants_processed, "variants,", variants_written, "passed filters")
                cat(" (", pass_rate, "% pass rate)\n")
            }
            
            # Reset chunk
            chunk_lines <- character(0)
        }
    }
    
    # Process any remaining variants in the last incomplete chunk
    if(length(chunk_lines) > 0) {
        if(num_cores > 1 && length(chunk_lines) >= num_cores) {
            # Split remaining chunk for parallel processing
            sub_chunk_size <- ceiling(length(chunk_lines) / num_cores)
            sub_chunks <- split(chunk_lines, ceiling(seq_along(chunk_lines) / sub_chunk_size))
            
            chunk_results <- mclapply(sub_chunks, function(sub_chunk) {
                process_variant_chunk(
                    sub_chunk, 
                    sample_names, 
                    arg$STRmode, 
                    pheno_info, 
                    arg$missing_cutoff, 
                    arg$minimal_length, 
                    arg$outcometype, 
                    arg$quiet
                )
            }, mc.cores = num_cores)
            
            all_results <- unlist(chunk_results, recursive = FALSE)
        } else {
            # Process remaining variants sequentially
            all_results <- process_variant_chunk(
                chunk_lines, 
                sample_names, 
                arg$STRmode, 
                pheno_info, 
                arg$missing_cutoff, 
                arg$minimal_length, 
                arg$outcometype, 
                arg$quiet,
                input_type
            )
        }
        
        # Write final results
        if(length(all_results) > 0) {
            result_lines <- character(length(all_results))
            for(i in 1:length(all_results)) {
                result_values <- as.character(unlist(all_results[[i]]))
                result_lines[i] <- paste(result_values, collapse = "\t")
            }
            write(paste(result_lines, collapse = "\n"), file = arg$out, append = TRUE)
            variants_written <- variants_written + length(all_results)
        }
    }
    
    close(con)
    
    # Calculate summary statistics
    pass_rate <- if(variants_processed > 0) round(variants_written / variants_processed * 100, 1) else 0
    
    if(!arg$quiet) {
        cat("Analysis complete!\n")
        cat("Total variants processed:", variants_processed, "\n")
        cat("Total variants written:", variants_written, "\n")
        cat("Pass rate:", pass_rate, "%\n")
        cat("Results written to:", arg$out, "\n")
        
        # Report filtering statistics
        cat("\nFiltering summary:\n")
        cat("- Variants failing missing cutoff:", variants_processed - variants_written, "\n")
        if(!is.na(arg$minimal_length)) {
            cat("- Minimal length filter:", arg$minimal_length, "\n")
        }
    }
    
    # Add multiple testing correction
    if(variants_written > 0) {
        if(!arg$quiet) {
            cat("\nAdding multiple testing correction...\n")
        }
        add_corrected_pvalues(arg$out, variants_written, arg$quiet, arg$sort)
    }
}

# Add Bonferroni corrected p-values to output file using streaming
# If sort=TRUE, loads full results into memory to sort by p-value
add_corrected_pvalues <- function(output_file, n_tests, quiet = FALSE, sort = FALSE) {
    
    # If sorting is requested, load into memory and sort
    if(sort) {
        if(!quiet) {
            cat("Loading results into memory for sorting...\n")
        }
        
        results <- tryCatch({
            fread(output_file, header = TRUE)
        }, error = function(e) {
            if(!quiet) {
                cat("Warning: Could not read output file: ", e$message, "\n", sep = "", file = stderr())
            }
            return(NULL)
        })
        
        if(is.null(results) || nrow(results) == 0 || !"Pvalue" %in% colnames(results)) {
            if(!quiet) {
                cat("Warning: Cannot sort results\n", file = stderr())
            }
            return(invisible(NULL))
        }
        
        # Add Bonferroni correction
        results$Pvalue_bonf <- pmin(results$Pvalue * n_tests, 1.0)
        
        # Sort by raw p-value
        results <- results[order(results$Pvalue), ]
        
        # Count significant results
        bonf_sig <- sum(results$Pvalue_bonf < 0.05, na.rm = TRUE)
        
        # Write sorted results
        tryCatch({
            fwrite(results, output_file, sep = "\t", quote = FALSE)
            if(!quiet) {
                cat("Added Bonferroni corrected p-values and sorted by significance\n")
                cat("- Total tests:", n_tests, "\n")
                cat("- Bonferroni significant (p < 0.05):", bonf_sig, "\n")
            }
        }, error = function(e) {
            if(!quiet) {
                cat("Warning: Could not write sorted results: ", e$message, "\n", sep = "", file = stderr())
            }
        })
        
        return(invisible(NULL))
    }
    
    # Otherwise use streaming approach
    tryCatch({
        con_in <- file(output_file, "r")
        temp_file <- tempfile(fileext = ".tsv")
        con_out <- file(temp_file, "w")
        
        # Read and modify header
        header <- readLines(con_in, n = 1)
        header_cols <- strsplit(header, "\t")[[1]]
        pvalue_col <- which(header_cols == "Pvalue")
        
        if(length(pvalue_col) == 0) {
            close(con_in)
            close(con_out)
            if(!quiet) {
                cat("Warning: Pvalue column not found, skipping correction\n", file = stderr())
            }
            return(invisible(NULL))
        }
        
        writeLines(paste0(header, "\tPvalue_bonf"), con_out)
        
        # Stream through file and add Bonferroni correction
        bonf_sig <- 0
        
        while(length(line <- readLines(con_in, n = 1)) > 0) {
            fields <- strsplit(line, "\t")[[1]]
            pval <- as.numeric(fields[pvalue_col])
            
            if(!is.na(pval)) {
                # Bonferroni: multiply by number of tests (cap at 1.0)
                pval_bonf <- min(pval * n_tests, 1.0)
                
                # Count significant results
                if(pval_bonf < 0.05) bonf_sig <- bonf_sig + 1
                
                # Write line with corrected p-value
                writeLines(paste0(line, "\t", pval_bonf), con_out)
            } else {
                # Missing p-value, write NA
                writeLines(paste0(line, "\tNA"), con_out)
            }
        }
        
        close(con_in)
        close(con_out)
        
        # Replace original file with corrected version
        file.rename(temp_file, output_file)
        
        if(!quiet) {
            cat("Added Bonferroni corrected p-values\n")
            cat("- Total tests:", n_tests, "\n")
            cat("- Bonferroni significant (p < 0.05):", bonf_sig, "\n")
        }
    }, error = function(e) {
        if(!quiet) {
            cat("Warning: Could not write corrected p-values: ", e$message, "\n", sep = "", file = stderr())
        }
    })
}

# Generate QQ and Manhattan plots
generate_plots <- function(results_file, out_prefix, quiet) {
    # Check if qqman is available, install if needed
    if (!requireNamespace("qqman", quietly = TRUE)) {
        if(!quiet) {
            cat("qqman package not installed. Installing now...\n")
        }
        
        install_result <- tryCatch({
            install.packages("qqman", repos = "https://cran.rstudio.com/", quiet = quiet)
            if(!quiet) {
                cat("qqman installed successfully.\n")
            }
            TRUE
        }, error = function(e) {
            cat("Warning: Failed to install qqman package: ", e$message, "\n", sep = "", file = stderr())
            cat("Skipping plot generation. You can manually install with: install.packages('qqman')\n", file = stderr())
            FALSE
        })
        
        if (!install_result) {
            return(invisible(NULL))
        }
        
        # Verify installation succeeded
        if (!requireNamespace("qqman", quietly = TRUE)) {
            cat("Warning: qqman installation failed. Skipping plot generation.\n", file = stderr())
            return(invisible(NULL))
        }
    }
    
    suppressWarnings(suppressMessages(library(qqman)))
    
    if(!quiet) {
        cat("Generating plots...\n")
    }
    
    # Read results
    results <- tryCatch({
        fread(results_file, header = TRUE)
    }, error = function(e) {
        cat("Warning: Could not read results file for plotting: ", e$message, "\n", sep = "", file = stderr())
        NULL
    })
    
    if(is.null(results) || nrow(results) == 0) {
        cat("Warning: No results to plot\n", file = stderr())
        return(invisible(NULL))
    }
    
    # Check required columns
    if(!"Pvalue" %in% colnames(results)) {
        cat("Warning: Pvalue column not found in results\n", file = stderr())
        return(invisible(NULL))
    }
    
    # Parse VariantID to get CHR and BP
    # VariantID format is: chr_start_end (underscore-separated)
    results[, c("CHR", "BP") := tstrsplit(VariantID, "_", keep = 1:2)]
    results[, CHR := gsub("^chr", "", CHR)]  # Remove 'chr' prefix if present
    results[, BP := as.numeric(BP)]
    
    # Filter out rows with missing CHR or BP
    results <- results[!is.na(CHR) & !is.na(BP) & !is.na(Pvalue) & Pvalue > 0]
    
    if(nrow(results) == 0) {
        cat("Warning: No valid variants for plotting after filtering\n", file = stderr())
        return(invisible(NULL))
    }
    
    # Convert CHR to numeric (handle X, Y, MT)
    results[, CHR_numeric := CHR]
    results[CHR == "X", CHR_numeric := "23"]
    results[CHR == "Y", CHR_numeric := "24"]
    results[CHR == "MT" | CHR == "M", CHR_numeric := "25"]
    results[, CHR_numeric := as.numeric(CHR_numeric)]
    
    # Filter out chromosomes that couldn't be converted
    results <- results[!is.na(CHR_numeric)]
    
    if(nrow(results) == 0) {
        cat("Warning: No variants with valid chromosomes for plotting\n", file = stderr())
        return(invisible(NULL))
    }
    
    # Prepare data for qqman (needs CHR as numeric, BP, P, and SNP columns)
    plot_data <- data.frame(
        SNP = results$VariantID,
        CHR = results$CHR_numeric,
        BP = results$BP,
        P = results$Pvalue
    )
    
    # Generate Manhattan plot
    manhattan_file <- paste0(out_prefix, "_manhattan.png")
    tryCatch({
        png(manhattan_file, width = 1200, height = 600, res = 100)
        manhattan(plot_data, 
                 chr = "CHR", 
                 bp = "BP", 
                 p = "P", 
                 snp = "SNP",
                 main = "Manhattan Plot - STR Association",
                 col = c("blue4", "orange3"),
                 suggestiveline = -log10(1e-5),
                 genomewideline = -log10(5e-8),
                 chrlabs = NULL)
        dev.off()
        if(!quiet) {
            cat("Manhattan plot saved to:", manhattan_file, "\n")
        }
    }, error = function(e) {
        cat("Warning: Failed to generate Manhattan plot: ", e$message, "\n", sep = "", file = stderr())
        if(dev.cur() > 1) dev.off()
    })
    
    # Generate QQ plot
    qq_file <- paste0(out_prefix, "_qq.png")
    tryCatch({
        png(qq_file, width = 600, height = 600, res = 100)
        qq(plot_data$P, 
           main = "QQ Plot - STR Association")
        dev.off()
        if(!quiet) {
            cat("QQ plot saved to:", qq_file, "\n")
        }
    }, error = function(e) {
        cat("Warning: Failed to generate QQ plot: ", e$message, "\n", sep = "", file = stderr())
        if(dev.cur() > 1) dev.off()
    })
}

# Main execution
main <- function() {
    arg <- parse_arguments()
    
    if(!arg$quiet) {
        cat("inquiSTR - STR_regression Parallel Version 2.1, Dec 2025\n")
    }
    
    run_streaming_analysis(arg)
    
    # Generate plots if requested
    if(!is.na(arg$plot)) {
        # Use provided plot prefix, or default to output filename stem
        plot_prefix <- if(arg$plot != "") {
            arg$plot
        } else {
            sub("\\.[^.]*$", "", arg$out)
        }
        generate_plots(arg$out, plot_prefix, arg$quiet)
    }
}

# Run main function
main()