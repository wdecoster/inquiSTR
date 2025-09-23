#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(argparser)))
suppressWarnings(suppressMessages(library(parallel)))

# Chunk-based parallel version of STR_regression.R that processes variants in parallel chunks
# for improved memory efficiency and faster processing
# for improved memory efficiency

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
    
    arg <- argparser::parse_args(p)
    
    # Validation
    if (is.na(arg$input) || is.na(arg$phenocovar) || is.na(arg$phenotype) || 
        is.na(arg$out) || is.na(arg$STRmode) || is.na(arg$outcometype)) {
        stop("Error: Missing required arguments")
    }
    
    if (arg$outcometype == "binary" && is.na(arg$binaryOrder)) {
        stop("Error: --binaryOrder required for binary outcomes")
    }
    
    if (!arg$STRmode %in% c("MEAN", "MAX", "MIN")) {
        stop("Error: STRmode must be MEAN, MAX, or MIN")
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
    
    return(arg)
}

# Prepare phenotype data for analysis
prepare_phenotype_data <- function(phenocovar_file, phenotype_col, covariates_str, binaryOrder = NULL, quiet = FALSE) {
    pheno_data <- fread(phenocovar_file)
    
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
    }
    
    # Remove samples with missing phenotype
    pheno_data <- pheno_data[!is.na(get(phenotype_col))]
    
    return(list(
        data = pheno_data,
        sample_id_col = sample_id_col,
        phenotype_col = phenotype_col,
        covariates = covariates,
        binary_levels = if(!is.null(binaryOrder)) trimws(unlist(strsplit(binaryOrder, ","))) else NULL
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
    
    # Check missingness for this variant
    non_missing_rate <- sum(!is.na(str_values)) / length(str_values)
    if(non_missing_rate < missing_cutoff) {
        return(NULL)  # Skip this variant due to high missingness
    }
    
    # Check minimal length filter if specified
    if(!is.na(minimal_length)) {
        max_length <- max(str_values, na.rm = TRUE)
        if(is.na(max_length) || max_length <= minimal_length) {
            return(NULL)  # Skip this variant due to insufficient length
        }
    }
    
    # Create data frame for this variant
    variant_df <- data.frame(
        sample_id = names(str_values),
        variant_value = as.numeric(str_values),
        stringsAsFactors = FALSE
    )
    names(variant_df)[1] <- pheno_info$sample_id_col
    
    # Join with phenotype data
    analysis_data <- merge(variant_df, pheno_info$data, by = pheno_info$sample_id_col)
    
    # Remove samples with missing variant data
    analysis_data <- analysis_data[!is.na(analysis_data$variant_value), ]
    
    if(nrow(analysis_data) < 10) {
        return(NULL)  # Skip if too few samples
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
                    group_name <- paste0("Group", LETTERS[i])  # GroupA, GroupB, etc.
                    group_stats[[paste0(group_name, "_N")]] <- nrow(group_data)
                    group_stats[[paste0(group_name, "_AvgSize")]] <- round(mean(group_data$variant_value, na.rm = TRUE), 3)
                }
            }
            
            result <- data.frame(
                VariantID = paste(variant_info$chrom, variant_info$start, variant_info$end, sep = "_"),
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
            
            result <- data.frame(
                VariantID = paste(variant_info$chrom, variant_info$start, variant_info$end, sep = "_"),
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
process_variant_chunk <- function(variant_lines, sample_names, str_mode, pheno_info, missing_cutoff, minimal_length, outcometype, quiet = FALSE) {
    results <- list()
    
    for(i in 1:length(variant_lines)) {
        line_data <- strsplit(variant_lines[i], "\t")[[1]]
        
        # Parse variant data
        variant_data <- parse_variant_line(line_data, sample_names, str_mode)
        
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
    
    # Read header to get sample names
    header_line <- readLines(con, n = 1)
    header_cols <- strsplit(header_line, "\t")[[1]]
    
    # Extract sample names from header (every other column starting from 4th, removing _H1 suffix)
    h1_cols <- header_cols[seq(4, length(header_cols), by = 2)]
    sample_names <- gsub("_H1$", "", h1_cols)
    
    if(!arg$quiet) {
        cat("Found", length(sample_names), "samples in input file\n")
        cat("Starting variant processing...\n")
    }
    
    # Write output header
    if(arg$outcometype == "binary") {
        if(!is.null(pheno_info$binary_levels)) {
            output_header <- "VariantID\tOR\tOR_L95\tOR_U95\tOR_stdErr\tPvalue\tN\tAvgSize\tGroupA_N\tGroupB_N\tGroupA_AvgSize\tGroupB_AvgSize"
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
                        arg$quiet
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
                    arg$quiet
                )
            }
            
            # Write results
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
                arg$quiet
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
}

# Main execution
main <- function() {
    arg <- parse_arguments()
    
    if(!arg$quiet) {
        cat("inquiSTR - STR_regression Parallel Version 2.1, Dec 2025\n")
    }
    
    run_streaming_analysis(arg)
}

# Run main function
main()