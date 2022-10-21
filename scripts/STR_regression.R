#!~/miniconda3/bin/Rscript
## Run association testing genome-wide for STRs

# For now in R, to be converted into Rust

# TODOS TO PREVENT USER STUPIDITY
# TODO: check if mode is only SUM MIN or MAX, if not throw early error
# TODO: check that if users supply a begin they also supply an end argument
# TODO: either support running the whole file, or check that users at least supply one of chr, chr-begin-end or bed
# TODO: add a progress bar :-D

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(argparser, quietly = TRUE))

assoc_binary <- function(arg, calls_file, phenotype, no_cols, covariates) {
    binaryOrder <- gsub(",", " ", arg$binaryOrder)
    binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))
    calls_file_selected <- as.data.table(calls_file[calls_file[[phenotype]] == c(binaryOrder_prepared), ])
    calls_file_selected[[phenotype]] <- factor(calls_file_selected[[phenotype]], c(binaryOrder_prepared))
    results_calls_file_selected <- as.data.frame(matrix(0, 1, 11))
    colnames(results_calls_file_selected) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")
    for (i in seq(no_cols, ncol(calls_file_selected), 1)) {
        VariantToBeTested <- as.character(colnames(calls_file_selected)[i])
        selectedtable <- as.data.table(cbind(as.character(calls_file_selected[[phenotype]]), as.numeric(calls_file_selected[[VariantToBeTested]])))
        colnames(selectedtable) <- c(phenotype, VariantToBeTested)
        group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])
        group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])
        AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        if (!is.na(covariates)) {
            covlist <- gsub(",", " ", covariates)
            covlist_prepared <- unlist(strsplit(covlist, split = " "))
            formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")
        } else {
            formulax <- paste(phenotype, VariantToBeTested, sep = "~")
        }
        glm_result <- glm(formula = formulax, data = calls_file_selected, family = binomial(link = "logit"))
        Predictors <- names(glm_result$coefficients)
        VariantID <- names(glm_result$coefficients)[2]
        OR <- round(as.numeric(exp(glm_result$coefficients)), digits = 3)
        OR_L95 <- round(as.numeric(exp(confint.default(glm_result, level = 0.95)[, 1])), digits = 3)
        OR_U95 <- round(as.numeric(exp(confint.default(glm_result, level = 0.95)[, 2])), digits = 3)
        OR_stdErr <- round(as.numeric(coef(summary(glm_result))[, "Std. Error"]), digits = 3)
        Pvalue <- as.numeric(coef(summary(glm_result))[, "Pr(>|z|)"])
        N <- nobs(glm_result)
        model <- as.character(glm_result$formula)[1]
        tabular_result <- as.data.frame(cbind(Predictors, OR, OR_L95, OR_U95, OR_stdErr, Pvalue, N, AvgSize, Group1_AvgSize, Group2_AvgSize, model))
        tabular_result <- subset(tabular_result, Predictors == VariantID)
        colnames(tabular_result)[1] <- "VariantID"
        results_calls_file_selected <- rbind.data.frame(results_calls_file_selected, tabular_result)
    }
    results_calls_file_selected <- results_calls_file_selected[-1, ]
    sorted_results_calls_file_selected <- results_calls_file_selected[order(results_calls_file_selected$Pvalue), ]
    write.table(sorted_results_calls_file_selected, arg$out, sep = "\t", col.names = TRUE, quote = F, row.names = F)
}

assoc_continuous <- function(arg, calls_file, phenotype, no_cols, covariates) {
    results_calls_file <- as.data.frame(matrix(0, 1, 11))
    colnames(results_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")
    for (i in seq(no_cols, ncol(calls_file), 1)) {
        VariantToBeTested <- as.character(colnames(calls_file)[i])
        selectedtable <- as.data.table(cbind(as.character(calls_file[[phenotype]]), as.numeric(calls_file[[VariantToBeTested]])))
        colnames(selectedtable) <- c(phenotype, VariantToBeTested)
        AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        MaxSize <- round(max(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        MinSize <- round(min(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        if (!is.na(covariates)) {
            covlist <- gsub(",", " ", covariates)
            covlist_prepared <- unlist(strsplit(covlist, split = " "))
            formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")
        } else {
            formulax <- paste(phenotype, VariantToBeTested, sep = "~")
        }
        glm_result <- glm(formula = formulax, data = calls_file, family = gaussian(link = "identity"))
        Predictors <- names(glm_result$coefficients)
        VariantID <- names(glm_result$coefficients)[2]
        Beta <- round(as.numeric(exp(glm_result$coefficients)), digits = 3)
        Beta_L95 <- round(as.numeric(exp(confint.default(glm_result, level = 0.95)[, 1])), digits = 3)
        Beta_U95 <- round(as.numeric(exp(confint.default(glm_result, level = 0.95)[, 2])), digits = 3)
        Beta_stdErr <- round(as.numeric(coef(summary(glm_result))[, "Std. Error"]), digits = 3)
        Pvalue <- as.numeric(coef(summary(glm_result))[, "Pr(>|t|)"])
        N <- nobs(glm_result)
        model <- as.character(glm_result$formula)[1]
        tabular_result <- as.data.frame(cbind(Predictors, Beta, Beta_L95, Beta_U95, Beta_stdErr, Pvalue, N, AvgSize, MinSize, MaxSize, model))
        tabular_result <- subset(tabular_result, Predictors == VariantID)
        colnames(tabular_result)[1] <- "VariantID"
        results_calls_file <- rbind.data.frame(results_calls_file, tabular_result)
    }
    results_calls_file <- results_calls_file[-1, ]
    sorted_results_calls_file <- results_calls_file[order(results_calls_file$Pvalue), ]
    write.table(sorted_results_calls_file, arg$out, sep = "\t", col.names = TRUE, quote = F, row.names = F)
}

read_calls <- function(input, chrom) {
    calls_file <- fread(input, header = TRUE)
    calls_file <- subset(calls_file, chrom == chrom)
    strnames <- paste(calls_file$chrom, calls_file$begin, calls_file$end, sep = "_")
    rest <- calls_file[, -c(1:3)]
    col_index <- seq(1:ncol(rest))
    inqH1 <- as.data.table(rest %>% select(col_index[col_index %% 2 != 0]))
    inqH2 <- as.data.table(rest %>% select(col_index[col_index %% 2 == 0]))
    colnames(inqH1) <- gsub(".[^.]+$", "", colnames(inqH1))
    colnames(inqH2) <- gsub(".[^.]+$", "", colnames(inqH2))
    return(list("H1" = inqH1, "H2" = inqH2, "strnames" = strnames))
    rm(calls_file)
}

read_calls_chr_begin_end <- function(input, chr, chr_begin, chr_end) {
    calls_file <- fread(input, header = TRUE)
    calls_file <- subset(calls_file, ((chr == chr) & (chr_begin >= chr_begin) & (chr_end <= chr_end)))
    strnames <- paste(calls_file$chrom, calls_file$begin, calls_file$end, sep = "_")
    rest <- calls_file[, -c(1:3)]
    col_index <- seq(1:ncol(rest))
    inqH1 <- as.data.table(rest %>% select(col_index[col_index %% 2 != 0]))
    inqH2 <- as.data.table(rest %>% select(col_index[col_index %% 2 == 0]))
    colnames(inqH1) <- gsub(".[^.]+$", "", colnames(inqH1))
    colnames(inqH2) <- gsub(".[^.]+$", "", colnames(inqH2))
    return(list("H1" = inqH1, "H2" = inqH2, "strnames" = strnames))
    rm(calls_file)
}

read_calls_bed <- function(input, bed) {
    calls_file <- fread(input, header = TRUE)
    bedfile <- fread(bed, header = FALSE)
    colnames(bedfile) <- c("chrom", "start", "end")
    colnames(calls_file)[2] <- "start"
    intersecttable <- as.data.table(bed_intersect(calls_file, bedfile, suffix = c("", ".y")))
    intersecttable <- intersecttable[, 1:(length(intersecttable) - 3)]
    intersecttable <- subset(intersecttable, !is.na(chrom)) ## WDC here is something wrong - there is no chrom
    colnames(intersecttable)[2] <- "begin"
    intersect_strnames <- paste(intersecttable$chrom, intersecttable$begin, intersecttable$end, sep = "_")
    rest <- intersecttable[, -c(1:3)]
    col_index <- seq(1:ncol(rest))
    inqH1 <- as.data.table(rest %>% select(col_index[col_index %% 2 != 0]))
    inqH2 <- as.data.table(rest %>% select(col_index[col_index %% 2 == 0]))
    colnames(inqH1) <- gsub(".[^.]+$", "", colnames(inqH1))
    colnames(inqH2) <- gsub(".[^.]+$", "", colnames(inqH2))
    return(list("H1" = inqH1, "H2" = inqH2, "strnames" = intersect_strnames))
    rm(calls_file)
    rm(intersecttable)
    rm(bedfile)
}

prepare_phenotype <- function(phenofile, phenotype, sample_list) {
    phenocovar <- fread(phenofile, header = TRUE)
    colnames(sample_list) <- "sample_id"
    return(list(
        "sample_list" = left_join(sample_list, phenocovar, by = "sample_id"),
        "phenotype" = phenotype <- paste0(phenotype, ""),
        "no_cols" = ncol(phenocovar) + 1
    ))
}

prepare_calls <- function(calls, sample_list_wPheno, mode, missing_cutoff) {
    if (mode == "SUM") {
        calls2 <- calls$H1 + calls$H2
    } else if (mode == "MAX") {
        calls2 <- pmax(calls$H1, calls$H2)
    } else if (mode == "MIN") {
        calls2 <- pmin(calls$H1, calls$H2)
    }
    calls_file <- transpose(calls2)
    colnames(calls_file) <- calls$strnames
    calls_file <- cbind(sample_list_wPheno, calls_file)
    calls_file <- calls_file[, which(unlist(lapply(calls_file, function(x) !all(is.na(x))))), with = F]
    calls_file <- data.table(data.frame(calls_file)[, which(colMeans(!is.na(data.frame(calls_file))) >= missing_cutoff)])
    calls_file <- calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))
    return(calls_file)
}

parse_arguments <- function() {
    p <- arg_parser("Run association testing for STRs with different modes and options. Version 1.3, October 20, 2022")
    p <- add_argument(p, "--input", help = "Input STR file with a header, first 3 columns are chrom, begin, end, and rest are sample ids with inqH1 & inqH2 STR lengths", type = "character", nargs = 1)
    p <- add_argument(p, "--phenocovar", help = "Phenotype and covariate file with header, first column is sample_id", type = "character", nargs = 1)
    p <- add_argument(p, "--covnames", help = "Covariate names you want to use (optional), separated by comma", type = "character", nargs = "*")
    p <- add_argument(p, "--phenotype", help = "Column name of your phenotype of interest variable in the --phenocovar file", type = "character", nargs = 1)
    p <- add_argument(p, "--out", help = "Output file name", type = "character", nargs = 1)
    p <- add_argument(p, "--mode", help = "Select a mode name from following: SUM, MAX, MIN", type = "character", nargs = 1)
    p <- add_argument(p, "--missing_cutoff", help = "Defines the call rate cutoff for variants, default is 0.80 meaning that keeping all variants present in at least 80% of individuals", type = "numeric", default = "0.80")
    p <- add_argument(p, "--outcometype", help = "Select a outcome variable type: binary or continuous", type = "character", nargs = 1)
    p <- add_argument(p, "--binaryOrder", help = "Give the binary phenotype order, comma separated, e.g. Control, Patient will code Control as 0/Group1 and Patient as 1/Group2. This will also be used to further filter your data to only two groups (if you had more than >2 groups for your categorical data)", type = "character", nargs = "*")
    p <- add_argument(p, "--chr", help = "Indicate chromosome number to be analyzed (with chr prefix). Optional if bed file is provided.", type = "character", nargs = "*")
    p <- add_argument(p, "--chr_begin", help = "Define a begin position (inclusive) for a region of interest (optional, and should be combined with --chr_end)", type = "integer", nargs = "*")
    p <- add_argument(p, "--chr_end", help = "Define a end position (inclusive) for a region of interest (optional, and should be combined with --chr_begin)", type = "integer", nargs = "*")
    p <- add_argument(p, "--bed", help = "A bed file (without a header) with three columns: chromosome (with chr prefix), begin, and end positions for region(s) of interest (optional). valr Rpackage is required - can be installed with mamba install r-valr on conda environment", type = "character", nargs = "*")
    Version <- "Run association testing for STRs with different modes and options. Version 1.3, October 20, 2022"
    print(Version)

    return(parse_args(p))
}

arg <- parse_arguments()

# calls is a list with H1, H2 and strnames attributes
if ((!is.na(arg$chr) && !is.na(arg$chr_begin) && !is.na(arg$chr_end))) {
    calls <- read_calls_chr_begin_end(input = arg$input, chr = arg$chr, chr_begin = arg$chr_begin, chr_end = arg$chr_end)
} else if (!is.na(arg$bed)) {
    suppressPackageStartupMessages(library(valr))
    calls <- read_calls_bed(input = arg$input, bed = arg$bed)
} else if (!is.na(arg$chr)) {
    calls <- read_calls(input = arg$input, chrom = arg$chr)
}

# pheno_info is a list with sample_list, phenotype and no_cols attributes
pheno_info <- prepare_phenotype(
    phenofile = arg$phenocovar,
    phenotype = arg$phenotype,
    sample_list = as.data.table(colnames(calls$H1))
)

calls_file <- prepare_calls(
    calls = calls,
    sample_list_wPheno = pheno_info$sample_list,
    mode = arg$mode,
    missing_cutoff = arg$missing_cutoff
)

if (arg$outcometype == "binary") {
    assoc_binary(
        arg = arg,
        calls_file = calls_file,
        phenotype = pheno_info$phenotype,
        no_cols = pheno_info$no_cols,
        covariates = arg$covnames
    )
} else if (arg$outcometype == "continuous") {
    assoc_continuous(
        arg = arg,
        calls_file = calls_file,
        phenotype = pheno_info$phenotype,
        no_cols = pheno_info$no_cols,
        covariates = arg$covnames
    )
}
