#!~/miniconda3/bin/Rscript
## Written by Fahri Kucukali to run association testing genome-wide for STRs

# For now in R, to be converted into Rust

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(argparser, quietly = TRUE))

p <- arg_parser("Run association testing for STRs with different modes and options, written by Fahri Kucukali. Version 1.2, October 9, 2022")
p <- add_argument(p, "--input", help = "Input STR file with a header, first 3 columns are chrom, start, end, and rest are sample ids with inqH1 & inqH2 STR lengths", nargs = 1)
p <- add_argument(p, "--phenocovar", help = "Phenotype and covariate file with header, first column is sample_id", nargs = "*")
p <- add_argument(p, "--covnames", help = "Covariate names you want to use (if you want to), separated by comma", nargs = "*")
p <- add_argument(p, "--phenotype", help = "Column name of your phenotype of interest variable in the --phenocovar file", nargs = 1)
p <- add_argument(p, "--out", help = "Output file name", nargs = 1)
p <- add_argument(p, "--mode", help = "Select a mode name from following: SUMinqH1H2, MAXinqH1H2, MINinqH1H2", nargs = 1)
p <- add_argument(p, "--missing_cutoff", help = "Defines the call rate cutoff for variants, default is 0.80 meaning that keeping all variants present in at least 80% of individuals", default = "0.80")
p <- add_argument(p, "--outcometype", help = "Select a outcome variable type: binary or continuous", nargs = 1)
p <- add_argument(p, "--binaryOrder", help = "Give the binary phenotype order, comma separated, e.g. Control, Patient will code Control as 0/Group1 and Patient as 1/Group2", nargs = "*")
p <- add_argument(p, "--chr", help = "Indicate chromosome number to be analyzed", nargs = 1)

arg <- parse_args(p)

Version <- "Run association testing for STRs with different modes and options, written by Fahri Kucukali. Version 1.2, October 9, 2022"
print(Version)

calls_file <- fread(arg$input, header = TRUE)

calls_file <- subset(calls_file, chrom == arg$chr)

first_3col <- calls_file[, c(1:3)]

rest <- calls_file[, -c(1:3)]

first_3col$STR_chr_begin_end <- paste(first_3col$chrom, first_3col$begin, first_3col$end, sep = "_")

col_index <- seq(1:ncol(rest))

inqH1 <- as.data.table(rest %>% select(col_index[col_index %% 2 != 0]))

inqH2 <- as.data.table(rest %>% select(col_index[col_index %% 2 == 0]))

colnames(inqH1) <- gsub(".[^.]+$", "", colnames(inqH1))

colnames(inqH2) <- gsub(".[^.]+$", "", colnames(inqH2))

sample_list <- as.data.table(colnames(inqH1))

colnames(sample_list) <- "sample_id"

phenocovar <- fread(arg$phenocovar, header = TRUE)

phenotype <- paste0(arg$phenotype, "")

no_cols <- ncol(phenocovar) + 1

sample_list_wPheno <- left_join(sample_list, phenocovar, by = "sample_id")

outputname <- arg$out

prepare_calls <- function(inqH1, inqH2, strnames, sample_list_wPheno, mode) {
    if (mode == "SUMinqH1H2") {
        calls <- inqH1 + inqH2
    } else if (mode == "MAXinqH1H2") {
        calls <- pmax(inqH1, inqH2)
    } else if (mode == "MINinqH1H2") {
        calls <- pmin(inqH1, inqH2)
    }
    calls_file <- transpose(calls)
    colnames(calls_file) <- strnames
    calls_file <- cbind(sample_list_wPheno, calls_file)
    calls_file <- calls_file[, which(unlist(lapply(calls_file, function(x) !all(is.na(x))))), with = F]
    calls_file <- data.table(data.frame(calls_file)[, which(colMeans(!is.na(data.frame(calls_file))) >= arg$missing_cutoff)])
    calls_file <- calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))
    return(calls_file)
}

calls_file <- prepare_calls(inqH1, inqH2, first_3col$STR_chr_begin_end, sample_list_wPheno, arg$mode)

if (!is.na(arg$covnames)) {
    covlist <- gsub(",", " ", arg$covnames)
    covlist_prepared <- unlist(strsplit(covlist, split = " "))
    if (arg$outcometype == "binary") {
        binaryOrder <- gsub(",", " ", arg$binaryOrder)
        binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))
        calls_file[[phenotype]] <- factor(calls_file[[phenotype]], c(binaryOrder_prepared))
        results_calls_file <- as.data.frame(matrix(0, 1, 11))
        colnames(results_calls_file) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")
        for (i in seq(no_cols, ncol(calls_file), 1)) {
            VariantToBeTested <- as.character(colnames(calls_file)[i])
            selectedtable <- as.data.table(cbind(as.character(calls_file[[phenotype]]), as.numeric(calls_file[[VariantToBeTested]])))
            colnames(selectedtable) <- c(phenotype, VariantToBeTested)
            group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])
            group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])
            AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")
            glm_result <- glm(formula = formulax, data = SUMinqH1H2_calls_file, family = binomial(link = "logit"))
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
            results_calls_file <- rbind.data.frame(results_calls_file, tabular_result)
        }
        results_calls_file <- results_calls_file[-1, ]
        sorted_results_calls_file <- results_calls_file[order(results_calls_file$Pvalue), ]
        write.table(sorted_results_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
    } else if (arg$outcometype == "continuous") {
        results_calls_file <- as.data.frame(matrix(0, 1, 11))
        colnames(results_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")
        for (i in seq(no_cols, ncol(calls_file), 1)) {
            VariantToBeTested <- as.character(colnames(calls_file)[i])
            selectedtable <- as.data.table(cbind(as.character(calls_file[[phenotype]]), as.numeric(calls_file[[VariantToBeTested]])))
            colnames(selectedtable) <- c(phenotype, VariantToBeTested)
            group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])
            group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])
            AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            MaxSize <- round(max(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            MinSize <- round(min(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")
            glm_result <- glm(formula = formulax, data = SUMinqH1H2_calls_file, family = gaussian(link = "identity"))
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
        write.table(sorted_results_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
    }
} else if (is.na(arg$covnames)) {
    if (arg$outcometype == "binary") {
        binaryOrder <- gsub(",", " ", arg$binaryOrder)
        binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))
        calls_file[[phenotype]] <- factor(calls_file[[phenotype]], c(binaryOrder_prepared))
        results_calls_file <- as.data.frame(matrix(0, 1, 11))
        colnames(results_calls_file) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")
        for (i in seq(no_cols, ncol(calls_file), 1)) {
            VariantToBeTested <- as.character(colnames(calls_file)[i])
            selectedtable <- as.data.table(cbind(as.character(calls_file[[phenotype]]), as.numeric(calls_file[[VariantToBeTested]])))
            colnames(selectedtable) <- c(phenotype, VariantToBeTested)
            group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])
            group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])
            AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            formulax <- paste(phenotype, VariantToBeTested, sep = "~")
            glm_result <- glm(formula = formulax, data = SUMinqH1H2_calls_file, family = binomial(link = "logit"))
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
            results_calls_file <- rbind.data.frame(results_calls_file, tabular_result)
        }
        results_calls_file <- results_calls_file[-1, ]
        sorted_results_calls_file <- results_calls_file[order(results_calls_file$Pvalue), ]
        write.table(sorted_results_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
    } else if (arg$outcometype == "continuous") {
        results_calls_file <- as.data.frame(matrix(0, 1, 11))
        colnames(results_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")
        for (i in seq(no_cols, ncol(calls_file), 1)) {
            VariantToBeTested <- as.character(colnames(calls_file)[i])
            selectedtable <- as.data.table(cbind(as.character(calls_file[[phenotype]]), as.numeric(calls_file[[VariantToBeTested]])))
            colnames(selectedtable) <- c(phenotype, VariantToBeTested)
            group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])
            group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])
            AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            MaxSize <- round(max(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            MinSize <- round(min(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            formulax <- paste(phenotype, VariantToBeTested, sep = "~")
            glm_result <- glm(formula = formulax, data = SUMinqH1H2_calls_file, family = gaussian(link = "identity"))
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
        write.table(sorted_results_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
    }
}
