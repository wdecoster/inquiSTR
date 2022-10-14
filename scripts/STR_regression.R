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
p <- add_argument(p, "--mode", help = "Select a mode name from following: inqH1, inqH2, SUMinqH1H2, MEANinqH1H2, MAXinqH1H2, MINinqH1H2", nargs = 1)
p <- add_argument(p, "--missing_cutoff", help = "Defines the call rate cutoff for variants, default is 0.80 meaning that keeping all variants present in at least 80% of individuals", default = "0.80")
p <- add_argument(p, "--outcometype", help = "Select a outcome variable type: binary or continuous", nargs = 1)
p <- add_argument(p, "--binaryOrder", help = "Give the binary phenotype order, comma separated, e.g. Control, Patient will code Control as 0/Group1 and Patient as 1/Group2", nargs = "*")
p <- add_argument(p, "--chr", help = "Indicate chromosome number to be analyzed", nargs = 1)

arg <- parse_args(p)

Version <- "Run association testing for STRs with different modes and options, written by Fahri Kucukali. Version 1.2, October 9, 2022"
print(Version)

calls_file <- fread(arg$input, header = TRUE)

chr_to_analyze <- arg$chr

calls_file <- subset(calls_file, chrom == chr_to_analyze)

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

if (!is.na(arg$covnames)) {
    covlist <- gsub(",", " ", arg$covnames)

    covlist_prepared <- unlist(strsplit(covlist, split = " "))

    if (arg$mode == "inqH1") {
        inqH1_calls_file <- transpose(inqH1)

        colnames(inqH1_calls_file) <- first_3col$STR_chr_begin_end

        inqH1_calls_file <- cbind(sample_list_wPheno, inqH1_calls_file)

        inqH1_calls_file <- inqH1_calls_file[, which(unlist(lapply(inqH1_calls_file, function(x) !all(is.na(x))))), with = F]

        inqH1_calls_file <- data.table(data.frame(inqH1_calls_file)[, which(colMeans(!is.na(data.frame(inqH1_calls_file))) >= arg$missing_cutoff)])

        inqH1_calls_file <- inqH1_calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))

        if (arg$outcometype == "binary") {
            binaryOrder <- gsub(",", " ", arg$binaryOrder)

            binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))

            inqH1_calls_file[[phenotype]] <- factor(inqH1_calls_file[[phenotype]], c(binaryOrder_prepared))

            results_inqH1_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_inqH1_calls_file) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")

            for (i in seq(no_cols, ncol(inqH1_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(inqH1_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(inqH1_calls_file[[phenotype]]), as.numeric(inqH1_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")

                glm_result <- glm(formula = formulax, data = inqH1_calls_file, family = binomial(link = "logit"))

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

                results_inqH1_calls_file <- rbind.data.frame(results_inqH1_calls_file, tabular_result)
            }

            results_inqH1_calls_file <- results_inqH1_calls_file[-1, ]

            sorted_results_inqH1_calls_file <- results_inqH1_calls_file[order(results_inqH1_calls_file$Pvalue), ]

            write.table(sorted_results_inqH1_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        } else if (arg$outcometype == "continuous") {
            results_inqH1_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_inqH1_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")

            for (i in seq(no_cols, ncol(inqH1_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(inqH1_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(inqH1_calls_file[[phenotype]]), as.numeric(inqH1_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MaxSize <- round(max(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MinSize <- round(min(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")

                glm_result <- glm(formula = formulax, data = inqH1_calls_file, family = gaussian(link = "identity"))

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

                results_inqH1_calls_file <- rbind.data.frame(results_inqH1_calls_file, tabular_result)
            }

            results_inqH1_calls_file <- results_inqH1_calls_file[-1, ]

            sorted_results_inqH1_calls_file <- results_inqH1_calls_file[order(results_inqH1_calls_file$Pvalue), ]

            write.table(sorted_results_inqH1_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        }
    } else if (arg$mode == "inqH2") {
        inqH2_calls_file <- transpose(inqH2)

        colnames(inqH2_calls_file) <- first_3col$STR_chr_begin_end

        inqH2_calls_file <- cbind(sample_list_wPheno, inqH2_calls_file)

        inqH2_calls_file <- inqH2_calls_file[, which(unlist(lapply(inqH2_calls_file, function(x) !all(is.na(x))))), with = F]

        inqH2_calls_file <- data.table(data.frame(inqH2_calls_file)[, which(colMeans(!is.na(data.frame(inqH2_calls_file))) >= arg$missing_cutoff)])

        inqH2_calls_file <- inqH2_calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))

        if (arg$outcometype == "binary") {
            binaryOrder <- gsub(",", " ", arg$binaryOrder)

            binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))

            inqH2_calls_file[[phenotype]] <- factor(inqH2_calls_file[[phenotype]], c(binaryOrder_prepared))

            results_inqH2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_inqH2_calls_file) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")

            for (i in seq(no_cols, ncol(inqH2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(inqH2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(inqH2_calls_file[[phenotype]]), as.numeric(inqH2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")

                glm_result <- glm(formula = formulax, data = inqH2_calls_file, family = binomial(link = "logit"))

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

                results_inqH2_calls_file <- rbind.data.frame(results_inqH2_calls_file, tabular_result)
            }

            results_inqH2_calls_file <- results_inqH2_calls_file[-1, ]

            sorted_results_inqH2_calls_file <- results_inqH2_calls_file[order(results_inqH2_calls_file$Pvalue), ]

            write.table(sorted_results_inqH2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        } else if (arg$outcometype == "continuous") {
            results_inqH2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_inqH2_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")

            for (i in seq(no_cols, ncol(inqH2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(inqH2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(inqH2_calls_file[[phenotype]]), as.numeric(inqH2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MaxSize <- round(max(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MinSize <- round(min(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")

                glm_result <- glm(formula = formulax, data = inqH2_calls_file, family = gaussian(link = "identity"))

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

                results_inqH2_calls_file <- rbind.data.frame(results_inqH2_calls_file, tabular_result)
            }

            results_inqH2_calls_file <- results_inqH2_calls_file[-1, ]

            sorted_results_inqH2_calls_file <- results_inqH2_calls_file[order(results_inqH2_calls_file$Pvalue), ]

            write.table(sorted_results_inqH2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        }
    } else if (arg$mode == "SUMinqH1H2") {
        SUMinqH1H2 <- inqH1 + inqH2

        SUMinqH1H2_calls_file <- transpose(SUMinqH1H2)

        colnames(SUMinqH1H2_calls_file) <- first_3col$STR_chr_begin_end

        SUMinqH1H2_calls_file <- cbind(sample_list_wPheno, SUMinqH1H2_calls_file)

        SUMinqH1H2_calls_file <- SUMinqH1H2_calls_file[, which(unlist(lapply(SUMinqH1H2_calls_file, function(x) !all(is.na(x))))), with = F]

        SUMinqH1H2_calls_file <- data.table(data.frame(SUMinqH1H2_calls_file)[, which(colMeans(!is.na(data.frame(SUMinqH1H2_calls_file))) >= arg$missing_cutoff)])

        SUMinqH1H2_calls_file <- SUMinqH1H2_calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))

        if (arg$outcometype == "binary") {
            binaryOrder <- gsub(",", " ", arg$binaryOrder)

            binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))

            SUMinqH1H2_calls_file[[phenotype]] <- factor(SUMinqH1H2_calls_file[[phenotype]], c(binaryOrder_prepared))

            results_SUMinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_SUMinqH1H2_calls_file) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")

            for (i in seq(no_cols, ncol(SUMinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(SUMinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(SUMinqH1H2_calls_file[[phenotype]]), as.numeric(SUMinqH1H2_calls_file[[VariantToBeTested]])))

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

                results_SUMinqH1H2_calls_file <- rbind.data.frame(results_SUMinqH1H2_calls_file, tabular_result)
            }

            results_SUMinqH1H2_calls_file <- results_SUMinqH1H2_calls_file[-1, ]

            sorted_results_SUMinqH1H2_calls_file <- results_SUMinqH1H2_calls_file[order(results_SUMinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_SUMinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        } else if (arg$outcometype == "continuous") {
            results_SUMinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_SUMinqH1H2_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")

            for (i in seq(no_cols, ncol(SUMinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(SUMinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(SUMinqH1H2_calls_file[[phenotype]]), as.numeric(SUMinqH1H2_calls_file[[VariantToBeTested]])))

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

                results_SUMinqH1H2_calls_file <- rbind.data.frame(results_SUMinqH1H2_calls_file, tabular_result)
            }

            results_SUMinqH1H2_calls_file <- results_SUMinqH1H2_calls_file[-1, ]

            sorted_results_SUMinqH1H2_calls_file <- results_SUMinqH1H2_calls_file[order(results_SUMinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_SUMinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        }
    } else if (arg$mode == "MEANinqH1H2") {
        MEANinqH1H2 <- (inqH1 + inqH2) / 2

        MEANinqH1H2_calls_file <- transpose(MEANinqH1H2)

        colnames(MEANinqH1H2_calls_file) <- first_3col$STR_chr_begin_end

        MEANinqH1H2_calls_file <- cbind(sample_list_wPheno, MEANinqH1H2_calls_file)

        MEANinqH1H2_calls_file <- MEANinqH1H2_calls_file[, which(unlist(lapply(MEANinqH1H2_calls_file, function(x) !all(is.na(x))))), with = F]

        MEANinqH1H2_calls_file <- data.table(data.frame(MEANinqH1H2_calls_file)[, which(colMeans(!is.na(data.frame(MEANinqH1H2_calls_file))) >= arg$missing_cutoff)])

        MEANinqH1H2_calls_file <- MEANinqH1H2_calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))

        if (arg$outcometype == "binary") {
            binaryOrder <- gsub(",", " ", arg$binaryOrder)

            binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))

            MEANinqH1H2_calls_file[[phenotype]] <- factor(MEANinqH1H2_calls_file[[phenotype]], c(binaryOrder_prepared))

            results_MEANinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_MEANinqH1H2_calls_file) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")

            for (i in seq(no_cols, ncol(MEANinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(MEANinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(MEANinqH1H2_calls_file[[phenotype]]), as.numeric(MEANinqH1H2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")

                glm_result <- glm(formula = formulax, data = MEANinqH1H2_calls_file, family = binomial(link = "logit"))

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

                results_MEANinqH1H2_calls_file <- rbind.data.frame(results_MEANinqH1H2_calls_file, tabular_result)
            }

            results_MEANinqH1H2_calls_file <- results_MEANinqH1H2_calls_file[-1, ]

            sorted_results_MEANinqH1H2_calls_file <- results_MEANinqH1H2_calls_file[order(results_MEANinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_MEANinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        } else if (arg$outcometype == "continuous") {
            results_MEANinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_MEANinqH1H2_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")

            for (i in seq(no_cols, ncol(MEANinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(MEANinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(MEANinqH1H2_calls_file[[phenotype]]), as.numeric(MEANinqH1H2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MaxSize <- round(max(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MinSize <- round(min(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")

                glm_result <- glm(formula = formulax, data = MEANinqH1H2_calls_file, family = gaussian(link = "identity"))

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

                results_MEANinqH1H2_calls_file <- rbind.data.frame(results_MEANinqH1H2_calls_file, tabular_result)
            }

            results_MEANinqH1H2_calls_file <- results_MEANinqH1H2_calls_file[-1, ]

            sorted_results_MEANinqH1H2_calls_file <- results_MEANinqH1H2_calls_file[order(results_MEANinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_MEANinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        }
    } else if (arg$mode == "MAXinqH1H2") {
        MAXinqH1H2 <- pmax(inqH1, inqH2)

        MAXinqH1H2_calls_file <- transpose(MAXinqH1H2)

        colnames(MAXinqH1H2_calls_file) <- first_3col$STR_chr_begin_end

        MAXinqH1H2_calls_file <- cbind(sample_list_wPheno, MAXinqH1H2_calls_file)

        MAXinqH1H2_calls_file <- MAXinqH1H2_calls_file[, which(unlist(lapply(MAXinqH1H2_calls_file, function(x) !all(is.na(x))))), with = F]

        MAXinqH1H2_calls_file <- data.table(data.frame(MAXinqH1H2_calls_file)[, which(colMeans(!is.na(data.frame(MAXinqH1H2_calls_file))) >= arg$missing_cutoff)])

        MAXinqH1H2_calls_file <- MAXinqH1H2_calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))

        if (arg$outcometype == "binary") {
            binaryOrder <- gsub(",", " ", arg$binaryOrder)

            binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))

            MAXinqH1H2_calls_file[[phenotype]] <- factor(MAXinqH1H2_calls_file[[phenotype]], c(binaryOrder_prepared))

            results_MAXinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_MAXinqH1H2_calls_file) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")

            for (i in seq(no_cols, ncol(MAXinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(MAXinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(MAXinqH1H2_calls_file[[phenotype]]), as.numeric(MAXinqH1H2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")

                glm_result <- glm(formula = formulax, data = MAXinqH1H2_calls_file, family = binomial(link = "logit"))

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

                results_MAXinqH1H2_calls_file <- rbind.data.frame(results_MAXinqH1H2_calls_file, tabular_result)
            }

            results_MAXinqH1H2_calls_file <- results_MAXinqH1H2_calls_file[-1, ]

            sorted_results_MAXinqH1H2_calls_file <- results_MAXinqH1H2_calls_file[order(results_MAXinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_MAXinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        } else if (arg$outcometype == "continuous") {
            results_MAXinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_MAXinqH1H2_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")

            binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))

            for (i in seq(no_cols, ncol(MAXinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(MAXinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(MAXinqH1H2_calls_file[[phenotype]]), as.numeric(MAXinqH1H2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MaxSize <- round(max(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MinSize <- round(min(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")

                glm_result <- glm(formula = formulax, data = MAXinqH1H2_calls_file, family = gaussian(link = "identity"))

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

                results_MAXinqH1H2_calls_file <- rbind.data.frame(results_MAXinqH1H2_calls_file, tabular_result)
            }

            results_MAXinqH1H2_calls_file <- results_MAXinqH1H2_calls_file[-1, ]

            sorted_results_MAXinqH1H2_calls_file <- results_MAXinqH1H2_calls_file[order(results_MAXinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_MAXinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        }
    } else if (arg$mode == "MINinqH1H2") {
        MINinqH1H2 <- pmin(inqH1, inqH2)

        MINinqH1H2_calls_file <- transpose(MINinqH1H2)

        colnames(MINinqH1H2_calls_file) <- first_3col$STR_chr_begin_end

        MINinqH1H2_calls_file <- cbind(sample_list_wPheno, MINinqH1H2_calls_file)

        MINinqH1H2_calls_file <- MINinqH1H2_calls_file[, which(unlist(lapply(MINinqH1H2_calls_file, function(x) !all(is.na(x))))), with = F]

        MINinqH1H2_calls_file <- data.table(data.frame(MINinqH1H2_calls_file)[, which(colMeans(!is.na(data.frame(MINinqH1H2_calls_file))) >= arg$missing_cutoff)])

        MINinqH1H2_calls_file <- MINinqH1H2_calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))

        if (arg$outcometype == "binary") {
            binaryOrder <- gsub(",", " ", arg$binaryOrder)

            binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))

            MINinqH1H2_calls_file[[phenotype]] <- factor(MINinqH1H2_calls_file[[phenotype]], c(binaryOrder_prepared))

            results_MINinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_MINinqH1H2_calls_file) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")

            for (i in seq(no_cols, ncol(MINinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(MINinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(MINinqH1H2_calls_file[[phenotype]]), as.numeric(MINinqH1H2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")

                glm_result <- glm(formula = formulax, data = MINinqH1H2_calls_file, family = binomial(link = "logit"))

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

                results_MINinqH1H2_calls_file <- rbind.data.frame(results_MINinqH1H2_calls_file, tabular_result)
            }

            results_MINinqH1H2_calls_file <- results_MINinqH1H2_calls_file[-1, ]

            sorted_results_MINinqH1H2_calls_file <- results_MINinqH1H2_calls_file[order(results_MINinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_MINinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        } else if (arg$outcometype == "continuous") {
            results_MINinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_MINinqH1H2_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")

            for (i in seq(no_cols, ncol(MINinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(MINinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(MINinqH1H2_calls_file[[phenotype]]), as.numeric(MINinqH1H2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MaxSize <- round(max(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MinSize <- round(min(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")

                glm_result <- glm(formula = formulax, data = MINinqH1H2_calls_file, family = gaussian(link = "identity"))

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

                results_MINinqH1H2_calls_file <- rbind.data.frame(results_MINinqH1H2_calls_file, tabular_result)
            }

            results_MINinqH1H2_calls_file <- results_MINinqH1H2_calls_file[-1, ]

            sorted_results_MINinqH1H2_calls_file <- results_MINinqH1H2_calls_file[order(results_MINinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_MINinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        }
    }
} else if (is.na(arg$covnames)) {
    if (arg$mode == "inqH1") {
        inqH1_calls_file <- transpose(inqH1)

        colnames(inqH1_calls_file) <- first_3col$STR_chr_begin_end

        inqH1_calls_file <- cbind(sample_list_wPheno, inqH1_calls_file)

        inqH1_calls_file <- inqH1_calls_file[, which(unlist(lapply(inqH1_calls_file, function(x) !all(is.na(x))))), with = F]

        inqH1_calls_file <- data.table(data.frame(inqH1_calls_file)[, which(colMeans(!is.na(data.frame(inqH1_calls_file))) >= arg$missing_cutoff)])

        inqH1_calls_file <- inqH1_calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))

        if (arg$outcometype == "binary") {
            binaryOrder <- gsub(",", " ", arg$binaryOrder)

            binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))

            inqH1_calls_file[[phenotype]] <- factor(inqH1_calls_file[[phenotype]], c(binaryOrder_prepared))

            results_inqH1_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_inqH1_calls_file) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")

            for (i in seq(no_cols, ncol(inqH1_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(inqH1_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(inqH1_calls_file[[phenotype]]), as.numeric(inqH1_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, VariantToBeTested, sep = "~")

                glm_result <- glm(formula = formulax, data = inqH1_calls_file, family = binomial(link = "logit"))

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

                results_inqH1_calls_file <- rbind.data.frame(results_inqH1_calls_file, tabular_result)
            }

            results_inqH1_calls_file <- results_inqH1_calls_file[-1, ]

            sorted_results_inqH1_calls_file <- results_inqH1_calls_file[order(results_inqH1_calls_file$Pvalue), ]

            write.table(sorted_results_inqH1_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        } else if (arg$outcometype == "continuous") {
            results_inqH1_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_inqH1_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")

            for (i in seq(no_cols, ncol(inqH1_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(inqH1_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(inqH1_calls_file[[phenotype]]), as.numeric(inqH1_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MaxSize <- round(max(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MinSize <- round(min(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, VariantToBeTested, sep = "~")

                glm_result <- glm(formula = formulax, data = inqH1_calls_file, family = gaussian(link = "identity"))

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

                results_inqH1_calls_file <- rbind.data.frame(results_inqH1_calls_file, tabular_result)
            }

            results_inqH1_calls_file <- results_inqH1_calls_file[-1, ]

            sorted_results_inqH1_calls_file <- results_inqH1_calls_file[order(results_inqH1_calls_file$Pvalue), ]

            write.table(sorted_results_inqH1_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        }
    } else if (arg$mode == "inqH2") {
        inqH2_calls_file <- transpose(inqH2)

        colnames(inqH2_calls_file) <- first_3col$STR_chr_begin_end

        inqH2_calls_file <- cbind(sample_list_wPheno, inqH2_calls_file)

        inqH2_calls_file <- inqH2_calls_file[, which(unlist(lapply(inqH2_calls_file, function(x) !all(is.na(x))))), with = F]

        inqH2_calls_file <- data.table(data.frame(inqH2_calls_file)[, which(colMeans(!is.na(data.frame(inqH2_calls_file))) >= arg$missing_cutoff)])

        inqH2_calls_file <- inqH2_calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))

        if (arg$outcometype == "binary") {
            binaryOrder <- gsub(",", " ", arg$binaryOrder)

            binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))

            inqH2_calls_file[[phenotype]] <- factor(inqH2_calls_file[[phenotype]], c(binaryOrder_prepared))

            results_inqH2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_inqH2_calls_file) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")

            for (i in seq(no_cols, ncol(inqH2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(inqH2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(inqH2_calls_file[[phenotype]]), as.numeric(inqH2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, VariantToBeTested, sep = "~")

                glm_result <- glm(formula = formulax, data = inqH2_calls_file, family = binomial(link = "logit"))

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

                results_inqH2_calls_file <- rbind.data.frame(results_inqH2_calls_file, tabular_result)
            }

            results_inqH2_calls_file <- results_inqH2_calls_file[-1, ]

            sorted_results_inqH2_calls_file <- results_inqH2_calls_file[order(results_inqH2_calls_file$Pvalue), ]

            write.table(sorted_results_inqH2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        } else if (arg$outcometype == "continuous") {
            results_inqH2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_inqH2_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")

            for (i in seq(no_cols, ncol(inqH2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(inqH2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(inqH2_calls_file[[phenotype]]), as.numeric(inqH2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MaxSize <- round(max(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MinSize <- round(min(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, VariantToBeTested, sep = "~")

                glm_result <- glm(formula = formulax, data = inqH2_calls_file, family = gaussian(link = "identity"))

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

                results_inqH2_calls_file <- rbind.data.frame(results_inqH2_calls_file, tabular_result)
            }

            results_inqH2_calls_file <- results_inqH2_calls_file[-1, ]

            sorted_results_inqH2_calls_file <- results_inqH2_calls_file[order(results_inqH2_calls_file$Pvalue), ]

            write.table(sorted_results_inqH2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        }
    } else if (arg$mode == "SUMinqH1H2") {
        SUMinqH1H2 <- inqH1 + inqH2

        SUMinqH1H2_calls_file <- transpose(SUMinqH1H2)

        colnames(SUMinqH1H2_calls_file) <- first_3col$STR_chr_begin_end

        SUMinqH1H2_calls_file <- cbind(sample_list_wPheno, SUMinqH1H2_calls_file)

        SUMinqH1H2_calls_file <- SUMinqH1H2_calls_file[, which(unlist(lapply(SUMinqH1H2_calls_file, function(x) !all(is.na(x))))), with = F]

        SUMinqH1H2_calls_file <- data.table(data.frame(SUMinqH1H2_calls_file)[, which(colMeans(!is.na(data.frame(SUMinqH1H2_calls_file))) >= arg$missing_cutoff)])

        SUMinqH1H2_calls_file <- SUMinqH1H2_calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))

        if (arg$outcometype == "binary") {
            binaryOrder <- gsub(",", " ", arg$binaryOrder)

            binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))

            SUMinqH1H2_calls_file[[phenotype]] <- factor(SUMinqH1H2_calls_file[[phenotype]], c(binaryOrder_prepared))

            results_SUMinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_SUMinqH1H2_calls_file) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")

            for (i in seq(no_cols, ncol(SUMinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(SUMinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(SUMinqH1H2_calls_file[[phenotype]]), as.numeric(SUMinqH1H2_calls_file[[VariantToBeTested]])))

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

                results_SUMinqH1H2_calls_file <- rbind.data.frame(results_SUMinqH1H2_calls_file, tabular_result)
            }

            results_SUMinqH1H2_calls_file <- results_SUMinqH1H2_calls_file[-1, ]

            sorted_results_SUMinqH1H2_calls_file <- results_SUMinqH1H2_calls_file[order(results_SUMinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_SUMinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        } else if (arg$outcometype == "continuous") {
            results_SUMinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_SUMinqH1H2_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")

            for (i in seq(no_cols, ncol(SUMinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(SUMinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(SUMinqH1H2_calls_file[[phenotype]]), as.numeric(SUMinqH1H2_calls_file[[VariantToBeTested]])))

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

                results_SUMinqH1H2_calls_file <- rbind.data.frame(results_SUMinqH1H2_calls_file, tabular_result)
            }

            results_SUMinqH1H2_calls_file <- results_SUMinqH1H2_calls_file[-1, ]

            sorted_results_SUMinqH1H2_calls_file <- results_SUMinqH1H2_calls_file[order(results_SUMinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_SUMinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        }
    } else if (arg$mode == "MEANinqH1H2") {
        MEANinqH1H2 <- (inqH1 + inqH2) / 2

        MEANinqH1H2_calls_file <- transpose(MEANinqH1H2)

        colnames(MEANinqH1H2_calls_file) <- first_3col$STR_chr_begin_end

        MEANinqH1H2_calls_file <- cbind(sample_list_wPheno, MEANinqH1H2_calls_file)

        MEANinqH1H2_calls_file <- MEANinqH1H2_calls_file[, which(unlist(lapply(MEANinqH1H2_calls_file, function(x) !all(is.na(x))))), with = F]

        MEANinqH1H2_calls_file <- data.table(data.frame(MEANinqH1H2_calls_file)[, which(colMeans(!is.na(data.frame(MEANinqH1H2_calls_file))) >= arg$missing_cutoff)])

        MEANinqH1H2_calls_file <- MEANinqH1H2_calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))

        if (arg$outcometype == "binary") {
            binaryOrder <- gsub(",", " ", arg$binaryOrder)

            binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))

            MEANinqH1H2_calls_file[[phenotype]] <- factor(MEANinqH1H2_calls_file[[phenotype]], c(binaryOrder_prepared))

            results_MEANinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_MEANinqH1H2_calls_file) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")

            for (i in seq(no_cols, ncol(MEANinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(MEANinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(MEANinqH1H2_calls_file[[phenotype]]), as.numeric(MEANinqH1H2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, VariantToBeTested, sep = "~")

                glm_result <- glm(formula = formulax, data = MEANinqH1H2_calls_file, family = binomial(link = "logit"))

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

                results_MEANinqH1H2_calls_file <- rbind.data.frame(results_MEANinqH1H2_calls_file, tabular_result)
            }

            results_MEANinqH1H2_calls_file <- results_MEANinqH1H2_calls_file[-1, ]

            sorted_results_MEANinqH1H2_calls_file <- results_MEANinqH1H2_calls_file[order(results_MEANinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_MEANinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        } else if (arg$outcometype == "continuous") {
            results_MEANinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_MEANinqH1H2_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")

            for (i in seq(no_cols, ncol(MEANinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(MEANinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(MEANinqH1H2_calls_file[[phenotype]]), as.numeric(MEANinqH1H2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MaxSize <- round(max(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MinSize <- round(min(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, VariantToBeTested, sep = "~")

                glm_result <- glm(formula = formulax, data = MEANinqH1H2_calls_file, family = gaussian(link = "identity"))

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

                results_MEANinqH1H2_calls_file <- rbind.data.frame(results_MEANinqH1H2_calls_file, tabular_result)
            }

            results_MEANinqH1H2_calls_file <- results_MEANinqH1H2_calls_file[-1, ]

            sorted_results_MEANinqH1H2_calls_file <- results_MEANinqH1H2_calls_file[order(results_MEANinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_MEANinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        }
    } else if (arg$mode == "MAXinqH1H2") {
        MAXinqH1H2 <- pmax(inqH1, inqH2)

        MAXinqH1H2_calls_file <- transpose(MAXinqH1H2)

        colnames(MAXinqH1H2_calls_file) <- first_3col$STR_chr_begin_end

        MAXinqH1H2_calls_file <- cbind(sample_list_wPheno, MAXinqH1H2_calls_file)

        MAXinqH1H2_calls_file <- MAXinqH1H2_calls_file[, which(unlist(lapply(MAXinqH1H2_calls_file, function(x) !all(is.na(x))))), with = F]

        MAXinqH1H2_calls_file <- data.table(data.frame(MAXinqH1H2_calls_file)[, which(colMeans(!is.na(data.frame(MAXinqH1H2_calls_file))) >= arg$missing_cutoff)])

        MAXinqH1H2_calls_file <- MAXinqH1H2_calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))

        if (arg$outcometype == "binary") {
            binaryOrder <- gsub(",", " ", arg$binaryOrder)

            binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))

            MAXinqH1H2_calls_file[[phenotype]] <- factor(MAXinqH1H2_calls_file[[phenotype]], c(binaryOrder_prepared))

            results_MAXinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_MAXinqH1H2_calls_file) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")

            for (i in seq(no_cols, ncol(MAXinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(MAXinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(MAXinqH1H2_calls_file[[phenotype]]), as.numeric(MAXinqH1H2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, VariantToBeTested, sep = "~")

                glm_result <- glm(formula = formulax, data = MAXinqH1H2_calls_file, family = binomial(link = "logit"))

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

                results_MAXinqH1H2_calls_file <- rbind.data.frame(results_MAXinqH1H2_calls_file, tabular_result)
            }

            results_MAXinqH1H2_calls_file <- results_MAXinqH1H2_calls_file[-1, ]

            sorted_results_MAXinqH1H2_calls_file <- results_MAXinqH1H2_calls_file[order(results_MAXinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_MAXinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        } else if (arg$outcometype == "continuous") {
            results_MAXinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_MAXinqH1H2_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")

            for (i in seq(no_cols, ncol(MAXinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(MAXinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(MAXinqH1H2_calls_file[[phenotype]]), as.numeric(MAXinqH1H2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MaxSize <- round(max(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MinSize <- round(min(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, VariantToBeTested, sep = "~")

                glm_result <- glm(formula = formulax, data = MAXinqH1H2_calls_file, family = gaussian(link = "identity"))

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

                results_MAXinqH1H2_calls_file <- rbind.data.frame(results_MAXinqH1H2_calls_file, tabular_result)
            }

            results_MAXinqH1H2_calls_file <- results_MAXinqH1H2_calls_file[-1, ]

            sorted_results_MAXinqH1H2_calls_file <- results_MAXinqH1H2_calls_file[order(results_MAXinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_MAXinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        }
    } else if (arg$mode == "MINinqH1H2") {
        MINinqH1H2 <- pmin(inqH1, inqH2)

        MINinqH1H2_calls_file <- transpose(MINinqH1H2)

        colnames(MINinqH1H2_calls_file) <- first_3col$STR_chr_begin_end

        MINinqH1H2_calls_file <- cbind(sample_list_wPheno, MINinqH1H2_calls_file)

        MINinqH1H2_calls_file <- MINinqH1H2_calls_file[, which(unlist(lapply(MINinqH1H2_calls_file, function(x) !all(is.na(x))))), with = F]

        MINinqH1H2_calls_file <- data.table(data.frame(MINinqH1H2_calls_file)[, which(colMeans(!is.na(data.frame(MINinqH1H2_calls_file))) >= arg$missing_cutoff)])

        MINinqH1H2_calls_file <- MINinqH1H2_calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))

        if (arg$outcometype == "binary") {
            binaryOrder <- gsub(",", " ", arg$binaryOrder)

            binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))

            MINinqH1H2_calls_file[[phenotype]] <- factor(MINinqH1H2_calls_file[[phenotype]], c(binaryOrder_prepared))

            results_MINinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_MINinqH1H2_calls_file) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "model")

            for (i in seq(no_cols, ncol(MINinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(MINinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(MINinqH1H2_calls_file[[phenotype]]), as.numeric(MINinqH1H2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, VariantToBeTested, sep = "~")

                glm_result <- glm(formula = formulax, data = MINinqH1H2_calls_file, family = binomial(link = "logit"))

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

                results_MINinqH1H2_calls_file <- rbind.data.frame(results_MINinqH1H2_calls_file, tabular_result)
            }

            results_MINinqH1H2_calls_file <- results_MINinqH1H2_calls_file[-1, ]

            sorted_results_MINinqH1H2_calls_file <- results_MINinqH1H2_calls_file[order(results_MINinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_MINinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        } else if (arg$outcometype == "continuous") {
            results_MINinqH1H2_calls_file <- as.data.frame(matrix(0, 1, 11))

            colnames(results_MINinqH1H2_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "model")

            for (i in seq(no_cols, ncol(MINinqH1H2_calls_file), 1)) {
                VariantToBeTested <- as.character(colnames(MINinqH1H2_calls_file)[i])

                selectedtable <- as.data.table(cbind(as.character(MINinqH1H2_calls_file[[phenotype]]), as.numeric(MINinqH1H2_calls_file[[VariantToBeTested]])))

                colnames(selectedtable) <- c(phenotype, VariantToBeTested)

                group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])

                group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])

                AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MaxSize <- round(max(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                MinSize <- round(min(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)

                formulax <- paste(phenotype, VariantToBeTested, sep = "~")

                glm_result <- glm(formula = formulax, data = MINinqH1H2_calls_file, family = gaussian(link = "identity"))

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

                results_MINinqH1H2_calls_file <- rbind.data.frame(results_MINinqH1H2_calls_file, tabular_result)
            }

            results_MINinqH1H2_calls_file <- results_MINinqH1H2_calls_file[-1, ]

            sorted_results_MINinqH1H2_calls_file <- results_MINinqH1H2_calls_file[order(results_MINinqH1H2_calls_file$Pvalue), ]

            write.table(sorted_results_MINinqH1H2_calls_file, outputname, sep = "\t", col.names = TRUE, quote = F, row.names = F)
        }
    }
}
