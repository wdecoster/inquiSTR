#!~/miniconda3/bin/Rscript
## Run association testing for STRs with different modes and options

# For now in R, to be converted into Rust (maybe)

# TODOS TO PREVENT USER STUPIDITY - FK: language! :D
# TODO: check if STRmode is only MEAN, MIN or MAX, if not throw early error - FK: done
# TODO: check that if users supply a begin they also supply an end argument - FK: (at least partially) addressed, but it is difficult to account for all combinations
# TODO: either support running the whole file, or check that users at least supply one of chr, chr-begin-end or bed - FK: done, and many more other run options added
# TODO: add a progress bar :-D - FK: done! Using https://github.com/SciViews/svMisc

# FK: Many other functionalities are added, this is a major update

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(argparser)))
suppressWarnings(suppressMessages(library(svMisc)))
suppressWarnings(suppressMessages(library(valr)))

assoc_binary <- function(arg, calls_file, phenotype, no_cols, covariates) {
    binaryOrder <- gsub(",", " ", arg$binaryOrder)
    binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))
    calls_file_selected <- as.data.table(calls_file[calls_file[[phenotype]] %in% c(binaryOrder_prepared), ])
    calls_file_selected[[phenotype]] <- factor(calls_file_selected[[phenotype]], c(binaryOrder_prepared))
    results_calls_file_selected <- as.data.frame(matrix(0, 1, 15))
    colnames(results_calls_file_selected) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "Group1_N", "Group2_N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "Group2_Group1_absAvgSizeDiff", "Group2_Group1_OR_for_absAvgSizeDiff", "model")
    message(paste0("Running association testing for ", (ncol(calls_file_selected) - no_cols) + 1, " qualifying variants..."))
    for (i in seq(no_cols, ncol(calls_file_selected), 1)) {
        VariantToBeTested <- as.character(colnames(calls_file_selected)[i])
        if (!is.na(covariates)) {
        covlist <- gsub(",", " ", covariates)
        covlist_prepared <- unlist(strsplit(covlist, split = " "))
        formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")
        selectedtable <- na.omit(as.data.table(cbind(as.character(calls_file_selected[[phenotype]]),as.numeric(calls_file_selected[[VariantToBeTested]]),  calls_file_selected[,..covlist_prepared])))
        colnames(selectedtable)[1:2] <- c(phenotype, VariantToBeTested)
        group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])
        group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])
        AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_Group1_absAvgSizeDiff <- round(abs(Group2_AvgSize - Group1_AvgSize), digits = 3)
        Group2_N <- nrow(subset(group2, group2[[VariantToBeTested]] != "NaN"))
        Group1_N <- nrow(subset(group1, group1[[VariantToBeTested]] != "NaN"))
        } else {
        formulax <- paste(phenotype, VariantToBeTested, sep = "~")
        selectedtable <- as.data.table(cbind(as.character(calls_file_selected[[phenotype]]), as.numeric(calls_file_selected[[VariantToBeTested]])))
        colnames(selectedtable) <- c(phenotype, VariantToBeTested)
        group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])
        group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])
        AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_Group1_absAvgSizeDiff <- round(abs(Group2_AvgSize - Group1_AvgSize), digits = 3)
        Group2_N <- nrow(subset(group2, group2[[VariantToBeTested]] != "NaN"))
        Group1_N <- nrow(subset(group1, group1[[VariantToBeTested]] != "NaN"))
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
        Group2_Group1_OR_for_absAvgSizeDiff <- round((exp(Group2_Group1_absAvgSizeDiff*log(OR))), digits=3)
        model <- as.character(glm_result$formula)[1]
        tabular_result <- as.data.frame(cbind(Predictors, OR, OR_L95, OR_U95, OR_stdErr, Pvalue, N, Group1_N, Group2_N, AvgSize, Group1_AvgSize, Group2_AvgSize, Group2_Group1_absAvgSizeDiff, Group2_Group1_OR_for_absAvgSizeDiff, model))
        tabular_result <- subset(tabular_result, Predictors == VariantID)
        colnames(tabular_result)[1] <- "VariantID"
        results_calls_file_selected <- rbind.data.frame(results_calls_file_selected, tabular_result)
        progress(i, init = TRUE, progress.bar = FALSE, console = TRUE, gui = FALSE)
        if (i == ncol(calls_file_selected)) {
        message("Done!")
        }}
    results_calls_file_selected <- as.data.table(results_calls_file_selected[-1, ])
    sorted_results_calls_file_selected <- results_calls_file_selected[order(as.numeric(results_calls_file_selected$Pvalue)), ]
    write.table(sorted_results_calls_file_selected, arg$out, sep = "\t", col.names = TRUE, quote = F, row.names = F)
}

assoc_continuous <- function(arg, calls_file, phenotype, no_cols, covariates) {
    results_calls_file <- as.data.frame(matrix(0, 1, 13))
    colnames(results_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "Max_Min_absSizeDiff", "Max_Min_Beta_for_absSizeDiff", "model")
    message(paste0("Running association testing for ", (ncol(calls_file) - no_cols) + 1, " qualifying variants..."))
    for (i in seq(no_cols, ncol(calls_file), 1)) {
        VariantToBeTested <- as.character(colnames(calls_file)[i])
        if (!is.na(covariates)) {
        covlist <- gsub(",", " ", covariates)
        covlist_prepared <- unlist(strsplit(covlist, split = " "))
        formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")
        selectedtable <- na.omit(as.data.table(cbind(as.character(calls_file_selected[[phenotype]]),as.numeric(calls_file_selected[[VariantToBeTested]]),  calls_file_selected[,..covlist_prepared])))
        colnames(selectedtable)[1:2] <- c(phenotype, VariantToBeTested)
        AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        MaxSize <- round(max(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        MinSize <- round(min(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Max_Min_absSizeDiff <- round(abs(max(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE) - min(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE)), digits = 3)
        } else {
        formulax <- paste(phenotype, VariantToBeTested, sep = "~")
        selectedtable <- as.data.table(cbind(as.character(calls_file[[phenotype]]), as.numeric(calls_file[[VariantToBeTested]])))
        colnames(selectedtable) <- c(phenotype, VariantToBeTested)
        AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        MaxSize <- round(max(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        MinSize <- round(min(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Max_Min_absSizeDiff <- round(abs(max(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE) - min(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE)), digits = 3)
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
        Max_Min_Beta_for_absSizeDiff <- round((Max_Min_absSizeDiff*Beta), digits=3)
        model <- as.character(glm_result$formula)[1]
        tabular_result <- as.data.frame(cbind(Predictors, Beta, Beta_L95, Beta_U95, Beta_stdErr, Pvalue, N, AvgSize, MinSize, MaxSize, Max_Min_absSizeDiff, Max_Min_Beta_for_absSizeDiff, model))
        tabular_result <- subset(tabular_result, Predictors == VariantID)
        colnames(tabular_result)[1] <- "VariantID"
        results_calls_file <- rbind.data.frame(results_calls_file, tabular_result)
        progress(i, init = TRUE, progress.bar = FALSE, console = TRUE, gui = FALSE)
        if (i == ncol(calls_file)) {
        message("Done!")
        }}
    results_calls_file <- results_calls_file[-1, ]
    sorted_results_calls_file <- results_calls_file[order(as.numeric(results_calls_file$Pvalue)), ]
    write.table(sorted_results_calls_file, arg$out, sep = "\t", col.names = TRUE, quote = F, row.names = F)
}

assoc_binary_expandedAllele <- function(arg, calls_file, phenotype, no_cols, covariates, expandedAllele) {
    expandedAllele <- as.numeric(arg$expandedAllele)
    binaryOrder <- gsub(",", " ", arg$binaryOrder)
    binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))
    calls_file_selected <- as.data.table(calls_file[calls_file[[phenotype]] %in% c(binaryOrder_prepared), ])
    calls_file_selected[[phenotype]] <- factor(calls_file_selected[[phenotype]], c(binaryOrder_prepared))
    results_calls_file_selected <- as.data.frame(matrix(0, 1, 15))
    colnames(results_calls_file_selected) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", "Group1_N", "Group2_N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "Group2_Group1_absAvgSizeDiff", "Group2_Group1_OR_for_absAvgSizeDiff", "model")
    message(paste0("Running association testing for ", (ncol(calls_file_selected) - no_cols) + 1, " qualifying variants..."))
    for (i in seq(no_cols, ncol(calls_file_selected), 1)) {
        VariantToBeTested <- as.character(colnames(calls_file_selected)[i])
        if (!is.na(covariates)) {
        covlist <- gsub(",", " ", covariates)
        covlist_prepared <- unlist(strsplit(covlist, split = " "))
        selectedtable <- na.omit(as.data.table(cbind(as.character(calls_file_selected[[phenotype]]),as.numeric(calls_file_selected[[VariantToBeTested]]),  calls_file_selected[,..covlist_prepared])))
        colnames(selectedtable)[1:2] <- c(phenotype, VariantToBeTested)
        group2 <- subset(selectedtable, selectedtable[[VariantToBeTested]] >= expandedAllele)
        group1 <- subset(selectedtable, selectedtable[[VariantToBeTested]] < expandedAllele)
        calls_file_selected$expanded_allele_group <- ifelse(calls_file_selected[[VariantToBeTested]] >= expandedAllele, "Expanded", "notExpanded")
        calls_file_selected$expanded_allele_group <- factor(calls_file_selected$expanded_allele_group, c("notExpanded","Expanded"))
        formulax <- paste(phenotype, paste(c("expanded_allele_group", covlist_prepared), collapse = "+"), sep = "~")
        AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_Group1_absAvgSizeDiff <- round(abs(Group2_AvgSize - Group1_AvgSize), digits = 3)
        Group2_N <- nrow(subset(group2, group2[[VariantToBeTested]] != "NaN"))
        Group1_N <- nrow(subset(group1, group1[[VariantToBeTested]] != "NaN"))
        } else {
        selectedtable <- as.data.table(cbind(as.character(calls_file_selected[[phenotype]]), as.numeric(calls_file_selected[[VariantToBeTested]])))
        colnames(selectedtable) <- c(phenotype, VariantToBeTested)
        group2 <- subset(selectedtable, selectedtable[[VariantToBeTested]] >= expandedAllele)
        group1 <- subset(selectedtable, selectedtable[[VariantToBeTested]] < expandedAllele)
        calls_file_selected$expanded_allele_group <- ifelse(calls_file_selected[[VariantToBeTested]] >= expandedAllele, "Expanded", "notExpanded")
        calls_file_selected$expanded_allele_group <- factor(calls_file_selected$expanded_allele_group, c("notExpanded","Expanded"))
        formulax <- paste(phenotype, "expanded_allele_group", sep = "~")
        AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_Group1_absAvgSizeDiff <- round(abs(Group2_AvgSize - Group1_AvgSize), digits = 3)
        Group2_N <- nrow(subset(group2, group2[[VariantToBeTested]] != "NaN"))
        Group1_N <- nrow(subset(group1, group1[[VariantToBeTested]] != "NaN"))
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
        Group2_Group1_OR_for_absAvgSizeDiff <- round((exp(Group2_Group1_absAvgSizeDiff*log(OR))), digits=3)
        model <- as.character(glm_result$formula)[1]
        tabular_result <- as.data.frame(cbind(Predictors, OR, OR_L95, OR_U95, OR_stdErr, Pvalue, N, Group1_N, Group2_N, AvgSize, Group1_AvgSize, Group2_AvgSize, Group2_Group1_absAvgSizeDiff, Group2_Group1_OR_for_absAvgSizeDiff, model))
        tabular_result <- subset(tabular_result, Predictors == VariantID)
        colnames(tabular_result)[1] <- "VariantID"
        tabular_result$VariantID <- paste0(as.character(arg$single_variant),"_ExpandedAllele")
        results_calls_file_selected <- rbind.data.frame(results_calls_file_selected, tabular_result)
        progress(i, init = TRUE, progress.bar = FALSE, console = TRUE, gui = FALSE)
        if (i == ncol(calls_file_selected)) {
        message("Done!")
        }}
    results_calls_file_selected <- as.data.table(results_calls_file_selected[-1, ])
    sorted_results_calls_file_selected <- results_calls_file_selected[order(as.numeric(results_calls_file_selected$Pvalue)), ]
    write.table(sorted_results_calls_file_selected, arg$out, sep = "\t", col.names = TRUE, quote = F, row.names = F)
}

assoc_continuous_expandedAllele <- function(arg, calls_file, phenotype, no_cols, covariates, expandedAllele) {
    expandedAllele <- as.numeric(arg$expandedAllele)
    results_calls_file <- as.data.frame(matrix(0, 1, 19))
    colnames(results_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "Group1_N", "Group2_N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "Group2_Group1_absAvgSizeDiff", "Group2_Group1_Beta_for_absAvgSizeDiff", "MinSize", "MaxSize", "Max_Min_absSizeDiff","Max_Min_Beta_for_absSizeDiff", "model")
    message(paste0("Running association testing for ", (ncol(calls_file) - no_cols) + 1, " qualifying variants..."))
    for (i in seq(no_cols, ncol(calls_file), 1)) {
        VariantToBeTested <- as.character(colnames(calls_file)[i])
        if (!is.na(covariates)) {
        covlist <- gsub(",", " ", covariates)
        covlist_prepared <- unlist(strsplit(covlist, split = " "))
        selectedtable <- na.omit(as.data.table(cbind(as.character(calls_file_selected[[phenotype]]),as.numeric(calls_file_selected[[VariantToBeTested]]),  calls_file_selected[,..covlist_prepared])))
        colnames(selectedtable)[1:2] <- c(phenotype, VariantToBeTested)
        group2 <- subset(selectedtable, selectedtable[[VariantToBeTested]] >= expandedAllele)
        group1 <- subset(selectedtable, selectedtable[[VariantToBeTested]] < expandedAllele)
        calls_file_selected$expanded_allele_group <- ifelse(calls_file_selected[[VariantToBeTested]] >= expandedAllele, "Expanded", "notExpanded")
        calls_file_selected$expanded_allele_group <- factor(calls_file_selected$expanded_allele_group, c("notExpanded","Expanded"))
        formulax <- paste(phenotype, paste(c("expanded_allele_group", covlist_prepared), collapse = "+"), sep = "~")
        AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_Group1_absAvgSizeDiff <- round(abs(Group2_AvgSize - Group1_AvgSize), digits = 3)
        MaxSize <- round(max(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        MinSize <- round(min(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Max_Min_absSizeDiff <- round(abs(max(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE) - min(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE)), digits = 3)
        } else {
        selectedtable <- as.data.table(cbind(as.character(calls_file_selected[[phenotype]]), as.numeric(calls_file_selected[[VariantToBeTested]])))
        colnames(selectedtable) <- c(phenotype, VariantToBeTested)
        group2 <- subset(selectedtable, selectedtable[[VariantToBeTested]] >= expandedAllele)
        group1 <- subset(selectedtable, selectedtable[[VariantToBeTested]] < expandedAllele)
        calls_file_selected$expanded_allele_group <- ifelse(calls_file_selected[[VariantToBeTested]] >= expandedAllele, "Expanded", "notExpanded")
        calls_file_selected$expanded_allele_group <- factor(calls_file_selected$expanded_allele_group, c("notExpanded","Expanded"))
        formulax <- paste(phenotype, "expanded_allele_group", sep = "~")
        AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_Group1_absAvgSizeDiff <- round(abs(Group2_AvgSize - Group1_AvgSize), digits = 3)
        MaxSize <- round(max(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        MinSize <- round(min(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Max_Min_absSizeDiff <- round(abs(max(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE) - min(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE)), digits = 3)
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
        Group2_N <- nrow(subset(group2, !is.na(group2[[VariantToBeTested]])))
        Group1_N <- nrow(subset(group1, !is.na(group1[[VariantToBeTested]])))
        Group2_Group1_Beta_for_absAvgSizeDiff <- round((Group2_Group1_absAvgSizeDiff*Beta), digits=3)
        Max_Min_Beta_for_absSizeDiff <- round((Max_Min_absSizeDiff*Beta), digits=3)
        model <- as.character(glm_result$formula)[1]
        tabular_result <- as.data.frame(cbind(Predictors, Beta, Beta_L95, Beta_U95, Beta_stdErr, Pvalue, N, Group1_N, Group2_N, AvgSize, Group1_AvgSize, Group2_AvgSize, Group2_Group1_absAvgSizeDiff, Group2_Group1_Beta_for_absAvgSizeDiff, MinSize, MaxSize, Max_Min_absSizeDiff, Max_Min_Beta_for_absSizeDiff, model))
        tabular_result <- subset(tabular_result, Predictors == VariantID)
        colnames(tabular_result)[1] <- "VariantID"
        tabular_result$VariantID <- paste0(as.character(arg$single_variant),"_ExpandedAllele")
        results_calls_file <- rbind.data.frame(results_calls_file, tabular_result)
        progress(i, init = TRUE, progress.bar = FALSE, console = TRUE, gui = FALSE)
        if (i == ncol(calls_file)) {
        message("Done!")
        }}
    results_calls_file <- results_calls_file[-1, ]
    sorted_results_calls_file <- results_calls_file[order(as.numeric(results_calls_file$Pvalue)), ]
    write.table(sorted_results_calls_file, arg$out, sep = "\t", col.names = TRUE, quote = F, row.names = F)
}

read_calls_full <- function(input) {
    message("Loading and processing the input file...")
    calls_file <- fread(input, header = TRUE)
    strnames <- paste(calls_file$chrom, calls_file$begin, calls_file$end, sep = "_")
    rest <- calls_file[, -c(1:3)]
    col_index <- seq(1:ncol(rest))
    inqH1 <- as.data.table(rest %>% select(col_index[col_index %% 2 != 0]))
    inqH2 <- as.data.table(rest %>% select(col_index[col_index %% 2 == 0]))
    colnames(inqH1) <- gsub(".[^.]+$", "", colnames(inqH1))
    colnames(inqH2) <- gsub(".[^.]+$", "", colnames(inqH2))
    return(list("H1" = inqH1, "H2" = inqH2, "strnames" = strnames))
}

read_calls_chr <- function(input, chr) {
    message("Loading and processing the input file...")
    calls_file <- fread(input, header = TRUE)
    calls_file <- subset(calls_file, chrom == chr)
    strnames <- paste(calls_file$chrom, calls_file$begin, calls_file$end, sep = "_")
    rest <- calls_file[, -c(1:3)]
    col_index <- seq(1:ncol(rest))
    inqH1 <- as.data.table(rest %>% select(col_index[col_index %% 2 != 0]))
    inqH2 <- as.data.table(rest %>% select(col_index[col_index %% 2 == 0]))
    colnames(inqH1) <- gsub(".[^.]+$", "", colnames(inqH1))
    colnames(inqH2) <- gsub(".[^.]+$", "", colnames(inqH2))
    return(list("H1" = inqH1, "H2" = inqH2, "strnames" = strnames))
}

read_calls_chr_begin_end <- function(input, chr, chr_begin, chr_end) {
    message("Loading and processing the input file...")
    calls_file <- fread(input, header = TRUE)
    calls_file <- subset(calls_file, ((chrom == chr) & (begin >= chr_begin) & (end <= chr_end)))
    strnames <- paste(calls_file$chrom, calls_file$begin, calls_file$end, sep = "_")
    rest <- calls_file[, -c(1:3)]
    col_index <- seq(1:ncol(rest))
    inqH1 <- as.data.table(rest %>% select(col_index[col_index %% 2 != 0]))
    inqH2 <- as.data.table(rest %>% select(col_index[col_index %% 2 == 0]))
    colnames(inqH1) <- gsub(".[^.]+$", "", colnames(inqH1))
    colnames(inqH2) <- gsub(".[^.]+$", "", colnames(inqH2))
    return(list("H1" = inqH1, "H2" = inqH2, "strnames" = strnames))
}

read_calls_bed <- function(input, bed) {
    message("Loading and processing the input file...")
    calls_file <- fread(input, header = TRUE)
    bedfile <- fread(bed, header = FALSE)
    colnames(bedfile) <- c("chrom", "start", "end")
    colnames(calls_file)[2] <- "start"
    intersecttable <- as.data.table(bed_intersect(calls_file, bedfile, suffix = c("", ".y")))
    intersecttable <- intersecttable[, 1:(length(intersecttable) - 3)]
    intersecttable <- subset(intersecttable, !is.na(chrom)) ## WDC here is something wrong - there is no chrom ## FK - no it is fine, we can keep it is. It is just to make sure that there are no "NA" entries in the datatable, based on a chrom column != NA condition (I did it because I've seen it sometimes after bed_intersect). (feel free to remove this comment)
    colnames(intersecttable)[2] <- "begin"
    intersect_strnames <- paste(intersecttable$chrom, intersecttable$begin, intersecttable$end, sep = "_")
    rest <- intersecttable[, -c(1:3)]
    col_index <- seq(1:ncol(rest))
    inqH1 <- as.data.table(rest %>% select(col_index[col_index %% 2 != 0]))
    inqH2 <- as.data.table(rest %>% select(col_index[col_index %% 2 == 0]))
    colnames(inqH1) <- gsub(".[^.]+$", "", colnames(inqH1))
    colnames(inqH2) <- gsub(".[^.]+$", "", colnames(inqH2))
    return(list("H1" = inqH1, "H2" = inqH2, "strnames" = intersect_strnames))
}

read_calls_singleVariant_expandedAllele <- function(input, single_variant) {
    message("Loading and processing the input file...")
    calls_file <- fread(input, header = TRUE)
    single_variant_toAnalyze <- unlist(strsplit(arg$single_variant, split = "_"))
    calls_file <- subset(calls_file, ((chrom == single_variant_toAnalyze[1]) & (begin == single_variant_toAnalyze[2]) & (end == single_variant_toAnalyze[3])))
    strnames <- paste(calls_file$chrom, calls_file$begin, calls_file$end, sep = "_")
    rest <- calls_file[, -c(1:3)]
    col_index <- seq(1:ncol(rest))
    inqH1 <- as.data.table(rest %>% select(col_index[col_index %% 2 != 0]))
    inqH2 <- as.data.table(rest %>% select(col_index[col_index %% 2 == 0]))
    colnames(inqH1) <- gsub(".[^.]+$", "", colnames(inqH1))
    colnames(inqH2) <- gsub(".[^.]+$", "", colnames(inqH2))
    return(list("H1" = inqH1, "H2" = inqH2, "strnames" = strnames))
}

prepare_phenotype <- function(phenofile, phenotype, sample_list) {
    message("Processing the phenotype file...")
    phenocovar <- fread(phenofile, header = TRUE)
    colnames(sample_list) <- "sample_id"
    return(list(
        "sample_list" = left_join(sample_list, phenocovar, by = "sample_id"),
        "phenotype" = phenotype <- paste0(phenotype, ""),
        "no_cols" = ncol(phenocovar) + 1
    ))
}

prepare_calls <- function(calls, sample_list_wPheno, STRmode, missing_cutoff) {
    message("Processing the input file based on the STRmode chosen...")
    if (STRmode == "MEAN") {
        calls2 <- (pmax(calls$H1, calls$H2,na.rm=TRUE) + pmin(inqH1[,2],inqH2[,2],na.rm=TRUE)) / 2
    } else if (STRmode == "MAX") {
        calls2 <- pmax(calls$H1, calls$H2, na.rm = TRUE)
    } else if (STRmode == "MIN") {
        calls2 <- pmin(calls$H1, calls$H2, na.rm = TRUE)
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
    p <- arg_parser("Run association testing for STRs with different modes and options. Version 1.4, October 28, 2022")
    p <- add_argument(p, "--input", help = "inquiSTR input STR file with a header, first 3 columns are chrom, begin, end, and rest are sample ids with inqH1 & inqH2 STR lengths", type = "character", nargs = 1)
    p <- add_argument(p, "--phenocovar", help = "Phenotype and covariate file with header, first column is sample_id", type = "character", nargs = 1)
    p <- add_argument(p, "--covnames", help = "Covariate names you want to use (optional), separated by comma", type = "character", nargs = "*")
    p <- add_argument(p, "--phenotype", help = "Column name of your phenotype of interest variable in the --phenocovar file", type = "character", nargs = 1)
    p <- add_argument(p, "--out", help = "Output file name", type = "character", nargs = 1)
    p <- add_argument(p, "--STRmode", help = "Choose a STRmode from following: MEAN, MAX, MIN; meaning H1+H2 alleles divided by two, maximum of two, or minimum of two. Missing alleles are not considered.", type = "character", nargs = 1)
    p <- add_argument(p, "--missing_cutoff", help = "Defines the call rate cutoff for variants, default is 0.80 meaning that keeping all variants present in at least 80% of individuals. Might mean different things in each of the MEAN, MAX, MIN modes, use carefully.", type = "numeric", default = "0.80")
    p <- add_argument(p, "--outcometype", help = "Select a outcome variable type: binary or continuous", type = "character", nargs = 1)
    p <- add_argument(p, "--binaryOrder", help = "Give the binary phenotype order, comma separated, e.g. Control, Patient will code Control as 0/Group1 and Patient as 1/Group2. This will also be used to further filter your data to only two groups (if you had more than >2 groups for your categorical data)", type = "character", nargs = "*")
    p <- add_argument(p, "--run", help = "Run mode, select among: full (entire file will be analyzed), chromosome (chromosome indicated by --chr will be analyzed), chr_interval (intervals between chr:begin-end defined by --chr --begin --end will be analyzed), bed_interval (an intersection with a bed file defined by --bed will be analyzed), single_variant (a single variant defined by --single_variant will be analyzed; should be combined with --expandedAllele)", type = "character", nargs = 1)
    p <- add_argument(p, "--chr", help = "Indicate chromosome number to be analyzed (with chr prefix). Optional if bed file is provided.", type = "character", nargs = "?")
    p <- add_argument(p, "--chr_begin", help = "Define a begin position (inclusive) for a region of interest (optional, and should be combined with --chr_end)", type = "integer", nargs = "?")
    p <- add_argument(p, "--chr_end", help = "Define a end position (inclusive) for a region of interest (optional, and should be combined with --chr_begin)", type = "integer", nargs = "?")
    p <- add_argument(p, "--bed", help = "A bed file (without a header) with three columns: chromosome (with chr prefix), begin, and end positions for region(s) of interest (optional). valr Rpackage is required - can be installed with mamba install r-valr on conda environment", type = "character", nargs = "?")
    p <- add_argument(p, "--single_variant", help = "A single variant to be analyzed, indicated as chr_begin_end, based on a cut-off value provided by the argument --expandedAllele", type = "character", nargs = "?")
    p <- add_argument(p, "--expandedAllele", help = "A cut-off number (integer or decimal) to define 2 groups with and without expanded allele. STR_length >= expandedAllele is Group 2, STR_length < expandedAllele is Group 1.", type = "integer", nargs = "?")
    Version <- "inquiSTR - STR_regression Rscript Version 1.4, October 28, 2022"
    message(Version)
    return(parse_args(p))
}

arg <- parse_arguments()

# argument checks
if (is.na(arg$input) || is.na(arg$phenocovar) || is.na(arg$phenotype) || is.na(arg$out) || is.na(arg$STRmode) || is.na(arg$outcometype) || is.na(arg$run)) {
  message("Error: exiting because at least one of the following required arguments is missing: --input, --phenocovar, --phenotype, --out, --STRmode, --outcometype, --run")
  q("no")
}

if ((arg$outcometype == "binary") && is.na(arg$binaryOrder)) {
  message("Error: exiting because --binaryOrder argument is missing, please provide it when you use --outcometype binary")
  q("no")
}

if ((arg$run == "chromosome") && is.na(arg$chr)) {
  message("Error: exiting because --chr argument is missing, please provide it when you use --run chromosome")
  q("no")
}

if ((arg$run == "chr_interval") && ((is.na(arg$chr)) || (is.na(arg$chr_begin)) || (is.na(arg$chr_end)))) {
  message("Error: exiting because At least one of the --chr, --chr_begin, or --chr_end arguments is missing, please provide these when you use --run chr_interval")
  q("no")
}

if ((arg$run == "bed_interval") && is.na(arg$bed)) {
  message("Error: exiting because --bed argument, therefore input bed file, is missing; please provide it when you use --run bed_interval")
  q("no")
}

if ((arg$run == "single_variant") && is.na(arg$expandedAllele)) {
  message("Error: exiting because At least one of the two following aruguments is missing: --single_variant --expandedAllele; please provide these when you use --run single_variant")
  q("no")
}

# calls is a list with H1, H2 and strnames attributes
if (arg$run == "chr_interval") {
    calls <- read_calls_chr_begin_end(input = arg$input, chr = arg$chr, chr_begin = arg$chr_begin, chr_end = arg$chr_end)
    message("chr_interval run mode is selected")
} else if (arg$run == "bed_interval") {
    calls <- read_calls_bed(input = arg$input, bed = arg$bed)
    message("bed_interval run mode is selected")
} else if (arg$run == "chromosome") {
    calls <- read_calls_chr(input = arg$input, chr = arg$chr)
    message("chromosome run mode is selected")
} else if (arg$run == "full") {
    calls <- read_calls_full(input = arg$input)
    message("full run mode is selected")
} else if (arg$run == "single_variant") {
    calls <- read_calls_singleVariant_expandedAllele(input = arg$input, single_variant = arg$single_variant)
    message("single_variant run mode is selected (with expandedAllele)")
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
    STRmode = arg$STRmode,
    missing_cutoff = arg$missing_cutoff
)

if (arg$outcometype == "binary" && arg$run != "single_variant") {
    assoc_binary(
        arg = arg,
        calls_file = calls_file,
        phenotype = pheno_info$phenotype,
        no_cols = pheno_info$no_cols,
        covariates = arg$covnames
    )
} else if (arg$outcometype == "continuous" && arg$run != "single_variant") {
    assoc_continuous(
        arg = arg,
        calls_file = calls_file,
        phenotype = pheno_info$phenotype,
        no_cols = pheno_info$no_cols,
        covariates = arg$covnames
    )
} else if (arg$outcometype == "binary" && arg$run == "single_variant") {
    assoc_binary_expandedAllele(
        arg = arg,
        calls_file = calls_file,
        phenotype = pheno_info$phenotype,
        no_cols = pheno_info$no_cols,
        covariates = arg$covnames,
        expandedAllele = arg$expandedAllele
    )
} else if (arg$outcometype == "continuous" && arg$run == "single_variant") {
    assoc_continuous_expandedAllele(
        arg = arg,
        calls_file = calls_file,
        phenotype = pheno_info$phenotype,
        no_cols = pheno_info$no_cols,
        covariates = arg$covnames,
        expandedAllele = arg$expandedAllele
    )
}
