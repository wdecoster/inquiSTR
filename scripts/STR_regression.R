#!~/miniconda3/bin/Rscript
## Run association testing for STRs with different modes and options
# For now in R, to be converted into Rust (maybe)

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(argparser)))
suppressWarnings(suppressMessages(library(svMisc)))
suppressWarnings(suppressMessages(library(valr)))

assoc_binary <- function(arg, calls_file, phenotype, no_cols, covariates, missing_cutoff) {
    binaryOrder_prepared <- unlist(strsplit(arg$binaryOrder, split = ","))
    # make sure that all values in binaryOrder_prepared are in the calls_file[[phenotype]] column
    for (i in binaryOrder_prepared) {
        if (!i %in% calls_file[[phenotype]]) {
            stop(paste0("The value ", i, " in binaryOrder is not present in the phenotype column of the input file."))
        }
    }
    calls_file_selected <- as.data.table(calls_file[calls_file[[phenotype]] %in% c(binaryOrder_prepared), ])
    calls_file_selected[[phenotype]] <- factor(calls_file_selected[[phenotype]], c(binaryOrder_prepared))
    calls_file_selected <- calls_file_selected[, which(unlist(lapply(calls_file_selected, function(x) !all(is.na(x))))), with = FALSE]
    calls_file_selected <- data.table(data.frame(calls_file_selected)[, which(colMeans(!is.na(data.frame(calls_file_selected))) >= missing_cutoff)])
    calls_file_selected <- calls_file_selected %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))
    results_calls_file_selected <- as.data.frame(matrix(0, 1, 16))
    colnames(results_calls_file_selected) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", paste0(binaryOrder_prepared[1], "_N"), paste0(binaryOrder_prepared[2], "_N"), "AvgSize", paste0(binaryOrder_prepared[1], "_AvgSize"), paste0(binaryOrder_prepared[2], "_AvgSize"), paste0(binaryOrder_prepared[2], "_", binaryOrder_prepared[1], "_absAvgSizeDiff"), paste0(binaryOrder_prepared[2], "_", binaryOrder_prepared[1], "_OR_for_absAvgSizeDiff"), "model", "binaryOrder")
    if (!arg$quiet) {
        message(paste0("Running association testing for ", (ncol(calls_file_selected) - no_cols) + 1, " qualifying variants..."))
    }
    for (i in seq(no_cols, ncol(calls_file_selected), 1)) {
        VariantToBeTested <- as.character(colnames(calls_file_selected)[i])
        if (!is.na(covariates)) {
            covlist <- gsub(",", " ", covariates)
            covlist_prepared <- unlist(strsplit(covlist, split = " "))
            formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")
            selectedtable <- na.omit(as.data.table(cbind(as.character(calls_file_selected[[phenotype]]), as.numeric(calls_file_selected[[VariantToBeTested]]), calls_file_selected[, ..covlist_prepared])))
        } else {
            formulax <- paste(phenotype, VariantToBeTested, sep = "~")
            selectedtable <- as.data.table(cbind(as.character(calls_file_selected[[phenotype]]), as.numeric(calls_file_selected[[VariantToBeTested]])))
        }
        colnames(selectedtable)[1:2] <- c(phenotype, VariantToBeTested)
        group2 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[2])
        group1 <- subset(selectedtable, selectedtable[[phenotype]] == binaryOrder_prepared[1])
        AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
        Group2_Group1_absAvgSizeDiff <- round(abs(Group2_AvgSize - Group1_AvgSize), digits = 3)
        Group2_N <- nrow(subset(group2, group2[[VariantToBeTested]] != "NaN"))
        Group1_N <- nrow(subset(group1, group1[[VariantToBeTested]] != "NaN"))
        binaryOrderInTable <- arg$binaryOrder
        glm_result <- glm(formula = formulax, data = calls_file_selected, family = binomial(link = "logit"))
        Predictors <- names(glm_result$coefficients)
        VariantID <- names(glm_result$coefficients)[2]
        OR <- round(as.numeric(exp(glm_result$coefficients)), digits = 3)
        OR_L95 <- round(as.numeric(exp(confint.default(glm_result, level = 0.95)[, 1])), digits = 3)
        OR_U95 <- round(as.numeric(exp(confint.default(glm_result, level = 0.95)[, 2])), digits = 3)
        OR_stdErr <- round(as.numeric(coef(summary(glm_result))[, "Std. Error"]), digits = 3)
        Pvalue <- as.numeric(coef(summary(glm_result))[, "Pr(>|z|)"])
        N <- nobs(glm_result)
        Group2_Group1_OR_for_absAvgSizeDiff <- round((exp(Group2_Group1_absAvgSizeDiff * log(OR))), digits = 3)
        model <- as.character(glm_result$formula)[1]
        # instead of creating data frames and rbind'ing them every iteration, it would be better to just print to stdout
        # we don't get sorting then, but that's not too bad
        tabular_result <- as.data.frame(cbind(Predictors, OR, OR_L95, OR_U95, OR_stdErr, Pvalue, N, Group1_N, Group2_N, AvgSize, Group1_AvgSize, Group2_AvgSize, Group2_Group1_absAvgSizeDiff, Group2_Group1_OR_for_absAvgSizeDiff, model, binaryOrderInTable))
        tabular_result <- subset(tabular_result, Predictors == VariantID)
        colnames(tabular_result) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", paste0(binaryOrder_prepared[1], "_N"), paste0(binaryOrder_prepared[2], "_N"), "AvgSize", paste0(binaryOrder_prepared[1], "_AvgSize"), paste0(binaryOrder_prepared[2], "_AvgSize"), paste0(binaryOrder_prepared[2], "_", binaryOrder_prepared[1], "_absAvgSizeDiff"), paste0(binaryOrder_prepared[2], "_", binaryOrder_prepared[1], "_OR_for_absAvgSizeDiff"), "model", "binaryOrder")
        results_calls_file_selected <- rbind.data.frame(results_calls_file_selected, tabular_result)
        svMisc::progress(i, init = TRUE, progress.bar = FALSE, console = TRUE, gui = FALSE)
        if (i == ncol(calls_file_selected)) {
            message("Done!")
        }
    }
    results_calls_file_selected <- as.data.table(results_calls_file_selected[-1, ])
    sorted_results_calls_file_selected <- results_calls_file_selected[order(as.numeric(results_calls_file_selected$Pvalue)), ]
    write.table(sorted_results_calls_file_selected, arg$out, sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
}

assoc_continuous <- function(arg, calls_file, phenotype, no_cols, covariates, missing_cutoff) {
    calls_file <- calls_file[, which(unlist(lapply(calls_file, function(x) !all(is.na(x))))), with = FALSE]
    calls_file <- data.table(data.frame(calls_file)[, which(colMeans(!is.na(data.frame(calls_file))) >= missing_cutoff)])
    calls_file <- calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))
    results_calls_file <- as.data.frame(matrix(0, 1, 13))
    colnames(results_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "AvgSize", "MinSize", "MaxSize", "Max_Min_absSizeDiff", "Max_Min_Beta_for_absSizeDiff", "model")
    if (!arg$quiet) {
        message(paste0("Running association testing for ", (ncol(calls_file) - no_cols) + 1, " qualifying variants..."))
    }
    for (i in seq(no_cols, ncol(calls_file), 1)) {
        VariantToBeTested <- as.character(colnames(calls_file)[i])
        if (!is.na(covariates)) {
            covlist <- gsub(",", " ", covariates)
            covlist_prepared <- unlist(strsplit(covlist, split = " "))
            formulax <- paste(phenotype, paste(c(VariantToBeTested, covlist_prepared), collapse = "+"), sep = "~")
            selectedtable <- na.omit(as.data.table(cbind(as.character(calls_file_selected[[phenotype]]), as.numeric(calls_file_selected[[VariantToBeTested]]), calls_file_selected[, ..covlist_prepared])))
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
        Max_Min_Beta_for_absSizeDiff <- round((Max_Min_absSizeDiff * Beta), digits = 3)
        model <- as.character(glm_result$formula)[1]
        tabular_result <- as.data.frame(cbind(Predictors, Beta, Beta_L95, Beta_U95, Beta_stdErr, Pvalue, N, AvgSize, MinSize, MaxSize, Max_Min_absSizeDiff, Max_Min_Beta_for_absSizeDiff, model))
        tabular_result <- subset(tabular_result, Predictors == VariantID)
        colnames(tabular_result)[1] <- "VariantID"
        results_calls_file <- rbind.data.frame(results_calls_file, tabular_result)
        svMisc::progress(i, init = TRUE, progress.bar = FALSE, console = TRUE, gui = FALSE)
        if (i == ncol(calls_file)) {
            message("Done!")
        }
    }
    results_calls_file <- results_calls_file[-1, ]
    sorted_results_calls_file <- results_calls_file[order(as.numeric(results_calls_file$Pvalue)), ]
    write.table(sorted_results_calls_file, arg$out, sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
}

assoc_binary_expandedAllele <- function(arg, calls_file, phenotype, no_cols, covariates, expandedAllele, missing_cutoff) {
    expandedAllele <- as.numeric(arg$expandedAllele)
    binaryOrder <- gsub(",", " ", arg$binaryOrder)
    binaryOrder_prepared <- unlist(strsplit(binaryOrder, split = " "))
    calls_file_selected <- as.data.table(calls_file[calls_file[[phenotype]] %in% c(binaryOrder_prepared), ])
    calls_file_selected[[phenotype]] <- factor(calls_file_selected[[phenotype]], c(binaryOrder_prepared))
    calls_file_selected <- calls_file_selected[, which(unlist(lapply(calls_file_selected, function(x) !all(is.na(x))))), with = F]
    calls_file_selected <- data.table(data.frame(calls_file_selected)[, which(colMeans(!is.na(data.frame(calls_file_selected))) >= missing_cutoff)])
    calls_file_selected <- calls_file_selected %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))
    results_calls_file_selected <- as.data.frame(matrix(0, 1, 16))
    colnames(results_calls_file_selected) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", paste0(binaryOrder_prepared[1], "_N"), paste0(binaryOrder_prepared[2], "_N"), "AvgSize", paste0(binaryOrder_prepared[1], "_AvgSize"), paste0(binaryOrder_prepared[2], "_AvgSize"), paste0(binaryOrder_prepared[2], "_", binaryOrder_prepared[1], "_absAvgSizeDiff"), paste0(binaryOrder_prepared[2], "_", binaryOrder_prepared[1], "_OR_for_absAvgSizeDiff"), "model", "binaryOrder")
    if (!arg$quiet) {
        message(paste0("Running association testing for ", (ncol(calls_file_selected) - no_cols) + 1, " qualifying variants..."))
    }
    for (i in seq(no_cols, ncol(calls_file_selected), 1)) {
        VariantToBeTested <- as.character(colnames(calls_file_selected)[i])
        if (!is.na(covariates)) {
            covlist <- gsub(",", " ", covariates)
            covlist_prepared <- unlist(strsplit(covlist, split = " "))
            selectedtable <- na.omit(as.data.table(cbind(as.character(calls_file_selected[[phenotype]]), as.numeric(calls_file_selected[[VariantToBeTested]]), calls_file_selected[, ..covlist_prepared])))
            colnames(selectedtable)[1:2] <- c(phenotype, VariantToBeTested)
            group2 <- subset(selectedtable, selectedtable[[VariantToBeTested]] >= expandedAllele)
            group1 <- subset(selectedtable, selectedtable[[VariantToBeTested]] < expandedAllele)
            calls_file_selected$expanded_allele_group <- ifelse(calls_file_selected[[VariantToBeTested]] >= expandedAllele, "Expanded", "notExpanded")
            calls_file_selected$expanded_allele_group <- factor(calls_file_selected$expanded_allele_group, c("notExpanded", "Expanded"))
            formulax <- paste(phenotype, paste(c("expanded_allele_group", covlist_prepared), collapse = "+"), sep = "~")
            AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            Group2_Group1_absAvgSizeDiff <- round(abs(Group2_AvgSize - Group1_AvgSize), digits = 3)
            Group2_N <- nrow(subset(group2, group2[[VariantToBeTested]] != "NaN"))
            Group1_N <- nrow(subset(group1, group1[[VariantToBeTested]] != "NaN"))
            binaryOrderInTable <- arg$binaryOrder
        } else {
            selectedtable <- as.data.table(cbind(as.character(calls_file_selected[[phenotype]]), as.numeric(calls_file_selected[[VariantToBeTested]])))
            colnames(selectedtable) <- c(phenotype, VariantToBeTested)
            group2 <- subset(selectedtable, selectedtable[[VariantToBeTested]] >= expandedAllele)
            group1 <- subset(selectedtable, selectedtable[[VariantToBeTested]] < expandedAllele)
            calls_file_selected$expanded_allele_group <- ifelse(calls_file_selected[[VariantToBeTested]] >= expandedAllele, "Expanded", "notExpanded")
            calls_file_selected$expanded_allele_group <- factor(calls_file_selected$expanded_allele_group, c("notExpanded", "Expanded"))
            formulax <- paste(phenotype, "expanded_allele_group", sep = "~")
            AvgSize <- round(mean(as.numeric(selectedtable[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            Group2_AvgSize <- round(mean(as.numeric(group2[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            Group1_AvgSize <- round(mean(as.numeric(group1[[VariantToBeTested]]), na.rm = TRUE), digits = 3)
            Group2_Group1_absAvgSizeDiff <- round(abs(Group2_AvgSize - Group1_AvgSize), digits = 3)
            Group2_N <- nrow(subset(group2, group2[[VariantToBeTested]] != "NaN"))
            Group1_N <- nrow(subset(group1, group1[[VariantToBeTested]] != "NaN"))
            binaryOrderInTable <- arg$binaryOrder
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
        Group2_Group1_OR_for_absAvgSizeDiff <- round((exp(Group2_Group1_absAvgSizeDiff * log(OR))), digits = 3)
        model <- as.character(glm_result$formula)[1]
        tabular_result <- as.data.frame(cbind(Predictors, OR, OR_L95, OR_U95, OR_stdErr, Pvalue, N, Group1_N, Group2_N, AvgSize, Group1_AvgSize, Group2_AvgSize, Group2_Group1_absAvgSizeDiff, Group2_Group1_OR_for_absAvgSizeDiff, model, binaryOrderInTable))
        tabular_result <- subset(tabular_result, Predictors == VariantID)
        colnames(tabular_result) <- c("VariantID", "OR", "OR_L95", "OR_U95", "OR_stdErr", "Pvalue", "N", paste0(binaryOrder_prepared[1], "_N"), paste0(binaryOrder_prepared[2], "_N"), "AvgSize", paste0(binaryOrder_prepared[1], "_AvgSize"), paste0(binaryOrder_prepared[2], "_AvgSize"), paste0(binaryOrder_prepared[2], "_", binaryOrder_prepared[1], "_absAvgSizeDiff"), paste0(binaryOrder_prepared[2], "_", binaryOrder_prepared[1], "_OR_for_absAvgSizeDiff"), "model", "binaryOrder")
        tabular_result$VariantID <- paste0(as.character(arg$single_variant), "_ExpandedAllele")
        results_calls_file_selected <- rbind.data.frame(results_calls_file_selected, tabular_result)
        svMisc::progress(i, init = TRUE, progress.bar = FALSE, console = TRUE, gui = FALSE)
        if (i == ncol(calls_file_selected)) {
            message("Done!")
        }
    }
    results_calls_file_selected <- as.data.table(results_calls_file_selected[-1, ])
    sorted_results_calls_file_selected <- results_calls_file_selected[order(as.numeric(results_calls_file_selected$Pvalue)), ]
    write.table(sorted_results_calls_file_selected, arg$out, sep = "\t", col.names = TRUE, quote = F, row.names = F)
}

assoc_continuous_expandedAllele <- function(arg, calls_file, phenotype, no_cols, covariates, expandedAllele, missing_cutoff) {
    expandedAllele <- as.numeric(arg$expandedAllele)
    calls_file <- calls_file[, which(unlist(lapply(calls_file, function(x) !all(is.na(x))))), with = F]
    calls_file <- data.table(data.frame(calls_file)[, which(colMeans(!is.na(data.frame(calls_file))) >= missing_cutoff)])
    calls_file <- calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))
    results_calls_file <- as.data.frame(matrix(0, 1, 19))
    colnames(results_calls_file) <- c("VariantID", "Beta", "Beta_L95", "Beta_U95", "Beta_stdErr", "Pvalue", "N", "Group1_N", "Group2_N", "AvgSize", "Group1_AvgSize", "Group2_AvgSize", "Group2_Group1_absAvgSizeDiff", "Group2_Group1_Beta_for_absAvgSizeDiff", "MinSize", "MaxSize", "Max_Min_absSizeDiff", "Max_Min_Beta_for_absSizeDiff", "model")
    if (!arg$quiet) {
        message(paste0("Running association testing for ", (ncol(calls_file) - no_cols) + 1, " qualifying variants..."))
    }
    for (i in seq(no_cols, ncol(calls_file), 1)) {
        VariantToBeTested <- as.character(colnames(calls_file)[i])
        if (!is.na(covariates)) {
            covlist <- gsub(",", " ", covariates)
            covlist_prepared <- unlist(strsplit(covlist, split = " "))
            selectedtable <- na.omit(as.data.table(cbind(as.character(calls_file_selected[[phenotype]]), as.numeric(calls_file_selected[[VariantToBeTested]]), calls_file_selected[, ..covlist_prepared])))
            colnames(selectedtable)[1:2] <- c(phenotype, VariantToBeTested)
            group2 <- subset(selectedtable, selectedtable[[VariantToBeTested]] >= expandedAllele)
            group1 <- subset(selectedtable, selectedtable[[VariantToBeTested]] < expandedAllele)
            calls_file_selected$expanded_allele_group <- ifelse(calls_file_selected[[VariantToBeTested]] >= expandedAllele, "Expanded", "notExpanded")
            calls_file_selected$expanded_allele_group <- factor(calls_file_selected$expanded_allele_group, c("notExpanded", "Expanded"))
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
            calls_file_selected$expanded_allele_group <- factor(calls_file_selected$expanded_allele_group, c("notExpanded", "Expanded"))
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
        Group2_Group1_Beta_for_absAvgSizeDiff <- round((Group2_Group1_absAvgSizeDiff * Beta), digits = 3)
        Max_Min_Beta_for_absSizeDiff <- round((Max_Min_absSizeDiff * Beta), digits = 3)
        model <- as.character(glm_result$formula)[1]
        tabular_result <- as.data.frame(cbind(Predictors, Beta, Beta_L95, Beta_U95, Beta_stdErr, Pvalue, N, Group1_N, Group2_N, AvgSize, Group1_AvgSize, Group2_AvgSize, Group2_Group1_absAvgSizeDiff, Group2_Group1_Beta_for_absAvgSizeDiff, MinSize, MaxSize, Max_Min_absSizeDiff, Max_Min_Beta_for_absSizeDiff, model))
        tabular_result <- subset(tabular_result, Predictors == VariantID)
        colnames(tabular_result)[1] <- "VariantID"
        tabular_result$VariantID <- paste0(as.character(arg$single_variant), "_ExpandedAllele")
        results_calls_file <- rbind.data.frame(results_calls_file, tabular_result)
        svMisc::progress(i, init = TRUE, progress.bar = FALSE, console = TRUE, gui = FALSE)
        if (i == ncol(calls_file)) {
            message("Done!")
        }
    }
    results_calls_file <- results_calls_file[-1, ]
    sorted_results_calls_file <- results_calls_file[order(as.numeric(results_calls_file$Pvalue)), ]
    write.table(sorted_results_calls_file, arg$out, sep = "\t", col.names = TRUE, quote = F, row.names = F)
}

prepare_phenotype <- function(sample_list, arg) {
    if (!arg$quiet) {
        message("Processing the phenotype file...")
    }
    phenocovar <- fread(arg$phenocovar, header = TRUE)
    # check if phenotype is a column in phenocovar
    if (!arg$phenotype %in% colnames(phenocovar)) {
        stop(paste0("The phenotype variable you provided is not a column in the phenotype file you provided. Please check the column names in your phenotype file."))
    }
    colnames(sample_list) <- "individual"
    return(list(
        "sample_list" = left_join(sample_list, phenocovar, by = "individual"),
        "phenotype" = arg$phenotype <- paste0(arg$phenotype, ""),
        "no_cols" = ncol(phenocovar) + 1
    ))
}

prepare_calls <- function(calls, sample_list_wPheno, arg) {
    if (!arg$quiet) {
        message("Processing the input file based on the STRmode chosen...")
    }
    if (arg$STRmode == "MEAN") {
        calls_file <- transpose((pmax(calls$H1, calls$H2, na.rm = TRUE) + pmin(calls$H1, calls$H2, na.rm = TRUE)) / 2)
    } else if (arg$STRmode == "MAX") {
        calls_file <- transpose(pmax(calls$H1, calls$H2, na.rm = TRUE))
    } else if (arg$STRmode == "MIN") {
        calls_file <- transpose(pmin(calls$H1, calls$H2, na.rm = TRUE))
    }
    if (all(is.na(calls_file))) {
        stop(paste0("The STRmode and run mode you chose resulted in all missing values. Aborting."))
    }

    colnames(calls_file) <- calls$strnames
    calls_file <- cbind(sample_list_wPheno, calls_file)
    calls_file <- calls_file[, which(unlist(lapply(calls_file, function(x) !all(is.na(x))))), with = FALSE]
    calls_file <- data.table(data.frame(calls_file)[, which(colMeans(!is.na(data.frame(calls_file))) >= arg$missing_cutoff)])
    calls_file <- calls_file %>% select(where(~ n_distinct(., na.rm = TRUE) > 1))
    return(calls_file)
}

parse_arguments <- function() {
    p <- argparser::arg_parser("Run association testing for STRs with different modes and options. Version 1.5, November 14, 2022")
    p <- argparser::add_argument(p, "--input", help = "inquiSTR input STR file with a header, first 3 columns are chr, begin, end, and rest are sample ids with inqH1 & inqH2 STR lengths", type = "character", nargs = 1)
    p <- argparser::add_argument(p, "--phenocovar", help = "Phenotype and covariate file with header, first column is individual", type = "character", nargs = 1)
    p <- argparser::add_argument(p, "--covnames", help = "Covariate names you want to use (optional), separated by comma", type = "character", nargs = "*")
    p <- argparser::add_argument(p, "--phenotype", help = "Column name of your phenotype of interest variable in the --phenocovar file", type = "character", nargs = 1)
    p <- argparser::add_argument(p, "--out", help = "Output file name", type = "character", nargs = 1)
    p <- argparser::add_argument(p, "--STRmode", help = "Choose a STRmode from following: MEAN, MAX, MIN; meaning H1+H2 alleles divided by two, maximum of two, or minimum of two. Missing alleles are not considered.", type = "character", nargs = 1)
    p <- argparser::add_argument(p, "--missing_cutoff", help = "Defines the call rate cutoff for variants, default is 0.80 meaning that keeping all variants present in at least 80% of individuals (importantly, both for the input file and for the subset groups you selected to include in association testing). Might mean different things in each of the MEAN, MAX, MIN modes, use carefully.", type = "numeric", default = "0.80")
    p <- argparser::add_argument(p, "--outcometype", help = "Select a outcome variable type: binary or continuous", type = "character", nargs = 1)
    p <- argparser::add_argument(p, "--binaryOrder", help = "Give the binary phenotype order, comma separated, e.g. Control, Patient will code Control as 0/Group1 and Patient as 1/Group2. This will also be used to further filter your data to only two groups (if you had more than >2 groups for your categorical data)", type = "character", nargs = "*")
    p <- argparser::add_argument(p, "--run", help = "Run mode, select among: full (entire file will be analyzed), chromosome (chromosome indicated by --chr will be analyzed), chr_interval (intervals between chr:begin-end defined by --chr --begin --end will be analyzed), bed_interval (an intersection with a bed file defined by --bed will be analyzed), single_variant (a single variant defined by --single_variant will be analyzed; should be combined with --expandedAllele)", type = "character", nargs = 1)
    p <- argparser::add_argument(p, "--chr", help = "Indicate chromosome number to be analyzed (with chr prefix). Optional if bed file is provided.", type = "character", nargs = "?")
    p <- argparser::add_argument(p, "--chr_begin", help = "Define a begin position (inclusive) for a region of interest (optional, and should be combined with --chr_end)", type = "integer", nargs = "?")
    p <- argparser::add_argument(p, "--chr_end", help = "Define a end position (inclusive) for a region of interest (optional, and should be combined with --chr_begin)", type = "integer", nargs = "?")
    p <- argparser::add_argument(p, "--bed", help = "A bed file (without a header) with three columns: chromosome (with chr prefix), begin, and end positions for region(s) of interest (optional). valr Rpackage is required - can be installed with mamba install r-valr on conda environment", type = "character", nargs = "?")
    p <- argparser::add_argument(p, "--single_variant", help = "A single variant to be analyzed, indicated as chr_begin_end OR chr:begin-end, based on a cut-off value provided by the argument --expandedAllele", type = "character", nargs = "?")
    p <- argparser::add_argument(p, "--expandedAllele", help = "A cut-off number (integer or decimal) to define 2 groups with and without expanded allele. STR_length >= expandedAllele is Group 2, STR_length < expandedAllele is Group 1.", type = "integer", nargs = "?")
    p <- argparser::add_argument(p, "--quiet", help = "Do not print progress messages", flag = TRUE)
    arg <- argparser::parse_args(p)
    if (!arg$quiet) {
        Version <- "inquiSTR - STR_regression Rscript Version 1.5.1, July 04, 2023"
        message(Version)
    }
    # argument checks
    if (is.na(arg$input) || is.na(arg$phenocovar) || is.na(arg$phenotype) || is.na(arg$out) || is.na(arg$STRmode) || is.na(arg$outcometype) || is.na(arg$run)) {
        stop("Error: exiting because at least one of the following required arguments is missing: --input, --phenocovar, --phenotype, --out, --STRmode, --outcometype, --run")
    }

    if ((arg$outcometype == "binary") && is.na(arg$binaryOrder)) {
        stop("Error: exiting because --binaryOrder argument is missing, please provide it when you use --outcometype binary")
    }

    if ((arg$run == "chromosome") && is.na(arg$chr)) {
        stop("Error: exiting because --chr argument is missing, please provide it when you use --run chromosome")
    }

    if ((arg$run == "chr_interval") && ((is.na(arg$chr)) || (is.na(arg$chr_begin)) || (is.na(arg$chr_end)))) {
        stop("Error: exiting because At least one of the --chr, --chr_begin, or --chr_end arguments is missing, please provide these when you use --run chr_interval")
    }

    if ((arg$run == "bed_interval") && is.na(arg$bed)) {
        stop("Error: exiting because --bed argument, therefore input bed file, is missing; please provide it when you use --run bed_interval")
    }

    if ((arg$run == "single_variant") && is.na(arg$expandedAllele)) {
        stop("Error: exiting because At least one of the two following aruguments is missing: --single_variant --expandedAllele; please provide these when you use --run single_variant")
    }
    return(arg)
}

arg <- parse_arguments()
# for interactive debugging:
# arg <- list(input = "all.inq.gz", out="test.tsv", chr = "chr22", run="chromosome", phenocovar = "sample_info_with_haplotype.tsv", phenotype = "fus", binaryOrder = "N,Y", covnames="Sex", STRmode = 'MAX', missing_cutoff = 0.8, quiet=FALSE)



# calls is a list with H1, H2 and strnames attributes
if (!arg$quiet) {
    message("Loading and processing the input file...")
}
calls_file <- fread(arg$input, header = TRUE)
if (arg$run == "chr_interval") {
    calls_file <- subset(calls_file, ((chr == arg$chr) & (begin >= arg$chr_begin) & (end <= arg$chr_end)))
    if (!arg$quiet) {
        message("chr_interval run mode is selected")
    }
} else if (arg$run == "bed_interval") {
    bedfile <- fread(arg$bed, header = FALSE)
    colnames(bedfile) <- c("chrom", "start", "end")
    colnames(calls_file)[2] <- "start"
    calls_file <- as.data.table(bed_intersect(calls_file, bedfile, suffix = c("", ".y")))
    calls_file <- intersecttable[, 1:(length(intersecttable) - 3)]
    calls_file <- subset(intersecttable, !is.na(chrom))
    colnames(calls_file)[2] <- "begin"
    if (!arg$quiet) {
        message("bed_interval run mode is selected")
    }
} else if (arg$run == "chromosome") {
    calls_file <- subset(calls_file, chr == arg$chr)
    if (!arg$quiet) {
        message("chromosome run mode is selected")
    }
} else if (arg$run == "full") {
    if (!arg$quiet) {
        message("full run mode is selected")
    }
} else if (arg$run == "single_variant") {
    single_variant_toAnalyze <- unlist(strsplit(arg$single_variant, split = "_"))
    single_variant_toAnalyze <- unlist(strsplit(unlist(strsplit(single_variant_toAnalyze, split = "-")), ":"))
    calls_file <- subset(calls_file, ((chrom == single_variant_toAnalyze[1]) & (begin == single_variant_toAnalyze[2]) & (end == single_variant_toAnalyze[3])))
    if (!arg$quiet) {
        message("single_variant run mode is selected (with expandedAllele)")
    }
}
strnames <- paste(paste(calls_file$chr, calls_file$begin, sep = ":"), calls_file$end, sep = "_")
rest <- calls_file[, -c(1:3)]
col_index <- seq_len(ncol(rest))
inqH1 <- as.data.table(rest %>% select(col_index[col_index %% 2 != 0]))
inqH2 <- as.data.table(rest %>% select(col_index[col_index %% 2 == 0]))
colnames(inqH1) <- gsub("_H1$", "", colnames(inqH1))
colnames(inqH2) <- gsub("_H2$", "", colnames(inqH2))
calls <- list("H1" = inqH1, "H2" = inqH2, "strnames" = strnames)

# pheno_info is a list with sample_list, phenotype and no_cols attributes
pheno_info <- prepare_phenotype(
    sample_list = as.data.table(colnames(calls$H1)),
    arg = arg
)


calls_file <- prepare_calls(
    calls = calls,
    sample_list_wPheno = pheno_info$sample_list,
    arg = arg
)

if (arg$outcometype == "binary" && arg$run != "single_variant") {
    assoc_binary(
        arg = arg,
        calls_file = calls_file,
        phenotype = pheno_info$phenotype,
        no_cols = pheno_info$no_cols,
        covariates = arg$covnames,
        missing_cutoff = arg$missing_cutoff
    )
} else if (arg$outcometype == "continuous" && arg$run != "single_variant") {
    assoc_continuous(
        arg = arg,
        calls_file = calls_file,
        phenotype = pheno_info$phenotype,
        no_cols = pheno_info$no_cols,
        covariates = arg$covnames,
        missing_cutoff = arg$missing_cutoff
    )
} else if (arg$outcometype == "binary" && arg$run == "single_variant") {
    assoc_binary_expandedAllele(
        arg = arg,
        calls_file = calls_file,
        phenotype = pheno_info$phenotype,
        no_cols = pheno_info$no_cols,
        covariates = arg$covnames,
        expandedAllele = arg$expandedAllele,
        missing_cutoff = arg$missing_cutoff
    )
} else if (arg$outcometype == "continuous" && arg$run == "single_variant") {
    assoc_continuous_expandedAllele(
        arg = arg,
        calls_file = calls_file,
        phenotype = pheno_info$phenotype,
        no_cols = pheno_info$no_cols,
        covariates = arg$covnames,
        expandedAllele = arg$expandedAllele,
        missing_cutoff = arg$missing_cutoff
    )
}
