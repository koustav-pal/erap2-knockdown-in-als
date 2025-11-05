data_dir <- normalizePath("/export2/home/koustav/")
project_dir <- file.path(data_dir, "Projects", "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "18-screening-for-erap2-mutations-in-postmortem")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_output_dir <- file.path(analysis_dir, "output_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
temp_dir <- file.path(analysis_dir, "tmp")
tmp_dir <- file.path(analysis_dir, "slurm")
patient_count_df <- readRDS(file.path("../11-postmortem-data-analysis/rdata_files/02-summarised-counts-for-patients.rds"))
mutation_screening_results <- file.path(analysis_output_dir, "FILTERED_MUTATION_SCREENING_RESULTS.vcf")
# ==========================================================================
# Packages
# ==========================================================================
setwd("/export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis")
require("HiCBricks")
require("GenomicRanges")
require("stringr")
# require("Biostrings")
# require("rtracklayer")
# require("GenomeInfoDb")
# require("GenomicFeatures")
require("ggplot2")
require("ggbeeswarm") 
library(rlang)
library(survival)
library(reshape2)
library(knitr)
library(rmarkdown)
require("RColorBrewer")
require("survminer")
require("survival")
require("randomForestSRC")
require("parallel")
require("matrixStats")
require("patchwork")
# ==========================================================================
# Analysis: Prepare the survival data.frame
# ==========================================================================
patient_count_df$type <- factor(patient_count_df$type, levels = c("control", "c9orf72-ALS", "other_mutation-fALS", "sporadic"))
patient_count_df$type_group <- as.character(patient_count_df$type)
patient_count_df$type_group[patient_count_df$type %in% c("other_mutation-fALS", "sporadic")] <- "sporadic+other_fALS"
patient_count_df$tissue_group <- "spinal_cord"
patient_count_df$tissue_group[patient_count_df$tissue %in% c("Cortex", "Cortex_Motor", "Cerebellum")] <- "brain"

patient_surv_df <- do.call(rbind, lapply(split(patient_count_df, paste(patient_count_df$subject_id, patient_count_df$tissue)), function(temp_df){
    # message(unique(temp_df$subject_id))
    data.frame(subject_id = unique(temp_df$subject_id),
        # tissue_group = unique(temp_df$tissue_group),
        tissue = unique(temp_df$tissue),
        age_at_onset = unique(temp_df$age_at_symptom_onset),
        age_at_death = unique(temp_df$age_at_death),
        disease_duration_in_months = unique(temp_df$disease_duration_in_months),
        sex_genotype = unique(temp_df$sex_genotype),
        mnd_with_ftd = unique(temp_df$mnd_with_ftd),
        condition = unique(temp_df$condition),
        mnd_with_dementia = unique(temp_df$mnd_with_dementia),
        type_group = unique(temp_df$type_group),
        erap2_expression = mean(temp_df$expr))
}))

patient_surv_df$age_at_death[patient_surv_df$age_at_death == "90 or Older"] <- "90"
patient_surv_df$age_at_death <- as.numeric(patient_surv_df$age_at_death)
patient_surv_df$age_at_onset <- as.numeric(patient_surv_df$age_at_onset)
patient_surv_df$disease_duration <- as.numeric(patient_surv_df$disease_duration_in_months)

# ==========================================================================
# Analysis: Select the spinal cord samples for ALS patients
# ==========================================================================

als_surv_df <- patient_surv_df[patient_surv_df$condition %in% c("als"),]
als_surv_df$sex_genotype <- factor(als_surv_df$sex_genotype, levels = c("XX", "XY"))

# als_surv_df$condition <- factor(als_surv_df$condition, levels = c("ctrl", "als"))
als_surv_df$status <- 1
als_surv_df$status[is.na(als_surv_df$age_at_death)] <- 0

spinal_cord_df <- als_surv_df[als_surv_df$tissue %in% c("Spinal_Cord_Cervical", "Spinal_Cord_Lumbar", "Spinal_Cord_Thoracic"),c("subject_id", "tissue", "age_at_onset", "age_at_death", "disease_duration", "status", "condition", "sex_genotype", "erap2_expression")]


# ==========================================================================
# Add information about mutations
# ==========================================================================
mutation_list <- readRDS("18-screening-for-erap2-mutations-in-postmortem/rdata_files/01-list-of-selected-mutations-for-each-subject.rds")

patient_specific_mutation_list <- mutation_list$patient_mutation_list
frequency_matrix <- mutation_list$frequency_matrix
selected_mutation_df <- mutation_list$mutation_list


#
# Select mutations that are present in at least 10% of ALS cases and 
# less that 5% of Control groups.
# --------------------------------------------------------------------------
candidate_mutations_df <- selected_mutation_df[selected_mutation_df$als > 0.1 & selected_mutation_df$ctrl < 0.05,]

keep_patient_specific_mutation_list <- patient_specific_mutation_list[names(patient_specific_mutation_list) %in% spinal_cord_df$subject_id]

patient_mutation_frequency_df <- t(do.call(cbind,lapply(keep_patient_specific_mutation_list, function(x){ x$occurence_vector[candidate_mutations_df$snp_id] })))
colnames(patient_mutation_frequency_df) <- paste("rs", candidate_mutations_df$snp_id, sep = "")

spinal_cord_with_mutation_df <- cbind(spinal_cord_df, patient_mutation_frequency_df[spinal_cord_df$subject_id,])

# ==========================================================================
# Prepare the random survival forests
# ==========================================================================
# Recycle the disease duration seeds
seed_file <- file.path(project_dir, "analysis", "11-postmortem-data-analysis/rdata_files/01-seed-file-for-comparing-disease-duration-RFS-survival-curves-between-als-and-control.rds")

if(!file.exists(seed_file)){
    seed_list <- sample.int(.Machine$integer.max, 1000)    
    saveRDS(object = seed_list, file = seed_file)
}else{
    seed_list <- readRDS(seed_file)
}

text_formula <- formula(paste("Surv(disease_duration, status)", paste("age_at_onset", "sex_genotype", "erap2_expression", paste(colnames(patient_mutation_frequency_df), collapse = " + "), sep = " + "), sep = " ~ "))

#
# RSF models for cervical samples
# --------------------------------------------------------------------------
cervical_df <- spinal_cord_with_mutation_df[spinal_cord_with_mutation_df$tissue == "Spinal_Cord_Cervical",]
message("Preparing Cervical RSF models...")
cervical_models <- lapply(seq_along(seed_list), function(x){
    set.seed(seed_list[x])
    rsf_model <- rfsrc(text_formula, 
        data = cervical_df, ntree = 1000, importance = TRUE)
})

#
# RSF models for lumbar samples
# --------------------------------------------------------------------------
lumbar_df <- spinal_cord_with_mutation_df[spinal_cord_with_mutation_df$tissue == "Spinal_Cord_Lumbar",]
message("Preparing Lumbar RSF models...")
lumbar_models <- lapply(seq_along(seed_list), function(x){
    set.seed(seed_list[x])
    rsf_model <- rfsrc(text_formula, 
        data = lumbar_df, ntree = 1000, importance = TRUE)
})

#
# RSF models for lumbar samples
# --------------------------------------------------------------------------
thoracic_df <- spinal_cord_with_mutation_df[spinal_cord_with_mutation_df$tissue == "Spinal_Cord_Thoracic",]
message("Preparing Thoracic RSF models...")
thoracic_models <- lapply(seq_along(seed_list), function(x){
    set.seed(seed_list[x])
    rsf_model <- rfsrc(Surv(disease_duration, status) ~ age_at_onset + sex_genotype + erap2_expression, 
        data = thoracic_df, ntree = 1000, importance = TRUE)
})


# ==========================================================================
# Fetch importance of all features
# ==========================================================================
a_list <- list("cervical" = cervical_models, "lumbar" = lumbar_models, "thoracic" = thoracic_models)
all_importance_df <- do.call(rbind, lapply(names(a_list), function(x){
    importance_df <- do.call(rbind, lapply(a_list[[x]], function(y){
        importance_vector <- y$importance[c("age_at_onset", "sex_genotype", 
        	colnames(patient_mutation_frequency_df), "erap2_expression")]
        data.frame(tissue = x,
            variable = names(importance_vector),
            values = as.vector(importance_vector))
    }))
}))
all_importance_df$tissue <- factor(all_importance_df$tissue, levels = c("cervical", "lumbar", "thoracic"))

median_values_df <- do.call(rbind, lapply(split(all_importance_df, all_importance_df$variable), function(x){
    return_df <- data.frame(variable = unique(x[,"variable"]), median = NA)
    return_df$median <- median(x$values)
    return(return_df)
}))

# ==========================================================================
# Plot importance violin plots
# ==========================================================================

current_theme <- theme(
    legend.position = "bottom",
    text = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "azure3",
    linetype = "dashed",
    size = 0.2))

all_importance_df$variable <- factor(all_importance_df$variable, 
	levels = median_values_df$variable[order(median_values_df$median, decreasing = T)])

colours <- c("high" = "#FF5666", "low" = "#FFCCC9")
The_plot <- ggplot(all_importance_df, aes(x = variable, y = values))
The_plot <- The_plot + geom_violin(fill = "#F1FAEE", colour = "#000000")
The_plot <- The_plot + geom_hline(yintercept = 0.01)
The_plot <- The_plot + geom_boxplot(fill = "#F1FAEE", colour = "#000000", 
    width = 0.1, lwd = 0.5, outlier.shape = NA)
The_plot <- The_plot + scale_x_discrete("model variables")
The_plot <- The_plot + scale_y_continuous("importance")
# The_plot <- The_plot + scale_fill_manual("Risk group", values = colours)
The_plot <- The_plot + labs(title = "Importance of model variables")
The_plot <- The_plot + facet_wrap(tissue ~ ., scales = "free")
The_plot <- The_plot + theme_bw() + current_theme
ggsave(width = 10, height = 15, file = file.path(analysis_plots_dir, "01-01-plot-importance-of-model-variables.pdf"), The_plot)


# ==========================================================================
# Fetch risk scores for cervical samples
# ==========================================================================
fetch_columns <- c("subject_id", "tissue", "disease_duration", 
	"age_at_onset", "sex_genotype", "status", "condition", 
	"erap2_expression", colnames(patient_mutation_frequency_df))
non_na_cervical_df <- cervical_df[
    !is.na(cervical_df$disease_duration) & 
    !is.na(cervical_df$age_at_onset), fetch_columns]

# get predicted risk 
message("Computing Cervical risk scores...")

cervical_risk_score_matrix <- do.call(cbind, lapply(cervical_models, function(x){
    rsf_pred <- predict(x, newdata = non_na_cervical_df)
    rsf_pred$predicted
}))

non_na_cervical_df$risk_score <- rowMedians(cervical_risk_score_matrix)
median_risk <- median(non_na_cervical_df$risk_score, na.rm = T)
non_na_cervical_df$risk_group <- ifelse(non_na_cervical_df$risk_score >= median_risk, "high", "low")


differential_survival <- survdiff(Surv(disease_duration, status) ~ risk_group, data = non_na_cervical_df)
hazard_ratio <- (differential_survival$obs[1] / differential_survival$exp[1]) / (differential_survival$obs[2] / differential_survival$exp[2])

message("Cervical disease duration HR:", hazard_ratio)

wilcox_test <- wilcox.test(non_na_cervical_df[non_na_cervical_df$risk_group == "high", "erap2_expression"],
    non_na_cervical_df[non_na_cervical_df$risk_group == "low", "erap2_expression"])

message("Wilcoxon test p-value for Cervical ERAP2 stratification:", wilcox_test$p.value)


# ==========================================================================
# Fetch risk scores for lumbar samples
# ==========================================================================
non_na_lumbar_df <- lumbar_df[
    !is.na(lumbar_df$disease_duration) & !is.na(lumbar_df$age_at_onset), fetch_columns]

message("Computing Lumbar risk scores...")
lumbar_risk_score_matrix <- do.call(cbind, lapply(lumbar_models, function(x){
    rsf_pred <- predict(x, newdata = non_na_lumbar_df)
    rsf_pred$predicted
}))

non_na_lumbar_df$risk_score <- rowMedians(lumbar_risk_score_matrix)
median_risk <- median(non_na_lumbar_df$risk_score, na.rm = T)
non_na_lumbar_df$risk_group <- ifelse(non_na_lumbar_df$risk_score >= median_risk, "high", "low")

differential_survival <- survdiff(Surv(disease_duration, status) ~ risk_group, data = non_na_lumbar_df)
hazard_ratio <- (differential_survival$obs[1] / differential_survival$exp[1]) / (differential_survival$obs[2] / differential_survival$exp[2])
message("Lumbar disease duration HR:", hazard_ratio)

wilcox_test <- wilcox.test(non_na_lumbar_df[non_na_lumbar_df$risk_group == "high", "erap2_expression"],
    non_na_lumbar_df[non_na_lumbar_df$risk_group == "low", "erap2_expression"])

message("Wilcoxon test p-value for Lumbar ERAP2 stratification:", wilcox_test$p.value)

# ==========================================================================
# Fetch risk scores for thoracic samples
# ==========================================================================
non_na_thoracic_df <- thoracic_df[
    !is.na(thoracic_df$disease_duration) & !is.na(thoracic_df$age_at_onset), fetch_columns]

message("Computing Thoracic risk scores...")
thoracic_risk_score_matrix <- do.call(cbind, lapply(thoracic_models, function(x){
    rsf_pred <- predict(x, newdata = non_na_thoracic_df)
    rsf_pred$predicted
}))

non_na_thoracic_df$risk_score <- rowMedians(thoracic_risk_score_matrix)
median_risk <- median(non_na_thoracic_df$risk_score, na.rm = T)
non_na_thoracic_df$risk_group <- ifelse(non_na_thoracic_df$risk_score >= median_risk, "high", "low")

differential_survival <- survdiff(Surv(disease_duration, status) ~ risk_group, data = non_na_thoracic_df)
hazard_ratio <- (differential_survival$obs[1] / differential_survival$exp[1]) / (differential_survival$obs[2] / differential_survival$exp[2])
message("Thoracic disease duration HR:", hazard_ratio)

wilcox_test <- wilcox.test(non_na_thoracic_df[non_na_thoracic_df$risk_group == "high", "erap2_expression"],
    non_na_thoracic_df[non_na_thoracic_df$risk_group == "low", "erap2_expression"])

message("Wilcoxon test p-value for Thoracic ERAP2 stratification:", wilcox_test$p.value)

non_na_cervical_df$tissue <- "cervical"
non_na_lumbar_df$tissue <- "lumbar"
non_na_thoracic_df$tissue <- "thoracic"
all_patient_risk_stratified_df <- rbind(non_na_cervical_df, non_na_lumbar_df, non_na_thoracic_df)
saveRDS(object = all_patient_risk_stratified_df, file = file.path(analysis_save_dir, "02-all-patient-disease-duration-risk-stratified-sampleinfo.rds"))

# ==========================================================================
# Plot cervical survival curves for disease duration
# ==========================================================================

#
# Survival plots
# --------------------------------------------------------------------------
cervical_fit <- survfit(Surv(disease_duration, status) ~ risk_group, data = non_na_cervical_df)

the_plot <- ggsurvplot(cervical_fit, 
    data = non_na_cervical_df,
    palette = c("#FF5666", "#FFCCC9"),
    pval = T, pval.method = T,
    conf.int = F,
    # surv.median.line = c("v"),
    risk.table = "absolute",
    ggtheme = theme_bw() + current_theme)$plot
cervical_plot <- the_plot + ggtitle(paste("RSF Cervical"))

#
# Expression plots for cervical region
# --------------------------------------------------------------------------
colours <- c("high" = "#FF5666", "low" = "#FFCCC9")
cervical_expr_plot <- ggplot(non_na_cervical_df, aes(x = risk_group, y = erap2_expression))
cervical_expr_plot <- cervical_expr_plot + geom_boxplot(aes(fill = risk_group), 
	colour = "#000000", width = 0.1, lwd = 0.5, outlier.shape = NA)
cervical_expr_plot <- cervical_expr_plot + geom_quasirandom(aes(fill = risk_group), shape = 21, colour = "#000000")
cervical_expr_plot <- cervical_expr_plot + scale_x_discrete("Risk group")
cervical_expr_plot <- cervical_expr_plot + scale_y_continuous("ERAP2 expression")
cervical_expr_plot <- cervical_expr_plot + scale_fill_manual("Risk group", values = colours)
cervical_expr_plot <- cervical_expr_plot + labs(title = "Cervical ERAP2 expr")
cervical_expr_plot <- cervical_expr_plot + theme_bw() + current_theme

total_cervical_plot <- cervical_plot | cervical_expr_plot

# ==========================================================================
# Plot lumbar survival curves for disease duration
# ==========================================================================
lumbar_fit <- survfit(Surv(disease_duration, status) ~ risk_group, data = non_na_lumbar_df)

the_plot <- ggsurvplot(lumbar_fit, 
    data = non_na_lumbar_df,
    palette = c("#FF5666", "#FFCCC9"),
    pval = T, pval.method = T,
    conf.int = F,
    # surv.median.line = c("v"),
    risk.table = "absolute",
    ggtheme = theme_bw() + current_theme)$plot
lumbar_plot <- the_plot + ggtitle(paste("RSF Lumbar"))

#
# Expression plots for lumbar region
# --------------------------------------------------------------------------
colours <- c("high" = "#FF5666", "low" = "#FFCCC9")
lumbar_expr_plot <- ggplot(non_na_lumbar_df, aes(x = risk_group, y = erap2_expression))
lumbar_expr_plot <- lumbar_expr_plot + geom_boxplot(aes(fill = risk_group), colour = "#000000", width = 0.1, lwd = 0.5, outlier.shape = NA)
lumbar_expr_plot <- lumbar_expr_plot + geom_quasirandom(aes(fill = risk_group), shape = 21, colour = "#000000")
lumbar_expr_plot <- lumbar_expr_plot + scale_x_discrete("Risk group")
lumbar_expr_plot <- lumbar_expr_plot + scale_y_continuous("ERAP2 expression")
lumbar_expr_plot <- lumbar_expr_plot + scale_fill_manual("Risk group", values = colours)
lumbar_expr_plot <- lumbar_expr_plot + labs(title = "Lumbar ERAP2 expr")
lumbar_expr_plot <- lumbar_expr_plot + theme_bw() + current_theme

total_lumbar_plot <- lumbar_plot | lumbar_expr_plot

# ==========================================================================
# Plot thoracic survival curves for disease duration
# ==========================================================================
thoracic_fit <- survfit(Surv(disease_duration, status) ~ risk_group, data = non_na_thoracic_df)

the_plot <- ggsurvplot(thoracic_fit, 
    data = non_na_thoracic_df,
    palette = c("#FF5666", "#FFCCC9"),
    pval = T, pval.method = T,
    conf.int = F,
    # surv.median.line = c("v"),
    risk.table = "absolute",
    ggtheme = theme_bw() + current_theme)$plot
thoracic_plot <- the_plot + ggtitle(paste("RSF Thoracic"))

#
# Expression plots for thoracic region
# --------------------------------------------------------------------------
colours <- c("high" = "#FF5666", "low" = "#FFCCC9")
thoracic_expr_plot <- ggplot(non_na_thoracic_df, aes(x = risk_group, y = erap2_expression))
thoracic_expr_plot <- thoracic_expr_plot + geom_boxplot(aes(fill = risk_group), colour = "#000000", width = 0.1, lwd = 0.5, outlier.shape = NA)
thoracic_expr_plot <- thoracic_expr_plot + geom_quasirandom(aes(fill = risk_group), shape = 21, colour = "#000000")
thoracic_expr_plot <- thoracic_expr_plot + scale_x_discrete("Risk group")
thoracic_expr_plot <- thoracic_expr_plot + scale_y_continuous("ERAP2 expression")
thoracic_expr_plot <- thoracic_expr_plot + scale_fill_manual("Risk group", values = colours)
thoracic_expr_plot <- thoracic_expr_plot + labs(title = "Thoracic ERAP2 expr")
thoracic_expr_plot <- thoracic_expr_plot + theme_bw() + current_theme

total_thoracic_plot <- thoracic_plot | thoracic_expr_plot

all_plots <- total_cervical_plot / total_lumbar_plot / total_thoracic_plot


ggsave(width = 10, height = 15, 
	file = file.path(analysis_plots_dir, "01-03-plot-survival-curves-for-patient-prognosis-using-disease-duration-and-erap2-expression.pdf"), 
	all_plots)

