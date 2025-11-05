data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects", "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "11-postmortem-data-analysis/")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
temp_dir <- file.path(analysis_dir, "tmp")
alsod_genes_df <- read.csv(file.path(project_dir, "data_files/alsod_genelist.csv"))
annotation_gtf <- "/nemo/lab/patanir/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
tmp_dir <- file.path(analysis_dir, "slurm")
genic_count_list <- readRDS(file.path(analysis_save_dir, "01-gene-level-counts-of-postmortem-samples.rds"))
# ==========================================================================
# Imports
# ==========================================================================
library("stringr")
library("GenomicRanges")
library("DESeq2")
library("rtracklayer")
library("ggplot2")
library("GenomicFeatures")
library("GenomeInfoDb")
require("Rsamtools")
require("STRINGdb")
require("survival")
source("scripts/src/functions.R")
# source("../src/rmats-analysis/src/source.R")
# ==========================================================================
# Analysis
# ==========================================================================
annotation_ranges <- import.gff2(annotation_gtf)
annotation_gene_ranges <- annotation_ranges[annotation_ranges$type == "gene"]

# feature_ranges_with_repeats <- significant_de_event_list$repeat_containing_ranges
counts_matrix <- genic_count_list$counts 
tpm_matrix <- counts_matrix / (width(annotation_gene_ranges[match(rownames(counts_matrix), annotation_gene_ranges$gene_id)]) / 1e3)
tpm_matrix <- t(t(tpm_matrix) / (colSums(tpm_matrix) / 1e6))

coldata_df <- genic_count_list$coldata
coldata_df$expr <- tpm_matrix[annotation_gene_ranges$gene_id[annotation_gene_ranges$gene_name == "ERAP2" & !is.na(annotation_gene_ranges$gene_name)],
    coldata_df$column_id]

coldata_split_by_tissue <- split(coldata_df, coldata_df$tissue)
patients_expr_zscore_df <- do.call(rbind, lapply(coldata_split_by_tissue, function(x){
    erap2_zscore <- (x$expr - mean(x$expr)) / sd(x$expr)
    x$erap2_zscore <- erap2_zscore
    return(x)
}))
rownames(patients_expr_zscore_df) <- NULL
patients_expr_zscore_df$erap2_group <- ifelse(patients_expr_zscore_df$erap2_zscore > 0, "high", "low")

patients_expr_zscore_df$patients <- paste(patients_expr_zscore_df$subject_id, patients_expr_zscore_df$tissue)


summary_df <- summarise_tissue_specific_parameter_by_patients(coldata_df)
summary_df$spinal_cord_lumbar_group <- patients_expr_zscore_df$erap2_group[match(paste(summary_df$subject_id, "Spinal_Cord_Lumbar"), patients_expr_zscore_df$patients)]
summary_df$spinal_cord_cervical_group <- patients_expr_zscore_df$erap2_group[match(paste(summary_df$subject_id, "Spinal_Cord_Cervical"), patients_expr_zscore_df$patients)]
summary_df$spinal_cord_thoracic_group <- patients_expr_zscore_df$erap2_group[match(paste(summary_df$subject_id, "Spinal_Cord_Thoracic"), patients_expr_zscore_df$patients)]
summary_df$cortex_group <- patients_expr_zscore_df$erap2_group[match(paste(summary_df$subject_id, "Cortex"), patients_expr_zscore_df$patients)]
summary_df$cortex_motor_group <- patients_expr_zscore_df$erap2_group[match(paste(summary_df$subject_id, "Cortex_Motor"), patients_expr_zscore_df$patients)]
summary_df$Cerebellum_group <- patients_expr_zscore_df$erap2_group[match(paste(summary_df$subject_id, "Cerebellum"), patients_expr_zscore_df$patients)]





perform_survival_df <- summary_df[summary_df$als_group %in% c("als", "als+ftd"),]
perform_survival_df$disease_duration[perform_survival_df$disease_duration == "~84"] <- "84"
# perform_survival_df <- perform_survival_df[!(perform_survival_df$disease_duration %in% c("Unknown", "Not Applicable")),]
perform_survival_df$disease_duration <- as.numeric(perform_survival_df$disease_duration)
perform_survival_df <- perform_survival_df[perform_survival_df$disease_duration <= 48,]
perform_survival_df$event <- ifelse(!is.na(perform_survival_df$disease_duration), 1, 0)
# perform_survival_df$event <- ifelse(perform_survival_df$als_group %in% c("als", "als+ftd"), 1, 0)


# no_ctrl_tissue_specific_df <- tissue_specific_df[tissue_specific_df$als_group %in% c("als", "als+ftd"),]
# no_ctrl_tissue_specific_df <- no_ctrl_tissue_specific_df[no_ctrl_tissue_specific_df$disease_duration > 12 & no_ctrl_tissue_specific_df$disease_duration <= 36,]
# no_ctrl_tissue_specific_df$event <- ifelse(!is.na(no_ctrl_tissue_specific_df$disease_duration), 1, 0)



require("survival")
require("survminer")
require("dplyr")
require("tidycmprsk")

tissues <- c(lumbar = "spinal_cord_lumbar", cervical = "spinal_cord_cervical", 
    thoracic = "spinal_cord_thoracic", cerebellum = "Cerebellum", 
    cortex = "cortex", motor_cortex = "cortex_motor")

survival_list <- lapply(names(tissues), function(x){
    a_tissue <- tissues[x]
    a_col <- paste(a_tissue, "group", sep = "_")

    survival_df <- perform_survival_df
    survival_df$erap2_group <- survival_df[,a_col]
    survival_df[,"erap2_group"] <- factor(survival_df[,"erap2_group"], levels = c("low", "high"))
    fit1 <- survfit(Surv(time = disease_duration, event = event, type = "right") ~ erap2_group, data = survival_df)
    file_path <- file.path(analysis_plots_dir, 
        paste("08-als-survival-analysis-of", x, "patient-survival-duration-vs-higher-erap2-expression.pdf", sep = "-"))
    the_plot <- ggsurvplot(fit = fit1, data = survival_df, pval = TRUE, color = "strata", palette = c("#fb9a99", "#e31a1c"))
    the_plot <- the_plot$plot + 
    labs(title = paste("Survival prob in", x, "versus", "onset", sep = " "), x = "Disease duration", y = "Survival probability")
    ggsave(file = file_path, plot = the_plot)
})
names(survival_list) <- tissues








perform_survival_df <- summary_df
perform_survival_df$age_at_onset <- as.numeric(perform_survival_df$age_at_onset)
perform_survival_df$event <- ifelse(perform_survival_df$als_group %in% c("als", "als+ftd"), 1, 0)

survival_list <- lapply(names(tissues), function(x){
    a_tissue <- tissues[x]
    a_col <- paste(a_tissue, "group", sep = "_")

    survival_df <- perform_survival_df
    survival_df$erap2_group <- survival_df[,a_col]
    survival_df[,"erap2_group"] <- factor(survival_df[,"erap2_group"], levels = c("low", "high"))
    fit1 <- survfit(Surv(time = age_at_death, event = event, type = "right") ~ erap2_group, data = survival_df)
    file_path <- file.path(analysis_plots_dir, 
        paste("08-als-survival-analysis-of", x, "patient-survival-age-at-death-vs-higher-erap2-expression.pdf", sep = "-"))
    the_plot <- ggsurvplot(fit = fit1, data = survival_df, pval = TRUE, color = "strata", palette = c("#fb9a99", "#e31a1c"))
    the_plot <- the_plot$plot + 
    labs(title = paste("Survival prob in", x, "versus", "duration", sep = " "), x = "Age at death", y = "Survival probability")
    ggsave(file = file_path, plot = the_plot)
})
names(survival_list) <- tissues



perform_survival_df <- summary_df
perform_survival_df$age_at_onset <- as.numeric(perform_survival_df$age_at_onset)
perform_survival_df$event <- ifelse(perform_survival_df$als_group %in% c("als", "als+ftd"), 1, 0)

survival_list <- lapply(names(tissues), function(x){
    a_tissue <- tissues[x]
    a_col <- paste(a_tissue, "group", sep = "_")

    survival_df <- perform_survival_df
    survival_df$erap2_group <- survival_df[,a_col]
    survival_df[,"erap2_group"] <- factor(survival_df[,"erap2_group"], levels = c("low", "high"))
    fit1 <- survfit(Surv(time = age_at_onset, event = event, type = "right") ~ erap2_group, data = survival_df)
    file_path <- file.path(analysis_plots_dir, 
        paste("09-als-survival-analysis-of", x, "patient-survival-age-at-onset-vs-higher-erap2-expression.pdf", sep = "-"))
    the_plot <- ggsurvplot(fit = fit1, data = survival_df, pval = TRUE, color = "strata", palette = c("#fb9a99", "#e31a1c"))
    the_plot <- the_plot$plot + 
    labs(title = paste("Survival prob in", x, "versus", "onset", sep = " "), x = "Age at onset", y = "Survival probability")
    ggsave(file = file_path, plot = the_plot)
})
names(survival_list) <- tissues





