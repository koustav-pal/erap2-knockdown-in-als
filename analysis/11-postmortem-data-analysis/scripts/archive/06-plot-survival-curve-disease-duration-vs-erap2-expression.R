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
erap2_kd_de_df <- readRDS(file.path(analysis_save_dir, "01-all-contrasts-differential-splicing-dataframe.rds"))
temp_dir <- file.path(analysis_dir, "tmp")
alsod_genes_df <- read.csv(file.path(project_dir, "data_files/alsod_genelist.csv"))
annotation_gtf <- "/nemo/lab/patanir/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
rmats_data_list <- readRDS(file.path(analysis_save_dir, "01-imported-rmats-object-list.rds"))
tmp_dir <- file.path(analysis_dir, "slurm")
genic_count_list <- readRDS(file.path(analysis_save_dir, "01-gene-level-counts-of-postmortem-samples.rds"))
umap_clustering_list <- readRDS(file.path(analysis_save_dir, "01-patient-latent-factor-umap-clustering.rds"))
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

clustering_df <- umap_clustering_list$umap
patient_matrix <- umap_clustering_list$counts

clustering_df$subject_id <- str_split(clustering_df$patients, pattern = " ", simplify = T)[,1]
clustering_df$sex_genotype <- coldata_df$sex_genotype[match(clustering_df$subject_id, coldata_df$subject_id)]
clustering_df$age_at_symptom_onset <- coldata_df$age_at_symptom_onset[match(clustering_df$subject_id, coldata_df$subject_id)]
clustering_df$age_at_death <- coldata_df$age_at_death[match(clustering_df$subject_id, coldata_df$subject_id)]
clustering_df$disease_duration_in_months <- coldata_df$disease_duration_in_months[match(clustering_df$subject_id, coldata_df$subject_id)]
clustering_df$group <- NA
clustering_df$group[clustering_df$cluster %in% c(1)] <- "cortex"
clustering_df$group[clustering_df$tissue %in% c("Spinal_Cord_Lumbar", "Spinal_Cord_Cervical", "Spinal_Cord_Thoracic")] <- "spinal_cord"
clustering_df$group[clustering_df$tissue %in% c("Cerebellum")] <- "cerebellum"

cols_split_by_group <- split(clustering_df$patients, clustering_df$tissue)
patients_expr_zscore_matrix <- do.call(cbind, lapply(cols_split_by_group, function(x){
    current_matrix <- patient_matrix[,x]
    zscore_matrix <- (current_matrix - rowMeans(current_matrix)) / rowSds(current_matrix)
    return(zscore_matrix)
}))
clustering_df$erap2_zscore <- patients_expr_zscore_matrix["ENSG00000164308",clustering_df$patients]
clustering_df$erap2_group <- ifelse(clustering_df$erap2_zscore > 0, "high", "low")

summary_df <- summarise_tissue_specific_parameter_by_patients(coldata_df)
summary_df$spinal_cord_lumbar_group <- clustering_df$erap2_group[match(paste(summary_df$subject_id, "Spinal_Cord_Lumbar"), clustering_df$patients)]
summary_df$spinal_cord_cervical_group <- clustering_df$erap2_group[match(paste(summary_df$subject_id, "Spinal_Cord_Cervical"), clustering_df$patients)]
summary_df$spinal_cord_thoracic_group <- clustering_df$erap2_group[match(paste(summary_df$subject_id, "Spinal_Cord_Thoracic"), clustering_df$patients)]
summary_df$cortex_group <- clustering_df$erap2_group[match(paste(summary_df$subject_id, "Cortex"), clustering_df$patients)]
summary_df$cortex_motor_group <- clustering_df$erap2_group[match(paste(summary_df$subject_id, "Cortex_Motor"), clustering_df$patients)]
summary_df$Cerebellum_group <- clustering_df$erap2_group[match(paste(summary_df$subject_id, "Cerebellum"), clustering_df$patients)]





perform_survival_df <- summary_df[summary_df$als_group %in% c("als", "als+ftd"),]
perform_survival_df$disease_duration[perform_survival_df$disease_duration == "~84"] <- "84"
# perform_survival_df <- perform_survival_df[!(perform_survival_df$disease_duration %in% c("Unknown", "Not Applicable")),]
perform_survival_df$disease_duration <- as.numeric(perform_survival_df$disease_duration)
# perform_survival_df <- perform_survival_df[perform_survival_df$disease_duration > 12 & perform_survival_df$disease_duration <= 36,]
perform_survival_df$event <- ifelse(!is.na(perform_survival_df$disease_duration), 1, 0)
# perform_survival_df$event <- ifelse(perform_survival_df$als_group %in% c("als", "als+ftd"), 1, 0)


# no_ctrl_tissue_specific_df <- tissue_specific_df[tissue_specific_df$als_group %in% c("als", "als+ftd"),]
# no_ctrl_tissue_specific_df <- no_ctrl_tissue_specific_df[no_ctrl_tissue_specific_df$disease_duration > 12 & no_ctrl_tissue_specific_df$disease_duration <= 36,]
# no_ctrl_tissue_specific_df$event <- ifelse(!is.na(no_ctrl_tissue_specific_df$disease_duration), 1, 0)



tissues <- c(lumbar = "spinal_cord_lumbar", cervical = "spinal_cord_cervical", 
    thoracic = "spinal_cord_thoracic", cerebellum = "Cerebellum", 
    cortex = "cortex", motor_cortex = "cortex_motor")


require("survival")
require("survminer")
require("dplyr")
require("tidycmprsk")

survival_list <- lapply(names(tissues), function(x){
    a_tissue <- tissues[x]
    a_col <- paste(a_tissue, "group", sep = "_")

    survival_df <- perform_survival_df
    survival_df$erap2_group <- survival_df[,a_col]
    survival_df[,"erap2_group"] <- factor(survival_df[,"erap2_group"], levels = c("low", "high"))
    fit1 <- survfit(Surv(time = disease_duration, event = event, type = "right") ~ erap2_group, data = survival_df)
    file_path <- file.path(analysis_plots_dir, 
        paste("07-als-survival-analysis-of", x, "patient-survival-duration-vs-higher-erap2-expression.pdf", sep = "-"))
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
    fit1 <- survfit(Surv(time = as.numeric(age_at_death), event = event, type = "right") ~ erap2_group, data = survival_df)
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
