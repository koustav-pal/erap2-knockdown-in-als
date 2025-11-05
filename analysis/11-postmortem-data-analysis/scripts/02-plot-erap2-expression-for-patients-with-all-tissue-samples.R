data_dir <- normalizePath("/nemo/lab/patanir/home/users/palk/")
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


coldata_df_split <- split(coldata_df, paste(coldata_df$subject_id, coldata_df$tissue))
patient_plot_df <- do.call(rbind, lapply(coldata_df_split, function(temp_df){
    return_df <- unique(temp_df[,!(colnames(temp_df) %in% 
        c("X", "sample_name", "alignment_file", "alignment_file_dir",
          "isPairedEnd", "strand_specificity", "samplesheet_path", 
          "column_id", "sample_id", "mutation", "site_of_motor_onset",
          "instrument", "rin", "expr"))])
    if(nrow(return_df) > 1){
        # message(paste(unique(return_df$subject_id), unique(return_df$tissue)))
        if("N/A" %in% return_df$mnd_with_ftd & "Yes" %in% return_df$mnd_with_ftd){
            return_df <- return_df[return_df$mnd_with_ftd == "Yes",]
        }else if("N/A" %in% return_df$mnd_with_ftd & "No" %in% return_df$mnd_with_ftd){
            return_df <- return_df[return_df$mnd_with_ftd == "No",]
        }else{
            stop(paste(return_df$subject_id, return_df$tissue))
        }
        # message(nrow(return_df))
    }
    return_df$expr <- mean(temp_df$expr)
    return(return_df)
}))
patient_plot_df$tissue <- droplevels(patient_plot_df$tissue)
saveRDS(object = patient_plot_df, file = file.path(analysis_save_dir, "02-summarised-counts-for-patients.rds"))
# temp_df <- coldata_df_split[["PF-UCL-111 Cortex"]]


patient_with_multiple_tissues <- table(patient_plot_df$subject_id, patient_plot_df$tissue)
total_num_tissues <- rowSums(patient_with_multiple_tissues)

multiple_patient_df <- patient_plot_df[patient_plot_df$subject_id %in% names(total_num_tissues)[total_num_tissues == 6],]
multiple_patient_df$tissue <- factor(multiple_patient_df$tissue, 
    levels = c("Spinal_Cord_Lumbar", "Spinal_Cord_Thoracic", "Spinal_Cord_Cervical", 
        "Cortex", "Cortex_Motor", "Cerebellum"))

plot_multiple_patient_df <- do.call(rbind, lapply(split(multiple_patient_df, multiple_patient_df$subject_id), function(x){
    x$erap2_expr_zscore <- (x$expr - mean(x$expr)) / sd(x$expr)
    return(x)
}))
MinMax <- min(abs(c(min(plot_multiple_patient_df$erap2_expr_zscore), max(plot_multiple_patient_df$erap2_expr_zscore))))
plot_multiple_patient_df$scaled_zscore <- plot_multiple_patient_df$erap2_expr_zscore
plot_multiple_patient_df$scaled_zscore[plot_multiple_patient_df$erap2_expr_zscore > MinMax] <- MinMax
plot_multiple_patient_df$scaled_zscore[plot_multiple_patient_df$erap2_expr_zscore < (MinMax * -1)] <- (MinMax * -1)



colours <- c("ctrl" = "#EFF7F6", "als" = "#F87575", "als+ftd" = "#FFA9A3")
ThePlot <- ggplot(plot_multiple_patient_df, aes(x = tissue, y = subject_id))
# ThePlot <- ThePlot + geom_violin(aes(fill = als_group), colour = "#000000", 
#     lwd = 0.5, position = position_dodge(width = 0.9), scale = "width")
ThePlot <- ThePlot + geom_tile(aes(fill = scaled_zscore))
# ThePlot <- ThePlot + geom_boxplot(aes(fill = als_group), colour = "#000000", 
#     width = 0.5, lwd = 0.5, position = position_dodge(width = 0.9))
# ThePlot <- ThePlot + geom_hline(aes(yintercept = 1), colour = "#000000", lty = "dashed")
ThePlot <- ThePlot + scale_x_discrete("tissues")
ThePlot <- ThePlot + scale_y_discrete("patients")
ThePlot <- ThePlot + scale_fill_distiller("erap2 expr zscore", palette = "RdBu")
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("erap2 expression violins")
ThePlot <- ThePlot + facet_wrap(. ~ condition, ncol = 1, nrow = 3, scales = "free")
ThePlot <- ThePlot + theme(
    legend.position = "bottom",
    text = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "azure3",
    linetype = "dashed",
    size = 0.2))
ggsave(width = 7, height = 15, file = file.path(analysis_plots_dir, "01-03-plot-erap2-expression-heatmap.pdf"), ThePlot)


plot_df <- plot_multiple_patient_df[plot_multiple_patient_df$condition %in% c("als", "als+ftd", "ftd") & 
    !(plot_multiple_patient_df$disease_duration_in_months %in% "Unknown"),]
plot_df$disease_duration_in_months <- as.numeric(plot_df$disease_duration_in_months)
colours <- c("ctrl" = "#EFF7F6", "als" = "#F87575", "als+ftd" = "#FFA9A3")
ThePlot <- ggplot(plot_df, aes(x = disease_duration_in_months, y = erap2_expr_zscore))
# ThePlot <- ThePlot + geom_violin(aes(fill = als_group), colour = "#000000", 
#     lwd = 0.5, position = position_dodge(width = 0.9), scale = "width")
ThePlot <- ThePlot + geom_point(aes(colour = condition))
ThePlot <- ThePlot + geom_smooth(method = "lm", formula = y ~ x)
# ThePlot <- ThePlot + geom_boxplot(aes(fill = als_group), colour = "#000000", 
#     width = 0.5, lwd = 0.5, position = position_dodge(width = 0.9))
# ThePlot <- ThePlot + geom_hline(aes(yintercept = 1), colour = "#000000", lty = "dashed")
ThePlot <- ThePlot + scale_x_continuous("disease duration in months")
ThePlot <- ThePlot + scale_y_continuous("ERAP2 expression zscore")
ThePlot <- ThePlot + scale_colour_brewer("ALS group", palette = "Set1")
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("erap2 expression violins")
ThePlot <- ThePlot + facet_wrap(tissue ~ ., nrow = 2, ncol = 3, scales = "free")
ThePlot <- ThePlot + theme(
    legend.position = "bottom",
    text = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "azure3",
    linetype = "dashed",
    size = 0.2))
ggsave(width = 15, height = 5, file = file.path(analysis_plots_dir, "01-04-plot-erap2-expression-versus-disease-duration-by-tissue.pdf"), ThePlot)




lapply(split(plot_df, plot_df$tissue), function(temp_df){
    summary(lm(data = temp_df, formula = erap2_expr_zscore ~ disease_duration_in_months))
})

# lapply(split(plot_df, plot_df$tissue), function(temp_df){
#     cor.test(x = temp_df$disease_duration_in_months, y = temp_df$erap2_expr_zscore)
# })














# rownames(patient_plot_df) <- NULL
# patient_plot_df$tissue <- droplevels(patient_plot_df$tissue)
# wilcox_pval_df <- do.call(rbind, lapply(split(patient_plot_df, patient_plot_df$tissue), function(temp_df){
#     message(unique(temp_df$tissue))
#     als_pval <- wilcox.test(temp_df$expr[temp_df$als_group == "als"], temp_df$expr[temp_df$als_group == "ctrl"])$p.value
#     ftd_pval <- wilcox.test(temp_df$expr[temp_df$als_group == "als+ftd"], temp_df$expr[temp_df$als_group == "ctrl"])$p.value
#     # ftd_vs_als_pval <- wilcox.test(temp_df$expr[temp_df$als_group == "als+ftd"], temp_df$expr[temp_df$als_group == "als"])
#     data.frame(tissue = rep(unique(temp_df$tissue),2), als_group = c("als", "als+ftd"), pval = c(als_pval, ftd_pval))
# }))
# rownames(wilcox_pval_df) <- NULL

#                  tissue als_group        pval
# 1    Spinal_Cord_Lumbar       als 0.006507642
# 2    Spinal_Cord_Lumbar   als+ftd 0.073338136
# 3  Spinal_Cord_Cervical       als 0.096854507
# 4  Spinal_Cord_Cervical   als+ftd 0.506705846
# 5            Cerebellum       als 0.052983128
# 6            Cerebellum   als+ftd 0.045207345
# 7          Cortex_Motor       als 0.668195319
# 8          Cortex_Motor   als+ftd 0.643884441
# 9                Cortex       als 0.136765892
# 10               Cortex   als+ftd 0.022310162
# 11 Spinal_Cord_Thoracic       als 0.264686443
# 12 Spinal_Cord_Thoracic   als+ftd 0.953046953

# colours <- c("ctrl" = "#EFF7F6", "als" = "#F87575", "als+ftd" = "#FFA9A3")
# ThePlot <- ggplot(patient_plot_df, aes(x = type, y = expr))
# ThePlot <- ThePlot + geom_boxplot(aes(fill = als_group), colour = "#000000", 
#     width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9))
# ThePlot <- ThePlot + scale_x_discrete("mutation")
# ThePlot <- ThePlot + scale_y_continuous("ERAP2 expression tpm")
# ThePlot <- ThePlot + scale_fill_manual("conditions", values = colours)
# ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + ggtitle("erap2 expression violins")
# ThePlot <- ThePlot + facet_grid(tissue ~ .)
# ThePlot <- ThePlot + theme(
#     legend.position = "bottom",
#     text = element_text(size = 12),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     legend.text = element_text(size = 10),
#     legend.title = element_text(size = 10),
#     plot.background = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(colour = "azure3",
#     linetype = "dashed",
#     size = 0.2))
# ggsave(width = 15, height = 15, file = file.path(analysis_plots_dir, "01-01-plot-erap2-expression-tpm-in-mutations.pdf"), ThePlot)

# colours <- c("ctrl" = "#EFF7F6", "als" = "#F87575", "als+ftd" = "#FFA9A3")
# ThePlot <- ggplot(patient_plot_df[patient_plot_df$als_group != "ctrl",], aes(x = as.numeric(disease_duration_in_months), y = expr))
# ThePlot <- ThePlot + geom_point(aes(fill = als_group), shape = 21, colour = "#000000")
# ThePlot <- ThePlot + scale_x_continuous("age at death")
# ThePlot <- ThePlot + scale_y_continuous("ERAP2 expression tpm")
# ThePlot <- ThePlot + scale_fill_manual("tissue", values = colours)
# ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + ggtitle("erap2 expression ~ age at death")
# ThePlot <- ThePlot + facet_grid(. ~ tissue)
# ThePlot <- ThePlot + theme(
#     legend.position = "bottom",
#     text = element_text(size = 12),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     legend.text = element_text(size = 10),
#     legend.title = element_text(size = 10),
#     plot.background = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(colour = "azure3",
#     linetype = "dashed",
#     size = 0.2))
# ggsave(width = 15, height = 5, file = file.path(analysis_plots_dir, "01-02-plot-age-associated-erap2-expression-tpm.pdf"), ThePlot)

# summary_df <- summarise_tissue_specific_parameter_by_patients(coldata_df)
# saveRDS(object = summary_df, file = file.path(analysis_save_dir, "02-postmortem-patient-counts-summarised.rds"))
