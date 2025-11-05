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
# rmats_data_list <- readRDS(file.path(analysis_save_dir, "01-imported-rmats-object-list.rds"))
tmp_dir <- file.path(analysis_dir, "slurm")
genic_count_list <- readRDS(file.path(analysis_save_dir, "01-gene-level-counts-of-postmortem-samples.rds"))
candidate_de_genes_affected_by_erap2_kd <- readRDS("../02-01-differential-expression-profiling/rdata_files/02-significantly-affected-genes-by-kd-in-als-vs-control.rds")
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
require("umap")
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

subset_matrix <- tpm_matrix[rownames(tpm_matrix) %in% annotation_gene_ranges$gene_id[(annotation_gene_ranges$gene_name == "ERAP2"  | 
    annotation_gene_ranges$gene_id %in% candidate_de_genes_affected_by_erap2_kd$gene_ids[candidate_de_genes_affected_by_erap2_kd$selected_by_significance]) & !is.na(annotation_gene_ranges$gene_name)],
    coldata_df$column_id]


cols_split_by_patients <- split(colnames(subset_matrix),paste(
        coldata_df$subject_id[match(colnames(subset_matrix), coldata_df$column_id)],
        coldata_df$tissue[match(colnames(subset_matrix), coldata_df$column_id)]))

patient_matrix <- do.call(cbind, lapply(cols_split_by_patients, function(x){
    if(length(x) == 1){
        return(subset_matrix[,x])
    }
    return(rowMeans(subset_matrix[,x]))
}))
patient_erap2_expr <- patient_matrix["ENSG00000164308",]
patient_coldata <- data.frame(column_id = colnames(patient_matrix),
    patient_id = str_split(colnames(patient_matrix), pattern = " ", simplify = T)[,1],
    tissue = str_split(colnames(patient_matrix), pattern = " ", simplify = T)[,2],
    erap2_expr = patient_erap2_expr[colnames(patient_matrix)])


pca_object <- prcomp(patient_matrix, scale = TRUE)
pca_summary <- summary(pca_object)
loadings <- pca_object$rotation
loadings_df <- as.data.frame(loadings[,c("PC1", "PC2", "PC3", "PC4", "PC5")])
loadings_df$erap2_expr <- patient_coldata$erap2_expr[match(rownames(loadings_df), patient_coldata$column_id)]
loadings_df$patient_id <- patient_coldata$patient_id[match(rownames(loadings_df), patient_coldata$column_id)]
loadings_df$mnd_with_ftd <- coldata_df$mnd_with_ftd[match(loadings_df$patient_id, coldata_df$subject_id)]
loadings_df$mnd_with_dementia <- coldata_df$mnd_with_dementia[match(loadings_df$patient_id, coldata_df$subject_id)]
loadings_df$tissue <- patient_coldata$tissue[match(rownames(loadings_df), patient_coldata$column_id)]
loadings_df$type <- coldata_df$type[match(loadings_df$patient_id, coldata_df$subject_id)]
loadings_df$als_group <- coldata_df$als_group[match(loadings_df$patient_id, coldata_df$subject_id)]


pca_object <- prcomp(patient_matrix, scale = TRUE)
pca_summary <- summary(pca_object)
variance_df <- data.frame(pc = colnames(pca_summary$importance), 
    variance_explained = pca_summary$importance[2,])

require(ggplot2)
current_theme <- theme(
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


ThePlot <- ggplot(loadings_df, aes(x = PC1, y = PC2, colour = tissue))
# ThePlot <- ThePlot + geom_vline(xintercept = 0, linetype = "dashed", colour = "#000000")
# ThePlot <- ThePlot + geom_hline(yintercept = 0, linetype = "dashed", colour = "#000000")
ThePlot <- ThePlot + geom_point()
ThePlot <- ThePlot + scale_colour_brewer("genotypes", palette = "Set1")
# ThePlot <- ThePlot + scale_colour_distiller("erap2 expr", palette = "Reds", direction = 1)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ tissue)
ThePlot <- ThePlot + ggtitle("principal component analysis of erap2 associated genes")
ThePlot <- ThePlot + current_theme
ggsave(height = 7, width = 7, file = file.path(analysis_plots_dir, "04-plot-pca-of-erap2-associated-genes.pdf"), ThePlot)




all_correlation_values <- do.call(c, lapply(seq_len(ncol(loadings)), function(x){
    cor(patient_matrix["ENSG00000164308",], loadings[,x])
}))
names(all_correlation_values) <- colnames(loadings)
all_correlation_values <- all_correlation_values[order(all_correlation_values, decreasing = T)]




latent_factors <- loadings[,1:10]
umap_object <- umap(latent_factors)
graph <- bluster::makeSNNGraph(latent_factors)
clust <- igraph::cluster_louvain(graph)


plot_data <- data.frame(patients = rownames(umap_object$layout), 
    xaxis = umap_object$layout[,1], 
    yaxis = umap_object$layout[,2],
    cluster = factor(clust$membership))
plot_data$erap2_expr <- patient_coldata$erap2_expr[match(rownames(plot_data), patient_coldata$column_id)]
plot_data$patient_id <- patient_coldata$patient_id[match(rownames(plot_data), patient_coldata$column_id)]
plot_data$mnd_with_ftd <- coldata_df$mnd_with_ftd[match(plot_data$patient_id, coldata_df$subject_id)]
plot_data$mnd_with_dementia <- coldata_df$mnd_with_dementia[match(plot_data$patient_id, coldata_df$subject_id)]
plot_data$tissue <- patient_coldata$tissue[match(rownames(plot_data), patient_coldata$column_id)]
plot_data$type <- coldata_df$type[match(plot_data$patient_id, coldata_df$subject_id)]
plot_data$als_group <- coldata_df$als_group[match(plot_data$patient_id, coldata_df$subject_id)]


require(ggplot2)
current_theme <- theme(
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

plot_data$plot_erap_expr <- plot_data$erap2_expr
plot_data$plot_erap_expr[plot_data$erap2_expr > 1] <- 1
ThePlot <- ggplot(plot_data, aes(x = xaxis, y = yaxis, colour = tissue))
# ThePlot <- ThePlot + geom_vline(xintercept = 0, linetype = "dashed", colour = "#000000")
# ThePlot <- ThePlot + geom_hline(yintercept = 0, linetype = "dashed", colour = "#000000")
ThePlot <- ThePlot + geom_point()
ThePlot <- ThePlot + scale_x_continuous("umap 1")
ThePlot <- ThePlot + scale_y_continuous("umap 2")
ThePlot <- ThePlot + scale_colour_brewer("genotypes", palette = "Set1")
# ThePlot <- ThePlot + scale_colour_distiller("erap2 expr", palette = "Reds", direction = 1)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ tissue)
ThePlot <- ThePlot + ggtitle("umap analysis of erap2 associated genes")
ThePlot <- ThePlot + current_theme
ggsave(height = 7, width = 7, file = file.path(analysis_plots_dir, "05-02-plot-umap-of-erap2-associated-factors-highlight-tissue.pdf"), ThePlot)


ThePlot <- ggplot(plot_data, aes(x = xaxis, y = yaxis, colour = as.character(cluster)))
# ThePlot <- ThePlot + geom_vline(xintercept = 0, linetype = "dashed", colour = "#000000")
# ThePlot <- ThePlot + geom_hline(yintercept = 0, linetype = "dashed", colour = "#000000")
ThePlot <- ThePlot + geom_point()
ThePlot <- ThePlot + scale_x_continuous("umap 1")
ThePlot <- ThePlot + scale_y_continuous("umap 2")
ThePlot <- ThePlot + scale_colour_brewer("genotypes", palette = "Set1")
# ThePlot <- ThePlot + scale_colour_distiller("erap2 expr", palette = "Reds", direction = 1)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ tissue)
ThePlot <- ThePlot + ggtitle("umap analysis of erap2 associated genes")
ThePlot <- ThePlot + current_theme
ggsave(height = 7, width = 7, file = file.path(analysis_plots_dir, "05-02-plot-umap-of-erap2-associated-factors-highlight-clusters.pdf"), ThePlot)


ThePlot <- ggplot(plot_data, aes(x = xaxis, y = yaxis, colour = plot_erap_expr))
# ThePlot <- ThePlot + geom_vline(xintercept = 0, linetype = "dashed", colour = "#000000")
# ThePlot <- ThePlot + geom_hline(yintercept = 0, linetype = "dashed", colour = "#000000")
ThePlot <- ThePlot + geom_point()
ThePlot <- ThePlot + scale_x_continuous("umap 1")
ThePlot <- ThePlot + scale_y_continuous("umap 2")
# ThePlot <- ThePlot + scale_colour_brewer("genotypes", palette = "Set1")
ThePlot <- ThePlot + scale_colour_distiller("erap2 expr", palette = "Reds", direction = 1)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ tissue)
ThePlot <- ThePlot + ggtitle("umap analysis of erap2 associated genes")
ThePlot <- ThePlot + current_theme
ggsave(height = 7, width = 7, file = file.path(analysis_plots_dir, "05-02-plot-umap-of-erap2-associated-factors-highlight-erap2-expression.pdf"), ThePlot)



saveRDS(object = list(umap = plot_data, counts = patient_matrix), file = file.path(analysis_save_dir, "01-patient-latent-factor-umap-clustering.rds"))
