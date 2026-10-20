data_dir <- normalizePath("/data/shared/lab_internal")
project_dir <- file.path(data_dir, "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "21-marker-gene-analysis")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
temp_dir <- file.path(analysis_dir, "tmp")
dir.create(analysis_save_dir, recursive = T)
dir.create(analysis_input_dir, recursive = T)
dir.create(analysis_plots_dir, recursive = T)
dir.create(temp_dir, recursive = T)
gene_level_count_list <- readRDS("nfcore/counts/gene_level_counts.rds")
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
require("ComplexHeatmap")
require("circlize")
require("viridis")
# ==========================================================================
# Analysis
# ==========================================================================

gene_count_matrix <- gene_level_count_list$counts
row_ranges <- gene_level_count_list$row_ranges
sampleinfo_df <- gene_level_count_list$sampleinfo
meta <- gene_level_count_list$meta
annotation_df <- gene_level_count_list$annotation


sampleinfo_df$groups <- "bradley_et_al_astrocytes"
sampleinfo_df$groups[grepl(x = sampleinfo_df$sample, pattern = "astrocyte_adult")] <- "adult_astrocytes"
sampleinfo_df$groups[sampleinfo_df$sample %in% 
	c("C1_CTRL", "C2_CTRL", "C6_CTRL")] <- "control_scramble_astrocytes"
sampleinfo_df$groups[sampleinfo_df$sample %in% 
	c("C1_ERAP2", "C2_ERAP2", "C6_ERAP2")] <- "control_erap2kd_astrocytes"
sampleinfo_df$groups[sampleinfo_df$sample %in% 
	c("GLiB_CTRL", "NCRME6_CTRL", "NCRM6_CTRL", "GLiA_CTRL")] <- "als_scramble_astrocytes"
sampleinfo_df$groups[sampleinfo_df$sample %in% 
	c("GLiB_ERAP2", "NCRME6_ERAP2", "NCRM6_ERAP2", "GLiA_ERAP2")] <- "als_erap2kd_astrocytes"
sampleinfo_df$groups <- factor(sampleinfo_df$groups, 
	levels = c("adult_astrocytes", "bradley_et_al_astrocytes", "control_scramble_astrocytes",
	"control_erap2kd_astrocytes", "als_scramble_astrocytes", "als_erap2kd_astrocytes"))


marker_genes <- c("SLC7A2", "MGST1", "NFIX", "VIM", "TNC",
"FGFR3", "CYBRD1", "CCDC80", "EZR", "MLC1", "CHRDL1",
"ACSBG1", "GLI3", "SLC9A3R1", "RFX4", "F3", "SOX9",
"GFAP", "EDNRB", "BMPR1B", "EGFR", "TMEM47", "GJA1",
"CBS", "S100B", "PBXIP1", "S1PR1", "AQP4", "ID4",
"FZD2", "EMP2")

deseq_object <- DESeqDataSetFromMatrix(countData = round(gene_count_matrix[,sampleinfo_df$sample]), 
    colData = sampleinfo_df, 
    design = ~ groups)
differential_object <- DESeq(deseq_object)
normalised_counts_matrix <- assay(vst(differential_object, blind = T))

sampleinfo_split <- split(sampleinfo_df$sample, sampleinfo_df$groups)
average_vst_matrix <- do.call(cbind, lapply(sampleinfo_split, function(x){
	rowMeans(normalised_counts_matrix[,x])
}))

selected_annotation_df <- annotation_df[annotation_df$gene_name %in% marker_genes,]
selected_average_vst_matrix <- average_vst_matrix[selected_annotation_df$gene_id,]
rownames(selected_average_vst_matrix) <- selected_annotation_df$gene_name[match(rownames(selected_average_vst_matrix), selected_annotation_df$gene_id)]


zscore_matrix <- (selected_average_vst_matrix - rowMeans(selected_average_vst_matrix)) / rowSds(selected_average_vst_matrix)

col_fun = colorRamp2(unique(c(seq(from = -1, to = 0, length.out = 5), seq(from = 0, to = 1, length.out = 5))), magma(n = 9))

# col_anno <- HeatmapAnnotation(condition = rep("ALS",ncol(predictor_zscore_matrix)),
#     col = list(condition = c(
#     	# "sALS" = "#D81E5B", 
#     	"ALS" = "#F0544F", 
#     	# "FUS-ALS" = "#F0544F", "SOD1-ALS" = "#F0544F",
#     	# "TARDBP-ALS" = "#F0544F", 
#     	"AD" = "#8CBCB9", 
#     	"PD" = "#DDA448",
#     	"control" = "#C5CBD3")))

# row_anno <- rowAnnotation(
# # predictor_group = predictor_ranges$predictor[match(rownames(predictor_zscore_matrix), predictor_ranges$hash_id)],
# predictor_group = predictor_df$predictor[match(rownames(predictor_zscore_matrix), predictor_df$gene_id)],
# col = list(predictor_group = c("strong_predictors" = "#F96E46", "moderate_predictors" = "#5688C7", "modest_predictors" = "#FFE3E3")))

# row_split <- predictor_df$predictor[match(rownames(predictor_zscore_matrix), predictor_df$gene_id)]
# row_split <- factor(row_split, levels = c("strong_predictors", "moderate_predictors", "modest_predictors"))


# pdf(file = file.path(analysis_plots_dir, "01-plot-marker-gene-expression-matrix.pdf"), width = 5, height = 15)
# 	Heatmap(zscore_matrix, name = "Marker expression heatmap", 
# 	top_annotation = NULL,
# 	left_annotation = NULL,
# 	# row_split = row_split,
# 	# column_split = col_split,
# 	col = col_fun,
# 	show_row_names = T,
# 	show_column_names = T,
# 	cluster_rows = T, 
# 	cluster_columns = T) 
# dev.off()

