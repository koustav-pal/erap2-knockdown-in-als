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
umap_clustering_list <- readRDS(file.path(analysis_save_dir, "01-patient-latent-factor-umap-clustering.rds"))
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
require("umap")
source("scripts/src/functions.R")
# source("../src/rmats-analysis/src/source.R")
# ==========================================================================
# Analysis
# ==========================================================================
annotation_ranges <- import.gff2(annotation_gtf)
annotation_gene_ranges <- annotation_ranges[annotation_ranges$type == "gene"]


coldata_df <- genic_count_list$coldata


clustering_df <- umap_clustering_list$umap
patient_matrix <- umap_clustering_list$counts

clustering_df$subject_id <- str_split(clustering_df$patients, pattern = " ", simplify = T)[,1]
clustering_df$sex_genotype <- coldata_df$sex_genotype[match(clustering_df$subject_id, coldata_df$subject_id)]
clustering_df$age_at_symptom_onset <- coldata_df$age_at_symptom_onset[match(clustering_df$subject_id, coldata_df$subject_id)]
clustering_df$age_at_death <- coldata_df$age_at_death[match(clustering_df$subject_id, coldata_df$subject_id)]
clustering_df$disease_duration_in_months <- coldata_df$disease_duration_in_months[match(clustering_df$subject_id, coldata_df$subject_id)]




tissue_to_cluster_df <- as.matrix(table(clustering_df$tissue, clustering_df$cluster))
tissue_to_cluster_prop_df <- as.data.frame(tissue_to_cluster_df / colSums(tissue_to_cluster_df))
colnames(tissue_to_cluster_prop_df) <- c("tissue", "cluster", "proportion")

tissue_to_cluster_prop_df$tissue <- factor(tissue_to_cluster_prop_df$tissue, 
	levels = c("Spinal_Cord_Lumbar", "Spinal_Cord_Thoracic", "Spinal_Cord_Cervical",
		"Cortex", "Cortex_Motor", "Cerebellum"))


ThePlot <- ggplot(tissue_to_cluster_prop_df, aes(x = cluster, y = tissue))
# ThePlot <- ThePlot + geom_vline(xintercept = 0, linetype = "dashed", colour = "#000000")
# ThePlot <- ThePlot + geom_hline(yintercept = 0, linetype = "dashed", colour = "#000000")
ThePlot <- ThePlot + geom_tile(aes(fill = proportion))
ThePlot <- ThePlot + scale_x_discrete("cluster")
ThePlot <- ThePlot + scale_y_discrete("tissue types")
# ThePlot <- ThePlot + scale_colour_brewer("genotypes", palette = "Set1")
ThePlot <- ThePlot + scale_fill_distiller("proportion of patient in cluster", palette = "Greys", direction = 1)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ tissue)
ThePlot <- ThePlot + ggtitle("Proportion of tissue type in each cluster")
ThePlot <- ThePlot + current_theme
ggsave(height = 7, width = 7, file = file.path(analysis_plots_dir, "06-plot-proportion-of-tissue-type-in-each-cluster.pdf"), ThePlot)
