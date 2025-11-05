data_dir <- normalizePath("/export2/home/koustav/")
project_dir <- file.path(data_dir, "Projects", "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "11-postmortem-data-analysis")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
temp_dir <- file.path(analysis_dir, "tmp")
alsod_genes_df <- read.csv(file.path(project_dir, "data_files/alsod_genelist.csv"))
annotation_gtf <- "/export2/home/koustav/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
tmp_dir <- file.path(analysis_dir, "slurm")
rsf_risk_stratified_patient_df <- file.path(save_dir, "04-01-save-patients-with-rsf-risk-based-on-disease-duration.rds")
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
source("../src/rmats-analysis/src/source.R")
# ==========================================================================
# Analysis
# ==========================================================================


# feature_ranges_with_repeats <- significant_de_event_list$repeat_containing_ranges
counts_matrix <- postmortem_list$counts


gene_level_counts_matrix <- do.call(rbind, lapply(split(feature_ranges$hash_id[feature_ranges$type == "exon"], feature_ranges$gene_id[feature_ranges$type == "exon"]), function(hash_id){
	if(length(hash_id) == 1){
		return(counts_matrix[hash_id,])
	}
	return(colSums(counts_matrix[hash_id,]))
}))

genic_count_list <- list(counts = gene_level_counts_matrix, coldata = sampleinfo_df)

saveRDS(object = genic_count_list, file = file.path(analysis_save_dir, "01-gene-level-counts-of-postmortem-samples.rds"))