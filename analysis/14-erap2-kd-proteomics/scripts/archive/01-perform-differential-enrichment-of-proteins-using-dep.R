data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects", "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "14-erap2-kd-proteomics")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
erap2_kd_df <- read.table(file.path(analysis_input_dir, "report.unique_genes_matrix.tsv"), header = T, sep = "\t")
annotation_df <- read.csv(file.path(analysis_input_dir, "annotation.csv"))
annotation_gtf <- "/nemo/lab/patanir/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
# ==========================================================================
# Imports
# ==========================================================================
library("stringr")
library("GenomicRanges")
library("DESeq2")
library("DEP")
library("stringr")
library("rtracklayer")
library("biomaRt")
library("org.Hs.eg.db")
# ==========================================================================
# Analysis
# ==========================================================================

annotation_ranges <- import.gff2(annotation_gtf)

# reformat proteomics data
melt_erap2_kd_df <- reshape2::melt(erap2_kd_df)
colnames(melt_erap2_kd_df) <-c("Genes", "variable", "value")

variable_split <- str_split(melt_erap2_kd_df$variable, pattern = "[_,.]", simplify = T)[,-c(1)]

annotation_df$sample_name <- paste(annotation_df[,2], annotation_df[,3], annotation_df[,4], annotation_df[,5], sep = "_")


rep_info <- unlist(lapply(split(variable_split[,7], variable_split[,5]), function(x){
    return_list <- seq_along(unique(x))
    names(return_list) <- unique(x)
    return(return_list)
}))
names(rep_info) <- str_split(names(rep_info), pattern = "\\.", simplify = T)[,2]




melt_erap2_kd_df$line <- annotation_df$Line[as.numeric(variable_split[,5])]
melt_erap2_kd_df$genotype <- annotation_df$Genotype[as.numeric(variable_split[,5])]
melt_erap2_kd_df$treatment <- annotation_df$Treatment[as.numeric(variable_split[,5])]
melt_erap2_kd_df$erap2 <- annotation_df$siRNA[as.numeric(variable_split[,5])]
melt_erap2_kd_df$replicate <- rep_info[variable_split[,7]]
melt_erap2_kd_df$sample_name <- annotation_df$sample_name[as.numeric(variable_split[,5])]
melt_erap2_kd_df$unique_sample_name <- paste(annotation_df$sample_name[as.numeric(variable_split[,5])], melt_erap2_kd_df$replicate, sep = "_")


keep_melt_erap2_kd_df <- melt_erap2_kd_df[melt_erap2_kd_df$treatment != "TIC" & melt_erap2_kd_df$erap2 != "Lipo",]
non_na_keep_melt_erap2_kd_df <-  keep_melt_erap2_kd_df[!is.na(keep_melt_erap2_kd_df$value) & keep_melt_erap2_kd_df$Genes != "NA",]

protein_matrix <- reshape2::acast(data = non_na_keep_melt_erap2_kd_df, formula = Genes ~ unique_sample_name, value.var = "value")
experiment_df <- unique(keep_melt_erap2_kd_df[,c("line", "genotype", "treatment", "erap2", "sample_name", "unique_sample_name")])


summarised_protein_matrix <- do.call(cbind, lapply(split(experiment_df, experiment_df$sample_name), function(temp_df){
    rowSums(protein_matrix[,temp_df$unique_sample_name], na.rm = T) / ncol(protein_matrix)
}))
experiment_df <- unique(experiment_df[,-6])
summarised_protein_df <- as.data.frame(summarised_protein_matrix)
summarised_protein_df$protein_name <- rownames(summarised_protein_matrix)
rownames(summarised_protein_df) <- NULL

protein_matrix_df <- as.data.frame(protein_matrix)
protein_matrix_df$protein_name <- rownames(protein_matrix)
rownames(protein_matrix_df) <- NULL



mart <- useEnsembl("ensembl","hsapiens_gene_ensembl")
get_ids_df <- getBM(c("entrezgene_id","ensembl_gene_id","uniprotswissprot","description","hgnc_symbol"), "hgnc_symbol", 
    summarised_protein_df$protein_name, mart)

protein_accessions <- vapply(split(get_ids_df$uniprotswissprot, get_ids_df$hgnc_symbol), function(x){
    paste(unique(x[x != ""]), collapse = ";")
},"1")
summarised_protein_df$protein_accessions <- protein_accessions[summarised_protein_df$protein_name]
protein_matrix_df$protein_accessions <- protein_accessions[protein_matrix_df$protein_name]

unique_df <- make_unique(proteins = protein_matrix_df,
    names = "protein_name",
    ids = "protein_accessions",
    delim = ";")

se_experiment_df <- unique(experiment_df[,c("unique_sample_name", "line", "genotype", "erap2")])
se_experiment_df$condition <- paste(se_experiment_df$genotype, se_experiment_df$erap2, sep = "_")
se_experiment_df <- se_experiment_df[,-c(3,4)]
colnames(se_experiment_df) <- c("label", "replicate", "condition")



summarised_experiment <- make_se(unique_df, as.integer(c(1:16)), se_experiment_df)

filtered_data <- filter_missval(summarised_experiment, thr = 0)

data_norm <- normalize_vsn(filtered_data)

pdf(file = file.path(analysis_plots_dir, "01-plot-data-normalisation.pdf"))
plot_normalization(filtered_data, data_norm)
plot_coverage(filtered_data)
dev.off()


pdf(file = file.path(analysis_plots_dir, "02-plot-missing-values.pdf"))
plot_missval(filtered_data)
plot_detect(filtered_data)
dev.off()


# The data was detected as an MNAR data, meaning missing not at random. 
# Apply both imputation methods and see what happens

# Apply MinProb
minprob_impute <- impute(data_norm, fun = "MinProb", q = 0.01)

# Apply QRILC
regressed_impute <- impute(data_norm, fun = "QRILC")


pdf(file = file.path(analysis_plots_dir, "03-plot-imputation-results.pdf"))
plot_imputation(data_norm, minprob_impute)
plot_imputation(data_norm, regressed_impute)
dev.off()


differential_enrichment <- test_diff(regressed_impute, type = "manual", 
    test = c("VCP_CTRL_vs_CTRL_CTRL", "CTRL_ERAP2_vs_CTRL_CTRL","VCP_ERAP2_vs_VCP_CTRL", "VCP_ERAP2_vs_CTRL_ERAP2"))

dep <- add_rejections(differential_enrichment, alpha = 0.05, lfc = log2(1.5))


# Make PCA plot

pdf(file = file.path(analysis_plots_dir, "04-plot-PCA.pdf"))
plot_pca(dep, x = 1, y = 2, n = 500, indicate = "condition", point_size = 4)
dev.off()


# Plot correlation

pdf(file = file.path(analysis_plots_dir, "05-plot-correlation.pdf"))
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
dev.off()


# Plot heatmap
pdf(file = file.path(analysis_plots_dir, "06-plot-heatmap.pdf"))
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 3, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition"))
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 3, col_limit = 3, show_row_names = FALSE)
dev.off()

pdf(file = file.path(analysis_plots_dir, "07-plot-volcano.pdf"))
plot_volcano(dep, contrast = "VCP_CTRL_vs_CTRL_CTRL", label_size = 1, adjusted = F, add_names = TRUE)
plot_volcano(dep, contrast = "CTRL_ERAP2_vs_CTRL_CTRL", label_size = 1, adjusted = F, add_names = TRUE)
plot_volcano(dep, contrast = "VCP_ERAP2_vs_VCP_CTRL", label_size = 1, adjusted = F, add_names = TRUE)
plot_volcano(dep, contrast = "VCP_ERAP2_vs_CTRL_ERAP2", label_size =1, adjusted = F, add_names = TRUE)
dev.off()



pdf(file = file.path(analysis_plots_dir, "08-plot-frequency.pdf"))
plot_cond(dep)
dev.off()


data_results <- get_results(dep)
saveRDS(object = data_results, file = file.path(analysis_save_dir, "01-save-differential-protein-enrichment-results.rds"))
write.csv(data_results, file = file.path(analysis_save_dir, "01-save-differential-protein-enrichment-results.csv"))

