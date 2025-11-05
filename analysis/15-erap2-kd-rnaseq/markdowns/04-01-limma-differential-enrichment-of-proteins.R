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
library("limma")
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


keep_melt_erap2_kd_df <- melt_erap2_kd_df[melt_erap2_kd_df$treatment != "TIC" & melt_erap2_kd_df$erap2 != "Lipo" & melt_erap2_kd_df$line != "CTRL4",]
non_na_keep_melt_erap2_kd_df <-  keep_melt_erap2_kd_df[!is.na(keep_melt_erap2_kd_df$value) & keep_melt_erap2_kd_df$Genes != "NA",]

protein_matrix <- reshape2::acast(data = non_na_keep_melt_erap2_kd_df, formula = Genes ~ unique_sample_name, value.var = "value")
experiment_df <- unique(keep_melt_erap2_kd_df[,c("line", "genotype", "treatment", "erap2", "sample_name", "unique_sample_name")])


protein_sparsity <- rowSums(protein_matrix == 0 | is.na(protein_matrix)) / ncol(protein_matrix)
average_signal <- rowMeans(protein_matrix, na.rm = T)
summarised_protein_matrix <- do.call(cbind, lapply(split(experiment_df, experiment_df$sample_name), function(temp_df){
    rowSums(protein_matrix[,temp_df$unique_sample_name], na.rm = T) / ncol(protein_matrix[,temp_df$unique_sample_name])
}))
# summarised_protein_matrix[is.na(summarised_protein_matrix)] <- 0
experiment_df <- unique(experiment_df[,-6])
summarised_protein_df <- as.data.frame(summarised_protein_matrix)
summarised_protein_df$protein_name <- rownames(summarised_protein_matrix)
summarised_protein_df$average_signal <- average_signal
summarised_protein_df$sparsity <- protein_sparsity
rownames(summarised_protein_df) <- NULL


# protein_matrix_df <- as.data.frame(protein_matrix)
# protein_matrix_df$protein_name <- rownames(protein_matrix)
# rownames(protein_matrix_df) <- NULL



mart <- useEnsembl("ensembl","hsapiens_gene_ensembl")
get_ids_df <- getBM(c("entrezgene_id","ensembl_gene_id","uniprotswissprot","description","hgnc_symbol", "peptide"), "hgnc_symbol", 
    summarised_protein_df$protein_name, mart)

protein_accessions <- vapply(split(get_ids_df$uniprotswissprot, get_ids_df$hgnc_symbol), function(x){
    paste(unique(x[x != ""]), collapse = ";")
},"1")
summarised_protein_df$protein_accessions <- protein_accessions[summarised_protein_df$protein_name]
# protein_matrix_df$protein_accessions <- protein_accessions[protein_matrix_df$protein_name]

unique_df <- make_unique(proteins = summarised_protein_df,
    names = "protein_name",
    ids = "protein_accessions",
    delim = ";")

se_experiment_df <- unique(experiment_df[,c("sample_name", "line", "genotype", "erap2")])
se_experiment_df$condition <- paste(se_experiment_df$genotype, se_experiment_df$erap2, sep = "_")
se_experiment_df <- se_experiment_df[,-c(3,4)]
colnames(se_experiment_df) <- c("label", "replicate", "condition")



summarised_experiment <- make_se(unique_df, as.integer(c(1:14)), se_experiment_df)


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
minprob_impute <- DEP::impute(data_norm, fun = "MinProb")

# Apply QRILC
regressed_impute <- DEP::impute(data_norm, fun = "QRILC")


pdf(file = file.path(analysis_plots_dir, "03-plot-imputation-results.pdf"))
plot_imputation(data_norm, minprob_impute)
plot_imputation(data_norm, regressed_impute)
dev.off()


current_coldata <- colData(filtered_data)
assay_data <- assay(filtered_data)

# assay_data


current_coldata$ERAP <- ifelse(str_split(current_coldata$condition, pattern = "_", simplify = T)[,2] == "ERAP2", "ERAP2", "CTRL")
current_coldata$mutation <- ifelse(str_split(current_coldata$condition, pattern = "_", simplify = T)[,1] == "CTRL", "CTRL", "VCP")


design <- model.matrix(~ 0 + replicate + mutation + ERAP + mutation:ERAP, current_coldata)
colnames(design) <- make.names(colnames(design))

fit <- lmFit(assay_data, design)
contrast_matrix <- makeContrasts(als_vs_control = mutationVCP - mutationCTRL, 
    erap2_effect_in_control = ERAPERAP2 - mutationCTRL,
    erap2_effect_in_als = mutationVCP.ERAPERAP2 - mutationVCP, 
    erap2_effect_on_als_vs_control = (mutationVCP - mutationCTRL) + (mutationVCP.ERAPERAP2 - mutationVCP), levels = design)


fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = FALSE)

als_vs_ctrl_df <- topTable(fit2, coef="als_vs_control", number=Inf, adjust.method="BH", p.value = 0.1)
erap2_effect_in_control_df <- topTable(fit2, coef="erap2_effect_in_control", number=Inf, adjust.method="BH", p.value = 0.01)
erap2_effect_in_als_df <- topTable(fit2, coef="erap2_effect_in_als", number=Inf, adjust.method="BH", p.value = 0.01)
erap2_effect_on_als_vs_control_df <- topTable(fit2, coef="erap2_effect_on_als_vs_control", number=Inf, adjust.method="BH", p.value = 0.01)

nrow(als_vs_ctrl_df)
nrow(erap2_effect_in_control_df)
nrow(erap2_effect_in_als_df)
nrow(erap2_effect_on_als_vs_control_df)

results <- nrow(topTable(fit2, coef="erap2_effect_in_control", number=Inf, adjust.method="BH"))




require("parallel")
comparisons <- c("VCP_CTRL_vs_CTRL_CTRL", "CTRL_ERAP2_vs_CTRL_CTRL", "VCP_ERAP2_vs_VCP_CTRL", "VCP_ERAP2_vs_CTRL_ERAP2")

saveRDS(object = diff_expr_quantile_df, file = file.path(analysis_save_dir, "02-01-limma-differential-enrichment-analysis-raw-data.rds"))



