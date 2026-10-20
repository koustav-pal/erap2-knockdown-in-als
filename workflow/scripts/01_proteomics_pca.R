# Adapted analysis; see provenance/source_ranges.json.
source(here::here("workflow", "scripts", "common.R"))
library(stringr)
library(reshape2)
library(DEP)
library(PCAtools)
library(SummarizedExperiment)
erap2_kd_df <- read.delim(input("proteomics_matrix.tsv"), check.names = TRUE)
annotation_df <- read.csv(input("proteomics_annotation.csv"))
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
experiment_df <- unique(experiment_df[,-6])
summarised_protein_df <- as.data.frame(summarised_protein_matrix)
summarised_protein_df$protein_name <- rownames(summarised_protein_matrix)
summarised_protein_df$average_signal <- average_signal
summarised_protein_df$sparsity <- protein_sparsity
rownames(summarised_protein_df) <- NULL

get_ids_df <- readRDS(input("proteomics_gene_mapping.rds"))
protein_accessions <- vapply(split(get_ids_df$uniprotswissprot, get_ids_df$hgnc_symbol), function(x){
    paste(unique(x[x != ""]), collapse = ";")
},"1")
summarised_protein_df$protein_accessions <- protein_accessions[summarised_protein_df$protein_name]

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

minprob_impute <- DEP::impute(data_norm, fun = "MinProb")

current_coldata <- colData(minprob_impute)
assay_data <- assay(minprob_impute)
current_coldata$mutation <- str_split(current_coldata$condition, pattern = "_", simplify = T)[,1]
current_coldata$treatment <- str_split(current_coldata$condition, pattern = "_", simplify = T)[,2]
current_coldata$treatment[current_coldata$treatment == "CTRL"] <- "Scr"
rownames(current_coldata) <- paste(current_coldata$replicate, current_coldata$treatment, sep = "_")

colnames(assay_data) <- rownames(current_coldata)[match(colnames(assay_data), current_coldata$ID)]

pca_object <- pca(assay_data, metadata = current_coldata, removeVar = 0.1)

pdf(file.path(analysis_plots_dir, "04-plot-PCA.pdf"))
  biplot(pca_object, showLoadings = F,
    labSize = 3,
    pointSize = 3,
    sizeLoadingsNames = 3,
    colby = "mutation",
    legendPosition = "bottom",
    hline = 0,
    vline = 0)

writeLines(capture.output(sessionInfo()), file.path(plot_dir, "sessionInfo.txt"))
