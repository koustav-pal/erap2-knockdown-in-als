data_dir <- normalizePath("/export2/home/koustav/")
project_dir <- file.path(data_dir, "Projects", "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "14-01-erap2-ov-proteomics")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
erap2_ov_df <- read.csv(file.path(analysis_input_dir, "ERAP2_overexpression_proteomics.csv"))
# annotation_df <- read.csv(file.path(analysis_input_dir, "annotation.csv"))
annotation_gtf <- "/export2/home/koustav/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
# ==========================================================================
# Imports
# ==========================================================================
library("stringr")
library("GenomicRanges")
library("DESeq2")
library("DEP")
library("stringr")
library("rtracklayer")
library("reshape2")
# library("biomaRt")
# library("org.Hs.eg.db")
library("ROTS")
# ==========================================================================
# Analysis
# ==========================================================================

# annotation_ranges <- import.gff2(annotation_gtf)
unknown_proteins <- c("AG10A_HUMAN","ANXA8_HUMAN","ARF1_HUMAN","AT5EL_HUMAN","CALM1_HUMAN"
"CCZ1B_HUMAN","CDK19_HUMAN","CLD3_HUMAN","CLK1_HUMAN","DVL1_HUMAN"
"FABP4_HUMAN","F133A_HUMAN","GAL3A_HUMAN","GGT1_HUMAN","GBG5B_HUMAN"
"GPHRA_HUMAN","GREM1_HUMAN","GST2_HUMAN","T2H2L_HUMAN","H2A1D_HUMAN"
"H2A1B_HUMAN","H2AV_HUMAN","H2B1C_HUMAN","H31T_HUMAN","HBB_HUMAN"
"HBE_HUMAN","HS71A_HUMAN","IGLC2_HUMAN","IQIP1_HUMAN","LIMS3_HUMAN"
"MGN2_HUMAN","MLP3B_HUMAN","MOB1A_HUMAN","MT1G_HUMAN","MZT2A_HUMAN"
"NKX61_HUMAN","R39L5_HUMAN","RUXGL_HUMAN","TCAL3_HUMAN","TPC2A_HUMAN"
"UB2D1_HUMAN","WASH2_HUMAN","ZFX_HUMAN","Z324A_HUMAN","ZNG1A_HUMAN"
"U2AFL_HUMAN")
# reformat proteomics data
annotation_df <- reshape2::melt(as.matrix(erap2_ov_df[-1,c(41, 42, 43)]))
genes_df <- annotation_df[annotation_df$Var2 == "T..Genes",]
protein_df <- annotation_df[annotation_df$Var2 == "T..Protein.Group",]
gene_to_protein_df <- reshape2::acast(annotation_df, formula = Var1 ~ Var2, value.var = "value")

melt_erap2_ov_df <- reshape2::melt(as.matrix(erap2_ov_df[-1,c(1:27)]))
colnames(melt_erap2_ov_df) <-c("row_x", "variable", "value")

variable_split <- str_split(melt_erap2_ov_df$variable, pattern = "[_,.]", simplify = T)[,-c(1)]


melt_erap2_ov_df$line <- variable_split[,1]
melt_erap2_ov_df$erap2 <- variable_split[,2]
melt_erap2_ov_df$replicate <- variable_split[,3]
melt_erap2_ov_df$sample_name <- paste(melt_erap2_ov_df$line, melt_erap2_ov_df$erap2, sep = "_")
melt_erap2_ov_df$unique_sample_name <- paste(melt_erap2_ov_df$line, melt_erap2_ov_df$erap2, melt_erap2_ov_df$replicate, sep = "_")
melt_erap2_ov_df$Genes <- genes_df[match(melt_erap2_ov_df$row_x, genes_df$Var1),"value"]
melt_erap2_ov_df$protein_id <- protein_df[match(melt_erap2_ov_df$row_x, protein_df$Var1),"value"]
melt_erap2_ov_df$value <- as.numeric(melt_erap2_ov_df$value)


non_na_keep_melt_erap2_kd_df <-  melt_erap2_ov_df[!is.na(melt_erap2_ov_df$value) & melt_erap2_ov_df$Genes != "NA" & melt_erap2_ov_df$Genes != "",]

protein_matrix <- reshape2::acast(data = non_na_keep_melt_erap2_kd_df, formula = Genes ~ unique_sample_name, value.var = "value")
experiment_df <- unique(non_na_keep_melt_erap2_kd_df[,c("line", "erap2", "sample_name", "unique_sample_name")])


protein_sparsity <- rowSums(protein_matrix == 0 | is.na(protein_matrix)) / ncol(protein_matrix)
average_signal <- rowMeans(protein_matrix, na.rm = T)
summarised_protein_matrix <- do.call(cbind, lapply(split(experiment_df, experiment_df$sample_name), function(temp_df){
    rowSums(protein_matrix[,temp_df$unique_sample_name], na.rm = T) / ncol(protein_matrix[,temp_df$unique_sample_name])
}))
# summarised_protein_matrix[is.na(summarised_protein_matrix)] <- 0
experiment_df <- unique(experiment_df[,-4])
summarised_protein_df <- as.data.frame(summarised_protein_matrix)
summarised_protein_df$protein_accessions <- gene_to_protein_df[match(rownames(summarised_protein_matrix), gene_to_protein_df[,3]),1]
summarised_protein_df$protein_name <- gene_to_protein_df[match(rownames(summarised_protein_matrix), gene_to_protein_df[,3]),2]
summarised_protein_df$average_signal <- average_signal
summarised_protein_df$sparsity <- protein_sparsity
rownames(summarised_protein_df) <- NULL


# protein_matrix_df <- as.data.frame(protein_matrix)
# protein_matrix_df$protein_name <- rownames(protein_matrix)
# rownames(protein_matrix_df) <- NULL

unique_df <- make_unique(proteins = summarised_protein_df,
    names = "protein_name",
    ids = "protein_accessions",
    delim = ";")

se_experiment_df <- unique(experiment_df[,c("sample_name", "line", "erap2")])
se_experiment_df$condition <- se_experiment_df$erap2
se_experiment_df <- se_experiment_df[,-c(3)]
colnames(se_experiment_df) <- c("label", "replicate", "condition")



summarised_experiment <- make_se(unique_df, as.integer(c(1:9)), se_experiment_df)


# filtered_data <- filter_missval(summarised_experiment, thr = 0)
# data_norm <- normalize_vsn(filtered_data)

# pdf(file = file.path(analysis_plots_dir, "01-plot-data-normalisation.pdf"))
# plot_normalization(filtered_data, data_norm)
# plot_coverage(filtered_data)
# dev.off()


# pdf(file = file.path(analysis_plots_dir, "02-plot-missing-values.pdf"))
# plot_missval(filtered_data)
# plot_detect(filtered_data)
# dev.off()


# The data was detected as an MNAR data, meaning missing not at random. 
# Apply both imputation methods and see what happens

# Apply MinProb
# minprob_impute <- DEP::impute(data_norm, fun = "MinProb")

# # Apply QRILC
# regressed_impute <- DEP::impute(data_norm, fun = "QRILC")


# pdf(file = file.path(analysis_plots_dir, "03-plot-imputation-results.pdf"))
# plot_imputation(data_norm, minprob_impute)
# plot_imputation(data_norm, regressed_impute)
# dev.off()


current_coldata <- colData(summarised_experiment)

# make_se also does log transformation, we don't want this happening twice!
assay_data <- 2^assay(summarised_experiment)



run_differential_protein_enrichment_using_rots <- function(comparison, coldata, paired = F, my_seed = 923572, protein_matrix){
    if(length(comparison) != 1){
        stop("comparison variable must be provided...")
    }
    comparison_split <- str_split(string = comparison, pattern = "_vs_", simplify = T)[1,]
    colnames_split <- split(coldata$ID[coldata$condition %in% comparison_split], coldata$condition[coldata$condition %in% comparison_split])

    if(paired){
        line_split <- split(coldata$replicate[coldata$condition %in% comparison_split], coldata$condition[coldata$condition %in% comparison_split])
        new_colnames_split <- lapply(names(colnames_split), function(x){
            colnames_split[[x]][match(line_split[[1]], line_split[[x]])]
        })
        names(new_colnames_split) <- names(colnames_split)
        colnames_split <- new_colnames_split
    }

    denominator <- comparison_split[2]
    numerator <- comparison_split[1]

    keep_matrix <- protein_matrix[,c(colnames_split[[denominator]],colnames_split[[numerator]])]

    numerator_means <- rowMeans(keep_matrix[,colnames_split[[numerator]]])
    denominator_means <- rowMeans(keep_matrix[,colnames_split[[denominator]]])
    all_means <- rowMeans(keep_matrix)
    log_ratios <- numerator_means - denominator_means

    rots.out <- ROTS(data = keep_matrix, 
        groups = c(rep(1, length(colnames_split[[denominator]])),
            rep(0, length(colnames_split[[numerator]]))),
        paired = paired,
        log = T, 
        seed = my_seed,
        K = ceiling(nrow(keep_matrix) * 0.7))
    # rots.summary <- summary(rots.out, fdr = fdr_threshold)

    return_df <- data.frame(comparison = comparison, 
        protein_name = names(rots.out$logfc),
        logfc = rots.out$logfc,
        pvalue = rots.out$pvalue, 
        fdr = rots.out$FDR)
    return(return_df)
}


require("parallel")
comparisons <- c("ERAP2_vs_SSCFP", "SSCFP_vs_UN", "ERAP2_vs_UN")
diff_expr_quantile_df_list <- lapply(c(0), function(x){
    message(x)
    set.seed(29359263)
    assay_data <- 2^assay(summarised_experiment)
    current_coldata <- colData(summarised_experiment)
    assay_data[is.na(assay_data)] <- 0
    differential_expression_results <- mclapply(comparisons, function(a_comparison){
        if(a_comparison %in% c("ERAP2_vs_SSCFP", "SSCFP_vs_UN", "ERAP2_vs_UN")){
            paired <- T
        }else{
            paired <- F
        }
        differential_df <- run_differential_protein_enrichment_using_rots(comparison = a_comparison,
            coldata = current_coldata, 
            paired = paired,
            protein_matrix = assay_data)
        differential_df$minprob_q_threshold <- x
        return(differential_df)
    })
    names(differential_expression_results) <- comparisons
    differential_expression_df <- do.call(rbind, differential_expression_results)
    rownames(differential_expression_df) <- NULL
    return(differential_expression_df)
})

diff_expr_quantile_df <- do.call(rbind, diff_expr_quantile_df_list)
diff_expr_quantile_df$sparsity <- protein_sparsity[diff_expr_quantile_df$protein_name]
diff_expr_quantile_df$average_signal <- average_signal[diff_expr_quantile_df$protein_name]


saveRDS(object = diff_expr_quantile_df, file = file.path(analysis_save_dir, "02-01-erap2-ov-rots-differential-enrichment-analysis-raw-data.rds"))