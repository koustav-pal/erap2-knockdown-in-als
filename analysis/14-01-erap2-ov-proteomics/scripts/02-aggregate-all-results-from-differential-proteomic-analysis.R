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
annotation_gtf <- "/nemo/lab/patanir/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
diff_expr_quantile_df <- readRDS(file.path(analysis_save_dir, "02-01-erap2-ov-rots-differential-enrichment-analysis-raw-data.rds"))
# ==========================================================================
# Imports
# ==========================================================================
library("stringr")
library("GenomicRanges")
library("DESeq2")
library("stringr")
library("rtracklayer")
# library("biomaRt")
# library("org.Hs.eg.db")
library("qvalue")
# library("swfdr")
# ==========================================================================
# Analysis
# ==========================================================================

split_df <- split(diff_expr_quantile_df, paste(diff_expr_quantile_df$comparison, diff_expr_quantile_df$minprob_q_threshold, sep = "_"))


# compute qvalues
diff_expr_df_with_qval <- do.call(rbind, lapply(split_df, function(temp_df){
    qvalues <- qvalue(temp_df$pvalue)
    temp_df$qvalue <- qvalues$qvalues
    return(temp_df)
}))


# Plot scatter plots of pvalue vs average intensity and sparsity
test_data <- split_df[[1]]


# summary(lm(data = test_data, formula = pvalue ~ logfc + sparsity + average_signal))



require("ggplot2")
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


unimputed_data_df <- diff_expr_df_with_qval[diff_expr_df_with_qval$minprob_q_threshold == 0,]
split_df2 <- split(unimputed_data_df, unimputed_data_df$comparison)

unimputed_effect_df <- lapply(split_df2, function(temp_df){
    return_df <- data.frame(comparison = unique(temp_df$comparison),
        protein_name = temp_df$protein_name,
        logfc = temp_df$logfc,
        pvalue = temp_df$pvalue,
        num_pval_sig = temp_df$pvalue < 0.05,
        pvalue_prop = temp_df$pvalue < 0.05,
        fdr = temp_df$fdr,
        qval  = temp_df$qvalue,
        num_fdr_sig = temp_df$fdr < 0.05,
        fdr_prop = temp_df$fdr < 0.05,
        # qval = rowMeans(qval_matrix),
        # num_qval_sig = rowSums(qval_matrix < 0.05),
        # qval_prop = rowSums(qval_matrix < 0.05) / ncol(qval_matrix),
        sparsity = temp_df$sparsity,
        average_signal = temp_df$average_signal)
    # covariate_matrix <- data.frame(sparsity = return_df$sparsity,
    #     average_signal = return_df$average_signal / max(return_df$average_signal),
    #     peptide_length = return_df$peptide_length)
    # covariate_matrix$peptide_length[is.na(covariate_matrix$peptide_length)] <- 0
    # qvalues <- lm_qvalue(return_df$pvalue, X = covariate_matrix)
    # return_df$covariate_qval <- qvalues$qvalues
    return(return_df)    
})
all_effect_df <- do.call(rbind, unimputed_effect_df)

saveRDS(object = list(unimputed_data = unimputed_data_df), 
    file = file.path(analysis_save_dir, "02-01-erap2-OV-rots-differential-enrichment-analysis-results.rds"))
