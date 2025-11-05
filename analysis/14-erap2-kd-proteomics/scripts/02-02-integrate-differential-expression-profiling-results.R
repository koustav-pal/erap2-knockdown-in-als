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
diff_expr_quantile_df <- readRDS(file.path(analysis_save_dir, "02-01-rots-differential-enrichment-analysis-raw-data.rds"))
# ==========================================================================
# Imports
# ==========================================================================
library("stringr")
library("GenomicRanges")
library("DESeq2")
library("stringr")
library("rtracklayer")
library("biomaRt")
library("org.Hs.eg.db")
library("qvalue")
library("swfdr")
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


ThePlot <- ggplot(test_data, aes(x = -log10(pvalue), y = sparsity))
# ThePlot <- ThePlot + geom_violin(aes(fill = als_group), colour = "#000000", 
#     lwd = 0.5, position = position_dodge(width = 0.9), scale = "width")
ThePlot <- ThePlot + geom_point(colour = "#F5EFED")
ThePlot <- ThePlot + scale_x_continuous("pvalue")
ThePlot <- ThePlot + scale_y_continuous("sparsity")
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("pvalue vs sparsity")
ThePlot <- ThePlot + current_theme
ggsave(width = 7, height = 7, file = file.path(analysis_plots_dir, "02-01-plot-pvalue-vs-sparsity.pdf"), ThePlot)

ThePlot <- ggplot(test_data, aes(x = -log10(pvalue), y = log10(average_signal)))
# ThePlot <- ThePlot + geom_violin(aes(fill = als_group), colour = "#000000", 
#     lwd = 0.5, position = position_dodge(width = 0.9), scale = "width")
ThePlot <- ThePlot + geom_point(colour = "#F5EFED")
ThePlot <- ThePlot + scale_x_continuous("pvalue")
ThePlot <- ThePlot + scale_y_continuous("sparsity")
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("pvalue vs sparsity")
ThePlot <- ThePlot + current_theme
ggsave(width = 7, height = 7, file = file.path(analysis_plots_dir, "02-02-plot-pvalue-vs-average-signal.pdf"), ThePlot)


# cor.test(test_data$pvalue, test_data$sparsity, method = "pearson")
# cor.test(test_data$pvalue, test_data$average_signal, method = "pearson")
# cor.test(test_data$pvalue, test_data$peptide_length, method = "pearson")

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
        num_fdr_sig = temp_df$fdr < 0.05,
        fdr_prop = temp_df$fdr < 0.05,
        # qval = rowMeans(qval_matrix),
        # num_qval_sig = rowSums(qval_matrix < 0.05),
        # qval_prop = rowSums(qval_matrix < 0.05) / ncol(qval_matrix),
        sparsity = temp_df$sparsity,
        average_signal = temp_df$average_signal,
        peptide_length = temp_df$peptide_length)
    covariate_matrix <- data.frame(sparsity = return_df$sparsity,
        average_signal = return_df$average_signal / max(return_df$average_signal),
        peptide_length = return_df$peptide_length)
    covariate_matrix$peptide_length[is.na(covariate_matrix$peptide_length)] <- 0
    qvalues <- lm_qvalue(return_df$pvalue, X = covariate_matrix)
    return_df$covariate_qval <- qvalues$qvalues
    return(return_df)    
})


t(do.call(rbind, lapply(unimputed_effect_df, function(x){
    c("pval<0.05" = nrow(x[x$pvalue < 0.05,]),
        "fdr<0.05" = nrow(x[x$fdr < 0.05,]), 
        "qval<0.05" = nrow(x[x$qval < 0.05,]),
        "cov_qval < 0.05" = nrow(x[x$covariate_qval < 0.05,]),
        "fdr_prop > 0.3" = nrow(x[x$fdr_prop > 0.3,]), 
        "qval_prop > 0.3" = nrow(x[x$qval_prop > 0.3,]))
})))
#                 CTRL_ERAP2_vs_CTRL_CTRL VCP_CTRL_vs_CTRL_CTRL
# pval<0.05                           309                  1393
# fdr<0.05                              2                   790
# qval<0.05                             0                     0
# cov_qval < 0.05                       3                  1807

#                 VCP_ERAP2_vs_CTRL_ERAP2 VCP_ERAP2_vs_VCP_CTRL
# pval<0.05                          1128                   414
# fdr<0.05                             81                     1
# qval<0.05                             0                     0
# cov_qval < 0.05                      71                     1


diff_expr_df_with_qval <- diff_expr_df_with_qval[diff_expr_df_with_qval$minprob_q_threshold != 0,]
split_df3 <- split(diff_expr_df_with_qval, diff_expr_df_with_qval$comparison)

mean_effect_df <- lapply(split_df3, function(temp_df){
    logfc_matrix <- reshape2::acast(data = temp_df, formula = protein_name ~ minprob_q_threshold, value.var = "logfc")
    pval_matrix <- reshape2::acast(data = temp_df, formula = protein_name ~ minprob_q_threshold, value.var = "pvalue")
    fdr_matrix <- reshape2::acast(data = temp_df, formula = protein_name ~ minprob_q_threshold, value.var = "fdr")
    qval_matrix <- reshape2::acast(data = temp_df, formula = protein_name ~ minprob_q_threshold, value.var = "qvalue")
    return_df <- data.frame(comparison = unique(temp_df$comparison),
        protein_name = rownames(logfc_matrix),
        logfc = rowMeans(logfc_matrix),
        pvalue = rowMeans(pval_matrix),
        num_pval_sig = rowSums(pval_matrix < 0.05),
        pvalue_prop = rowSums(pval_matrix < 0.05) / ncol(pval_matrix),
        fdr = rowMeans(fdr_matrix),
        num_fdr_sig = rowSums(pval_matrix < 0.05),
        fdr_prop = rowSums(fdr_matrix < 0.05) / ncol(fdr_matrix),
        # qval = rowMeans(qval_matrix),
        # num_qval_sig = rowSums(qval_matrix < 0.05),
        # qval_prop = rowSums(qval_matrix < 0.05) / ncol(qval_matrix),
        sparsity = temp_df$sparsity[match(rownames(logfc_matrix), temp_df$protein_name)],
        average_signal = temp_df$average_signal[match(rownames(logfc_matrix), temp_df$protein_name)],
        peptide_length = temp_df$peptide_length[match(rownames(logfc_matrix), temp_df$protein_name)])
    covariate_matrix <- data.frame(sparsity = return_df$sparsity,
        average_signal = return_df$average_signal / max(return_df$average_signal),
        peptide_length = return_df$peptide_length)
    covariate_matrix$peptide_length[is.na(covariate_matrix$peptide_length)] <- 0
    qvalues <- lm_qvalue(return_df$pvalue, X = covariate_matrix)
    return_df$covariate_qval <- qvalues$qvalues
    return(return_df)
})


t(do.call(rbind, lapply(mean_effect_df, function(x){
    c("pval<0.05" = nrow(x[x$pvalue < 0.05,]),
        "fdr<0.05" = nrow(x[x$fdr < 0.05,]), 
        "qval<0.05" = nrow(x[x$qval < 0.05,]),
        "cov_qval < 0.05" = nrow(x[x$covariate_qval < 0.05,]),
        "fdr_prop > 0.3" = nrow(x[x$fdr_prop > 0.3,]), 
        "qval_prop > 0.3" = nrow(x[x$qval_prop > 0.3,]))
})))

#                 CTRL_ERAP2_vs_CTRL_CTRL VCP_CTRL_vs_CTRL_CTRL
# pval<0.05                            35                  1209
# fdr<0.05                              0                    94
# qval<0.05                             0                     0
# cov_qval < 0.05                       0                  1921
# fdr_prop > 0.3                        3                   154
# qval_prop > 0.3                       0                     0
#                 VCP_ERAP2_vs_CTRL_ERAP2 VCP_ERAP2_vs_VCP_CTRL
# pval<0.05                           859                   348
# fdr<0.05                             10                     0
# qval<0.05                             0                     0
# cov_qval < 0.05                     812                   207
# fdr_prop > 0.3                       32                     1
# qval_prop > 0.3                       0                     0


minprob_default_imputed_data_df <- diff_expr_df_with_qval[diff_expr_df_with_qval$minprob_q_threshold == 0.01,]
split_df4 <- split(minprob_default_imputed_data_df, minprob_default_imputed_data_df$comparison)

minprob_default_imputed_effect_df <- lapply(split_df4, function(temp_df){
    return_df <- data.frame(comparison = unique(temp_df$comparison),
        protein_name = temp_df$protein_name,
        logfc = temp_df$logfc,
        pvalue = temp_df$pvalue,
        num_pval_sig = temp_df$pvalue < 0.05,
        pvalue_prop = temp_df$pvalue < 0.05,
        fdr = temp_df$fdr,
        num_fdr_sig = temp_df$fdr < 0.05,
        fdr_prop = temp_df$fdr < 0.05,
        # qval = rowMeans(qval_matrix),
        # num_qval_sig = rowSums(qval_matrix < 0.05),
        # qval_prop = rowSums(qval_matrix < 0.05) / ncol(qval_matrix),
        sparsity = temp_df$sparsity,
        average_signal = temp_df$average_signal,
        peptide_length = temp_df$peptide_length)
    covariate_matrix <- data.frame(sparsity = return_df$sparsity,
        average_signal = return_df$average_signal / max(return_df$average_signal),
        peptide_length = return_df$peptide_length)
    covariate_matrix$peptide_length[is.na(covariate_matrix$peptide_length)] <- 0
    qvalues <- lm_qvalue(return_df$pvalue, X = covariate_matrix)
    return_df$covariate_qval <- qvalues$qvalues
    return(return_df)    
})

t(do.call(rbind, lapply(minprob_default_imputed_effect_df, function(x){
    c("pval<0.05" = nrow(x[x$pvalue < 0.05,]),
        "fdr<0.05" = nrow(x[x$fdr < 0.05,]), 
        "qval<0.05" = nrow(x[x$qval < 0.05,]),
        "cov_qval < 0.05" = nrow(x[x$covariate_qval < 0.05,]),
        "fdr_prop > 0.3" = nrow(x[x$fdr_prop > 0.3,]), 
        "qval_prop > 0.3" = nrow(x[x$qval_prop > 0.3,]))
})))

#                 CTRL_ERAP2_vs_CTRL_CTRL VCP_CTRL_vs_CTRL_CTRL
# pval<0.05                           306                  1407
# fdr<0.05                              1                   792
# qval<0.05                             0                     0
# cov_qval < 0.05                       0                   874
# fdr_prop > 0.3                        1                   792
# qval_prop > 0.3                       0                     0
#                 VCP_ERAP2_vs_CTRL_ERAP2 VCP_ERAP2_vs_VCP_CTRL
# pval<0.05                           966                   425
# fdr<0.05                             71                     1
# qval<0.05                             0                     0
# cov_qval < 0.05                     259                     0
# fdr_prop > 0.3                       71                     1
# qval_prop > 0.3                       0                     0



saveRDS(object = list(imputed_data = mean_effect_df, 
        unimputed_data = unimputed_data_df, 
        minprob_default_imputed_data = minprob_default_imputed_effect_df), 
    file = file.path(analysis_save_dir, "02-02-rots-differential-enrichment-analysis-results.rds"))
