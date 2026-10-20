data_dir <- normalizePath("/export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis")
project_dir <- file.path(data_dir, "Projects", "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "11-postmortem-data-analysis/")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
temp_dir <- file.path(analysis_dir, "tmp")
alsod_genes_df <- read.csv(file.path(project_dir, "data_files/alsod_genelist.csv"))
annotation_gtf <- "/nemo/lab/patanir/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
tmp_dir <- file.path(analysis_dir, "slurm")
# genic_count_list <- readRDS(file.path(analysis_save_dir, "01-gene-level-counts-of-postmortem-samples.rds"))
patient_count_df <- readRDS(file.path(analysis_save_dir, "02-summarised-counts-for-patients.rds"))
# ==========================================================================
# Imports
# ==========================================================================
library("stringr")
library("GenomicRanges")
library("DESeq2")
library("rtracklayer")
library("ggplot2")
library("ggbeeswarm")
library("GenomicFeatures")
library("GenomeInfoDb")
require("Rsamtools")
require("STRINGdb")
source("scripts/src/functions.R")
# source("../src/rmats-analysis/src/source.R")
# ==========================================================================
# Analysis
# ==========================================================================

colours <- c("ctrl" = "#EFF7F6", "als" = "#FF1F4B", "ftd" = "#FFA9A3", "als+ftd" = "#F87575")
patient_count_df$condition <- factor(patient_count_df$condition, levels = c("ctrl",  "ftd", "als", "als+ftd"))
patient_count_df$tissue <- factor(patient_count_df$tissue, levels = c("Spinal_Cord_Lumbar", "Spinal_Cord_Thoracic", "Spinal_Cord_Cervical", "Cortex", "Cortex_Motor", "Cerebellum"))


expr_zscore_df <- do.call(rbind, lapply(split(patient_count_df, patient_count_df$tissue),function(temp_df){
    temp_df$expr_zscore <- (temp_df$expr - mean(temp_df$expr)) / sd(temp_df$expr)
    return(temp_df)
}))

expr_zscore_df$tissue_group <- "Brain"
expr_zscore_df$tissue_group[expr_zscore_df$tissue %in% c("Spinal_Cord_Lumbar", "Spinal_Cord_Thoracic", "Spinal_Cord_Cervical")] <- "Spinal_cord"
expr_zscore_df$super_group <- as.character(expr_zscore_df$type)
expr_zscore_df$super_group[expr_zscore_df$super_group %in% "other_mutation-fALS" | expr_zscore_df$super_group %in% "c9orf72-ALS"] <- "known_mutations+fALS"
expr_zscore_df$super_group <- factor(expr_zscore_df$super_group, levels = c("control", "known_mutations+fALS", "sporadic"))


# ThePlot <- ggplot(patient_count_df, aes(x = log10(expr)))
# # ThePlot <- ThePlot + geom_violin(aes(fill = als_group), colour = "#000000", 
# #     lwd = 0.5, position = position_dodge(width = 0.9), scale = "width")
# ThePlot <- ThePlot + geom_density(aes(fill = condition), colour = "#000000", alpha = 0.6)
# ThePlot <- ThePlot + geom_hline(aes(yintercept = 1), colour = "#000000", lty = "dashed")
# ThePlot <- ThePlot + scale_x_continuous("ERAP2 TPM")
# ThePlot <- ThePlot + scale_y_continuous("density")
# ThePlot <- ThePlot + scale_fill_manual("conditions", values = colours)
# ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + ggtitle("erap2 expression violins")
# ThePlot <- ThePlot + facet_wrap(. ~ tissue, ncol = 6, nrow = 1, scales = "free")
# ThePlot <- ThePlot + theme(
#     legend.position = "bottom",
#     text = element_text(size = 12),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     legend.text = element_text(size = 10),
#     legend.title = element_text(size = 10),
#     plot.background = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(colour = "azure3",
#     linetype = "dashed",
#     size = 0.2))
# ggsave(width = 15, height = 5, file = file.path(analysis_plots_dir, "01-00-plot-density-plot-of-erap2-expression-tpm.pdf"), ThePlot)

# boo <- unique(patient_count_df[,c("subject_id", "type", "condition")])
# table(boo$type, boo$condition)



ThePlot <- ggplot(expr_zscore_df, aes(x = condition, y = expr_zscore))
# ThePlot <- ThePlot + geom_violin(aes(fill = als_group), colour = "#000000", 
#     lwd = 0.5, position = position_dodge(width = 0.9), scale = "width")
ThePlot <- ThePlot + geom_boxplot(fill = "#F5EFED", colour = "#000000", 
    width = 0.5, lwd = 0.5, position = position_dodge(width = 0.9), outlier.shape = NA)
ThePlot <- ThePlot + geom_quasirandom(aes(fill = condition), shape = 21, colour = "#000000")
ThePlot <- ThePlot + geom_hline(aes(yintercept = 1), colour = "#000000", lty = "dashed")
ThePlot <- ThePlot + scale_x_discrete("condition")
ThePlot <- ThePlot + scale_y_continuous("ERAP2 expression tpm")
ThePlot <- ThePlot + scale_fill_manual("conditions", values = colours)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("erap2 expression violins")
ThePlot <- ThePlot + facet_wrap(. ~ tissue, ncol = 6, nrow = 1, scales = "free")
ThePlot <- ThePlot + theme(
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
ggsave(width = 15, height = 7, file = file.path(analysis_plots_dir, "01-01-plot-density-of-erap2-expression-tpm.pdf"), ThePlot)



rownames(patient_count_df) <- NULL
patient_count_df$tissue <- droplevels(patient_count_df$tissue)
wilcox_pval_df <- do.call(rbind, lapply(split(patient_count_df, patient_count_df$tissue), function(temp_df){
    als_pval <- wilcox.test(temp_df$expr[temp_df$condition == "als"], temp_df$expr[temp_df$condition == "ctrl"])$p.value
    ftd_pval <- wilcox.test(temp_df$expr[temp_df$condition == "als+ftd"], temp_df$expr[temp_df$condition == "ctrl"])$p.value
    ftd_only_pval <- NA
    if("ftd" %in% temp_df$condition){
        message(unique(temp_df$tissue))
        ftd_only_pval <- wilcox.test(temp_df$expr[temp_df$condition == "ftd"], temp_df$expr[temp_df$condition == "als"])$p.value
    }
    # ftd_vs_als_pval <- wilcox.test(temp_df$expr[temp_df$als_group == "als+ftd"], temp_df$expr[temp_df$als_group == "als"])
    data.frame(tissue = rep(unique(temp_df$tissue),3), als_group = c("als", "als+ftd", "ftd"), pval = c(als_pval, ftd_pval, ftd_only_pval))
}))
rownames(wilcox_pval_df) <- NULL

#                  tissue als_group        pval
# 1    Spinal_Cord_Lumbar       als 0.003452384
# 2    Spinal_Cord_Lumbar   als+ftd 0.073338136
# 3  Spinal_Cord_Thoracic       als 0.118192480
# 4  Spinal_Cord_Thoracic   als+ftd 0.953046953
# 5  Spinal_Cord_Cervical       als 0.057973023
# 6  Spinal_Cord_Cervical   als+ftd 0.506705846
# 7                Cortex       als 0.159612216
# 8                Cortex   als+ftd 0.022310162
# 9          Cortex_Motor       als 0.576438710
# 10         Cortex_Motor   als+ftd 0.643884441
# 11           Cerebellum       als 0.016409862
# 12           Cerebellum   als+ftd 0.045207345




colours <- c("ctrl" = "#EFF7F6", "als" = "#FF1F4B", "ftd" = "#FFA9A3", "als+ftd" = "#F87575")
genotype_colours <- c("XX" = "#CAA8F5", "XY" = "#F1D302")
ThePlot <- ggplot(expr_zscore_df, aes(x = sex_genotype, y = expr_zscore))
ThePlot <- ThePlot + geom_quasirandom(aes(fill = sex_genotype), shape = 21, colour = "#000000")
ThePlot <- ThePlot + geom_boxplot(aes(fill = sex_genotype), colour = "#000000", 
    width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9), outlier.shape = NA)
ThePlot <- ThePlot + scale_x_discrete("mutation")
ThePlot <- ThePlot + scale_y_continuous("Within tissue standardised ERAP2 expression TPM z-score")
ThePlot <- ThePlot + scale_fill_manual("condition", values = genotype_colours)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("erap2 expression violins")
ThePlot <- ThePlot + facet_grid(tissue_group ~ super_group)
ThePlot <- ThePlot + theme(
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
ggsave(width = 7, height = 7, file = file.path(analysis_plots_dir, "01-01-plot-erap2-expression-tpm-in-mutations.pdf"), ThePlot)



rownames(expr_zscore_df) <- NULL
wilcox_pval_df <- do.call(rbind, lapply(split(expr_zscore_df, paste(expr_zscore_df$tissue_group, expr_zscore_df$super_group)), function(temp_df){
    message(unique(temp_df$tissue_group, temp_df$super_group))
    geno_pval <- wilcox.test(temp_df$expr_zscore[temp_df$sex_genotype == "XY"], temp_df$expr_zscore[temp_df$sex_genotype == "XX"])$p.value
    data.frame(tissue_group = unique(temp_df$tissue_group), 
        super_group = unique(temp_df$super_group), 
        genotype = "XY", pval = geno_pval)
}))
rownames(wilcox_pval_df) <- NULL



all_cdf_df <- do.call(rbind, lapply(split(patient_count_df, patient_count_df$tissue),function(temp_df){
    quantiles <- quantile(log10(temp_df$expr), seq(0,1,0.01)) 
    tissue_df <- do.call(rbind, lapply(unique(temp_df$condition), function(a_condition){
        current_cond_df <- temp_df[temp_df$condition == a_condition,]
        cdf_df <- do.call(rbind, lapply(quantiles, function(x){
            data.frame(expr = x, 
                prop = nrow(current_cond_df[log10(current_cond_df$expr) <= x,]) / nrow(current_cond_df))
        }))
        cdf_df$condition <- a_condition
        return(cdf_df)
    }))
    tissue_df$tissue <- unique(temp_df$tissue)
    return(tissue_df)
}))
rownames(all_cdf_df) <- NULL


ctrl_thresholds_df <- do.call(rbind, lapply(split(all_cdf_df, paste(all_cdf_df$tissue, all_cdf_df$condition)), function(x){
    if(!("ctrl" %in% x$condition)){
        return(NULL)
    }
    x[min(which(x$prop >= 0.5)),]
}))

all_thresholds_df <- do.call(rbind, lapply(split(all_cdf_df, paste(all_cdf_df$tissue, all_cdf_df$condition)), function(x){
    x[min(which(x$prop >= 0.5)),]
}))

perc_change_df <- do.call(rbind, lapply(split(all_thresholds_df, all_thresholds_df$tissue), function(x){
    x$linear_expr <- 10^x$expr
    x$perc_change <- x$linear_expr - x$linear_expr[x$condition == "ctrl"]
    return(x)
}))

rownames(perc_change_df) <- NULL


colours <- c("ctrl" = "#000000", "als" = "#FF1F4B", "ftd" = "#FFA9A3", "als+ftd" = "#F87575")
ThePlot <- ggplot(all_cdf_df, aes(x = expr, y = prop))
ThePlot <- ThePlot + geom_line(aes(color = condition))
ThePlot <- ThePlot + geom_vline(data = ctrl_thresholds_df, aes(xintercept = expr), lty = "dashed")
ThePlot <- ThePlot + geom_hline(data = ctrl_thresholds_df, aes(yintercept = prop), lty = "dashed")
ThePlot <- ThePlot + scale_x_continuous("log10 ERAP2 expr")
ThePlot <- ThePlot + scale_y_continuous("CDF")
ThePlot <- ThePlot + scale_color_manual("conditions", values = colours)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("CDF of erap2 expression")
ThePlot <- ThePlot + facet_grid(. ~ tissue, scales = "free")
ThePlot <- ThePlot + theme(
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
ggsave(width = 15, height = 4, file = file.path(analysis_plots_dir, "01-02-plot-CDF-of-erap2-expression-in-tissues.pdf"), ThePlot)


ks_test_results <- do.call(rbind, lapply(split(all_cdf_df, all_cdf_df$tissue), function(temp_df){
    current_conditions <- unique(temp_df$condition)
    seq_condition <- seq_along(current_conditions)
    condition_pairs <- cbind(rep(seq_condition, times = length(seq_condition)), 
        rep(seq_condition, each = length(seq_condition)))
    # condition_pairs <- condition_pairs[condition_pairs[,2] != condition_pairs[,1],]
    pval_df <- do.call(rbind, lapply(seq_len(nrow(condition_pairs)), function(y){
        pval <- ks.test(temp_df[temp_df$condition == current_conditions[condition_pairs[y,1]],"prop"], temp_df[temp_df$condition == current_conditions[condition_pairs[y,2]],"prop"])$p.value
        data.frame(x = current_conditions[condition_pairs[y,1]], y = current_conditions[condition_pairs[y,2]], pval = pval)
    }))
    pval_df$tissue <- unique(temp_df$tissue)
    return(pval_df)
}))
ks_test_results$sig <- ifelse(ks_test_results$pval < 0.001, "***", ifelse(ks_test_results$pval < 0.01, "**", ifelse(ks_test_results$pval < 0.05, "*", "NS")))


ThePlot <- ggplot(ks_test_results, aes(x = x, y = y))
ThePlot <- ThePlot + geom_tile(aes(fill = -log10(pval)))
ThePlot <- ThePlot + geom_text(aes(label = sig))
ThePlot <- ThePlot + scale_x_discrete("conditions")
ThePlot <- ThePlot + scale_y_discrete("conditions")
ThePlot <- ThePlot + scale_fill_distiller("ks test significance", palette = "Reds", direction = 1)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("CDF of erap2 expression")
ThePlot <- ThePlot + facet_wrap(. ~ tissue, scales = "free")
ThePlot <- ThePlot + theme(
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
ggsave(width = 15, height = 4, file = file.path(analysis_plots_dir, "01-02-plot-ks-test-significance-of-erap2-expression-cdfs-in-tissues.pdf"), ThePlot)



# selected_patient_df <- do.call(rbind, lapply(split(patient_count_df, patient_count_df$tissue),function(temp_df){
#     current_threshold <- ctrl_thresholds_df[ctrl_thresholds_df$tissue == unique(temp_df$tissue), "expr"]
#     patients_df <- temp_df[temp_df$condition != "ctrl",]
#     patients_df$erap2_expr_group <- ifelse(log10(patients_df$expr) < current_threshold, "erap2_low", "erap2_high")
#     return(patients_df)
# }))

# selected_patient_df$erap2_expr_group <- factor(selected_patient_df$erap2_expr_group, levels = c("erap2_low", "erap2_high"))

# colours <- c("ctrl" = "#EFF7F6", "als" = "#FF1F4B", "ftd" = "#FFA9A3", "als+ftd" = "#F87575")
# ThePlot <- ggplot(selected_patient_df, aes(x = log10(as.numeric(age_at_symptom_onset)), y = log10(as.numeric(disease_duration_in_months)/12)))
# ThePlot <- ThePlot + geom_point(aes(fill = erap2_expr_group), shape = 21, colour = "#000000")
# ThePlot <- ThePlot + stat_smooth(method = "lm", formula = y ~ x)
# ThePlot <- ThePlot + scale_x_continuous("age at symptom onset")
# ThePlot <- ThePlot + scale_y_continuous("ERAP2 expression tpm")
# ThePlot <- ThePlot + scale_fill_brewer("mutation groups", values = "Set3")
# ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + ggtitle("erap2 expression ~ age at death")
# ThePlot <- ThePlot + facet_grid(. ~ tissue, scales = "free")
# ThePlot <- ThePlot + theme(
#     legend.position = "bottom",
#     text = element_text(size = 12),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     legend.text = element_text(size = 10),
#     legend.title = element_text(size = 10),
#     plot.background = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(colour = "azure3",
#     linetype = "dashed",
#     size = 0.2))
# ggsave(width = 15, height = 5, file = file.path(analysis_plots_dir, "01-03-plot-age-at-symptom-onset-vs-disease-duration.pdf"), ThePlot)


# colours <- c("ctrl" = "#EFF7F6", "als" = "#FF1F4B", "ftd" = "#FFA9A3", "als+ftd" = "#F87575")
# erap2_colours <- c("erap2_low" = "#0496FF","erap2_high" = "#D81159")
# ThePlot <- ggplot(selected_patient_df, aes(x = log10(as.numeric(age_at_symptom_onset)), y = log10(as.numeric(disease_duration_in_months)/12), group = erap2_expr_group))
# ThePlot <- ThePlot + geom_point(aes(fill = erap2_expr_group), shape = 21, colour = "#000000")
# ThePlot <- ThePlot + stat_smooth(data = selected_patient_df, aes(colour = erap2_expr_group), method = "lm", formula = y ~ x)
# ThePlot <- ThePlot + scale_x_continuous("log10 age at symptom onset")
# ThePlot <- ThePlot + scale_y_continuous("log10 disease duration in years")
# ThePlot <- ThePlot + scale_fill_manual("ERAP2 expression", values = erap2_colours)
# ThePlot <- ThePlot + scale_color_manual("ERAP2 expression", values = erap2_colours)
# ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + ggtitle("erap2 expression ~ age at death")
# ThePlot <- ThePlot + facet_grid(. ~ tissue, scales = "free")
# ThePlot <- ThePlot + theme(
#     legend.position = "bottom",
#     text = element_text(size = 12),
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     legend.text = element_text(size = 10),
#     legend.title = element_text(size = 10),
#     plot.background = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(colour = "azure3",
#     linetype = "dashed",
#     size = 0.2))
# ggsave(width = 15, height = 5, file = file.path(analysis_plots_dir, "01-03-plot-disease-duration-associated-erap2-expression-tpm.pdf"), ThePlot)


# lm_list <- lapply(split(selected_patient_df, paste(selected_patient_df$tissue)), function(temp_df){
#     current_df <- temp_df
#     current_df$duration <- as.numeric(current_df$disease_duration_in_months)/12
#     summary(lm(data = current_df, 
#         formula = log10(duration) ~ log10(as.numeric(age_at_symptom_onset)) + erap2_expr_group:log10(as.numeric(age_at_symptom_onset))))
# })



# boo_df <- selected_patient_df[selected_patient_df$tissue == "Cortex",]
# boo_df$duration <- as.numeric(boo_df$disease_duration_in_months)/12


# dist_fit <- fitdistr(x = boo_df$duration[!is.na(boo_df$duration)], densfun = "exponential")
# sim_y <- rexp(n = length(boo_df$duration[!is.na(boo_df$duration)]), rate = dist_fit$estimate[1])

# # ks.test(sim_y, boo_df$duration[!is.na(boo_df$duration)])


# patient_with_multiple_tissues <- table(patient_count_df$subject_id, patient_count_df$tissue)

# summary_df <- summarise_tissue_specific_parameter_by_patients(patient_count_df)
# saveRDS(object = summary_df, file = file.path(analysis_save_dir, "02-postmortem-patient-counts-summarised.rds"))


