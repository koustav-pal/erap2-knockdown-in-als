# Adapted analysis; see provenance/source_ranges.json.
source(here::here("workflow", "scripts", "common.R"))
library(stringr)
library(ggbeeswarm)
patient_count_df <- readRDS(input("patient_expression.rds"))
patient_count_df$tissue_group <- "spinal_cord"
patient_count_df$tissue_group[patient_count_df$tissue %in% c("Cortex", "Cortex_Motor", "Cerebellum")] <- "brain"
patient_count_df$type <- factor(patient_count_df$type, levels = c("control", "c9orf72-ALS", "other_mutation-fALS", "sporadic"))
patient_count_df$type_group <- as.character(patient_count_df$type)
patient_count_df$type_group[patient_count_df$type %in% c("other_mutation-fALS", "sporadic")] <- "sporadic+other_fALS"
patient_count_df$tissue_group <- "spinal_cord"
patient_count_df$tissue_group[patient_count_df$tissue %in% c("Cortex", "Cortex_Motor", "Cerebellum")] <- "brain"

patient_surv_df <- do.call(rbind, lapply(split(patient_count_df, paste(patient_count_df$subject_id, patient_count_df$tissue)), function(temp_df){
    data.frame(subject_id = unique(temp_df$subject_id),
        tissue_group = unique(temp_df$tissue_group),
        age_at_onset = unique(temp_df$age_at_symptom_onset),
        age_at_death = unique(temp_df$age_at_death),
        disease_duration_in_months = unique(temp_df$disease_duration_in_months),
        sex_genotype = unique(temp_df$sex_genotype),
        mnd_with_ftd = unique(temp_df$mnd_with_ftd),
        condition = unique(temp_df$condition),
        mnd_with_dementia = unique(temp_df$mnd_with_dementia),
        type_group = unique(temp_df$type_group),
        expr = mean(temp_df$expr))
}))

expr_zscore_df <- do.call(rbind, lapply(split(patient_surv_df, patient_surv_df$tissue_group), function(temp_df){
    temp_df$expr_zscore <- (temp_df$expr - mean(temp_df$expr)) / sd(temp_df$expr)
    return(temp_df)
}))

expr_zscore_df$type_group <- factor(expr_zscore_df$type_group, levels = c("control", "c9orf72-ALS", "sporadic+other_fALS"))

wilcox.test(expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "brain" & expr_zscore_df$type_group == "control"],
    expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "brain" & expr_zscore_df$type_group == "sporadic+other_fALS"])
wilcox.test(expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "brain" & expr_zscore_df$type_group == "control"],
    expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "brain" & expr_zscore_df$type_group == "c9orf72-ALS"])

wilcox.test(expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "spinal_cord" & expr_zscore_df$type_group == "control"],
    expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "spinal_cord" & expr_zscore_df$type_group == "sporadic+other_fALS"])
wilcox.test(expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "spinal_cord" & expr_zscore_df$type_group == "control"],
    expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "spinal_cord" & expr_zscore_df$type_group == "c9orf72-ALS"])

wilcox.test(expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "brain" & expr_zscore_df$type_group == "c9orf72-ALS"],
    expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "brain" & expr_zscore_df$type_group == "sporadic+other_fALS"])
wilcox.test(expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "spinal_cord" & expr_zscore_df$type_group == "c9orf72-ALS"],
    expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "spinal_cord" & expr_zscore_df$type_group == "sporadic+other_fALS"])

fill_colours <- c("control" = "#F6F4F3",
    "c9orf72-ALS" = "#169873", "sporadic" = "#F49FBC",
    "sporadic+other_fALS" = "#F49FBC",
     "other_mutation-fALS" = "#9EBD6E")

ThePlot <- ggplot(expr_zscore_df, aes(x = type_group, y = expr_zscore))
ThePlot <- ThePlot + geom_quasirandom(aes(fill = type_group), shape = 21, colour = "#000000")
ThePlot <- ThePlot + geom_boxplot(aes(fill = type_group), colour = "#000000",
    width = 0.2, lwd = 0.5, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.6)
ThePlot <- ThePlot + geom_hline(aes(yintercept = 0), colour = "#000000", lty = "dashed")
ThePlot <- ThePlot + scale_x_discrete("mutation status")
ThePlot <- ThePlot + scale_y_continuous("ERAP2 expression z-score")
ThePlot <- ThePlot + scale_fill_manual("conditions", values = fill_colours)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("erap2 expression violins")
ThePlot <- ThePlot + facet_wrap(. ~ tissue_group, ncol = 6, nrow = 1, scales = "free")
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
ggsave(width = 7, height = 7, file = file.path(plot_dir, "01-01-plot-expression-zscore-of-erap2-expression-tpm.pdf"), ThePlot)

wilcox.test(expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "brain" & expr_zscore_df$condition == "ctrl" & expr_zscore_df$sex_genotype == "XX"],
    expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "brain" & expr_zscore_df$condition == "ctrl" & expr_zscore_df$sex_genotype == "XY"])
wilcox.test(expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "brain" & expr_zscore_df$condition != "ctrl" & expr_zscore_df$sex_genotype == "XX"],
    expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "brain" & expr_zscore_df$condition != "ctrl" & expr_zscore_df$sex_genotype == "XY"])

wilcox.test(expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "spinal_cord" & expr_zscore_df$condition == "ctrl" & expr_zscore_df$sex_genotype == "XX"],
    expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "spinal_cord" & expr_zscore_df$condition == "ctrl" & expr_zscore_df$sex_genotype == "XY"])
wilcox.test(expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "spinal_cord" & expr_zscore_df$condition != "ctrl" & expr_zscore_df$sex_genotype == "XX"],
    expr_zscore_df$expr_zscore[expr_zscore_df$tissue_group == "spinal_cord" & expr_zscore_df$condition != "ctrl" & expr_zscore_df$sex_genotype == "XY"])

fill_colours <- c("XX" = "#EF9CDA", "XY" = "#00A5E0")

expr_zscore_df$disease_status <- "control"
expr_zscore_df$disease_status[expr_zscore_df$condition != "ctrl"] <- "ALS/ALS+FTD/FTD"

ThePlot <- ggplot(expr_zscore_df, aes(x = sex_genotype, y = expr_zscore))
ThePlot <- ThePlot + geom_quasirandom(aes(fill = sex_genotype), shape = 21, colour = "#000000")
ThePlot <- ThePlot + geom_boxplot(aes(fill = sex_genotype), colour = "#000000",
    width = 0.2, lwd = 0.5, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.6)
ThePlot <- ThePlot + geom_hline(aes(yintercept = 0), colour = "#000000", lty = "dashed")
ThePlot <- ThePlot + scale_x_discrete("mutation status")
ThePlot <- ThePlot + scale_y_continuous("ERAP2 expression z-score")
ThePlot <- ThePlot + scale_fill_manual("conditions", values = fill_colours)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("erap2 expression violins")
ThePlot <- ThePlot + facet_grid(tissue_group ~ disease_status, scales = "free_y")
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
ggsave(width = 7, height = 7, file = file.path(plot_dir, "01-02-plot-expression-zscore-of-erap2-expression-tpm-stratified-by-genotypic-sex.pdf"), ThePlot)
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

writeLines(capture.output(sessionInfo()), file.path(plot_dir, "sessionInfo.txt"))
