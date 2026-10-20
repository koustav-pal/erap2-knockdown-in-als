data_dir <- normalizePath("/data/shared/lab_internal/")
project_dir <- file.path(data_dir, "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "99-01-bodipy-erap2-kd")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
bodipy_cell_level_list <- readRDS(file.path(analysis_save_dir, "all-CLO-metrics.rds"))
# ==========================================================================
# Imports
# ==========================================================================
library("stringr")
library("ggplot2")
library("lme4")
# source("../src/bodipy_cdfs/source.R")
# source("../src/rmats-analysis/src/source.R")
# ==========================================================================
# Analysis
# ==========================================================================
# 1. Build the membership matrix
# Rows = cells in X
# Columns = membership vectors in K
bodipy_cell_level_list$well_uid <- with(
  bodipy_cell_level_list,
  interaction(experiment, well_id, drop = TRUE)
)

model_formula <- CLO_mean_SLO_spot_area_um2 ~ condition * treatment +
  (1 | cell_line) +
  (1 | experiment) +
  (1 | cell_line:experiment) +
  (1 | well_uid)


fit_lme <- lme4::lmer(
  model_formula,
  data = bodipy_cell_level_list,
  REML = TRUE
)
summary(fit_lme)
# Linear mixed model fit by REML ['lmerMod']
# Formula:
# CLO_mean_SLO_spot_area_um2 ~ condition * treatment + (1 | cell_line) +
#     (1 | experiment) + (1 | cell_line:experiment) + (1 | well_uid)
#    Data: bodipy_cell_level_list

# REML criterion at convergence: 148087

# Scaled residuals:
#     Min      1Q  Median      3Q     Max
# -2.3130 -0.5499 -0.1544  0.2990 16.7956

# Random effects:
#  Groups               Name        Variance Std.Dev.
#  well_uid             (Intercept) 0.008982 0.09478
#  cell_line:experiment (Intercept) 0.102374 0.31996
#  cell_line            (Intercept) 0.004921 0.07015
#  experiment           (Intercept) 0.295086 0.54322
#  Residual                         1.446747 1.20281
# Number of obs: 46096, groups:
# well_uid, 104; cell_line:experiment, 26; cell_line, 8; experiment, 4

# Fixed effects:
#                             Estimate Std. Error t value
# (Intercept)                  3.13736    0.29059  10.796
# conditionvcp                -0.05626    0.14155  -0.397
# treatmentERAP2               0.02818    0.03202   0.880
# conditionvcp:treatmentERAP2  0.03521    0.04393   0.801

# Correlation of Fixed Effects:
#             (Intr) cndtnv tERAP2
# conditinvcp -0.259
# trtmntERAP2 -0.054  0.110
# cndtn:ERAP2  0.039 -0.150 -0.729

emm <- emmeans::emmeans(
  fit_lme,
  ~ treatment | condition
)

# Check treatment-level order before interpreting the contrast sign.
emmeans::contrast(
  emm,
  method = "revpairwise"
)
# condition = control:
#  contrast         estimate     SE  df z.ratio p.value
#  ERAP2 - Scramble   0.0282 0.0320 Inf   0.880  0.3789

# condition = vcp:
#  contrast         estimate     SE  df z.ratio p.value
#  ERAP2 - Scramble   0.0634 0.0301 Inf   2.107  0.0351

# Degrees-of-freedom method: asymptotic

# predict(fit_lme)

treatment_means <- emmeans::emmeans(
  fit_lme,
  ~ treatment | condition,
  at = list(
    treatment = c("Scramble", "ERAP2"),
    condition = c("control", "vcp")
  )
)

# Treatment response within each condition
treatment_responses <- emmeans::contrast(
  treatment_means,
  method = list("ERAP2 - Scramble" = c(-1, 1))
)

# Compare the two responses
interaction_comparison <- emmeans::contrast(
  treatment_responses,
  method = list("VCP response - control response" = c(-1, 1)),
  by = NULL
)

summary(
  interaction_comparison,
  infer = c(TRUE, TRUE),
  adjust = "none"
)

cell_data <- model.frame(fit_lme)

cell_data$adjusted <- model.response(cell_data) -
  predict(fit_lme, re.form = NULL) +
  predict(fit_lme, re.form = NA)

cell_data_split <- split(cell_data, cell_data$well_uid)
well_average_df <- do.call(rbind, lapply(cell_data_split, function(temp_df){
    return_df <- unique(cell_data[,c("well_uid", "condition", "treatment", "cell_line", "experiment", "label")])
    return_df$adjusted <- mean(temp_df$adjusted)
    return(return_df)
}))


# group_fitted <- emmeans::emmeans(
#   fit_lme,
#   ~ treatment | condition
# )

# as.data.frame(confint(group_fitted))

# condition = control:
#  treatment   emmean        SE  df asymp.LCL asymp.UCL
#  Scramble  3.137363 0.2905911 Inf  2.567815  3.706911
#  ERAP2     3.165541 0.2906391 Inf  2.595899  3.735183

# condition = vcp:
#  treatment   emmean        SE  df asymp.LCL asymp.UCL
#  Scramble  3.081099 0.2883304 Inf  2.515982  3.646217
#  ERAP2     3.144487 0.2883817 Inf  2.579269  3.709705
# Degrees-of-freedom method: asymptotic
# Confidence level used: 0.95


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

condition_colours <- c("control_Scramble" = "#a6cee3", "control_ERAP2" = "#1f78b4", "vcp_Scramble" = "#fb9a99", "vcp_ERAP2" =  "#e31a1c")
# condition_colours <- c("control" = "#1f78b4", "vcp" = "#e31a1c")
cell_data$label <- paste(cell_data$condition, cell_data$treatment, sep = "_")
# average_spot_size_of_cell_lines_df$label <- factor(average_spot_size_of_cell_lines_df$label, levels = names(condition_colours))
cell_data$treatment <- factor(cell_data$treatment, levels = c("Scramble", "ERAP2"))
cell_data <- cell_data[order(cell_data$treatment, decreasing = T),]
cell_data$condition <- factor(cell_data$condition, levels = c("control", "vcp"))
cell_data$label <- factor(cell_data$label, levels = c("control_Scramble", "control_ERAP2", "vcp_Scramble", "vcp_ERAP2"))

ThePlot <- ggplot(cell_data, aes(x = condition, y = log10(adjusted)))
ThePlot <- ThePlot + geom_violin(aes(fill = label), colour = "#000000", position = position_dodge(width = 0.9))
ThePlot <- ThePlot + geom_boxplot(aes(fill = label), colour = "#000000", 
    width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9))
ThePlot <- ThePlot + scale_x_discrete("condition")
ThePlot <- ThePlot + scale_y_continuous("log10 covariate adjusted mean spot size")
# ThePlot <- ThePlot + scale_colour_manual("condition", values = condition_colours)
ThePlot <- ThePlot + scale_fill_manual("condition", values = condition_colours)
# ThePlot <- ThePlot + scale_colour_manual("status", values = Colours)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ treatment, scales = "free")
ThePlot <- ThePlot + ggtitle("Violin plot of spot sizes")
ggsave(width = 5, height = 10, file = file.path(analysis_plots_dir, "07-01-plot-covariate-adjusted-mean-spot-size-per-cell.pdf"), ThePlot)

