# Adapted analysis; see provenance/source_ranges.json.
source(here::here("workflow", "scripts", "common.R"))
library(lme4)
library(emmeans)
bodipy_cell_level_list <- readRDS(input("bodipy_cell_metrics.rds"))
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

emm <- emmeans::emmeans(
  fit_lme,
  ~ treatment | condition
)

emmeans::contrast(
  emm,
  method = "revpairwise"
)

treatment_means <- emmeans::emmeans(
  fit_lme,
  ~ treatment | condition,
  at = list(
    treatment = c("Scramble", "ERAP2"),
    condition = c("control", "vcp")
  )
)

treatment_responses <- emmeans::contrast(
  treatment_means,
  method = list("ERAP2 - Scramble" = c(-1, 1))
)

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
cell_data$label <- paste(cell_data$condition, cell_data$treatment, sep = "_")
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
ThePlot <- ThePlot + scale_fill_manual("condition", values = condition_colours)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("Violin plot of spot sizes")
ggsave(width = 5, height = 10, file = file.path(analysis_plots_dir, "07-01-plot-covariate-adjusted-mean-spot-size-per-cell.pdf"), ThePlot)

write.csv(as.data.frame(treatment_responses), file.path(plot_dir, "treatment_contrasts.csv"), row.names = FALSE)
saveRDS(fit_lme, file.path(plot_dir, "bodipy_lmer.rds"))

writeLines(capture.output(sessionInfo()), file.path(plot_dir, "sessionInfo.txt"))
