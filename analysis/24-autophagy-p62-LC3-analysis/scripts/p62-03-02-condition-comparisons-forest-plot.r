data_dir <- normalizePath("/data/shared/lab_internal/")
project_dir <- file.path(data_dir, "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "99-02-autophagy-erap2-kd")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")

all_models_v1 <- readRDS(file.path(analysis_save_dir, "p62/models-v1-all-nbd2.rds"))

# ==========================================================================
# Imports
# ==========================================================================

library("stringr")
library("ggplot2")
library("lme4")
library("dplyr")
library("tidyr")
library("ggbeeswarm")

# ==========================================================================
# Obtain contrast significance values
# ==========================================================================

m_v1_06 <- all_models_v1$models[[6]]

m_v1_06_adj_means <- lapply(m_v1_06, function(model_group) {
  emmeans::emmeans(
    model_group,
    ~ bafa1_timepoint | treatment * condition,
    at = list(
      treatment = c("Scramble", "ERAP2"),
      condition = c("control", "vcp"),
      bafa1_timepoint = c("0h", "4h", "8h")
    ),
    type = "response"
  )
})

m_v1_06_contrast_bafa1_timepoints <- lapply(m_v1_06_adj_means, function(emm) {
  emmeans::contrast(
    emm,
    method = list(
      "4h / 0h" = c(-1,  1, 0),
      "8h / 0h" = c(-1,  0, 1),
      "8h / 4h" = c( 0, -1, 1)
    ),
    by = c("treatment", "condition"),
    adjust = "tukey"
  )
})

m_v1_06_contrast_treatments <- lapply(
  m_v1_06_contrast_bafa1_timepoints,
  function(time_change) {
    emmeans::contrast(
      time_change,
      method = list(
        "ERAP2 change / Scramble change" = c(-1, 1)
      ),
      by = c("condition", "contrast"),
      adjust = "none"
    )
  }
)

m_v1_06_contrast_conditions <- lapply(
  m_v1_06_contrast_treatments,
  function(treatment_change) {
    emmeans::contrast(
      treatment_change,
      method = list(
        "VCP effect / control effect" = c(-1, 1)
      ),
      by = "contrast",
      adjust = "none"
    )
  }
)

response_names <- names(m_v1_06_contrast_conditions)

if (is.null(response_names)) {
  response_names <- paste0(
    "Area bin ",
    seq_along(m_v1_06_contrast_conditions)
  )
}

dod_df <- dplyr::bind_rows(
  lapply(seq_along(m_v1_06_contrast_conditions), function(i) {
    result <- as.data.frame(
      summary(
        m_v1_06_contrast_conditions[[i]],
        infer = c(TRUE, TRUE)
      )
    )

    result$response <- response_names[i]
    result
  })
)

dod_df$p_label <- sprintf("p = %.3g", dod_df$p.value)

dod_df$ratio <- dod_df$ratio - 1
dod_df$asymp.LCL <- dod_df$asymp.LCL - 1
dod_df$asymp.UCL <- dod_df$asymp.UCL - 1

contrast_levels <- c(
  "4h / 0h",
  "8h / 0h",
  "8h / 4h"
)

# Reverse this if you want the current response order flipped
response_levels <- unique(as.character(dod_df$response))

response_spacing <- 1.4
contrast_offsets <- c(-0.25, 0, 0.25)

dod_plot_df <- dod_df |>
  dplyr::mutate(
    response = factor(response, levels = response_levels),
    contrast = factor(contrast, levels = contrast_levels),
    response_centre =
      match(as.character(response), response_levels) *
      response_spacing,
    contrast_offset =
      contrast_offsets[match(as.character(contrast), contrast_levels)],
    y_position = response_centre + contrast_offset
  )

response_breaks <- seq_along(response_levels) * response_spacing

ThePlot <- ggplot(
  dod_plot_df,
  aes(
    x = ratio,
    y = y_position,
    colour = contrast
  )
) +
  geom_vline(
    xintercept = 0,
    colour = "grey60",
    linetype = "dashed"
  ) +
  geom_segment(
    aes(
      x = asymp.LCL,
      xend = asymp.UCL,
      yend = y_position
    ),
    linewidth = 0.7
  ) +
  geom_point(
    aes(shape = contrast),
    size = 2.8
  ) +
  geom_text(
    aes(
      x = asymp.UCL,
      label = p_label
    ),
    hjust = -0.15,
    size = 5,
    show.legend = FALSE
  ) +
  scale_y_continuous(
    breaks = response_breaks,
    labels = response_levels,
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.08, 0.15))
  ) +
  labs(
    x = paste(
      "Difference in treatment-related change:",
      "VCP - control"
    ),
    y = NULL,
    colour = "Timepoint comparison",
    shape = "Timepoint comparison"
  ) +
  theme_bw(base_size = 17) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )

ThePlot

ggsave(
    plot = ThePlot,
    filename = file.path(analysis_plots_dir, "p62", "02-01-plot-p62-pairwise-condition-comparison-pvalues-m_v1_06.pdf"),
    width = 10,
    height = 7,
    units = "in"
)
