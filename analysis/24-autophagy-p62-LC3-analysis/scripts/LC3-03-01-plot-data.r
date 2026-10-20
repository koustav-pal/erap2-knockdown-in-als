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
LC3_cell_level_list <- readRDS(file.path(analysis_save_dir, "LC3/all-CLO-metrics-binned.rds"))

# ==========================================================================
# Imports
# ==========================================================================

library("stringr")
library("ggplot2")
library("dplyr")
library("tidyr")
library("ggbeeswarm")
library("ggpubr")

# ==========================================================================
# Analysis - grouped by BafA1 timepoint
# ==========================================================================

best_model_v1 <- readRDS(file.path(analysis_save_dir, "LC3", "models-v1-all-nbd2.rds"))
m_v1_06 <- best_model_v1$models[[6]]

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
    adjust = "tukey",
    type = "response"
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
      adjust = "none",
      type = "response"
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
      adjust = "none",
      type = "response"
    )
  }
)

# ==========================================================================
# Extract significant values
# ==========================================================================

# List elements correspond to bins 1, 2, ... unless bin_labels is supplied.
extract_significant_p_values <- function(
    contrast_list,
    alpha = 0.05,
    treatment = NULL,
    bin_labels = as.character(seq_along(contrast_list))
) {
    if (!is.list(contrast_list)) {
        stop("contrast_list must be a list of emmeans contrast objects.")
    }
    if (length(alpha) != 1L || !is.numeric(alpha) ||
        !is.finite(alpha) || alpha <= 0 || alpha > 1) {
        stop("alpha must be a number greater than 0 and at most 1.")
    }
    if (length(bin_labels) != length(contrast_list) ||
        anyNA(bin_labels) || anyDuplicated(bin_labels)) {
        stop("Provide one unique, non-missing bin label per list element.")
    }
    if (!is.null(treatment) &&
        (length(treatment) != 1L || is.na(treatment))) {
        stop("Select a single treatment, e.g. treatment = 'ERAP2'.")
    }

    empty_result <- data.frame(
        condition = character(),
        bin = character(),
        group1 = character(),
        group2 = character(),
        p_value = numeric(),
        stringsAsFactors = FALSE
    )
    if (length(contrast_list) == 0L) {
        return(empty_result)
    }
    if (!requireNamespace("emmeans", quietly = TRUE)) {
        stop("The emmeans package is required.")
    }

    results <- lapply(seq_along(contrast_list), function(i) {
        # Summarize all contrasts before filtering to preserve test families.
        tab <- as.data.frame(summary(
            contrast_list[[i]],
            infer = c(FALSE, TRUE)
        ))
        required <- c("condition", "contrast", "p.value")
        if (!all(required %in% names(tab))) {
            stop("List element ", i, " must contain condition, contrast, and p.value.")
        }

        if ("treatment" %in% names(tab)) {
            available <- unique(as.character(tab$treatment))
            selected <- treatment
            if (is.null(selected)) {
                if (length(available) != 1L) {
                    stop(
                        "Timepoint contrasts contain multiple treatments. ",
                        "Call separately with treatment = 'Scramble' or 'ERAP2'."
                    )
                }
                selected <- available
            }
            if (!selected %in% available) {
                stop("Treatment '", selected, "' is absent from list element ", i, ".")
            }
            tab <- tab[as.character(tab$treatment) %in% selected, , drop = FALSE]
        } else if (!is.null(treatment)) {
            stop("This contrast list has no treatment column; omit treatment.")
        }

        tab <- tab[
            is.finite(tab$p.value) & tab$p.value < alpha,
            ,
            drop = FALSE
        ]
        if (nrow(tab) == 0L) {
            return(empty_result)
        }

        groups <- strsplit(trimws(as.character(tab$contrast)), "\\s+/\\s+", perl = TRUE)
        if (any(lengths(groups) != 2L) ||
            any(vapply(groups, function(x) any(!nzchar(trimws(x))), logical(1)))) {
            stop("Expected timepoint contrast labels such as '8h / 4h'.")
        }

        data.frame(
            condition = as.character(tab$condition),
            bin = rep(as.character(bin_labels[i]), nrow(tab)),
            group1 = vapply(groups, function(x) trimws(x[2L]), character(1)),
            group2 = vapply(groups, function(x) trimws(x[1L]), character(1)),
            p_value = tab$p.value,
            stringsAsFactors = FALSE
        )
    })

    result <- do.call(rbind, results)
    rownames(result) <- NULL
    result
}

# Change this prefix for future models
model_prefix <- "m_v1_06"

treatments_p_values <- extract_significant_p_values(
    get(paste0(model_prefix, "_contrast_treatments")),
    alpha = 0.05
)

# BafA1 contrasts (not depicted in final plot)
timepoint_contrasts <- get(
    paste0(model_prefix, "_contrast_bafa1_timepoints")
)
p_values_ERAP2 <- extract_significant_p_values(
    timepoint_contrasts,
    treatment = "ERAP2"
)
p_values_Scramble <- extract_significant_p_values(
    timepoint_contrasts,
    treatment = "Scramble"
)

combined_bafa1_timepoints_p_values <- rbind(
    cbind(
        p_values_Scramble,
        treatment = rep("Scramble", nrow(p_values_Scramble))
    ),
    cbind(
        p_values_ERAP2,
        treatment = rep("ERAP2", nrow(p_values_ERAP2))
    )
)

combined_bafa1_timepoints_p_values$treatment <- factor(
    combined_bafa1_timepoints_p_values$treatment,
    levels = c("Scramble", "ERAP2")
)
combined_bafa1_timepoints_p_values$condition <- factor(
    combined_bafa1_timepoints_p_values$condition,
    levels = c("control", "vcp")
)
timepoint_levels <- c("0h", "4h", "8h")
for (column in c("group1", "group2")) {
    combined_bafa1_timepoints_p_values[[column]] <- factor(
        combined_bafa1_timepoints_p_values[[column]],
        levels = timepoint_levels
    )
}
rownames(combined_bafa1_timepoints_p_values) <- NULL

# ==========================================================================
# Create model frame for covariate-adjusted data
# ==========================================================================

m_v1_06_model_frame <- lapply(m_v1_06, function(model_group) {
  mf <- model.frame(model_group)

  mf$data_adj <-
    model.response(mf) -
    predict(model_group, type = "response", re.form = NULL) +
    predict(model_group, type = "response", re.form = NA)

  mf
})

condition_colours <- c("control_Scramble" = "#a6cee3", "control_ERAP2" = "#1f78b4", "vcp_Scramble" = "#fb9a99", "vcp_ERAP2" =  "#e31a1c")

m_v1_06_model_frame <- lapply(m_v1_06_model_frame, function(bin_df) {
  bin_df$label <- paste(bin_df$condition, bin_df$treatment, sep = "_")
  bin_df$label <- factor(bin_df$label, levels = names(condition_colours))
  bin_df
})

m_v1_06_model_frame_combined <- m_v1_06_model_frame[[1]]

for (x in seq_along(m_v1_06_model_frame)) {
  col_name <- paste0("adj_num_spots_bin_", x)

  m_v1_06_model_frame_combined[[col_name]] <- m_v1_06_model_frame[[x]]$data_adj
}

m_v1_06_model_frame_combined$data_adj <- NULL

# ==========================================================================
# Plot final data
# ==========================================================================

plot_combined_bin_v3 <- function(
    model_frame_df,
    treatments_p_values = NULL,
    file_name = NA_character_,
    timepoint_levels = c("0h", "4h", "8h"),
    p_label_digits = 3,
    combined_bafa1_timepoints_p_values = NULL,
    dodge_width = 1,
    box_width = 0.75,
    y_limits = NULL,
    plot_title = "LC3 spots",
    output_dir = file.path(analysis_plots_dir, "LC3"),
    annotation_step = NULL,
    annotation_gap = NULL
) {
    metric_cols <- paste0("adj_num_spots_bin_", seq_len(4))
    bin_labels <- as.character(seq_len(4))
    condition_levels <- c("control", "vcp")
    treatment_levels <- c("Scramble", "ERAP2")
    condition_colours <- c(
        "control_Scramble" = "#a6cee3",
        "control_ERAP2" = "#1f78b4",
        "vcp_Scramble" = "#fb9a99",
        "vcp_ERAP2" = "#e31a1c"
    )
    label_levels <- names(condition_colours)
    legend_labels <- c("control / Scr", "control / ERAP2", "vcp / Scr", "vcp / ERAP2")
    required_cols <- c(metric_cols, "condition", "treatment", "bafa1_timepoint")
    missing_cols <- setdiff(required_cols, names(model_frame_df))
    if (length(missing_cols) > 0) {
        stop("model_frame_df is missing: ", paste(missing_cols, collapse = ", "))
    }
    if (length(timepoint_levels) == 0 || anyNA(timepoint_levels) ||
        anyDuplicated(timepoint_levels)) {
        stop("timepoint_levels must contain unique, non-missing timepoints.")
    }
    if (length(dodge_width) != 1L || !is.numeric(dodge_width) ||
        !is.finite(dodge_width) || dodge_width <= 0) {
        stop("dodge_width must be a finite, positive number.")
    }
    if (length(box_width) != 1L || !is.numeric(box_width) ||
        !is.finite(box_width) || box_width <= 0) {
        stop("box_width must be a finite, positive number.")
    }
    if (!is.null(y_limits) &&
        (length(y_limits) != 2L || !is.numeric(y_limits) ||
         any(!is.na(y_limits) & !is.finite(y_limits)))) {
        stop("y_limits must be NULL or two numbers; NA selects an automatic endpoint.")
    }
    for (argument in c("annotation_step", "annotation_gap")) {
        value <- get(argument, inherits = FALSE)
        if (!is.null(value) &&
            (length(value) != 1L || !is.numeric(value) || !is.finite(value) || value <= 0)) {
            stop(argument, " must be NULL or a finite, positive number.")
        }
    }

    # Match the four centres used by position_dodge(), in legend order.
    label_offsets <- setNames(
        dodge_width * ((seq_along(label_levels) - 0.5) / length(label_levels) - 0.5),
        label_levels
    )
    condition_offsets <- setNames(vapply(condition_levels, function(condition) {
        mean(label_offsets[paste(condition, treatment_levels, sep = "_")])
    }, numeric(1)), condition_levels)
    point_width <- min(0.12, 0.4 * dodge_width / length(label_levels))
    x_padding <- max(
        0.5,
        max(abs(label_offsets)) + max(point_width, box_width / (2 * length(label_levels))) + 0.1
    )
    x_display_limits <- c(1 - x_padding, length(timepoint_levels) + x_padding)

    plot_df <- do.call(rbind, lapply(seq_along(metric_cols), function(i) {
        data.frame(
            condition = as.character(model_frame_df$condition),
            treatment = as.character(model_frame_df$treatment),
            bafa1_timepoint = as.character(model_frame_df$bafa1_timepoint),
            bin = bin_labels[i],
            plot_value = model_frame_df[[metric_cols[i]]],
            stringsAsFactors = FALSE
        )
    }))
    plot_df <- plot_df[
        is.finite(plot_df$plot_value) &
            plot_df$condition %in% condition_levels &
            plot_df$treatment %in% treatment_levels &
            plot_df$bafa1_timepoint %in% timepoint_levels,
        ,
        drop = FALSE
    ]
    if (nrow(plot_df) == 0) {
        stop("No finite rows available for combined bin plots.")
    }

    expected_layout <- expand.grid(
        condition = condition_levels,
        bin = bin_labels,
        bafa1_timepoint = timepoint_levels,
        treatment = treatment_levels,
        stringsAsFactors = FALSE
    )
    expected_keys <- with(expected_layout, paste(condition, bin, bafa1_timepoint, treatment, sep = "|"))
    observed_keys <- unique(with(plot_df, paste(condition, bin, bafa1_timepoint, treatment, sep = "|")))
    missing_layout <- expected_layout[!expected_keys %in% observed_keys, , drop = FALSE]
    if (nrow(missing_layout) > 0) {
        stop(
            "The plot data are missing condition/bin/timepoint/treatment combinations. ",
            "First missing combination: ",
            paste(missing_layout[1, ], collapse = " / ")
        )
    }

    plot_df$condition <- factor(plot_df$condition, levels = condition_levels)
    plot_df$treatment <- factor(plot_df$treatment, levels = treatment_levels)
    plot_df$bin <- factor(plot_df$bin, levels = bin_labels)
    plot_df$bafa1_timepoint <- factor(plot_df$bafa1_timepoint, levels = timepoint_levels)
    plot_df$label <- factor(
        paste(plot_df$condition, plot_df$treatment, sep = "_"),
        levels = label_levels
    )

    prepare_p_values <- function(p_values, source_name, by_treatment = FALSE) {
        if (is.null(p_values)) {
            return(NULL)
        }
        required <- c("condition", "bin", "group1", "group2", "p_value")
        if (by_treatment) {
            required <- c(required, "treatment")
        }
        missing <- setdiff(required, names(p_values))
        if (length(missing) > 0) {
            stop(source_name, " is missing: ", paste(missing, collapse = ", "))
        }
        p_values <- as.data.frame(p_values, stringsAsFactors = FALSE)
        if (nrow(p_values) == 0) {
            return(NULL)
        }
        for (column in c("condition", "bin", "group1", "group2")) {
            p_values[[column]] <- as.character(p_values[[column]])
        }
        p_values$p_value <- as.numeric(as.character(p_values$p_value))
        if (any(!is.finite(p_values$p_value) |
                p_values$p_value < 0 | p_values$p_value > 1)) {
            stop(source_name, ": p_value must be finite and between 0 and 1.")
        }
        if (any(!p_values$condition %in% condition_levels)) {
            stop("Unknown condition in ", source_name, ".")
        }
        if (any(!p_values$bin %in% bin_labels)) {
            stop("Unknown bin in ", source_name, ".")
        }
        if (any(!p_values$group1 %in% timepoint_levels) ||
            any(!p_values$group2 %in% timepoint_levels)) {
            stop("Unknown timepoint in ", source_name, ".")
        }
        if (any(p_values$group1 == p_values$group2)) {
            stop(source_name, ": group1 and group2 must be different timepoints.")
        }
        if (by_treatment) {
            p_values$treatment <- as.character(p_values$treatment)
            if (any(!p_values$treatment %in% treatment_levels)) {
                stop("Unknown treatment in ", source_name, ".")
            }
        } else {
            p_values$treatment <- rep(NA_character_, nrow(p_values))
        }
        group1_index <- match(p_values$group1, timepoint_levels)
        group2_index <- match(p_values$group2, timepoint_levels)
        keys <- data.frame(
            condition = p_values$condition,
            bin = p_values$bin,
            treatment = p_values$treatment,
            left = pmin(group1_index, group2_index),
            right = pmax(group1_index, group2_index)
        )
        if (anyDuplicated(keys)) {
            stop(source_name, " contains a duplicated condition/bin/treatment/comparison.")
        }
        if (!"y_position" %in% names(p_values)) {
            p_values$y_position <- rep(NA_real_, nrow(p_values))
        } else {
            numeric_y <- as.numeric(as.character(p_values$y_position))
            if (any(!is.na(p_values$y_position) & !is.finite(numeric_y))) {
                stop(source_name, ": y_position must be finite or NA for automatic placement.")
            }
            p_values$y_position <- numeric_y
        }
        if (by_treatment) {
            p_values <- p_values[p_values$p_value < 0.05, , drop = FALSE]
            if (nrow(p_values) == 0) {
                return(NULL)
            }
        }

        significance <- as.character(cut(
            p_values$p_value,
            breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
            labels = c("***", "**", "*", "ns"),
            right = FALSE
        ))
        p_label <- significance
        annotation_colour <- rep("black", nrow(p_values))
        x_offset <- unname(condition_offsets[p_values$condition])
        if (by_treatment) {
            label_ids <- paste(p_values$condition, p_values$treatment, sep = "_")
            x_offset <- unname(label_offsets[label_ids])
            annotation_colour <- unname(condition_colours[label_ids])
        } else {
            formatted_p <- trimws(format.pval(
                p_values$p_value,
                digits = p_label_digits,
                eps = 10^(-p_label_digits)
            ))
            p_label <- paste(
                ifelse(startsWith(formatted_p, "<"), "p", "p ="),
                formatted_p,
                significance
            )
        }
        data.frame(
            condition = factor(p_values$condition, levels = condition_levels),
            bin = factor(p_values$bin, levels = bin_labels),
            by_treatment = rep(by_treatment, nrow(p_values)),
            y_position = p_values$y_position,
            p_label = p_label,
            annotation_colour = annotation_colour,
            x_start = match(p_values$group1, timepoint_levels) + x_offset,
            x_end = match(p_values$group2, timepoint_levels) + x_offset,
            stringsAsFactors = FALSE
        )
    }

    annotation_df <- rbind(
        prepare_p_values(
            combined_bafa1_timepoints_p_values,
            "combined_bafa1_timepoints_p_values",
            by_treatment = TRUE
        ),
        prepare_p_values(treatments_p_values, "treatments_p_values")
    )

    current_plot <- ggplot(
        plot_df,
        aes(
            x = bafa1_timepoint,
            y = plot_value,
            fill = label,
            group = interaction(bafa1_timepoint, label, drop = TRUE)
        )
    ) +
    geom_boxplot(
        colour = "black",
        alpha = 1,
        width = box_width,
        linewidth = 0.3,
        outlier.shape = NA,
        position = position_dodge(width = dodge_width)
    ) +
    facet_wrap(~ bin, nrow = 1, drop = FALSE) +
    scale_x_discrete("BafA1 treatment timepoint", limits = timepoint_levels, drop = FALSE) +
    scale_y_continuous(
        "Number of spots per cell (covariate-adjusted)",
        expand = expansion(mult = c(0.05, 0.08)),
        breaks = function(limits) {
            ticks <- pretty(data_limits, n = 5)
            ticks[ticks >= data_limits[1] & ticks <= data_limits[2]]
        }
    ) +
    scale_fill_manual(
        "Condition / Treatment",
        values = condition_colours,
        breaks = label_levels,
        labels = legend_labels,
        drop = FALSE
    ) +
    labs(title = plot_title) +
    theme_bw(base_size = 15) +
    theme(
        legend.position = "top",
        text = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90", linetype = "dashed", linewidth = 0.2),
        plot.title = element_text(size = 16, hjust = 0.5),
        strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 1),
        strip.text = element_text(size = 15)
    )

    # Coordinate limits leave all observations available to boxplot statistics.
    bounds <- if (is.null(y_limits)) c(NA_real_, NA_real_) else y_limits
    visible_values <- plot_df$plot_value[
        (is.na(bounds[1]) | plot_df$plot_value >= bounds[1]) &
            (is.na(bounds[2]) | plot_df$plot_value <= bounds[2])
    ]
    if (length(visible_values) == 0) {
        stop("No plot data fall within y_limits.")
    }
    visible_range <- range(visible_values)
    range_hint <- c(
        if (is.na(bounds[1])) visible_range[1] else bounds[1],
        if (is.na(bounds[2])) visible_range[2] else bounds[2]
    )
    visible_span <- max(diff(range_hint), 1)
    data_limits <- c(
        if (is.na(bounds[1])) visible_range[1] - 0.05 * visible_span else bounds[1],
        if (is.na(bounds[2])) visible_range[2] + 0.08 * visible_span else bounds[2]
    )
    if (data_limits[1] >= data_limits[2]) {
        stop("The lower y-limit must be smaller than the upper y-limit.")
    }
    display_limits <- data_limits

    if (!is.null(annotation_df)) {
        # Keep the data viewport separate from the annotation band above it.
        data_ceiling <- data_limits[2]
        data_span <- diff(data_limits)
        stack_step <- if (is.null(annotation_step)) max(0.1 * data_span, 0.5) else annotation_step
        band_gap <- if (is.null(annotation_gap)) max(0.04 * data_span, 0.3) else annotation_gap
        tip_depth <- 0.3 * stack_step
        label_offset <- 0.15 * stack_step
        label_padding <- 1 * stack_step
        annotation_df$branch_height <- ifelse(annotation_df$by_treatment, 0, 0.65 * stack_step)
        for (current_bin in unique(as.character(annotation_df$bin))) {
            facet_rows <- which(as.character(annotation_df$bin) == current_bin)
            automatic_rows <- facet_rows[is.na(annotation_df$y_position[facet_rows])]
            manual_heights <- annotation_df$y_position[facet_rows]
            manual_heights <- manual_heights[is.finite(manual_heights)]
            if (length(automatic_rows) > 0) {
                next_child_height <- max(
                    data_ceiling + band_gap + tip_depth,
                    manual_heights + stack_step
                )
                for (i in automatic_rows) {
                    annotation_df$y_position[i] <- next_child_height + annotation_df$branch_height[i]
                    next_child_height <- annotation_df$y_position[i] + stack_step
                }
            }
        }
        annotation_df$x_middle <- (annotation_df$x_start + annotation_df$x_end) / 2
        annotation_df$child_y <- annotation_df$y_position - annotation_df$branch_height
        annotation_df$tip_bottom <- annotation_df$child_y - tip_depth
        annotation_df$label_y <- annotation_df$y_position + label_offset
        if (any(annotation_df$tip_bottom <= data_ceiling)) {
            stop("Manual y_position values must place all bracket tips above ", data_ceiling,
                 ". Remove y_position for automatic placement.")
        }
        display_limits[2] <- max(annotation_df$label_y) + label_padding

        child_rows <- parent_rows <- vector("list", nrow(annotation_df))
        for (i in seq_len(nrow(annotation_df))) {
            current <- annotation_df[i, , drop = FALSE]
            if (current$by_treatment) {
                left <- current$x_start
                right <- current$x_end
            } else {
                # Each small bracket compares Scr and ERAP2 at one timepoint.
                members <- unname(label_offsets[paste(current$condition, treatment_levels, sep = "_")])
                pair_half_width <- diff(members) / 2
                centres <- c(current$x_start, current$x_end)
                left <- centres - pair_half_width
                right <- centres + pair_half_width
                parent_rows[[i]] <- data.frame(
                    bin = rep(current$bin, 3),
                    x = c(centres[1], centres[1], centres[2]),
                    xend = c(centres[2], centres[1], centres[2]),
                    y = c(current$y_position, current$child_y, current$child_y),
                    yend = rep(current$y_position, 3)
                )
            }
            child_rows[[i]] <- data.frame(
                bin = rep(current$bin, length(left)),
                condition = rep(current$condition, length(left)),
                x_start = left,
                x_end = right,
                y_position = rep(current$child_y, length(left)),
                tip_bottom = rep(current$tip_bottom, length(left)),
                annotation_colour = rep(current$annotation_colour, length(left))
            )
        }
        child_brackets <- do.call(rbind, child_rows)
        parent_brackets <- do.call(rbind, parent_rows)

        # Paint over data above the requested ceiling without changing statistics.
        annotation_band <- data.frame(ymin = data_ceiling, ymax = Inf)
        current_plot <- current_plot +
            geom_rect(
                data = annotation_band,
                aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = "white", colour = NA, show.legend = FALSE
            ) +
            geom_hline(yintercept = data_ceiling, colour = "grey80", linewidth = 0.3) +
            geom_segment(
                data = child_brackets,
                aes(x = x_start, xend = x_end, y = y_position, yend = y_position, colour = annotation_colour),
                inherit.aes = FALSE, linewidth = 0.45
            ) +
            geom_segment(
                data = child_brackets,
                aes(x = x_start, xend = x_start, y = tip_bottom, yend = y_position, colour = annotation_colour),
                inherit.aes = FALSE, linewidth = 0.45
            ) +
            geom_segment(
                data = child_brackets,
                aes(x = x_end, xend = x_end, y = tip_bottom, yend = y_position, colour = annotation_colour),
                inherit.aes = FALSE, linewidth = 0.45
            ) +
            geom_text(
                data = annotation_df,
                aes(x = x_middle, y = label_y + 0.1, label = p_label),
                inherit.aes = FALSE, colour = "black", size = 5, vjust = 0
            ) +
            scale_colour_identity()
        if (!is.null(parent_brackets)) {
            current_plot <- current_plot + geom_segment(
                data = parent_brackets,
                aes(x = x, xend = xend, y = y, yend = yend),
                inherit.aes = FALSE, colour = "black", linewidth = 0.45
            )
        }
    }
    current_plot <- current_plot + coord_cartesian(
        xlim = x_display_limits, ylim = display_limits, expand = FALSE
    )

    if (length(file_name) == 1L && !is.na(file_name)) {
        output_file <- file.path(output_dir, file_name)
        dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
        ggsave(filename = output_file, plot = current_plot, width = 14, height = 8,
               units = "in", bg = "white")
    }
    current_plot
}

p <- plot_combined_bin_v3(
    model_frame_df = m_v1_06_model_frame_combined,
    treatments_p_values = treatments_p_values,
    # combined_bafa1_timepoints_p_values = combined_bafa1_timepoints_p_values,
    file_name = "01-08-plot-LC3-m_v1_06-covariate-adjusted-num-spots-per-cell-boxplot.pdf",
    y_limits = c(0, 6),
    dodge_width = 0.75,
    box_width = 0.55,
    annotation_step = 0.5
)

p
