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
bodipy_stats_list <- readRDS(file.path(analysis_save_dir, "2026_09_14-sampled_cell_metrics_summary.rds"))
# ==========================================================================
# Imports
# ==========================================================================
library("stringr")
library("GenomicRanges")
library("DESeq2")
library("rtracklayer")
library("ggplot2")
library("GenomicFeatures")
library("GenomeInfoDb")
require("Rsamtools")
require("STRINGdb")
require("umap")
require("geosphere")
require("TSP")
require("ggbeeswarm")
# source("../src/bodipy_cdfs/source.R")
# source("../src/rmats-analysis/src/source.R")
# ==========================================================================
# Analysis
# ==========================================================================
# 1. Build the membership matrix
# Rows = cells in X
# Columns = membership vectors in K

fetch_membership_matrix <- function(K, X) {
  membership_matrix <- do.call(
    cbind,
    lapply(K, function(k) {
      X %in% k
    })
  )

  rownames(membership_matrix) <- X
  membership_matrix
}

# 5. Prepare the reference ECDF of all M

make_cdf_reference <- function(M) {
  values <- sort(unique(M))
  reference_ecdf <- stats::ecdf(M)

  list(
    values = values,
    probabilities = reference_ecdf(values),
    ecdf = reference_ecdf
  )
}


# 6. Calculate the KS distance between selected M and all M
score_selection <- function(M, selected, reference) {
  if (!length(selected)) {
    stop("At least one mean must be selected.")
  }

  selected_ecdf <- stats::ecdf(M[selected])

  ks_D <- max(abs(
    reference$probabilities -
      selected_ecdf(reference$values)
  ))

  list(
    ecdf = selected_ecdf,
    ks_D = ks_D
  )
}

# 1. Initial selection: one central observation per rank group

select_quantiles <- function(M, q) {
  stopifnot(
    length(q) == 1L,
    is.finite(q),
    q >= 1L,
    q <= length(M),
    q == floor(q)
  )

  ordered_indices <- order(M)

  groups <- floor(
    (seq_along(M) - 1L) * q / length(M)
  ) + 1L

  selected <- integer(q)

  for (j in seq_len(q)) {
    positions <- which(groups == j)
    centre <- mean(range(positions))

    selected[j] <- ordered_indices[
      positions[which.min(abs(positions - centre))]
    ]
  }

  selected
}


# 2. Calculate membership statistics

membership_stats <- function(counts, ks_D) {
  represented <- sum(counts > 0L)

  c(
    total_cells = length(counts),
    selected_cells = represented,
    proportion_unique = if (represented > 0L) {
      sum(counts == 1L) / represented
    } else {
      NA_real_
    },
    proportion_shared = if (represented > 0L) {
      sum(counts > 1L) / represented
    } else {
      NA_real_
    },
    ks_D = ks_D
  )
}


# 3. Find an initial selection meeting the KS tolerance
# Full coverage is optional and is NOT required during optimisation.

quantile_sweep <- function(
    M,
    membership_matrix,
    reference,
    start_quantiles = 10L,
    max_quantiles = length(M),
    shape_tolerance = 0.05,
    initial_full_coverage = FALSE
) {
  stopifnot(
    start_quantiles >= 1L,
    start_quantiles == floor(start_quantiles),
    max_quantiles >= start_quantiles,
    max_quantiles <= length(M),
    max_quantiles == floor(max_quantiles)
  )

  if (initial_full_coverage &&
      any(rowSums(membership_matrix) == 0L)) {
    stop("Full initial coverage is impossible.")
  }

  history <- list()

  for (q in seq.int(start_quantiles, max_quantiles)) {
    selected <- select_quantiles(M, q)

    counts <- rowSums(
      membership_matrix[, selected, drop = FALSE]
    )

    fit <- score_selection(M, selected, reference)
    stats <- membership_stats(counts, fit$ks_D)

    history[[length(history) + 1L]] <- data.frame(
      quantiles = q,
      selected_cells = unname(stats["selected_cells"]),
      proportion_shared = unname(stats["proportion_shared"]),
      ks_D = fit$ks_D
    )

    coverage_ok <- !initial_full_coverage || all(counts > 0L)

    if (coverage_ok &&
        sum(counts) > 0L &&
        fit$ks_D <= shape_tolerance) {
      return(list(
        selected = selected,
        history = do.call(rbind, history)
      ))
    }
  }

  stop("No initial selection met the requirements within max_quantiles.")
}


# 4. Evaluate one removal or swap

evaluate_move <- function(
    M,
    selected,
    counts,
    remove,
    add,
    membership_matrix,
    reference,
    shape_tolerance
) {
  candidate_selected <- selected[selected != remove]
  candidate_counts <- counts - membership_matrix[, remove]

  if (!is.na(add)) {
    candidate_selected <- c(candidate_selected, add)
    candidate_counts <- candidate_counts +
      membership_matrix[, add]
  }

  if (!length(candidate_selected) ||
      !any(candidate_counts > 0L)) {
    return(NULL)
  }

  fit <- score_selection(
    M,
    candidate_selected,
    reference
  )

  if (fit$ks_D > shape_tolerance) {
    return(NULL)
  }

  stats <- membership_stats(candidate_counts, fit$ks_D)

  list(
    selected = candidate_selected,
    counts = candidate_counts,
    ecdf = fit$ecdf,
    ks_D = fit$ks_D,
    shared = unname(stats["proportion_shared"]),
    removed = remove,
    added = add
  )
}


# 5. Find the best strictly overlap-reducing removal or swap

best_membership_move <- function(
    M,
    selected,
    counts,
    membership_matrix,
    reference,
    shape_tolerance,
    allow_swaps = TRUE,
    improvement_tolerance = 1e-12
) {
  current_shared <- sum(counts > 1L) / sum(counts > 0L)

  unselected <- setdiff(seq_along(M), selected)

  # NA means remove without replacement.
  additions <- NA_integer_

  if (allow_swaps) {
    additions <- c(additions, unselected)
  }

  best <- NULL

  for (remove in selected) {
    for (add in additions) {
      candidate <- evaluate_move(
        M = M,
        selected = selected,
        counts = counts,
        remove = remove,
        add = add,
        membership_matrix = membership_matrix,
        reference = reference,
        shape_tolerance = shape_tolerance
      )

      if (is.null(candidate)) next

      # Require an actual reduction in shared-cell proportion.
      if (candidate$shared >=
          current_shared - improvement_tolerance) next

      if (is.null(best)) {
        best <- candidate
        next
      }

      better_overlap <- candidate$shared <
        best$shared - improvement_tolerance

      tied_overlap <- abs(candidate$shared - best$shared) <=
        improvement_tolerance

      # Among overlap ties, favour a better ECDF match.
      if (better_overlap ||
          (tied_overlap && candidate$ks_D < best$ks_D)) {
        best <- candidate
      }
    }
  }

  best
}


# 6. Record optimisation progress

optimisation_record <- function(
    step,
    action,
    removed,
    added,
    selected,
    counts,
    ks_D
) {
  stats <- membership_stats(counts, ks_D)

  data.frame(
    step = step,
    action = action,
    removed = removed,
    added = added,
    n_selected = length(selected),
    selected_cells = unname(stats["selected_cells"]),
    proportion_unique = unname(stats["proportion_unique"]),
    proportion_shared = unname(stats["proportion_shared"]),
    ks_D = ks_D
  )
}


# 7. Improve membership diversity subject to the KS constraint

optimise_selection <- function(
    M,
    selected,
    membership_matrix,
    reference,
    shape_tolerance = 0.05,
    allow_swaps = TRUE,
    max_steps = 100L
) {
  stopifnot(
    length(max_steps) == 1L,
    is.finite(max_steps),
    max_steps >= 0L,
    max_steps == floor(max_steps)
  )

  selected <- unique(selected)
  counts <- rowSums(
    membership_matrix[, selected, drop = FALSE]
  )

  current <- score_selection(M, selected, reference)

  if (current$ks_D > shape_tolerance) {
    stop("Initial selection exceeds shape_tolerance.")
  }

  if (!any(counts > 0L)) {
    stop("Initial selection contains no represented cells.")
  }

  history <- list(optimisation_record(
    step = 0L,
    action = "initial",
    removed = NA_integer_,
    added = NA_integer_,
    selected = selected,
    counts = counts,
    ks_D = current$ks_D
  ))

  stop_reason <- "Maximum optimisation steps reached."

  for (step in seq_len(max_steps)) {
    best <- best_membership_move(
      M = M,
      selected = selected,
      counts = counts,
      membership_matrix = membership_matrix,
      reference = reference,
      shape_tolerance = shape_tolerance,
      allow_swaps = allow_swaps
    )

    if (is.null(best)) {
      stop_reason <- paste(
        "No tested move reduces sharing",
        "while satisfying the KS tolerance."
      )
      break
    }

    selected <- best$selected
    counts <- best$counts
    current <- best

    history[[length(history) + 1L]] <- optimisation_record(
      step = step,
      action = if (is.na(best$added)) "remove" else "swap",
      removed = best$removed,
      added = best$added,
      selected = selected,
      counts = counts,
      ks_D = best$ks_D
    )
  }

  list(
    selected = selected,
    counts = counts,
    ecdf = current$ecdf,
    ks_D = current$ks_D,
    history = do.call(rbind, history),
    stop_reason = stop_reason
  )
}


# 8. Coordinate initialisation and optimisation

select_and_fit <- function(
    M,
    membership_matrix,
    start_quantiles = 10L,
    max_quantiles = length(M),
    shape_tolerance = 0.05,
    initial_full_coverage = FALSE,
    allow_swaps = TRUE,
    max_steps = 100L
) {
  stopifnot(
    is.numeric(M),
    length(M) > 0L,
    all(is.finite(M)),
    is.matrix(membership_matrix),
    nrow(membership_matrix) > 0L,
    ncol(membership_matrix) == length(M),
    !anyNA(membership_matrix),
    all(membership_matrix %in% c(0, 1)),
    length(shape_tolerance) == 1L,
    is.finite(shape_tolerance),
    shape_tolerance >= 0,
    shape_tolerance <= 1
  )

  reference <- make_cdf_reference(M)

  initial <- quantile_sweep(
    M = M,
    membership_matrix = membership_matrix,
    reference = reference,
    start_quantiles = start_quantiles,
    max_quantiles = max_quantiles,
    shape_tolerance = shape_tolerance,
    initial_full_coverage = initial_full_coverage
  )

  optimised <- optimise_selection(
    M = M,
    selected = initial$selected,
    membership_matrix = membership_matrix,
    reference = reference,
    shape_tolerance = shape_tolerance,
    allow_swaps = allow_swaps,
    max_steps = max_steps
  )

  list(
    selected = optimised$selected,
    stats = membership_stats(
      optimised$counts,
      optimised$ks_D
    ),
    means = M[optimised$selected],
    ecdf = optimised$ecdf,
    reference_ecdf = reference$ecdf,
    ks_D = optimised$ks_D,
    within_tolerance = optimised$ks_D <= shape_tolerance,
    membership_counts = optimised$counts,
    initial_selected = initial$selected,
    quantile_sweep = initial$history,
    optimisation_history = optimised$history,
    stop_reason = optimised$stop_reason
  )
}

# Evaluate one proposed addition or swap.

evaluate_coverage_move <- function(
    M,
    selected,
    counts,
    add,
    membership_matrix,
    reference,
    shape_tolerance,
    max_shared,
    remove = NA_integer_
) {
  candidate_selected <- selected
  candidate_counts <- counts

  if (!is.na(remove)) {
    candidate_selected <- selected[selected != remove]
    candidate_counts <- candidate_counts -
      membership_matrix[, remove]
  }

  candidate_selected <- c(candidate_selected, add)
  candidate_counts <- candidate_counts +
    membership_matrix[, add]

  n_covered <- sum(candidate_counts > 0L)
  coverage_gain <- n_covered - sum(counts > 0L)

  # Only consider moves that increase distinct-cell coverage.
  if (coverage_gain <= 0L) return(NULL)

  proportion_shared <- sum(candidate_counts > 1L) / n_covered

  if (proportion_shared > max_shared + 1e-12) return(NULL)

  fit <- score_selection(
    M,
    candidate_selected,
    reference
  )

  if (fit$ks_D > shape_tolerance) return(NULL)

  list(
    selected = candidate_selected,
    counts = candidate_counts,
    ecdf = fit$ecdf,
    ks_D = fit$ks_D,
    coverage_gain = coverage_gain,
    proportion_shared = proportion_shared,
    added = add,
    removed = remove
  )
}


# Compare eligible moves:
# greatest coverage gain, then lowest sharing, then lowest KS distance.

better_coverage_move <- function(candidate, best) {
  if (is.null(candidate)) return(FALSE)
  if (is.null(best)) return(TRUE)

  if (candidate$coverage_gain != best$coverage_gain) {
    return(candidate$coverage_gain > best$coverage_gain)
  }

  if (candidate$proportion_shared != best$proportion_shared) {
    return(candidate$proportion_shared < best$proportion_shared)
  }

  candidate$ks_D < best$ks_D
}


# Find the best addition or swap.

best_coverage_move <- function(
    M,
    selected,
    counts,
    membership_matrix,
    reference,
    shape_tolerance,
    max_shared,
    allow_swaps = TRUE
) {
  unselected <- setdiff(seq_along(M), selected)
  best <- NULL

  removals <- NA_integer_

  if (allow_swaps) {
    removals <- c(removals, selected)
  }

  for (add in unselected) {
    # A candidate with no currently uncovered cells cannot
    # increase coverage through either addition or swapping.
    adds_new_cells <- any(
      membership_matrix[, add] & (counts == 0L)
    )

    if (!adds_new_cells) next

    for (remove in removals) {
      candidate <- evaluate_coverage_move(
        M = M,
        selected = selected,
        counts = counts,
        add = add,
        membership_matrix = membership_matrix,
        reference = reference,
        shape_tolerance = shape_tolerance,
        max_shared = max_shared,
        remove = remove
      )

      if (better_coverage_move(candidate, best)) {
        best <- candidate
      }
    }
  }

  best
}


# Record coverage and shape after each accepted move.

coverage_record <- function(
    step,
    action,
    added,
    removed,
    selected,
    counts,
    ks_D
) {
  represented <- sum(counts > 0L)

  data.frame(
    step = step,
    action = action,
    added = added,
    removed = removed,
    n_selected = length(selected),
    selected_cells = represented,
    proportion_covered = represented / length(counts),
    proportion_unique = sum(counts == 1L) / represented,
    proportion_shared = sum(counts > 1L) / represented,
    ks_D = ks_D
  )
}


# Increase coverage from an existing selection.
# By default, sharing cannot exceed its starting proportion.

increase_coverage <- function(
    M,
    membership_matrix,
    selected,
    shape_tolerance = 0.05,
    max_shared = NULL,
    allow_swaps = TRUE,
    max_steps = 100L
) {
  stopifnot(
    is.numeric(M),
    length(M) > 0L,
    all(is.finite(M)),
    is.matrix(membership_matrix),
    nrow(membership_matrix) > 0L,
    ncol(membership_matrix) == length(M),
    !anyNA(membership_matrix),
    all(membership_matrix %in% c(0, 1)),
    is.numeric(selected),
    length(selected) > 0L,
    all(is.finite(selected)),
    all(selected == floor(selected)),
    all(selected >= 1L & selected <= length(M)),
    length(shape_tolerance) == 1L,
    is.finite(shape_tolerance),
    shape_tolerance >= 0,
    shape_tolerance <= 1,
    length(allow_swaps) == 1L,
    is.logical(allow_swaps),
    !is.na(allow_swaps),
    length(max_steps) == 1L,
    is.finite(max_steps),
    max_steps >= 0L,
    max_steps == floor(max_steps)
  )

  selected <- unique(selected)
  initial_selected <- selected

  reference <- make_cdf_reference(M)

  counts <- rowSums(
    membership_matrix[, selected, drop = FALSE]
  )

  represented <- sum(counts > 0L)

  if (represented == 0L) {
    stop("The initial selection contains no represented cells.")
  }

  current <- score_selection(M, selected, reference)
  current_shared <- sum(counts > 1L) / represented

  if (is.null(max_shared)) {
    max_shared <- current_shared
  }

  stopifnot(
    length(max_shared) == 1L,
    is.finite(max_shared),
    max_shared >= 0,
    max_shared <= 1
  )

  if (current$ks_D > shape_tolerance) {
    stop("Initial selection exceeds shape_tolerance.")
  }

  if (current_shared > max_shared + 1e-12) {
    stop("Initial selection exceeds max_shared.")
  }

  history <- list(coverage_record(
    step = 0L,
    action = "initial",
    added = NA_integer_,
    removed = NA_integer_,
    selected = selected,
    counts = counts,
    ks_D = current$ks_D
  ))

  stop_reason <- "Maximum coverage steps reached."

  for (step in seq_len(max_steps)) {
    if (all(counts > 0L)) {
      stop_reason <- "Full coverage achieved."
      break
    }

    best <- best_coverage_move(
      M = M,
      selected = selected,
      counts = counts,
      membership_matrix = membership_matrix,
      reference = reference,
      shape_tolerance = shape_tolerance,
      max_shared = max_shared,
      allow_swaps = allow_swaps
    )

    if (is.null(best)) {
      stop_reason <- paste(
        "No tested addition or swap increases coverage",
        "within the sharing and KS limits."
      )
      break
    }

    selected <- best$selected
    counts <- best$counts
    current <- best

    history[[length(history) + 1L]] <- coverage_record(
      step = step,
      action = if (is.na(best$removed)) "add" else "swap",
      added = best$added,
      removed = best$removed,
      selected = selected,
      counts = counts,
      ks_D = current$ks_D
    )
  }

  if (all(counts > 0L)) {
    stop_reason <- "Full coverage achieved."
  }

  represented <- sum(counts > 0L)

  stats <- c(
    total_cells = nrow(membership_matrix),
    selected_cells = represented,
    proportion_unique = sum(counts == 1L) / represented,
    proportion_shared = sum(counts > 1L) / represented,
    ks_D = current$ks_D
  )

  list(
    selected = selected,
    stats = stats,
    means = M[selected],
    ecdf = current$ecdf,
    reference_ecdf = reference$ecdf,
    membership_counts = counts,
    initial_selected = initial_selected,
    coverage_history = do.call(rbind, history),
    shape_tolerance = shape_tolerance,
    max_shared = max_shared,
    stop_reason = stop_reason
  )
}



bodipy_stats_split_by_label <- split(bodipy_stats_list, bodipy_stats_list$label)
a_metric <- "mean_absolute_spot_area_um2"
q <- 10
lapply(bodipy_stats_split_by_label, function(a_bodipy_df){
    all_cells <- unique(do.call(c, a_bodipy_df$cell_id_list))

    membership_matrix <- fetch_membership_matrix(
    K = a_bodipy_df$cell_id_list,
    X = all_cells
    )

    # M must follow the same row order as a_bodipy_df.
    M <- a_bodipy_df[,a_metric]
    out <- select_and_fit(
        M = M,
        membership_matrix = membership_matrix,
        start_quantiles = 10,
        shape_tolerance = 0.05,
        initial_full_coverage = FALSE,
        allow_swaps = TRUE
    )

    expanded <- increase_coverage(
        M = M,
        membership_matrix = membership_matrix,
        selected = out$selected,
        shape_tolerance = 0.05,
        max_shared = 0.25
    )

    rbind(
    before = out$stats,
    after = expanded$stats
    )

    expanded$coverage_history
    expanded$stop_reason


    selected_rows <- a_bodipy_df[out$selected, , drop = FALSE]
    out$quantile_sweep
    out$pruning_history
    out$within_tolerance
    out$stop_reason
    # cur_quantile <- select_quantiles(M = a_bodipy_df[,a_metric], q)
})

lapply(bodipy_stats_split_by_label, function(a_bodipy_df){
    all_cells <- unique(do.call(c, a_bodipy_df$cell_id_list))

    membership_matrix <- fetch_membership_matrix(
    K = a_bodipy_df$cell_id_list,
    X = all_cells
    )

    # M must follow the same row order as a_bodipy_df.
    out <- select_and_fit(
    M = M,
    membership_matrix = membership_matrix,
    densfun = "normal",
    cdf = stats::pnorm,
    shape_tolerance = 0.05
    )

    selected_rows <- a_bodipy_df[out$selected, , drop = FALSE]

    out$pruning_history
    out$fit
    out$stop_reason
})


