# Adapted analysis; see provenance/source_ranges.json.
here::i_am("workflow/scripts/common.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
data_root <- normalizePath(args[[1]], mustWork = TRUE)
plot_dir <- args[[2]]
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- normalizePath(plot_dir, mustWork = TRUE)
analysis_plots_dir <- plot_dir
output_dir <- plot_dir
analysis_output_dir <- plot_dir
save_dir <- plot_dir
analysis_save_dir <- plot_dir
input <- function(path) {
  p <- file.path(data_root, path)
  if (!file.exists(p) || is.na(file.size(p)) || file.size(p) == 0) stop("Missing/empty input: ", p)
  p
}
set.seed(20261020)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))
my_theme <- theme_bw(base_size = 12) + theme(legend.position = "bottom")
current_theme <- my_theme
