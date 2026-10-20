data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects", "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "99-01-bodipy-erap2-kd")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
bodypy_df <- read.csv(file.path(analysis_input_dir, "test_spots_object.csv"))
process_files <- c("150324_for_Koustav.txt", "190424_for_Koustav.txt", "210624_for_Koustav.txt", "280424_for_Koustav.txt")
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
source("../src/bodipy_cdfs/source.R")
# source("../src/rmats-analysis/src/source.R")
# ==========================================================================
# Analysis
# ==========================================================================
row_annotation_list <- list(
    "150324_for_Koustav.txt" = data.frame(
    Row = seq(1,5),
    cell_line = c("CTRL1", "CTRL2", "CTRL6", "GliB", "NCRMC2"),
    condition = c("control", "control", "control", "vcp", "vcp")),

    "190424_for_Koustav.txt" = data.frame(
    Row = seq(1,8),
    treatment = c("Untx", "Untx", "Lipo", "Lipo", "Scramble", "Scramble", "ERAP2", "ERAP2"),
    replicate = c(1,2,1,2,1,2,1,2)),

    "280424_for_Koustav.txt" = data.frame(
    Row = seq(1,6),
    cell_line = c("CTRL1", "CTRL2", "CTRL6", "GliB", "NCRMC2", "NCRME6"),
    condition = c("control", "control", "control", "vcp", "vcp", "vcp")),

    "210624_for_Koustav.txt" = data.frame(
    Row = seq(1,8),
    treatment = c("Untx", "Untx", "Lipo", "Lipo", "Scramble", "Scramble", "ERAP2", "ERAP2"),
    replicate = c(1,2,1,2,1,2,1,2)))

col_annotation_list <- list(
    "150324_for_Koustav.txt" = data.frame(
    Column = seq(1,10),
    treatment = c("Untx", "Untx", "Lipo", "Lipo", "Scramble", "Scramble", "ERAP2", "ERAP2", "TIC", "TIC"),
    replicate = c(1,2,1,2,1,2,1,2,1,2)),

    "190424_for_Koustav.txt" = data.frame(
    Column = seq(1,7),
    cell_line = c("CTRL1", "CTRL2", "CTRL6", "GliA", "GliB", "NCRMC2", "NCRME6"),
    condition = c("control", "control", "control", "vcp", "vcp", "vcp", "vcp")),

    "280424_for_Koustav.txt" = data.frame(
    Column = seq(1,10),
    treatment = c("Untx", "Untx", "Lipo", "Lipo", "Scramble", "Scramble", "ERAP2", "ERAP2", "TIC", "TIC"),
    replicate = c(1,2,1,2,1,2,1,2,1,2)),

    "210624_for_Koustav.txt" = data.frame(
    Column = seq(1,9),
    cell_line = c("CTRL1", "CTRL2", "CTRL3", "CTRL6", "GliA", "GliB", "CB1E", "NCRMC2", "NCRME6"),
    condition = c("control", "control", "control", "control", "vcp", "vcp", "vcp", "vcp", "vcp")))

cell_columns <- lapply(process_files, function(x){
    bodypy_df <- read.table(file.path(analysis_input_dir, x), sep= "\t", skip = 9, header = T, nrow = 100)
    colnames(bodypy_df)[23]
})
names(cell_columns) <- process_files

spot_area_columns <- lapply(process_files, function(x){
    bodypy_df <- read.table(file.path(analysis_input_dir, x), sep= "\t", skip = 9, header = T, nrow = 100)
    colnames(bodypy_df)[20]
})
names(spot_area_columns) <- process_files

all_colnames <- lapply(process_files, function(x){
    bodypy_df <- read.table(file.path(analysis_input_dir, x), sep= "\t", skip = 9, header = T, nrow = 100)
    colnames(bodypy_df)
})[[1]]

data_object_list <- lapply(process_files, function(x){
    row_info_df <- row_annotation_list[[x]]
    col_info_df <- col_annotation_list[[x]]
    bodypy_df <- read.table(file.path(analysis_input_dir, x), sep= "\t", skip = 9, header = T)
    colnames(bodypy_df) <- all_colnames
    message(x)
    if(x %in% c("150324_for_Koustav.txt", "280424_for_Koustav.txt")){
        x_with_treatment <- assign_labels_to_columns(x = bodypy_df, 
            experiment_df = col_info_df, 
            map_col = "Column",
            value.var = c("treatment", "replicate"))
        x_with_cond <- assign_labels_to_columns(x = x_with_treatment, 
            experiment_df = row_info_df, 
            map_col = "Row",
            value.var = c("cell_line", "condition"))
    }else{
        x_with_treatment <- assign_labels_to_columns(x = bodypy_df, 
            experiment_df = row_info_df, 
            map_col = "Row",
            value.var = c("treatment", "replicate"))
        x_with_cond <- assign_labels_to_columns(x = x_with_treatment, 
            experiment_df = col_info_df, 
            map_col = "Column",
            value.var = c("cell_line", "condition"))        
    }
    rownames(x_with_cond) <- NULL
    x_with_cond <- x_with_cond[!(x_with_cond$cell_line %in% c("CTRL3", "CB1E")) & 
        x_with_cond$treatment != "TIC",]
    total_spot_area <- do.call(rbind, lapply(split(x_with_cond, paste(x_with_cond$treatment, x_with_cond$replicate, x_with_cond$cell_line, x_with_cond$condition)), function(temp_df){
        sum(temp_df[,spot_area_columns[[1]]]) 
    }))
    split_names <- str_split(rownames(total_spot_area), pattern = " ", simplify = T)
    spot_area_df <- data.frame(
        spot_area = total_spot_area, 
        batch = x, 
        treatment = split_names[,1],
        replicate = split_names[,2],
        line = split_names[,3],
        condition = split_names[,4])
    return(list(spot_area_df = spot_area_df, data_df = x_with_cond))
})

corrected_data_object_list <- lapply(data_object_list, function(x){
    x_with_cond <- x$data_df 
    spot_area_df <- x$spot_area_df
    scaled_signal_x_df <- do.call(rbind, lapply(split(x_with_cond, paste(x_with_cond$treatment, x_with_cond$replicate, x_with_cond$cell_line, x_with_cond$condition)), function(temp_df){
        # current_spot_area_df <- all_spot_area_df[paste(all_spot_area_df$batch, 
        #     all_spot_area_df$treatment, all_spot_area_df$replicate, 
        #     all_spot_area_df$line, all_spot_area_df$condition) == 
        #     paste(unique(spot_area_df$batch), unique(temp_df$treatment), unique(temp_df$replicate), unique(temp_df$cell_line), unique(temp_df$condition)),]
            # message(current_spot_area_df$adjusted_spot_area)
        # temp_df$normalised_spot_area <- temp_df[,spot_area_columns[[1]]] / (current_spot_area_df$adjusted_spot_area/1e6)
        temp_df$max_normalised_spot_area <- temp_df[,spot_area_columns[[1]]] / (max(temp_df[,spot_area_columns[[1]]]))
        # temp_df$scaled_spot_area <- scale(temp_df[,spot_area_columns[[1]]])
        return(temp_df)
    }))
    scaled_signal_x_df$experiment <- unique(spot_area_df$batch)
    return(scaled_signal_x_df)
})

data_object_df <- do.call(rbind, corrected_data_object_list)
rownames(data_object_df) <- NULL


data_object <- spot_analysR(data_object_df, 
experiment_block = 1, 
vars = c("experiment", "cell_line", "condition", "treatment", "replicate"),
var_order = c("experiment", "cell_line", "condition", "treatment", "replicate"),
field_col = "Field",
cell_col = cell_columns[[1]])


density_data_df <- plot_data_density(object = data_object,
    plot_col = "max_normalised_spot_area",
    log = T, 
    return_data = T,
    file = file.path(analysis_plots_dir, paste("01-plot-spot-area-density.pdf", sep = "-")))



quantiles <- quantile(density_data_df$x, c(0, 0.25, 0.5, 0.75, 0.99))
cell_count_data_df <- do.call(rbind, lapply(seq_along(quantiles[-length(quantiles)]), function(y){
    message(y)
    use_data <- data_object_df[data_object_df[,"max_normalised_spot_area"] > quantiles[y] & data_object_df[,"max_normalised_spot_area"] <= quantiles[y+1],]
    rownames(use_data) <- NULL

    data_object <- spot_analysR(use_data, 
    experiment_block = 1, 
    vars = c("experiment", "cell_line", "condition", "treatment", "replicate"),
    var_order = c("experiment", "cell_line", "condition", "treatment", "replicate"),
    field_col = "Field",
    cell_col = cell_columns[[1]])
    cell_count_df <- fetch_cell_count(data_object)
    cell_count_df$spot_threshold <- names(quantiles)[y]
    return(cell_count_df)
}))
cell_count_df_split <- split(cell_count_data_df, paste(cell_count_data_df$spot_threshold, paste(cell_count_data_df$experiment, cell_count_data_df$cell_line, cell_count_data_df$condition, cell_count_data_df$treatment)))

summarised_cell_count_df <- do.call(rbind, lapply(cell_count_df_split, function(x){
    return_df <- data.frame(sample = unique(paste(x$experiment, x$cell_line, x$condition, x$treatment, sep = "_")),
        experiment = unique(x$experiment),
        cell_line = unique(x$cell_line),
        condition = unique(x$condition),
        treatment = unique(x$treatment),
        spot_threshold = unique(x$spot_threshold),
        cell_count = sum(x$num_cells))
    return(return_df)
}))
rownames(summarised_cell_count_df) <- NULL
cell_count_df_split <- split(summarised_cell_count_df, summarised_cell_count_df$spot_threshold)
all_samples <- unique(summarised_cell_count_df$sample)

cell_count_matrix <- do.call(rbind, lapply(cell_count_df_split, function(x){
    return_vector <- rep(F, length(all_samples))
    names(return_vector) <- all_samples
    a_vector <- x$cell_count > 100
    names(a_vector) <- x$sample
    return_vector[names(a_vector)] <- T
    return(return_vector)
}))

cell_count_matrix[,colSums(cell_count_matrix) != 4]
> cell_count_matrix[,colSums(cell_count_matrix) != 4]

#     190424_for_Koustav.txt_CTRL1_control_Scramble
# 0%                                          FALSE
# 25%                                          TRUE
# 50%                                          TRUE
# 75%                                          TRUE
#     190424_for_Koustav.txt_GliA_vcp_Scramble
# 0%                                     FALSE
# 25%                                     TRUE
# 50%                                     TRUE
# 75%                                     TRUE
#     190424_for_Koustav.txt_GliA_vcp_Untx
# 0%                                 FALSE
# 25%                                 TRUE
# 50%                                 TRUE
# 75%                                 TRUE
#     280424_for_Koustav.txt_CTRL1_control_Untx
# 0%                                      FALSE
# 25%                                      TRUE
# 50%                                      TRUE
# 75%                                      TRUE
#     280424_for_Koustav.txt_CTRL2_control_Untx
# 0%                                      FALSE
# 25%                                      TRUE
# 50%                                      TRUE
# 75%                                      TRUE




quantiles <- quantile(density_data_df$x, c(0, 0.25, 0.5, 0.75, 0.999))
all_threshold_cdf_df <- do.call(rbind, lapply(seq_along(quantiles[-length(quantiles)]), function(y){
    message(y)
    use_data <- data_object_df[data_object_df[,"max_normalised_spot_area"] > quantiles[y] & data_object_df[,"max_normalised_spot_area"] <= quantiles[y+1],]
    rownames(use_data) <- NULL

    data_object <- spot_analysR(use_data, 
    experiment_block = 1, 
    vars = c("experiment", "cell_line", "condition", "treatment", "replicate"),
    var_order = c("experiment", "cell_line", "condition", "treatment", "replicate"),
    field_col = "Field",
    cell_col = cell_columns[[1]])
    # return(fetch_cell_count(data_object))

    data_object_with_summary <- summarise_data(data_object, 
        summarise_col = "max_normalised_spot_area",
        col_name = "spot_area")

    data_object_with_params <- prep_params(data_object_with_summary, select_cells = 100, iteration_size = 300, duplication_threshold = 1)
    data_object_with_seeds <- generate_seeds(data_object_with_params)

    data_object_with_sampled_cells <- sample_cells_across_replicates(data_object_with_seeds)
    data_object_with_range <- get_data_range(data_object_with_sampled_cells)
    data_object_with_quantile <- get_data_quantiles(data_object_with_range, fends = 30)

    data_object_with_cdf <- compute_cdfs(data_object_with_quantile)

    cdf_df <- fetch_metadata(data_object_with_cdf, "cdf_function")[[1]]
    total_spot_area_cdf <- summarise_cdf_by_col(cdf_df, fun = "total", col_names = c("condition", "treatment"))
    total_spot_area_cdf$FUN <- "total"
    total_spot_area_cdf$spot_threshold <- names(quantiles)[y]
    
    mean_spot_area_cdf <- summarise_cdf_by_col(cdf_df, fun = "mean", col_names = c("condition", "treatment"))
    mean_spot_area_cdf$FUN <- "mean"
    mean_spot_area_cdf$spot_threshold <- names(quantiles)[y]

    number_of_spot_area_cdf <- summarise_cdf_by_col(cdf_df, fun = "n", col_names = c("condition", "treatment"))
    number_of_spot_area_cdf$FUN <- "n"
    number_of_spot_area_cdf$spot_threshold <- names(quantiles)[y]

    summarised_cdf <- rbind(total_spot_area_cdf, mean_spot_area_cdf, number_of_spot_area_cdf)
    rownames(summarised_cdf) <- NULL
    return(summarised_cdf)
}))
all_threshold_cdf_df$group <- paste(all_threshold_cdf_df$condition, all_threshold_cdf_df$treatment, sep = "_")


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
Colours <- c("control_ERAP2" = "#1f78b4", "control_Scramble" = "#a6cee3", "vcp_Scramble" = "#fb9a99", "vcp_ERAP2" =  "#e31a1c")
ThePlot <- ggplot(all_threshold_cdf_df[
    all_threshold_cdf_df$FUN == "total" & 
    all_threshold_cdf_df$treatment %in% c("ERAP2", "Scramble"),], aes(x = log10(quantiles), y = mean, group = group))
# ThePlot <- ThePlot + geom_ribbon(aes(ymin = lower, ymax = upper, fill = condition), alpha = 0.6)
ThePlot <- ThePlot + geom_step(aes(colour = group))
ThePlot <- ThePlot + scale_x_continuous("log10 mean total spot area per cell")
ThePlot <- ThePlot + scale_y_continuous("CDF")
ThePlot <- ThePlot + scale_colour_manual("condition", values = Colours)
ThePlot <- ThePlot + scale_fill_brewer("condition", palette = "Set1")
# ThePlot <- ThePlot + scale_colour_manual("status", values = Colours)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + facet_grid(condition ~ spot_threshold, scales = "free")
ThePlot <- ThePlot + ggtitle("CDF of mean total spot area per cell")
ggsave(width = 10, height = 5, file = file.path(analysis_plots_dir, "04-01-plot-cdf-of-average-total-max-normalised-spot-area-in-bodipy-experiments-facet-condition.pdf"), ThePlot)


ThePlot <- ggplot(all_threshold_cdf_df[
    all_threshold_cdf_df$FUN == "mean" & 
    all_threshold_cdf_df$treatment %in% c("ERAP2", "Scramble"),], aes(x = log10(quantiles), y = mean, group = group))
# ThePlot <- ThePlot + geom_ribbon(aes(ymin = lower, ymax = upper, fill = condition), alpha = 0.6)
ThePlot <- ThePlot + geom_step(aes(colour = group))
ThePlot <- ThePlot + scale_x_continuous("log10 mean spot area per cell")
ThePlot <- ThePlot + scale_y_continuous("CDF")
ThePlot <- ThePlot + scale_colour_manual("condition", values = Colours)
ThePlot <- ThePlot + scale_fill_brewer("condition", palette = "Set1")
# ThePlot <- ThePlot + scale_colour_manual("status", values = Colours)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + facet_grid(condition ~ spot_threshold, scales = "free")
ThePlot <- ThePlot + ggtitle("CDF of mean spot area per cell")
ggsave(width = 10, height = 5, file = file.path(analysis_plots_dir, "04-02-plot-cdf-of-average-max-normalised-spot-area-in-bodipy-experiments-facet-condition.pdf"), ThePlot)



ThePlot <- ggplot(all_threshold_cdf_df[
    all_threshold_cdf_df$FUN == "n" & 
    all_threshold_cdf_df$treatment %in% c("ERAP2", "Scramble"),], aes(x = log10(quantiles), y = mean, group = group))
# ThePlot <- ThePlot + geom_ribbon(aes(ymin = lower, ymax = upper, fill = condition), alpha = 0.6)
ThePlot <- ThePlot + geom_step(aes(colour = group))
ThePlot <- ThePlot + scale_x_continuous("log10 number of spots per cell")
ThePlot <- ThePlot + scale_y_continuous("CDF")
ThePlot <- ThePlot + scale_colour_manual("condition", values = Colours)
ThePlot <- ThePlot + scale_fill_brewer("condition", palette = "Set1")
# ThePlot <- ThePlot + scale_colour_manual("status", values = Colours)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + facet_grid(condition ~ spot_threshold, scales = "free")
ThePlot <- ThePlot + ggtitle("CDF of number of spots per cell")
ggsave(width = 10, height = 5, file = file.path(analysis_plots_dir, "04-02-plot-cdf-of-average-of-number-of-spots-in-bodipy-experiments-facet-condition.pdf"), ThePlot)



keep_cdf_df <- all_threshold_cdf_df[all_threshold_cdf_df$treatment %in% c("ERAP2", "Scramble"),]
cdf_split_df <- split(keep_cdf_df, paste(keep_cdf_df$FUN, keep_cdf_df$spot_threshold))
treatment_effect_df <- do.call(rbind, lapply(cdf_split_df, function(temp_df){
    condition_split <- split(temp_df, temp_df$condition)
    return_df <- do.call(rbind, lapply(condition_split, function(cond_df){
        pval <- ks.test(cond_df$mean[cond_df$treatment == "ERAP2"], cond_df$mean[cond_df$treatment == "Scramble"])$p.value
        data.frame(condition = unique(cond_df$condition), pval = pval)
    }))
    return_df$FUN <- unique(temp_df$FUN)
    return_df$spot_threshold <- unique(temp_df$spot_threshold)
    return(return_df)
}))
#                   condition       pval   FUN spot_threshold
# mean 0%.control     control 0.82345352  mean             0%
# mean 0%.vcp             vcp 0.96342971  mean             0%
# mean 25%.control    control 0.96342971  mean            25%
# mean 25%.vcp            vcp 1.00000000  mean            25%
# mean 50%.control    control 1.00000000  mean            50%
# mean 50%.vcp            vcp 0.99911526  mean            50%
# mean 75%.control    control 1.00000000  mean            75%
# mean 75%.vcp            vcp 1.00000000  mean            75%
# n 0%.control        control 0.24066829     n             0%
# n 0%.vcp                vcp 0.24066829     n             0%
# n 25%.control       control 0.10076790     n            25%
# n 25%.vcp               vcp 0.09643439     n            25%
# n 50%.control       control 0.21783088     n            50%
# n 50%.vcp               vcp 0.21417439     n            50%
# n 75%.control       control 0.20497003     n            75%
# n 75%.vcp               vcp 0.20457236     n            75%
# total 0%.control    control 1.00000000 total             0%
# total 0%.vcp            vcp 0.82345352 total             0%
# total 25%.control   control 1.00000000 total            25%
# total 25%.vcp           vcp 1.00000000 total            25%
# total 50%.control   control 1.00000000 total            50%
# total 50%.vcp           vcp 1.00000000 total            50%
# total 75%.control   control 1.00000000 total            75%
# total 75%.vcp           vcp 0.99911526 total            75%

condition_effect_df <- do.call(rbind, lapply(cdf_split_df, function(temp_df){
    condition_split <- split(temp_df, temp_df$condition)
    compare_df <- do.call(rbind, lapply(condition_split, function(cond_df){
        delta <- cond_df$mean[cond_df$treatment == "ERAP2"] - cond_df$mean[cond_df$treatment == "Scramble"]
        data.frame(condition = unique(cond_df$condition), deltas = delta)
    }))
    pval <- ks.test(cumsum(compare_df$deltas[compare_df$condition == "vcp"]), cumsum(compare_df$deltas[compare_df$condition == "control"]))$p.value
    greater <- ks.test(cumsum(compare_df$deltas[compare_df$condition == "vcp"]), cumsum(compare_df$deltas[compare_df$condition == "control"]), alternative = "greater")$p.value
    less <- ks.test(cumsum(compare_df$deltas[compare_df$condition == "vcp"]), cumsum(compare_df$deltas[compare_df$condition == "control"]), alternative = "less")$p.value
    return_df <- data.frame(FUN = unique(temp_df$FUN), 
        spot_threshold = unique(temp_df$spot_threshold),
        pval = pval,
        greater = greater,
        less = less)
    return(return_df)
}))

condition_effect_df$sig <- condition_effect_df$pval < 0.05
condition_effect_df$greater_sig <- condition_effect_df$greater < 0.05
condition_effect_df$less_sig <- condition_effect_df$less < 0.05



