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
process_spot_files <- c("150324_for_Koustav.txt", "190424_for_Koustav.txt", "210624_for_Koustav.txt", "280424_for_Koustav.txt")
process_cyto_files <- c("15_03_24_for_Koustav_cell_info.txt", "19_04_24_for_Koustav_cell_info.txt", "21_06_24_for_Koustav_cell_info.txt", "28_04_24_for_Koustav_cell_info.txt")
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

cell_columns <- lapply(process_spot_files, function(x){
    bodypy_df <- read.table(file.path(analysis_input_dir, x), sep= "\t", skip = 9, header = T, nrow = 100)
    colnames(bodypy_df)[23]
})
names(cell_columns) <- process_spot_files

spot_area_columns <- lapply(process_spot_files, function(x){
    bodypy_df <- read.table(file.path(analysis_input_dir, x), sep= "\t", skip = 9, header = T, nrow = 100)
    colnames(bodypy_df)[20]
})
names(spot_area_columns) <- process_spot_files

cytoplasm_area_columns <- lapply(process_cyto_files, function(x){
    bodypy_df <- read.table(file.path(analysis_input_dir, x), sep= "\t", skip = 9, header = T, nrow = 100)
    colnames(bodypy_df)[21]
})
names(cytoplasm_area_columns) <- process_cyto_files

cytoplasm_num_spots_per_area_column <- lapply(process_cyto_files, function(x){
    bodypy_df <- read.table(file.path(analysis_input_dir, x), sep= "\t", skip = 9, header = T, nrow = 100)
    colnames(bodypy_df)[grepl(x = colnames(bodypy_df), pattern = "Number.of.Spots.per.Area.of.Cytoplasm")]
})
names(cytoplasm_num_spots_per_area_column) <- process_cyto_files


spot_area_um_column <- lapply(process_cyto_files, function(x){
    bodypy_df <- read.table(file.path(analysis_input_dir, x), sep= "\t", skip = 9, header = T, nrow = 100)
    colnames(bodypy_df)[grepl(x = colnames(bodypy_df), pattern = "Spots.Area..µm..")]
})
names(spot_area_um_column) <- process_cyto_files


cell_info_cell_columns <- lapply(process_cyto_files, function(x){
    bodypy_df <- read.table(file.path(analysis_input_dir, x), sep= "\t", skip = 9, header = T, nrow = 100)
    colnames(bodypy_df)[20]
})
names(cell_info_cell_columns) <- process_cyto_files


all_colnames <- lapply(process_spot_files, function(x){
    bodypy_df <- read.table(file.path(analysis_input_dir, x), sep= "\t", skip = 9, header = T, nrow = 100)
    colnames(bodypy_df)
})[[1]]

data_object_list <- lapply(seq_along(process_spot_files), function(x){
    row_info_df <- row_annotation_list[[ process_spot_files[x] ]]
    col_info_df <- col_annotation_list[[ process_spot_files[x] ]]
    bodypy_df <- read.table(file.path(analysis_input_dir, process_spot_files[x]), sep= "\t", skip = 9, header = T)
    bodypy_df$temp_id <- paste(bodypy_df$Row, bodypy_df$Column, bodypy_df$Timepoint, bodypy_df$Field, bodypy_df[,cell_columns[[ process_spot_files[x] ]] ], sep = "_")
    bodypy_df$experiment_label <- paste(bodypy_df$Row, bodypy_df$Column, bodypy_df$Timepoint, sep = "_")


    # compute total spot area in px
    total_spot_area <- vapply(split(bodypy_df[,spot_area_columns[[ process_spot_files[x] ]]], bodypy_df[,"temp_id"]), sum, 1)
    # compute total number of spots
    total_n_of_spots <- vapply(split(bodypy_df[,spot_area_columns[[ process_spot_files[x] ]]], bodypy_df[,"temp_id"]), length, 1)
    # sd within cell


    cell_info_df <- read.table(file.path(analysis_input_dir, process_cyto_files[x]), sep= "\t", skip = 9, header = T)
    cell_info_df$temp_id <- paste(cell_info_df$Row, cell_info_df$Column, cell_info_df$Timepoint, cell_info_df$Field, cell_info_df$Object.No, sep = "_")
    cell_info_df$experiment_label <- paste(cell_info_df$Row, cell_info_df$Column, cell_info_df$Timepoint, sep = "_")
    
    match_bodypy_to_cell <- match(bodypy_df$temp_id, cell_info_df$temp_id)

    # Fetch total cytoplasmic area in cell in um
    cytoplasmic_area <- cell_info_df[ match_bodypy_to_cell , cytoplasm_area_columns[[ process_cyto_files[x] ]] ]
    # Fetch total cytoplasmic area in field
    total_cytoplasmic_area <- vapply(split(cell_info_df[, cytoplasm_area_columns[[ process_cyto_files[x] ]]], cell_info_df[,"experiment_label"]), sum, 1)
    # Fetch total spot area in cell in um
    spot_area_in_um <- cell_info_df[ match_bodypy_to_cell , spot_area_um_column[[ process_cyto_files[x] ]] ]    
    # Rescale spot area to 100um
    scaled_spot_area <- (spot_area_in_um / cytoplasmic_area) * 100


    summarised_cell_df <- unique(data.frame(experiment = process_spot_files[x],
            temp_id = bodypy_df$temp_id,
            experiment_label = bodypy_df$experiment_label,
            Row = bodypy_df$Row,
            Column = bodypy_df$Column,
            spot_area_in_um = spot_area_in_um, 
            total_cytoplasmic_area = total_cytoplasmic_area[bodypy_df$experiment_label],
            cytoplasmic_area_in_um = cytoplasmic_area, 
            spot_area_scaled_to_100um = scaled_spot_area,
            number_of_spots = total_n_of_spots[bodypy_df$temp_id]))
    # compute expected size of a spot per cell
    summarised_cell_df$mean_size_of_a_spot <- summarised_cell_df$spot_area_scaled_to_100um / summarised_cell_df$number_of_spots
    average_spot_size_in_experiment <- vapply(split(summarised_cell_df$mean_size_of_a_spot, summarised_cell_df$experiment_label), mean, 1)
    sd_spot_size_in_experiment <- vapply(split(summarised_cell_df$mean_size_of_a_spot, summarised_cell_df$experiment_label), sd, 1)

    summarised_cell_df$zscore <- (summarised_cell_df$mean_size_of_a_spot - average_spot_size_in_experiment[summarised_cell_df$experiment_label]) / sd_spot_size_in_experiment[summarised_cell_df$experiment_label]

    bodypy_df$experiment <- process_spot_files[x]
    bodypy_df$spot_area_in_um <- (bodypy_df[, spot_area_columns[[ process_spot_files[x] ]] ] / total_spot_area[bodypy_df[,"temp_id"]]) * spot_area_in_um
    bodypy_df$spot_area_scaled_to_100um <- (bodypy_df$spot_area_in_um / cytoplasmic_area) * 100

    std_dev_of_scaled_spot_area <- vapply(split(bodypy_df[,"spot_area_scaled_to_100um"], bodypy_df[,"temp_id"]), sd, 1)
    bodypy_df$relative_spot_size <- bodypy_df$spot_area_scaled_to_100um / summarised_cell_df$mean_size_of_a_spot[match(bodypy_df$temp_id, summarised_cell_df$temp_id)]
    # bodypy_df$spot_size_zscore <- (bodypy_df$spot_area_scaled_to_100um - summarised_cell_df$mean_size_of_a_spot[match(bodypy_df$temp_id, summarised_cell_df$temp_id)]) / std_dev_of_scaled_spot_area[bodypy_df[,"temp_id"]]


    if(process_spot_files[x] %in% c("150324_for_Koustav.txt", "280424_for_Koustav.txt")){
        x_with_treatment <- assign_labels_to_columns(x = summarised_cell_df, 
            experiment_df = col_info_df, 
            map_col = "Column",
            value.var = c("treatment", "replicate"))
        x_with_cond <- assign_labels_to_columns(x = x_with_treatment, 
            experiment_df = row_info_df, 
            map_col = "Row",
            value.var = c("cell_line", "condition"))
    }else{
        x_with_treatment <- assign_labels_to_columns(x = summarised_cell_df, 
            experiment_df = row_info_df, 
            map_col = "Row",
            value.var = c("treatment", "replicate"))
        x_with_cond <- assign_labels_to_columns(x = x_with_treatment, 
            experiment_df = col_info_df, 
            map_col = "Column",
            value.var = c("cell_line", "condition"))
    }
    if(process_spot_files[x] %in% c("150324_for_Koustav.txt", "280424_for_Koustav.txt")){
        y_with_treatment <- assign_labels_to_columns(x = bodypy_df, 
            experiment_df = col_info_df, 
            map_col = "Column",
            value.var = c("treatment", "replicate"))
        y_with_cond <- assign_labels_to_columns(x = y_with_treatment, 
            experiment_df = row_info_df, 
            map_col = "Row",
            value.var = c("cell_line", "condition"))
    }else{
        y_with_treatment <- assign_labels_to_columns(x = bodypy_df, 
            experiment_df = row_info_df, 
            map_col = "Row",
            value.var = c("treatment", "replicate"))
        y_with_cond <- assign_labels_to_columns(x = y_with_treatment, 
            experiment_df = col_info_df, 
            map_col = "Column",
            value.var = c("cell_line", "condition"))        
    }
    return(list(summarised_cell_df = x_with_cond, bodypy_df = y_with_cond))
})

# fetch all summarised cell data

all_summarised_cell_df <- do.call(rbind, lapply(data_object_list, function(x){
    x$summarised_cell_df
}))
all_bodypy_df <- do.call(rbind, lapply(data_object_list, function(x){
    colnames(x$bodypy_df) <- colnames(data_object_list[[1]]$bodypy_df)
    return(x$bodypy_df)
}))
all_summarised_cell_df <- all_summarised_cell_df[all_summarised_cell_df$treatment %in% c("ERAP2", "Scramble", "Untx") &
    !(all_summarised_cell_df$cell_line %in% c("CTRL3")),]

max_spot_size_per_cell <- do.call(rbind, lapply(split(all_bodypy_df, paste(all_bodypy_df$experiment, all_bodypy_df$temp_id)), function(x){
    data.frame(experiment = unique(x$experiment),
        temp_id = unique(x$temp_id),
        spot_area_scaled_to_100um = median(x$spot_area_scaled_to_100um),
        relative_spot_size = max(x$relative_spot_size))
}))


# Investigate the average effect using all cells 
# The simplest procedure possible
average_spot_size_of_cell_lines_df <- do.call(rbind, lapply(split(all_summarised_cell_df, paste(all_summarised_cell_df$experiment, all_summarised_cell_df$cell_line, all_summarised_cell_df$treatment)), function(temp_df){
    data.frame(
        experiment = unique(temp_df$experiment), 
        cell_line = unique(temp_df$cell_line),
        treatment = unique(temp_df$treatment),
        condition = unique(temp_df$condition),
        mean_size_of_a_spot = median(temp_df$mean_size_of_a_spot))
}))

rownames(average_spot_size_of_cell_lines_df) <- NULL
average_spot_size_of_cell_lines_df$label <- paste(average_spot_size_of_cell_lines_df$condition, average_spot_size_of_cell_lines_df$treatment, sep = "_")


lapply(split(average_spot_size_of_cell_lines_df, average_spot_size_of_cell_lines_df$treatment), function(temp_df){
    wilcox.test(temp_df$mean_size_of_a_spot[temp_df$condition == "vcp"], temp_df$mean_size_of_a_spot[temp_df$condition == "control"])$p.value
})

lapply(split(average_spot_size_of_cell_lines_df, average_spot_size_of_cell_lines_df$condition), function(temp_df){
    wilcox.test(temp_df$mean_size_of_a_spot[temp_df$treatment == "ERAP2"], temp_df$mean_size_of_a_spot[temp_df$treatment == "Scramble"])$p.value
})



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

average_spot_size_of_cell_lines_df$label <- factor(average_spot_size_of_cell_lines_df$label, levels = names(condition_colours))
ThePlot <- ggplot(average_spot_size_of_cell_lines_df[average_spot_size_of_cell_lines_df$treatment != "Untx",], aes(x = condition, y = mean_size_of_a_spot))
ThePlot <- ThePlot + geom_boxplot(aes(fill = label), colour = "#000000", 
    width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9))
ThePlot <- ThePlot + geom_jitter(aes(fill = label), colour = "#000000", position = position_jitterdodge(dodge.width = 0.9), shape = 21)
ThePlot <- ThePlot + scale_x_discrete("condition")
ThePlot <- ThePlot + scale_y_continuous("average (median) spot size scaled to 100 um")
ThePlot <- ThePlot + scale_colour_manual("condition", values = condition_colours)
ThePlot <- ThePlot + scale_fill_manual("condition", values = condition_colours)
# ThePlot <- ThePlot + scale_colour_manual("status", values = Colours)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ treatment, scales = "free")
ThePlot <- ThePlot + ggtitle("CDF of mean total spot area per cell")
ggsave(width = 5, height = 10, file = file.path(analysis_plots_dir, "05-01-plot-median-spot-size-of-a-spot.pdf"), ThePlot)





# all_bodipy_df <- do.call(rbind, lapply(data_object_list, function(x){
#     current_bodipy_df <- x$bodypy_df
#     summarised_cell_df <- x$summarised_cell_df

#     prop_cell_df <- do.call(rbind, lapply(zscore_thresholds, function(x){
#         keep_bodipy_df <- current_bodipy_df[current_bodipy_df$zscore > x,]
#         keep_summarised_cell_df <- summarised_cell_df[summarised_cell_df$temp_id %in% keep_bodipy_df$temp_id,]
#         keep_num_cells <- table(keep_summarised_cell_df$cell_line, keep_summarised_cell_df$treatment)
#         prop_cells <- keep_num_cells / total_num_cells
#         ret_prop_cells <- reshape2::melt(prop_cells)
#         ret_prop_cells$zscore <- x
#         return(ret_prop_cells)
#     }))

#     num_cells

# }))



unique_experiment_cytoarea_df <- unique(data.frame(total_cytoplasmic_area = all_summarised_cell_df$total_cytoplasmic_area,
            experiment_label = paste(all_summarised_cell_df$experiment, all_summarised_cell_df$experiment_label, sep = "_")))
total_cyto_area <- sum(unique_experiment_cytoarea_df$total_cytoplasmic_area)
unique_experiment_cytoarea_df$proportion_over_total <- unique_experiment_cytoarea_df$total_cytoplasmic_area / total_cyto_area

sample_cell_proportion <- table(paste(all_summarised_cell_df$experiment, all_summarised_cell_df$experiment_label, sep = "_")) / length(paste(all_summarised_cell_df$experiment, all_summarised_cell_df$experiment_label))

all_summarised_cell_df$sample_total_cyto_area_probability <- unique_experiment_cytoarea_df$proportion_over_total[match(paste(all_summarised_cell_df$experiment, all_summarised_cell_df$experiment_label, sep = "_"), unique_experiment_cytoarea_df$experiment_label)]
all_summarised_cell_df$cell_selection_probability <- all_summarised_cell_df$cytoplasmic_area_in_um / (total_cyto_area * sample_cell_proportion[paste(all_summarised_cell_df$experiment, all_summarised_cell_df$experiment_label, sep = "_")])

all_summarised_cell_df$relative_spot_size <- max_spot_size_per_cell[match(paste(all_summarised_cell_df$experiment, all_summarised_cell_df$temp_id), paste(max_spot_size_per_cell$experiment, max_spot_size_per_cell$temp_id)),"relative_spot_size"]
all_summarised_cell_df$spot_area_scaled_to_100um <- max_spot_size_per_cell[match(paste(all_summarised_cell_df$experiment, all_summarised_cell_df$temp_id), paste(max_spot_size_per_cell$experiment, max_spot_size_per_cell$temp_id)),"spot_area_scaled_to_100um"]






seed_list <- lapply(seq(1,1000), function(y){
    treatment_split <- split(all_summarised_cell_df, all_summarised_cell_df$treatment)
    temp_list <- lapply(names(treatment_split), function(a_treatment){
        treatment_df <- treatment_split[[a_treatment]]
        pval <- 0
        while(pval < 0.05){
            message(y)
            current_vector <- sample(x = seq_along(treatment_df$sample_total_cyto_area_probability), size = 1000, prob = treatment_df$sample_total_cyto_area_probability)
            current_treatment_df <- treatment_df[current_vector,]
            current_proportions <- table(paste(current_treatment_df$experiment, current_treatment_df$experiment_label, sep = "_")) / length(paste(current_treatment_df$experiment, current_treatment_df$experiment_label))
            pval <- wilcox.test(rank(current_proportions), rank(sample_cell_proportion[names(current_proportions)]))$p.value
        }
        return(current_vector)
    })
    names(temp_list) <- names(treatment_split)
    return(temp_list)
})

# how many cells were not selected
treatment_split <- split(all_summarised_cell_df, all_summarised_cell_df$treatment)
prop_coverage <- lapply(names(treatment_split), function(a_treatment){
    treatment_df <- treatment_split[[a_treatment]]
    selected_cells <- do.call(c, lapply(seq(1,1000), function(y){
            seed_list[[y]][[a_treatment]]
        }))
    all_cell_seq <- seq_len(nrow(treatment_df))
    return(length(all_cell_seq[all_cell_seq %in% selected_cells]) / length(all_cell_seq))
})



average_spot_size_df <- do.call(rbind, lapply(seq(1,1000), function(y){
    message(y)
    treatment_split <- split(all_summarised_cell_df, all_summarised_cell_df$treatment)
    do.call(rbind, lapply(names(treatment_split), function(a_treatment){
        treatment_df <- treatment_split[[a_treatment]]
        current_vector <- seed_list[[y]][[a_treatment]]
        current_selected_df <- treatment_df[current_vector,]

        mean_relative_spot_size <- vapply(split(current_selected_df$relative_spot_size, current_selected_df$condition), median, 1)
        spot_area_scaled_to_100um <- vapply(split(current_selected_df$spot_area_scaled_to_100um, current_selected_df$condition), median, 1)
        mean_spot_size <- vapply(split(current_selected_df$mean_size_of_a_spot, current_selected_df$condition), median, 1)
        data.frame(iteration = y, treatment = a_treatment, 
            condition = names(mean_spot_size), 
            median_spot_size = mean_spot_size,
            spot_area_scaled_to_100um = spot_area_scaled_to_100um,
            median_relative_spot_size = mean_relative_spot_size)
    }))
}))
average_spot_size_df$label <- paste(average_spot_size_df$condition, average_spot_size_df$treatment, sep = "_")


lapply(split(average_spot_size_df, average_spot_size_df$treatment), function(temp_df){
    wilcox.test(temp_df$median_relative_spot_size[temp_df$condition == "vcp"], temp_df$median_relative_spot_size[temp_df$condition == "control"])$p.value
})

lapply(split(average_spot_size_df, average_spot_size_df$condition), function(temp_df){
    wilcox.test(temp_df$median_relative_spot_size[temp_df$treatment == "ERAP2"], temp_df$median_relative_spot_size[temp_df$treatment == "Scramble"])$p.value
})


lapply(split(average_spot_size_df, average_spot_size_df$treatment), function(temp_df){
    wilcox.test(temp_df$spot_area_scaled_to_100um[temp_df$condition == "vcp"], temp_df$spot_area_scaled_to_100um[temp_df$condition == "control"])$p.value
})

lapply(split(average_spot_size_df, average_spot_size_df$condition), function(temp_df){
    wilcox.test(temp_df$spot_area_scaled_to_100um[temp_df$treatment == "ERAP2"], temp_df$spot_area_scaled_to_100um[temp_df$treatment == "Scramble"])$p.value
})


condition_colours <- c("control_Scramble" = "#a6cee3", "control_ERAP2" = "#1f78b4", "vcp_Scramble" = "#fb9a99", "vcp_ERAP2" =  "#e31a1c")
average_spot_size_df$label <- factor(average_spot_size_df$label, levels = names(condition_colours))
ThePlot <- ggplot(average_spot_size_df[average_spot_size_df$treatment != "Untx",], aes(x = condition, y = median_spot_size))
ThePlot <- ThePlot + geom_quasirandom(fill = "#F2FDFF", aes(group = label), colour = "#000000", shape = 21, dodge.width = 0.9, alpha = 0.6)
ThePlot <- ThePlot + geom_boxplot(aes(fill = label), colour = "#000000", 
    width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9))
ThePlot <- ThePlot + scale_x_discrete("condition")
ThePlot <- ThePlot + scale_y_continuous("average estimated spot size scaled to 100um")
ThePlot <- ThePlot + scale_colour_manual("condition", values = condition_colours)
ThePlot <- ThePlot + scale_fill_manual("condition", values = condition_colours)
# ThePlot <- ThePlot + scale_colour_manual("status", values = Colours)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ treatment, scales = "free")
ThePlot <- ThePlot + ggtitle("CDF of mean total spot area per cell")
ggsave(width = 5, height = 10, file = file.path(analysis_plots_dir, "05-02-plot-sampled-median-spot-size.pdf"), ThePlot)


condition_colours <- c("control_Scramble" = "#a6cee3", "control_ERAP2" = "#1f78b4", "vcp_Scramble" = "#fb9a99", "vcp_ERAP2" =  "#e31a1c")
average_spot_size_df$label <- factor(average_spot_size_df$label, levels = names(condition_colours))
ThePlot <- ggplot(average_spot_size_df[average_spot_size_df$treatment != "Untx",], aes(x = condition, y = spot_area_scaled_to_100um))
ThePlot <- ThePlot + geom_quasirandom(fill = "#F2FDFF", aes(group = label), colour = "#000000", shape = 21, dodge.width = 0.9, alpha = 0.6)
ThePlot <- ThePlot + geom_boxplot(aes(fill = label), colour = "#000000", 
    width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9))
ThePlot <- ThePlot + scale_x_discrete("condition")
ThePlot <- ThePlot + scale_y_continuous("average (median) spot size scaled to 100 um")
ThePlot <- ThePlot + scale_colour_manual("condition", values = condition_colours)
ThePlot <- ThePlot + scale_fill_manual("condition", values = condition_colours)
# ThePlot <- ThePlot + scale_colour_manual("status", values = Colours)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ treatment, scales = "free")
ThePlot <- ThePlot + ggtitle("CDF of mean total spot area per cell")
ggsave(width = 5, height = 10, file = file.path(analysis_plots_dir, "05-02-plot-sampled-spot_area_scaled_to_100um.pdf"), ThePlot)


condition_colours <- c("control_Scramble" = "#a6cee3", "control_ERAP2" = "#1f78b4", "vcp_Scramble" = "#fb9a99", "vcp_ERAP2" =  "#e31a1c")
average_spot_size_df$label <- factor(average_spot_size_df$label, levels = names(condition_colours))
ThePlot <- ggplot(average_spot_size_df[average_spot_size_df$treatment != "Untx",], aes(x = condition, y = median_relative_spot_size))
ThePlot <- ThePlot + geom_quasirandom(fill = "#F2FDFF", aes(group = label), colour = "#000000", shape = 21, dodge.width = 0.9, alpha = 0.6)
ThePlot <- ThePlot + geom_boxplot(aes(fill = label), colour = "#000000", 
    width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9))
ThePlot <- ThePlot + scale_x_discrete("condition")
ThePlot <- ThePlot + scale_y_continuous("median relative spot size")
ThePlot <- ThePlot + scale_colour_manual("condition", values = condition_colours)
ThePlot <- ThePlot + scale_fill_manual("condition", values = condition_colours)
# ThePlot <- ThePlot + scale_colour_manual("status", values = Colours)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ treatment, scales = "free")
ThePlot <- ThePlot + ggtitle("CDF of mean total spot area per cell")
ggsave(width = 5, height = 10, file = file.path(analysis_plots_dir, "05-03-plot-sampled-relative-spot-size.pdf"), ThePlot)




# Plot video of true mean

require("gganimate")
require("av")

condition_colours <- c("control_Scramble" = "#a6cee3", "control_ERAP2" = "#1f78b4", "vcp_Scramble" = "#fb9a99", "vcp_ERAP2" =  "#e31a1c")
condition_colours <- c("control" = "#1f78b4", "vcp" =  "#e31a1c")
average_spot_size_df$label <- factor(average_spot_size_df$label, levels = names(condition_colours))
vline_df <- do.call(rbind, lapply(split(all_summarised_cell_df, paste(all_summarised_cell_df$condition, all_summarised_cell_df$treatment)), function(x){
    data.frame(condition = unique(x$condition),
        treatment = unique(x$treatment),
        spot_area_scaled_to_100um = median(x$spot_area_scaled_to_100um),
        animate_spot_area_scaled_to_100um = average_spot_size_df$spot_area_scaled_to_100um,
        relative_spot_size = median(x$relative_spot_size),
        animate_relative_spot_size = average_spot_size_df$median_relative_spot_size)
}))

ThePlot <- ggplot(all_summarised_cell_df[all_summarised_cell_df$treatment != "Untx",], aes(x = log10(spot_area_scaled_to_100um)))
ThePlot <- ThePlot + geom_density(aes(fill = condition), colour = "#000000", alpha = 0.6)
ThePlot <- ThePlot + geom_vline(data = vline_df[vline_df$treatment != "Untx",], aes(xintercept = log10(spot_area_scaled_to_100um)), colour = "#000000")
ThePlot <- ThePlot + transition_states(states = log10(animate_spot_area_scaled_to_100um), transition_length = 2, state_length = 1) + enter_drift() + exit_drift() + ease_aes("elastic-in-out")
ThePlot <- ThePlot + scale_x_continuous("estimated spot size scaled to 100um")
ThePlot <- ThePlot + scale_y_continuous("density")
ThePlot <- ThePlot + scale_colour_manual("condition", values = condition_colours)
ThePlot <- ThePlot + scale_fill_manual("condition", values = condition_colours)
# ThePlot <- ThePlot + scale_colour_manual("status", values = Colours)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + facet_grid(condition ~ treatment, scales = "free")
ThePlot <- ThePlot + ggtitle("density of spot area scaled to 100um")
anim_plot <- animate(ThePlot, width = 1525, height = 972, fps = 5,
    renderer = av_renderer())
anim_save(file.path(analysis_plots_dir, "plot-transition-of-mean.mp4"), anim_plot)

# ggsave(width = 10, height = 10, file = file.path(analysis_plots_dir, "05-04-plot-density-distribution-of-spot-area-scaled-to-100um.pdf"), ThePlot)








