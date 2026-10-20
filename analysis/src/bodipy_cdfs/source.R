.spot_analysR <- setClass("spot_analysR", 
    slots = representation(
        experiment_block = "numeric",
        variables = "character",
        order = "character",
        field_col = "character",
        cell_col = "character",
        summarised_columns = "character"),
    contains="SummarizedExperiment")

setGeneric("fetch_data", function(object) standardGeneric("fetch_data"))
setMethod("fetch_data", "spot_analysR", function(object) {
    return(assay(object))
})

setGeneric("fetch_metadata", function(object, col_name) standardGeneric("fetch_metadata"))
setMethod("fetch_metadata", "spot_analysR", function(object, col_name) {
    if(!(col_name %in% names(object@metadata))){
        stop(paste(col_name, "not found in object..."))
    }
    return(object@metadata[[col_name]])
})

setGeneric("fetch_cell_count", function(object) standardGeneric("fetch_cell_count"))
setMethod("fetch_cell_count", "spot_analysR", function(object) {
    return(fetch_metadata(object, "cell_matrix"))
})

setGeneric("count_cells", function(object) standardGeneric("count_cells"))
setMethod("count_cells", "spot_analysR", function(object) {
    sample_labels <- fetch_assay_labels(object)
    sampleinfo_df <- fetch_sampleinfo(object)
    experiment_labels <- fetch_experiment_labels(object)
    assay_data <- fetch_data(object)
    cell_df <- do.call(rbind, lapply(seq_along(experiment_labels), function(x){
        current_sample <- experiment_labels[x]
        current_df <- assay_data[sample_labels == current_sample,]
        per_field_per_cell_df <- split(current_df, paste(current_df[,object@field_col], current_df[,object@cell_col]))
        return_rows_df <- sampleinfo_df[experiment_labels == current_sample,]
        return_rows_df$sample <- current_sample
        return_rows_df$num_cells <- length(per_field_per_cell_df)
        return(return_rows_df)
    }))
    object@metadata[["cell_matrix"]] <- cell_df
    return(object)
})


.make_labels <- function(data_df, vars){
    labels <- vapply(seq_len(nrow(data_df[,vars])), function(x){
        paste(data_df[x,vars], collapse = "_")  
    },"1")
    return(labels)
}

spot_analysR <- function(data_df, experiment_block, field_col, cell_col, vars, var_order, ...){
    sampleinfo_df <- unique(data_df[,vars])
    assay_labels <- .make_labels(data_df, vars)
    experiment_labels <- .make_labels(sampleinfo_df, vars)
    summarized_experiment <- SummarizedExperiment(
    assays = list(data = data_df),
    metadata = list(data = sampleinfo_df, 
        assay_labels = assay_labels,
        experiment_labels = experiment_labels))
    return_object <- .spot_analysR(summarized_experiment)
    return_object@experiment_block <- experiment_block
    return_object@variables <- vars
    return_object@order <- var_order
    return_object@field_col <- field_col
    return_object@cell_col <- cell_col
    object_with_cells_n <- count_cells(return_object)
    return(object_with_cells_n)
}

assign_labels_to_columns <- function(x, experiment_df, map_col, value.var){
    if(length(map_col) > 1){
        stop("map_col cannot have length greater than 1")
    }
    if(!(map_col %in% colnames(x)) | !(map_col %in% colnames(experiment_df))){
        stop(paste(map_col, " not found in x or experiment_df"))
    }
    if(any(!(value.var %in% colnames(experiment_df)))){
        stop(paste(paste(value.var[which(!(value.var %in% colnames(experiment_df)))], collapse = ","), " not found in experiment_df"))
    }
    unique_expt_df <- unique((experiment_df[,c(map_col, value.var)]))
    if(any(duplicated(unique_expt_df[,map_col]))){
        stop(paste("Found duplicated entries for", map_col, "in", value.var, "column of experiment_df...\n", map_col, "cannot map to multiple entries..."))
    }
    x[,value.var] <- unique_expt_df[match(x[,map_col], unique_expt_df[,map_col]),value.var]
    return(x)
}

setGeneric("fetch_sampleinfo", function(object) standardGeneric("fetch_sampleinfo"))
setMethod("fetch_sampleinfo", "spot_analysR", function(object) {
    return(fetch_metadata(object, "data"))
})
setGeneric("fetch_assay_labels", function(object) standardGeneric("fetch_assay_labels"))
setMethod("fetch_assay_labels", "spot_analysR", function(object) {
    return(fetch_metadata(object, "assay_labels"))
})

setGeneric("fetch_experiment_labels", function(object) standardGeneric("fetch_experiment_labels"))
setMethod("fetch_experiment_labels", "spot_analysR", function(object) {
    return(fetch_metadata(object, "experiment_labels"))
})

setGeneric("fetch_summarised_data", function(object, col_name) standardGeneric("fetch_summarised_data"))
setMethod("fetch_summarised_data", "spot_analysR", function(object, col_name) {
    if(!(col_name %in% object@summarised_columns)){
        stop(paste(col_name, "not found in object..."))
    }
    return(fetch_metadata(object, col_name))
})

# setGeneric("plot_data_density", function(object, summarise_col, col_name) standardGeneric("plot_data_density"))

setGeneric("summarise_data", function(object, summarise_col, col_name) standardGeneric("summarise_data"))
setMethod("summarise_data", "spot_analysR", function(object, summarise_col, col_name) {
    sample_labels <- fetch_assay_labels(object)
    sampleinfo_df <- fetch_sampleinfo(object)
    experiment_labels <- fetch_experiment_labels(object)
    assay_data <- fetch_data(object)
    summary_df <- do.call(rbind, lapply(seq_along(experiment_labels), function(x){
        current_sample <- experiment_labels[x]
        current_df <- assay_data[sample_labels == current_sample,]
        per_field_per_cell_df <- split(current_df, paste(current_df[,object@field_col], current_df[,object@cell_col]))
        mean_val <- as.vector(vapply(per_field_per_cell_df, function(y){
            mean(y[,summarise_col])}, 1))
        std_dev <- as.vector(vapply(per_field_per_cell_df, function(y){
            sd(y[,summarise_col])}, 1))
        median_val <- as.vector(vapply(per_field_per_cell_df, function(y){
            median(y[,summarise_col])}, 1))
        total <- as.vector(vapply(per_field_per_cell_df, function(y){
            sum(y[,summarise_col])}, 1))
        max <- as.vector(vapply(per_field_per_cell_df, function(y){
            max(y[,summarise_col])}, 1))
        n <- as.vector(vapply(per_field_per_cell_df, function(y){
            length(y[,summarise_col])}, 1))
        return_rows_df <- sampleinfo_df[experiment_labels == current_sample,]
        return_rows_df <- return_rows_df[rep(1, length(mean_val)),]
        return_rows_df$sample <- current_sample
        return_rows_df$mean <- mean_val
        return_rows_df$std_dev <- std_dev
        return_rows_df$median <- median_val
        return_rows_df$total <- total
        return_rows_df$max <- max
        return_rows_df$n <- n
        return(return_rows_df)
    }))
    object@metadata[[col_name]] <- summary_df
    object@summarised_columns <- c(object@summarised_columns, col_name)
    return(object)
})

setGeneric("prep_params", function(object, select_cells = 100, iteration_size = 10000, duplication_threshold = 1) standardGeneric("prep_params"))
setMethod("prep_params", "spot_analysR", function(object, select_cells = 100, iteration_size = 10000, duplication_threshold = 1) {
    object@metadata[["params"]] <- c(total_cells_to_select = select_cells, 
        num_iter = iteration_size,
        prop_dup = duplication_threshold)
    return(object)
})

setGeneric("fetch_params", function(object) standardGeneric("fetch_params"))
setMethod("fetch_params", "spot_analysR", function(object) {
    return(fetch_metadata(object, "params"))
})

setGeneric("generate_seeds", function(object) standardGeneric("generate_seeds"))
setMethod("generate_seeds", "spot_analysR", function(object) {
    params <- fetch_params(object)
    experiment_labels <- fetch_experiment_labels(object)
    sample_seed_matrix <- do.call(cbind, lapply(seq_along(experiment_labels), function(x){
        seed_list <- sample(x = seq(1, as.numeric(Sys.time())), size = params["num_iter"])
        return(seed_list)
    }))
    colnames(sample_seed_matrix) <- experiment_labels
    object@metadata[["seed_list"]] <- sample_seed_matrix
    return(object)
})

setGeneric("fetch_seeds", function(object) standardGeneric("fetch_seeds"))
setMethod("fetch_seeds", "spot_analysR", function(object) {
    fetch_metadata(object, "seed_list")
})

setGeneric("sample_cells_across_replicates", function(object, summarise_by = "replicate") standardGeneric("sample_cells_across_replicates"))
setMethod("sample_cells_across_replicates", "spot_analysR", function(object, summarise_by = "replicate") {
    params <- fetch_params(object)
    summarised_cols <- object@summarised_columns

    all_summarised_data <- object@metadata[summarised_cols]

    sampleinfo_df <- fetch_sampleinfo(object)
    all_vars <- object@variables
    summary_vars <- all_vars[!(all_vars %in% summarise_by)]

    experiment_labels <- fetch_experiment_labels(object)
    summarised_experiment_labels <- .make_labels(sampleinfo_df, summary_vars)

    split_experiment_labels <- split(experiment_labels, summarised_experiment_labels)

    seed_matrix <- fetch_seeds(object)
    cell_count <- fetch_cell_count(object)

    sampled_cell_list <- lapply(all_summarised_data, function(summary_df){
        all_sampled_cells_df <- do.call(rbind, lapply(names(split_experiment_labels), function(x){
            current_samples <- split_experiment_labels[[x]]

            current_summary_sample_df <- summary_df[summary_df$sample %in% current_samples,]
            current_summary_labels <- current_summary_sample_df$sample

            current_cell_count <- cell_count[cell_count$sample %in% current_samples,]
            current_cell_count$prop <- current_cell_count$num_cells / sum(current_cell_count$num_cells)
            if(params["total_cells_to_select"] > sum(current_cell_count$num_cells)){
                warning("Skipping ",x," due to low number of cells...")
                return(NULL)
            }
            current_summary_df_split <- split(current_summary_sample_df, current_summary_labels)

            selected_cell_df_list <- lapply(names(current_summary_df_split), function(y){
                temp_df <- current_summary_df_split[[y]]
                select_cells <- round(current_cell_count$prop[current_cell_count$sample == y] * params["total_cells_to_select"])
                sequence_over_cells <- seq_len(nrow(temp_df))

                cell_list <- lapply(seq_along(seed_matrix[,y]), function(y2){
                    set.seed(seed_matrix[y2,y])
                    sample(x = sequence_over_cells, size = select_cells)
                })

                dup_matrix <- do.call(rbind, lapply(cell_list, function(y2){
                    sequence_over_cells %in% y2
                }))

                prop_values_duplicated <- vapply(seq_len(nrow(dup_matrix))[-1], function(y2){
                    sum(colSums(dup_matrix[seq(1,y2),]) > 1) / current_cell_count$num_cells[current_cell_count$sample == y]
                },1)
                
                prop_values_duplicated <- c(0, prop_values_duplicated)
                selected_num_iter <- min(which(prop_values_duplicated == 1))

                sampled_iter_df <- do.call(rbind, lapply(seq(1, selected_num_iter), function(y2){
                    return_df <- temp_df[cell_list[[y2]],]
                    return_df$iter <- y2
                    return(return_df)
                }))
                return(sampled_iter_df)
            })
            selected_cell_df <- do.call(rbind, selected_cell_df_list)
            selected_cell_df$summary_group <- x
            return(selected_cell_df)
        }))
        return(all_sampled_cells_df)
    })
    object@metadata[["sampled_cells"]] <- sampled_cell_list
    return(object)
})


cdf_fun <- function(a_vector, quantiles){
    proportions <- vapply(quantiles, function(y){
        length(a_vector[a_vector <= y]) / length(a_vector)
    },1)
    return_df <- data.frame(quantile = quantiles, cdf = proportions)
    return(return_df)
}

setGeneric("get_data_range", function(object, summarise_by = "replicate") standardGeneric("get_data_range"))
setMethod("get_data_range", "spot_analysR", function(object){
    summarised_cols <- object@summarised_columns
    range_list <- lapply(summarised_cols, function(x){
        fun_list <- c("mean", "median", "total", "max", "n")
        min_list <- vapply(fun_list, function(y){
            min(object@metadata[[x]][,y])
        },1)
        max_list <- vapply(fun_list, function(y){
            max(object@metadata[[x]][,y])
        },1)
        cbind(min_list, max_list)
    })
    names(range_list) <- summarised_cols
    object@metadata[["value_ranges"]] <- range_list
    return(object)
})

setGeneric("get_data_quantiles", function(object, fends = 30) standardGeneric("get_data_quantiles"))
setMethod("get_data_quantiles", "spot_analysR", function(object, fends = 30){
    summarised_cols <- object@summarised_columns
    quantile_matrix <- lapply(summarised_cols, function(x){
        fun_list <- c("mean", "median", "total", "max", "n")
        quantile_list <- do.call(rbind, lapply(fun_list, function(y){
            quantile(object@metadata[[x]][,y], seq(0,1,1/fends))
        }))
        rownames(quantile_list) <- fun_list
        return(quantile_list)
    })
    names(quantile_matrix) <- summarised_cols
    object@metadata[["value_quantiles"]] <- quantile_matrix
    return(object)
})

setGeneric("compute_cdfs", function(object, summarise_by = "replicate") standardGeneric("compute_cdfs"))
setMethod("compute_cdfs", "spot_analysR", function(object, summarise_by = "replicate") {
    summarised_cols <- object@summarised_columns
    all_summarised_data <- object@metadata[summarised_cols]

    # data_value_range_list <- fetch_metadata(object, "value_ranges")
    data_value_quantile_list <- fetch_metadata(object, "value_quantiles")

    # sampled_cell_list <- fetch_metadata(object, "sampled_cells")
    fun_list <- c("mean", "median", "total", "max", "n")

    all_vars <- object@variables
    summary_vars <- all_vars[!(all_vars %in% summarise_by)]

    all_sampled_cdf_df <- lapply(summarised_cols, function(summary_col){
        summary_df <- sampled_cell_list[[summary_col]]
        rownames(summary_df) <- NULL
        current_quantiles <- data_value_quantile_list[[summary_col]]
        summary_df_split <- split(summary_df, summary_df$summary_group)
        cdf_df <- lapply(summary_df_split, function(temp_df){
            temp_df_iter_split <- split(temp_df, temp_df$iter)
            iter_df <- do.call(rbind, lapply(temp_df_iter_split, function(temp_df_2){
                fun_df <- do.call(rbind, lapply(fun_list, function(x){
                    return_df <- cdf_fun(temp_df_2[,x], current_quantiles[x,])
                    return_df$q <- seq_len(ncol(current_quantiles))
                    return_df$fun <- x
                    return(return_df)
                }))
                fun_df$iter <- unique(temp_df_2$iter)
                return(fun_df)
            }))

            iter_df_split <- split(iter_df, iter_df$fun)
            n <- length(iter_df_split)

            mean_values <- lapply(iter_df_split, function(temp_iter_df){
                fun_mat <- reshape2::acast(data = temp_iter_df, formula = iter ~ q, value.var = "cdf")
                colMeans(fun_mat)
            })
            col_sd <- lapply(iter_df_split, function(temp_iter_df){
                fun_mat <- reshape2::acast(data = temp_iter_df, formula = iter ~ q, value.var = "cdf")
                colSds(fun_mat)
            })

            cdf_df <- do.call(rbind, lapply(fun_list, function(x){
                current_fun_means <- mean_values[[x]]
                current_sds <- col_sd[[x]]
                margin <- qt(0.975, df = n - 1) * current_sds / sqrt(n)
                return_df <- data.frame(
                    fun = x,
                    quantiles = current_quantiles[x,], 
                    quantiles_scaled = current_quantiles[x,] / max(current_quantiles[x,]),
                    q = seq_len(ncol(current_quantiles)),
                    cdf = current_fun_means,
                    margin = margin,
                    lower = current_fun_means - margin,
                    upper = current_fun_means + margin)
                return(return_df)
            }))
            rownames(cdf_df) <- NULL
            return(cbind(unique(temp_df[,summary_vars]), cdf_df))
        })
        do.call(rbind, cdf_df)
    })
    names(all_sampled_cdf_df) <- names(sampled_cell_list)
    object@metadata[["cdf_function"]] <- all_sampled_cdf_df
    return(object)
})

return_basic_theme <- function(type = "basic"){
    if(type == "bw"){
        current_theme <- theme_bw()
    }
    if(type == "basic"){
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
    }

    return(current_theme)
}

setGeneric("plot_correlations", function(object, file, variables, cor.method = "pearson", summarise = T) standardGeneric("plot_correlations"))
setMethod("plot_correlations", "spot_analysR", function(object, file, variables, cor.method = "pearson", summarise = T) {
    seq_over_vars <- seq_along(variables)
    seq_rep_mat <- cbind(rep(seq_over_vars, each = length(seq_over_vars)),
        rep(seq_over_vars, times = length(seq_over_vars)))

    sample_labels <- fetch_assay_labels(object)
    sampleinfo_df <- fetch_sampleinfo(object)
    experiment_labels <- fetch_experiment_labels(object)
    assay_data <- fetch_data(object)
    assay_vars <- object@variables
    return_df <- unique(assay_data[,assay_vars])
    all_cor_df <- do.call(rbind, lapply(seq_along(experiment_labels), function(x){
        current_sample <- experiment_labels[x]
        current_df <- assay_data[sample_labels == current_sample,]
        cor_df <- do.call(rbind, lapply(seq_len(nrow(seq_rep_mat)), function(x){
            cor_val <- cor(current_df[,variables[seq_rep_mat[x,1]]], current_df[,variables[seq_rep_mat[x,2]]], method = cor.method)
            cbind(return_df, data.frame(x = variables[seq_rep_mat[x,1]], 
                            y = variables[seq_rep_mat[x,2]],
                            cor = cor_val))
        }))
        cor_df$sample <- current_sample
        return(cor_df)
    }))
    
    plot_df <- all_cor_df
    if(summarise){
        summary_df_split <- split(all_cor_df, paste(all_cor_df$x, all_cor_df$y, sep = "_"))
        summary_cor_df <- do.call(rbind, lapply(summary_df_split, function(temp_df){
            return_df <- unique(temp_df[,c("x", "y")])
            return_df$cor <- mean(temp_df$cor)
            return(return_df)
        }))
        plot_df <- summary_cor_df
    }
    if(!is.null(names(variables))){
        plot_df$x <- names(variables)[match(plot_df$x, variables)]
        plot_df$y <- names(variables)[match(plot_df$y, variables)]
    }

    current_theme <- return_basic_theme(type = "bw")
    plot_df$cor[plot_df$x == plot_df$y] <- NA
    MinMax <- abs(min(c(max(plot_df$cor, na.rm = T), min(plot_df$cor, na.rm = T))))
    plot_df$scaled_cor <- plot_df$cor
    plot_df$scaled_cor[plot_df$scaled_cor > MinMax] <- MinMax
    plot_df$scaled_cor[plot_df$scaled_cor < (MinMax*-1)] <- MinMax*-1
    ThePlot <- ggplot(plot_df, aes(x = x, y = y))
    ThePlot <- ThePlot + geom_tile(aes(fill = scaled_cor))
    ThePlot <- ThePlot + geom_text(aes(label = round(scaled_cor,2)))
    ThePlot <- ThePlot + scale_x_discrete("variables")
    ThePlot <- ThePlot + scale_y_discrete("variables")
    ThePlot <- ThePlot + current_theme
    ThePlot <- ThePlot + scale_fill_distiller(palette = "RdBu", direction = 1)
    if(!summarise){
        ThePlot <- ThePlot + facet_wrap(. ~ sample)
    }
    ThePlot <- ThePlot + ggtitle("correlation between variables")
    ggsave(width = 7, height = 7, file = file, ThePlot)
})

setGeneric("plot_model_cell_count_vs_var", function(object, file, var, FUN = sum, colour_by = NA) standardGeneric("plot_model_cell_count_vs_var"))
setMethod("plot_model_cell_count_vs_var", "spot_analysR", function(object, file, var, FUN = sum, colour_by = NA) {
    sample_labels <- fetch_assay_labels(object)
    sampleinfo_df <- fetch_sampleinfo(object)
    experiment_labels <- fetch_experiment_labels(object)
    assay_data <- fetch_data(object)
    cell_count_df <- fetch_cell_count(object)
    all_summarised_df <- do.call(rbind, lapply(seq_along(experiment_labels), function(x){
        current_sample <- experiment_labels[x]
        current_df <- assay_data[sample_labels == current_sample,]
        return_df <- cell_count_df[cell_count_df[,"sample"] == current_sample,]
        return_df[,"plot_var"] <- FUN(current_df[,var])
        return(return_df)
    }))
    if(!is.na(colour_by)){
        all_summarised_df$colour_var <- all_summarised_df[,colour_by]
    }

    current_theme <- return_basic_theme(type = "bw")
    ThePlot <- ggplot(all_summarised_df, aes(x = num_cells, y = plot_var))
    if(!is.na(colour_by)){
        ThePlot <- ThePlot + geom_point(aes(fill = colour_var), shape = 21, colour = "#000000")
        ThePlot <- ThePlot + geom_smooth(method = "lm", aes(colour = colour_var), formula = y ~ x)
        ThePlot <- ThePlot + scale_fill_brewer(palette = "Set1")
        ThePlot <- ThePlot + scale_colour_brewer(palette = "Set1")
    }else{
        ThePlot <- ThePlot + geom_point(fill = "#F6F4F3", shape = 21, colour = "#000000")
        ThePlot <- ThePlot + geom_smooth(method = "lm", formula = y ~ x)
    }
    
    ThePlot <- ThePlot + scale_x_continuous("N of cells in experiment")
    ThePlot <- ThePlot + scale_y_continuous(var)
    ThePlot <- ThePlot + current_theme
    ThePlot <- ThePlot + ggtitle(paste("Association between number of cells and ", var))
    ggsave(width = 7, height = 7, file = file, ThePlot)
})

setGeneric("plot_violin_per_sample", function(object, file, var, log = T, colour_by = NA, width = 7, height = 7) standardGeneric("plot_violin_per_sample"))
setMethod("plot_violin_per_sample", "spot_analysR", function(object, file, var, log = T, colour_by = NA, width = 7, height = 7) {
    sample_labels <- fetch_assay_labels(object)
    sampleinfo_df <- fetch_sampleinfo(object)
    experiment_labels <- fetch_experiment_labels(object)
    assay_data <- fetch_data(object)
    assay_data$plot_var <- assay_data[,var]
    assay_data$sample <- sample_labels
    
    if(!is.na(colour_by)){
        assay_data$colour_var <- assay_data[,colour_by]
    }
    if(log){
        assay_data$plot_var <- log10(assay_data$plot_var)
    }

    current_theme <- return_basic_theme(type = "bw")
    ThePlot <- ggplot(assay_data, aes(x = sample, y = plot_var))
    if(!is.na(colour_by)){
        ThePlot <- ThePlot + geom_violin(aes(fill = colour_var), colour = "#000000")
        ThePlot <- ThePlot + geom_boxplot(aes(fill = colour_var), colour = "#000000", 
            width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9))
        ThePlot <- ThePlot + scale_fill_brewer(palette = "Set1")
        ThePlot <- ThePlot + scale_colour_brewer(palette = "Set1")
    }else{
        ThePlot <- ThePlot + geom_violin(fill = "#F6F4F3", colour = "#000000")
        ThePlot <- ThePlot + geom_boxplot(fill = "#F6F4F3", colour = "#000000", 
            width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9))
    }    
    ThePlot <- ThePlot + scale_x_discrete("Sample")
    ThePlot <- ThePlot + scale_y_continuous(var)
    ThePlot <- ThePlot + current_theme
    ThePlot <- ThePlot + ggtitle(paste("Distribution of", var))
    ggsave(width = width, height = height, file = file, ThePlot)
})

setGeneric("plot_model_var1_vs_var2", function(object, file, var1, var2, log = T, colour_by = NA) standardGeneric("plot_model_var1_vs_var2"))
setMethod("plot_model_var1_vs_var2", "spot_analysR", function(object, file, var1, var2, log = T, colour_by = NA) {
    sample_labels <- fetch_assay_labels(object)
    sampleinfo_df <- fetch_sampleinfo(object)
    experiment_labels <- fetch_experiment_labels(object)
    assay_data <- fetch_data(object)
    assay_data[,"plot_var1"] <- assay_data[,var1]
    assay_data[,"plot_var2"] <- assay_data[,var2]

    if(log){
        assay_data[,"plot_var1"] <- log10(assay_data[,var1])
        assay_data[,"plot_var2"] <- log10(assay_data[,var2])
    }
    if(!is.na(colour_by)){
        assay_data$colour_var <- assay_data[,colour_by]
    }

    current_theme <- return_basic_theme(type = "bw")

    ThePlot <- ggplot(assay_data, aes(x = plot_var1, y = plot_var2))
    if(!is.na(colour_by)){
        ThePlot <- ThePlot + geom_point(aes(fill = colour_var), shape = 21, colour = "#000000")
        ThePlot <- ThePlot + geom_smooth(method = "lm", aes(colour = colour_var), formula = y ~ x)
        ThePlot <- ThePlot + scale_fill_brewer(palette = "Set1")
        ThePlot <- ThePlot + scale_colour_brewer(palette = "Set1")
    }else{
        ThePlot <- ThePlot + geom_point(fill = "#F6F4F3", shape = 21, colour = "#000000")
        ThePlot <- ThePlot + geom_smooth(method = "lm", formula = y ~ x)
    }

    ThePlot <- ThePlot + geom_smooth(method = "lm", formula = y ~ x)
    ThePlot <- ThePlot + scale_x_continuous(var1)
    ThePlot <- ThePlot + scale_y_continuous(var2)
    ThePlot <- ThePlot + current_theme
    ThePlot <- ThePlot + ggtitle(paste("Association between",var1,"and",var2))
    ggsave(width = 7, height = 7, file = file, ThePlot)
})


setGeneric("fetch_spot_distances", function(object, x_col = "X", y_col = "Y") standardGeneric("fetch_spot_distances"))
setMethod("fetch_spot_distances", "spot_analysR", function(object, x_col = "X", y_col = "Y") {
    sample_labels <- fetch_assay_labels(object)
    sampleinfo_df <- fetch_sampleinfo(object)
    experiment_labels <- fetch_experiment_labels(object)
    assay_data <- fetch_data(object)
    summary_df <- do.call(rbind, lapply(seq_along(experiment_labels), function(x){
        current_sample <- experiment_labels[x]
        current_df <- assay_data[sample_labels == current_sample,]
        per_field_per_cell_df <- split(current_df, paste(current_df[,object@field_col], current_df[,object@cell_col]))
        distances <- vapply(per_field_per_cell_df, function(temp_df){
            x_values <- temp_df[,x_col]
            y_values <- temp_df[,y_col]
            if(length(x_values) == 1){
                return(NA)
            }
            polar_coordscart2pol(x_values, y_values, degrees = F)
            pair_matrix <- cbind(x_values, y_values)
            
            seq_over_val <- seq_len(nrow(pair_matrix))
            
            x_rep <- rep(seq_over_val, each = length(pair_matrix))
            y_rep <- rep(seq_over_val, times = length(seq_over_val))

            distance <- sqrt((pair_matrix[x_rep, 1] - pair_matrix[y_rep, 1])^2 + (pair_matrix[x_rep, 2] - pair_matrix[y_rep, 2])^2)



            }, 1)

        return_rows_df <- sampleinfo_df[experiment_labels == current_sample,]
        return_rows_df <- return_rows_df[rep(1, length(mean_val)),]
        return_rows_df$sample <- current_sample
        return_rows_df$mean <- mean_val
        return_rows_df$std_dev <- std_dev
        return_rows_df$median <- median_val
        return_rows_df$total <- total
        return_rows_df$max <- max
        return_rows_df$n <- n
        return(return_rows_df)
    }))
    object@metadata[[col_name]] <- summary_df
    object@summarised_columns <- c(object@summarised_columns, col_name)
    return(object)
})





summarise_cdf_by_col <- function(cdf_df, fun, q_col = "quantiles", col_names){
    fun_cdf_df <- cdf_df[cdf_df$fun == fun,]
    all_labels <- .make_labels(fun_cdf_df, col_names)
    fun_cdf_df$label <- all_labels
    split_cdf_df <- split(fun_cdf_df, fun_cdf_df$label)
    use_lab_cols <- colnames(cdf_df)
    use_lab_cols <- use_lab_cols[!(use_lab_cols %in% c(col_names, "fun", "quantiles", "quantiles_scaled", "q", "cdf", "margin", "lower", "upper"))]
    return_df <- do.call(rbind, lapply(split_cdf_df, function(temp_df){
        temp_df$xaxis <- apply(temp_df[,use_lab_cols], 1, function(x){paste(x, collapse = "_") })
        a_cell_matrix <- (reshape2::acast(data = temp_df, formula = xaxis ~ q, value.var = "cdf"))
        quantile_matrix <- (reshape2::acast(data = temp_df, formula = xaxis ~ q, value.var = "quantiles"))
        mean_values <- colMeans(a_cell_matrix)
        std_dev <- colSds(a_cell_matrix)
        n <- nrow(a_cell_matrix) 
        margin <- qt(0.975, df = n - 1) * std_dev / sqrt(n)
        cbind(unique(temp_df[,col_names]), data.frame(
            q = as.numeric(colnames(a_cell_matrix)),
            quantiles = colMeans(quantile_matrix),
            lower = mean_values - margin, 
            mean = mean_values,
            upper = mean_values + margin))
    }))
    rownames(return_df) <- NULL
    return(return_df)
}

plot_data_density <- function(object, plot_col, log = F, q = c(0.1, 0.25, 0.5, 0.75, 0.9), return_data = F, file) {
    sample_labels <- fetch_assay_labels(object)
    sampleinfo_df <- fetch_sampleinfo(object)
    experiment_labels <- fetch_experiment_labels(object)
    assay_data <- fetch_data(object)
    summary_df <- do.call(rbind, lapply(seq_along(experiment_labels), function(x){
        current_sample <- experiment_labels[x]
        current_df <- assay_data[sample_labels == current_sample,]
        return_df <- data.frame(cell = current_df[,object@cell_col], 
            field = current_df[,object@field_col], 
            sample = current_sample, 
            x = current_df[,plot_col])
        return(return_df)
    }))

    value.var <- "x"
    if(log){
        summary_df$logT <- log10(summary_df$x)
        value.var <- "logT"
    }

    density_prof <- density(summary_df[,value.var])
    density_df <- data.frame(x = density_prof$x, y = density_prof$y)

    quantiles <- quantile(summary_df[,value.var], q)
    quantile_df <- data.frame(q = quantiles, y = max(density_prof$y)/2, labs = names(quantiles))

    current_theme <- return_basic_theme(type = "bw")

    ThePlot <- ggplot(density_df, aes(x = x, y = y))
    ThePlot <- ThePlot + geom_area(fill = "#F6F4F3")
    ThePlot <- ThePlot + geom_line(colour = "#000000")
    ThePlot <- ThePlot + geom_vline(data = quantile_df, aes(xintercept = q), lty = "dashed")
    ThePlot <- ThePlot + geom_text(data = quantile_df, aes(x = q, y = y, label = labs))
    ThePlot <- ThePlot + scale_x_continuous("log10 spot area")
    ThePlot <- ThePlot + scale_y_continuous("density")
    ThePlot <- ThePlot + current_theme
    ThePlot <- ThePlot + ggtitle("spot area density")
    ggsave(width = 7, height = 7, file = file, ThePlot)
    if(return_data){
        return(summary_df)
    }
}