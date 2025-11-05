
cols_split_by_patients <- split(colnames(keep_matrix),paste(
        coldata_df$subject_id[match(colnames(keep_matrix), coldata_df$column_id)],
        coldata_df$tissue[match(colnames(keep_matrix), coldata_df$column_id)]))

patient_matrix <- do.call(cbind, lapply(cols_split_by_patients, function(x){
    if(length(x) == 1){
        return(keep_matrix[,x])
    }
    return(rowMeans(keep_matrix[,x]))
}))
patient_matrix <- patient_matrix[genes_in_same_module,]
patient_erap2_expr <- patient_matrix["ENSG00000164308",]
patient_coldata <- data.frame(column_id = colnames(patient_matrix),
    patient_id = str_split(colnames(patient_matrix), pattern = " ", simplify = T)[,1],
    tissue = str_split(colnames(patient_matrix), pattern = " ", simplify = T)[,2],
    erap2_expr = patient_erap2_expr[colnames(patient_matrix)])


pca_object <- prcomp(patient_matrix, scale = TRUE)
pca_summary <- summary(pca_object)
loadings <- pca_object$rotation
loadings_df <- as.data.frame(loadings[,c("PC1", "PC2", "PC3", "PC4", "PC5")])
loadings_df$erap2_expr <- patient_coldata$erap2_expr[match(rownames(loadings_df), patient_coldata$column_id)]
loadings_df$patient_id <- patient_coldata$patient_id[match(rownames(loadings_df), patient_coldata$column_id)]
loadings_df$mnd_with_ftd <- coldata_df$mnd_with_ftd[match(loadings_df$patient_id, coldata_df$subject_id)]
loadings_df$mnd_with_dementia <- coldata_df$mnd_with_dementia[match(loadings_df$patient_id, coldata_df$subject_id)]
loadings_df$tissue <- patient_coldata$tissue[match(rownames(loadings_df), patient_coldata$column_id)]
loadings_df$type <- coldata_df$type[match(loadings_df$patient_id, coldata_df$subject_id)]
loadings_df$condition <- coldata_df$condition[match(loadings_df$patient_id, coldata_df$subject_id)]


pca_object <- prcomp(patient_matrix, scale = TRUE)
pca_summary <- summary(pca_object)
variance_df <- data.frame(pc = colnames(pca_summary$importance), 
    variance_explained = pca_summary$importance[2,])

require(ggplot2)
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


ThePlot <- ggplot(loadings_df, aes(x = PC1, y = PC2, colour = tissue))
# ThePlot <- ThePlot + geom_vline(xintercept = 0, linetype = "dashed", colour = "#000000")
# ThePlot <- ThePlot + geom_hline(yintercept = 0, linetype = "dashed", colour = "#000000")
ThePlot <- ThePlot + geom_point()
ThePlot <- ThePlot + scale_colour_brewer("genotypes", palette = "Set1")
# ThePlot <- ThePlot + scale_colour_distiller("erap2 expr", palette = "Reds", direction = 1)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ tissue)
ThePlot <- ThePlot + ggtitle("principal component analysis of erap2 associated genes")
ThePlot <- ThePlot + current_theme
ggsave(height = 7, width = 7, file = file.path(analysis_plots_dir, "04-plot-pca-of-erap2-correlated-genes-in-postmortem.pdf"), ThePlot)




latent_factors <- loadings[,1:10]
umap_object <- umap(latent_factors)
graph <- bluster::makeSNNGraph(latent_factors)
clust <- igraph::cluster_louvain(graph)


plot_data <- data.frame(patients = rownames(umap_object$layout), 
    xaxis = umap_object$layout[,1], 
    yaxis = umap_object$layout[,2],
    cluster = factor(clust$membership))
plot_data$erap2_expr <- patient_coldata$erap2_expr[match(rownames(plot_data), patient_coldata$column_id)]
plot_data$patient_id <- patient_coldata$patient_id[match(rownames(plot_data), patient_coldata$column_id)]
plot_data$mnd_with_ftd <- coldata_df$mnd_with_ftd[match(plot_data$patient_id, coldata_df$subject_id)]
plot_data$mnd_with_dementia <- coldata_df$mnd_with_dementia[match(plot_data$patient_id, coldata_df$subject_id)]
plot_data$tissue <- patient_coldata$tissue[match(rownames(plot_data), patient_coldata$column_id)]
plot_data$type <- coldata_df$type[match(plot_data$patient_id, coldata_df$subject_id)]
plot_data$condition <- coldata_df$condition[match(plot_data$patient_id, coldata_df$subject_id)]


require(ggplot2)
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

plot_data$plot_erap_expr <- plot_data$erap2_expr
plot_data$plot_erap_expr <- log10(plot_data$plot_erap_expr)
ThePlot <- ggplot(plot_data, aes(x = xaxis, y = yaxis, colour = tissue))
# ThePlot <- ThePlot + geom_vline(xintercept = 0, linetype = "dashed", colour = "#000000")
# ThePlot <- ThePlot + geom_hline(yintercept = 0, linetype = "dashed", colour = "#000000")
ThePlot <- ThePlot + geom_point()
ThePlot <- ThePlot + scale_x_continuous("umap 1")
ThePlot <- ThePlot + scale_y_continuous("umap 2")
ThePlot <- ThePlot + scale_colour_brewer("genotypes", palette = "Set1")
# ThePlot <- ThePlot + scale_colour_distiller("erap2 expr", palette = "Reds", direction = 1)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ tissue)
ThePlot <- ThePlot + ggtitle("umap analysis of erap2 associated genes")
ThePlot <- ThePlot + current_theme
ggsave(height = 7, width = 7, file = file.path(analysis_plots_dir, "05-02-plot-umap-of-erap2-associated-factors-highlight-tissue.pdf"), ThePlot)


ThePlot <- ggplot(plot_data, aes(x = xaxis, y = yaxis, colour = as.character(cluster)))
# ThePlot <- ThePlot + geom_vline(xintercept = 0, linetype = "dashed", colour = "#000000")
# ThePlot <- ThePlot + geom_hline(yintercept = 0, linetype = "dashed", colour = "#000000")
ThePlot <- ThePlot + geom_point()
ThePlot <- ThePlot + scale_x_continuous("umap 1")
ThePlot <- ThePlot + scale_y_continuous("umap 2")
ThePlot <- ThePlot + scale_colour_brewer("genotypes", palette = "Set1")
# ThePlot <- ThePlot + scale_colour_distiller("erap2 expr", palette = "Reds", direction = 1)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ tissue)
ThePlot <- ThePlot + ggtitle("umap analysis of erap2 associated genes")
ThePlot <- ThePlot + current_theme
ggsave(height = 7, width = 7, file = file.path(analysis_plots_dir, "05-02-plot-umap-of-erap2-associated-factors-highlight-clusters.pdf"), ThePlot)


ThePlot <- ggplot(plot_data, aes(x = xaxis, y = yaxis, colour = as.character(condition)))
# ThePlot <- ThePlot + geom_vline(xintercept = 0, linetype = "dashed", colour = "#000000")
# ThePlot <- ThePlot + geom_hline(yintercept = 0, linetype = "dashed", colour = "#000000")
ThePlot <- ThePlot + geom_point()
ThePlot <- ThePlot + scale_x_continuous("umap 1")
ThePlot <- ThePlot + scale_y_continuous("umap 2")
ThePlot <- ThePlot + scale_colour_brewer("genotypes", palette = "Set1")
# ThePlot <- ThePlot + scale_colour_distiller("erap2 expr", palette = "Reds", direction = 1)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ tissue)
ThePlot <- ThePlot + ggtitle("umap analysis of erap2 associated genes")
ThePlot <- ThePlot + current_theme
ggsave(height = 7, width = 7, file = file.path(analysis_plots_dir, "05-02-plot-umap-of-erap2-associated-factors-highlight-condition.pdf"), ThePlot)


ThePlot <- ggplot(plot_data, aes(x = xaxis, y = yaxis, colour = plot_erap_expr))
# ThePlot <- ThePlot + geom_vline(xintercept = 0, linetype = "dashed", colour = "#000000")
# ThePlot <- ThePlot + geom_hline(yintercept = 0, linetype = "dashed", colour = "#000000")
ThePlot <- ThePlot + geom_point()
ThePlot <- ThePlot + scale_x_continuous("umap 1")
ThePlot <- ThePlot + scale_y_continuous("umap 2")
# ThePlot <- ThePlot + scale_colour_brewer("genotypes", palette = "Set1")
ThePlot <- ThePlot + scale_colour_distiller("erap2 expr", palette = "Reds", direction = 1)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_grid(. ~ tissue)
ThePlot <- ThePlot + ggtitle("umap analysis of erap2 associated genes")
ThePlot <- ThePlot + current_theme
ggsave(height = 7, width = 7, file = file.path(analysis_plots_dir, "05-02-plot-umap-of-erap2-associated-factors-highlight-erap2-expression.pdf"), ThePlot)

saveRDS(object = list(umap = plot_data, counts = patient_matrix), file = file.path(analysis_save_dir, "01-patient-latent-factor-umap-clustering.rds"))