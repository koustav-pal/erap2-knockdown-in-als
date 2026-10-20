# Adapted analysis; see provenance/source_ranges.json.
source(here::here("workflow", "scripts", "common.R"))
library(stringr)
library(fgsea)
go_bp_lines <- readLines(input("GO_Biological_Process_2023.txt"))
go_bp_list <- str_split(go_bp_lines, pattern = "\t")
go_bp_list_names <- vapply(go_bp_list, function(x){
	x[1]
},"1")
names(go_bp_list) <- go_bp_list_names
go_bp_list_final <- lapply(go_bp_list, function(x){
	x[-c(1,2,length(x))]
})

proteomics_lfc_df <- read.csv(input("protein_background.csv"))
protein_columns <- c("vcp_ctrl_vs_ctrl_ctrl", "vcp_erap2_vs_vcp_ctrl")
all_enrichment_results <- lapply(protein_columns, function(x){
	logfc_values <- proteomics_lfc_df[,paste(x, "logfc", sep = "_")]
	names(logfc_values) <- proteomics_lfc_df$protein_name
	logfc_values <- logfc_values[order(logfc_values, decreasing = T)]
	enrichment_df <- fgsea(pathways = go_bp_list_final,
		stats = logfc_values,
		minSize = 10,
        nproc = 1,
		eps = 0.0,
		nPermSimple = 10000)
	enrichment_df$contrast <- x
	return(enrichment_df)
})
names(all_enrichment_results) <- protein_columns

significant_pathways <- do.call(c, lapply(all_enrichment_results[c("vcp_ctrl_vs_ctrl_ctrl", "vcp_erap2_vs_vcp_ctrl")], function(temp_df){
	temp_df$pathway[temp_df$pval < 0.1]
}))

pathway_frequency <- table(significant_pathways)
present_in_both <- names(pathway_frequency[pathway_frequency == 2])

selected_enrichment_matrix <- do.call(cbind, lapply(all_enrichment_results[c("vcp_ctrl_vs_ctrl_ctrl", "vcp_erap2_vs_vcp_ctrl")], function(temp_df){
	current_df <- temp_df[match(present_in_both, temp_df$pathway),]
	a_vector <- current_df$NES
	return(a_vector)
}))
selected_enrichment_df <- as.data.frame(selected_enrichment_matrix)
selected_enrichment_df$pathway <- present_in_both

fisher.test(table(selected_enrichment_df$vcp_ctrl_vs_ctrl_ctrl < 0, selected_enrichment_df$vcp_erap2_vs_vcp_ctrl > 0))

require("ggplot2")
model <- lm(data = selected_enrichment_df,
    formula = vcp_erap2_vs_vcp_ctrl ~ vcp_ctrl_vs_ctrl_ctrl)
summary(model)

my_theme <- theme(
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

ThePlot <- ggplot(selected_enrichment_df, aes(x = vcp_ctrl_vs_ctrl_ctrl, y = vcp_erap2_vs_vcp_ctrl))
ThePlot <- ThePlot + geom_hline(yintercept = 0, linetype = "dashed")
ThePlot <- ThePlot + geom_vline(xintercept = 0, linetype = "dashed")
ThePlot <- ThePlot + geom_point(colour = "#FFC43D")
ThePlot <- ThePlot + geom_smooth(method = "lm", formula = y ~ x, lty = "dashed", colour = "#000000")
ThePlot <- ThePlot + scale_x_continuous("VCP vs Control effect")
ThePlot <- ThePlot + scale_y_continuous("Effect of ERAP2 KD in VCP")

ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("VCP vs Control effect versus ERAP2 KD in VCP")
ThePlot <- ThePlot + my_theme
ggsave(width = 10, height = 10, file = file.path(analysis_plots_dir, "01-03-plot-scatterplot-of-gsea-nes-scores-for-proteins.pdf"), ThePlot)

write.csv(selected_enrichment_df, file.path(plot_dir, "gsea_NES.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(plot_dir, "sessionInfo.txt"))
