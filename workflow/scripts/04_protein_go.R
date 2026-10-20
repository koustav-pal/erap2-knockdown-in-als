# Adapted analysis; see provenance/source_ranges.json.
source(here::here("workflow", "scripts", "common.R"))
library(stringr)
library(ggbeeswarm)
library(patchwork)
annotated_gene_set_enrichment_df <- read.csv(input("protein_go_annotated.csv"))

upregulated_terms_df <- annotated_gene_set_enrichment_df[annotated_gene_set_enrichment_df$contrasts == "VCP_CTRL_vs_CTRL_CTRL" & annotated_gene_set_enrichment_df$direction == "up",]

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

median_gene_ratio <- vapply(split(upregulated_terms_df$gene_ratio, upregulated_terms_df$Group), length, 1)
upregulated_terms_df$Group <- factor(upregulated_terms_df$Group,
    levels = names(median_gene_ratio[order(median_gene_ratio, decreasing = F)]))
keep_levels <- levels(upregulated_terms_df$Group)
keep_levels <- keep_levels[keep_levels != "Other Processes"]

ThePlot_1 <- ggplot(upregulated_terms_df[upregulated_terms_df$Group %in% rev(keep_levels)[1:10],], aes(x = -log10(P.value), y = Group))
ThePlot_1 <- ThePlot_1 + geom_quasirandom(aes(fill = direction), colour = "#000000",
     position = position_jitter(), shape = 21)
ThePlot_1 <- ThePlot_1 + scale_x_continuous("-log10 pvalue")
ThePlot_1 <- ThePlot_1 + scale_y_discrete("Upregulated terms from GO:BP")
ThePlot_1 <- ThePlot_1 + scale_fill_manual("direction", values = c("up" = "#E0777D", "down" = "#4C86A8"))
ThePlot_1 <- ThePlot_1 + theme_bw()
ThePlot_1 <- ThePlot_1 + ggtitle("Upregulated genes")
ThePlot_1 <- ThePlot_1 + my_theme

downregulated_terms_df <- annotated_gene_set_enrichment_df[annotated_gene_set_enrichment_df$contrasts == "VCP_CTRL_vs_CTRL_CTRL" & annotated_gene_set_enrichment_df$direction == "down",]
median_gene_ratio <- vapply(split(downregulated_terms_df$gene_ratio, downregulated_terms_df$Group), length, 1)
downregulated_terms_df$Group <- factor(downregulated_terms_df$Group,
    levels = names(median_gene_ratio[order(median_gene_ratio, decreasing = F)]))
keep_levels <- levels(downregulated_terms_df$Group)
keep_levels <- keep_levels[keep_levels != "Other Processes"]

ThePlot_2 <- ggplot(downregulated_terms_df[downregulated_terms_df$Group %in% rev(keep_levels)[1:10],], aes(x = -log10(P.value), y = Group))
ThePlot_2 <- ThePlot_2 + geom_quasirandom(aes(fill = direction), colour = "#000000", shape = 21)
ThePlot_2 <- ThePlot_2 + scale_x_continuous("-log10 pvalue")
ThePlot_2 <- ThePlot_2 + scale_y_discrete("Downregulated terms from GO:BP")
ThePlot_2 <- ThePlot_2 + scale_fill_manual("direction", values = c("up" = "#E0777D", "down" = "#4C86A8"))
ThePlot_2 <- ThePlot_2 + theme_bw()
ThePlot_2 <- ThePlot_2 + ggtitle("Downregulated genes")
ThePlot_2 <- ThePlot_2 + my_theme

require("patchwork")
ThePlot <- ThePlot_1 | ThePlot_2

ggsave(width = 10, height = 12, file = file.path(plot_dir, "03-04-plot-gene-set-enrichment-using-gobp.pdf"), ThePlot)

annotated_gene_set_enrichment_df <- read.csv(input("pairwise_go_annotated.csv"))

Colours <- c("down->up" = "#FB3640", "up->down" = "#247BA0",
    "down->down" = "#72A3C0", "up->up" = "#E0777D")

down_up_terms_df <- annotated_gene_set_enrichment_df[annotated_gene_set_enrichment_df$direction_group == "down->up",]

median_gene_ratio <- vapply(split(down_up_terms_df$gene_ratio, down_up_terms_df$Group), length, 1)
down_up_terms_df$Group <- factor(down_up_terms_df$Group,
    levels = names(median_gene_ratio[order(median_gene_ratio, decreasing = F)]))
keep_levels <- levels(down_up_terms_df$Group)
keep_levels <- keep_levels[keep_levels != "Other Processes"]

ThePlot_1 <- ggplot(down_up_terms_df[down_up_terms_df$Group %in% rev(keep_levels)[1:10],],
    aes(x = -log10(P.value), y = Group))
ThePlot_1 <- ThePlot_1 + geom_quasirandom(aes(fill = direction_group), colour = "#000000",
     position = position_jitter(), shape = 21)
ThePlot_1 <- ThePlot_1 + scale_x_continuous("-log10 pvalue")
ThePlot_1 <- ThePlot_1 + scale_y_discrete("Downregulated -> upregulated terms from GO:BP")
ThePlot_1 <- ThePlot_1 + scale_fill_manual("direction", values = Colours)
ThePlot_1 <- ThePlot_1 + theme_bw()
ThePlot_1 <- ThePlot_1 + ggtitle("Down -> up proteins")
ThePlot_1 <- ThePlot_1 + my_theme

up_down_terms_df <- annotated_gene_set_enrichment_df[annotated_gene_set_enrichment_df$direction_group == "up->down",]

median_gene_ratio <- vapply(split(up_down_terms_df$gene_ratio, up_down_terms_df$Group), length, 1)
up_down_terms_df$Group <- factor(up_down_terms_df$Group,
    levels = names(median_gene_ratio[order(median_gene_ratio, decreasing = F)]))
keep_levels <- levels(up_down_terms_df$Group)
keep_levels <- keep_levels[keep_levels != "Other Processes"]

ThePlot_2 <- ggplot(up_down_terms_df[up_down_terms_df$Group %in% rev(keep_levels)[1:10],], aes(x = -log10(P.value), y = Group))
ThePlot_2 <- ThePlot_2 + geom_quasirandom(aes(fill = direction_group), colour = "#000000", shape = 21)
ThePlot_2 <- ThePlot_2 + scale_x_continuous("-log10 pvalue")
ThePlot_2 <- ThePlot_2 + scale_y_discrete("Upregulated -> downregulated terms from GO:BP")
ThePlot_2 <- ThePlot_2 + scale_fill_manual("direction", values = Colours)
ThePlot_2 <- ThePlot_2 + theme_bw()
ThePlot_2 <- ThePlot_2 + ggtitle("Up -> down proteins")
ThePlot_2 <- ThePlot_2 + my_theme

down_down_terms_df <- annotated_gene_set_enrichment_df[annotated_gene_set_enrichment_df$direction_group == "down->down",]

median_gene_ratio <- vapply(split(down_down_terms_df$gene_ratio, down_down_terms_df$Group), length, 1)
down_down_terms_df$Group <- factor(down_down_terms_df$Group,
    levels = names(median_gene_ratio[order(median_gene_ratio, decreasing = F)]))
keep_levels <- levels(down_down_terms_df$Group)
keep_levels <- keep_levels[keep_levels != "Other Processes"]

ThePlot_3 <- ggplot(down_down_terms_df[down_down_terms_df$Group %in% rev(keep_levels)[1:10],], aes(x = -log10(P.value), y = Group))
ThePlot_3 <- ThePlot_3 + geom_quasirandom(aes(fill = direction_group), colour = "#000000", shape = 21)
ThePlot_3 <- ThePlot_3 + scale_x_continuous("-log10 pvalue")
ThePlot_3 <- ThePlot_3 + scale_y_discrete("Downregulated -> downregulated terms from GO:BP")
ThePlot_3 <- ThePlot_3 + scale_fill_manual("direction", values = Colours)
ThePlot_3 <- ThePlot_3 + theme_bw()
ThePlot_3 <- ThePlot_3 + ggtitle("Down -> down proteins")
ThePlot_3 <- ThePlot_3 + my_theme

up_up_terms_df <- annotated_gene_set_enrichment_df[annotated_gene_set_enrichment_df$direction_group == "up->up",]

median_gene_ratio <- vapply(split(up_up_terms_df$gene_ratio, up_up_terms_df$Group), length, 1)
up_up_terms_df$Group <- factor(up_up_terms_df$Group,
    levels = names(median_gene_ratio[order(median_gene_ratio, decreasing = F)]))
keep_levels <- levels(up_up_terms_df$Group)
keep_levels <- keep_levels[keep_levels != "Other Processes"]

ThePlot_4 <- ggplot(up_up_terms_df[up_up_terms_df$Group %in% rev(keep_levels)[1:10],], aes(x = -log10(P.value), y = Group))
ThePlot_4 <- ThePlot_4 + geom_quasirandom(aes(fill = direction_group), colour = "#000000", shape = 21)
ThePlot_4 <- ThePlot_4 + scale_x_continuous("-log10 pvalue")
ThePlot_4 <- ThePlot_4 + scale_y_discrete("Upregulated -> upregulated terms from GO:BP")
ThePlot_4 <- ThePlot_4 + scale_fill_manual("direction", values = Colours)
ThePlot_4 <- ThePlot_4 + theme_bw()
ThePlot_4 <- ThePlot_4 + ggtitle("Up -> up proteins")
ThePlot_4 <- ThePlot_4 + my_theme

require("patchwork")
main_plot_1 <- ThePlot_1 | ThePlot_4
main_plot_2 <- ThePlot_3 | ThePlot_2
main_plot <- main_plot_1 / main_plot_2

ggsave(width = 12, height = 12, file = file.path(plot_dir, "03-05-plot-gene-set-enrichment-using-gobp-for-ameliorated-and-abrogated-proteins.pdf"), main_plot)

ggsave(file.path(plot_dir, "Figure_1D.pdf"), ThePlot_1, width = 6, height = 6)
ggsave(file.path(plot_dir, "Figure_1E.pdf"), ThePlot_2, width = 6, height = 6)
ggsave(file.path(plot_dir, "Figure_S2C_down.pdf"), ThePlot_3, width = 6, height = 6)
ggsave(file.path(plot_dir, "Figure_S2C_up.pdf"), ThePlot_4, width = 6, height = 6)

writeLines(capture.output(sessionInfo()), file.path(plot_dir, "sessionInfo.txt"))
