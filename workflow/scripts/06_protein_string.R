# Adapted analysis; see provenance/source_ranges.json.
source(here::here("workflow", "scripts", "common.R"))
library(ggbeeswarm)
significant_enrichment_df <- read.csv(here::here("workflow", "resources", "protein_string_grouped.csv"))
compartment_df <- subset(significant_enrichment_df, fdr < 0.05 & category == "COMPARTMENTS")
keyword_df <- subset(significant_enrichment_df, fdr < 0.05 & category == "Keyword")

require(ggbeeswarm)

my_theme <- theme(
    legend.position = "bottom",
    text = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "azure3",
    linetype = "dashed",
    size = 0.2))

max_gene_ratio <- vapply(split(compartment_df$gene_ratio, compartment_df$description), max, 1)

Colours <- c("down->up" = "#FB3640", "up->down" = "#247BA0",
    "down->down" = "#72A3C0", "up->up" = "#E0777D")

compartment_df$description <- factor(compartment_df$description, levels = names(max_gene_ratio)[order(max_gene_ratio, decreasing = F)])

ThePlot <- ggplot(compartment_df, aes(x = gene_ratio, y = description))
ThePlot <- ThePlot + geom_vline(xintercept = 0.1, linetype = "dashed")
ThePlot <- ThePlot + geom_quasirandom(aes(fill = direction_group), colour = "#000000", shape = 21, size = 1.5)
ThePlot <- ThePlot + scale_x_continuous("gene ratio")
ThePlot <- ThePlot + scale_y_discrete("Enriched compartments (fdr < 0.05)")
ThePlot <- ThePlot + scale_fill_manual("Direction groups", values = Colours)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("Compartment enrichment from STRING")
ThePlot <- ThePlot + my_theme
ggsave(width = 10, height = 10, file = file.path(plot_dir, "03-06-plot-enriched-compartments-from-STRING.pdf"), ThePlot)

my_theme <- theme(
    legend.position = "bottom",
    text = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "azure3",
    linetype = "dashed",
    size = 0.2))

max_gene_ratio <- vapply(split(keyword_df$gene_ratio, keyword_df$description), max, 1)

Colours <- c("down->up" = "#FB3640", "up->down" = "#247BA0",
    "down->down" = "#72A3C0", "up->up" = "#E0777D")

keyword_df$description <- factor(keyword_df$description, levels = names(max_gene_ratio)[order(max_gene_ratio, decreasing = F)])

ThePlot <- ggplot(keyword_df, aes(x = gene_ratio, y = description))
ThePlot <- ThePlot + geom_vline(xintercept = 0.1, linetype = "dashed")
ThePlot <- ThePlot + geom_quasirandom(aes(fill = direction_group), colour = "#000000", shape = 21, size = 1.5)
ThePlot <- ThePlot + scale_x_continuous("gene ratio")
ThePlot <- ThePlot + scale_y_discrete("Enriched compartments (fdr < 0.05)")
ThePlot <- ThePlot + scale_fill_manual("Direction groups", values = Colours)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("Compartment enrichment from STRING")
ThePlot <- ThePlot + my_theme
ggsave(width = 5, height = 10, file = file.path(plot_dir, "03-07-plot-enriched-keywords-from-STRING.pdf"), ThePlot)

writeLines(capture.output(sessionInfo()), file.path(plot_dir, "sessionInfo.txt"))
