# Adapted analysis; see provenance/source_ranges.json.
source(here::here("workflow", "scripts", "common.R"))
significant_events_df <- readRDS(input("significant_proteins.rds"))
frequency_df <- as.data.frame(table(significant_events_df$comparison, significant_events_df$direction))
colnames(frequency_df) <- c("comparison", "direction", "frequency")
frequency_df <- frequency_df[!(frequency_df$comparison %in% c("VCP_ERAP2_vs_CTRL_ERAP2")),]
frequency_df$comparison <- factor(frequency_df$comparison, levels = c("VCP_CTRL_vs_CTRL_CTRL", "CTRL_ERAP2_vs_CTRL_CTRL", "VCP_ERAP2_vs_VCP_CTRL"))
prop_df <- do.call(rbind, lapply(split(frequency_df, frequency_df$comparison), function(x){
    x <- x[rev(order(x$direction)),]
    x$prop <- x$frequency / sum(x$frequency)
    x$lab_y <-  cumsum(x$prop) - (0.5 * x$prop)
    return(x)
}))

require("ggplot2")

my_theme <- theme(legend.position = "bottom",
    text = element_text(size = 12),
    axis.text.x = element_text(size = 8, hjust = 1, angle = 45),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "azure3",
    linetype = "dashed", size = 0.2))

Colours <- c("up" = "#D81E5B", "down" = "#3993DD")
the_plot <- ggplot(prop_df, aes(x = comparison, y = prop))
the_plot <- the_plot + geom_col(aes(fill = direction), colour = "#000000")
the_plot <- the_plot + geom_text(data = prop_df, aes(x = comparison, y = lab_y, label = frequency), colour = "#000000")
the_plot <- the_plot + scale_x_discrete("comparisons")
the_plot <- the_plot + scale_y_continuous("Proportion of differentially expressed proteins")
the_plot <- the_plot + scale_fill_manual("lfc direction", values = Colours)
the_plot <- the_plot + theme_bw()
the_plot <- the_plot + my_theme
ggsave(the_plot, file = file.path(plot_dir, "02-06-frequency-of-differential-protein-enrichment.pdf"), height = 9, width = 3)

writeLines(capture.output(sessionInfo()), file.path(plot_dir, "sessionInfo.txt"))
