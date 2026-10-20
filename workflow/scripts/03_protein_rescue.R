# Adapted analysis; see provenance/source_ranges.json.
source(here::here("workflow", "scripts", "common.R"))
pairwise_group_sig_df <- readRDS(input("pairwise_proteins.rds"))
background_groups_df <- read.csv(input("protein_background.csv"))
significantly_affected_proteins_df <- pairwise_group_sig_df[pairwise_group_sig_df$vcp_ctrl_vs_ctrl_ctrl_sig,]

Colours <- c("down->up" = "#FB3640", "up->down" = "#247BA0",
    "down->down" = "#72A3C0", "up->up" = "#E0777D")

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

ThePlot <- ggplot(background_groups_df[!(background_groups_df$protein_name %in% significantly_affected_proteins_df$protein_name),], aes(x = vcp_ctrl_vs_ctrl_ctrl_logfc, y = vcp_erap2_vs_vcp_ctrl_logfc))
ThePlot <- ThePlot + geom_hline(yintercept = 0, linetype = "dashed")
ThePlot <- ThePlot + geom_vline(xintercept = 0, linetype = "dashed")
ThePlot <- ThePlot + geom_point(colour = "azure3", alpha = 0.4)
ThePlot <- ThePlot + geom_smooth(data = background_groups_df, method = "lm", formula = y ~ x, lty = "dashed", colour = "#000000")
ThePlot <- ThePlot + geom_point(data = significantly_affected_proteins_df,
    aes(x = vcp_ctrl_vs_ctrl_ctrl_logfc, y = vcp_erap2_vs_vcp_ctrl_logfc, colour = direction_group))
ThePlot <- ThePlot + scale_x_continuous("VCP vs Control effect")
ThePlot <- ThePlot + scale_y_continuous("Effect of ERAP2 KD in VCP")
ThePlot <- ThePlot + scale_colour_manual("direction groups", values = Colours)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + ggtitle("VCP vs Control effect versus ERAP2 KD in VCP")
ThePlot <- ThePlot + my_theme
ggsave(width = 10, height = 10, file = file.path(plot_dir, "03-05-plot-vcp-vs-control-versus-effect-of-erap2-kd-in-vcp.pdf"), ThePlot)

writeLines(capture.output(sessionInfo()), file.path(plot_dir, "sessionInfo.txt"))
