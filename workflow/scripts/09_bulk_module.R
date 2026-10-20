# Adapted analysis; see provenance/source_ranges.json.
source(here::here("workflow", "scripts", "common.R"))
library(stringr)
library(ggbeeswarm)
module_df <- readRDS(input("bulk_modules.rds"))
erap2_module <- unique(module_df$modules[!is.na(module_df$gene_name) & module_df$gene_name == "ERAP2"])
stopifnot(length(erap2_module) == 1)
enrichment <- read.csv(input("bulk_module_go_annotated.csv"))
writeLines(erap2_module, file.path(plot_dir, "ERAP2_module.txt"))
module_col <- intersect(c("modules", "module"), names(enrichment))
stopifnot(length(module_col) == 1)
enrichment <- enrichment[enrichment[[module_col]] == erap2_module,]
p_col <- intersect(c("Adjusted.P.value", "adj_pval"), names(enrichment))[1]
stopifnot(!is.na(p_col))
enrichment <- enrichment[!is.na(enrichment[[p_col]]) & enrichment[[p_col]] < 0.05,]

label_col <- intersect(c("Group", "group", "Term"), names(enrichment))[1]
stopifnot(nrow(enrichment) > 0, !is.na(label_col))
enrichment$label <- enrichment[[label_col]]
enrichment$label <- reorder(enrichment$label, -log10(enrichment$P.value), FUN = median)
p <- ggplot(enrichment, aes(-log10(P.value), label)) + geom_quasirandom(shape = 21, fill = "#F4AC45") +
  my_theme + labs(x = "-log10 p-value", y = "ERAP2 module GO:BP", title = erap2_module)
ggsave(file.path(plot_dir, "Figure_3C.pdf"), p, width = 8, height = 7)
write.csv(enrichment, file.path(plot_dir, "Figure_3C_terms.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(plot_dir, "sessionInfo.txt"))
