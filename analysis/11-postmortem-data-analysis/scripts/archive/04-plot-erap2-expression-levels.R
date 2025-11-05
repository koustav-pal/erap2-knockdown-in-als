data_dir <- normalizePath("/camp/lab/patanir/home/users/palk/")
project_dir <- file.path(data_dir, "Projects", "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "11-postmortem-data-analysis/")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
erap2_kd_de_df <- readRDS(file.path(analysis_save_dir, "01-all-contrasts-differential-splicing-dataframe.rds"))
temp_dir <- file.path(analysis_dir, "tmp")
alsod_genes_df <- read.csv(file.path(project_dir, "data_files/alsod_genelist.csv"))
annotation_gtf <- "/nemo/lab/patanir/home/users/palk/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
umap_clustering_list <- readRDS(file.path(analysis_save_dir, "01-patient-latent-factor-umap-clustering.rds"))
genic_count_list <- readRDS(file.path(analysis_save_dir, "01-gene-level-counts-of-postmortem-samples.rds"))
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
source("scripts/src/functions.R")
# source("../src/rmats-analysis/src/source.R")
# ==========================================================================
# Analysis
# ==========================================================================
annotation_ranges <- import.gff2(annotation_gtf)
annotation_gene_ranges <- annotation_ranges[annotation_ranges$type == "gene"]


coldata_df <- genic_count_list$coldata


clustering_df <- umap_clustering_list$umap
patient_matrix <- umap_clustering_list$counts

clustering_df$subject_id <- str_split(clustering_df$patients, pattern = " ", simplify = T)[,1]
clustering_df$sex_genotype <- coldata_df$sex_genotype[match(clustering_df$subject_id, coldata_df$subject_id)]
clustering_df$age_at_symptom_onset <- coldata_df$age_at_symptom_onset[match(clustering_df$subject_id, coldata_df$subject_id)]
clustering_df$age_at_death <- coldata_df$age_at_death[match(clustering_df$subject_id, coldata_df$subject_id)]
clustering_df$disease_duration_in_months <- coldata_df$disease_duration_in_months[match(clustering_df$subject_id, coldata_df$subject_id)]
clustering_df$group <- NA
clustering_df$group[clustering_df$cluster %in% c(1,2,5)] <- "cortex"
clustering_df$group[clustering_df$cluster %in% c(4,3)] <- "spinal_cord"
clustering_df$group[clustering_df$cluster %in% c(6,7)] <- "cerebellum"


cols_split_by_group <- split(clustering_df$patients, clustering_df$group)

patients_expr_zscore_matrix <- do.call(cbind, lapply(cols_split_by_group, function(x){
    current_matrix <- patient_matrix[,x]
    zscore_matrix <- (current_matrix - rowMeans(current_matrix)) / rowSds(current_matrix)
    return(zscore_matrix)
}))

clustering_df$erap2_zscore <- patients_expr_zscore_matrix["ENSG00000164308",clustering_df$patients]

clustering_df <- clustering_df[clustering_df$type != "none",]
clustering_df$type <- factor(clustering_df$type, levels = c("ctrl", "c9orf72", "tardbp", "mutation", "sporadic"))


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



colours <- c("ctrl" = "#EFF7F6", "als" = "#F87575", "als+ftd" = "#FFA9A3")
ThePlot <- ggplot(clustering_df, aes(x = type, y = erap2_zscore))
ThePlot <- ThePlot + geom_violin(aes(fill = als_group), colour = "#000000", 
    lwd = 0.5, position = position_dodge(width = 0.9), scale = "width")
ThePlot <- ThePlot + geom_boxplot(aes(fill = als_group), colour = "#000000", 
    width = 0.1, lwd = 0.5, position = position_dodge(width = 0.9))
ThePlot <- ThePlot + scale_x_discrete("mutation groups")
ThePlot <- ThePlot + scale_y_continuous("erap2 expression zscore")
ThePlot <- ThePlot + scale_fill_manual("condition", values = colours)
# ThePlot <- ThePlot + scale_colour_distiller("overall expr", palette = "RdBu", direction = -1)
ThePlot <- ThePlot + theme_bw()
# ThePlot <- ThePlot + facet_wrap(. ~ als_group, nrow = 1)
ThePlot <- ThePlot + ggtitle("expression landscape across mutations")
ThePlot <- ThePlot + current_theme
ggsave(height = 7, width = 7, file = file.path(analysis_plots_dir, "05-08-plot-erap2-expression-in-mutations.pdf"), ThePlot)


require("dplyr")
frequency_df <- do.call(rbind, lapply(unique(clustering_df$als_group), function(x){
    temp_df <- clustering_df[clustering_df$als_group == x,]
    temp_df$type <- droplevels(temp_df$type)
    frequency_table <- table(temp_df$type, temp_df$tissue)
    proportion_table <- frequency_table / rowSums(frequency_table)
    freq_df <- as.data.frame(frequency_table)
    prop_df <- as.data.frame(proportion_table)
    colnames(freq_df) <- c("type", "tissue", "frequency")
    freq_df$proportion <- prop_df[,3]
    freq_df$als_group <- x
    # return_df <- freq_df %>% group_by(type) %>% arrange(tissue, desc(proportion)) %>% mutate(prop_y = cumsum(proportion) - 0.5 * proportion)
    # freq_df$proportion <- freq_df$frequency / sum(freq_df$frequency)
    return(return_df)
}))

frequency_df <- frequency_df %>% 
group_by(als_group, type) %>% 
arrange(tissue, desc(proportion*-1)) %>% 
mutate(prop_y = cumsum(proportion) - 0.5 * proportion) %>% 
as.data.frame()

colours <- c("ctrl" = "#EFF7F6", "als" = "#F87575", "als+ftd" = "#FFA9A3")
ThePlot <- ggplot(frequency_df, aes(x = type, y = frequency, label = frequency))
ThePlot <- ThePlot + geom_bar(aes(fill = tissue), colour = "#000000", stat = "identity")
ThePlot <- ThePlot + geom_text(size = 3, position = position_stack(vjust = 0.5))
ThePlot <- ThePlot + scale_x_discrete("mutation groups")
ThePlot <- ThePlot + scale_y_continuous("proportion of patients")
ThePlot <- ThePlot + scale_fill_brewer("tissue", values = "Set1")
# ThePlot <- ThePlot + scale_colour_distiller("overall expr", palette = "RdBu", direction = -1)
ThePlot <- ThePlot + theme_bw()
ThePlot <- ThePlot + facet_wrap(. ~ als_group, scales = "free_x", nrow = 1)
ThePlot <- ThePlot + ggtitle("expression landscape across mutations")
ThePlot <- ThePlot + current_theme
ggsave(height = 7, width = 7, file = file.path(analysis_plots_dir, "05-08-plot-patient-proportions.pdf"), ThePlot)