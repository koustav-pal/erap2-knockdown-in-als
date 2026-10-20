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
bodypy_df <- read.csv(file.path(analysis_input_dir, "ERAP2_AC_18_03_24.csv"))
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


bodypy_df$sample_name <- paste(bodypy_df[,1], bodypy_df[,3], sep = "_")
bodypy_df$Treatment <- stringr::str_split(string = bodypy_df$Treatment, pattern = " ", simplify = T)[,1]
adjusted_mean_df <- do.call(rbind, lapply(split(bodypy_df, bodypy_df$sample_name), function(x){
    error <- sd(x[,"raw.spots.per.area.cyto"]) / sqrt(nrow(x))
    return_df <- unique(x[,-c(4,5,6)])
    return_df$adjusted_value <- mean(x[,"raw.spots.per.area.cyto"]) / error
    return(return_df)
}))

keep_df <- adjusted_mean_df[adjusted_mean_df$Line != "CTRL3",]


bodypy_df$Genotype <- factor(bodypy_df$Genotype, levels = c("CTRL", "VCP"))
bodypy_df$Treatment <- factor(bodypy_df$Treatment, levels = c("Scr", "ERAP2"))

summary(lm(data = bodypy_df[bodypy_df$Genotype != "VCP",], formula = raw.spots.per.area.cyto ~ Treatment))


wilcox.test(bodypy_df$Fold.change.spots.per.area.of.cyto[bodypy_df$Genotype == "VCP" & 
    bodypy_df$Treatment == "Scr siRNA"],
    bodypy_df$Fold.change.spots.per.area.of.cyto[bodypy_df$Genotype == "CTRL" & 
    bodypy_df$Treatment == "Scr siRNA"])


wilcox.test(bodypy_df$Fold.change.spots.per.area.of.cyto[bodypy_df$Genotype == "VCP" & 
    bodypy_df$Treatment == "Scr siRNA"],
    bodypy_df$Fold.change.spots.per.area.of.cyto[bodypy_df$Genotype == "VCP" & 
    bodypy_df$Treatment == "ERAP2  siRNA"])



summary(lm(bodypy_df, formula = Fold.change.spots.per.area.of.cyto ~ Genotype + Genotype:Treatment))



