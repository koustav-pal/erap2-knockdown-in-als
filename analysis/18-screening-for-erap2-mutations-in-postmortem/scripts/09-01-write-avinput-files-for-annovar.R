data_dir <- normalizePath("/export2/home/koustav/")
project_dir <- file.path(data_dir, "Projects", "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "18-screening-for-erap2-mutations-in-postmortem")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_output_dir <- file.path(analysis_dir, "output_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
temp_dir <- file.path(analysis_dir, "tmp")
tmp_dir <- file.path(analysis_dir, "slurm")
patient_count_df <- readRDS(file.path("../11-postmortem-data-analysis/rdata_files/02-summarised-counts-for-patients.rds"))
mutation_screening_results <- file.path(analysis_output_dir, "FILTERED_MUTATION_SCREENING_RESULTS.vcf")
# ==========================================================================
# Packages
# ==========================================================================
setwd("/export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/analysis")
require("HiCBricks")
require("GenomicRanges")
require("stringr")
# require("Biostrings")
# require("rtracklayer")
# require("GenomeInfoDb")
# require("GenomicFeatures")
require("ggplot2")
library(rlang)
library(survival)
library(reshape2)
library(knitr)
library(rmarkdown)
require("RColorBrewer")
require("survminer")
require("survival")
require("randomForestSRC")
require("parallel")
require("matrixStats")
require("patchwork")
# ==========================================================================
# Add information about mutations
# ==========================================================================
mutation_list <- readRDS("18-screening-for-erap2-mutations-in-postmortem/rdata_files/01-list-of-selected-mutations-for-each-subject.rds")

patient_specific_mutation_list <- mutation_list$patient_mutation_list
frequency_matrix <- mutation_list$frequency_matrix
selected_mutation_df <- mutation_list$mutation_list

#
# Select mutations that are present in at least 10% of ALS cases and 
# less that 5% of Control groups.
# --------------------------------------------------------------------------
candidate_mutations_df <- selected_mutation_df[selected_mutation_df$als > 0.1 & selected_mutation_df$ctrl < 0.05,]

candidate_mutations_df <- candidate_mutations_df[,c(1,2,2,3,4,5)]
colnames(candidate_mutations_df) <- c("chr", "start", "end", "ref", "alt", "snp_id")

# write.table(x = candidate_mutations_df, 
#     sep = "\t",
#     col.names = T,
#     row.names = F,
#     file = file.path(analysis_input_dir, "02-erap2-inferred-mutation-candidates-with-als-frequency-gt-10perc-ctrl-frequency-lt-5perc.avinput"))


all_selected_mutations <- selected_mutation_df[,c(1,2,2,3,4,5)]
colnames(all_selected_mutations) <- c("chr", "start", "end", "ref", "alt", "snp_id")
all_selected_mutations$mutation <- paste(all_selected_mutations$ref, all_selected_mutations$alt, sep = "/")
all_selected_mutations$strand <- "."


write.table(x = all_selected_mutations[,c("chr", "start", "end", "mutation", "strand", "snp_id")], 
    sep = "\t",
    col.names = F,
    row.names = F,
    quote = F, 
    file = file.path(analysis_input_dir, "03-all-erap2-inferred-mutation-candidates-ensembl-format.txt"))
