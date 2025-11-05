# ==========================================================================
# Intial script setup
# ==========================================================================
fetch_file_from_nas <- function(x, method){
    tempfile <- tempfile()
    system(paste(paste("rsync patanilab@10.246.164.99",x,sep = ":"), tempfile))
    return(method(tempfile))
}
data_dir <- normalizePath("/export2/home/koustav")
project_dir <- file.path(data_dir, "Projects",
    "0006-phase-condensation-of-repeats-in-neurons/")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "00-manuscript", "13-nygc-consortium")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plot_dir <- file.path(analysis_dir, "plots")
rslurm_dir <- file.path(analysis_dir, "rslurm")
alignment_dir <- file.path(project_dir, "analysis", "Alignments")
object_list <- fetch_file_from_nas("/var/services/homes/patanilab/2025-crick-user-folders/palk/Projects/0006-phase-condensation-of-repeats-in-neurons/analysis/00-manuscript/13-nygc-consortium/rdata_files/03-02-list_of_counts_stat_samplesheet_from_nygc.rds", method = readRDS)
isPairedEnd <- T
strand_specificity <- 2
sampleinfo_df <- read.csv(file.path(analysis_input_dir, "00_filtered_sample_gsm_list.txt"))
clinical_metadata_df <- fetch_file_from_nas("/var/services/homes/patanilab/2025-crick-projects/proj-patanir-seq/palk/0103-NYGC-consortium/input_files/samplesheet.csv", method = read.csv)
# ==========================================================================
# Imports
# ==========================================================================
require("GenomicRanges")
require("rslurm")
require("stringr")
require("matrixStats")
# require("sva")
require("DESeq2")
# library("bladderdata")
source("../09-01-QC-compare-fractional-counts-vs-abs-counts/scripts/src/deseq_functions.R")
source("../../useful-functions/plot-utils.R")

clinical_metadata_df$type <- clinical_metadata_df$mutation
clinical_metadata_df$mutation[!grepl(x = tolower(clinical_metadata_df$genetic_test), pattern = "positive for") 
    & !is.na(clinical_metadata_df$genetic_test)
    & clinical_metadata_df$family_history_of_als_ftd %in% c("Unknown", "No", "")
    & clinical_metadata_df$mutation != "ctrl"] <- "sporadic"

clinical_metadata_df[clinical_metadata_df$mutation == "tardbp",]

clinical_metadata_df$type[clinical_metadata_df$mutation == "c9orf72"] <- "c9orf72-ALS"
clinical_metadata_df$type[clinical_metadata_df$mutation != "c9orf72" & clinical_metadata_df$family_history_of_als_ftd == "Yes"] <- "other_mutation-fALS"
clinical_metadata_df$type[clinical_metadata_df$mutation != "c9orf72" & clinical_metadata_df$family_history_of_als_ftd != "Yes"] <- "sporadic"
clinical_metadata_df$type[clinical_metadata_df$mutation == "ctrl"] <- "control"

counts_matrix <- object_list$counts

sample_df <- object_list$sample_df
sample_df <- sample_df[sample_df$sample_name %in% sampleinfo_df$Experiment,]

sample_df$subject_id <- clinical_metadata_df$external_subject_id[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$sample_id <- clinical_metadata_df$external_sample_id[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$age_at_symptom_onset <- clinical_metadata_df$age_at_symptom_onset[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$age_at_death <- clinical_metadata_df$age_at_death[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$disease_duration_in_months <- clinical_metadata_df$disease_duration_in_months[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$mutation <- clinical_metadata_df$mutation[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$site_of_motor_onset <- clinical_metadata_df$site_of_motor_onset[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$family_history_of_als_ftd <- clinical_metadata_df$family_history_of_als_ftd[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]

sample_df$tissue <- clinical_metadata_df$sample_source[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$instrument <- sampleinfo_df$Instrument[match(sample_df$sample_name, sampleinfo_df$Experiment)]
sample_df$sex_genotype <- clinical_metadata_df$sex_genotype[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$rin <- clinical_metadata_df$rin[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$family_history_of_als_ftd <- clinical_metadata_df$family_history_of_als_ftd[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$mnd_with_ftd <- clinical_metadata_df$mnd_with_ftd[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$condition <- clinical_metadata_df$condition[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$condition[sample_df$mnd_with_ftd == "Yes"] <- "als+ftd"

sample_df$mnd_with_ftd <- clinical_metadata_df$mnd_with_ftd[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$mnd_with_ftd[sample_df$mnd_with_ftd == ""] <- "N/A"
sample_df$mnd_with_ftd[sample_df$mnd_with_ftd == "Not Applicable"] <- "N/A"
sample_df$mnd_with_dementia <- clinical_metadata_df$mnd_with_dementia[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]
sample_df$mnd_with_dementia[sample_df$mnd_with_dementia == ""] <- "N/A"
sample_df$mnd_with_dementia[sample_df$mnd_with_dementia == "Not Applicable"] <- "N/A"
sample_df$type <- clinical_metadata_df$type[match(sample_df$sample_name, clinical_metadata_df$experiment_accession)]

no_low_rin_sample_df <- sample_df[sample_df$rin > 4 & !is.na(sample_df$rin) & sample_df$condition %in% c("als", "ctrl", "als+ftd"),]

low_counts_removed_matrix <- counts_matrix[rowMeans(counts_matrix) > 1,]

current_sampleinfo <- no_low_rin_sample_df[no_low_rin_sample_df$column_id %in% colnames(low_counts_removed_matrix),]
current_sampleinfo$tissue[current_sampleinfo$tissue %in% c("Cortex_Frontal", "Cortex_Temporal", "Cortex_Occipital")] <- "Cortex"
current_sampleinfo$tissue[current_sampleinfo$tissue %in% c("Cortex_Motor_Lateral", "Cortex_Motor_Medial", "Cortex_Motor_Unspecified")] <- "Cortex_Motor"
current_sampleinfo <- current_sampleinfo[current_sampleinfo$tissue != "Hippocampus",]


unique(current_sampleinfo[,c("subject_id", "mutation")])$mutation %>% table



saveRDS(object = current_sampleinfo, file = file.path("/nemo/lab/patanir/home/users/palk/Projects/012-Ben-Clarke-ERAP1-KD/analysis/11-postmortem-data-analysis/rdata_files/00-nygc-patient-sampleinfo-datasheet.rds"))