data_dir <- normalizePath("/export2/home/koustav/")
project_dir <- file.path(data_dir, "Projects", "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "18-screening-for-erap2-mutations-in-postmortem/")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
temp_dir <- file.path(analysis_dir, "tmp")
tmp_dir <- file.path(analysis_dir, "slurm")
genic_count_list <- readRDS(file.path("../11-postmortem-data-analysis/rdata_files/01-gene-level-counts-of-postmortem-samples.rds"))
# ==========================================================================
# Analysis
# ==========================================================================
coldata_df <- genic_count_list$coldata
write.table(x = coldata_df, 
	file = file.path(analysis_input_dir, "02-list-of-patient-samples.txt"),
	sep = "\t", quote = F, col.names = T, row.names = F)