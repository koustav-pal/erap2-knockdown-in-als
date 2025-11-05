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
# rmats_data_list <- readRDS(file.path(analysis_save_dir, "01-imported-rmats-object-list.rds"))
tmp_dir <- file.path(analysis_dir, "slurm")
genic_count_list <- readRDS(file.path(analysis_save_dir, "01-gene-level-counts-of-postmortem-samples.rds"))
# candidate_de_genes_affected_by_erap2_kd <- readRDS("../02-01-differential-expression-profiling/rdata_files/02-significantly-affected-genes-by-kd-in-als-vs-control.rds")
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
require("CEMiTool")
source("scripts/src/functions.R")
# source("../src/rmats-analysis/src/source.R")
# ==========================================================================
# Analysis
# ==========================================================================
annotation_ranges <- import.gff2(annotation_gtf)
annotation_gene_ranges <- annotation_ranges[annotation_ranges$type == "gene"]

# feature_ranges_with_repeats <- significant_de_event_list$repeat_containing_ranges
counts_matrix <- genic_count_list$counts 
tpm_matrix <- counts_matrix / (width(annotation_gene_ranges[match(rownames(counts_matrix), annotation_gene_ranges$gene_id)]) / 1e3)
tpm_matrix <- t(t(tpm_matrix) / (colSums(tpm_matrix) / 1e6))


coldata_df <- genic_count_list$coldata
coldata_df$expr <- tpm_matrix[annotation_gene_ranges$gene_id[annotation_gene_ranges$gene_name == "ERAP2" & !is.na(annotation_gene_ranges$gene_name)],
    coldata_df$column_id]
coldata_df$tissue_group <- "spinal_cord"
coldata_df$tissue_group[coldata_df$tissue %in% c("Cortex_Motor", "Cortex", "Cerebellum")] <- "brain"

subset_matrix <- as.data.frame(counts_matrix[rownames(counts_matrix) %in% 
    annotation_gene_ranges$gene_id[annotation_gene_ranges$gene_biotype == "protein_coding"],])
keep_matrix <- subset_matrix[(rowSums(subset_matrix == 0) / ncol(subset_matrix)) < 0.10,]


coldata_df$matrix_row <- match(coldata_df$column_id, colnames(keep_matrix))

cols_split_by_patient <- split(coldata_df, paste(coldata_df$subject_id, coldata_df$tissue_group))

patient_matrix <- do.call(cbind, lapply(cols_split_by_patient, function(x){
    if(nrow(x) == 1){
        return(keep_matrix[,x$matrix_row])
    }
    return(rowMeans(keep_matrix[,x$matrix_row]))
}))

patient_df <- do.call(rbind, lapply(cols_split_by_patient, function(temp_df){
    return_df <- unique(temp_df[,!(colnames(temp_df) %in% 
        c("X", "sample_name", "alignment_file", "alignment_file_dir",
          "isPairedEnd", "strand_specificity", "samplesheet_path", 
          "column_id", "sample_id", "mutation", "site_of_motor_onset",
          "instrument", "rin", "expr", "tissue", "matrix_row"))])
    if(nrow(return_df) > 1){
        # message(paste(unique(return_df$subject_id), unique(return_df$tissue)))
        if("N/A" %in% return_df$mnd_with_ftd & "Yes" %in% return_df$mnd_with_ftd){
            return_df <- return_df[return_df$mnd_with_ftd == "Yes",]
        }else if("N/A" %in% return_df$mnd_with_ftd & "No" %in% return_df$mnd_with_ftd){
            return_df <- return_df[return_df$mnd_with_ftd == "No",]
        }else{
            stop(paste(return_df$subject_id, return_df$tissue))
        }
        # message(nrow(return_df))
    }
    return_df$expr <- mean(temp_df$expr)
    return_df$column_id <- paste(return_df$subject_id, return_df$tissue_group)
    return(return_df)
}))
patient_df$age_at_death[patient_df$age_at_death == "90 or Older"] <- "90"
patient_df$age_at_death <- as.numeric(patient_df$age_at_death)
patient_df$age_at_onset <- as.numeric(patient_df$age_at_onset)
patient_df$disease_duration <- as.numeric(patient_df$disease_duration_in_months)

patient_df$disease_duration_in_months <- as.numeric(patient_df$disease_duration_in_months)

# select patients

selected_patients <- patient_df[patient_df$tissue_group == "spinal_cord" &
    patient_df$disease_duration_in_months <= 12 & 
    !is.na(patient_df$disease_duration_in_months),]






tissue_level_connectivity_list <- lapply(unique(patient_df$tissue_group), function(a_tissue){
    current_patients_df <- patient_df[patient_df$tissue_group == a_tissue &
        patient_df$disease_duration_in_months <= 12 & 
        !is.na(patient_df$disease_duration_in_months),]
        message(nrow(current_patients_df))
    # als_patient_expr_quantiles <- quantile(current_patients_df$expr[current_patients_df$condition != "ctrl"], c(0.9))
    anno_df <- current_patients_df
    # [
    #     current_patients_df$disease_duration_in_months <= 12 & 
    #     !is.na(current_patients_df$disease_duration_in_months) &
    #     current_patients_df$expr > als_patient_expr_quantiles,]
    keep_matrix_df <- as.data.frame(patient_matrix[,anno_df$column_id])
    cem <- cemitool(keep_matrix_df, annot = anno_df, 
        apply_vst = T, 
        sample_name_column = "column_id", class_column = "condition")
    return(cem)
})
names(tissue_level_connectivity_list) <- unique(patient_df$tissue_group)

# Fetch modules
module_df <- do.call(rbind, lapply(names(tissue_level_connectivity_list), function(x){
    message(x)
    temp_object <- tissue_level_connectivity_list[[x]]
    module_list <- module_genes(temp_object)
    if(!("ENSG00000164308" %in% module_list[,1])){
        return(NULL)
    }
    if(module_list[module_list[,1] == "ENSG00000164308",2] == "Not.Correlated"){
        return(NULL)
    }
    return_df <- module_list
    return(return_df)
}))

genes_in_same_module <- unique(module_df$genes)

saveRDS(object = module_df, file = file.path(analysis_save_dir, "04-postmortem-WGCNA-gene-modules.rds"))
