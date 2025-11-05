data_dir <- normalizePath("/export2/home/koustav/")
project_dir <- file.path(data_dir, "Projects", "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "14-01-erap2-ov-proteomics")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
# erap2_kd_df <- read.table(file.path(analysis_input_dir, "report.unique_genes_matrix.tsv"), header = T, sep = "\t")
# annotation_df <- read.csv(file.path(analysis_input_dir, "annotation.csv"))
annotation_gtf <- "/export2/home/koustav/Projects/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
# ==========================================================================
# Imports
# ==========================================================================
require("stringr")
library("biomaRt")
require("org.Hs.eg.db")
# source("useful_functions/plot_utils.R")
# ==========================================================================
# Analysis
# ==========================================================================

differential_protein_enrichment_list <- readRDS(file.path(analysis_save_dir, "02-01-erap2-OV-rots-differential-enrichment-analysis-results.rds"))
erap2_ov_df <- read.csv(file.path(analysis_input_dir, "ERAP2_overexpression_proteomics.csv"))
unimputed_effect_df <- differential_protein_enrichment_list$unimputed_data

gene_name_indexes <- vapply(unimputed_effect_df$protein_name, function(x){
    if(x %in% erap2_ov_df[,"T..Protein.Names"]){
        return(which(erap2_ov_df[,"T..Protein.Names"] == x))
    }else{
        return(which(grepl(x = erap2_ov_df[,"T..Protein.Names"], pattern = x) & grepl(x = erap2_ov_df[,"T..Protein.Names"], pattern = ";")))
    }
},1)
unimputed_effect_df$gene_names <- erap2_ov_df[gene_name_indexes,"T..Genes"]

all_gene_names <- do.call(c, str_split(unimputed_effect_df$gene_names, pattern = ";"))
unique_gene_names <- unique(all_gene_names)

mart <- useEnsembl("ensembl","hsapiens_gene_ensembl")
all_info_list <- lapply(split(unique_gene_names, ceiling(seq_along(unique_gene_names)/1000)), function(x){
    message("Done")
    all_attributes <- listAttributes(mart)
    # all_attributes[grepl(x = all_attributes$name, pattern = "synonym"),]
    get_ids_df <- getBM(c("entrezgene_id","ensembl_gene_id", "ensembl_gene_id_version", "uniprotswissprot","description","hgnc_symbol", "peptide"), "hgnc_symbol", 
        x, mart)

    all_ids <- get_ids_df$hgnc_symbol

    # which genes were not found
    unidentified_genes <- x[!(x %in% get_ids_df$hgnc_symbol)]
    if(length(unidentified_genes) > 0){
        get_unknown_ids_df <- getBM(c("entrezgene_id", "ensembl_gene_id", "ensembl_gene_id_version", "uniprotswissprot","description","hgnc_symbol", "external_synonym", "peptide"), 
            "external_synonym", unidentified_genes, mart)
        all_ids <- c(all_ids,get_unknown_ids_df$hgnc_symbol)
    }

    synonyms_df <- getBM(c("hgnc_symbol", "external_synonym"), "hgnc_symbol", 
        unique(all_ids), mart)
    intersect_cols <- intersect(colnames(get_unknown_ids_df), colnames(get_ids_df))
    all_biomart_info_df <- rbind(get_ids_df[,intersect_cols],get_unknown_ids_df[,intersect_cols])
    return(list(biomart_info = all_biomart_info_df, synonyms = synonyms_df))
})

all_biomart_info_df <- do.call(rbind, lapply(all_info_list, function(x){
    x$biomart_info
}))

synonyms_df <- do.call(rbind, lapply(all_info_list, function(x){
    x$synonyms
}))
message(all(unique_gene_names %in% all_biomart_info_df$hgnc_symbol | unique_gene_names %in% all_biomart_info_df$external_synonym))
saveRDS(object = list(biomart_info = all_biomart_info_df, synonyms = synonyms_df), file = file.path(analysis_save_dir, "00-01-erap2-ov-biomart-metadata.rds"))