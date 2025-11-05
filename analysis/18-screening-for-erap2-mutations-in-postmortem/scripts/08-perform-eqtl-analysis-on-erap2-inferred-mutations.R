data_dir <- normalizePath("/export2/home/koustav/")
project_dir <- file.path(data_dir, "Projects", "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "18-screening-for-erap2-mutations-in-postmortem/")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_output_dir <- file.path(analysis_dir, "output_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
temp_dir <- file.path(analysis_dir, "tmp")
tmp_dir <- file.path(analysis_dir, "slurm")
genic_count_list <- readRDS(file.path("../11-postmortem-data-analysis/rdata_files/01-gene-level-counts-of-postmortem-samples.rds"))
annotation_gtf <- "/export2/home/koustav/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
# ==========================================================================
# Packages
# ==========================================================================
require("BSgenome.Hsapiens.UCSC.hg38")
require("stringr")
require("edgeR")
require("MatrixEQTL")
source("scripts/src/source.R")
# ==========================================================================
# Analysis: Prepare the list of SNPs in ERAP1/2
# ==========================================================================
mutation_df <- data.table::fread(file.path(analysis_input_dir, "snp_result.txt"), data.table = F)
erap2_location_df <- data.frame(chr = 5, start = 96875984, end = 96919716)
mutation_df <- mutation_df[mutation_df$pos >= erap2_location_df$start & mutation_df$pos <= erap2_location_df$end,]

hg38 <- BSgenome.Hsapiens.UCSC.hg38
all_mutation_df <- do.call(rbind, lapply(seq_len(nrow(mutation_df)), function(x){
    message(x)
    all_alleles <- str_split(string = mutation_df$variation[x], 
        pattern = ">", simplify = T)
    ref_allele <- all_alleles[1,1]
    alt_alleles <- all_alleles[1,2]
    if(grepl(x = all_alleles[1,2], pattern = ",")){
        alt_alleles <- str_split(string = all_alleles[1,2], 
        pattern = ",", simplify = T)[1,]
    }
    return_mutation_df <- data.frame(chromosome = mapSeqlevels(paste("chr",mutation_df[x,"#chr"], sep = ""), "ENSEMBL"), 
        pos = mutation_df[x,"pos"],
        reference = ref_allele, 
        alternative = alt_alleles,
        name = mutation_df[x,"snp_id"])
    return(return_mutation_df)
}))

# ==========================================================================
# Analysis: Prepare the coldata
# ==========================================================================
coldata_df <- genic_count_list$coldata
coldata_df$sample_row <- seq_len(nrow(coldata_df))
coldata_df$sample_name <- paste(coldata_df$condition, coldata_df$sample_row, sep = "_")

coldata_df$tissue_group <- "spinal_cord"
coldata_df$tissue_group[coldata_df$tissue %in% c("Cortex", "Cortex_Motor", "Cerebellum")] <- "brain"


# ==========================================================================
# Analysis: Import the screening results and prepare per patient mutation profile
# ==========================================================================

mutation_screening_df <- data.table::fread(file = "tmp/ALL_MUTATION_SCREENING_RESULTS.txt", header = T, fill = T, data.table = F)

columns <- c(reference_fwd = "reff", reference_rev = "refv", alternate_fwd = "altf", alternate_rev = "altv")
mutation_screening_df$snp_id <- all_mutation_df$name[match(
    paste(mutation_screening_df[,1], mutation_screening_df[,2], mutation_screening_df[,3], mutation_screening_df[,4], sep = "_"),
    paste(all_mutation_df[,1], all_mutation_df[,2], all_mutation_df[,3], all_mutation_df[,4], sep = "_"))]

per_patient_mutation_profile <- do.call(cbind, lapply(split(coldata_df, coldata_df$subject_id), function(a_row){
    message(unique(a_row$subject_id))
    ids <- paste(a_row$condition, a_row$sample_row, sep = "_")
    column_rep <- rep(columns, times = length(ids))
    ids_rep <- rep(ids, each = length(columns))
    all_columns <- paste(ids_rep, column_rep, sep = ".")
    column_split <- split(all_columns, column_rep)
    a_snp_freq_matrix <- do.call(cbind, lapply(column_split, function(x){
        if(length(x) == 1){
            return(mutation_screening_df[,x])
        }
        return(rowSums(mutation_screening_df[,x]))
    }))
    a_snp_freq_matrix <- cbind(a_snp_freq_matrix, "total" = rowSums(a_snp_freq_matrix))
    a_snp_freq_matrix <- cbind(a_snp_freq_matrix, "vaf" = rowSums(a_snp_freq_matrix[,c("altf", "altv")]) / a_snp_freq_matrix[,"total"])

    strand_bias <- vapply(seq_len(nrow(a_snp_freq_matrix)), function(x){
        fisher.test(matrix(data = c(a_snp_freq_matrix[x,3],a_snp_freq_matrix[x,1],a_snp_freq_matrix[x,4],a_snp_freq_matrix[x,2]),nrow = 2, ncol = 2))$p.value
    },1)
    a_snp_freq_matrix <- cbind(a_snp_freq_matrix, strand_bias)
    colnames(a_snp_freq_matrix) <- paste(unique(a_row$subject_id), colnames(a_snp_freq_matrix), sep = "_")
    return(a_snp_freq_matrix)
}))


# ==========================================================================
# Analysis: Filter the patient mutation profile keeping only 
# those mutations with a depth > 10
# ==========================================================================

summarised_mutation_df <- cbind(mutation_screening_df[,c(1,2,3,4,ncol(mutation_screening_df))], per_patient_mutation_profile)
unique_coldata_df <- unique(coldata_df[, c("subject_id", "disease_duration_in_months", "age_at_symptom_onset", "age_at_death", "sex_genotype", "condition")])

unique_coldata_df$age_at_death[unique_coldata_df$age_at_death == "90 or Older"] <- "90"
unique_coldata_df$age_at_death <- as.numeric(unique_coldata_df$age_at_death)
unique_coldata_df$age_at_onset <- as.numeric(unique_coldata_df$age_at_symptom_onset)
unique_coldata_df$disease_duration <- as.numeric(unique_coldata_df$disease_duration_in_months)

all_columns <- c(columns, c("total", "vaf", "strand_bias"))
min_reads <- 10

patient_specific_mutation_list <- lapply(seq_len(nrow(unique_coldata_df)), function(x){
    current_row <- unique_coldata_df[x,]
    col_names <- paste(current_row$subject_id, all_columns, sep = "_")
    patient_specific_mutations <- summarised_mutation_df[,c(1:5, which(colnames(summarised_mutation_df) %in% col_names))]
    colnames(patient_specific_mutations) <- c(colnames(patient_specific_mutations)[1:5], all_columns)
    filtered_patient_specific_mutations <- filter_mutations_for_eqtl(patient_specific_mutations, depth = 10)
    a_vector <- rep(0, nrow(summarised_mutation_df))
    names(a_vector) <- summarised_mutation_df$snp_id
    genotype_vector <- a_vector
    vaf_vector <- a_vector
    vaf_vector[filtered_patient_specific_mutations$snp_id] <- filtered_patient_specific_mutations$vaf
    a_vector[filtered_patient_specific_mutations$snp_id] <- 1
    genotype_vector[filtered_patient_specific_mutations$snp_id[filtered_patient_specific_mutations$vaf > 0.1]] <- 1
    list(filtered_mutations = filtered_patient_specific_mutations, 
        occurence_vector = a_vector,
        genotype_vector = genotype_vector,
        vaf_vector = vaf_vector)
})
names(patient_specific_mutation_list) <- unique_coldata_df$subject_id

# ==========================================================================
# Analysis: Based on occurence rate of 50% in either ALS or Control we 
# filter for the SNPs that we want to test. Patients have a genotype label
# of 1 if they have VAF > 0.2 and occurence = 1 (depth > 10 reads). However, 
# if patients have an occurrence = 1 but VAF < 0.2 they are set to NA.
# ==========================================================================

occurence_matrix <- do.call(cbind, lapply(patient_specific_mutation_list, function(x){
    x$occurence_vector
}))
colnames(occurence_matrix) <- unique_coldata_df$subject_id

als_occurrence_filter <- rowMeans(occurence_matrix[,unique_coldata_df$subject_id[unique_coldata_df$condition == "als"]]) > 0.1
ctrl_occurrence_filter <- rowMeans(occurence_matrix[,unique_coldata_df$subject_id[unique_coldata_df$condition == "ctrl"]]) < 0.05
across_all_occurrence_filter <- rowSums(occurence_matrix[,unique_coldata_df$subject_id[unique_coldata_df$condition %in% c("ctrl", "als")]]) >= 3

genotype_matrix <- do.call(cbind, lapply(patient_specific_mutation_list, function(x){
    x$genotype_vector
}))
colnames(genotype_matrix) <- unique_coldata_df$subject_id


vaf_matrix <- do.call(cbind, lapply(patient_specific_mutation_list, function(x){
    x$vaf_vector
}))
colnames(vaf_matrix) <- unique_coldata_df$subject_id

genotype_matrix[occurence_matrix == 1 & vaf_matrix <= 0.1] <- NA
genotype_matrix <- genotype_matrix[als_occurrence_filter & ctrl_occurrence_filter & across_all_occurrence_filter,]

# ==========================================================================
# Analysis: Summarise all counts per gene for cervical, 
# lumbar and thoracic spinal cord segments
# ==========================================================================

gene_counts_matrix <- genic_count_list$counts

cervical_patient_df <- coldata_df[coldata_df$tissue %in% "Spinal_Cord_Cervical",]
cervical_counts_matrix <- do.call(cbind, lapply(split(cervical_patient_df$column_id, cervical_patient_df$subject_id), function(x){
    if(length(x) == 1){
        return(gene_counts_matrix[,x])
    }
    return(rowMeans(gene_counts_matrix[,x]))
}))
cervical_subject_ids <- intersect(intersect(colnames(cervical_counts_matrix), colnames(genotype_matrix)), unique_coldata_df$subject_id[unique_coldata_df$condition %in% c("als")])

cervical_list <- list(
    coldata = unique_coldata_df[match(cervical_subject_ids,unique_coldata_df$subject_id),],
    counts = cervical_counts_matrix[,cervical_subject_ids],
    genotype_matrix = genotype_matrix[,cervical_subject_ids])


lumbar_patient_df <- coldata_df[coldata_df$tissue %in% "Spinal_Cord_Lumbar",]
lumbar_counts_matrix <- do.call(cbind, lapply(split(lumbar_patient_df$column_id, lumbar_patient_df$subject_id), function(x){
    if(length(x) == 1){
        return(gene_counts_matrix[,x])
    }
    return(rowMeans(gene_counts_matrix[,x]))
}))
lumbar_subject_ids <- intersect(intersect(colnames(lumbar_counts_matrix), colnames(genotype_matrix)), unique_coldata_df$subject_id[unique_coldata_df$condition %in% c("als")])

lumbar_list <- list(
    coldata = unique_coldata_df[match(lumbar_subject_ids,unique_coldata_df$subject_id),],
    counts = lumbar_counts_matrix[,lumbar_subject_ids],
    genotype_matrix = genotype_matrix[,lumbar_subject_ids])


thoracic_patient_df <- coldata_df[coldata_df$tissue %in% "Spinal_Cord_Thoracic",]
thoracic_counts_matrix <- do.call(cbind, lapply(split(thoracic_patient_df$column_id, thoracic_patient_df$subject_id), function(x){
    if(length(x) == 1){
        return(gene_counts_matrix[,x])
    }
    return(rowMeans(gene_counts_matrix[,x]))
}))
thoracic_subject_ids <- intersect(intersect(colnames(thoracic_counts_matrix), colnames(genotype_matrix)), unique_coldata_df$subject_id[unique_coldata_df$condition %in% c("als")])

thoracic_list <- list(
    coldata = unique_coldata_df[match(thoracic_subject_ids,unique_coldata_df$subject_id),],
    counts = thoracic_counts_matrix[,thoracic_subject_ids],
    genotype_matrix = genotype_matrix[,thoracic_subject_ids])

all_tissue_list <- list(cervical = cervical_list, 
    lumbar = lumbar_list,
    thoracic = thoracic_list)

# ==========================================================================
# Analysis: 
# ==========================================================================

annotation_ranges <- import.gff2(annotation_gtf)
gene_ranges <- annotation_ranges[annotation_ranges$type == "gene"]


all_eqtl_results <- lapply(all_tissue_list, function(a_list){
    # Genotype data
    snps <- SlicedData$new()
    snps$CreateFromMatrix(as.matrix(a_list$genotype_matrix))
    snps$ResliceCombined(sliceSize = 100)

    snp_pos_df <- all_mutation_df[match(rownames(a_list$genotype_matrix), all_mutation_df$name),]
    snp_pos_df <- snp_pos_df[,c("name", "chromosome", "pos")]
    snp_pos_df$pos <- as.numeric(snp_pos_df$pos)
    colnames(snp_pos_df) <- c("snp_id", "chr", "pos")


    dge <- DGEList(counts = a_list$counts)
    dge <- calcNormFactors(dge, method = "TMM")
    log2_expr <- cpm(dge, log = TRUE, prior.count = 1)

    keep_gene_ids <- intersect(rownames(log2_expr), gene_ranges$gene_id)
    log2_expr <- log2_expr[keep_gene_ids,]
    # Expression data
    gene <- SlicedData$new()
    gene$CreateFromMatrix(as.matrix(log2_expr))
    gene$ResliceCombined(sliceSize = 100)

    covariates <- a_list$coldata
    covariate_df <- covariates[
        match(colnames(a_list$genotype_matrix), covariates$subject_id),
        c("subject_id", "age_at_death", "sex_genotype", "disease_duration")]
    covariate_df$sex_genotype <- ifelse(covariate_df$sex_genotype == "XY", 1, 0)
    # covariate_df$condition <- ifelse(covariate_df$condition == "als", 1, 0)
    covariate_matrix <- covariate_df[,-1]
    rownames(covariate_matrix) <- covariate_df$subject_id
    covariate_matrix <- t(covariate_matrix)

    covariates_object <- SlicedData$new()
    covariates_object$CreateFromMatrix(covariate_matrix)

    gene_pos_ranges <- annotation_ranges[match(keep_gene_ids, gene_ranges$gene_id)]
    gene_pos_df <- data.frame(
      geneid = keep_gene_ids,
      chr = as.vector(seqnames(gene_pos_ranges)),
      start = start(gene_pos_ranges),
      end = end(gene_pos_ranges))

    results <- Matrix_eQTL_main(
      snps = snps,
      gene = gene,
      cvrt = covariates_object,  # covariates; can be empty for now
      output_file_name = NULL,
      pvOutputThreshold = 1e-5,
      useModel = modelLINEAR,
      errorCovariance = numeric(),  # identity covariance
      verbose = TRUE,
      output_file_name.cis = NULL,
      pvOutputThreshold.cis = 5e-2,
      snpspos = snp_pos_df,
      genepos = gene_pos_df,
      cisDist = 1e6)
})
