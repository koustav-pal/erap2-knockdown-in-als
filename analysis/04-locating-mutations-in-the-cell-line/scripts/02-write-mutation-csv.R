# ==========================================================================
# Intial script setup
# ==========================================================================
data_dir <- normalizePath("/camp/lab/luscomben/home/users/palk/")
project_dir <- file.path(data_dir, "Projects","012-Ben-Clarke-ERAP1-KD/analysis/")
setwd(project_dir)
# disable devtools errors from warnings mode
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# ==========================================================================
# Imports
# ==========================================================================
require("GenomicRanges")
require("rtracklayer")
require("stringr")
require("BSgenome.Hsapiens.UCSC.hg38")
require("GenomeInfoDb")

mutations_bed <- import.bed("/camp/home/palk/home/users/palk/Genomes/als_mutation_list.bed")
seqlevels(mutations_bed) <- mapSeqlevels(seqlevels(mutations_bed), "UCSC")
experiment_matrix_df <- data.table::fread("01-alignment/input_files/nfcore_sampleinfo.csv", data.table = F)
experiment_matrix_df$ERAP <- str_split(experiment_matrix_df$sample, pattern = "_", simplify = T)[,2]
experiment_matrix_df$replicate <- str_split(experiment_matrix_df$sample, pattern = "_", simplify = T)[,1]
experiment_matrix_df <- experiment_matrix_df[,-c(2,3,4)]
experiment_matrix_df$condition <- ifelse(experiment_matrix_df$replicate %in% c("C2", "C6"), "ctrl", "vcp")
experiment_matrix_df$fraction <- "whole"
# ==========================================================================
# Fetch genes
# ==========================================================================
hg38 <- BSgenome.Hsapiens.UCSC.hg38


alt_allele <- str_split(string = str_split(string = mutations_bed$name, 
    pattern = ":", simplify = T)[,5], pattern = ">", simplify = T)[,2]

fetch_left_seq_pos <- getSeq(hg38, names = as.vector(seqnames(mutations_bed)), 
    start = start(mutations_bed) - 10,
    end = start(mutations_bed) - 1)
fetch_left_seq_neg <- getSeq(hg38, names = as.vector(seqnames(mutations_bed)), 
    start = start(mutations_bed) - 10,
    end = start(mutations_bed) - 1)
fetch_right_seq <- getSeq(hg38, names = as.vector(seqnames(mutations_bed)), 
    start = start(mutations_bed) + 1,
    end = start(mutations_bed) + 10)

mutation_df <- data.frame(names = mutations_bed$name, 
    left_seq = fetch_left_seq, 
    alt_allele = alt_allele, 
    right_seq = fetch_right_seq, 
    chromosome = mapSeqlevels(as.vector(seqnames(mutations_bed)), "ENSEMBL"))

write.csv(x = mutation_df, row.names = F, quote = F, 
    file = "/camp/home/palk/home/users/palk/Genomes/als_mutation_list_mutscan.csv")