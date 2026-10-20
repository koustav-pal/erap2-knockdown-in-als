data_dir <- normalizePath("/home/koustav/")
project_dir <- file.path(data_dir, "Projects", "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "23-astrocyte-metaanalysis")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
genic_count_list <- readRDS(file.path(analysis_save_dir, "gene_level_counts.rds"))
annotation_gtf <- "/data/shared/public_datasets/Genomes/ENSEMBL_hg38/Homo_sapiens.GRCh38.104.chr_patch_hapl_scaff.gtf"
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
require("rtracklayer")
require("STRINGdb")
require("umap")
require("CEMiTool")
# source("scripts/src/functions.R")
# source("../src/rmats-analysis/src/source.R")
# ==========================================================================
# Analysis
# ==========================================================================
annotation_ranges <- import.gff2(annotation_gtf)
annotation_gene_ranges <- annotation_ranges[annotation_ranges$type == "gene"]

# feature_ranges_with_repeats <- significant_de_event_list$repeat_containing_ranges

counts_matrix <- genic_count_list$counts 
coldata_df <- genic_count_list$sampleinfo
coldata_df <- coldata_df[coldata_df$study != "ben_astrocyte_sfpqkd",]
coldata_df <- coldata_df[!(coldata_df$treatment %in% c("a0", "a1")),]
coldata_df$condition[coldata_df$condition != "control"] <- "als"



subset_matrix <- as.data.frame(counts_matrix[rownames(counts_matrix) %in% 
    annotation_gene_ranges$gene_id[annotation_gene_ranges$gene_biotype == "protein_coding"],])
keep_matrix <- subset_matrix[(rowSums(subset_matrix == 0) / ncol(subset_matrix)) < 0.10,]
keep_matrix <- keep_matrix[,coldata_df$sample]
keep_matrix_df <- as.data.frame(keep_matrix)



cem <- cemitool(keep_matrix_df, annot = coldata_df, 
    apply_vst = T, 
    sample_name_column = "sample", class_column = "condition")

tissue_level_connectivity_list <- list("cell_line" = cem)
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


module_df$gene_name <- annotation_gene_ranges$gene_name[match(module_df$genes, annotation_gene_ranges$gene_id)]

saveRDS(object = module_df, file = file.path(analysis_save_dir, "01-iPSC-astrocyte-WGCNA-gene-modules.rds"))
postmortem_module_df <- readRDS(file.path("../11-postmortem-data-analysis/rdata_files/04-postmortem-WGCNA-gene-modules.rds"))


selected_postmortem_module <- postmortem_module_df[postmortem_module_df$modules == postmortem_module_df$modules[postmortem_module_df$genes == "ENSG00000164308"],]
selected_astrocyte_module <- module_df[module_df$modules == module_df$modules[module_df$genes == "ENSG00000164308"],]

# length(intersect(selected_astrocyte_module$genes , selected_postmortem_module$genes)) / length(union(selected_astrocyte_module$genes , selected_postmortem_module$genes))
# union(selected_astrocyte_module$genes , selected_postmortem_module$genes)
#   [1] "ENSG00000108821" "ENSG00000115461" "ENSG00000149591" "ENSG00000137801"
#   [5] "ENSG00000060718" "ENSG00000107796" "ENSG00000130635" "ENSG00000120708"
#   [9] "ENSG00000146674" "ENSG00000166147" "ENSG00000091986" "ENSG00000198467"
#  [13] "ENSG00000144810" "ENSG00000106366" "ENSG00000163453" "ENSG00000147604"
#  [17] "ENSG00000148848" "ENSG00000134013" "ENSG00000177954" "ENSG00000105974"
#  [21] "ENSG00000166923" "ENSG00000075223" "ENSG00000159176" "ENSG00000101335"
#  [25] "ENSG00000104332" "ENSG00000137809" "ENSG00000092758" "ENSG00000161638"
#  [29] "ENSG00000122026" "ENSG00000198910" "ENSG00000186340" "ENSG00000158258"
#  [33] "ENSG00000103888" "ENSG00000139926" "ENSG00000113739" "ENSG00000173546"
#  [37] "ENSG00000130176" "ENSG00000164106" "ENSG00000197956" "ENSG00000138061"
#  [41] "ENSG00000148677" "ENSG00000166741" "ENSG00000128510" "ENSG00000133816"
#  [45] "ENSG00000134954" "ENSG00000116016" "ENSG00000135636" "ENSG00000115232"
#  [49] "ENSG00000113389" "ENSG00000130203" "ENSG00000138080" "ENSG00000168386"
#  [53] "ENSG00000095752" "ENSG00000136542" "ENSG00000198682" "ENSG00000141526"
#  [57] "ENSG00000117228" "ENSG00000185070" "ENSG00000157613" "ENSG00000096433"
#  [61] "ENSG00000135269" "ENSG00000149596" "ENSG00000019549" "ENSG00000213949"
#  [65] "ENSG00000118515" "ENSG00000151572" "ENSG00000116774" "ENSG00000105976"
#  [69] "ENSG00000172348" "ENSG00000102760" "ENSG00000132031" "ENSG00000167772"
#  [73] "ENSG00000074416" "ENSG00000164308" "ENSG00000067066" "ENSG00000087510"
#  [77] "ENSG00000104324" "ENSG00000134668" "ENSG00000164687" "ENSG00000169129"
#  [81] "ENSG00000198597" "ENSG00000054356" "ENSG00000005243" "ENSG00000050344"
#  [85] "ENSG00000148926" "ENSG00000081923" "ENSG00000118785" "ENSG00000026508"
#  [89] "ENSG00000112096" "ENSG00000151929" "ENSG00000090104" "ENSG00000173369"
#  [93] "ENSG00000177575" "ENSG00000128016" "ENSG00000133048" "ENSG00000105835"
#  [97] "ENSG00000168209" "ENSG00000142089" "ENSG00000275395" "ENSG00000221869"
# [101] "ENSG00000125148" "ENSG00000117525" "ENSG00000047457" "ENSG00000096060"
# [105] "ENSG00000124762" "ENSG00000112715" "ENSG00000165507" "ENSG00000075426"
# [109] "ENSG00000164764" "ENSG00000072110" "ENSG00000122862" "ENSG00000099860"
# [113] "ENSG00000121966" "ENSG00000143226" "ENSG00000064886" "ENSG00000060138"
# [117] "ENSG00000018280" "ENSG00000059804" "ENSG00000154262" "ENSG00000135046"
# [121] "ENSG00000134531" "ENSG00000028137" "ENSG00000184557" "ENSG00000108551"
# [125] "ENSG00000092969" "ENSG00000203747" "ENSG00000038945" "ENSG00000135245"
# [129] "ENSG00000064932" "ENSG00000185201" "ENSG00000196136" "ENSG00000145623"
# [133] "ENSG00000125810" "ENSG00000197249" "ENSG00000104635" "ENSG00000124491"
# [137] "ENSG00000108846" "ENSG00000162645" "ENSG00000187193" "ENSG00000116741"
# [141] "ENSG00000197405" "ENSG00000137462" "ENSG00000155659" "ENSG00000172578"
# [145] "ENSG00000108691" "ENSG00000025708" "ENSG00000171051" "ENSG00000147872"
# [149] "ENSG00000065911" "ENSG00000116717" "ENSG00000169715" "ENSG00000165030"
# [153] "ENSG00000086544" "ENSG00000121005" "ENSG00000169508" "ENSG00000115594"
# [157] "ENSG00000185022" "ENSG00000106066" "ENSG00000196843" "ENSG00000162772"
# [161] "ENSG00000188483" "ENSG00000187479" "ENSG00000173846" "ENSG00000132965"
# [165] "ENSG00000105639" "ENSG00000143333" "ENSG00000158869" "ENSG00000163220"
# [169] "ENSG00000000938" "ENSG00000117226" "ENSG00000235750" "ENSG00000136997"
# [173] "ENSG00000173221" "ENSG00000178726" "ENSG00000179593" "ENSG00000127507"
# [177] "ENSG00000105697" "ENSG00000152049" "ENSG00000165929" "ENSG00000205364"
# [181] "ENSG00000106236" "ENSG00000138135" "ENSG00000169413" "ENSG00000139832"
# [185] "ENSG00000157766" "ENSG00000069399" "ENSG00000110079" "ENSG00000164674"
# [189] "ENSG00000198814" "ENSG00000222009" "ENSG00000143382" "ENSG00000136286"
# [193] "ENSG00000163814" "ENSG00000204577" "ENSG00000196639" "ENSG00000011422"
# [197] "ENSG00000258227" "ENSG00000113749" "ENSG00000166106" "ENSG00000188766"
# [201] "ENSG00000158055" "ENSG00000158714" "ENSG00000164761" "ENSG00000121933"
# [205] "ENSG00000115602" "ENSG00000151651" "ENSG00000137558" "ENSG00000273259"
# [209] "ENSG00000254415" "ENSG00000115919" "ENSG00000143546" "ENSG00000203710"
# [213] "ENSG00000100336" "ENSG00000188488" "ENSG00000020633" "ENSG00000110848"
# [217] "ENSG00000101916" "ENSG00000175793" "ENSG00000159167"
# intersect(selected_astrocyte_module$genes , selected_postmortem_module$genes)
#  [1] "ENSG00000137801" "ENSG00000120708" "ENSG00000106366" "ENSG00000134013"
#  [5] "ENSG00000113739" "ENSG00000166741" "ENSG00000141526" "ENSG00000149596"
#  [9] "ENSG00000167772" "ENSG00000164308" "ENSG00000134668" "ENSG00000164687"
# [13] "ENSG00000148926"