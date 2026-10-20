# Adapted analysis; see provenance/source_ranges.json.
source(here::here("workflow", "scripts", "common.R"))
library(dplyr)
library(tidyr)
library(stringr)
library(Seurat)
library(patchwork)
library(viridis)
library(lme4)
library(emmeans)
library(rtracklayer)
library(broom)
library(reshape2)
annotated_cells <- readRDS(input("annotated_cells.rds"))
DefaultAssay(annotated_cells) <- "RNA"
gene_ranges <- rtracklayer::import(input("genes.gtf"))
gene_ranges <- gene_ranges[gene_ranges$type == "gene"]
cell_type_colours <- c(Astrocytes = "#E69F00", Endothelial = "#56B4E9", Lymphocytes = "#009E73", Microglia = "#F0E442",
    Neurons = "#0072B2", Oligodendrocytes = "#D55E00", OPC = "#CC79A7", Meninges = "#999999")
annotated_cells$predicted.cell_type <-
    as.character(annotated_cells$predicted.cell_type)

annotated_cells$predicted.cell_type[
    annotated_cells$predicted.cell_type == "Astrcytes"
] <- "Astrocytes"

annotated_cells$predicted.cell_type <-
    factor(annotated_cells$predicted.cell_type)

Idents(annotated_cells) <- "seurat_clusters"

annotated_cells$patient_id <- sub(
    "^g?([0-9]+).*",
    "\\1",
    as.character(annotated_cells$sample)
)

annotated_cells$condition <- ifelse(
    grepl(
        "(^|[-_])ALS([-_]|$)",
        annotated_cells$sample,
        ignore.case = TRUE
    ),
    "ALS",
    "control"
)

annotated_cells$condition <- factor(
    annotated_cells$condition,
    levels = c("control", "ALS")
)

table(unique(annotated_cells[[]][,c("patient_id", "condition")])$condition)

library(dplyr)

cluster_composition <- annotated_cells[[]] |>
    count(
        seurat_clusters,
        predicted.cell_type,
        name = "n"
    ) |>
    group_by(seurat_clusters) |>
    mutate(
        total = sum(n),
        proportion = n / total
    ) |>
    arrange(
        as.numeric(as.character(seurat_clusters)),
        desc(proportion)
    ) |>
    ungroup()

cluster_composition

cluster_top2 <- cluster_composition |>
    group_by(seurat_clusters) |>
    slice_max(
        order_by = proportion,
        n = 2,
        with_ties = FALSE
    ) |>
    ungroup()

cluster_annotation <- cluster_composition |>
    group_by(seurat_clusters) |>
    slice_max(
        order_by = proportion,
        n = 1,
        with_ties = FALSE
    ) |>
    ungroup() |>
    select(
        seurat_clusters,
        cluster_cell_type = predicted.cell_type,
        cluster_purity = proportion,
        n,
        total
    )

cluster_annotation

cluster_to_celltype <- setNames(
    cluster_annotation$cluster_cell_type,
    cluster_annotation$seurat_clusters
)

annotated_cells$cluster_cell_type <- unname(
    cluster_to_celltype[
        as.character(annotated_cells$seurat_clusters)
    ]
)
colnames(patient_cell_type_df) <- c("condition", "cell_type", "frequency")
patient_condition_df <- unique(annotated_cells[[]][,c("patient_id", "condition")])
rownames(patient_condition_df) <- NULL
condition_freq <- table(patient_condition_df$condition)
patient_cell_type_df$per_patient_frequency <- patient_cell_type_df$frequency / condition_freq[patient_cell_type_df$condition]
patient_cell_type_df <- patient_cell_type_df[
    patient_cell_type_df$cell_type %in%
    c( "Neurons", "Astrocytes", "Oligodendrocytes",
    "Endothelial", "Lymphocytes", "OPC", "Microglia",
    "Meninges"),]

require(RColorBrewer)
require(viridis)
require(ggplot2)
condition_colours <- c("control" = "#F0EFF4", "ALS" = "#BD1E1E")

p_frequency <- ggplot(
    patient_cell_type_df,
    aes(
        x = reorder(
            cell_type,
            frequency
        ),
        y = log10(per_patient_frequency)
    )
) +
    geom_col(aes(fill = condition), color = "#000000", position = position_dodge(width = 0.9)) +
    scale_fill_manual("cell type", values = condition_colours) +
    geom_hline(
        yintercept = 0,
        linetype = 2
    ) +
    coord_flip() +
    labs(
        x = NULL,
        y = "log10 N of cells detected per patient",
        title = "Frequency of cell annotations"
    ) +
    theme_classic()

ggsave(
    filename = file.path(
        plot_dir,
        "20-01-plot-frequency-of-predicted-cell-types.pdf"
    ),
    plot = p_frequency,
    width = 10,
    height = 5,
    units = "in",
    dpi = 300
)

require(ggplot2)

a_plot <- DimPlot(
    annotated_cells,
    reduction = "umap",
    group.by = "cluster_cell_type",
    cols = cell_type_colours,
    label = TRUE,
    repel = TRUE,
    raster = TRUE
) + NoLegend()

ggsave(
    filename = file.path(
        plot_dir,
        "20-03-plot-unsupervised-clusters-from-seurat-with-cleaned-up-reference-annotation.pdf"
    ),
    plot = a_plot,
    width = 15,
    height = 15,
    units = "in",
    dpi = 300
)

library(dplyr)
library(ggplot2)

counts <- LayerData(
    annotated_cells,
    assay = "RNA",
    layer = "counts"
)

meta <- annotated_cells[[]]
meta$ERAP2_count <- as.numeric(
    counts["ERAP2", rownames(meta)]
)

meta$total_counts <- Matrix::colSums(
    counts[, rownames(meta), drop = FALSE]
)

erap2_patient_cpm <- meta |>
    group_by(
        patient_id,
        condition,
        cluster_cell_type
    ) |>
    summarise(
        ERAP2_counts = sum(ERAP2_count),
        total_counts = sum(total_counts),

        n_cells = n(),

        n_ERAP2_positive = sum(ERAP2_count > 0),
        pct_expression = 100 * mean(ERAP2_count > 0),

        .groups = "drop"
    ) |>
    mutate(
        ERAP2_CPM =
            ERAP2_counts / total_counts * 1e6
    )

erap2_summary <- erap2_patient_cpm |>
    group_by(
        cluster_cell_type,
        condition
    ) |>
    summarise(
        mean_CPM = mean(ERAP2_CPM),
        median_CPM = median(ERAP2_CPM),

        mean_pct_expression = mean(pct_expression),

        n_patients = n_distinct(patient_id),

        .groups = "drop"
    ) |>
    mutate(
        ERAP2_condition = paste("ERAP2", condition)
    )

erap2_summary$ERAP2_condition <- factor(
    erap2_summary$ERAP2_condition,
    levels = c(
        "ERAP2 ALS",
        "ERAP2 control"
    )
)

p_erap2_condition <- ggplot(
    erap2_summary,
    aes(
        x = cluster_cell_type,
        y = ERAP2_condition,
        size = mean_pct_expression,
        colour = mean_CPM
    )
) +
    geom_point() +
    scale_size(
        name = "Percent expressed",
        range = c(1, 8)
    ) +
    scale_colour_viridis(
        option = "A",
        name = "ERAP2 CPM"
    ) +
    labs(
        x = NULL,
        y = NULL,
        title = "ERAP2 expression by cell type and condition"
    ) +
    theme_classic() +
    theme(
        axis.text.x = element_text(
            angle = 90,
            hjust = 0.5,
            vjust = 1
        )
    )

ggsave(
    filename = file.path(
        plot_dir,
        "20-06-plot-interaction-effect-between-ERAP2-and-condition-in-cleaned-up-reference-annotation.pdf"
    ),
    plot = p_erap2_condition,
    width = 10,
    height = 5,
    units = "in",
    dpi = 300
)
DefaultAssay(annotated_cells) <- "RNA"

counts <- LayerData(
    annotated_cells,
    assay = "RNA",
    layer = "counts"
)

meta <- annotated_cells[[]]

meta$ERAP2_positive <- as.numeric(
    counts["ERAP2", rownames(meta)]
) > 0

library(dplyr)

erap2_patient <- meta |>
    group_by(
        patient_id,
        condition,
        cluster_cell_type
    ) |>
    summarise(
        n_cells = n(),
        n_ERAP2_positive = sum(ERAP2_positive),
        frac_ERAP2_positive = mean(ERAP2_positive),
        .groups = "drop"
    )

erap2_patient |>
    arrange(cluster_cell_type, patient_id)

library(tidyr)

erap2_enrichment <- erap2_patient |>
    group_by(
        cluster_cell_type,
        condition
    ) |>
    summarise(
        mean_frac = mean(frac_ERAP2_positive),
        median_frac = median(frac_ERAP2_positive),
        n_patients = n_distinct(patient_id),
        .groups = "drop"
    ) |>
    pivot_wider(
        names_from = condition,
        values_from = c(
            mean_frac,
            median_frac,
            n_patients
        )
    ) |>
    mutate(

        delta_percent =
            100 * (mean_frac_ALS - mean_frac_control),

        fold_enrichment =
            mean_frac_ALS / mean_frac_control,

        log2_enrichment =
            log2(fold_enrichment)
    )

p_enrichment <- ggplot(
    erap2_enrichment,
    aes(
        x = reorder(
            cluster_cell_type,
            log2_enrichment
        ),
        y = log2_enrichment
    )
) +
    geom_col(aes(fill = cluster_cell_type)) +
    scale_fill_manual("cell type", values = cell_type_colours) +
    geom_hline(
        yintercept = 0,
        linetype = 2
    ) +
    coord_flip() +
    labs(
        x = NULL,
        y = "ERAP2 detection enrichment\nlog2(ALS / control)",
        title = "Relative enrichment of ERAP2-positive cells"
    ) +
    theme_classic()

ggsave(
    filename = file.path(
        plot_dir,
        "20-07-plot-enrichment-of-ERAP2-expressing-cells-cleaned-up-reference-annotation.pdf"
    ),
    plot = p_enrichment,
    width = 10,
    height = 5,
    units = "in",
    dpi = 300
)

library(lme4)

meta$condition <- factor(
    meta$condition,
    levels = c("control", "ALS")
)

fit_erap2_detection <- glmer(
    ERAP2_positive ~
        condition * cluster_cell_type +
        (1 | patient_id),
    data = meta,
    family = binomial
)

library(emmeans)

emm <- emmeans(
    fit_erap2_detection,
    ~ condition | cluster_cell_type
)

erap2_OR <- contrast(
    emm,
    method = list(
        ALS_vs_control = c(-1, 1)
    )
)

summary(
    erap2_OR,
    infer = TRUE,
    type = "response"
)

erap2_OR_df <- as.data.frame(
    summary(
        erap2_OR,
        infer = TRUE,
        type = "response"
    )
)

library(ggplot2)

p_or_plot <- ggplot(
    erap2_OR_df,
    aes(
        y = reorder(
            cluster_cell_type,
            odds.ratio
        ),
        x = odds.ratio
    )
) +
    geom_errorbarh(
        aes(
            xmin = asymp.LCL,
            xmax = asymp.UCL
        ),
        height = 0.15
    ) +
    geom_point(size = 3, shape = 21, color = "#000000", aes(fill = cluster_cell_type)) +
    geom_vline(
        xintercept = 1,
        linetype = 2
    ) +
    scale_x_log10() +
    scale_fill_manual("cell type", values = cell_type_colours) +
    labs(
        x = "ALS / control odds ratio for ERAP2 detection",
        y = NULL,
        title = "ERAP2-positive cell enrichment in ALS"
    ) +
    theme_classic()

ggsave(
    filename = file.path(
        plot_dir,
        "20-08-plot-odds-ratio-of-ERAP2-expressing-cells-cleaned-up-reference-annotation.pdf"
    ),
    plot = p_or_plot,
    width = 5,
    height = 10,
    units = "in",
    dpi = 300
)
genic_count_list <- readRDS(input("bulk_counts.rds"))
wgcna_modules_df <- readRDS(input("bulk_modules.rds"))

erap2_module_df <- wgcna_modules_df[wgcna_modules_df$modules == "M3",]
gene_list <- erap2_module_df$gene_name

gene_count_matrix <- genic_count_list$counts
cpm_matrix <- t(t(gene_count_matrix) / (colSums(gene_count_matrix)/1e6))

coldata_df <- genic_count_list$coldata
coldata_df$expr <- cpm_matrix[gene_ranges$gene_id[gene_ranges$gene_name == "ERAP2" & !is.na(gene_ranges$gene_name)],
    coldata_df$column_id]
coldata_df$tissue_group <- "spinal_cord"
coldata_df$tissue_group[coldata_df$tissue %in% c("Cortex_Motor", "Cortex", "Cerebellum")] <- "brain"

subset_matrix <- as.data.frame(cpm_matrix[rownames(cpm_matrix) %in%
    gene_ranges$gene_id[gene_ranges$gene_biotype == "protein_coding"],])
keep_matrix <- subset_matrix[(rowSums(subset_matrix == 0) / ncol(subset_matrix)) < 0.10,]

coldata_df$matrix_row <- match(coldata_df$column_id, colnames(keep_matrix))

cols_split_by_patient <- split(coldata_df, paste(coldata_df$subject_id, coldata_df$tissue))

patient_matrix <- do.call(cbind, lapply(cols_split_by_patient, function(x){
    if(nrow(x) == 1){
        return(keep_matrix[,x$matrix_row])
    }
    return(rowMeans(keep_matrix[,x$matrix_row]))
}))
selected_patient_matrix <- patient_matrix[erap2_module_df$genes,]

patient_df <- do.call(rbind, lapply(cols_split_by_patient, function(temp_df){
    return_df <- unique(temp_df[,!(colnames(temp_df) %in%
        c("X", "sample_name", "alignment_file", "alignment_file_dir",
          "isPairedEnd", "strand_specificity", "samplesheet_path",
          "column_id", "sample_id", "mutation", "site_of_motor_onset",
          "instrument", "rin", "expr", "tissue_group", "matrix_row"))])
    if(nrow(return_df) > 1){
        if("N/A" %in% return_df$mnd_with_ftd & "Yes" %in% return_df$mnd_with_ftd){
            return_df <- return_df[return_df$mnd_with_ftd == "Yes",]
        }else if("N/A" %in% return_df$mnd_with_ftd & "No" %in% return_df$mnd_with_ftd){
            return_df <- return_df[return_df$mnd_with_ftd == "No",]
        }else{
            stop(paste(return_df$subject_id, return_df$tissue))
        }
    }
    return_df$column_id <- paste(return_df$subject_id, return_df$tissue)
    return(return_df)
}))
patient_df$age_at_death[patient_df$age_at_death == "90 or Older"] <- "90"
patient_df$age_at_death <- as.numeric(patient_df$age_at_death)
patient_df$age_at_onset <- as.numeric(patient_df$age_at_symptom_onset)
patient_df$disease_duration <- as.numeric(patient_df$disease_duration_in_months)
patient_df$disease_duration_in_months <- as.numeric(patient_df$disease_duration_in_months)

cervical_patient_df <- patient_df[patient_df$tissue == "Spinal_Cord_Cervical",]
cervical_selected_matrix <- selected_patient_matrix[,cervical_patient_df$column_id]
rownames(cervical_selected_matrix) <- erap2_module_df$gene_name
cervical_selected_matrix <- cervical_selected_matrix[!is.na(rownames(cervical_selected_matrix)),]

X_bulk <- t(cervical_selected_matrix)

X_bulk_z <- scale(X_bulk)

train_center <- attr(X_bulk_z, "scaled:center")
train_scale  <- attr(X_bulk_z, "scaled:scale")

pca_module <- prcomp(
    X_bulk_z,
    center = FALSE,
    scale. = FALSE
)
module_weights <- pca_module$rotation[, 1]

sort(
    module_weights,
    decreasing = TRUE
)

nygc_module_score <- pca_module$x[, 1]

erap2_bulk <- cervical_selected_matrix[
    "ERAP2",
    rownames(X_bulk)
]

if (
    cor(
        nygc_module_score,
        erap2_bulk,
        method = "spearman"
    ) < 0
) {
    module_weights <- -module_weights
    nygc_module_score <- -nygc_module_score
}

cell_types <- c("Astrocytes", "Lymphocytes", "Microglia", "Neurons", "Oligodendrocytes")
patient_cell_type_scores_df_list <- lapply(cell_types, function(a_cell_type){
    selected_cells <- subset(
        annotated_cells,
        subset = cluster_cell_type == a_cell_type
    )
    patient_meta <- selected_cells[[]] |>
        distinct(
            patient_id,
            condition
        )
    patient_meta
    pb_counts <- AggregateExpression(
        selected_cells,
        assays = "RNA",
        group.by = "patient_id",
        return.seurat = FALSE
    )$RNA
    pb_cpm <- sweep(
        pb_counts,
        2,
        Matrix::colSums(pb_counts),
        "/"
    ) * 1e6

    genes_pb <- intersect(
        names(module_weights),
        rownames(pb_cpm)
    )

    length(module_weights)
    length(genes_pb)

    "ERAP2" %in% genes_pb

    X_pb <- t(
        pb_cpm[
            genes_pb,
            ,
            drop = FALSE
        ]
    )
    pb_sd <- apply(
        X_pb,
        2,
        sd,
        na.rm = TRUE
    )

    keep_pb <- is.finite(pb_sd) & pb_sd > 0

    X_pb <- X_pb[
        ,
        keep_pb,
        drop = FALSE
    ]

    X_pb_z <- scale(X_pb)

    weights_pb <- module_weights[
        colnames(X_pb_z)
    ]

    stopifnot(
        identical(
            names(weights_pb),
            colnames(X_pb_z)
        )
    )

    weights_pb <- weights_pb /
        sqrt(sum(weights_pb^2))

    pb_module_score <- as.numeric(
        X_pb_z %*% weights_pb
    )

    names(pb_module_score) <- rownames(X_pb_z)

    patient_scores_df <- data.frame(
        cell_type = a_cell_type,
        patient_id = stringr::str_replace_all(string = names(pb_module_score), pattern = "g", replacement = ""),
        ERAP2_module_score = pb_module_score
    ) |>
        left_join(
            patient_meta,
            by = "patient_id"
        )

    gene_contrib <- sweep(
        X_pb_z,
        2,
        weights_pb,
        "*"
    )

    condition_vector <- patient_scores_df$condition[
        match(
            rownames(gene_contrib),
            paste("g",patient_scores_df$patient_id, sep = "")
        )
    ]

    delta_contrib <- apply(
        gene_contrib,
        2,
        function(x) {

            mean(
                x[condition_vector == "ALS"],
                na.rm = TRUE
            ) -
            mean(
                x[condition_vector == "control"],
                na.rm = TRUE
            )
        }
    )

    delta_df <- data.frame(
        cell_type = a_cell_type,
        gene = names(delta_contrib),
        delta_contribution = delta_contrib)

    delta_df <- delta_df[order(abs(delta_contrib), decreasing = F),]
    delta_df$rank <- seq_len(nrow(delta_df))

    return(list(deltas = delta_df, scores = patient_scores_df))
})
names(patient_cell_type_scores_df_list) <- cell_types

patient_cell_type_scores_df <- do.call(rbind, lapply(patient_cell_type_scores_df_list, function(x){
    x$scores
}))
rownames(patient_cell_type_scores_df) <- NULL

require(dplyr)
require(broom)

wilcox_test_tib <- patient_cell_type_scores_df |>
group_by(cell_type) |>
group_modify(~ broom::tidy(
    wilcox.test(ERAP2_module_score ~ condition, data = .x)
)) |>
ungroup() |>
mutate(padj = p.adjust(p.value, method = "BH"))

as.data.frame(wilcox_test_tib)

patient_cell_type_scores_df$cell_type <- factor(patient_cell_type_scores_df$cell_type,
levels = wilcox_test_tib$cell_type[order(wilcox_test_tib$padj, decreasing = F)])

library(ggplot2)

p_module <- ggplot(
    patient_cell_type_scores_df,
    aes(
        x = condition,
        y = ERAP2_module_score
    )
) +
    geom_boxplot(
        aes(fill = cell_type),
        color = "#000000",
        outlier.shape = NA,
        width = 0.5
    ) +
    geom_jitter(
        aes(fill = cell_type),
        color = "#000000",
        shape = 21,
        width = 0.08,
        size = 3
    ) +
    scale_fill_manual("cell types", values = cell_type_colours[unique(as.character(patient_cell_type_scores_df$cell_type))]) +
    facet_wrap(cell_type ~ ., scales = "free", nrow = 1, ncol = 5) +
    labs(
        x = NULL,
        y = "ERAP2-associated module score",
        title = "ERAP2-associated programme in astrocytes"
    ) +
    theme_classic()

ggsave(
    filename = file.path(
        plot_dir,
        "20-09-plot-module-score-of-ERAP2-module-M3-in-scrnaseq-data.pdf"
    ),
    plot = p_module,
    width = 10,
    height = 5,
    units = "in",
    dpi = 300
)

writeLines(capture.output(sessionInfo()), file.path(plot_dir, "sessionInfo.txt"))
