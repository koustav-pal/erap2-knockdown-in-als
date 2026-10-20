# Adapted analysis; see provenance/source_ranges.json.
source(here::here("workflow", "scripts", "common.R"))
library(stringr)
library(survival)
library(survminer)
library(ggbeeswarm)
library(patchwork)
library(randomForestSRC)
library(matrixStats)
patient_count_df <- readRDS(input("patient_expression.rds"))
patient_count_df$type <- factor(patient_count_df$type, levels = c("control", "c9orf72-ALS", "other_mutation-fALS", "sporadic"))
patient_count_df$type_group <- as.character(patient_count_df$type)
patient_count_df$type_group[patient_count_df$type %in% c("other_mutation-fALS", "sporadic")] <- "sporadic+other_fALS"
patient_count_df$tissue_group <- "spinal_cord"
patient_count_df$tissue_group[patient_count_df$tissue %in% c("Cortex", "Cortex_Motor", "Cerebellum")] <- "brain"

patient_surv_df <- do.call(rbind, lapply(split(patient_count_df, paste(patient_count_df$subject_id, patient_count_df$tissue)), function(temp_df){
    data.frame(subject_id = unique(temp_df$subject_id),
        tissue = unique(temp_df$tissue),
        age_at_onset = unique(temp_df$age_at_symptom_onset),
        age_at_death = unique(temp_df$age_at_death),
        disease_duration_in_months = unique(temp_df$disease_duration_in_months),
        sex_genotype = unique(temp_df$sex_genotype),
        mnd_with_ftd = unique(temp_df$mnd_with_ftd),
        condition = unique(temp_df$condition),
        mnd_with_dementia = unique(temp_df$mnd_with_dementia),
        type_group = unique(temp_df$type_group),
        erap2_expression = mean(temp_df$expr))
}))

patient_surv_df$age_at_death[patient_surv_df$age_at_death == "90 or Older"] <- "90"
patient_surv_df$age_at_death <- as.numeric(patient_surv_df$age_at_death)
patient_surv_df$age_at_onset <- as.numeric(patient_surv_df$age_at_onset)
patient_surv_df$disease_duration <- as.numeric(patient_surv_df$disease_duration_in_months)

als_surv_df <- patient_surv_df[patient_surv_df$condition %in% c("als"),]
als_surv_df$sex_genotype <- factor(als_surv_df$sex_genotype, levels = c("XX", "XY"))

als_surv_df$status <- 1
als_surv_df$status[is.na(als_surv_df$age_at_death)] <- 0

spinal_cord_df <- als_surv_df[als_surv_df$tissue %in% c("Spinal_Cord_Cervical", "Spinal_Cord_Lumbar", "Spinal_Cord_Thoracic"),c("subject_id", "tissue", "age_at_onset", "age_at_death", "disease_duration", "status", "condition", "sex_genotype", "erap2_expression")]

mutation_list <- readRDS(input("patient_mutations.rds"))

patient_specific_mutation_list <- mutation_list$patient_mutation_list
frequency_matrix <- mutation_list$frequency_matrix
selected_mutation_df <- mutation_list$mutation_list

candidate_mutations_df <- selected_mutation_df[selected_mutation_df$als > 0.1 & selected_mutation_df$ctrl < 0.05,]

keep_patient_specific_mutation_list <- patient_specific_mutation_list[names(patient_specific_mutation_list) %in% spinal_cord_df$subject_id]

patient_mutation_frequency_df <- t(do.call(cbind,lapply(keep_patient_specific_mutation_list, function(x){ x$occurence_vector[candidate_mutations_df$snp_id] })))
colnames(patient_mutation_frequency_df) <- paste("rs", candidate_mutations_df$snp_id, sep = "")

spinal_cord_with_mutation_df <- cbind(spinal_cord_df, patient_mutation_frequency_df[spinal_cord_df$subject_id,])

seed_list <- readRDS(input("rsf_seeds.rds"))
stopifnot(length(seed_list) == 1000L, !anyDuplicated(seed_list))
text_formula <- formula(paste("Surv(disease_duration, status)", paste("age_at_onset", "sex_genotype", "erap2_expression", paste(colnames(patient_mutation_frequency_df), collapse = " + "), sep = " + "), sep = " ~ "))

all_importance <- list()
risk_data <- list()
for (region in c("Cervical", "Lumbar")) {
    d <- spinal_cord_with_mutation_df[spinal_cord_with_mutation_df$tissue == paste0("Spinal_Cord_", region),]
    predict_data <- d[!is.na(d$disease_duration) & !is.na(d$age_at_onset), ]
    risks <- matrix(NA_real_, nrow(predict_data), length(seed_list))
    imp <- vector("list", length(seed_list))
    for (i in seq_along(seed_list)) {
        set.seed(seed_list[[i]])
        model <- rfsrc(text_formula, data = d, ntree = 1000, importance = TRUE, nthread = 1)
        risks[,i] <- predict(model, newdata = predict_data)$predicted
        imp[[i]] <- data.frame(tissue = tolower(region), seed = seed_list[[i]], variable = names(model$importance), values = as.numeric(model$importance))
    }
    d <- predict_data
    d$risk_score <- matrixStats::rowMedians(risks)
    d$risk_group <- ifelse(d$risk_score >= median(d$risk_score, na.rm = TRUE), "high", "low")
    risk_data[[tolower(region)]] <- d
    all_importance[[tolower(region)]] <- do.call(rbind, imp)
}
all_importance_df <- do.call(rbind, all_importance)
saveRDS(risk_data, file.path(plot_dir, "risk_stratification.rds"))
write.csv(all_importance_df, file.path(plot_dir, "variable_importance.csv"), row.names = FALSE)
for (region in names(risk_data)) {
    d <- risk_data[[region]]
    fit <- survfit(Surv(disease_duration, status) ~ risk_group, data = d)
    surv <- ggsurvplot(fit, data = d, palette = c("#FF5666", "#FFCCC9"), pval = TRUE, pval.method = TRUE,
                       conf.int = FALSE, ggtheme = theme_bw())$plot + ggtitle(region)
    expr <- ggplot(d, aes(risk_group, erap2_expression, fill = risk_group)) + geom_boxplot(width = 0.1, outlier.shape = NA) +
      geom_quasirandom(shape = 21) + scale_fill_manual(values = c(high = "#FF5666", low = "#FFCCC9")) + my_theme
    ggsave(file.path(plot_dir, paste0("Figure_3D_", region, ".pdf")), surv | expr, width = 10, height = 5)
    v <- subset(all_importance_df, tissue == region)
    v$variable <- reorder(v$variable, v$values, median)
    importance_plot <- ggplot(v, aes(variable, values)) + geom_violin(fill = "#F1FAEE") +
        geom_boxplot(width = 0.1, outlier.shape = NA) + my_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(file.path(plot_dir, paste0("Figure_S6_", region, "_importance.pdf")), importance_plot, width = 12, height = 6)
}

writeLines(capture.output(sessionInfo()), file.path(plot_dir, "sessionInfo.txt"))
