summarise_tissue_specific_parameter_by_patients <- function(patient_df){
	main_return_df <- do.call(rbind, lapply(split(patient_df, patient_df$subject_id), function(temp_df){
	    message(unique(temp_df$subject_id))
	    return_df <- data.frame(
	        subject_id = unique(temp_df$subject_id),
	        age_at_onset = unique(temp_df$age_at_symptom_onset),
	        disease_duration = unique(temp_df$disease_duration),
	        age_at_death = unique(temp_df$age_at_death),
	        als_group = unique(temp_df$condition),
	        sex_genotype = unique(temp_df$sex_genotype))
	    	all_tissues <- c("Spinal_Cord_Lumbar", "Spinal_Cord_Cervical", "Spinal_Cord_Thoracic", "Cerebellum", "Cortex", "Cortex_Motor")
	    	tissue_expr_df <- do.call(cbind, lapply(all_tissues, function(x){
	    		return_df <- data.frame(expr = NA)
	    		if(x %in% temp_df$tissue){
		    		return_df <- data.frame(expr = mean(temp_df$expr[temp_df$tissue == x]))
	    		}
	    		colnames(return_df) <- paste(x, "expr", sep = "_")
	    		return(return_df)
	    	}))
	    	return(cbind(return_df, tissue_expr_df))
	}))
	return(main_return_df)
}


prepare_data_dose_response_curve <- function(x, value_col, n_combs = 10){
	value_range <- seq(min(x[,value_col]), max(x[,value_col]), length.out = n_combs)
	return_df <- do.call(rbind, lapply(value_range, function(y){
	    prop <- nrow(x[x[,value_col] < y,]) / nrow(x)
	    temp_df <- data.frame(value = y, prop = prop)
	    colnames(temp_df) <- c("value", "proportions")
	    return(temp_df)
	}))
    return(return_df)
}

# return_tissue_specific_survival_df <- function(patient_df, tissue, value_col, value){
# 	current_tissue_df <- patient_df[patient_df[,"tissue"] == tissue,]
# 	current_tissue_df$event <- ifelse(current_tissue_df[,value_col] < value)

# }

summarise_parameter_by_patients <- function(patient_df, params, fun = mean){
	main_return_df <- do.call(rbind, lapply(split(patient_df, patient_df$subject_id), function(temp_df){
	    summary_matrix <- do.call(cbind, lapply(params, function(x){
	    	fun(temp_df[,x])
	    }))
	    summary_df <- as.data.frame(summary_matrix)
	    colnames(summary_df) <- params
	    return_df <- data.frame(subject_id = unique(temp_df$subject_id),
	        age_at_onset = unique(temp_df$age_at_onset),
	        disease_duration = unique(temp_df$disease_duration),
	        age_at_death = unique(temp_df$age_at_death),
	        als_group = unique(temp_df$als_group))
	    return_df <- cbind(return_df, summary_df)
	    return(return_df)
	}))
	return(main_return_df)
}

summarise_pvalues_parameter_by_patients <- function(patient_df){
	main_return_df <- do.call(rbind, lapply(split(patient_df, paste(patient_df$subject_id, patient_df$tissue)), function(temp_df){
	    # message(unique(temp_df$subject_id))
	    data.frame(
	    	subject_id = unique(temp_df$subject_id),
	        age_at_onset = unique(temp_df$age_at_onset),
	        disease_duration = unique(temp_df$disease_duration),
	        age_at_death = unique(temp_df$age_at_death),
	        tissue = unique(temp_df$tissue),
	        gene_ids = names(temp_df$sig_rep_vs_exon[[1]]),
	        exonic_pvalues = temp_df$sig_rep_vs_exon[[1]],
	        intronic_pvalues = temp_df$sig_rep_vs_intron[[1]],
	        als_group = unique(temp_df$als_group))
	}))
	return(main_return_df)
}



return_tissue_specific_event_df <- function(patient_df, tissue, region, group_val){
	prop_cols <- c("exon" = "exonic_prop",
		"intron" = "intronic_prop",
		"both" = "both_prop")
	tissue_cols <- c("Spinal_Cord_Lumbar" = "lumbar",
		"Spinal_Cord_Cervical" = "cervical",
		"Spinal_Cord_Thoracic" = "thoracic",
		"Cerebellum" = "cerebellum",
		"Cortex" = "cortex",
		"Cortex_Motor" = "motor_cortex")
	group_cols <- c("exon" = "e_r", "intron" = "i_r", "both" = "b_r")
	group_labs <- c("exon" = "e", "intron" = "i", "both" = "b")
	current_col <- paste(tissue_cols[tissue], prop_cols[region], sep = "_")
	patient_df$event <- ifelse(!is.na(patient_df[,current_col]), 1, 0)
	patient_df[,group_cols[region]] <- paste(group_labs[region], 
		ifelse(!is.na(patient_df[,current_col]) & patient_df[,current_col] > group_val, "+", "-"), sep = "")
	patient_df[is.na(patient_df[,current_col]),group_cols[region]] <- "missing"
	return(patient_df)
}


return_tissue_specific_landmark_event_df <- function(patient_df, tissue, region, group_val, landmark_duration = 30){
	prop_cols <- c("exon" = "exonic_prop",
		"intron" = "intronic_prop",
		"both" = "both_prop")
	tissue_cols <- c("Spinal_Cord_Lumbar" = "lumbar",
		"Spinal_Cord_Cervical" = "cervical",
		"Spinal_Cord_Thoracic" = "thoracic",
		"Cerebellum" = "cerebellum",
		"Cortex" = "cortex",
		"Cortex_Motor" = "motor_cortex")
	group_cols <- c("exon" = "e_r", "intron" = "i_r", "both" = "b_r")
	group_labs <- c("exon" = "e", "intron" = "i", "both" = "b")
	current_col <- paste(tissue_cols[tissue], prop_cols[region], sep = "_")
	current_df <- patient_df[!is.na(patient_df[,"disease_duration"]) & patient_df[,"disease_duration"] <= landmark_duration,]
	current_df$event <- ifelse(!is.na(current_df[,current_col]), 1, 0)
	current_df[,group_cols[region]] <- paste(group_labs[region], 
		ifelse(!is.na(current_df[,current_col]) & current_df[,current_col] > group_val, "+", "-"), sep = "")
	current_df[is.na(current_df[,current_col]),group_cols[region]] <- "missing"
	return(current_df)
}




return_tissue_specific_interval_event_df <- function(patient_df, tissue, region, group_val){
	prop_cols <- c("exon" = "exonic_prop",
		"intron" = "intronic_prop",
		"both" = "both_prop")
	tissue_cols <- c("Spinal_Cord_Lumbar" = "lumbar",
		"Spinal_Cord_Cervical" = "cervical",
		"Spinal_Cord_Thoracic" = "thoracic",
		"Cerebellum" = "cerebellum",
		"Cortex" = "cortex",
		"Cortex_Motor" = "motor_cortex")
	group_cols <- c("exon" = "e_r", "intron" = "i_r", "both" = "b_r")
	group_labs <- c("exon" = "e", "intron" = "i", "both" = "b")
	current_col <- paste(tissue_cols[tissue], prop_cols[region], sep = "_")
	patient_df$event <- ifelse(!is.na(patient_df[,"age_at_onset"]), 1, 2)
	patient_df$event[is.na(patient_df[,"disease_duration"])] <- 3
	patient_df$event[is.na(patient_df[,"disease_duration"]) & is.na(patient_df[,"age_at_onset"])] <- 2
	patient_df[,group_cols[region]] <- paste(group_labs[region], 
		ifelse(!is.na(patient_df[,current_col]) & patient_df[,current_col] > group_val, "+", "-"), sep = "")
	patient_df[is.na(patient_df[,current_col]),group_cols[region]] <- "missing"
	return(patient_df)
}