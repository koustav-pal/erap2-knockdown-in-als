data_dir <- normalizePath("/data/shared/lab_internal/")
project_dir <- file.path(data_dir, "012-Ben-Clarke-ERAP1-KD")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
# save_dir <- file.path(project_dir, "rdata_files")
sampleinfo_dir <- file.path(project_dir, "sampleinfo")
analysis_dir <- file.path(project_dir, "analysis", "99-02-autophagy-erap2-kd")
setwd(analysis_dir)
analysis_save_dir <- file.path(analysis_dir, "rdata_files")
analysis_input_dir <- file.path(analysis_dir, "input_files")
analysis_plots_dir <- file.path(analysis_dir, "plots")
LC3_cell_level_list <- readRDS(file.path(analysis_save_dir, "LC3", "all-CLO-metrics-binned.rds"))

# ==========================================================================
# Imports
# ==========================================================================

library("stringr")
library("dplyr")
library("tidyr")
library("glmmTMB")

# ==========================================================================
# Setting variables as factors with levels
# ==========================================================================

LC3_cell_level_list$well_uid <- with(
  LC3_cell_level_list,
  interaction(experiment, well_id, drop = TRUE)
)

LC3_cell_level_list$cell_line <- factor(
    LC3_cell_level_list$cell_line,
    levels = unique(LC3_cell_level_list$cell_line)
)

# ==========================================================================
# Model formulae and design construction
# ==========================================================================

bin_vars <- c(
  "CLO_spots_bin_1st",
  "CLO_spots_bin_2nd",
  "CLO_spots_bin_3rd",
  "CLO_spots_bin_4th"
)

# Model v1_01
m_v1_01_model_formula <- setNames(
  lapply(bin_vars, function(var) {
    model_formula <- reformulate(
      termlabels = c(
        "condition * treatment * bafa1_timepoint",
        "(1 | cell_line)",
        "(1 | experiment)",
        "(1 | cell_line:experiment)",
        "(1 | well_uid)"
      ),
      response = var
    )
  }),
  bin_vars
)

# Model v1_02 (incl. coverage_label as random effect)
LC3_cell_level_list$coverage_label <- "high"
LC3_cell_level_list$coverage_label[LC3_cell_level_list$CLO_spot_count < 50] <- "medium"
LC3_cell_level_list$coverage_label[LC3_cell_level_list$CLO_spot_count < 6] <- "low"
LC3_cell_level_list$coverage_label <- factor(LC3_cell_level_list$coverage_label, levels = c("low", "medium", "high"))

m_v1_02_model_formula <- setNames(
  lapply(bin_vars, function(var) {
    model_formula <- reformulate(
      termlabels = c(
        "condition * treatment * bafa1_timepoint",
        "(1 | cell_line)",
        "(1 | experiment)",
        "(1 | cell_line:experiment)",
        "(1 | well_uid)",
        "(1 | coverage_label)"
      ),
      response = var
    )
  }),
  bin_vars
)

# Model v1_03
m_v1_03_model_formula <- setNames(
  lapply(bin_vars, function(var) {
    model_formula <- reformulate(
      termlabels = c(
        "condition * treatment",
        "(1 | cell_line)",
        "(1 | experiment)",
        "(1 | cell_line:experiment)",
        "(1 | well_uid)",
        "(1 | bafa1_timepoint)"
      ),
      response = var
    )
  }),
  bin_vars
)

# Model v1_04
m_v1_04_model_formula <- setNames(
  lapply(bin_vars, function(var) {
    model_formula <- reformulate(
      termlabels = c(
        "condition * treatment",
        "(1 | cell_line)",
        "(1 | experiment)",
        "(1 | cell_line:experiment)",
        "(1 | well_uid)",
        "(1 | coverage_label)",
        "(1 | bafa1_timepoint)"
      ),
      response = var
    )
  }),
  bin_vars
)

# Model v1_05
m_v1_05_model_formula <- setNames(
    lapply(bin_vars, function(var) {
        model_formula <- reformulate(
            termlabels = c(
                "condition * treatment * bafa1_timepoint",
                "scale(log(CLO_cyto_area_um2))",
                "(1 | cell_line)",
                "(1 | experiment)",
                "(1 | cell_line:experiment)",
                "(1 | well_uid)"
            ),
            response = var
        )
    }),
    bin_vars
)

# Model v1_06
m_v1_06_model_formula <- setNames(
    lapply(bin_vars, function(var) {
        model_formula <- reformulate(
            termlabels = c(
                "condition * treatment * bafa1_timepoint",
                "scale(log(CLO_cyto_area_um2))",
                "(1 | cell_line)",
                "(1 | experiment)",
                "(1 | cell_line:experiment)",
                "(1 | well_uid)",
                "(1 | coverage_label)"
            ),
            response = var
        )
    }),
    bin_vars
)

model_formulae_v1 <- list(
    m_v1_01_model_formula,
    m_v1_02_model_formula,
    m_v1_03_model_formula,
    m_v1_04_model_formula,
    m_v1_05_model_formula,
    m_v1_06_model_formula
)
model_names <- c(
        "m_v1_01",
        "m_v1_02",
        "m_v1_03",
        "m_v1_04",
        "m_v1_05",
        "m_v1_06"
)
model_formulae_v1 <- setNames(
    object = model_formulae_v1,
    model_names
)

# ==========================================================================
# Model fitting
# ==========================================================================

models_v1_nbd2 <- lapply(seq_along(model_names), function(x) {
    message("Fitting negative binomial v1 models. Current model: ", x)

    lapply(bin_vars, function(var) {
        glmmTMB::glmmTMB(
            formula = model_formulae_v1[[x]][[var]],
            data = LC3_cell_level_list,
            family = glmmTMB::nbinom2(link = "log")
        )
    })
})
models_v1_nbd2 <- setNames(
    object = models_v1_nbd2,
    model_names
)

models_v1_genpois <- lapply(seq_along(model_names), function(x) {
    message("Fitting genpois v1 models. Current model: ", x)

    lapply(bin_vars, function(var) {
        glmmTMB::glmmTMB(
            formula = model_formulae_v1[[x]][[var]],
            data = LC3_cell_level_list,
            family = glmmTMB::genpois(link = "log")
        )
    })
})
models_v1_genpois <- setNames(
    object = models_v1_genpois,
    model_names
)

models_v1_poisson <- lapply(seq_along(model_names), function(x) {
    message("Fitting poisson v1 models. Current model: ", x)

    lapply(bin_vars, function(var) {
        glmmTMB::glmmTMB(
            model_formulae_v1[[x]][[var]],
            data = LC3_cell_level_list,
            family = poisson(link = "log")
        )
    })
})
models_v1_poisson <- setNames(
    object = models_v1_poisson,
    model_names
)

# ==========================================================================
# Model evaluation (AIC/BIC/logLik)
# ==========================================================================

# Run eval for every model type
models_v1_nbd_eval <- lapply(
  models_v1_nbd,
  function(model_group) {

    model_names <- names(model_group)

    if (is.null(model_names)) {
      model_names <- paste0(
        "model_",
        seq_along(model_group)
      )
    }

    aic_table <- data.frame(
      model = model_names,

      n = vapply(
        model_group,
        nobs,
        numeric(1)
      ),

      parameters = vapply(
        model_group,
        function(m) attr(logLik(m), "df"),
        numeric(1)
      ),

      logLik = vapply(
        model_group,
        function(m) as.numeric(logLik(m)),
        numeric(1)
      ),

      AIC = vapply(
        model_group,
        AIC,
        numeric(1)
      ),

      BIC = vapply(
        model_group,
        BIC,
        numeric(1)
      ),

      row.names = NULL
    )

    # Avoid min(..., na.rm = TRUE) returning Inf when all AICs are NA
    if (all(is.na(aic_table$AIC))) {
      aic_table$delta_AIC <- NA_real_
    } else {
      aic_table$delta_AIC <-
        aic_table$AIC -
        min(aic_table$AIC, na.rm = TRUE)
    }

    aic_table
  }
)

models_v1_nbd2 <- list(
    formulae = model_formulae_v1,
    models = models_v1_nbd2,
    evaluation = models_v1_nbd2_eval
)

models_v1_genpois_eval <- lapply(
  models_v1_genpois,
  function(model_group) {

    model_names <- names(model_group)

    if (is.null(model_names)) {
      model_names <- paste0(
        "model_",
        seq_along(model_group)
      )
    }

    aic_table <- data.frame(
      model = model_names,

      n = vapply(
        model_group,
        nobs,
        numeric(1)
      ),

      parameters = vapply(
        model_group,
        function(m) attr(logLik(m), "df"),
        numeric(1)
      ),

      logLik = vapply(
        model_group,
        function(m) as.numeric(logLik(m)),
        numeric(1)
      ),

      AIC = vapply(
        model_group,
        AIC,
        numeric(1)
      ),

      BIC = vapply(
        model_group,
        BIC,
        numeric(1)
      ),

      row.names = NULL
    )

    # Avoid min(..., na.rm = TRUE) returning Inf when all AICs are NA
    if (all(is.na(aic_table$AIC))) {
      aic_table$delta_AIC <- NA_real_
    } else {
      aic_table$delta_AIC <-
        aic_table$AIC -
        min(aic_table$AIC, na.rm = TRUE)
    }

    aic_table
  }
)

models_v1_genpois <- list(
    formulae = model_formulae_v1,
    models = models_v1_genpois,
    evaluation = models_v1_genpois_eval
)

models_v1_poisson_eval <- lapply(
  models_v1_poisson,
  function(model_group) {

    model_names <- names(model_group)

    if (is.null(model_names)) {
      model_names <- paste0(
        "model_",
        seq_along(model_group)
      )
    }

    aic_table <- data.frame(
      model = model_names,

      n = vapply(
        model_group,
        nobs,
        numeric(1)
      ),

      parameters = vapply(
        model_group,
        function(m) attr(logLik(m), "df"),
        numeric(1)
      ),

      logLik = vapply(
        model_group,
        function(m) as.numeric(logLik(m)),
        numeric(1)
      ),

      AIC = vapply(
        model_group,
        AIC,
        numeric(1)
      ),

      BIC = vapply(
        model_group,
        BIC,
        numeric(1)
      ),

      row.names = NULL
    )

    # Avoid min(..., na.rm = TRUE) returning Inf when all AICs are NA
    if (all(is.na(aic_table$AIC))) {
      aic_table$delta_AIC <- NA_real_
    } else {
      aic_table$delta_AIC <-
        aic_table$AIC -
        min(aic_table$AIC, na.rm = TRUE)
    }

    aic_table
  }
)

models_v1_poisson <- list(
    formulae = model_formulae_v1,
    models = models_v1_poisson,
    evaluation = models_v1_poisson
)

# saveRDS(models_v1_nbd, file.path(analysis_save_dir, "LC3", "models-v1-all-nbd2.rds"))
# saveRDS(models_v1_genpois, file.path(analysis_save_dir, "LC3", "models-v1-all-genpois.rds"))
# saveRDS(models_v1_poisson, file.path(analysis_save_dir, "LC3", "models-v1-all-poisson.rds"))

# ==========================================================================
# Model evaluation (DHARMa simulated residuals plots)
# ==========================================================================

simulation_output <- DHARMa::simulateResiduals(
  fittedModel = m_v1_06[[1]],
  n = 1000
)

plot(simulation_output)

DHARMa::testUniformity(simulation_output)
DHARMa::testDispersion(simulation_output)
DHARMa::testZeroInflation(simulation_output)
DHARMa::testOutliers(simulation_output)