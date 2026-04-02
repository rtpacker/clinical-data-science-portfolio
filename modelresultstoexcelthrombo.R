# =============================================================================
# Residential Segregation & Thrombo-Inflammatory Biomarkers
# Regression Tables — Unweighted Linear Models
# =============================================================================
# Description: Linear regression models (lm) examining associations between
#              residential segregation indices and thrombo-inflammatory
#              biomarkers in the REGARDS cohort.
#
# Exposures:   Dissimilarity index, Interaction index, Isolation index
#              (each standardized by SD)
#
# Outcomes:    log_ddimer, log_crp, log_ifn, log_tnf, log_il6 (log-transformed),
#              eselectin, fix (untransformed)
#
# Models:      Unadjusted + 11 sequentially adjusted models
#
# NOTE: Update data.path and output.path before running.
# =============================================================================


# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------
library(dplyr)
library(data.table)
library(tidyr)
library(tidyverse)
library(haven)
library(openxlsx)


# -----------------------------------------------------------------------------
# Paths — update before running
# -----------------------------------------------------------------------------
data.path   <- "path/to/data/"
output.path <- "path/to/output/"


# =============================================================================
# Load Data & Derive Exposure Variables
# =============================================================================
analyticsamplethr <- readRDS(paste0(data.path, "analyticsamplethr.rds"))

# Standardize residential segregation indices (divide by SD)
analyticsamplethr$dissimilarity_index_std <- analyticsamplethr$rs_dissimilarityb / sd(analyticsamplethr$rs_dissimilarityb)
analyticsamplethr$interaction_index_std   <- analyticsamplethr$rs_interactionb   / sd(analyticsamplethr$rs_interactionb)
analyticsamplethr$isolation_index_std     <- analyticsamplethr$rs_isolationb     / sd(analyticsamplethr$rs_isolationb)

# IL-1beta tertiles
analyticsamplethr <- analyticsamplethr %>%
  mutate(IL_1beta_tertiles = ntile(log_il1, 3))


# =============================================================================
# Fit Linear Models — All Exposures x Biomarkers
# =============================================================================
exposures  <- c("dissimilarity_index_std", "interaction_index_std", "isolation_index_std")
biomarkers <- c("log_ddimer", "log_crp", "log_ifn", "log_tnf", "log_il6", "eselectin", "fix")

model_results <- list()

for (exposure in exposures) {
  for (biomarker in biomarkers) {
    key <- paste0(exposure, "_", biomarker)

    model_results[[key]] <- list(
      Unadjusted                           = summary(lm(as.formula(paste(biomarker, "~", exposure)), data = analyticsamplethr)),
      `Age and gender adjusted`            = summary(lm(as.formula(paste(biomarker, "~", exposure, "+ age + gender")), data = analyticsamplethr)),
      `age, gender, income`                = summary(lm(as.formula(paste(biomarker, "~", exposure, "+ age + gender + income")), data = analyticsamplethr)),
      `age, gender, education`             = summary(lm(as.formula(paste(biomarker, "~", exposure, "+ age + gender + ed_cat")), data = analyticsamplethr)),
      `age, gender, nSES`                  = summary(lm(as.formula(paste(biomarker, "~", exposure, "+ age + gender + nses_tract")), data = analyticsamplethr)),
      `age, gender, physical activity`     = summary(lm(as.formula(paste(biomarker, "~", exposure, "+ age + gender + exercise_cat")), data = analyticsamplethr)),
      `age, gender, diet`                  = summary(lm(as.formula(paste(biomarker, "~", exposure, "+ age + gender + diet7")), data = analyticsamplethr)),
      `age, gender, stroke/TIA`            = summary(lm(as.formula(paste(biomarker, "~", exposure, "+ age + gender + tia_sr")), data = analyticsamplethr)),
      `age, gender, coronary artery disease` = summary(lm(as.formula(paste(biomarker, "~", exposure, "+ age + gender + cad_sr_ecg")), data = analyticsamplethr)),
      `smoking indvidual`                  = summary(lm(as.formula(paste(biomarker, "~", exposure, "+ smoke + age + gender")), data = analyticsamplethr)),
      `alcohol consumption indvidual`      = summary(lm(as.formula(paste(biomarker, "~", exposure, "+ alc_niaaa + age + gender")), data = analyticsamplethr)),
      `diabetes  indvidual`                = summary(lm(as.formula(paste(biomarker, "~", exposure, "+ diab_srmed_glu + age + gender")), data = analyticsamplethr))
    )
  }
}


# =============================================================================
# Helper Functions — Excel Export
# =============================================================================

# Model names used in every export block
models_of_interest <- c(
  "Unadjusted", "Age and gender adjusted",
  "age, gender, income", "age, gender, education", "age, gender, nSES",
  "age, gender, physical activity", "age, gender, diet",
  "age, gender, stroke/TIA", "age, gender, coronary artery disease",
  "smoking indvidual", "alcohol consumption indvidual", "diabetes  indvidual"
)

# -----------------------------------------------------------------------------
# extract_lm_row()
# Pulls beta, SE, p-value, reverse-transformed % change, and 95% CI
# from one lm summary for a given exposure coefficient.
# -----------------------------------------------------------------------------
extract_lm_row <- function(model, exposure, model_name) {
  coef_mat  <- model$coefficients
  beta      <- coef_mat[exposure, "Estimate"]
  std_err   <- coef_mat[exposure, "Std. Error"]
  p_raw     <- coef_mat[exposure, "Pr(>|t|)"]
  p_fmt     <- ifelse(round(p_raw, 3) == 0.000, "<0.001", sprintf("%.3f", p_raw))
  rev_t     <- (exp(abs(beta)) - 1) * 100
  ci_lo     <- (exp(abs(beta) - 1.96 * std_err) - 1) * 100
  ci_hi     <- (exp(abs(beta) + 1.96 * std_err) - 1) * 100

  data.frame(
    Model         = model_name,
    Beta          = sprintf("%.3f", beta),
    Std_Error     = sprintf("%.3f", std_err),
    P_Value       = p_fmt,
    Reverse_T     = sprintf("%.2f%%", rev_t),
    ReverseT_95CI = sprintf("(%.2f%%, %.2f%%)", ci_lo, ci_hi),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# export_biomarker_wb()
# Builds a 3-sheet workbook (Dissimilarity / Isolation / Interaction)
# for one biomarker and saves it to output.path.
# -----------------------------------------------------------------------------
export_biomarker_wb <- function(biomarker_key, filename) {
  exposure_map <- list(
    Dissimilarity_Results = "dissimilarity_index_std",
    Isolation_Results     = "isolation_index_std",
    Interaction_Results   = "interaction_index_std"
  )

  wb <- createWorkbook()

  for (sheet_name in names(exposure_map)) {
    exposure   <- exposure_map[[sheet_name]]
    result_key <- paste0(exposure, "_", biomarker_key)
    models     <- model_results[[result_key]]

    rows <- lapply(models_of_interest, function(mn) {
      extract_lm_row(models[[mn]], exposure, mn)
    })
    df <- do.call(rbind, rows)

    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, df, colNames = TRUE)
  }

  saveWorkbook(wb, paste0(output.path, filename), overwrite = TRUE)
  message("Saved: ", filename)
}


# =============================================================================
# Export Results — One Workbook per Biomarker
# =============================================================================

# FIX (Factor IX — untransformed, linear CI)
export_biomarker_wb("fix",        "fix_results_regression_index.xlsx")

# D-Dimer (log-transformed)
export_biomarker_wb("log_ddimer", "ddimer_results_regression_index.xlsx")

# C-Reactive Protein (log-transformed)
export_biomarker_wb("log_crp",    "crp_results_regression_index.xlsx")

# Interferon-gamma (log-transformed)
export_biomarker_wb("log_ifn",    "ifn_results_regression_index.xlsx")

# TNF-alpha (log-transformed)
export_biomarker_wb("log_tnf",    "tnf_results_regression_index.xlsx")

# IL-6 (log-transformed)
export_biomarker_wb("log_il6",    "il6_results_regression_index.xlsx")

# E-Selectin (untransformed, linear CI)
export_biomarker_wb("eselectin",  "eselectin_results_regression_index.xlsx")
