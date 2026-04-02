# =============================================================================
# Residential Segregation & Thrombo-Inflammatory Biomarkers
# Regression Tables
# =============================================================================
# Description: Survey-weighted linear and logistic regression models examining
#              associations between residential segregation indices and
#              thrombo-inflammatory biomarkers in the REGARDS cohort.
#
# Exposures:   Dissimilarity index, Interaction index, Isolation index
#              (each standardized by SD and split into tertiles)
#
# Outcomes:    log_ddimer, log_crp, log_ifn, log_tnf, log_il6 (log-transformed),
#              eselectin, fix (untransformed),
#              IL-1b (tertile comparisons: T1 vs T2, T1 vs T3)
#
# Models:      Unadjusted + 11 sequentially adjusted models (age/gender,
#              income, education, neighborhood SES, physical activity, diet,
#              stroke/TIA, CAD, smoking, alcohol, diabetes)
# =============================================================================




###############
library(dplyr)
library(data.table)
library(zoo) 
library(tidyr)
library(lubridate) 
library(tidyverse)
library(knitr)
library(kableExtra)
library(sas7bdat)
library(haven)
library("nimble")
library(xlsx) 
library(openxlsx)
library(survey)
library(broom)
library(purrr)
#read in data..... IL-1b needs to be broken into quartiles (see analytic plan)



#analytic sample data




analyticsamplethr=readRDS('path/to/analyticsamplethr.rds')
head(analyticsamplethr)#4362 before exclusions....2534 after exclusions


analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr$rs_dissimilarityb), ]#none
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr$rs_interactionb), ]#none
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr$rs_isolationb), ]#none

analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["log_ddimer"]]), ]#148
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["log_crp"]]), ]#48
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["log_ifn"]]), ]#37
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["log_tnf"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["log_il6"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["eselectin"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["fix"]]), ]#6
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["il1"]]), ]#6


analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["age"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["gender"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.na(analyticsamplethr[["race"]]), ]#0


analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["log_ddimer"]]), ]#148
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["log_crp"]]), ]#48
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["log_ifn"]]), ]#37
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["log_tnf"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["log_il6"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["eselectin"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["fix"]]), ]#6
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["il1"]]), ]#6

analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["age"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["gender"]]), ]#0
analyticsamplethr <- analyticsamplethr[!is.nan(analyticsamplethr[["race"]]), ]#0

# Create new standardized variables  "rs_dissimilarityb", "rs_interactionb", "rs_isolationb"
analyticsamplethr$dissimilarity_index_std <- analyticsamplethr$rs_dissimilarityb/ sd(analyticsamplethr$rs_dissimilarityb)
analyticsamplethr$interaction_index_std <- analyticsamplethr$rs_interactionb / sd(analyticsamplethr$rs_interactionb)
analyticsamplethr$isolation_index_std <- analyticsamplethr$rs_isolationb / sd(analyticsamplethr$rs_isolationb)


#make tertitles
analyticsamplethr<-analyticsamplethr %>%mutate(dissimilarity_index_tertiles = ntile(rs_dissimilarityb, 3))
analyticsamplethr<-analyticsamplethr %>%mutate(interaction_index_tertiles = ntile(rs_interactionb, 3))
analyticsamplethr<-analyticsamplethr %>%mutate(isolation_index_tertiles = ntile(rs_isolationb, 3))




#TERTILES FOR THE INDINCES
analyticsamplethr$dissimilarity_index_tertiles <- factor(analyticsamplethr$dissimilarity_index_tertiles, levels = c("1", "2","3"))
analyticsamplethr$dissimilarity_index_tertiles <- relevel(analyticsamplethr$dissimilarity_index_tertiles, ref = "1")
analyticsamplethr$interaction_index_tertiles <- factor(analyticsamplethr$interaction_index_tertiles, levels = c("1", "2","3"))
analyticsamplethr$interaction_index_tertiles <- relevel(analyticsamplethr$interaction_index_tertiles, ref = "1")
analyticsamplethr$isolation_index_tertiles <- factor(analyticsamplethr$isolation_index_tertiles, levels = c("1", "2","3"))
analyticsamplethr$isolation_index_tertiles <- relevel(analyticsamplethr$isolation_index_tertiles, ref = "1")



############Need to run stratification for the groups first...
#then to work with the IL_1B biomarker we need to run 2 logisitc model 2v1 and 3v1 for the regression tbales

#####

# List of exposures
exposures <- c("dissimilarity_index_std", "interaction_index_std", "isolation_index_std")  # exposures are the 3 index cats
b <- c("log_ddimer", "log_crp", "log_ifn", "log_tnf", "log_il6", "eselectin","fix")#makesure fix isnt log
# List of models
models <- c("Unadjusted", "Age and gender adjusted", "age, gender, income", "age, gender, education",
            "age, gender, nSES", "age, gender, physical activity", "age, gender, diet",
            "age, gender, stroke/TIA", "age, gender, coronary artery disease")












# List to store model results
model_results <- list()

# Iterate over exposures
for (exposure in exposures) {
  # Iterate over biomarkers
  for (biomarker in b) {
    # Formulate the formula for each exposure and biomarker
    formula_unadjusted <- as.formula(paste(biomarker, "~", exposure))
    formula_age_gender <- as.formula(paste(biomarker, "~", exposure, "+ age + gender"))
    formula_income <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + income"))
    formula_education <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + ed_cat"))
    formula_nses <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + nses_tract"))
    formula_physical_activity <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + exercise_cat"))
    formula_diet <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + diet7"))
    formula_stroke_tia <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + tia_sr"))
    formula_cad <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + cad_sr_ecg"))
    formula_smokingindiv <- as.formula(paste(biomarker, "~", exposure, "+ smoke+ age + gender"))
    formula_alcindivd <- as.formula(paste(biomarker, "~", exposure, "+ alc_niaaa+ age + gender"))
    formula_dmindivid <- as.formula(paste(biomarker, "~", exposure, "+ diab_srmed_glu + age + gender"))
    
    
    #survey weights
    # Remove missing values in strata variables
    dgn= svydesign(~1, prob = ~prob, strata=~race*gender,data = analyticsamplethr)
    
    # Fit the models 
    model_unadjusted <- svyglm(formula_unadjusted,design=dgn, data = analyticsamplethr)
    model_age_gender <- svyglm(formula_age_gender,design=dgn, data = analyticsamplethr)
    model_income <- svyglm(formula_income,design=dgn, data = analyticsamplethr)
    model_education <- svyglm(formula_education,design=dgn, data = analyticsamplethr)
    model_nses <- svyglm(formula_nses,design=dgn, data = analyticsamplethr)
    model_physical_activity <- svyglm(formula_physical_activity,design=dgn, data = analyticsamplethr)
    model_diet <- svyglm(formula_diet,design=dgn, data = analyticsamplethr)
    model_stroke_tia <- svyglm(formula_stroke_tia,design=dgn, data = analyticsamplethr)
    model_cad <- svyglm(formula_cad,design=dgn, data = analyticsamplethr)
    model_smokingindiv <- svyglm(formula_smokingindiv,design=dgn, data = analyticsamplethr)
    model_alcindivd <- svyglm(formula_alcindivd,design=dgn, data = analyticsamplethr)
    model_dmindivid <- svyglm(formula_dmindivid,design=dgn, data = analyticsamplethr)
    
    # Store the model results
    model_results[[paste(exposure, "_", biomarker, sep = "")]] <- list(
      Unadjusted = summary(model_unadjusted),
      `Age and gender adjusted` = summary(model_age_gender),
      `age, gender, income` = summary(model_income),
      `age, gender, education` = summary(model_education),
      `age, gender, nSES` = summary(model_nses),
      `age, gender, physical activity` = summary(model_physical_activity),
      `age, gender, diet` = summary(model_diet),
      `age, gender, stroke/TIA` = summary(model_stroke_tia),
      `age, gender, coronary artery disease` = summary(model_cad),
      `smoking indvidual` = summary(model_smokingindiv),
      `alcohol consumption indvidual` = summary(model_alcindivd),
      `diabetes  indvidual` = summary(model_dmindivid)
      
      
    )
  }
}

#LOOP WORKS
#Need to do: standardize...Fix model results output table(kat c code)....log the biomarkers

# Calculate weighted standard deviation for blineFVIII



models_of_interest <- c("Unadjusted", "Age and gender adjusted", "age, gender, income", "age, gender, education", "age, gender, nSES",
                        "age, gender, physical activity", "age, gender, diet", "age, gender, stroke/TIA", "age, gender, coronary artery disease",
                        "smoking indvidual","alcohol consumption indvidual","diabetes  indvidual"
)



# Create a workbook
wb <- createWorkbook()

# Create the first sheet for dissimilarity results
addWorksheet(wb, sheetName = "Dissimilarity_Results")

# Create the second sheet for isolation results
addWorksheet(wb, sheetName = "Isolation_Results")

# Create the third sheet for interaction results
addWorksheet(wb, sheetName = "Interaction_Results")

# Initialize empty data frames for each set of results
results_df <- data.frame()
results_df_2 <- data.frame()
results_df_3 <- data.frame()

# Iterate over each model of interest for dissimilarity results
for (model_name in models_of_interest) {
  model <- model_results$dissimilarity_index_std_fix[[model_name]]
  
  beta_value=model$coefficients["dissimilarity_index_std", "Estimate"]
  std_err=model$coefficients["dissimilarity_index_std", "Std. Error"]
  
  p_value = model$coefficients["dissimilarity_index_std", "Pr(>|t|)"]
  p_value_formatted = ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  p_value = model$coefficients
  
  # Calculate the margin of error
  margin_of_error <- qnorm(0.975) * std_err
  
  # Calculate the lower and upper bounds of the CI
  lower_bound <- beta_value - margin_of_error
  upper_bound <- beta_value + margin_of_error
  
  # Extract the relevant information
  result <- data.frame(
    Model = model_name,
    SD=sd(analyticsamplethr$fix),
    Beta = sprintf("%.3f",beta_value),
    Std_Error = sprintf("%.3f",std_err) ,
    P_Value = p_value_formatted,
    CI_95 = sprintf("(%.3f, %.3f)", lower_bound, upper_bound),
    
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df data frame
  results_df <- rbind(results_df, result)
}
colnames(results_df) <- c("Model","sd biomarker", "Beta", "Std_Error", "P_Value", "95CI")

# Write dissimilarity results to the first sheet
writeData(wb, sheet = "Dissimilarity_Results", results_df, colNames = TRUE)





# Iterate over each model of interest for isolation results
for (model_name in models_of_interest) {
  model <- model_results$isolation_index_std_fix[[model_name]]
  
  beta_value=model$coefficients["isolation_index_std", "Estimate"]
  std_err=model$coefficients["isolation_index_std", "Std. Error"]
  
  p_value = model$coefficients["isolation_index_std", "Pr(>|t|)"]
  p_value_formatted = ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate the margin of error
  margin_of_error <- qnorm(0.975) * std_err
  
  # Calculate the lower and upper bounds of the CI
  lower_bound <- beta_value - margin_of_error
  upper_bound <- beta_value + margin_of_error
  
  # Extract the relevant information
  result <- data.frame(
    Model = model_name,
    SD=sd(analyticsamplethr$fix),
    Beta = sprintf("%.3f",beta_value),
    Std_Error = sprintf("%.3f",std_err) ,
    P_Value = p_value_formatted,
    CI_95 = sprintf("(%.3f, %.3f)", lower_bound, upper_bound),
    
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_2 data frame
  results_df_2 <- rbind(results_df_2, result)
}
colnames(results_df_2) <- c("Model","sd biomarker", "Beta", "Std_Error", "P_Value", "95CI")

# Write isolation results to the second sheet
writeData(wb, sheet = "Isolation_Results", results_df_2,colNames = TRUE)

# Iterate over each model of interest for interaction results
for (model_name in models_of_interest) {
  model <- model_results$interaction_index_std_fix[[model_name]]
  
  beta_value=model$coefficients["interaction_index_std", "Estimate"]
  std_err=model$coefficients["interaction_index_std", "Std. Error"]
  
  p_value = model$coefficients["interaction_index_std", "Pr(>|t|)"]
  p_value_formatted = ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate the margin of error
  margin_of_error <- qnorm(0.975) * std_err
  
  # Calculate the lower and upper bounds of the CI
  lower_bound <- beta_value - margin_of_error
  upper_bound <- beta_value + margin_of_error
  
  # Extract the relevant information
  result <- data.frame(
    Model = model_name,
    SD=sd(analyticsamplethr$fix),
    Beta = sprintf("%.3f",beta_value),
    Std_Error = sprintf("%.3f",std_err) ,
    P_Value = p_value_formatted,
    CI_95 = sprintf("(%.3f, %.3f)", lower_bound, upper_bound),
    
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_3 data frame
  results_df_3 <- rbind(results_df_3, result)
}
colnames(results_df_3) <- c("Model","sd biomarker", "Beta", "Std_Error", "P_Value", "95CI")


# Write interaction results to the third sheet
writeData(wb, sheet = "Interaction_Results", results_df_3,colNames = TRUE)

# Save the Excel workbook
saveWorkbook(wb, "path/to/output/fixresultsregressindexweighted.xlsx")









#DDIMER

library(openxlsx)

models_of_interest <- c("Unadjusted", "Age and gender adjusted", "age, gender, income", "age, gender, education", "age, gender, nSES",
                        "age, gender, physical activity", "age, gender, diet", "age, gender, stroke/TIA", "age, gender, coronary artery disease",
                        "smoking indvidual","alcohol consumption indvidual","diabetes  indvidual"
)

# Create a new Excel workbook
wb <- createWorkbook()

# Create the first sheet for dissimilarity results
addWorksheet(wb, sheetName = "Dissimilarity_Results")

# Create the second sheet for isolation results
addWorksheet(wb, sheetName = "Isolation_Results")

addWorksheet(wb, sheetName = "Interaction_Results")

# Initialize empty data frames for each set of results
results_df <- data.frame()
results_df_2 <- data.frame()
results_df_3 <- data.frame()

# Iterate over each model of interest for dissimilarity results
for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$dissimilarity_index_std_log_ddimer[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["dissimilarity_index_std", "Estimate"]
  std_err <- model$coefficients["dissimilarity_index_std", "Std. Error"]
  p_value <- model$coefficients["dissimilarity_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  # Ensure lcl_p is the smaller value and ucl_p is the larger value
  if (lcl_p > ucl_p) {
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df dataframe
  results_df <- rbind(results_df, result)
}

# Write dissimilarity results to the second sheet
writeData(wb, sheet = "Dissimilarity_Results", results_df, colNames = TRUE)

# Print or further process results_df as needed
print(results_df)



for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$isolation_index_std_log_ddimer[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["isolation_index_std", "Estimate"]
  std_err <- model$coefficients["isolation_index_std", "Std. Error"]
  p_value <- model$coefficients["isolation_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  if (lcl_p > ucl_p) {
    # Swap values if lcl_p is greater than ucl_p
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_2 dataframe
  results_df_2 <- rbind(results_df_2, result)
}

# Write isolation results to the second sheet
writeData(wb, sheet = "Isolation_Results", results_df_2, colNames = TRUE)



# Iterate over each model of interest for interaction results
for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$interaction_index_std_log_ddimer[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["interaction_index_std", "Estimate"]
  std_err <- model$coefficients["interaction_index_std", "Std. Error"]
  p_value <- model$coefficients["interaction_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  # Ensure lcl_p is the smaller value and ucl_p is the larger value
  if (lcl_p > ucl_p) {
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_3 dataframe
  results_df_3 <- rbind(results_df_3, result)
}

# Write interaction results to the second sheet
writeData(wb, sheet = "Interaction_Results", results_df_3, colNames = TRUE)

# Print or further process results_df_3 as needed
print(results_df_3)

# Save the Excel workbook
saveWorkbook(wb, "path/to/output/ddimerresultsregressindexweighted.xlsx")




#crp

library(openxlsx)

models_of_interest <- c("Unadjusted", "Age and gender adjusted", "age, gender, income", "age, gender, education", "age, gender, nSES",
                        "age, gender, physical activity", "age, gender, diet", "age, gender, stroke/TIA", "age, gender, coronary artery disease",
                        "smoking indvidual","alcohol consumption indvidual","diabetes  indvidual"
)

# Create a new Excel workbook
wb <- createWorkbook()

# Create the first sheet for dissimilarity results
addWorksheet(wb, sheetName = "Dissimilarity_Results")

# Create the second sheet for isolation results
addWorksheet(wb, sheetName = "Isolation_Results")

addWorksheet(wb, sheetName = "Interaction_Results")

# Initialize empty data frames for each set of results
results_df <- data.frame()
results_df_2 <- data.frame()
results_df_3 <- data.frame()

# Iterate over each model of interest for dissimilarity results
for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$dissimilarity_index_std_log_crp[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["dissimilarity_index_std", "Estimate"]
  std_err <- model$coefficients["dissimilarity_index_std", "Std. Error"]
  p_value <- model$coefficients["dissimilarity_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  # Ensure lcl_p is the smaller value and ucl_p is the larger value
  if (lcl_p > ucl_p) {
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df dataframe
  results_df <- rbind(results_df, result)
}

# Write dissimilarity results to the second sheet
writeData(wb, sheet = "Dissimilarity_Results", results_df, colNames = TRUE)

# Print or further process results_df as needed
print(results_df)



for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$isolation_index_std_log_crp[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["isolation_index_std", "Estimate"]
  std_err <- model$coefficients["isolation_index_std", "Std. Error"]
  p_value <- model$coefficients["isolation_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  if (lcl_p > ucl_p) {
    # Swap values if lcl_p is greater than ucl_p
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_2 dataframe
  results_df_2 <- rbind(results_df_2, result)
}

# Write isolation results to the second sheet
writeData(wb, sheet = "Isolation_Results", results_df_2, colNames = TRUE)



# Iterate over each model of interest for interaction results
for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$interaction_index_std_log_crp[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["interaction_index_std", "Estimate"]
  std_err <- model$coefficients["interaction_index_std", "Std. Error"]
  p_value <- model$coefficients["interaction_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  # Ensure lcl_p is the smaller value and ucl_p is the larger value
  if (lcl_p > ucl_p) {
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_3 dataframe
  results_df_3 <- rbind(results_df_3, result)
}

# Write interaction results to the second sheet
writeData(wb, sheet = "Interaction_Results", results_df_3, colNames = TRUE)

# Print or further process results_df_3 as needed
print(results_df_3)

# Save the Excel workbook
saveWorkbook(wb, "path/to/output/crpresultsregressindexweighted.xlsx")



#ifn

library(openxlsx)

models_of_interest <- c("Unadjusted", "Age and gender adjusted", "age, gender, income", "age, gender, education", "age, gender, nSES",
                        "age, gender, physical activity", "age, gender, diet", "age, gender, stroke/TIA", "age, gender, coronary artery disease",
                        "smoking indvidual","alcohol consumption indvidual","diabetes  indvidual"
)

# Create a new Excel workbook
wb <- createWorkbook()

# Create the first sheet for dissimilarity results
addWorksheet(wb, sheetName = "Dissimilarity_Results")

# Create the second sheet for isolation results
addWorksheet(wb, sheetName = "Isolation_Results")

addWorksheet(wb, sheetName = "Interaction_Results")

# Initialize empty data frames for each set of results
results_df <- data.frame()
results_df_2 <- data.frame()
results_df_3 <- data.frame()

# Iterate over each model of interest for dissimilarity results
for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$dissimilarity_index_std_log_ifn[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["dissimilarity_index_std", "Estimate"]
  std_err <- model$coefficients["dissimilarity_index_std", "Std. Error"]
  p_value <- model$coefficients["dissimilarity_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  # Ensure lcl_p is the smaller value and ucl_p is the larger value
  if (lcl_p > ucl_p) {
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df dataframe
  results_df <- rbind(results_df, result)
}

# Write dissimilarity results to the second sheet
writeData(wb, sheet = "Dissimilarity_Results", results_df, colNames = TRUE)

# Print or further process results_df as needed
print(results_df)



for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$isolation_index_std_log_ifn[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["isolation_index_std", "Estimate"]
  std_err <- model$coefficients["isolation_index_std", "Std. Error"]
  p_value <- model$coefficients["isolation_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  if (lcl_p > ucl_p) {
    # Swap values if lcl_p is greater than ucl_p
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_2 dataframe
  results_df_2 <- rbind(results_df_2, result)
}

# Write isolation results to the second sheet
writeData(wb, sheet = "Isolation_Results", results_df_2, colNames = TRUE)



# Iterate over each model of interest for interaction results
for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$interaction_index_std_log_ifn[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["interaction_index_std", "Estimate"]
  std_err <- model$coefficients["interaction_index_std", "Std. Error"]
  p_value <- model$coefficients["interaction_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  # Ensure lcl_p is the smaller value and ucl_p is the larger value
  if (lcl_p > ucl_p) {
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_3 dataframe
  results_df_3 <- rbind(results_df_3, result)
}

# Write interaction results to the second sheet
writeData(wb, sheet = "Interaction_Results", results_df_3, colNames = TRUE)

# Print or further process results_df_3 as needed
print(results_df_3)

# Save the Excel workbook
saveWorkbook(wb, "path/to/output/ifnresultsregressindexweighted.xlsx")



#tnf

library(openxlsx)

models_of_interest <- c("Unadjusted", "Age and gender adjusted", "age, gender, income", "age, gender, education", "age, gender, nSES",
                        "age, gender, physical activity", "age, gender, diet", "age, gender, stroke/TIA", "age, gender, coronary artery disease",
                        "smoking indvidual","alcohol consumption indvidual","diabetes  indvidual"
)

# Create a new Excel workbook
wb <- createWorkbook()

# Create the first sheet for dissimilarity results
addWorksheet(wb, sheetName = "Dissimilarity_Results")

# Create the second sheet for isolation results
addWorksheet(wb, sheetName = "Isolation_Results")

addWorksheet(wb, sheetName = "Interaction_Results")

# Initialize empty data frames for each set of results
results_df <- data.frame()
results_df_2 <- data.frame()
results_df_3 <- data.frame()

# Iterate over each model of interest for dissimilarity results
for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$dissimilarity_index_std_log_tnf[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["dissimilarity_index_std", "Estimate"]
  std_err <- model$coefficients["dissimilarity_index_std", "Std. Error"]
  p_value <- model$coefficients["dissimilarity_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  # Ensure lcl_p is the smaller value and ucl_p is the larger value
  if (lcl_p > ucl_p) {
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df dataframe
  results_df <- rbind(results_df, result)
}

# Write dissimilarity results to the second sheet
writeData(wb, sheet = "Dissimilarity_Results", results_df, colNames = TRUE)

# Print or further process results_df as needed
print(results_df)



for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$isolation_index_std_log_tnf[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["isolation_index_std", "Estimate"]
  std_err <- model$coefficients["isolation_index_std", "Std. Error"]
  p_value <- model$coefficients["isolation_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  if (lcl_p > ucl_p) {
    # Swap values if lcl_p is greater than ucl_p
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_2 dataframe
  results_df_2 <- rbind(results_df_2, result)
}

# Write isolation results to the second sheet
writeData(wb, sheet = "Isolation_Results", results_df_2, colNames = TRUE)



# Iterate over each model of interest for interaction results
for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$interaction_index_std_log_tnf[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["interaction_index_std", "Estimate"]
  std_err <- model$coefficients["interaction_index_std", "Std. Error"]
  p_value <- model$coefficients["interaction_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  # Ensure lcl_p is the smaller value and ucl_p is the larger value
  if (lcl_p > ucl_p) {
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_3 dataframe
  results_df_3 <- rbind(results_df_3, result)
}

# Write interaction results to the second sheet
writeData(wb, sheet = "Interaction_Results", results_df_3, colNames = TRUE)

# Print or further process results_df_3 as needed
print(results_df_3)

# Save the Excel workbook
saveWorkbook(wb, "path/to/output/tnfresultsregressindexweighted.xlsx")



      
#il6

library(openxlsx)

models_of_interest <- c("Unadjusted", "Age and gender adjusted", "age, gender, income", "age, gender, education", "age, gender, nSES",
                        "age, gender, physical activity", "age, gender, diet", "age, gender, stroke/TIA", "age, gender, coronary artery disease",
                        "smoking indvidual","alcohol consumption indvidual","diabetes  indvidual"
)

# Create a new Excel workbook
wb <- createWorkbook()

# Create the first sheet for dissimilarity results
addWorksheet(wb, sheetName = "Dissimilarity_Results")

# Create the second sheet for isolation results
addWorksheet(wb, sheetName = "Isolation_Results")

addWorksheet(wb, sheetName = "Interaction_Results")

# Initialize empty data frames for each set of results
results_df <- data.frame()
results_df_2 <- data.frame()
results_df_3 <- data.frame()

# Iterate over each model of interest for dissimilarity results
for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$dissimilarity_index_std_log_il6[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["dissimilarity_index_std", "Estimate"]
  std_err <- model$coefficients["dissimilarity_index_std", "Std. Error"]
  p_value <- model$coefficients["dissimilarity_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  # Ensure lcl_p is the smaller value and ucl_p is the larger value
  if (lcl_p > ucl_p) {
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df dataframe
  results_df <- rbind(results_df, result)
}

# Write dissimilarity results to the second sheet
writeData(wb, sheet = "Dissimilarity_Results", results_df, colNames = TRUE)

# Print or further process results_df as needed
print(results_df)



for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$isolation_index_std_log_il6[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["isolation_index_std", "Estimate"]
  std_err <- model$coefficients["isolation_index_std", "Std. Error"]
  p_value <- model$coefficients["isolation_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  if (lcl_p > ucl_p) {
    # Swap values if lcl_p is greater than ucl_p
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_2 dataframe
  results_df_2 <- rbind(results_df_2, result)
}

# Write isolation results to the second sheet
writeData(wb, sheet = "Isolation_Results", results_df_2, colNames = TRUE)



# Iterate over each model of interest for interaction results
for (model_name in models_of_interest) {
  # Extract the model from model_results using the model_name
  model <- model_results$interaction_index_std_log_il6[[model_name]]
  
  # Extract coefficients and statistics of interest from the model
  beta_value <- model$coefficients["interaction_index_std", "Estimate"]
  std_err <- model$coefficients["interaction_index_std", "Std. Error"]
  p_value <- model$coefficients["interaction_index_std", "Pr(>|t|)"]
  p_value_formatted <- ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate sign of beta
  sign_beta <- ifelse(beta_value >= 0, 1, -1)
  
  # Calculate reverse transformation percentage
  reverse_t <- (exp(abs(beta_value)) - 1) * 100
  reverse_t_final <- sign_beta * reverse_t  # Multiply reverse_t by sign_beta
  
  # Calculate confidence interval (CI) and related statistics
  ci_stat <- 1.96 * std_err
  b_ci <- abs(beta_value)
  lcits <- b_ci - ci_stat
  ucits <- b_ci + ci_stat
  
  lcit <- (exp(lcits) - 1) * 100
  ucit <- (exp(ucits) - 1) * 100
  
  lcl_p <- sign_beta * lcit
  ucl_p <- sign_beta * ucit
  
  # Ensure lcl_p is the smaller value and ucl_p is the larger value
  if (lcl_p > ucl_p) {
    temp <- lcl_p
    lcl_p <- ucl_p
    ucl_p <- temp
  }
  
  # Format CI as a string
  ci <- sprintf("(%0.2f%%, %0.2f%%)", lcl_p, ucl_p)
  
  # Extract relevant information and create a result dataframe
  result <- data.frame(
    Model = model_name,
    sd = sd(analyticsamplethr$tnf),  # Assuming 'analyticsamplethr$tnf' is your data for sd calculation
    Beta = sprintf("%.3f", beta_value),
    Std_Error = sprintf("%.3f", std_err),
    P_Value = p_value_formatted,
    Reverse_T = sprintf("%.2f%%", reverse_t_final),  # Format reverse_t_final as percentage
    CI = ci,
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_3 dataframe
  results_df_3 <- rbind(results_df_3, result)
}

# Write interaction results to the second sheet
writeData(wb, sheet = "Interaction_Results", results_df_3, colNames = TRUE)

# Print or further process results_df_3 as needed
print(results_df_3)

# Save the Excel workbook
saveWorkbook(wb, "path/to/output/il6resultsregressindexweighted.xlsx")




#eselectin.............Reverse trans form this is wonky because not logged

library(openxlsx)

models_of_interest <- c("Unadjusted", "Age and gender adjusted", "age, gender, income", "age, gender, education", "age, gender, nSES",
                        "age, gender, physical activity", "age, gender, diet", "age, gender, stroke/TIA", "age, gender, coronary artery disease",
                        "smoking indvidual","alcohol consumption indvidual","diabetes  indvidual"
)

# Create a new Excel workbook
wb <- createWorkbook()

# Create the first sheet for dissimilarity results
addWorksheet(wb, sheetName = "Dissimilarity_Results")

# Create the second sheet for isolation results
addWorksheet(wb, sheetName = "Isolation_Results")

addWorksheet(wb, sheetName = "Interaction_Results")

# Initialize empty data frames for each set of results
results_df <- data.frame()
results_df_2 <- data.frame()
results_df_3 <- data.frame()

# Iterate over each model of interest for dissimilarity results
for (model_name in models_of_interest) {
  model <- model_results$dissimilarity_index_std_eselectin[[model_name]]
  
  beta_value=model$coefficients["dissimilarity_index_std", "Estimate"]
  std_err=model$coefficients["dissimilarity_index_std", "Std. Error"]
  
  p_value = model$coefficients["dissimilarity_index_std", "Pr(>|t|)"]
  p_value_formatted = ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate the margin of error
  margin_of_error <- qnorm(0.975) * std_err
  
  # Calculate the lower and upper bounds of the CI
  lower_bound <- beta_value - margin_of_error
  upper_bound <- beta_value + margin_of_error
  
  # Extract the relevant information
  result <- data.frame(
    Model = model_name,
    sd=sd(analyticsamplethr$eselectin),
    Beta = sprintf("%.3f",beta_value),
    Std_Error = sprintf("%.3f",std_err) ,
    P_Value = p_value_formatted,
    CI_95 = sprintf("(%.3f, %.3f)", lower_bound, upper_bound),
    
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df data frame
  results_df <- rbind(results_df, result)
}
colnames(results_df) <- c("Model","sd biomarker", "Beta", "Std_Error", "P_Value", "95CI")

# Write dissimilarity results to the first sheet
writeData(wb, sheet = "Dissimilarity_Results", results_df, colNames = TRUE)





# Iterate over each model of interest for isolation results
for (model_name in models_of_interest) {
  model <- model_results$isolation_index_std_eselectin[[model_name]]
  
  beta_value=model$coefficients["isolation_index_std", "Estimate"]
  std_err=model$coefficients["isolation_index_std", "Std. Error"]
  
  p_value = model$coefficients["isolation_index_std", "Pr(>|t|)"]
  p_value_formatted = ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate the margin of error
  margin_of_error <- qnorm(0.975) * std_err
  
  # Calculate the lower and upper bounds of the CI
  lower_bound <- beta_value - margin_of_error
  upper_bound <- beta_value + margin_of_error
  
  # Extract the relevant information
  result <- data.frame(
    Model = model_name,
    sd=sd(analyticsamplethr$eselectin),
    Beta = sprintf("%.3f",beta_value),
    Std_Error = sprintf("%.3f",std_err) ,
    P_Value = p_value_formatted,
    CI_95 = sprintf("(%.3f, %.3f)", lower_bound, upper_bound),
    
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_2 data frame
  results_df_2 <- rbind(results_df_2, result)
}
colnames(results_df_2) <- c("Model","sd biomarker", "Beta", "Std_Error", "P_Value", "95CI")

# Write isolation results to the second sheet
writeData(wb, sheet = "Isolation_Results", results_df_2,colNames = TRUE)

# Iterate over each model of interest for interaction results
for (model_name in models_of_interest) {
  model <- model_results$interaction_index_std_eselectin[[model_name]]
  
  beta_value=model$coefficients["interaction_index_std", "Estimate"]
  std_err=model$coefficients["interaction_index_std", "Std. Error"]
  
  p_value = model$coefficients["interaction_index_std", "Pr(>|t|)"]
  p_value_formatted = ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  # Calculate the margin of error
  margin_of_error <- qnorm(0.975) * std_err
  
  # Calculate the lower and upper bounds of the CI
  lower_bound <- beta_value - margin_of_error
  upper_bound <- beta_value + margin_of_error
  
  # Extract the relevant information
  result <- data.frame(
    Model = model_name,
    sd=sd(analyticsamplethr$eselectin),
    Beta = sprintf("%.3f",beta_value),
    Std_Error = sprintf("%.3f",std_err) ,
    P_Value = p_value_formatted,
    CI_95 = sprintf("(%.3f, %.3f)", lower_bound, upper_bound),
    
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_3 data frame
  results_df_3 <- rbind(results_df_3, result)
}
colnames(results_df_3) <- c("Model","sd biomarker", "Beta", "Std_Error", "P_Value", "95CI")


# Write interaction results to the third sheet
writeData(wb, sheet = "Interaction_Results", results_df_3,colNames = TRUE)

# Save the Excel workbook
saveWorkbook(wb, "path/to/output/eselectinresultsregressindexweighted.xlsx")




# 
# 
# #il-1b for this bio marker re run below



#for il1 if less than 0.01 replace with .009.....create this in a new var il1b2 us the new var for analysis
#to do ths section of anlysis.. make il1b2 into tertiles....then make 2 subsets of the data 1 that will have tertile 1 and tertile 2, and
# another that will have tertile 1 and tertile 3
analyticsamplethr$il1b2 <- ifelse(analyticsamplethr$il1 < 0.01, 0.009, analyticsamplethr$il1)


# Step 1: Create Tertiles for il1b2
analyticsamplethr$il1b2_tertiles <- cut(analyticsamplethr$il1b2, breaks = quantile(analyticsamplethr$il1b2, probs = c(0, 1/3, 2/3, 1)), labels = c("Tertile 1", "Tertile 2", "Tertile 3"))

# Step 2: Subset the data into two groups based on tertiles
# Subset 1: Tertile 1 and Tertile 2
il1b1v2dat <- subset(analyticsamplethr, il1b2_tertiles %in% c("Tertile 1", "Tertile 2"))#n=2889
il1b1v2dat$highest_tertile=ifelse(il1b1v2dat$il1b2_tertiles=="Tertile 2",1,0)


# Subset 2: Tertile 1 and Tertile 3
il1b1v3dat <- subset(analyticsamplethr, il1b2_tertiles %in% c("Tertile 1", "Tertile 3"))#n=2607
il1b1v3dat$highest_tertile=ifelse(il1b1v3dat$il1b2_tertiles=="Tertile 3",1,0)

##


############T1 Vs. T2

# List of exposures
exposures <- c("dissimilarity_index_std", "interaction_index_std", "isolation_index_std")  # exposures are the 3 index cats
b <- c("highest_tertile")#makesure fix isnt log
# List of models
models <- c("Unadjusted", "Age and gender adjusted", "age, gender, income", "age, gender, education",
            "age, gender, nSES", "age, gender, physical activity", "age, gender, diet",
            "age, gender, stroke/TIA", "age, gender, coronary artery disease")












# List to store model results
model_results_il1b1v2dat <- list()

# Iterate over exposures
for (exposure in exposures) {
  # Iterate over biomarkers
  for (biomarker in b) {
    # Formulate the formula for each exposure and biomarker
    formula_unadjusted <- as.formula(paste(biomarker, "~", exposure))
    formula_age_gender <- as.formula(paste(biomarker, "~", exposure, "+ age + gender"))
    formula_income <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + income"))
    formula_education <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + ed_cat"))
    formula_nses <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + nses_tract"))
    formula_physical_activity <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + exercise_cat"))
    formula_diet <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + diet7"))
    formula_stroke_tia <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + tia_sr"))
    formula_cad <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + cad_sr_ecg"))
    formula_smokingindiv <- as.formula(paste(biomarker, "~", exposure, "+ smoke+ age + gender"))
    formula_alcindivd <- as.formula(paste(biomarker, "~", exposure, "+ alc_niaaa+ age + gender"))
    formula_dmindivid <- as.formula(paste(biomarker, "~", exposure, "+ diab_srmed_glu + age + gender"))
    
    
    #survey weights
    # Remove missing values in strata variables
    dgn= svydesign(~1, prob = ~prob, strata=~race*gender,data = il1b1v2dat)
    
    # Fit the models 
    model_unadjusted <- svyglm(formula_unadjusted,design=dgn,family=binomial,data = il1b1v2dat)
    model_age_gender <- svyglm(formula_age_gender,design=dgn,family=binomial, data = il1b1v2dat)
    model_income <- svyglm(formula_income,design=dgn,family=binomial, data = il1b1v2dat)
    model_education <- svyglm(formula_education,design=dgn,family=binomial, data = il1b1v2dat)
    model_nses <- svyglm(formula_nses,design=dgn, family=binomial,data = il1b1v2dat)
    model_physical_activity <- svyglm(formula_physical_activity,design=dgn,family=binomial, data = il1b1v2dat)
    model_diet <- svyglm(formula_diet,design=dgn,family=binomial, data = il1b1v2dat)
    model_stroke_tia <- svyglm(formula_stroke_tia,design=dgn,family=binomial, data = il1b1v2dat)
    model_cad <- svyglm(formula_cad,design=dgn,family=binomial, data = il1b1v2dat)
    model_smokingindiv <- svyglm(formula_smokingindiv,design=dgn,family=binomial, data = il1b1v2dat)
    model_alcindivd <- svyglm(formula_alcindivd,design=dgn,family=binomial, data = il1b1v2dat)
    model_dmindivid <- svyglm(formula_dmindivid,design=dgn,family=binomial, data = il1b1v2dat)
    
    # Store the model results
    model_results_il1b1v2dat[[paste(exposure, "_", biomarker, sep = "")]] <- list(
      Unadjusted = summary(model_unadjusted),
      `Age and gender adjusted` = summary(model_age_gender),
      `age, gender, income` = summary(model_income),
      `age, gender, education` = summary(model_education),
      `age, gender, nSES` = summary(model_nses),
      `age, gender, physical activity` = summary(model_physical_activity),
      `age, gender, diet` = summary(model_diet),
      `age, gender, stroke/TIA` = summary(model_stroke_tia),
      `age, gender, coronary artery disease` = summary(model_cad),
      `smoking indvidual` = summary(model_smokingindiv),
      `alcohol consumption indvidual` = summary(model_alcindivd),
      `diabetes  indvidual` = summary(model_dmindivid)
      
      
      
    )
  }
}

library(openxlsx)

models_of_interest <- c("Unadjusted", "Age and gender adjusted", "age, gender, income", "age, gender, education", "age, gender, nSES",
                        "age, gender, physical activity", "age, gender, diet", "age, gender, stroke/TIA", "age, gender, coronary artery disease",
                        "smoking indvidual","alcohol consumption indvidual","diabetes  indvidual"
)

# Create a new Excel workbook
wb <- createWorkbook()

# Create the first sheet for dissimilarity results
addWorksheet(wb, sheetName = "Dissimilarity_Results")

# Create the second sheet for isolation results
addWorksheet(wb, sheetName = "Isolation_Results")

addWorksheet(wb, sheetName = "Interaction_Results")

# Initialize empty data frames for each set of results
results_df <- data.frame()
results_df_2 <- data.frame()
results_df_3 <- data.frame()

# Iterate over each model of interest for dissimilarity results
for (model_name in models_of_interest) {
  model <- model_results_il1b1v2dat$dissimilarity_index_std_highest_tertile[[model_name]]
  
  beta_value=model$coefficients["dissimilarity_index_std", "Estimate"]
  std_err=model$coefficients["dissimilarity_index_std", "Std. Error"]
  
  p_value = model$coefficients["dissimilarity_index_std", "Pr(>|t|)"]
  p_value_formatted = ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  odds_ratio <- exp(beta_value)
  
  # Calculate margin of error for OR
  or_margin_of_error <- qnorm(0.975) * std_err
  
  # Calculate lower and upper bounds of the CI for OR
  or_ci_lower <- exp(beta_value - or_margin_of_error)
  or_ci_upper <- exp(beta_value + or_margin_of_error)
  
  # Extract the relevant information
  result <- data.frame(
    Model = model_name,
    sd=sd(analyticsamplethr$eselectin),
    Beta = sprintf("%.3f",beta_value),
    Odds_Ratio = sprintf("%.3f", odds_ratio),
    Std_Error = sprintf("%.3f",std_err) ,
    P_Value = p_value_formatted,
    OR_CI_95 = sprintf("(%.3f, %.3f)", or_ci_lower, or_ci_upper),    
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df data frame
  results_df <- rbind(results_df, result)
}


# Write dissimilarity results to the first sheet
writeData(wb, sheet = "Dissimilarity_Results", results_df, colNames = TRUE)





# Iterate over each model of interest for isolation results
for (model_name in models_of_interest) {
  model <- model_results_il1b1v2dat$isolation_index_std_highest_tertile[[model_name]]
  
  beta_value=model$coefficients["isolation_index_std", "Estimate"]
  std_err=model$coefficients["isolation_index_std", "Std. Error"]
  
  p_value = model$coefficients["isolation_index_std", "Pr(>|t|)"]
  p_value_formatted = ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  odds_ratio <- exp(beta_value)
  
  # Calculate margin of error for OR
  or_margin_of_error <- qnorm(0.975) * std_err
  
  # Calculate lower and upper bounds of the CI for OR
  or_ci_lower <- exp(beta_value - or_margin_of_error)
  or_ci_upper <- exp(beta_value + or_margin_of_error)
  
  # Extract the relevant information
  result <- data.frame(
    Model = model_name,
    sd=sd(analyticsamplethr$eselectin),
    Beta = sprintf("%.3f",beta_value),
    Odds_Ratio = sprintf("%.3f", odds_ratio),
    Std_Error = sprintf("%.3f",std_err) ,
    P_Value = p_value_formatted,
    OR_CI_95 = sprintf("(%.3f, %.3f)", or_ci_lower, or_ci_upper),    
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_2 data frame
  results_df_2 <- rbind(results_df_2, result)
}


# Write isolation results to the second sheet
writeData(wb, sheet = "Isolation_Results", results_df_2,colNames = TRUE)

# Iterate over each model of interest for interaction results
for (model_name in models_of_interest) {
  model <- model_results_il1b1v2dat$interaction_index_std_highest_tertile[[model_name]]
  
  beta_value=model$coefficients["interaction_index_std", "Estimate"]
  std_err=model$coefficients["interaction_index_std", "Std. Error"]
  
  p_value = model$coefficients["interaction_index_std", "Pr(>|t|)"]
  p_value_formatted = ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  odds_ratio <- exp(beta_value)
  
  # Calculate margin of error for OR
  or_margin_of_error <- qnorm(0.975) * std_err
  
  # Calculate lower and upper bounds of the CI for OR
  or_ci_lower <- exp(beta_value - or_margin_of_error)
  or_ci_upper <- exp(beta_value + or_margin_of_error)
  
  # Extract the relevant information
  result <- data.frame(
    Model = model_name,
    sd=sd(analyticsamplethr$eselectin),
    Beta = sprintf("%.3f",beta_value),
    Odds_Ratio = sprintf("%.3f", odds_ratio),
    Std_Error = sprintf("%.3f",std_err) ,
    P_Value = p_value_formatted,
    OR_CI_95 = sprintf("(%.3f, %.3f)", or_ci_lower, or_ci_upper),    
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_3 data frame
  results_df_3 <- rbind(results_df_3, result)
}



# Write interaction results to the third sheet
writeData(wb, sheet = "Interaction_Results", results_df_3,colNames = TRUE)

# Save the Excel workbook
saveWorkbook(wb, "path/to/output/il1bt1vt2resultsregressindexweighted.xlsx")






############T1 Vs. T3

# List of exposures
exposures <- c("dissimilarity_index_std", "interaction_index_std", "isolation_index_std")  # exposures are the 3 index cats
b <- c("highest_tertile")#makesure fix isnt log
# List of models
models <- c("Unadjusted", "Age and gender adjusted", "age, gender, income", "age, gender, education",
            "age, gender, nSES", "age, gender, physical activity", "age, gender, diet",
            "age, gender, stroke/TIA", "age, gender, coronary artery disease")












# List to store model results
model_results_il1b1v3dat <- list()

# Iterate over exposures
for (exposure in exposures) {
  # Iterate over biomarkers
  for (biomarker in b) {
    # Formulate the formula for each exposure and biomarker
    formula_unadjusted <- as.formula(paste(biomarker, "~", exposure))
    formula_age_gender <- as.formula(paste(biomarker, "~", exposure, "+ age + gender"))
    formula_income <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + income"))
    formula_education <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + ed_cat"))
    formula_nses <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + nses_tract"))
    formula_physical_activity <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + exercise_cat"))
    formula_diet <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + diet7"))
    formula_stroke_tia <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + tia_sr"))
    formula_cad <- as.formula(paste(biomarker, "~", exposure, "+ age + gender + cad_sr_ecg"))
    formula_smokingindiv <- as.formula(paste(biomarker, "~", exposure, "+ smoke+ age + gender"))
    formula_alcindivd <- as.formula(paste(biomarker, "~", exposure, "+ alc_niaaa+ age + gender"))
    formula_dmindivid <- as.formula(paste(biomarker, "~", exposure, "+ diab_srmed_glu + age + gender"))
    
    
    #survey weights
    # Remove missing values in strata variables
    dgn= svydesign(~1, prob = ~prob, strata=~race*gender,data = il1b1v3dat)
    
    # Fit the models 
    model_unadjusted <- svyglm(formula_unadjusted,design=dgn,family=binomial,data = il1b1v3dat)
    model_age_gender <- svyglm(formula_age_gender,design=dgn,family=binomial, data = il1b1v3dat)
    model_income <- svyglm(formula_income,design=dgn,family=binomial, data = il1b1v3dat)
    model_education <- svyglm(formula_education,design=dgn,family=binomial, data = il1b1v3dat)
    model_nses <- svyglm(formula_nses,design=dgn, family=binomial,data = il1b1v3dat)
    model_physical_activity <- svyglm(formula_physical_activity,design=dgn,family=binomial, data = il1b1v3dat)
    model_diet <- svyglm(formula_diet,design=dgn,family=binomial, data = il1b1v3dat)
    model_stroke_tia <- svyglm(formula_stroke_tia,design=dgn,family=binomial, data = il1b1v3dat)
    model_cad <- svyglm(formula_cad,design=dgn,family=binomial, data = il1b1v3dat)
    model_smokingindiv <- svyglm(formula_smokingindiv,design=dgn,family=binomial, data = il1b1v3dat)
    model_alcindivd <- svyglm(formula_alcindivd,design=dgn,family=binomial, data = il1b1v3dat)
    model_dmindivid <- svyglm(formula_dmindivid,design=dgn,family=binomial, data = il1b1v3dat)
    
    # Store the model results
    model_results_il1b1v3dat[[paste(exposure, "_", biomarker, sep = "")]] <- list(
      Unadjusted = summary(model_unadjusted),
      `Age and gender adjusted` = summary(model_age_gender),
      `age, gender, income` = summary(model_income),
      `age, gender, education` = summary(model_education),
      `age, gender, nSES` = summary(model_nses),
      `age, gender, physical activity` = summary(model_physical_activity),
      `age, gender, diet` = summary(model_diet),
      `age, gender, stroke/TIA` = summary(model_stroke_tia),
      `age, gender, coronary artery disease` = summary(model_cad),
      `smoking indvidual` = summary(model_smokingindiv),
      `alcohol consumption indvidual` = summary(model_alcindivd),
      `diabetes  indvidual` = summary(model_dmindivid)
      
      
      
    )
  }
}

library(openxlsx)

models_of_interest <- c("Unadjusted", "Age and gender adjusted", "age, gender, income", "age, gender, education", "age, gender, nSES",
                        "age, gender, physical activity", "age, gender, diet", "age, gender, stroke/TIA", "age, gender, coronary artery disease",
                        "smoking indvidual","alcohol consumption indvidual","diabetes  indvidual"
)

# Create a new Excel workbook
wb <- createWorkbook()

# Create the first sheet for dissimilarity results
addWorksheet(wb, sheetName = "Dissimilarity_Results")

# Create the second sheet for isolation results
addWorksheet(wb, sheetName = "Isolation_Results")

addWorksheet(wb, sheetName = "Interaction_Results")

# Initialize empty data frames for each set of results
results_df <- data.frame()
results_df_2 <- data.frame()
results_df_3 <- data.frame()

# Iterate over each model of interest for dissimilarity results
for (model_name in models_of_interest) {
  model <- model_results_il1b1v3dat$dissimilarity_index_std_highest_tertile[[model_name]]
  
  beta_value=model$coefficients["dissimilarity_index_std", "Estimate"]
  std_err=model$coefficients["dissimilarity_index_std", "Std. Error"]
  
  p_value = model$coefficients["dissimilarity_index_std", "Pr(>|t|)"]
  p_value_formatted = ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  odds_ratio <- exp(beta_value)
  
  # Calculate margin of error for OR
  or_margin_of_error <- qnorm(0.975) * std_err
  
  # Calculate lower and upper bounds of the CI for OR
  or_ci_lower <- exp(beta_value - or_margin_of_error)
  or_ci_upper <- exp(beta_value + or_margin_of_error)
  
  # Extract the relevant information
  result <- data.frame(
    Model = model_name,
    sd=sd(analyticsamplethr$eselectin),
    Beta = sprintf("%.3f",beta_value),
    Odds_Ratio = sprintf("%.3f", odds_ratio),
    Std_Error = sprintf("%.3f",std_err) ,
    P_Value = p_value_formatted,
    OR_CI_95 = sprintf("(%.3f, %.3f)", or_ci_lower, or_ci_upper),    
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df data frame
  results_df <- rbind(results_df, result)
}


# Write dissimilarity results to the first sheet
writeData(wb, sheet = "Dissimilarity_Results", results_df, colNames = TRUE)





# Iterate over each model of interest for isolation results
for (model_name in models_of_interest) {
  model <- model_results_il1b1v3dat$isolation_index_std_highest_tertile[[model_name]]
  
  beta_value=model$coefficients["isolation_index_std", "Estimate"]
  std_err=model$coefficients["isolation_index_std", "Std. Error"]
  
  p_value = model$coefficients["isolation_index_std", "Pr(>|t|)"]
  p_value_formatted = ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  odds_ratio <- exp(beta_value)
  
  # Calculate margin of error for OR
  or_margin_of_error <- qnorm(0.975) * std_err
  
  # Calculate lower and upper bounds of the CI for OR
  or_ci_lower <- exp(beta_value - or_margin_of_error)
  or_ci_upper <- exp(beta_value + or_margin_of_error)
  
  # Extract the relevant information
  result <- data.frame(
    Model = model_name,
    sd=sd(analyticsamplethr$eselectin),
    Beta = sprintf("%.3f",beta_value),
    Odds_Ratio = sprintf("%.3f", odds_ratio),
    Std_Error = sprintf("%.3f",std_err) ,
    P_Value = p_value_formatted,
    OR_CI_95 = sprintf("(%.3f, %.3f)", or_ci_lower, or_ci_upper),    
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_2 data frame
  results_df_2 <- rbind(results_df_2, result)
}


# Write isolation results to the second sheet
writeData(wb, sheet = "Isolation_Results", results_df_2,colNames = TRUE)

# Iterate over each model of interest for interaction results
for (model_name in models_of_interest) {
  model <- model_results_il1b1v3dat$interaction_index_std_highest_tertile[[model_name]]
  
  beta_value=model$coefficients["interaction_index_std", "Estimate"]
  std_err=model$coefficients["interaction_index_std", "Std. Error"]
  
  p_value = model$coefficients["interaction_index_std", "Pr(>|t|)"]
  p_value_formatted = ifelse(round(p_value, 3) == 0.000, "<0.001", sprintf("%.3f", p_value))
  
  odds_ratio <- exp(beta_value)
  
  # Calculate margin of error for OR
  or_margin_of_error <- qnorm(0.975) * std_err
  
  # Calculate lower and upper bounds of the CI for OR
  or_ci_lower <- exp(beta_value - or_margin_of_error)
  or_ci_upper <- exp(beta_value + or_margin_of_error)
  
  # Extract the relevant information
  result <- data.frame(
    Model = model_name,
    sd=sd(analyticsamplethr$eselectin),
    Beta = sprintf("%.3f",beta_value),
    Odds_Ratio = sprintf("%.3f", odds_ratio),
    Std_Error = sprintf("%.3f",std_err) ,
    P_Value = p_value_formatted,
    OR_CI_95 = sprintf("(%.3f, %.3f)", or_ci_lower, or_ci_upper),    
    stringsAsFactors = FALSE
  )
  
  # Append the result to the results_df_3 data frame
  results_df_3 <- rbind(results_df_3, result)
}



# Write interaction results to the third sheet
writeData(wb, sheet = "Interaction_Results", results_df_3,colNames = TRUE)

# Save the Excel workbook
saveWorkbook(wb, "path/to/output/il1bt1vt3resultsregressindexweighted.xlsx")























