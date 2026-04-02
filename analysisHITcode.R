# =============================================================================
# HIT (Heparin-Induced Thrombocytopenia) Analysis
# Multi-Site Case-Cohort Study
# =============================================================================
# Description: Univariate logistic regression analysis examining patient-level
#              factors associated with anticoagulation agent selection in the
#              day following a positive HIT antibody lab result.
#
# Outcome:     outcome_dayafterhitlabacagent (binary: 1/0)
# Population:  Patients with a confirmed HIT antibody lab order across
#              multiple health system sites
#
# Sections covered:
#   1.  Data load & variable cleaning (race, gender, creatinine, HIT result)
#   2.  Baseline characteristics table
#   3.  Lab result turnaround time (overall and by site)
#   4.  Frequency checks
#   5.  Univariate logistic regression models:
#         - Age (continuous and categorical: <65 vs 65+)
#         - Gender
#         - Site
#         - Admission service
#         - Platelet count
#         - Creatinine category
#         - Anticoagulation intensity (day before HIT lab)
#         - Race
#         - HIT antibody test result category
#
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
library(stringi)
library(readxl)
library(openxlsx) 
library(data.table)
library(dplyr)
library(broom)







########bline char table

combined_hit_analytic=fread('path/to/hitanalytic.csv')
head(combined_hit_analytic)

#fixing race
combined_hit_analytic <- combined_hit_analytic %>%
  mutate(patient_race_category = if_else(patient_race_category == "hispanic", "other", patient_race_category))
combined_hit_analytic <- combined_hit_analytic %>%
  mutate(patient_race_category = if_else(patient_race_category %in% c("other", "missing"), "other/missing", patient_race_category))


combined_hit_analytic <- combined_hit_analytic %>%
  mutate(
    daybeforehitlabacagent  = if_else(is.na(daybeforehitlabacagent)  | daybeforehitlabacagent  == "", "None", daybeforehitlabacagent),
    dayafterhitlabacagent   = if_else(is.na(dayafterhitlabacagent)   | dayafterhitlabacagent   == "", "None", dayafterhitlabacagent)
  )


combined_hit_analytic <- combined_hit_analytic %>%
  mutate(patient_gender = case_when(
    patient_gender %in% c("M", "MALE", "Male") ~ "Male",
    patient_gender %in% c("F", "FEMALE", "Female") ~ "Female",
    TRUE ~ NA_character_  # optional: make anything unexpected NA
  ))


combined_hit_analytic <- combined_hit_analytic %>%
  mutate(
    creatinine_cat = case_when(
      HD == 1 ~ "HD",
      creatinine_cat == "2.0+" & (is.na(HD) | HD == 0) ~ "Cr > 2, NOT on HD",
      creatinine_cat %in% c("<1.2", "1.2-2.0", "<2.0") & (is.na(HD) | HD == 0) ~ "Cr < 2",
      TRUE ~ NA_character_
    )
  )

table(combined_hit_analytic$creatinine_cat)


combined_hit_analytic <- combined_hit_analytic %>%
  mutate(
    lab_result_clean = str_trim(lab_result),
    lab_result_num = suppressWarnings(as.numeric(lab_result_clean)),
    
    hit_ab_result_cat = case_when(
      # Numeric cases
      !is.na(lab_result_num) & lab_result_num < 0.5 ~ "<0.5 OR NEGATIVE",
      !is.na(lab_result_num) & lab_result_num >= 0.5 & lab_result_num < 1 ~ "0.5-1",
      !is.na(lab_result_num) & lab_result_num >= 1 & lab_result_num <= 1.8 ~ "1-1.8",
      !is.na(lab_result_num) & lab_result_num > 1.8 ~ ">1.8",
      
      # Explicit string inequalities like "<0.075"
      str_detect(lab_result_clean, "^<\\s*0\\.\\d+") ~ "<0.5 OR NEGATIVE",
      
      # Text-based labels
      str_detect(tolower(lab_result_clean), "negative") ~ "<0.5 OR NEGATIVE",
      str_detect(tolower(lab_result_clean), "positive") ~ "positive WITHOUT numerical value",
      
      TRUE ~ NA_character_
    )
  )

table(combined_hit_analytic$hit_ab_result_cat)


# combined_hit_analytic <- combined_hit_analytic %>%
#   mutate(lab_result_cathit = case_when(
#     # Exact text matches
#     str_to_upper(str_trim(lab_result)) %in% c("NEGATIVE", "<0.5 OR NEGATIVE") ~ "<0.5 OR NEGATIVE",
#     str_detect(str_to_upper(str_trim(lab_result)), "POSITIVE") ~ "positive WITHOUT numerical value",
#     str_detect(str_to_upper(str_trim(lab_result)), "INDETERM|AGGREGATION|NO AGGREGATION") ~ "positive WITHOUT numerical value",
#     
#     # Numeric conversion and classification
#     suppressWarnings(as.numeric(lab_result)) < 0.5 ~ "<0.5 OR NEGATIVE",
#     suppressWarnings(as.numeric(lab_result)) >= 0.5 & suppressWarnings(as.numeric(lab_result)) < 1 ~ "0.5-1",
#     suppressWarnings(as.numeric(lab_result)) >= 1 & suppressWarnings(as.numeric(lab_result)) <= 1.8 ~ "1-1.8",
#     suppressWarnings(as.numeric(lab_result)) > 1.8 ~ ">1.8",
#     
#     TRUE ~ NA_character_
#   ))

# Load necessary libraries
library(dplyr)
head(combined_hit_analytic)
combined_hit_analytic <- as_tibble(combined_hit_analytic)

# Calculate statistics for the requested table
summary_table <- list(
  age = list(
    mean_sd = combined_hit_analytic %>%
      dplyr::summarize(mean = mean(combined_hit_analytic$patient_age, na.rm = TRUE), sd = sd(combined_hit_analytic$patient_age, na.rm = TRUE)),
    range = range(combined_hit_analytic$patient_age, na.rm = TRUE)
  ),
  sex = combined_hit_analytic %>%
    count(patient_gender) %>%
    mutate(percent = n / sum(n) * 100),
  race = combined_hit_analytic %>%
    count(patient_race_category) %>%
    mutate(percent = n / sum(n) * 100),
  ethnicity = combined_hit_analytic %>%
    count(patient_ethnicity_category) %>%
    mutate(percent = n / sum(n) * 100),
  insurance_category = combined_hit_analytic %>%
    count(primary_payer_type_category) %>%
    mutate(percent = n / sum(n) * 100),
  healthcare_system = combined_hit_analytic %>%
    count(source) %>%
    mutate(percent = n / sum(n) * 100),
  admission_service = combined_hit_analytic %>%
    count(admission_service_category) %>%
    mutate(percent = n / sum(n) * 100)
)

# For the remaining non-yellow variables
remaining_stats <- list(
  platelet_count_day0 = combined_hit_analytic %>%
    count(platelet_count_cat) %>%
    mutate(percent = n / sum(n) * 100),
  creatinine_day1 = combined_hit_analytic %>%
    count(creatinine_cat) %>%
    mutate(percent = n / sum(n) * 100),
  anticoagulation_d1_intensity = combined_hit_analytic %>%
    count(daybeforehitlabaclvl) %>%
    mutate(percent = n / sum(n) * 100),
  anticoagulation_d1_agent = combined_hit_analytic %>%
    count(dayafterhitlabaclevel) %>%
    mutate(percent = n / sum(n) * 100)
)

# Print the results
summary_table
remaining_stats


combined_hit_analytic <- combined_hit_analytic %>%
  mutate(hit_result_lag_days = as.numeric(as.Date(hit_lab_result_date) - as.Date(hit_lab_date)))

# Median and IQR
summary_stats <- combined_hit_analytic %>%
  summarise(
    median = median(hit_result_lag_days, na.rm = TRUE),
    IQR_low = quantile(hit_result_lag_days, 0.25, na.rm = TRUE),
    IQR_high = quantile(hit_result_lag_days, 0.75, na.rm = TRUE),
    min = min(hit_result_lag_days, na.rm = TRUE),
    max = max(hit_result_lag_days, na.rm = TRUE)
  )
print(summary_stats)


summary_stats_by_site <- combined_hit_analytic %>%
  group_by(source) %>%
  summarise(
    median = median(hit_result_lag_days, na.rm = TRUE),
    IQR_low = quantile(hit_result_lag_days, 0.25, na.rm = TRUE),
    IQR_high = quantile(hit_result_lag_days, 0.75, na.rm = TRUE),
    min = min(hit_result_lag_days, na.rm = TRUE),
    max = max(hit_result_lag_days, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats_by_site)

table(combined_hit_analytic$HD)


table(combined_hit_analytic$daybeforehitlabaclvl)


table(combined_hit_analytic$daybeforehitlabacagent)

table(combined_hit_analytic$dayafterhitlabaclevel)


table(combined_hit_analytic$dayafterhitlabacagent)


table(combined_hit_analytic$accatdayprior)



table(combined_hit_analytic$creatinine_cat)

table(combined_hit_analytic$platelet_count_cat)


table(combined_hit_analytic$admission_cat_hit)


combined_hit_analytic <- combined_hit_analytic %>%
  mutate(days_after_admit_to_hitlab = ifelse(days_after_admit_to_hitlab < 0, NA, days_after_admit_to_hitlab))

median_value <- median(combined_hit_analytic$days_after_admit_to_hitlab, na.rm = TRUE)

print(median_value)

table(combined_hit_analytic$admission_cat_hit)


table(combined_hit_analytic$lab_result)

ktableQ<-combined_hit_analytic%>%select(bundle_id,dayafterhitlabacagent,dayafterhitlabaclevel,low_dose_hep)

#define out come.....1/0 if 

# Load required libraries

# combined_hit_analytic=fread('path/to/hitanalytic.csv')
# head(combined_hit_analytic)




### AGE model
age_model <- glm(outcome_dayafterhitlabacagent ~ patient_age, data = combined_hit_analytic, family = binomial)

# Extract OR + CI
age_coef <- summary(age_model)$coefficients
age_or <- exp(age_coef[, "Estimate"])
age_ci_lower <- exp(age_coef[, "Estimate"] - 1.96 * age_coef[, "Std. Error"])
age_ci_upper <- exp(age_coef[, "Estimate"] + 1.96 * age_coef[, "Std. Error"])
age_results <- data.frame(
  Variable = rownames(age_coef),
  OR = age_or,
  CI_Lower = age_ci_lower,
  CI_Upper = age_ci_upper,
  P_Value = age_coef[, "Pr(>|z|)"]
)
print(age_results)


# Create age_cat variable
combined_hit_analytic$age_cat <- ifelse(combined_hit_analytic$patient_age < 65, "<65", "65+")
combined_hit_analytic$age_cat <- factor(combined_hit_analytic$age_cat)
combined_hit_analytic$age_cat <- relevel(combined_hit_analytic$age_cat, ref = "<65")

# Run model
agecat_model <- glm(outcome_dayafterhitlabacagent ~ age_cat, data = combined_hit_analytic, family = binomial)

# Extract OR + CI
agecat_coef <- summary(agecat_model)$coefficients
agecat_or <- exp(agecat_coef[, "Estimate"])
agecat_ci_lower <- exp(agecat_coef[, "Estimate"] - 1.96 * agecat_coef[, "Std. Error"])
agecat_ci_upper <- exp(agecat_coef[, "Estimate"] + 1.96 * agecat_coef[, "Std. Error"])
agecat_results <- data.frame(
  Variable = rownames(agecat_coef),
  OR = agecat_or,
  CI_Lower = agecat_ci_lower,
  CI_Upper = agecat_ci_upper,
  P_Value = agecat_coef[, "Pr(>|z|)"]
)
print(agecat_results)






### GENDER model
combined_hit_analytic$patient_gender <- factor(combined_hit_analytic$patient_gender)
combined_hit_analytic$patient_gender <- relevel(combined_hit_analytic$patient_gender, ref = "Male")

gender_model <- glm(outcome_dayafterhitlabacagent ~ patient_gender, data = combined_hit_analytic, family = binomial)

gender_coef <- summary(gender_model)$coefficients
gender_or <- exp(gender_coef[, "Estimate"])
gender_ci_lower <- exp(gender_coef[, "Estimate"] - 1.96 * gender_coef[, "Std. Error"])
gender_ci_upper <- exp(gender_coef[, "Estimate"] + 1.96 * gender_coef[, "Std. Error"])
gender_results <- data.frame(
  Variable = rownames(gender_coef),
  OR = gender_or,
  CI_Lower = gender_ci_lower,
  CI_Upper = gender_ci_upper,
  P_Value = gender_coef[, "Pr(>|z|)"]
)
print(gender_results)


### SITE model   ####CHANGE TO 1-5


combined_hit_analytic$site <- factor(combined_hit_analytic$source)
table(combined_hit_analytic$site)
combined_hit_analytic$site <- relevel(combined_hit_analytic$site, ref = "Site_1")




site_model <- glm(outcome_dayafterhitlabacagent ~ site, data = combined_hit_analytic, family = binomial)

site_coef <- summary(site_model)$coefficients
site_or <- exp(site_coef[, "Estimate"])
site_ci_lower <- exp(site_coef[, "Estimate"] - 1.96 * site_coef[, "Std. Error"])
site_ci_upper <- exp(site_coef[, "Estimate"] + 1.96 * site_coef[, "Std. Error"])
site_results <- data.frame(
  Variable = rownames(site_coef),
  OR = site_or,
  CI_Lower = site_ci_lower,
  CI_Upper = site_ci_upper,
  P_Value = site_coef[, "Pr(>|z|)"]
)
print(site_results)


### ADMISSION SERVICE model
combined_hit_analytic$admission_cat_hit <- factor(combined_hit_analytic$admission_cat_hit)
combined_hit_analytic$admission_cat_hit <- relevel(combined_hit_analytic$admission_cat_hit, ref = "Non ICU Admission")

adm_model <- glm(outcome_dayafterhitlabacagent ~ admission_cat_hit, data = combined_hit_analytic, family = binomial)

adm_coef <- summary(adm_model)$coefficients
adm_or <- exp(adm_coef[, "Estimate"])
adm_ci_lower <- exp(adm_coef[, "Estimate"] - 1.96 * adm_coef[, "Std. Error"])
adm_ci_upper <- exp(adm_coef[, "Estimate"] + 1.96 * adm_coef[, "Std. Error"])
adm_results <- data.frame(
  Variable = rownames(adm_coef),
  OR = adm_or,
  CI_Lower = adm_ci_lower,
  CI_Upper = adm_ci_upper,
  P_Value = adm_coef[, "Pr(>|z|)"]
)
print(adm_results)


### PLATELET COUNT model
combined_hit_analytic$platelet_count_cat <- factor(combined_hit_analytic$platelet_count_cat)
combined_hit_analytic$platelet_count_cat <- relevel(combined_hit_analytic$platelet_count_cat, ref = "150+")

plt_model <- glm(outcome_dayafterhitlabacagent ~ platelet_count_cat, data = combined_hit_analytic, family = binomial)

plt_coef <- summary(plt_model)$coefficients
plt_or <- exp(plt_coef[, "Estimate"])
plt_ci_lower <- exp(plt_coef[, "Estimate"] - 1.96 * plt_coef[, "Std. Error"])
plt_ci_upper <- exp(plt_coef[, "Estimate"] + 1.96 * plt_coef[, "Std. Error"])
plt_results <- data.frame(
  Variable = rownames(plt_coef),
  OR = plt_or,
  CI_Lower = plt_ci_lower,
  CI_Upper = plt_ci_upper,
  P_Value = plt_coef[, "Pr(>|z|)"]
)
print(plt_results)


### CREATININE model
combined_hit_analytic$creatinine_cat <- factor(combined_hit_analytic$creatinine_cat)
combined_hit_analytic$creatinine_cat <- relevel(combined_hit_analytic$creatinine_cat, ref = "<2.0")

cr_model <- glm(outcome_dayafterhitlabacagent ~ creatinine_cat, data = combined_hit_analytic, family = binomial)

cr_coef <- summary(cr_model)$coefficients
cr_or <- exp(cr_coef[, "Estimate"])
cr_ci_lower <- exp(cr_coef[, "Estimate"] - 1.96 * cr_coef[, "Std. Error"])
cr_ci_upper <- exp(cr_coef[, "Estimate"] + 1.96 * cr_coef[, "Std. Error"])
cr_results <- data.frame(
  Variable = rownames(cr_coef),
  OR = cr_or,
  CI_Lower = cr_ci_lower,
  CI_Upper = cr_ci_upper,
  P_Value = cr_coef[, "Pr(>|z|)"]
)
print(cr_results)


### ANTICOAGULATION INTENSITY model
combined_hit_analytic$daybeforehitlabaclvl <- factor(combined_hit_analytic$daybeforehitlabaclvl)
combined_hit_analytic$daybeforehitlabaclvl <- relevel(combined_hit_analytic$daybeforehitlabaclvl, ref = "Prophylactic")

ac_model <- glm(outcome_dayafterhitlabacagent ~ daybeforehitlabaclvl, data = combined_hit_analytic, family = binomial)

ac_coef <- summary(ac_model)$coefficients
ac_or <- exp(ac_coef[, "Estimate"])
ac_ci_lower <- exp(ac_coef[, "Estimate"] - 1.96 * ac_coef[, "Std. Error"])
ac_ci_upper <- exp(ac_coef[, "Estimate"] + 1.96 * ac_coef[, "Std. Error"])
ac_results <- data.frame(
  Variable = rownames(ac_coef),
  OR = ac_or,
  CI_Lower = ac_ci_lower,
  CI_Upper = ac_ci_upper,
  P_Value = ac_coef[, "Pr(>|z|)"]
)
print(ac_results)


### RACE model
combined_hit_analytic$patient_race_category <- factor(combined_hit_analytic$patient_race_category)
combined_hit_analytic$patient_race_category <- relevel(combined_hit_analytic$patient_race_category, ref = "white")

race_model <- glm(outcome_dayafterhitlabacagent ~ patient_race_category, data = combined_hit_analytic, family = binomial)

race_coef <- summary(race_model)$coefficients
race_or <- exp(race_coef[, "Estimate"])
race_ci_lower <- exp(race_coef[, "Estimate"] - 1.96 * race_coef[, "Std. Error"])
race_ci_upper <- exp(race_coef[, "Estimate"] + 1.96 * race_coef[, "Std. Error"])
race_results <- data.frame(
  Variable = rownames(race_coef),
  OR = race_or,
  CI_Lower = race_ci_lower,
  CI_Upper = race_ci_upper,
  P_Value = race_coef[, "Pr(>|z|)"]
)
print(race_results)


### HIT TEST RESULT model
combined_hit_analytic$hit_ab_result_cat <- factor(combined_hit_analytic$hit_ab_result_cat)
combined_hit_analytic$hit_ab_result_cat <- relevel(combined_hit_analytic$hit_ab_result_cat, ref = "0.5-1")

hit_model <- glm(outcome_dayafterhitlabacagent ~ hit_ab_result_cat, data = combined_hit_analytic, family = binomial)

hit_coef <- summary(hit_model)$coefficients
hit_or <- exp(hit_coef[, "Estimate"])
hit_ci_lower <- exp(hit_coef[, "Estimate"] - 1.96 * hit_coef[, "Std. Error"])
hit_ci_upper <- exp(hit_coef[, "Estimate"] + 1.96 * hit_coef[, "Std. Error"])
hit_results <- data.frame(
  Variable = rownames(hit_coef),
  OR = hit_or,
  CI_Lower = hit_ci_lower,
  CI_Upper = hit_ci_upper,
  P_Value = hit_coef[, "Pr(>|z|)"]
)
print(hit_results)


#platelet count
# Set reference for platelet count category
combined_hit_analytic$platelet_count_cat <- factor(combined_hit_analytic$platelet_count_cat)
table(combined_hit_analytic$platelet_count_cat)
combined_hit_analytic$platelet_count_cat <- relevel(combined_hit_analytic$platelet_count_cat, ref = "150+")

# Run logistic regression
platelet_model <- glm(outcome_dayafterhitlabacagent ~ platelet_count_cat, data = combined_hit_analytic, family = binomial)

# Extract coefficients and calculate ORs and CIs
platelet_coef <- summary(platelet_model)$coefficients
platelet_or <- exp(platelet_coef[, "Estimate"])
platelet_ci_lower <- exp(platelet_coef[, "Estimate"] - 1.96 * platelet_coef[, "Std. Error"])
platelet_ci_upper <- exp(platelet_coef[, "Estimate"] + 1.96 * platelet_coef[, "Std. Error"])
platelet_results <- data.frame(
  Variable = rownames(platelet_coef),
  OR = round(platelet_or, 2),
  CI_Lower = round(platelet_ci_lower, 2),
  CI_Upper = round(platelet_ci_upper, 2),
  P_Value = signif(platelet_coef[, "Pr(>|z|)"], 3)
)

# Print formatted table
print(platelet_results)









































































