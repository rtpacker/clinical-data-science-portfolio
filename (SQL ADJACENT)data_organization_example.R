# =============================================================================
# Level 2 Data Organization — Hospital Encounter Pipeline
# =============================================================================
# Description: Processes Level 1 raw hospital data into a clean Level 2
#              analytic dataset. Each section loads a specific data domain,
#              applies filtering and date alignment to the study population,
#              and saves an output file.
#
# Data domains covered:
#   1.  Encounter information
#   2.  Vital signs
#   3.  Procedure codes
#   4.  Diagnosis (ICD)
#   5.  Laboratory results (with reference-range categorization)
#   6.  Transfusions (flagged for ≥2 units/48hrs)
#   7.  Hemoglobin drop (flagged for ≥2 g/dL drop within 24 or 48hrs)
#   8.  Height & weight (BMI calculation)
#   9.  Problem list
#   10. Medical history
#   11. Respiratory support
#   12. Inpatient medications
#   13. ICU stays
#
#
# =============================================================================


# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------
library(dplyr)
library(data.table)
library(zoo)
library(tidyr)
library(lubridate)
library(tidyverse)
library(knitr)
library(kableExtra)
library(haven)
library(nimble)
library(stringi)
library(readxl)


# -----------------------------------------------------------------------------
# Paths  — update before running
# -----------------------------------------------------------------------------
data.path  <- "path/to/data/"         # root data directory
look.path  <- "path/to/lookup_data/"  # lookup/reference tables directory

l1 <- function(f) sprintf('%sClean Data/Level 1/%s', data.path, f)
l2 <- function(f) sprintf('%sClean Data/Level 2/%s',  data.path, f)


# =============================================================================
# Helper Functions
# =============================================================================

# -----------------------------------------------------------------------------
# flag_sum_hour()
# Flags admissions where a lab value drops or accumulates beyond a threshold
# within a specified rolling time window.
#
# Args:
#   result_value : numeric vector of measurements (e.g., transfusion units, hgb)
#   time         : corresponding POSIXct or numeric timestamps
#   unit_sum     : minimum cumulative sum to flag (e.g., 2 units transfused)
#   unit_drop    : minimum drop in value to flag (e.g., 2 g/dL hgb drop)
#   time_days    : rolling window width in days
#   unit_vol     : volume per unit for sum conversion (default 250 mL)
#   bundle_id    : optional, used for debugging
#
# Returns: data.frame with columns sum_result_value, flag, time_flag
# -----------------------------------------------------------------------------
flag_sum_hour <- function(result_value, time,
                          unit_sum = NA, unit_drop = NA, time_days,
                          unit_vol = 250, bundle_id = NA) {
  if (!is.na(bundle_id)) print(bundle_id)

  flag             <- rep(0, length(time))
  sum_result_value <- rep(0, length(time))
  time_flag        <- rep(0, length(time))

  # Compute lower-triangular time-difference matrix
  t24 <- outer(time, time, `-`)
  t24[upper.tri(t24)] <- NA
  if (class(time)[[1]] == 'POSIXct') t24 <- (t24 / 3600) / 24

  which24 <- as.data.table(which(t24 <= time_days, arr.ind = TRUE))

  if (!is.na(unit_sum)) {
    # Cumulative sum branch (e.g., transfusions)
    which24$sum24 <- apply(which24, 1, function(x) sum(result_value[x[2]:x[1]]))
    which24max <- which24 %>%
      group_by(col) %>% mutate(max24 = max(sum24)) %>%
      select(col, max24) %>% distinct()
    which24max$flag_max <- as.numeric(which24max$max24 >= (unit_sum * unit_vol))
    first <- which(which24max$flag_max == 1)[1]
    which24max$time_first <- if (!is.na(first)) time[first] else NA
  }

  if (!is.na(unit_drop)) {
    # Value-drop branch (e.g., hemoglobin)
    which24$diff24 <- apply(which24, 1, function(x) result_value[x[2]] - result_value[x[1]])
    which24max <- which24 %>%
      group_by(col) %>% mutate(max24 = max(diff24)) %>%
      select(col, max24) %>% distinct()
    which24max$flag_max <- as.numeric(which24max$max24 >= unit_drop)
    first <- which(which24max$flag_max == 1)[1]
    if (!is.na(first))                             which24max$time_first <- time[first]
    else if (nrow(which24max) == 0)                which24max$time_first <- numeric()
    else                                           which24max$time_first <- NA
  }

  flag[which24max$col]             <- which24max$flag_max
  sum_result_value[which24max$col] <- which24max$max24
  time_flag[which24max$col]        <- which24max$time_first
  if (class(time)[[1]] == 'POSIXct') time_flag <- as_datetime(time_flag)

  return(data.frame(sum_result_value, flag, time_flag))
}


# -----------------------------------------------------------------------------
# sum_drop_summarise()
# Summarises the flag_sum_hour() output to one row per admission, capturing
# the first event time and whether the event occurred within 24hrs of admission.
#
# Args:
#   result_value   : numeric vector of measurements
#   time           : corresponding timestamps
#   admission_date : admission date/time for the bundle
#   flag, sum_lab, time_flag : outputs from flag_sum_hour()
#   unit_sum / unit_drop     : thresholds (one should be non-NA)
#   time_days, unit_vol      : passed through to label output columns
#
# Returns: one-row data.table with summary statistics
# -----------------------------------------------------------------------------
sum_drop_summarise <- function(result_value, time,
                               admission_date = NULL,
                               flag = NULL, sum_lab = NULL, time_flag = NULL,
                               unit_sum = NA, unit_drop = NA,
                               time_days = 2, unit_vol = 250) {
  if (is.null(admission_date)) stop('admission_date required')
  if (is.null(flag) | is.null(sum_lab) | is.null(time_flag)) {
    stop('flags not provided — run flag_sum_hour() first')
  }

  sum_drop      <- data.table(flag, sum_lab, time_flag)
  first_meas      <- result_value[1]
  first_meas_time <- time[1]
  any_flag        <- as.numeric(any(sum_drop$flag == 1))
  max_meas        <- max(sum_drop$sum_lab)
  min_meas        <- min(sum_drop$sum_lab[sum_lab > 0])

  admission_date <- unique(admission_date)
  if (length(admission_date) > 1) stop('bundle_id has multiple admission dates')
  time_first <- unique(sum_drop$time_flag)
  if (length(time_first) > 1) stop('bundle_id has multiple first-event times')

  if (class(time_first)[[1]] == 'numeric' & class(admission_date)[[1]] == 'numeric') {
    flag_first24hr <- as.numeric(time_first - admission_date <= 1)
  }
  if (class(time_first)[[1]] == 'POSIXct' & class(admission_date)[[1]] == 'POSIXct') {
    flag_first24hr <- as.numeric(difftime(time_first, admission_date, units = 'hours') <= 24)
  }

  ret <- data.table(first_meas, first_meas_time, any_flag, max_meas, min_meas, flag_first24hr, time_first)
  ret[min_meas == Inf, min_meas := 0]

  if (!is.na(unit_drop)) {
    setnames(ret, c('first_meas', 'first_meas_time', 'any_drop', 'max_drop',
                    'min_drop', 'drop_first24hr', 'start_time_first_drop'))
    ret[, unit_drop := unit_drop]
  }
  if (!is.na(unit_sum)) {
    setnames(ret, c('first_meas', 'first_meas_time', 'any_sum', 'max_sum',
                    'min_sum', 'sum_first24hr', 'start_time_first_sum'))
    ret[, unit_sum := unit_sum]
    ret[, unit_vol := unit_vol]
  }
  ret[, rolling_time_days := time_days]
  ret[, admission_date := admission_date]
  return(ret)
}


# -----------------------------------------------------------------------------
# labs_cat_func()
# Categorizes a lab result into clinically meaningful bins based on lab type
# and (for hemoglobin) patient gender.
#
# Args:
#   lab_cat    : character, one of the lab categories in labs_keep
#   lab_result : character or numeric result value
#   gender     : character vector ('Male' / 'Female')
#
# Returns: factor with labelled cut categories, or "" if lab_cat not recognized
# -----------------------------------------------------------------------------
labs_cat_func <- function(lab_cat, lab_result, gender) {
  print(lab_cat)
  result <- as.numeric(lab_result)

  cuts <- switch(lab_cat,
    'White Cell Count' = list(
      breaks = c(0, 11.1, Inf), labels = c('<11.1', '11.1+')
    ),
    'Platelet Count' = list(
      breaks = c(0, 10, 25, 50, 100, 150, 350, 450, 750, Inf),
      labels = c('<10', '10-24', '25-49', '50-99', '100-149', '150-349', '350-449', '450-749', '750+')
    ),
    'Red Cell Distribution Width Coefficient of Variation' = list(
      breaks = c(0, 14.7, Inf), labels = c('<14.7%', '14.7%+')
    ),
    'Mean Corpuscular Volume' = list(
      breaks = c(0, 70, 82, 86, 99, Inf), labels = c('<70', '70-81', '82-85', '86-98', '99+')
    ),
    'Serum Sodium'            = list(breaks = c(0, 136, 145, Inf), labels = c('<136', '136-144', '145+')),
    'Serum Potassium'         = list(breaks = c(0, 3.5, 5.0, Inf), labels = c('<3.5', '3.5-5.0', '5.1+')),
    'Serum Bicarbonate'       = list(breaks = c(0, 22, 32, Inf),   labels = c('<22', '22-31', '32+')),
    'Serum Blood Urea Nitrogen' = list(breaks = c(0, 5, 20, Inf),  labels = c('<5', '5-20', '20+')),
    'Serum Creatinine'        = list(breaks = c(0, 2, Inf),        labels = c('<2.0', '2.0+')),
    NULL
  )

  # Hemoglobin uses gender-specific thresholds
  if (lab_cat == 'Hemoglobin') {
    if (all(gender == 'Female')) {
      return(cut(result, breaks = c(0, 12.0, Inf), include.lowest = TRUE,
                 right = FALSE, labels = c('low', 'high')))
    } else if (all(gender == 'Male')) {
      return(cut(result, breaks = c(0, 13.6, Inf), include.lowest = TRUE,
                 right = FALSE, labels = c('low', 'high')))
    } else {
      return(cut(result, breaks = c(0, 12, 13.6, Inf), include.lowest = TRUE,
                 right = FALSE, labels = c('low', 'No or unclear Gender', 'high')))
    }
  }

  if (!is.null(cuts)) {
    return(cut(result, breaks = cuts$breaks, labels = cuts$labels,
               include.lowest = TRUE, right = FALSE))
  }
  return("")
}

labs_keep <- c('White Cell Count', 'Hemoglobin', 'Platelet Count',
               'Red Cell Distribution Width Standard Deviation',
               'Red Cell Distribution Width Coefficient of Variation',
               'Mean Corpuscular Volume', 'Serum Sodium', 'Serum Potassium',
               'Serum Bicarbonate', 'Serum Blood Urea Nitrogen', 'Serum Creatinine')


# =============================================================================
# SECTION 1: Encounter Information
# =============================================================================
ref       <- fread(l1('level1_reference.csv'))
ref_mith  <- ref[MITH_med_population == 1, .(patid, encounter_id, emergency_icu_admit, bundle_id)]

enc_l1 <- as.data.table(readRDS(l1('Encounter_Info.rds')))
enc_l1 <- enc_l1 %>% rename(patid = PATID)

# Merge to study population and sort by admission date
enc_mith <- unique(merge(enc_l1, ref_mith, by = c('patid', 'encounter_id', 'bundle_id'))[order(bundle_id, admission_date)])

# Recode emergency medicine service to use discharge service when not ICU
enc_mith[admission_service == 'EMERGENCY MEDICINE' & admission_service_category != 'Intensive Care Unit',
         admission_service := discharge_service]

# Collapse multi-encounter bundles to one row, capturing first/last values
enc_mith[, `:=`(
  adate_bundle                     = min(admission_date),
  ddate_bundle                     = max(discharge_date),
  first_encounter_type             = head(encounter_type, 1),
  first_admission_service          = head(admission_service, 1),
  first_admission_service_category = head(admission_service_category, 1),
  patient_age                      = head(patient_age, 1)
), .(bundle_id)]

enc_mith[, c('last_encounter_type', 'last_discharge_service',
             'last_discharge_service_category', 'last_discharge_disposition') :=
           lapply(.SD, last),
         .SDcols = c('encounter_type', 'discharge_service',
                     'discharge_service_category', 'discharge_disposition'),
         .(bundle_id)]

# Flag in-hospital death
enc_mith[, discharge_deceased := ifelse(last_discharge_disposition == 'EXPIRED', 1, 0), .(bundle_id)]

# Collapse to one row per bundle
enc_uniq <- enc_mith[order(bundle_id, admission_date, discharge_date), .SD[1], .(bundle_id)]
enc_uniq <- enc_uniq[, .(
  patid, bundle_id, adate_bundle, ddate_bundle,
  first_encounter_type, last_encounter_type,
  first_admission_service, last_discharge_service,
  admission_service_category  = first_admission_service_category,
  discharge_service_category  = last_discharge_service_category,
  last_discharge_disposition, discharge_deceased,
  patient_age, primary_payer_type,
  primary_payer_type_category = NA,
  zip_code
)]

# Keep only admissions with at least one overnight stay (LOS >= 1 day)
enc_uniq[, los_days := as.Date(ddate_bundle) - as.Date(adate_bundle)]
enc <- enc_uniq[los_days >= 1]


# =============================================================================
# SECTION 2: Vital Signs
# =============================================================================
vitals_l1 <- readRDS(l1('Vital_Signs_Flowsheet.rds'))

vitals_l2_func <- function(vitals_raw, enc_data) {
  # Drop records without a vital category
  vt <- vitals_raw[!is.na(vital_category) & vital_category != '']

  # Merge to encounter cohort and keep vitals within admission window
  vt <- merge(vt[, bundle_id := NULL],
              enc_data[, .(patid, adate_bundle, ddate_bundle, bundle_id)],
              by = 'patid', allow.cartesian = TRUE)
  vt <- vt[between(date(recorded_date), date(adate_bundle) - days(1), date(ddate_bundle))]

  vt <- vt[, .(patid, encounter_id, bundle_id,
               flowsheet_name      = template_flowsheet_name,
               vital_name          = template_flow_measure_name,
               vital_category,
               vital_result        = flow_measure_result,
               vital_result_units  = measure_result_unit,
               recorded_date, adate_bundle, ddate_bundle)]

  # Identify time of first core vital per bundle
  core_vital_cats <- c('Blood Pressure (Both Systolic and Diastolic)', 'Temperature',
                       'Heart Rate', 'Respiratory Rate',
                       'Systolic Blood Pressure', 'Diastolic Blood Pressure')

  vitals_first <- vt[vital_result != '' & !is.na(vital_result) & vital_category %in% core_vital_cats] %>%
    .[, min_vital := min(recorded_date), .(bundle_id)] %>%
    .[recorded_date == min_vital]
  vitals_first <- unique(vitals_first[, .(patid, bundle_id, adate_bundle, ddate_bundle,
                                          first_vital_time = min_vital)])

  saveRDS(vt, l2('MITH_vitals.rds'))
  return(vitals_first)
}

vitals_first <- vitals_l2_func(vitals_l1, enc)
fwrite(vitals_first, l2('MITH_firstvitals.csv'))

# Adjust admission date to first vital time if earlier than recorded admission
enc_fv <- merge(enc, vitals_first[, .(bundle_id, first_vital_time)], by = 'bundle_id', all.x = TRUE)
enc_fv[, adate_bundle_original := adate_bundle]
enc_fv[, adate_bundle := ifelse(is.na(first_vital_time), adate_bundle_original,
                                min(adate_bundle_original, first_vital_time)), .(bundle_id)]
enc_fv[, adate_bundle := as_datetime(adate_bundle)]
enc_fv[, los_days := difftime(date(ddate_bundle), date(adate_bundle), units = 'days')]

saveRDS(enc_fv, l2('MITH_medicine_encounters.rds'))

# Summary of admission services
admission_service_counts <- enc_fv %>%
  group_by(admission_service_category) %>%
  summarise(count = n(), percentage = n() / nrow(enc_fv) * 100) %>%
  arrange(desc(count))
print(admission_service_counts)


# =============================================================================
# SECTION 3: Procedure Information
# =============================================================================
ref       <- fread(l1('level1_reference.csv'))
ref_mith  <- subset(ref, MITH_med_population == 1)
enc       <- readRDS(l2('MITH_medicine_encounters.rds'))

proc_raw <- readRDS(l1('Procedure_info.rds'))
proc_raw <- proc_raw[PATID %in% ref_mith$patid]
proc     <- select(proc_raw, -bundle_id) %>% rename(patid = PATID)

# Align procedures to encounter windows
proc <- merge(proc, enc, by = 'patid', allow.cartesian = TRUE)
proc <- subset(proc, as.Date(adate_bundle) <= procedure_date & procedure_date <= as.Date(ddate_bundle))

# Determine procedure code type (CPT > ICD > HCPCS > other)
proc$procedure_code_type <- ifelse(!is.na(proc$cpt_procedure_code), 'CPT',
                              ifelse(!is.na(proc$icd_procedure_code), 'ICD',
                                ifelse(!is.na(proc$HCPCS_procedure_code), 'HCPCS',
                                  'other/Internal procedure code')))

proc$procedure_code <- NA
proc$procedure_code[!is.na(proc$cpt_procedure_code)  & proc$cpt_procedure_code  != 'NULL'] <- proc$cpt_procedure_code[!is.na(proc$cpt_procedure_code)  & proc$cpt_procedure_code  != 'NULL']
proc$procedure_code[is.na(proc$procedure_code) & !is.na(proc$icd_procedure_code)  & proc$icd_procedure_code  != 'NULL'] <- proc$icd_procedure_code[!is.na(proc$icd_procedure_code)  & proc$icd_procedure_code  != 'NULL']
proc$procedure_code[is.na(proc$procedure_code) & !is.na(proc$other_procedure_code) & proc$other_procedure_code != 'NULL'] <- proc$other_procedure_code[!is.na(proc$other_procedure_code) & proc$other_procedure_code != 'NULL']

proc_keep <- select(proc, PATID = patid, bundle_id, procedure_date, procedure_name,
                    procedure_code, procedure_code_type)
saveRDS(proc_keep, l2('MITH_procedure_info.rds'))


# =============================================================================
# SECTION 4: Diagnosis (ICD)
# =============================================================================
ref      <- fread(l1('level1_reference.csv'))
ref_mith <- subset(ref, MITH_med_population == 1)

dx <- readRDS(l1('diagnosis_icd_info.rds'))
dx <- dx[bundle_id %in% ref_mith$bundle_id] %>% select(-encounter_id)

saveRDS(dx, l2('diagnosis_icd_info.rds'))


# =============================================================================
# SECTION 5: Laboratory Results
# =============================================================================
ref      <- fread(l1('level1_reference.csv'))
ref_mith <- subset(ref, mith_med_population == 1)
enc      <- readRDS(l2('MITH_medicine_encounters.rds'))

# Combine annual lab files
lab_files <- paste0('Lab_info_', 2019:2022, '.rds')
labs <- rbindlist(lapply(lab_files, function(f) readRDS(l1(f)))) %>%
  rename(patid = PATID)

# Filter to study population and non-blank categories
labs <- subset(labs, lab_category != 'NA' & lab_category != '' & patid %in% ref_mith$patid)

# Align labs to encounter windows (up to 1.5 days before admission)
labs <- select(labs, -bundle_id)
labs <- merge(labs, enc, by = 'patid', allow.cartesian = TRUE)
labs$lab_drawn_date <- labs$lab_result_date
labs <- subset(labs, lab_drawn_date <= ddate_bundle & adate_bundle - 129600 <= lab_drawn_date)

# Standardize RDW category name
labs[lab_category == 'Red Cell Distribution Width Standard Deviation',
     lab_category := 'Red Cell Distribution Width Coefficient of Variation']

# Select output columns
labs <- select(labs, PATID = patid, bundle_id, lab_name, lab_category,
               lab_ordering_date, lab_drawn_date, lab_result_date,
               lab_result, lab_result_unit, lab_reference_value, lab_comments)

# Merge patient gender for hemoglobin categorization
dem         <- readRDS(l1('Patient_demographics.rds'))
enc_merge   <- readRDS(l2('MITH_medicine_encounters.rds'))[, .(patid, bundle_id, adate_bundle, ddate_bundle)]
gender      <- unique(dem[PATID %in% unique(enc_merge$patid), .(PATID, patient_gender)])
labs        <- merge(labs, gender, by = 'PATID', all.x = TRUE)

# Apply clinical categorization
labs[lab_category %in% labs_keep,
     lab_result_cat := labs_cat_func(lab_cat = lab_category,
                                      lab_result = lab_result,
                                      gender = patient_gender),
     .(lab_category, bundle_id)]

labs_final <- select(labs, PATID, bundle_id, lab_name, lab_category,
                     lab_ordering_date, lab_drawn_date, lab_result_date,
                     lab_result, lab_result_unit, lab_result_cat,
                     lab_reference_value, lab_comments)

saveRDS(labs_final, l2('MITH_lab_info_with_result_cat.rds'))


# =============================================================================
# SECTION 6: Transfusions
# =============================================================================
tran  <- readRDS(l1('Transfusions.rds')) %>% rename(patid = PATID)
enc   <- as.data.table(readRDS(l2('MITH_medicine_encounters.rds')))
tran  <- tran[patid %in% enc$patid]

# Set standard volume for "TRANSFUSE RED BLOOD CELLS" entries
tran_rbc <- subset(tran, TEMPLATE_FLOW_MEASURE_NAME == 'TRANSFUSE RED BLOOD CELLS')
tran_rbc$FLOW_MEASURE_RESULT  <- 250
tran_rbc$MEASURE_RESULT_UNIT  <- 'mL'
tran_rbc$result_n             <- 250

# Merge to encounter windows
tran_rbc$FLOW_MEASURE_RESULT  <- as.numeric(tran_rbc$FLOW_MEASURE_RESULT)
tran_rbc$transfusion_in_units <- tran_rbc$FLOW_MEASURE_RESULT / 250

tran_rbc <- merge(tran_rbc, enc, by = 'patid', allow.cartesian = TRUE)
tran_rbc$transfuseenddatetime <- as.POSIXct(tran_rbc$blood_end_instant)

# Keep transfusions within 24hrs before through end of admission
tf_data <- tran_rbc %>%
  filter(transfuseenddatetime >= adate_bundle - hours(24) &
           transfuseenddatetime <= ddate_bundle)

# Flag admissions with >= 2 units transfused within 48hrs
tf_sum <- tf_data[order(patid, bundle_id.y, transfuseenddatetime)] %>%
  .[, c('vol_48hr', 'flag_48hr', 'time_tf_48hr') :=
      flag_sum_hour(result_value = transfusion_in_units,
                    time         = transfuseenddatetime,
                    unit_sum     = 2, unit_vol = 1, time_days = 2),
    .(patid, bundle_id.y)]

tf_summary <- tf_sum[order(bundle_id.y, transfuseenddatetime)] %>%
  .[, sum_drop_summarise(result_value  = transfusion_in_units,
                         time          = transfuseenddatetime,
                         unit_sum      = 2, time_days = 2, unit_vol = 1,
                         admission_date = adate_bundle,
                         flag = flag_48hr, sum_lab = vol_48hr, time_flag = time_tf_48hr),
    .(patid, bundle_id.y)]

saveRDS(tf_summary, l2('MITH_transfusions_flagged.rds'))


# =============================================================================
# SECTION 7: Hemoglobin Drop
# =============================================================================
mithlab  <- readRDS(l2('MITH_lab_info_with_result_cat.rds'))
enc      <- as.data.table(readRDS(l2('MITH_medicine_encounters.rds')))

hgb_data <- subset(mithlab, lab_category == 'Hemoglobin')
hgb_data$lab_result <- as.numeric(hgb_data$lab_result)
hgb_data <- merge(hgb_data, enc[, .(bundle_id, adate_bundle, ddate_bundle)], by = 'bundle_id')
hgb_data <- hgb_data[adate_bundle - hours(24) <= lab_result_date & lab_result_date <= ddate_bundle]

# Flag >= 2 g/dL drop within 48hrs and within 24hrs
hgb_drop <- hgb_data[order(PATID, bundle_id, lab_result_date)] %>%
  .[, c('drop_48hr', 'flag_48hr', 'time_drop_48hr') :=
      flag_sum_hour(result_value = lab_result, time = lab_result_date,
                    unit_drop = 2, time_days = 2),
    .(PATID, bundle_id)]

hgb_drop24 <- hgb_data[order(PATID, bundle_id, lab_result_date)] %>%
  .[, c('drop_24hr', 'flag_24hr', 'time_drop_24hr') :=
      flag_sum_hour(result_value = lab_result, time = lab_result_date,
                    unit_drop = 2, time_days = 1),
    .(PATID, bundle_id)]

# Merge 24hr and 48hr flags and summarise
enc_dates <- enc[, .(bundle_id, adate_bundle)]
hgb_both  <- merge(hgb_drop, hgb_drop24,
                   by = c('PATID', 'bundle_id', 'lab_name', 'lab_category',
                          'lab_ordering_date', 'lab_drawn_date', 'lab_result_date',
                          'lab_result', 'lab_result_unit', 'lab_reference_value', 'lab_comments'))
hgb_both  <- merge(hgb_both, enc_dates, by = 'bundle_id')

hgb_summary <- hgb_both[order(bundle_id, lab_result_date)] %>%
  .[, sum_drop_summarise(result_value   = lab_result,
                         time           = lab_result_date,
                         unit_drop      = 2, time_days = 2,
                         admission_date = adate_bundle,
                         flag = flag_48hr, sum_lab = drop_48hr, time_flag = time_drop_48hr),
    .(PATID, bundle_id)]

saveRDS(hgb_summary, l2('MITH_hgb_drop_flagged.rds'))


# =============================================================================
# SECTION 8: Height & Weight (BMI)
# =============================================================================
enc      <- readRDS(l2('MITH_medicine_encounters.rds')) %>% rename(PATID = patid)
hw_raw   <- readRDS(l1('HeightWeight.rds'))
hw_merge <- merge(hw_raw, enc, by = 'PATID', allow.cartesian = TRUE)

# Height: filter to plausible values and closest measurement to admission
Height <- hw_merge %>%
  select(PATID, encounter_id, bundle_id = bundle_id.y, height, height_units,
         height_recorded_date, admission_date = adate_bundle, discharge_date = ddate_bundle)
Height$height[Height$height < 50 | Height$height > 250] <- NA
Height <- na.omit(Height)
Height$height_recorded_date <- ymd_hms(Height$height_recorded_date)
Height <- Height[Height$height_recorded_date <= Height$admission_date + years(5) &
                   Height$height_recorded_date >= Height$discharge_date - years(5), ]
Height$date_diff <- abs(Height$height_recorded_date - Height$admission_date)
Height <- Height %>%
  arrange(bundle_id, admission_date, height_recorded_date) %>%
  group_by(bundle_id) %>% filter(!is.na(height)) %>%
  slice(which.min(date_diff)) %>% ungroup() %>% select(-date_diff)

# Weight: filter to plausible values and closest measurement to admission
Weight <- hw_merge %>%
  select(PATID, encounter_id, bundle_id = bundle_id.y, weight, weight_units,
         weight_recorded_date, admission_date = adate_bundle, discharge_date = ddate_bundle)
Weight$weight[Weight$weight < 30 | Weight$weight > 250] <- NA
Weight <- na.omit(Weight)
Weight$weight_recorded_date <- ymd_hms(Weight$weight_recorded_date)
Weight <- Weight[Weight$weight_recorded_date <= Weight$discharge_date &
                   Weight$weight_recorded_date >= Weight$admission_date - months(3), ]
Weight$date_diff <- abs(Weight$weight_recorded_date - Weight$admission_date)
Weight <- Weight %>%
  arrange(bundle_id, admission_date, weight_recorded_date) %>%
  group_by(bundle_id) %>% filter(!is.na(weight)) %>%
  slice(which.min(date_diff)) %>% ungroup() %>% select(-date_diff)

# Merge and calculate BMI
HW <- merge(Height, Weight, by = 'bundle_id', allow.cartesian = TRUE) %>%
  mutate(bmi = weight / ((height / 100)^2)) %>%
  select(PATID = PATID.x, bundle_id, height, weight, bmi,
         height_recorded_date, weight_recorded_date)

saveRDS(HW, l2('MITH_height_weight.rds'))


# =============================================================================
# SECTION 9: Problem List
# =============================================================================
enc      <- readRDS(l2('MITH_medicine_encounters.rds'))
prob_raw <- readRDS(l1('Problem_List.rds')) %>% rename(patid = PATID)

# Expand comma-separated ICD10 codes to one row per code
prob_raw <- prob_raw %>% separate_rows(problem_list_icd10, sep = ',')

# Align to encounters and keep only problems documented before admission
prob_raw <- select(prob_raw, -bundle_id)
prob_raw <- merge(prob_raw, enc, by = 'patid', allow.cartesian = TRUE)
prob_raw <- subset(prob_raw, problem_list_entry_date < adate_bundle)
prob_raw <- subset(prob_raw, problem_list_status %in% c('ACTIVE', 'RESOLVED'))

prob_final <- prob_raw %>%
  select(PATID = patid, bundle_id, problem_list_icd10, problem_list_status,
         chronic_yn, problem_list_class, problem_list_entry_date,
         problem_list_date, problem_list_noted_date)

saveRDS(prob_final, l2('MITH_problem_list.rds'))


# =============================================================================
# SECTION 10: Medical History
# =============================================================================
enc        <- readRDS(l2('MITH_medicine_encounters.rds'))
enc_merge  <- enc[, .(PATID = patid, adate_bundle, ddate_bundle)]
med_hist   <- readRDS(l1('Medical_Diagnosis_History.rds'))
med_hist   <- merge(med_hist, enc_merge, by = 'PATID', allow.cartesian = TRUE)

# Keep only history recorded before or at admission (+12hr buffer)
vitals_first <- fread(l2('MITH_firstvitals.csv'))
names(med_hist) <- tolower(names(med_hist))
med_hist_merged <- merge(med_hist, vitals_first, by = 'patid', allow.cartesian = TRUE)
med_hist_f <- med_hist_merged %>%
  filter(medical_history_entry_date <= (first_vital_time + hours(12)))

saveRDS(med_hist, l2('MITH_medical_history.rds'))


# =============================================================================
# SECTION 11: Respiratory Support
# =============================================================================
enc       <- readRDS(l2('MITH_medicine_encounters.rds'))
enc_merge <- as.data.table(enc)[, .(patid, adate_bundle, ddate_bundle)]

resp_raw <- readRDS(l1('Respiratory_Support_flowsheet.rds'))
resp_raw$flow_measure_result_units <- NA
resp_raw$flowsheet_name            <- NA

resp_raw  <- merge(resp_raw, enc_merge, by = 'patid', allow.cartesian = TRUE)
resp_raw  <- resp_raw %>%
  filter(recorded_time <= ddate_bundle & recorded_time >= (adate_bundle - days(1)))

resp_final <- resp_raw %>%
  select(PATID = patid, bundle_id, adate_bundle, ddate_bundle,
         recorded_date           = recorded_time,
         flowsheet_name,
         flow_measure_name       = row_fs,
         flow_measure_result     = meas_value,
         flow_measure_result_units)

saveRDS(resp_final, l2('MITH_resp_support.rds'))


# =============================================================================
# SECTION 12: Inpatient Medications
# =============================================================================
enc        <- readRDS(l2('MITH_medicine_encounters.rds'))
enc_merge  <- as.data.table(enc)[, .(patid, bundle_id, adate_bundle, ddate_bundle)]

# Load medication name lookup
med_names      <- fread(sprintf('%sMedication Names Long.csv', look.path))
med_names_merge <- unique(med_names[, .(drug_category, parent_name_lower, brand_name_lower)])

# Build regex pattern per drug category
med_grepl_list <- sapply(unique(med_names$drug_category), function(cat) {
  parent <- paste(unique(med_names_merge[drug_category == cat]$parent_name_lower), collapse = '|')
  brand  <- paste(unique(med_names_merge[drug_category == cat]$brand_name_lower),  collapse = '|')
  paste(c(parent, brand), collapse = '|')
}, simplify = FALSE, USE.NAMES = TRUE)

ip_files <- paste0('Inpatient_Medications', 1:4, '.rds')

# Process each medication file: filter to cohort, align to encounter windows,
# extract drug category using regex patterns
ip_med_list <- lapply(seq_along(ip_files), function(x) {
  message(sprintf('Processing file %s', x))
  temp <- readRDS(l1(ip_files[[x]]))
  temp <- temp[PATID %in% unique(enc_merge$patid)]
  temp[, inpatient_med_administered_date := as.Date(inpatient_med_administered_date)]
  enc_merge[, adate_bundle := as.Date(adate_bundle)]
  enc_merge[, ddate_bundle := as.Date(ddate_bundle)]
  temp[, med_name_lower := tolower(inpatient_med_order_name)]
  temp <- merge(temp[, bundle_id := NULL], enc_merge, by.x = 'PATID', by.y = 'patid', allow.cartesian = TRUE)
  temp <- temp[between(inpatient_med_administered_date, adate_bundle, ddate_bundle)]

  lapply(names(med_grepl_list), function(cat) {
    temp[, med_name_simple := stri_extract_first_regex(med_name_lower, med_grepl_list[[cat]])]
    matched <- temp[!is.na(med_name_simple)]
    matched[, med_category := cat]
    matched
  }) %>% setNames(names(med_grepl_list))
})

# Combine across file chunks and save one file per drug category
med_df_list <- lapply(names(med_grepl_list), function(cat) {
  message(cat)
  df <- rbindlist(lapply(seq_along(ip_files), function(i) ip_med_list[[i]][[cat]]))
  df <- df[, .(patid = PATID, bundle_id, adate_bundle, ddate_bundle,
               inpatient_med_order_name,
               inpatient_med_administered_dose, inpatient_med_administered_dose_units,
               inpatient_med_administered_volume, inpatient_med_administered_volume_units,
               inpatient_med_administered_action, inpatient_med_administered_date,
               med_name_simple, med_category)]
  saveRDS(df, l2(sprintf('MITH_Inpatient_%s.rds', cat)))
  df
})

med_df_all <- rbindlist(med_df_list)
saveRDS(med_df_all, l2('MITH_Inpatient_Medications_all.rds'))


# =============================================================================
# SECTION 13: ICU Stays
# =============================================================================
enc       <- readRDS(l2('MITH_medicine_encounters.rds'))
proc_l2   <- readRDS(l2('MITH_procedure_info.rds'))

# Critical care CPT codes: 99291 (first 30-74 min) and 99292 (additional 30 min)
icu <- proc_l2[procedure_code %in% c('99291', '99292')][order(bundle_id, procedure_date)]
icu <- unique(icu[, .(bundle_id, procedure_date)])

# Build consecutive ICU stay groups
icu[, first_icu  := min(procedure_date), .(bundle_id)]
icu[, next_icu   := lead(procedure_date), .(bundle_id)]
icu[, next_dif   := as.numeric(difftime(next_icu, procedure_date, units = 'days'))]
icu[, last_dif   := as.numeric(difftime(procedure_date, lag(procedure_date), units = 'days')), .(bundle_id)]

icu[, icu_group  := ifelse((next_dif == 1 & last_dif > 1) | is.na(last_dif) |
                             (next_dif > 1 & last_dif > 1), 1, 0)]
icu[, icu_stay   := cumsum(icu_group), .(bundle_id)]
icu[, last_icu_stay := lag(icu_stay), .(bundle_id)]
icu[is.na(icu_stay), icu_stay := last_icu_stay + 1]

icu[, `:=`(icu_start = min(procedure_date),
           icu_end   = max(procedure_date)), .(bundle_id, icu_stay)]

icu_keep <- unique(icu[, .(bundle_id, icu_start, icu_end, first_icu)])
icu_keep <- unique(merge(icu_keep, enc[, .(bundle_id, adate_bundle, ddate_bundle)], by = 'bundle_id'))

# Flag whether ICU stay started on the day of admission
icu_keep[, icu_admission := ifelse(date(icu_start) == date(adate_bundle), 1, 0)]

message(sprintf('ICU admissions: %d of %d bundles (%.1f%%)',
                icu_keep[icu_admission == 1, uniqueN(bundle_id)],
                enc[, uniqueN(bundle_id)],
                100 * icu_keep[icu_admission == 1, uniqueN(bundle_id)] / enc[, uniqueN(bundle_id)]))

saveRDS(icu_keep, l2('MITH_ICU_stays.rds'))
