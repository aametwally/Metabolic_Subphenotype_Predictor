# Required Libraries
library(dplyr)
library(reshape2)
library(MESS)
library(imputeTS)
library(ggplot2)
library(goeveg)


# Set working directories and paths for data and output
BASE_DIR = "PATH-to-BASE-DIRECTORY"
setwd(BASE_DIR)
DATA_DIR <- file.path(BASE_DIR, "data")
OUTPUT_FIGURE_DIR <- file.path(BASE_DIR, "figures")


# Read OGTT glucose data
ogtt_glucose_all_melted <- read.csv(file.path(DATA_DIR, "all_cohort_metabolicsubphenotyping_ogtt_glucose_09102023.csv"))

# Check data structure
head(ogtt_glucose_all_melted)
dim(ogtt_glucose_all_melted)

#' Cast OGTT data to wide format for each subject
#' @param df Dataframe containing OGTT glucose values
#' @return Dataframe with timepoints as rows and subjects as columns
dcast_ogtt_timeseries <- function(df = NULL) {
  dt <- reshape2::dcast(df[, c("subject_id", "timepoint", "glucose")], timepoint ~ subject_id)
  rownames(dt) <- dt$timepoint
  dt <- dt[, -1]
  return(dt)
}


# Extract glucose timeseries for various experimental setups
# Venous without matching cgm
ogtt_glucose_venous_exp_without_matching = subset(ogtt_glucose_all_melted, sample_location_extraction_method == "CTRU_Venous" & exp_type == "venous_without_matching_cgm_and_without_planned_athome_cgm")
ogtt_glucose_venous_exp_without_matching = dcast_ogtt_timeseries(df = ogtt_glucose_venous_exp_without_matching)
dim(ogtt_glucose_venous_exp_without_matching)

# Venous with matching cgm
ogtt_glucose_venous_exp_with_matching = subset(ogtt_glucose_all_melted, sample_location_extraction_method == "CTRU_Venous" & exp_type == "venous_with_matching_cgm_and_with_planned_athome_cgm")
ogtt_glucose_venous_exp_with_matching = dcast_ogtt_timeseries(df = ogtt_glucose_venous_exp_with_matching)
dim(ogtt_glucose_venous_exp_with_matching)

# CTRU CGM
ogtt_glucose_ctru_cgm = subset(ogtt_glucose_all_melted, sample_location_extraction_method == "CTRU_CGM")
ogtt_glucose_ctru_cgm = dcast_ogtt_timeseries(df = ogtt_glucose_ctru_cgm)
dim(ogtt_glucose_ctru_cgm)

# Home CGM 1
ogtt_glucose_home_cgm_1 = subset(ogtt_glucose_all_melted, sample_location_extraction_method == "Home_CGM_1")
ogtt_glucose_home_cgm_1 = dcast_ogtt_timeseries(df = ogtt_glucose_home_cgm_1)
dim(ogtt_glucose_home_cgm_1)

# Home CGM 2
ogtt_glucose_home_cgm_2 = subset(ogtt_glucose_all_melted, sample_location_extraction_method == "Home_CGM_2")
ogtt_glucose_home_cgm_2 = dcast_ogtt_timeseries(df = ogtt_glucose_home_cgm_2)
dim(ogtt_glucose_home_cgm_2)



# Impute missing timepoints
impute_glucose_data <- function(df) {
  return(as.data.frame(apply(df, 2, na_interpolation)))
}

ogtt_glucose_venous_exp_without_matching_imputed <- impute_glucose_data(ogtt_glucose_venous_exp_without_matching)
ogtt_glucose_venous_exp_with_matching_imputed <- impute_glucose_data(ogtt_glucose_venous_exp_with_matching)
ogtt_glucose_ctru_cgm_imputed <- impute_glucose_data(ogtt_glucose_ctru_cgm)
ogtt_glucose_home_cgm_1_imputed <- impute_glucose_data(ogtt_glucose_home_cgm_1)
ogtt_glucose_home_cgm_2_imputed <- impute_glucose_data(ogtt_glucose_home_cgm_2)

# Save features as CSV
write.csv(ogtt_glucose_venous_exp_without_matching_imputed, file = file.path(DATA_DIR, "ogtt_glucose_venous_exp_without_matching_imputed.csv"))
write.csv(ogtt_glucose_venous_exp_with_matching_imputed, file = file.path(DATA_DIR, "ogtt_glucose_venous_exp_with_matching_imputed.csv"))
write.csv(ogtt_glucose_ctru_cgm_imputed, file = file.path(DATA_DIR, "ogtt_glucose_ctru_cgm_imputed.csv"))
write.csv(ogtt_glucose_home_cgm_1_imputed, file = file.path(DATA_DIR, "ogtt_glucose_home_cgm_1_imputed.csv"))
write.csv(ogtt_glucose_home_cgm_2_imputed, file = file.path(DATA_DIR, "ogtt_glucose_home_cgm_2_imputed.csv"))


# Normalize timeseries
znorm <- function(ts) {
  ts.mean <- mean(ts)
  ts.dev <- sd(ts)
  (ts - ts.mean) / ts.dev
}

normalize_timeseries <- function(df) {
  apply(df, 2, znorm)
}

# Applying normalization to each dataset
ogtt_glucose_venous_exp_without_matching_imputed_normalized <- normalize_timeseries(ogtt_glucose_venous_exp_without_matching_imputed)
ogtt_glucose_venous_exp_with_matching_imputed_normalized <- normalize_timeseries(ogtt_glucose_venous_exp_with_matching_imputed)
ogtt_glucose_ctru_cgm_imputed_normalized <- normalize_timeseries(ogtt_glucose_ctru_cgm_imputed)
ogtt_glucose_home_cgm_1_imputed_normalized <- normalize_timeseries(ogtt_glucose_home_cgm_1_imputed)
ogtt_glucose_home_cgm_2_imputed_normalized <- normalize_timeseries(ogtt_glucose_home_cgm_2_imputed)

# Save features as CSV
write.csv(ogtt_glucose_venous_exp_without_matching_imputed_normalized, file = file.path(DATA_DIR, "ogtt_glucose_venous_exp_without_matching_imputed_normalized.csv"))
write.csv(ogtt_glucose_venous_exp_with_matching_imputed_normalized, file = file.path(DATA_DIR, "ogtt_glucose_venous_exp_with_matching_imputed_normalized.csv"))
write.csv(ogtt_glucose_ctru_cgm_imputed_normalized, file = file.path(DATA_DIR, "ogtt_glucose_ctru_cgm_imputed_normalized.csv"))
write.csv(ogtt_glucose_home_cgm_1_imputed_normalized, file = file.path(DATA_DIR, "ogtt_glucose_home_cgm_1_imputed_normalized.csv"))
write.csv(ogtt_glucose_home_cgm_2_imputed_normalized, file = file.path(DATA_DIR, "ogtt_glucose_home_cgm_2_imputed_normalized.csv"))



# Smooth timeseries using spline smoothing
smooth_timeseries <- function(x) {
  smoothingSpline <- smooth.spline(x, spar = 0.35)
  return(smoothingSpline$y)
}

# Apply smoothing
apply_smoothing <- function(df) {
  smoothed_data <- apply(df, 2, smooth_timeseries)
  return(smoothed_data)
}

# Applying smoothing to each dataset
ogtt_glucose_venous_exp_without_matching_imputed_normalized_smoothed <- apply_smoothing(ogtt_glucose_venous_exp_without_matching_imputed_normalized)
ogtt_glucose_venous_exp_with_matching_imputed_normalized_smoothed <- apply_smoothing(ogtt_glucose_venous_exp_with_matching_imputed_normalized)
ogtt_glucose_ctru_cgm_imputed_normalized_smoothed <- apply_smoothing(ogtt_glucose_ctru_cgm_imputed_normalized)
ogtt_glucose_home_cgm_1_imputed_normalized_smoothed <- apply_smoothing(ogtt_glucose_home_cgm_1_imputed_normalized)
ogtt_glucose_home_cgm_2_imputed_normalized_smoothed <- apply_smoothing(ogtt_glucose_home_cgm_2_imputed_normalized)


# Save features as CSV
write.csv(ogtt_glucose_venous_exp_without_matching_imputed_normalized_smoothed, file = file.path(DATA_DIR, "ogtt_glucose_venous_exp_without_matching_imputed_normalized_smoothed.csv"))
write.csv(ogtt_glucose_venous_exp_with_matching_imputed_normalized_smoothed, file = file.path(DATA_DIR, "ogtt_glucose_venous_exp_with_matching_imputed_normalized_smoothed.csv"))
write.csv(ogtt_glucose_ctru_cgm_imputed_normalized_smoothed, file = file.path(DATA_DIR, "ogtt_glucose_ctru_cgm_imputed_normalized_smoothed.csv"))
write.csv(ogtt_glucose_home_cgm_1_imputed_normalized_smoothed, file = file.path(DATA_DIR, "ogtt_glucose_home_cgm_1_imputed_normalized_smoothed.csv"))
write.csv(ogtt_glucose_home_cgm_2_imputed_normalized_smoothed, file = file.path(DATA_DIR, "ogtt_glucose_home_cgm_2_imputed_normalized_smoothed.csv"))


### Feature Extraction from the OGTT curves
extract_curve_features = function(ogtt_imputed = NULL){
  ## remove timepoint -10 from imputed ogtt
  ogtt_imputed = ogtt_imputed[-1, ]

  ###### Total AUC (tAUC)
  ogtt_auc = apply(ogtt_imputed, 2, function(x){
    MESS::auc(as.numeric(rownames(ogtt_imputed)), x, from = 0, to = 180,  type = "linear")
  })

  ###### Incremental AUC (iAUC).
  ogtt_iauc = unlist(ogtt_auc - 180*ogtt_imputed["0", ])

  ###### pAUC: positive AUC
  qq = apply(ogtt_imputed, 2, function(x){
    x-x["0"]
  })
  qq[which(qq<0)]=0
  ogtt_pauc = apply(qq, 2, function(x){
    MESS::auc(as.numeric(rownames(qq)), x, from = 0, to = 180,  type = "linear")
  })

  ###### nAUC: negative AUC
  ogtt_nauc = ogtt_pauc - ogtt_iauc

  ###### OGTT Peak
  ogtt_max = apply(ogtt_imputed, 2, function(x){
    max(x)
  })

  ###### FPG (Glucose at 0)
  ogtt_fpg = unlist(ogtt_imputed["0", ])

  ###### Glucose at 60
  ogtt_60 = unlist(ogtt_imputed["60", ])

  ###### Glucose at 120
  ogtt_120 = unlist(ogtt_imputed["120", ])

  ###### Glucose at 180
  ogtt_180 = unlist(ogtt_imputed["180", ])

  ###### CV
  ogtt_cv = unlist(apply(ogtt_imputed, 2, cv))

  ###### Curve size
  ogtt_curve_size = apply(diff(as.matrix(ogtt_imputed)), 2, function(x){
    sum(abs(x))
  })

  ###### Time to peak
  time = as.numeric(rownames(ogtt_imputed))
  ogtt_time_baseline_peak = apply(ogtt_imputed, 2, function(x){
    mx = max(x)
    idx = which(x == mx)
    t = time[idx[1]]
  })
  ogtt_time_baseline_peak = unlist(ogtt_time_baseline_peak)

  ###### Time it takes to arrive to first sample below baseline after peak
  time = as.numeric(rownames(ogtt_imputed))
  ogtt_time_peak_baseline = apply(ogtt_imputed, 2, function(x){
    mx = max(x)
    idx_max = which(x == mx)
    idx_less_baseline = which(x<x[1])
    idx_less_baseline_more_max = which(idx_less_baseline > idx_max[1])
    idx = idx_less_baseline[idx_less_baseline_more_max[1]]
    t_peak_baseline = time[idx[1]] - time[idx_max[1]]
  })
  ogtt_time_peak_baseline = unlist(ogtt_time_peak_baseline)

  ###### Slope from baseline to peak
  time = as.numeric(rownames(ogtt_imputed))
  ogtt_slope_baseline_peak = apply(ogtt_imputed, 2, function(x){
    mx = max(x)
    idx = which(x == mx)
    slope = as.numeric((x[idx[1]]-x[1])/(time[idx[1]]-time[1]))
  })
  ogtt_slope_baseline_peak = unlist(ogtt_slope_baseline_peak)

  ###### Slope from peak to last time point
  time = as.numeric(rownames(ogtt_imputed))
  ogtt_slope_peak_last = apply(ogtt_imputed, 2, function(x){
    mx = max(x)
    idx = which(x == mx)
    slope = as.numeric(x[length(time)] - x[idx[1]])/(time[length(time)] - time[idx[1]])
  })
  ogtt_slope_peak_last = unlist(ogtt_slope_peak_last)

  ###### Mark subjects that went below baseline
  ogtt_time_below_basline = apply(ogtt_imputed, 2, function(x){
    mx = max(x)
    idx_max = which(x == mx)
    idx_less_baseline = which(x<x[1])
    idx_less_baseline_more_max = which(idx_less_baseline > idx_max[1])
    length(idx_less_baseline_more_max)>0
  })

  ogtt_features = data.frame(ogtt_fpg = ogtt_fpg, ogtt_60 = ogtt_60, ogtt_120 = ogtt_120,
                             ogtt_180 = (ogtt_180), ogtt_auc = ogtt_auc, ogtt_iauc = (ogtt_iauc),
                             ogtt_pauc = (ogtt_pauc), ogtt_nauc = (ogtt_nauc), ogtt_max = ogtt_max,
                             ogtt_curve_size = ogtt_curve_size, ogtt_cv = ogtt_cv,
                             ogtt_time_baseline_peak = ogtt_time_baseline_peak, ogtt_time_peak_baseline = ogtt_time_peak_baseline,
                             ogtt_slope_baseline_peak = ogtt_slope_baseline_peak, ogtt_slope_peak_last = ogtt_slope_peak_last,
                             ogtt_time_below_basline = ogtt_time_below_basline)
  
  
  ogtt_features$subject_id = rownames(ogtt_features)
  ogtt_features <- ogtt_features %>%
    select(subject_id, everything())
  rownames(ogtt_features) = NULL
  return(ogtt_features)
}

# Extract features for each dataset 
ogtt_glucose_venous_exp_without_matching_features <- extract_curve_features(ogtt_glucose_venous_exp_without_matching_imputed)
ogtt_glucose_venous_exp_with_matching_features <- extract_curve_features(ogtt_glucose_venous_exp_with_matching_imputed)
ogtt_glucose_ctru_cgm_features <- extract_curve_features(ogtt_glucose_ctru_cgm_imputed)
ogtt_glucose_home_cgm_1_features <- extract_curve_features(ogtt_glucose_home_cgm_1_imputed)
ogtt_glucose_home_cgm_2_features <- extract_curve_features(ogtt_glucose_home_cgm_2_imputed)

# Save features as CSV
write.csv(ogtt_glucose_venous_exp_without_matching_features, file = file.path(DATA_DIR, "ogtt_glucose_venous_exp_without_matching_features.csv"))
write.csv(ogtt_glucose_venous_exp_with_matching_features, file = file.path(DATA_DIR, "ogtt_glucose_venous_exp_with_matching_features.csv"))
write.csv(ogtt_glucose_ctru_cgm_features, file = file.path(DATA_DIR, "ogtt_glucose_ctru_cgm_features.csv"))
write.csv(ogtt_glucose_home_cgm_1_features, file = file.path(DATA_DIR, "ogtt_glucose_home_cgm_1_features.csv"))
write.csv(ogtt_glucose_home_cgm_2_features, file = file.path(DATA_DIR, "ogtt_glucose_home_cgm_2_features.csv"))
