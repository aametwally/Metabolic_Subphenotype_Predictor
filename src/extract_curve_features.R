library(corrr)
library(corrplot)
library(dplyr)
library(ggdendro)
library(dendextend)
library(grid)
library(goeveg)
library(glmnet)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(ggnetwork)
library(ggparty)
library(Hmisc)
library(hrbrthemes)
library(imputeTS)
library(intergraph)
library(jtools)
library(MESS)
library(Metrics)
library(naniar)
library(pheatmap)
library(plyr)
library(PerformanceAnalytics)
library(RColorBrewer)
library(scales)
library(viridis)
library(zoo)
library(reshape2)
library(data.table)

setwd("/Users/ahmedm/Dropbox/Diabetes_v4_resubmission/Metabolic_Subphenotype_Predictor_07102022/")

## Read glucose values from different OGTT experiment setup (Venous/CGM, CTRU/Home)
# ogtt_glucose_all_melted = read.csv(file = "data/ogtt_glucose_all_cohorts_09072023.csv")
ogtt_glucose_all_melted = read.csv("data/all_cohort_metabolicsubphenotyping_ogtt_glucose_09102023.csv")

head(ogtt_glucose_all_melted)
dim(ogtt_glucose_all_melted)


### Extract dataframes for each data type
dcast_ogtt_timeseries <- function(df = NULL){
  dt = reshape2::dcast(df[,c("subject_id", "timepoint", "glucose")], timepoint ~ subject_id)
  rownames(dt) = dt$timepoint
  dt = dt[,-1]
  return(dt)
}

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


######## Preprocess OGTT glucose timeseries
## Impute missing timepoints
ogtt_glucose_venous_exp_without_matching_imputed = as.data.frame(apply(ogtt_glucose_venous_exp_without_matching, 2, na_interpolation))
ogtt_glucose_venous_exp_with_matching_imputed = as.data.frame(apply(ogtt_glucose_venous_exp_with_matching, 2, na_interpolation)) 
ogtt_glucose_ctru_cgm_imputed = as.data.frame(apply(ogtt_glucose_ctru_cgm, 2, na_interpolation))
ogtt_glucose_home_cgm_1_imputed = as.data.frame(apply(ogtt_glucose_home_cgm_1, 2, na_interpolation))
ogtt_glucose_home_cgm_2_imputed = as.data.frame(apply(ogtt_glucose_home_cgm_2, 2, na_interpolation))


## Normalize imputed timeseries 
znorm <- function(ts){
  ts.mean <- mean(ts)
  ts.dev <- sd(ts)
  (ts - ts.mean)/ts.dev
}
ogtt_glucose_venous_exp_without_matching_imputed_normalized = apply(ogtt_glucose_venous_exp_without_matching_imputed, 2, znorm)
ogtt_glucose_venous_exp_with_matching_imputed_normalized = apply(ogtt_glucose_venous_exp_with_matching_imputed, 2, znorm)
ogtt_glucose_ctru_cgm_imputed_normalized = apply(ogtt_glucose_ctru_cgm_imputed, 2, znorm)
ogtt_glucose_home_cgm_1_imputed_normalized = apply(ogtt_glucose_home_cgm_1_imputed, 2, znorm)
ogtt_glucose_home_cgm_2_imputed_normalized = apply(ogtt_glucose_home_cgm_2_imputed, 2, znorm)


## Smooth normalized and imputed OGTT timeseries
smooth_timeseries = function(x){
  smoothingSpline = smooth.spline(x, spar=0.35)
  names(smoothingSpline$y) = names(x)
  smoothingSpline$y
}
ogtt_glucose_venous_exp_without_matching_imputed_normalized_smoothed = apply(ogtt_glucose_venous_exp_without_matching_imputed_normalized, 2, smooth_timeseries)
ogtt_glucose_venous_exp_with_matching_imputed_normalized_smoothed = apply(ogtt_glucose_venous_exp_with_matching_imputed_normalized, 2, smooth_timeseries)
ogtt_glucose_ctru_cgm_imputed_normalized_smoothed = apply(ogtt_glucose_ctru_cgm_imputed_normalized, 2, smooth_timeseries)
ogtt_glucose_home_cgm_1_imputed_normalized_smoothed = apply(ogtt_glucose_home_cgm_1_imputed_normalized, 2, smooth_timeseries)
ogtt_glucose_home_cgm_2_imputed_normalized_smoothed = apply(ogtt_glucose_home_cgm_2_imputed_normalized, 2, smooth_timeseries)


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
  return(ogtt_features)
}


ogtt_glucose_venous_exp_without_matching_features = extract_curve_features(ogtt_imputed = ogtt_glucose_venous_exp_without_matching_imputed)
ogtt_glucose_venous_exp_without_matching_features$subject_id = rownames(ogtt_glucose_venous_exp_without_matching_features)
ogtt_glucose_venous_exp_without_matching_features = relocate(ogtt_glucose_venous_exp_without_matching_features, column_name = 'subject_id', position = 1)
rownames(ogtt_glucose_venous_exp_without_matching_features) = NULL
write.csv(ogtt_glucose_venous_exp_without_matching_features, 
          file = "data/metabolic_subphenotyping_ogtt_glucose_venous_exp_without_matching_features_09102023.csv")

ogtt_glucose_venous_exp_with_matching_features = extract_curve_features(ogtt_imputed = ogtt_glucose_venous_exp_with_matching_imputed)
ogtt_glucose_venous_exp_with_matching_features$subject_id = rownames(ogtt_glucose_venous_exp_with_matching_features)
ogtt_glucose_venous_exp_with_matching_features = relocate(ogtt_glucose_venous_exp_with_matching_features, column_name = 'subject_id', position = 1)
rownames(ogtt_glucose_venous_exp_with_matching_features) = NULL
write.csv(ogtt_glucose_venous_exp_with_matching_features, 
          file = "data/metabolic_subphenotyping_ogtt_glucose_venous_exp_with_matching_features_09102023.csv")


ogtt_glucose_ctru_cgm_features = extract_curve_features(ogtt_imputed = ogtt_glucose_ctru_cgm_imputed)
ogtt_glucose_ctru_cgm_features$subject_id = rownames(ogtt_glucose_ctru_cgm_features)
ogtt_glucose_ctru_cgm_features = relocate(ogtt_glucose_ctru_cgm_features, column_name = 'subject_id', position = 1)
rownames(ogtt_glucose_ctru_cgm_features) = NULL
write.csv(ogtt_glucose_ctru_cgm_features, 
          file = "data/metabolic_subphenotyping_ogtt_glucose_ctru_cgm_features_09102023.csv")


ogtt_glucose_home_cgm_1_features = extract_curve_features(ogtt_imputed = ogtt_glucose_home_cgm_1_imputed)
ogtt_glucose_home_cgm_1_features$subject_id = rownames(ogtt_glucose_home_cgm_1_features)
ogtt_glucose_home_cgm_1_features = relocate(ogtt_glucose_home_cgm_1_features, column_name = 'subject_id', position = 1)
rownames(ogtt_glucose_home_cgm_1_features) = NULL
write.csv(ogtt_glucose_home_cgm_1_features, 
          file = "data/metabolic_subphenotyping_ogtt_glucose_home_cgm_1_matching_features_09102023.csv")

ogtt_glucose_home_cgm_2_matching_features = extract_curve_features(ogtt_imputed = ogtt_glucose_home_cgm_2_imputed)
ogtt_glucose_home_cgm_2_matching_features$subject_id = rownames(ogtt_glucose_home_cgm_2_matching_features)
ogtt_glucose_home_cgm_2_matching_features = relocate(ogtt_glucose_home_cgm_2_matching_features, column_name = 'subject_id', position = 1)
rownames(ogtt_glucose_home_cgm_2_matching_features) = NULL
write.csv(ogtt_glucose_home_cgm_2_matching_features, 
          file = "data/metabolic_subphenotyping_ogtt_glucose_home_cgm_2_matching_features_09102023.csv")





