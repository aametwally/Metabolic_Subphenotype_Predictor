# Libraries needed for this script (unused libraries removed)
library(dplyr)        # Data manipulation
library(ggplot2)      # Plotting
library(MESS)         # AUC calculation

# Set working directories and paths for data and output
BASE_DIR = "PATH-to-BASE-DIRECTORY"
setwd(BASE_DIR)
DATA_DIR <- file.path(BASE_DIR, "data")
OUTPUT_FIGURE_DIR <- file.path(BASE_DIR, "figures")

#########################
###### Loading all data
#########################
## laod demographics
demographics <- read.csv(file.path(DATA_DIR, "all_cohort_metabolicsubphenotyping_metadata_10072024.csv"))
demographics = demographics[, c("subject_id", "exp_type", "age", "bmi", "sex", "ethnicity", "t2d_pred_history_binary", "ogtt_2h", "a1c", "fpg", "sspg")]

# Load PRS and filter by SubjectID in demographics
prs <- read.csv(file.path(DATA_DIR, "t2d_prs.csv"))
prs$SubjectID = gsub("43883-0","S", prs$SubjectID)

# Load DI
di <- read.csv(file.path(DATA_DIR, "disposition_index_calculations.csv"))

# Load IE
ie <- read.csv(file.path(DATA_DIR, "ie_cpeptide_based_calculations.csv"))

# Load Hepatic IR
hepatic_ir <- read.csv(file.path(DATA_DIR, "hepatic_ir_calculations.csv"))

# Load HOMA and Matsuda
homa_and_matsuda <- read.csv(file.path(DATA_DIR, "homa_and_matsuda_calculations.csv"))


#############################################
######## Aggregate Metabolic indicator ######
#############################################
# Split demographics based on exp_type
demographics_main_cohort <- demographics %>% filter(exp_type == "venous_without_matching_cgm_and_without_planned_athome_cgm")
demographics_not_main <- demographics %>% filter(exp_type != "venous_without_matching_cgm_and_without_planned_athome_cgm")

aggregated_metabolic_indicators = merge(demographics_main_cohort, prs, by.x = "subject_id", by.y = "SubjectID", ,all.x = TRUE)
aggregated_metabolic_indicators = merge(aggregated_metabolic_indicators, di, by.x = "subject_id", by.y = "SubjectID", all.x = TRUE)
aggregated_metabolic_indicators = merge(aggregated_metabolic_indicators, ie, by.x = "subject_id", by.y = "SubjectID", all.x = TRUE)
aggregated_metabolic_indicators = merge(aggregated_metabolic_indicators, hepatic_ir, by.x = "subject_id", by.y = "subject_id", all.x = TRUE)
aggregated_metabolic_indicators = merge(aggregated_metabolic_indicators, homa_and_matsuda, by.x = "subject_id", by.y = "SubjectID", all.x = TRUE)

# Combine the two results back
aggregated_metabolic_indicators <- bind_rows(aggregated_metabolic_indicators, demographics_not_main)

# Round numbers
aggregated_metabolic_indicators$a1c = round(aggregated_metabolic_indicators$a1c, 3)
aggregated_metabolic_indicators$fpg = round(aggregated_metabolic_indicators$fpg, 0)
aggregated_metabolic_indicators$bmi = round(aggregated_metabolic_indicators$bmi, 1)
aggregated_metabolic_indicators$sspg = round(aggregated_metabolic_indicators$sspg, 0)
aggregated_metabolic_indicators$T2D_PRS = round(aggregated_metabolic_indicators$T2D_PRS, 2)
aggregated_metabolic_indicators$di = round(aggregated_metabolic_indicators$di, 3)
aggregated_metabolic_indicators$ie = round(aggregated_metabolic_indicators$ie, 3)
aggregated_metabolic_indicators$hepatic_ir = round(aggregated_metabolic_indicators$hepatic_ir, 3)
aggregated_metabolic_indicators$matsuda_index = round(aggregated_metabolic_indicators$matsuda_index, 2)


# raname columns to make them lowercase
colnames(aggregated_metabolic_indicators) = tolower(colnames(aggregated_metabolic_indicators))


## Thresholding on a1c/fpg/ogtt_2h and metabolic phenotypes
aggregated_metabolic_indicators <- aggregated_metabolic_indicators %>%
  mutate(
    glycemic_status_ogtt_2h = case_when(
      ogtt_2h < 140 ~ "Normoglycemic",
      ogtt_2h >= 140 & ogtt_2h < 200 ~ "PreDM",
      ogtt_2h >= 200 ~ "T2DM"
    ),
    glycemic_status_a1c = case_when(
      a1c < 5.7 ~ "Normoglycemic",
      a1c >= 5.7 & a1c < 6.5 ~ "PreDM",
      a1c >= 6.5 ~ "T2DM"
    ),
    glycemic_status_fpg = case_when(
      fpg < 100 ~ "Normoglycemic",
      fpg >= 100 & fpg < 126 ~ "PreDM",
      fpg >= 126 ~ "T2DM"
    )
  )

## Thresholding on metabolic phenotypes

# Thresholding for MuscleIR
aggregated_metabolic_indicators <- aggregated_metabolic_indicators %>%
  mutate(sspg_2_classes = case_when(
    sspg <= 120 ~ "IS",  # Insulin Sensitive
    sspg > 120 ~ "IR"    # Insulin Resistant
  ))

# Thresholding for BetaCell
aggregated_metabolic_indicators <- aggregated_metabolic_indicators %>%
  mutate(di_2_classes = case_when(
    di < 1.58 ~ "Dysfunction",  
    di >= 1.58 ~ "Normal"    
  ))

# Thresholding for IncretInEffect
aggregated_metabolic_indicators <- aggregated_metabolic_indicators %>%
  mutate(ie_2_classes = case_when(
    ie < 53.38 ~ "Dysfunction",  
    ie >= 53.38 ~ "Normal"    
  ))


# Thresholding for HepaticIR
aggregated_metabolic_indicators <- aggregated_metabolic_indicators %>%
  mutate(hepatic_ir_2_classes = case_when(
    hepatic_ir < 4.35 ~ "IS",  # Insulin Sensitive
    hepatic_ir >= 4.35 ~ "IR"    # Insulin Resistant
  ))


#### Save aggregated metabolic indicators
write.csv(aggregated_metabolic_indicators, file.path(DATA_DIR, "aggregated_metabolic_indicators.csv"), row.names = FALSE)
write.csv(aggregated_metabolic_indicators[,c("subject_id", "exp_type", "age", 
                                             "bmi", "sex", "ethnicity", "t2d_pred_history_binary", 
                                             "ogtt_2h", "a1c", "fpg", "sspg")], 
          file.path(DATA_DIR, "all_cohort_metabolicsubphenotyping_metadata_10082024.csv"), row.names = FALSE)


