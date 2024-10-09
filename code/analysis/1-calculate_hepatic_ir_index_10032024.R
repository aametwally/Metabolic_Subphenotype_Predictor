# Libraries needed for this script
library(dplyr)        # Data manipulation
library(ggplot2)      # Plotting
library(reshape2)     # Reshaping data
library(MESS)         # AUC calculation
library(data.table)   # Fast reading of large data
library(plyr)

# Set working directories and paths for data and output
BASE_DIR = "PATH-to-BASE-DIRECTORY"
setwd(BASE_DIR)
DATA_DIR <- file.path(BASE_DIR, "data")
OUTPUT_FIGURE_DIR <- file.path(BASE_DIR, "figures")

#########################
###### Loading all data
#########################

# Load metadata
metadata <- read.csv(file.path(DATA_DIR, "all_cohort_metabolicsubphenotyping_metadata_10072024.csv"))

# Filter metadata for the initial cohort
metadata_initial_cohort <- metadata %>%
  filter(exp_type == "venous_without_matching_cgm_and_without_planned_athome_cgm")

# Reorder metadata based on A1C levels in decreasing order
metabolic_indicators_orderd <- metadata_initial_cohort %>%
  arrange(desc(a1c))

# Lock the order of subject IDs based on the reordered data
metabolic_indicators_orderd$subject_id <- factor(metabolic_indicators_orderd$subject_id,
                                                 levels = metabolic_indicators_orderd$subject_id)

demographics = metadata_initial_cohort
a1c_order <- metabolic_indicators_orderd %>% pull(subject_id)


### Load and clean combined analytes
washu_combined_batches <- read.csv(file.path(DATA_DIR, "washu_combined_batches_10272020.csv")) %>%
  arrange(SubjectID, Assay, TimePoint) %>%
  mutate(SubjectID = gsub("43883-0", "S", SubjectID)) %>%
  filter(SubjectID %in% demographics$subject_id)

# Subsetting analytes based on assay type
ogtt_analytes <- subset(washu_combined_batches, Assay == "OGTT")
iigi_analytes <- subset(washu_combined_batches, Assay == "IIGI")

### Load OGTT and IIGI glucose measurements
ogtt <- fread(file.path(DATA_DIR, "ogtt.tsv")) 
colnames(ogtt)[which(colnames(ogtt) == "userID")] = "SubjectID"
ogtt <- ogtt %>% # rename(SubjectID = userID) %>%
  mutate(SubjectID = gsub("43883-0", "S", SubjectID)) %>%
  filter(SubjectID %in% demographics$subject_id) %>%
  mutate(ogtt_value = as.numeric(ogtt_value))


iigi <- fread(file.path(DATA_DIR, "iigi.tsv"))
colnames(iigi)[which(colnames(iigi) == "userID")] = "SubjectID"
iigi <- iigi %>% 
  mutate(SubjectID = gsub("43883-0", "S", SubjectID)) %>%
  filter(SubjectID %in% demographics$subject_id) %>%
  mutate(iigi_value = as.numeric(iigi_value))



### Load and clean metabolic panel data
metabolic_df <- fread(file.path(DATA_DIR, "blood.tsv"), sep = "\t") %>%
  na.omit() %>%
  filter(value != "Exact value cannot be determined") %>%
  mutate(SubjectID = gsub("43883-0", "S", userID)) %>%
  filter(SubjectID %in% demographics$subject_id)



#########################
###### Hepatic IR Calculations
#########################
# Merge and clean insulin data from OGTT and IIGI
ogtt_insulin <- ogtt_analytes %>%
  select(SubjectID, TimePoint, Insulin) %>%
  na.omit() %>%
  mutate(test = "OGTT")

### Calculate AUC for Insulin OGTT (0-120 minutes)
ogtt_insulin_wide <- dcast(ogtt_insulin, SubjectID ~ TimePoint, value.var = "Insulin")
ogtt_insulin_wide <- apply(ogtt_insulin_wide[, -1], 1, na_interpolation)


### Calculate AUC for Insulin OGTT 0-120
ogtt_insulin_dcast = reshape2::dcast(ogtt_insulin, SubjectID ~ TimePoint, value.var = "Insulin")
rownames(ogtt_insulin_dcast) = ogtt_insulin_dcast$SubjectID
ogtt_insulin_dcast = ogtt_insulin_dcast[,-1]
ogtt_insulin_dcast_imputed = t(apply(ogtt_insulin_dcast, 1, na_interpolation))

# Save imputed insulin levels
write.csv(ogtt_insulin_dcast_imputed, file = file.path(DATA_DIR, "ogtt_insulin_dcast_imputed.csv"))

# Calculation of AUC insulin curve for 0-120 min
timepoints = c(0, 15, 30, 60, 90, 120)
insulin_ogtt_auc_0_120 = apply(ogtt_insulin_dcast_imputed[,c("0", "15", "30", "60", "90", "120")], 1, function(x){
  MESS::auc(timepoints, x)
})
insulin_ogtt_auc_0_120 = as.data.frame(insulin_ogtt_auc_0_120)

### Add body fat percentage to the data
demographics <- demographics %>%
  mutate(
    sex_01 = ifelse(sex == "M", 1, 0),
    bf_percentage = 1.2 * bmi + 0.23 * age - 10.8 * sex_01 - 5.4
  )

### Merge insulin AUC and demographic data
demographics_insulinAUC <- merge(demographics, as.data.frame(insulin_ogtt_auc_0_120), 
                                 by.x = "subject_id", by.y = "row.names", all.x = TRUE)


# Extract HDL 
hdl = subset(metabolic_df, blood == "hdl")
hdl = hdl[,c("SubjectID", "value")]
colnames(hdl) = c("SubjectID", "hdl")
hdl$hdl = as.numeric(as.character(hdl$hdl))
hdl = ddply(hdl, ~SubjectID, summarise, hdl = mean(hdl))


### Calculate Hepatic IR Index
hepatic_ir_df <- merge(demographics_insulinAUC, hdl, by.x = "subject_id", by.y = "SubjectID", all.x = TRUE) %>%
  mutate(
    hepatic_ir = -0.091 + 0.4 * log(insulin_ogtt_auc_0_120) +
      0.346 * log(bf_percentage) - 0.408 * log(hdl) +
      0.435 * log(bmi)
  )

# Export Hepatic IR results
write.csv(hepatic_ir_df[,c("subject_id", "hepatic_ir")], file = file.path(DATA_DIR, "hepatic_ir_calculations.csv"), row.names = FALSE)
write.csv(hepatic_ir_df, file = file.path(DATA_DIR, "hepatic_ir_demographics.csv"), row.names = FALSE)
