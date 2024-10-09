# Required Libraries
library(reshape2) # Used for reshaping data
library(MESS)     # Used for calculating AUC
library(ggplot2)  # Used for data visualization
library(dplyr)    # Replaces ddply and helps simplify data manipulation
library(data.table)
library(imputeTS)

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



############################
####### Incretin Effect Calculations
############################

# 1. OGTT AUC for C-Peptide
# Reshape OGTT C-peptide data for AUC calculation
ogtt_cpeptide <- ogtt_analytes %>% select(SubjectID, TimePoint, C_peptide) %>% na.omit()
dt <- reshape2::dcast(ogtt_cpeptide, TimePoint ~ SubjectID) %>% t()
colnames(dt) <- dt[1,]
dt <- dt[-1,] %>% as.data.frame()
dt_imputed_ogtt <- apply(dt, 1, na_interpolation)

# AUC Calculation for OGTT
timepoints <- as.numeric(rownames(dt_imputed_ogtt))
cpeptide_ogtt_auc <- apply(dt_imputed_ogtt, 2, function(x) MESS::auc(timepoints, x))

# 2. IIGI AUC for C-Peptide (Similar steps to OGTT)
iigi_cpeptide <- iigi_analytes %>% select(SubjectID, TimePoint, C_peptide) %>% na.omit()
dt <- reshape2::dcast(iigi_cpeptide, TimePoint ~ SubjectID) %>% t()
colnames(dt) <- dt[1,]
dt <- dt[-1,] %>% as.data.frame()
dt_imputed_iigi <- apply(dt, 1, na_interpolation)

# AUC Calculation for IIGI
timepoints <- as.numeric(rownames(dt_imputed_iigi))
cpeptide_iigi_auc <- apply(dt_imputed_iigi, 2, function(x) MESS::auc(timepoints, x))

# 3. Incretin Effect Calculation (C-peptide)
# IE = % difference between OGTT and IIGI AUCs
ie <- 100 * (cpeptide_ogtt_auc - cpeptide_iigi_auc) / cpeptide_ogtt_auc
ie = data.frame(SubjectID = names(ie), ie = unlist(ie))


# Save results to CSV
write.csv(ie, file = file.path(DATA_DIR, "ie_cpeptide_based_calculations.csv"), row.names = FALSE)



############################
####### Visualization for C-Peptide (OGTT vs IIGI)
############################

# Prepare C-peptide data for OGTT and IIGI
cpeptide_ogtt_melted <- reshape2::melt(as.matrix(dt_imputed_ogtt))
cpeptide_ogtt_melted$test <- "OGTT"
cpeptide_iigi_melted <- reshape2::melt(as.matrix(dt_imputed_iigi))
cpeptide_iigi_melted$test <- "IIGI"
cpeptide_combined <- rbind(cpeptide_ogtt_melted, cpeptide_iigi_melted)
colnames(cpeptide_combined) <- c("time", "SubjectID", "c_pepdite", "test")

# Filter and reorder data by A1C
cpeptide_combined <- cpeptide_combined %>%
  filter(SubjectID %in% demographics$subject_id) %>%
  mutate(SubjectID = factor(SubjectID, levels = a1c_order))

# Plot C-peptide concentrations over time (OGTT vs IIGI)
gp <- ggplot(cpeptide_combined, aes(time, c_pepdite, col = test)) + 
  geom_line(aes(group = test), alpha = 0.8, size = 1.5) + 
  scale_color_manual(values = c("#56B4E9", "#E69F00", "#999999")) + 
  facet_wrap(~SubjectID, nrow = 4) + 
  theme_bw() + 
  theme(legend.position = "top", 
        axis.text.x = element_text(size = 10, angle = 45, hjust = 0.5),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 15), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) + 
  labs(y = "C-Peptide concentration (ng/mL)", x = "Time (mins)", 
       title = "C-Peptide concentration (OGTT vs IIGI)")

# Save Plot
ggsave(file.path(OUTPUT_FIGURE_DIR, "cpeptide_OGTT_IIGI_ordered.pdf"), gp, height = 6, width = 10)

############################
####### GIP and GLP-1 Features Calculation
############################

# 1. AUC for GLP-1 (OGTT and IIGI)
# Reshape OGTT GLP-1 data
ogtt_glp1 <- ogtt_analytes %>% select(SubjectID, TimePoint, GLP1) %>% na.omit()
dt_glp1 <- reshape2::dcast(ogtt_glp1, TimePoint ~ SubjectID) %>% t()
colnames(dt_glp1) <- dt_glp1[1,]
dt_glp1 <- dt_glp1[-1,] %>% as.data.frame()
dt_imputed_ogtt_glp1 <- apply(dt_glp1, 1, na_interpolation)

# AUC Calculation for OGTT GLP-1
timepoints <- as.numeric(rownames(dt_imputed_ogtt_glp1))
glp1_ogtt_auc <- apply(dt_imputed_ogtt_glp1, 2, function(x) MESS::auc(timepoints, x))

# Repeat for IIGI GLP-1
iigi_glp1 <- iigi_analytes %>% select(SubjectID, TimePoint, GLP1) %>% na.omit()
dt_glp1 <- reshape2::dcast(iigi_glp1, TimePoint ~ SubjectID) %>% t()
colnames(dt_glp1) <- dt_glp1[1,]
dt_glp1 <- dt_glp1[-1,] %>% as.data.frame()
dt_imputed_iigi_glp1 <- apply(dt_glp1, 1, na_interpolation)

# AUC Calculation for IIGI GLP-1
timepoints <- as.numeric(rownames(dt_imputed_iigi_glp1))
glp1_iigi_auc <- apply(dt_imputed_iigi_glp1, 2, function(x) MESS::auc(timepoints, x))

# 2. AUC for GIP (OGTT and IIGI)
# Reshape OGTT GIP data
ogtt_gip <- ogtt_analytes %>% select(SubjectID, TimePoint, GIP) %>% na.omit()
dt_gip <- reshape2::dcast(ogtt_gip, TimePoint ~ SubjectID) %>% t()
colnames(dt_gip) <- dt_gip[1,]
dt_gip <- dt_gip[-1,] %>% as.data.frame()
dt_imputed_ogtt_gip <- apply(dt_gip, 1, na_interpolation)

# AUC Calculation for OGTT GIP
timepoints <- as.numeric(rownames(dt_imputed_ogtt_gip))
gip_ogtt_auc <- apply(dt_imputed_ogtt_gip, 2, function(x) MESS::auc(timepoints, x))

# Repeat for IIGI GIP
iigi_gip <- iigi_analytes %>% select(SubjectID, TimePoint, GIP) %>% na.omit()
dt_gip <- reshape2::dcast(iigi_gip, TimePoint ~ SubjectID) %>% t()
colnames(dt_gip) <- dt_gip[1,]
dt_gip <- dt_gip[-1,] %>% as.data.frame()
dt_imputed_iigi_gip <- apply(dt_gip, 1, na_interpolation)

# AUC Calculation for IIGI GIP
timepoints <- as.numeric(rownames(dt_imputed_iigi_gip))
gip_iigi_auc <- apply(dt_imputed_iigi_gip, 2, function(x) MESS::auc(timepoints, x))


# 1. GLP-1 Maximum Concentration (OGTT and IIGI)
# Max GLP-1 for each subject (OGTT)
glp1_ogtt_max <- apply(dt_imputed_ogtt_glp1, 2, max, na.rm = TRUE)

# Max GLP-1 for each subject (IIGI)
glp1_iigi_max <- apply(dt_imputed_iigi_glp1, 2, max, na.rm = TRUE)

# 2. GIP Maximum Concentration (OGTT and IIGI)
# Max GIP for each subject (OGTT)
gip_ogtt_max <- apply(dt_imputed_ogtt_gip, 2, max, na.rm = TRUE)

# Max GIP for each subject (IIGI)
gip_iigi_max <- apply(dt_imputed_iigi_gip, 2, max, na.rm = TRUE)

# 3. GLP-1 and GIP at 30 min
glp1_ogtt_30_df = ogtt_glp1[which(ogtt_glp1$TimePoint == 30), c("SubjectID", "GLP1")]
colnames(glp1_ogtt_30_df) = c("SubjectID", "glp1_30min")
gip_ogtt_30_df = ogtt_gip[which(ogtt_gip$TimePoint == 30), c("SubjectID", "GIP")]
colnames(gip_ogtt_30_df) = c("SubjectID", "gip_30min")

# 4.  GLP-1 and GIP at 120 min
glp1_ogtt_120_df = ogtt_glp1[which(ogtt_glp1$TimePoint == 120), c("SubjectID", "GLP1")]
colnames(glp1_ogtt_120_df) = c("SubjectID", "glp1_120min")
gip_ogtt_120_df = ogtt_gip[which(ogtt_gip$TimePoint == 120), c("SubjectID", "GIP")]
colnames(gip_ogtt_120_df) = c("SubjectID", "gip_120min")


# Convert named lists to data frames first
glp1_ogtt_auc_df <- data.frame(SubjectID = names(glp1_ogtt_auc), glp1_ogtt_auc = unlist(glp1_ogtt_auc))
glp1_iigi_auc_df <- data.frame(SubjectID = names(glp1_iigi_auc), glp1_iigi_auc = unlist(glp1_iigi_auc))
gip_ogtt_auc_df <- data.frame(SubjectID = names(gip_ogtt_auc), gip_ogtt_auc = unlist(gip_ogtt_auc))
gip_iigi_auc_df <- data.frame(SubjectID = names(gip_iigi_auc), gip_iigi_auc = unlist(gip_iigi_auc))
glp1_ogtt_max_df <- data.frame(SubjectID = names(glp1_ogtt_max), glp1_ogtt_max = unlist(glp1_ogtt_max))
glp1_iigi_max_df <- data.frame(SubjectID = names(glp1_iigi_max), glp1_iigi_max = unlist(glp1_iigi_max))
gip_ogtt_max_df <- data.frame(SubjectID = names(gip_ogtt_max), gip_ogtt_max = unlist(gip_ogtt_max))
gip_iigi_max_df <- data.frame(SubjectID = names(gip_iigi_max), gip_iigi_max = unlist(gip_iigi_max))

# Combine ie with GLP1 and GIP features
ie_df <- ie #data.frame(SubjectID = names(ie), ie = unlist(ie))
ie_df = merge(ie_df, glp1_ogtt_auc_df, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, glp1_iigi_auc_df, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, gip_ogtt_auc_df, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, gip_iigi_auc_df, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, glp1_ogtt_max_df, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, glp1_iigi_max_df, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, gip_ogtt_max_df, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, gip_iigi_max_df, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, glp1_ogtt_30_df, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, gip_ogtt_30_df, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, glp1_ogtt_120_df, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, gip_ogtt_120_df, by = "SubjectID", all.x = TRUE)


# Save the combined results to CSV
write.csv(ie_df, file = file.path(DATA_DIR, "ie_glp1_gip_features_calculations.csv"))
