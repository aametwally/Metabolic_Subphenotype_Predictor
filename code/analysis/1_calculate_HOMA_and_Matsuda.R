# Libraries needed for this script (unused libraries removed)
library(dplyr)        # Data manipulation
library(ggplot2)      # Plotting
library(MESS)         # AUC calculation
library(data.table)
library(imputeTS)
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

### Load OGTT and IIGI glucose measurements
ogtt <- fread(file.path(DATA_DIR, "ogtt.tsv")) 
colnames(ogtt)[which(colnames(ogtt) == "userID")] = "SubjectID"
ogtt <- ogtt %>% # rename(SubjectID = userID) %>%
  mutate(SubjectID = gsub("43883-0", "S", SubjectID)) %>%
  filter(SubjectID %in% demographics$subject_id) %>%
  mutate(ogtt_value = as.numeric(ogtt_value))


############################
####### HOMA (IR/B/S)
############################

# Load baseline glucose and insulin data from OGTT
ogtt_glucose_baseline <- subset(ogtt, ogtt_min == 0)[, c("SubjectID", "ogtt_min", "ogtt_value")]
colnames(ogtt_glucose_baseline) <- c("SubjectID", "TimePoint", "Glucose")

ogtt_insulin_baseline <- subset(ogtt_analytes, TimePoint == 0 & Assay == "OGTT")[, c("SubjectID", "TimePoint", "Insulin")]
ogtt_insulin_baseline <- na.omit(ogtt_insulin_baseline)

# Merge glucose and insulin data for HOMA calculations
homa_df <- merge(ogtt_glucose_baseline, ogtt_insulin_baseline)

# Convert glucose from mg/dL to mmol/L (divide by 18)
homa_df$Glucose_mmol_L <- round(homa_df$Glucose / 18, 2)

# HOMA-B and HOMA-IR calculations
homa_df$HOMA_B <- round((20 * homa_df$Insulin) / (homa_df$Glucose_mmol_L - 3.5), 2)
homa_df$HOMA_IR <- round((homa_df$Insulin * homa_df$Glucose_mmol_L) / 22.5, 2)

# Calculate HOMA-S (insulin sensitivity)
homa_df$HOMA_S <- round(1 / homa_df$HOMA_IR, 2)


##################################
####### Plot HOMA-IR and HOMA-B
##################################

# Merge with demographic data and order by A1C levels
# homa_df_ordered <- merge(homa_df, demographics, by = "SubjectID")
homa_df_ordered <- merge(demographics, homa_df, by.x = "subject_id", by.y = "SubjectID")
homa_df_ordered <- homa_df_ordered[order(homa_df_ordered$a1c, decreasing = TRUE),]
homa_df_ordered$subject_id <- factor(homa_df_ordered$subject_id, levels = homa_df_ordered$subject_id)

# Plot HOMA-IR
gp <- ggplot(homa_df_ordered, aes(x = subject_id, y = HOMA_IR)) +
  geom_segment(aes(xend = subject_id, y = 0, yend = HOMA_IR), color = "pink") +
  geom_point(size = 2) +
  coord_flip() +
  theme_light() +
  theme(panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("HOMA-IR") +
  xlab("") +
  ylab("HOMA-IR")
gp
ggsave(file.path(OUTPUT_FIGURE_DIR, "HOMA_IR.pdf"), h = 6, w = 3)

# Plot HOMA-B
gp <- ggplot(homa_df_ordered, aes(x = subject_id, y = HOMA_B)) +
  geom_segment(aes(xend = subject_id, y = 0, yend = HOMA_B), color = "green") +
  geom_point(size = 2) +
  coord_flip() +
  theme_light() +
  theme(panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("HOMA-B") +
  xlab("") +
  ylab("HOMA-B")
gp
ggsave(file.path(OUTPUT_FIGURE_DIR, "HOMA_B.pdf"), h = 6, w = 3)

##################################
####### Matsuda Index  
##################################

# Calculate mean glucose from 0 to 120 minutes
ogtt_glucose_mean <- ddply(subset(ogtt, ogtt_min >= 0 & ogtt_min <= 120), ~SubjectID, 
                           summarise, ogtt_glucose_mean = mean(ogtt_value, na.rm = TRUE))
ogtt_glucose_mean$ogtt_glucose_mean_mmol_L <- ogtt_glucose_mean$ogtt_glucose_mean / 18

# Calculate mean insulin from 0 to 120 minutes
mean_insulin_df <- subset(ogtt_analytes, TimePoint >= 0 & TimePoint <= 120)[, c("SubjectID", "TimePoint", "Insulin")]
ogtt_insulin_mean <- ddply(na.omit(mean_insulin_df), ~SubjectID, summarise, ogtt_insulin_mean = mean(Insulin, na.rm = TRUE))

# Matsuda index calculation
matsuda_df <- merge(ogtt_glucose_mean, ogtt_insulin_mean)
matsuda_df <- merge(matsuda_df, homa_df)
matsuda_df$matsuda_index <- 10000 / sqrt(matsuda_df$ogtt_glucose_mean_mmol_L * 
                                           matsuda_df$ogtt_insulin_mean * 
                                           matsuda_df$Glucose_mmol_L * 
                                           matsuda_df$Insulin)


# Plot Matsuda index
matsuda_ordered <- merge(demographics, matsuda_df, by.x = "subject_id", by.y = "SubjectID")
# matsuda_ordered <- merge(matsuda_df, demographics, by = "SubjectID")
matsuda_ordered <- matsuda_ordered[order(matsuda_ordered$a1c, decreasing = TRUE),]
matsuda_ordered$subject_id <- factor(matsuda_ordered$subject_id, levels = matsuda_ordered$subject_id)

gp <- ggplot(matsuda_ordered, aes(x = subject_id, y = matsuda_index)) +
  geom_segment(aes(xend = subject_id, y = 0, yend = matsuda_index), color = "brown") +
  geom_point(size = 2) +
  coord_flip() +
  theme_light() +
  theme(panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Matsuda-Index") +
  xlab("") +
  ylab("Matsuda-Index")
gp
ggsave(file.path(OUTPUT_FIGURE_DIR, "Matsuda_Index.pdf"), h = 6, w = 3)



# Save Matsuda index results
write.csv(matsuda_df[,c("SubjectID",	"HOMA_B",	"HOMA_IR",	"HOMA_S",	"matsuda_index")], file = file.path(DATA_DIR, "homa_and_matsuda_calculations.csv"), row.names = FALSE)
