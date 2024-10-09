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


################################
###### Disposition Index
################################
#' Process insulin secretion curves from given OGTT/IIGI files
#' @param dir_path Directory path containing secretion result files
#' @return A data frame with insulin secretion rates for each subject
process_insulin_secretion <- function(dir_path) {
  files <- list.files(dir_path, recursive = TRUE, pattern = "RESULTS.TXT$")
  insulin_rate <- data.frame()
  s_name <- vector()
  
  for (f in files) {
    xx <- read.table(file.path(dir_path, f), fill = TRUE)
    s <- strsplit(f, "/")[[1]][1]
    s_name <- c(s_name, s)
    
    # Extract insulin rates from relevant rows and columns
    insulin_rate_tmp <- c(as.numeric(as.character(xx[11, 2])), as.numeric(as.character(xx[12:23, 3])))
    insulin_rate <- rbind(insulin_rate, insulin_rate_tmp)
  }
  
  rownames(insulin_rate) <- paste("43883-", s_name, sep = "")
  colnames(insulin_rate) <- c("before_0", "Interval_0_15", "Interval_15_30", "Interval_30_45", "Interval_45_60",
                              "Interval_60_75", "Interval_75_90", "Interval_90_105", "Interval_105_120",
                              "Interval_120_135", "Interval_135_150", "Interval_150_165", "Interval_165_180")
  return(insulin_rate)
}

# Process insulin secretion from OGTT and IIGI
insulin_rate_ogtt <- process_insulin_secretion(file.path(DATA_DIR, "isec_results/ogtt"))
insulin_rate_iigi <- process_insulin_secretion(file.path(DATA_DIR, "isec_results/iigi"))

# Remove outlier row based on row names
insulin_rate_ogtt <- insulin_rate_ogtt[rownames(insulin_rate_ogtt) != "43883-018", ]
insulin_rate_iigi <- insulin_rate_iigi[rownames(insulin_rate_iigi) != "43883-018", ]

# Clean row names
rownames(insulin_rate_ogtt) <- gsub("43883-0", "S", rownames(insulin_rate_ogtt))
rownames(insulin_rate_iigi) <- gsub("43883-0", "S", rownames(insulin_rate_iigi))

# Merge insulin rates with demographics data
ogtt_insulin_rate_clinical <- merge(insulin_rate_ogtt, demographics, by.x = "row.names", by.y = "subject_id", all.x = TRUE)
rownames(ogtt_insulin_rate_clinical) <- ogtt_insulin_rate_clinical$Row.names
ogtt_insulin_rate_clinical <- ogtt_insulin_rate_clinical[, -1]

# Export insulin secretion rate combined with clinical demographics
write.csv(ogtt_insulin_rate_clinical, file = file.path(DATA_DIR, "insulin_rate_clinical_demographics_annotation_OGTT.csv"))

### Visualize and compare insulin secretion from OGTT and IIGI
# Prepare combined data for plotting
insulin_rate_combined <- bind_rows(
  reshape2::melt(as.matrix(insulin_rate_ogtt)) %>% mutate(test = "OGTT"),
  reshape2::melt(as.matrix(insulin_rate_iigi)) %>% mutate(test = "IIGI")
)
colnames(insulin_rate_combined) <- c("SubjectID", "interval", "insulin_secretion", "test")
insulin_rate_combined$interval <- as.numeric(gsub(".*_(\\d+)$", "\\1", insulin_rate_combined$interval))
insulin_rate_combined$SubjectID <- factor(insulin_rate_combined$SubjectID, levels = a1c_order)

# Plot insulin secretion comparison
gp <- ggplot(insulin_rate_combined, aes(interval, insulin_secretion, col = test)) +
  geom_line(aes(group = test), alpha = 0.8, size = 1.5) +
  facet_wrap(~SubjectID, nrow = 4) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 10, angle = 45, hjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold")) +
  labs(y = "Insulin Secretion (pmol/kg/min)", x = "Time (mins)", title = "Insulin Secretion (OGTT vs IIGI)")
ggsave(file.path(OUTPUT_FIGURE_DIR, "Insulin_Secretion_OGTT_IIGI.pdf"), plot = gp, height = 6, width = 10)

### Calculate Insulin Secretion / SSPG
ogtt_insulinsecretion_over_sspg <- sweep(insulin_rate_ogtt, 1, demographics$sspg[match(rownames(insulin_rate_ogtt), demographics$subject_id)], `/`)

### Combine ogtt_insulinsecretion_over_sspg  Index with demographics
ogtt_insulinsecretion_over_sspg_demographics_annotation <- merge(ogtt_insulinsecretion_over_sspg, demographics, by.x = "row.names", by.y = "subject_id", all.x = TRUE)
ogtt_insulinsecretion_over_sspg_demographics_annotation$SubjectID = ogtt_insulinsecretion_over_sspg_demographics_annotation$Row.names

# Save insulin secretion rate divided by SSPG
write.csv(ogtt_insulinsecretion_over_sspg_demographics_annotation, file = file.path(DATA_DIR, "ogtt_insulinsecretion_over_sspg_demographics_annotation.csv"))



### Calculate Disposition Index (DI) using AUC over specified timepoints
timepoints <- seq(0, 30, by = 15)
disposition_index <- data.frame(
  di = apply(ogtt_insulinsecretion_over_sspg[, c("before_0", "Interval_0_15", "Interval_15_30")], 1, function(x) {
    MESS::auc(timepoints, x)
  }),
  SubjectID = rownames(ogtt_insulinsecretion_over_sspg)
)
disposition_index = disposition_index[, c("SubjectID", "di")]

# Save DI to file
write.csv(disposition_index, file = file.path(DATA_DIR, "disposition_index_calculations.csv"), row.names = FALSE)
