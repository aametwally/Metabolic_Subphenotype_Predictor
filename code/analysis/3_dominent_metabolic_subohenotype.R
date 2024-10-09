# Load necessary libraries
library(ggplot2) # for plotting
library(dplyr)   # for data manipulation
library(pheatmap) # for heatmap visualization
library(corrplot) # for correlation matrix visualization
library(Hmisc)   # for rcorr function
library(corrr)
library(PerformanceAnalytics)


# Set working directories and paths for data and output
BASE_DIR = "PATH-to-BASE-DIRECTORY"
setwd(BASE_DIR)
DATA_DIR <- file.path(BASE_DIR, "data")
OUTPUT_FIGURE_DIR <- file.path(BASE_DIR, "figures")

# Load metadata
metadata <- read.csv(file.path(DATA_DIR, "aggregated_metabolic_indicators.csv"))

# Filter metadata for the initial cohort
metadata_initial_cohort <- metadata %>%
  filter(exp_type == "venous_without_matching_cgm_and_without_planned_athome_cgm")

# Reorder metadata based on A1C levels in decreasing order
metabolic_indicators_orderd <- metadata_initial_cohort %>%
  arrange(desc(a1c))

# Lock the order of subject IDs based on the reordered data
metabolic_indicators_orderd$subject_id <- factor(metabolic_indicators_orderd$subject_id, 
                                                 levels = metabolic_indicators_orderd$subject_id)


## Determining Dominant Metabolic phenotype ##############
# Calculate deviance from average cohort
calculate_deviance <- function(data, variable) {
  (data[[variable]] - mean(data[[variable]], na.rm = TRUE)) / sd(data[[variable]], na.rm = TRUE)
}

# Apply the deviance calculation for each metabolic indicator
metabolic_indicators_orderd <- metabolic_indicators_orderd %>%
  mutate(
    MuscleIR_Deviance = calculate_deviance(., "sspg"),
    BetaCell_Deviance = -calculate_deviance(., "di"),
    Incretin_Deviance = -calculate_deviance(., "ie"),
    HepaticIR_Deviance = calculate_deviance(., "hepatic_ir")
  )

# Set row names to subject IDs
rownames(metabolic_indicators_orderd) <- metabolic_indicators_orderd$subject_id

# Write results to CSV
write.csv(metabolic_indicators_orderd, file.path(DATA_DIR, "metabolic_indicators_risk_score_10272023.csv"))

# Create a heatmap based on deviance of metabolic phenotypes
pdf(file.path(OUTPUT_FIGURE_DIR, "DominantMetabolicPhenotype_with_Total_MSP_Deviance.pdf"), width = 5, height = 7)
pheatmap(metabolic_indicators_orderd[, c("MuscleIR_Deviance", "BetaCell_Deviance", 
                                         "Incretin_Deviance", "HepaticIR_Deviance")],
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         color = RColorBrewer::brewer.pal(100, "YlOrRd"),
         fontsize_row = 10,
         fontsize_col = 12,
         main = "",
         angle_col = 45)
dev.off()

# Summarize metabolic subphenotype variables
summary_stats <- apply(metabolic_indicators_orderd[, c("sspg", "ie", "di", "hepatic_ir", "a1c", "fpg", "ogtt_2h")], 2, summary)

# Correlation between the deviance of metabolic subphenotypes
metabolic_indicators_tstat <- metabolic_indicators_orderd[, c("MuscleIR_Deviance", "BetaCell_Deviance", 
                                                              "Incretin_Deviance", "HepaticIR_Deviance")]

# Plot correlation network
pdf(file.path(OUTPUT_FIGURE_DIR, "DominantMetabolicPhenotype_correlation.pdf"), width = 6, height = 6)
metabolic_indicators_tstat %>%
  correlate() %>%
  network_plot(min_cor = 0.1, repel = TRUE, legend = "full", curved = TRUE,
               colors = c("#00AFBB", "white", "#FC4E07"))
dev.off()

# Create a correlation plot with p-values
pdf(file.path(OUTPUT_FIGURE_DIR, "DominantMetabolicPhenotype_correlation_pvalue.pdf"), width = 7, height = 7)
chart.Correlation(metabolic_indicators_tstat, histogram = TRUE, pch = 19)
dev.off()

# Calculate correlation coefficients with p-values
correlation_results <- rcorr(as.matrix(metabolic_indicators_tstat))

