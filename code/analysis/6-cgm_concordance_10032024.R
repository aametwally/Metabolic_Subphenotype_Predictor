# Required Libraries
library(ggplot2)          # For creating visualizations
library(reshape2)         # For reshaping data frames
library(dplyr)            # For data manipulation
library(plyr)             # For data aggregation (used with ddply)
library(ggrepel)          # For text annotations in ggplot2
library(Hmisc)            # For statistical tests


# Set working directories and paths for data and output
BASE_DIR = "PATH-to-BASE-DIRECTORY"
setwd(BASE_DIR)
DATA_DIR <- file.path(BASE_DIR, "data") 
OUTPUT_FIGURE_DIR <- file.path(BASE_DIR, "figures") 

# Load glucose data
glucose <- read.csv(file.path(DATA_DIR, "all_cohort_metabolicsubphenotyping_ogtt_glucose_09102023.csv"))


# Filter for validation cohort with matching CGM data
glucose_validation_cgm <- glucose %>%
  filter(exp_type == "venous_with_matching_cgm_and_with_planned_athome_cgm") %>%
  select(-exp_type)

# Remove subjects based on exclusion criteria
subjects_tobe_removed <- c("S061", "S105", "S108", "S109", "S112", "S115")
glucose_validation_cgm_filtered <- glucose_validation_cgm %>%
  filter(!(subject_id %in% subjects_tobe_removed))

# Plot glucose levels across time for each subject
gp <- ggplot(glucose_validation_cgm_filtered, aes(timepoint, glucose)) + 
  geom_line(aes(group = sample_location_extraction_method, colour = sample_location_extraction_method), alpha = 0.6, size = 1) + 
  # geom_hline(yintercept = 50, linetype = "dashed", color = "yellow", alpha = 0.6) +
  geom_hline(yintercept = 140, linetype = "dashed", color = "orange", alpha = 0.6) +  
  geom_hline(yintercept = 200, linetype = "dashed", color = "red", alpha = 0.6) +
  geom_vline(xintercept = 120, linetype = "dashed", alpha = 0.6) +  
  facet_wrap(~subject_id, nrow = 6) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 12, angle = 45, hjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.text=element_text(size=15, face="plain"), legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text.x = element_text(size = 12, face = "bold")
  ) +
  labs(y = "Glucose (mg/dL)", x = "Time (mins)")

# Save the plot
ggsave(file.path(OUTPUT_FIGURE_DIR, "concordance_metrics_ogtt_plasma_and_cgm.pdf"), plot = gp, height = 9, width = 12)

# Function to calculate correlation between two tests
#' @param df Data frame containing the test data
#' @param test_1 First test to compare
#' @param test_2 Second test to compare
#' @return A vector of correlation coefficients
get_corr <- function(df, test_1, test_2) {
  subjects <- unique(df$subject_id)
  correlations <- sapply(subjects, function(subj) {
    tmp <- df %>%
      filter(subject_id == subj) %>%
      select(all_of(c(test_1, test_2))) %>%
      na.omit()
    test_1_list <- as.numeric(tmp[[test_1]])
    test_2_list <- as.numeric(tmp[[test_2]])
    cor(test_1_list, test_2_list)
  })
  return(correlations)
}

# Apply exclusion criteria for additional subjects
subjects_tobe_removed_2 <- c("S101", "S104")
glucose_validation_cgm_filtered_2 <- glucose_validation_cgm_filtered %>%
  filter(!(subject_id %in% subjects_tobe_removed_2))

# Reshape data for further analysis
glucose_validation_cgm_filtered_dcast <- dcast(glucose_validation_cgm_filtered_2, subject_id + timepoint ~ sample_location_extraction_method)
glucose_validation_cgm_filtered_dcast$Home_CGM_mean <- rowMeans(glucose_validation_cgm_filtered_dcast[c("Home_CGM_1", "Home_CGM_2")], na.rm = TRUE)

# Calculate correlation between tests
ogtt_ctru_venous_cgm_cor <- get_corr(glucose_validation_cgm_filtered_dcast, "CTRU_Venous", "CTRU_CGM")
ogtt_ctru_venous_athome_cgm_cor <- get_corr(glucose_validation_cgm_filtered_dcast, "CTRU_Venous", "Home_CGM_1")
ogtt_ctru_cgm_athome_cgm_cor <- get_corr(glucose_validation_cgm_filtered_dcast, "CTRU_CGM", "Home_CGM_1")
ogtt_at_home_cgm_1_2 <- get_corr(glucose_validation_cgm_filtered_dcast, "Home_CGM_1", "Home_CGM_2")
ogtt_ctru_venous_athome_cgm_mean_cor <- get_corr(glucose_validation_cgm_filtered_dcast, "CTRU_Venous", "Home_CGM_mean")
ogtt_ctru_cgm_athome_cgm_mean_cor <- get_corr(glucose_validation_cgm_filtered_dcast, "CTRU_CGM", "Home_CGM_mean")

# Combine results into a data frame
ogtt_concordance <- data_frame(
  subjects = unique(glucose_validation_cgm_filtered_dcast$subject_id),
  CTRU_Venous_vs_CGM = ogtt_ctru_venous_cgm_cor,
  AtHome_CGM1_vs_CGM2 = ogtt_at_home_cgm_1_2,
  CTRU_CGM_vs_AtHome_CGM = ogtt_ctru_cgm_athome_cgm_mean_cor
)

# Export concordance metrics
write.csv(ogtt_concordance, file.path(DATA_DIR, "concordance_metrics_ogtt_plasma_and_cgm.csv"))

# Reshape data for visualization
ogtt_concordance_melt <- melt(ogtt_concordance)

# Summarize the concordance data
cdata_ogtt_concordance_melt <- ddply(ogtt_concordance_melt, "variable", summarise,
                                     N = length(value),
                                     median = median(value, na.rm = TRUE),
                                     mean = mean(value, na.rm = TRUE),
                                     sd = sd(value, na.rm = TRUE),
                                     se = sd / sqrt(N))

# Boxplot of correlation values
gp <- ggplot(ogtt_concordance_melt, aes(x = variable, y = value, color = variable)) +
  geom_boxplot() +
  geom_point() +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1.1)) +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=10, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=10, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 14, color = "black", face = "bold")) +
  labs(y = "Pearson Corr. Coef.", x = "") +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Correlation among different OGTT test settings")

# Save the boxplot
ggsave(file.path(OUTPUT_FIGURE_DIR, "correlations_between_ogtt_plasma_and_cgm.pdf"), plot = gp, height = 7, width = 7)


# Statistical analysis (Wilcoxon test)
wilcox.test(ogtt_concordance$CTRU_Venous_vs_CGM, ogtt_concordance$AtHome_CGM1_vs_CGM2, alternative = "two.sided")
wilcox.test(ogtt_concordance$CTRU_Venous_vs_CGM, ogtt_concordance$CTRU_CGM_vs_AtHome_CGM, alternative = "two.sided")
wilcox.test(ogtt_concordance$AtHome_CGM1_vs_CGM2, ogtt_concordance$CTRU_CGM_vs_AtHome_CGM, alternative = "two.sided")

