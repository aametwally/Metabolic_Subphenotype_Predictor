# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

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

# Create a custom theme for consistency across plots
custom_theme <- function() {
  theme_light() +
    theme(panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text.x = element_text(colour = "black", size = 14),
          axis.text.y = element_text(colour = "black", size = 14),
          axis.title.x = element_text(colour = "black", size = 14, face = "bold"),
          axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
          legend.text = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

# Plot HbA1c levels
group.color <- c(Normoglycemic = "#00AFBB", PreDM = "#E7B800", T2DM = "#FC4E07")
gp_a1c <- ggplot(metabolic_indicators_orderd, aes(x = subject_id, y = a1c, color = glycemic_status_a1c)) +
  geom_segment(aes(xend = subject_id, y = 0, yend = a1c, color = glycemic_status_a1c)) +
  geom_point(size = 3) +
  coord_flip() +
  scale_color_manual(values = group.color) +
  custom_theme() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  ggtitle("HbA1c") +
  ylab("HbA1c %") +
  xlab("")
ggsave(file.path(OUTPUT_FIGURE_DIR, "metabolicindicator_a1c.pdf"), plot = gp_a1c, height = 7, width = 3.5)

# Plot SSPG levels
group.color_sspg <- c(IS = "#CC6600", IR = "#333BFF")
metabolic_indicators_orderd$sspg_2_classes <- factor(metabolic_indicators_orderd$sspg_2_classes, 
                                                   levels = c("IS", "IR"))
gp_sspg <- ggplot(metabolic_indicators_orderd, aes(x = subject_id, y = sspg, color = sspg_2_classes)) +
  geom_segment(aes(xend = subject_id, y = 0, yend = sspg, color = sspg_2_classes)) +
  geom_point(size = 3) +
  coord_flip() +
  scale_color_manual(values = group.color_sspg) +
  custom_theme() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  ggtitle("Muscle IR") +
  ylab("SSPG mg/dL") +
  xlab("")
ggsave(file.path(OUTPUT_FIGURE_DIR, "metabolicindicator_sspg.pdf"), plot = gp_sspg, height = 7, width = 3.5)

# Plot DI levels
group.color_di <- c(Normal = "green", Dysfunction = "red")
metabolic_indicators_orderd$di_2_classes <- factor(metabolic_indicators_orderd$di_2_classes, 
                                                   levels = c("Normal", "Dysfunction"))
gp_di <- ggplot(metabolic_indicators_orderd, aes(x = subject_id, y = di, color = di_2_classes)) +
  geom_segment(aes(xend = subject_id, y = 0, yend = di, color = di_2_classes)) +
  geom_point(size = 3) +
  coord_flip() +
  scale_color_manual(values = group.color_di) +
  custom_theme() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  ggtitle("Disposition Index") +
  ylab("Disposition Index") +
  xlab("")
ggsave(file.path(OUTPUT_FIGURE_DIR, "metabolicindicator_di.pdf"), plot = gp_di, height = 7, width = 3.5)


# Plot IE levels
group.color_ie <- c(Normal = "#66c2a5", Dysfunction = "#fc8d62")
metabolic_indicators_orderd$ie_2_classes <- factor(metabolic_indicators_orderd$ie_2_classes, 
                                                   levels = c("Normal", "Dysfunction"))
gp_ie <- ggplot(metabolic_indicators_orderd, aes(x = subject_id, y = ie, color = ie_2_classes)) +
  geom_segment(aes(xend = subject_id, y = 0, yend = ie, color = ie_2_classes)) +
  geom_point(size = 3) +
  coord_flip() +
  scale_color_manual(values = group.color_ie) +
  custom_theme() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  ggtitle("Incretin Effect") +
  ylab("Incretin Effect %") +
  xlab("")
ggsave(file.path(OUTPUT_FIGURE_DIR, "metabolicindicator_ie.pdf"), plot = gp_ie, height = 7, width = 3.5)


# Plot HepaticIR levels
group.color_hepatic_ir <- c(IS = "#3288BD", IR = "#D53E4F")
metabolic_indicators_orderd$hepatic_ir_2_classes <- factor(metabolic_indicators_orderd$hepatic_ir_2_classes, 
                                                     levels = c("IS", "IR"))
gp_hepatic_ir <- ggplot(metabolic_indicators_orderd, aes(x = subject_id, y = hepatic_ir, color = hepatic_ir_2_classes)) +
  geom_segment(aes(xend = subject_id, y = 3, yend = hepatic_ir, color = hepatic_ir_2_classes)) +
  geom_point(size = 3) +
  coord_flip() +
  scale_color_manual(values = group.color_hepatic_ir) +
  custom_theme() +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  ggtitle("Hepatic IR") +
  ylim(3, 5.25) +
  ylab("Hepatic IR Index") +
  xlab("")
ggsave(file.path(OUTPUT_FIGURE_DIR, "metabolicindicator_hepatic_ir.pdf"), plot = gp_hepatic_ir, height = 7, width = 3.5)

