# Load necessary libraries
library(ggplot2) # for plotting
library(dplyr)   # for data manipulation

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

# Function to create histograms for various metabolic indicators
create_histogram <- function(data, x_var, fill_var, title, x_label, output_filename, fill_colors) {
  # Convert the fill variable to a factor with specified levels
  data[[fill_var]] <- factor(data[[fill_var]], levels = unique(data[[fill_var]]))
  
  # Generate histogram
  plot <- ggplot(data = data) +
    geom_histogram(aes_string(x = x_var, fill = fill_var), color = "#e9ecef", alpha = 0.6, position = 'identity') +
    scale_fill_manual(values = fill_colors) + 
    theme_bw() +
    theme(legend.position = "top") +
    ggtitle(title) +
    theme(axis.text.x = element_text(colour = "black", size = 10, face = "bold"),
          axis.text.y = element_text(colour = "black", size = 10),
          axis.title.x = element_text(colour = "black", size = 12, face = "bold"),
          axis.title.y = element_text(colour = "black", size = 12, face = "bold"),
          legend.text = element_text(size = 15),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    xlab(x_label) +
    ylab("Count") +
    labs(fill = "")
  
  # Save the plot
  ggsave(file.path(OUTPUT_FIGURE_DIR, output_filename), height = 3, width = 5)
}

# Histogram for Muscle IR
create_histogram(metabolic_indicators_orderd, "sspg", "sspg_2_classes", "Muscle IR",
                 "SSPG (mg/dL)", "metabolicindicator_histogram_sspg.pdf", 
                 c(IS = "#CC6600", IR = "#333BFF"))

# Histogram for Beta-cell Dysfunction
create_histogram(metabolic_indicators_orderd, "di", "di_2_classes", "Beta-cell Dysfunction",
                 "Disposition Index (pmol*dL)/(kg*ml)", "metabolicindicator_histogram_di.pdf", 
                 c(Normal = "green", Dysfunction = "red"))

# Histogram for Impaired Incretin Effect
create_histogram(metabolic_indicators_orderd, "ie", "ie_2_classes", "Impaired Incretin Effect",
                 "Incretin Effect %", "metabolicindicator_histogram_ie.pdf", 
                 c(Normal = "#66c2a5", Dysfunction = "#fc8d62"))

# Histogram for Hepatic IR
create_histogram(metabolic_indicators_orderd, "hepatic_ir", "hepatic_ir_2_classes", "Hepatic IR",
                 "Hepatic IR Index", "metabolicindicator_histogram_hepatic_ir_index.pdf", 
                 c(IS = "#3288BD", IR = "#D53E4F"))
