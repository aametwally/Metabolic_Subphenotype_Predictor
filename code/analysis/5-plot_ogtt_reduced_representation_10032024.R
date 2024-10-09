# Necessary libraries
library(dplyr)        # For data manipulation
library(ggplot2)      # For plotting
library(ggrepel)      # For text labels in ggplot
library(scales)       # For customizing scales in plots

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



# ######################
# ### PCA of OGTT TS
# ######################
ogtt_imputed_normalized = read.csv(file.path(DATA_DIR, "ogtt_glucose_venous_exp_without_matching_imputed_normalized.csv"), row.names = 1) 
df_pca = prcomp(t(ogtt_imputed_normalized))
write.csv(df_pca$x, file = file.path(DATA_DIR, "ogtt_glucose_venous_exp_without_matching_pcs.csv"))
df_pca_annotation = merge(df_pca$x, metabolic_indicators_orderd, by.x = "row.names", by.y = "subject_id", all.x = TRUE)
colnames(df_pca_annotation)[which(colnames(df_pca_annotation) == "Row.names")] = "SubjectID"
write.csv(df_pca_annotation, file.path(DATA_DIR, "ogtt_glucose_venous_exp_without_matching_pca_annotation.csv"), row.names = FALSE)


### Plotting PCA with metabolic phenotype
plot_pca_metabolic_subphenotype = function(df = NULL, subphenotype = NULL, text = NULL, group_color = NULL) {
  gp = ggplot(df , aes_string(x="PC1",y="PC2", color = subphenotype, label = "SubjectID" )) +
    geom_text_repel() +
    geom_point() +
    scale_color_manual(values=group_color) +
    theme_bw() +
    labs(title = "") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(title = text) +
    theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="plain"),
          axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="plain"),
          axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
          axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
          legend.text = element_text(size=11, face="plain"),
          legend.title = element_blank(),
          legend.position="top",
          plot.title = element_text(size=11)) +
    scale_x_continuous(breaks = waiver())
  gp
  ggsave(file.path(OUTPUT_FIGURE_DIR, paste(text, '.pdf', sep = "")) , h=4, w=6)
}


### PCA colored based on Muscle IR
df_pca_annotation$sspg_2_classes <- factor(df_pca_annotation$sspg_2_classes, levels = c("IS", "IR"))
plot_pca_metabolic_subphenotype(df = df_pca_annotation, subphenotype = "sspg_2_classes", 
                                group_color = c(IR = "#333BFF", IS = "#CC6600"), 
                                text = "Reduced representation of glucose-series colored by Muscle IR")

### PCA colored based on Beta Cell
df_pca_annotation$di_2_classes <- factor(df_pca_annotation$di_2_classes, 
                                                levels = c("Normal", "Dysfunction"))
plot_pca_metabolic_subphenotype(df = df_pca_annotation, subphenotype = "di_2_classes",
                                text = "Reduced representation of glucose-series colored by Beta-cell Dysfunction",
                                group_color = c(Normal = "green", Dysfunction = "red"))

### PCA colored based on IE
df_pca_annotation$ie_2_classes <- factor(df_pca_annotation$ie_2_classes, 
                                                levels = c("Normal", "Dysfunction"))
plot_pca_metabolic_subphenotype(df = df_pca_annotation, subphenotype = "ie_2_classes", 
                                group_color = c(Normal = "#66c2a5", Dysfunction = "#fc8d62"),
                                text = "Reduced representation of glucose-series colored by Incretin Effect")

### PCA colored based on Hepatic IR
df_pca_annotation$hepatic_ir_2_classes <- factor(df_pca_annotation$hepatic_ir_2_classes, 
                                                        levels = c("IS", "IR"))
plot_pca_metabolic_subphenotype(df = df_pca_annotation, subphenotype = "hepatic_ir_2_classes",
                                text = "Reduced representation of glucose-series colored by Hepatic IR",
                                group_color = c(IS = "#3288BD", IR = "#D53E4F"))
