library(corrr)
library(corrplot)
library(dplyr)
library(ggdendro)
library(dendextend)
library(grid)
library(goeveg)
library(glmnet)
library(ggplot2)
library(ggfortify)
library(ggrepel)
library(ggnetwork)
library(ggparty)
library(Hmisc)
library(hrbrthemes)
library(imputeTS)
library(intergraph)
library(jtools)
library(MESS)
library(Metrics)
library(naniar)
library(pheatmap)
library(plyr)
library(PerformanceAnalytics)
library(RColorBrewer)
library(scales)
library(viridis)
library(zoo)
library(reshape2)
library(data.table)
library(dtw)
setwd("/Users/ahmedm/Dropbox/Diabetes_v4_resubmission/Metabolic_Subphenotype_Predictor_07102022/")


### Prepare dataframes of validation/CGM cohort
glucose = read.csv("/Users/ahmedm/Dropbox/Diabetes_v4_resubmission/Metabolic_Subphenotype_Predictor_07102022/data/ogtt_glucose_all_cohorts_09072023.csv")
glucose_validation_cgm = subset(glucose, exp_type == "venous_with_matching_cgm_and_with_planned_athome_cgm")
glucose_validation_cgm = glucose_validation_cgm[,-which(colnames(glucose_validation_cgm) == "exp_type")]
head(glucose_validation_cgm)

## Remove subjects from study because they are not eligible based on our exculsion criteria
subjects_tobe_removed = c("S061", "S105", "S108", "S109", "S112", "S115")
glucose_validation_cgm_filtered <- glucose_validation_cgm[!(glucose_validation_cgm$subject_id %in% subjects_tobe_removed), ]



#### Visualize multiple tests for each subjects
gp = ggplot(glucose_validation_cgm_filtered, aes(timepoint, glucose)) + 
  geom_line(aes(group=sample_location_extraction_method, colour = sample_location_extraction_method), alpha=0.6, size = 1) + 
  geom_hline(yintercept=50, linetype="dashed", color = "yellow", alpha = 0.6) +
  geom_hline(yintercept=140, linetype="dashed", color = "orange", alpha = 0.6) +  
  geom_hline(yintercept=200, linetype="dashed", color = "red", alpha = 0.6) +
  geom_vline(xintercept=120, linetype="dashed", alpha = 0.6) +  
  facet_wrap(~subject_id, nrow = 6) +
  # scale_fill_manual(values = group.color) + 
  # scale_color_manual(values = group.color) + 
  theme_bw() +
  theme(legend.position="top") + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=45, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  labs(y = "Glucose (mg/dL)", x = "Time (mins)") +
  ggtitle("")
gp
ggsave('/Users/ahmedm/Dropbox/Diabetes_v4_resubmission/Metabolic_Subphenotype_Predictor_07102022/figures/validation_ogtt_11112023.pdf', h=9, w=12);
getwd()



###### Concordance analysis
get_corr = function(df = NULL, test_1 = NULL, test_2 = NULL){
  subjects = as.character(unique(df$subject_id))
  # print(subjects)
  ogtt_cor = vector()
  ogtt_dist = vector()
  ogtt_rmse = vector()
  ogtt_diff = vector()
  
  for(subj in subjects){
    tmp = subset(df, subject_id == subj)
    tmp = na.omit(tmp[,c(test_1, test_2)])
    test_1_list = as.numeric(unlist(tmp[test_1]))
    test_2_list = as.numeric(unlist(tmp[test_2]))
    print(subj)
    ogtt_cor = c(ogtt_cor, cor(test_1_list, test_2_list))
    print(ogtt_cor)
  }
  return(ogtt_cor)
}

#### Raw glucose values
subjects_tobe_removed_2 = c("S101", "S104") ## removed because of non-adherence
glucose_validation_cgm_filtered_2 <- glucose_validation_cgm_filtered[!(glucose_validation_cgm_filtered$subject_id %in% subjects_tobe_removed_2), ]
glucose_validation_cgm_filtered_dcast = reshape2::dcast(glucose_validation_cgm_filtered_2[,c("subject_id", "timepoint", "sample_location_extraction_method","glucose")], subject_id + timepoint ~ sample_location_extraction_method)
glucose_validation_cgm_filtered_dcast$Home_CGM_mean = rowMeans(glucose_validation_cgm_filtered_dcast[c("Home_CGM_1", "Home_CGM_2")], na.rm=TRUE)
glucose_validation_cgm_filtered_dcast$subject_id

ogtt_ctru_venous_cgm_cor = get_corr(df = glucose_validation_cgm_filtered_dcast, test_1 = "CTRU_Venous", test_2 = "CTRU_CGM")
ogtt_ctru_venous_athome_cgm_cor = get_corr(df = glucose_validation_cgm_filtered_dcast, test_1 = "CTRU_Venous", test_2 = "Home_CGM_1")
ogtt_ctru_cgm_athome_cgm_cor = get_corr(df = glucose_validation_cgm_filtered_dcast, test_1 = "CTRU_CGM", test_2 = "Home_CGM_1")
ogtt_at_home_cgm_1_2 = get_corr(df = glucose_validation_cgm_filtered_dcast, test_1 = "Home_CGM_1", test_2 = "Home_CGM_2")
ogtt_ctru_venous_athome_cgm_mean_cor = get_corr(df = glucose_validation_cgm_filtered_dcast, test_1 = "CTRU_Venous", test_2 = "Home_CGM_mean")
ogtt_ctru_cgm_athome_cgm_mean_cor = get_corr(df = glucose_validation_cgm_filtered_dcast, test_1 = "CTRU_CGM", test_2 = "Home_CGM_mean")
subjects = as.character(unique(glucose_validation_cgm_filtered_dcast$subject_id))
ogtt_concordance = data_frame(subjects, 
                              CTRU_Venous_vs_CGM = ogtt_ctru_venous_cgm_cor, 
                              AtHome_CGM1_vs_CGM2 = ogtt_at_home_cgm_1_2, 
                              CTRU_CGM_vs_AtHome_CGM = ogtt_ctru_cgm_athome_cgm_mean_cor)

head(ogtt_concordance)
write.csv(ogtt_concordance, "data/ogtt_concordance_among_tests_subset_12162023_new.csv")
print(ogtt_concordance)
ogtt_concordance_melt = reshape2::melt(ogtt_concordance)
cdata_ogtt_concordance_melt <- ddply(ogtt_concordance_melt, c("variable"), summarise,
                                     N    = length(value),
                                     median = median(value, na.rm=TRUE),
                                     mean = mean(value, na.rm=TRUE),
                                     sd   = sd(value, na.rm=TRUE),
                                     se   = sd / sqrt(N)
)

cdata_ogtt_concordance_melt

### Box plot
cbPalette <- c("#756bb1", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#c51b8a")
gp = ggplot(ogtt_concordance_melt, aes(x=variable, y=value, color=variable)) +
  geom_boxplot() + 
  geom_point() +
  scale_fill_manual(values=cbPalette) + 
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
  # theme(legend.position="top") + 
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Correlation among different OGTT test settings")
gp
ggsave(sprintf('/Users/ahmedm/Dropbox/Diabetes_v4_resubmission/Metabolic_Subphenotype_Predictor_07102022/figures/OGTT_tests_correlations_12162023.pdf'), 
       h=7, w=7)


wilcox.test(ogtt_concordance$CTRU_Venous_vs_CGM, ogtt_concordance$AtHome_CGM1_vs_CGM2, alternative = "two.sided")
wilcox.test(ogtt_concordance$CTRU_Venous_vs_CGM, ogtt_concordance$CTRU_CGM_vs_AtHome_CGM, alternative = "two.sided")
wilcox.test(ogtt_concordance$AtHome_CGM1_vs_CGM2, ogtt_concordance$CTRU_CGM_vs_AtHome_CGM, alternative = "two.sided")


### Correlation across all samples
glucose_lm = lm(home_CGM_mean ~ CTRU_Venous, data = validation_cohort_processed_glucose)
summary(glucose_lm)
glucose_lm = lm(CTRU_Venous ~ home_CGM_mean, data = validation_cohort_processed_glucose)
glucose_lm = lm(home_CGM_mean ~ CTRU_Venous, data = validation_cohort_processed_glucose)
summary(glucose_lm)
summ(glucose_lm, scale = TRUE)
