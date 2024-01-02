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

setwd("/Users/ahmedm/Dropbox/Diabetes_v4_resubmission/Metabolic_Subphenotype_Predictor_07102022/")

#########################
###### Loading all data
#########################
### Demographics
demographics = read.csv("data/demographics_09302021_AM.csv")

# Rename Normal to Normoglycmic:
demographics[which(demographics["a1c_t2d_status"] == "Normal"), "a1c_t2d_status"] = "Normoglycemic"
demographics[which(demographics["ogtt_t2d_status"] == "Normal"), "ogtt_t2d_status"] = "Normoglycemic"
demographics[which(demographics["fpg_t2d_status"] == "Normal"), "fpg_t2d_status"] = "Normoglycemic"


#### Load PRS
prs = read.csv("data/t2d_prs.csv")
prs = subset(prs, SubjectID %in% demographics$SubjectID)

### Load lipid panel
lipid_panel = read.csv("data/cgm_metadata_dnasamples_all_subj.csv")
colnames(lipid_panel)[which(colnames(lipid_panel) == "subject_id")] = "SubjectID"
lipid_panel$SubjectID = gsub("43883-0","S", lipid_panel$SubjectID)
lipid_panel = subset(lipid_panel, SubjectID %in% demographics$SubjectID)

### Load metabolic panel 
metabolic_df = read.table("data/blood.tsv", sep = "\t", header = TRUE)
metabolic_df = na.omit(metabolic_df)
metabolic_df = metabolic_df[-which(metabolic_df$value == "Exact value cannot be determined"),]
colnames(metabolic_df)[which(colnames(metabolic_df) == "userID")] = "SubjectID"
metabolic_df$SubjectID = gsub("43883-0","S", metabolic_df$SubjectID)
metabolic_df = subset(metabolic_df, SubjectID %in% demographics$SubjectID)

## Load analytes 
washu_combined_batches = read.csv("data/washu_combined_batches_10272020.csv")
washu_combined_batches = washu_combined_batches[with(washu_combined_batches, order(SubjectID, Assay, TimePoint)),]
washu_combined_batches$SubjectID = gsub("43883-0","S", washu_combined_batches$SubjectID)
washu_combined_batches = subset(washu_combined_batches, SubjectID %in% demographics$SubjectID)
ogtt_analytes = subset(washu_combined_batches, Assay == "OGTT")
iigi_analytes = subset(washu_combined_batches, Assay == "IIGI")
sspg_analytes = subset(washu_combined_batches, Assay == "SSPG")

## Load glucose measurements from OGTT/IIGI tests
ogtt = fread('data/ogtt.tsv')
colnames(ogtt)[which(colnames(ogtt) == "userID")] = "SubjectID"
ogtt$SubjectID = gsub("43883-0","S", ogtt$SubjectID)
ogtt = subset(ogtt, SubjectID %in% demographics$SubjectID)
ogtt$ogtt_value = as.numeric(ogtt$ogtt_value)

iigi = fread('data/iigi.tsv')
colnames(iigi)[which(colnames(iigi) == "userID")] = "SubjectID"
iigi$SubjectID = gsub("43883-0","S", iigi$SubjectID)
iigi = subset(iigi, SubjectID %in% demographics$SubjectID)
dim(iigi)
iigi$iigi_value = as.numeric(iigi$iigi_value)


##############################################################
######## Calculate/Infer Metabolic subphenotypes  ############
##############################################################

############################
####### Incretin Effect
############################
## OGTT AUC
ogtt_cpeptide = ogtt_analytes[,c("SubjectID", "TimePoint","C_peptide")]
ogtt_cpeptide_propcessed =  na.omit(ogtt_cpeptide)
dt = t(reshape2::dcast(ogtt_cpeptide_propcessed, TimePoint ~ SubjectID))
colnames(dt) = dt[1,]
dt = dt[-1,]
dt = as.data.frame(dt)
dt_imputed_ogtt = apply(dt, 1, na_interpolation)
timepoints = as.numeric(rownames(dt_imputed_ogtt))
cpeptide_ogtt_auc = apply(dt_imputed_ogtt, 2, function(x){
  MESS::auc(timepoints, x)
})

## IIGI AUC
iigi_cpeptide = iigi_analytes[,c("SubjectID", "TimePoint","C_peptide")]
iigi_cpeptide_propcessed =  na.omit(iigi_cpeptide)
dt = t(reshape2::dcast(iigi_cpeptide_propcessed, TimePoint ~ SubjectID))
colnames(dt) = dt[1,]
dt = dt[-1,]
dt = as.data.frame(dt)
dt_imputed_iigi = apply(dt, 1, na_interpolation)
timepoints = as.numeric(rownames(dt_imputed_iigi)) #seq(from = 0, to = 180, by = 15)
cpeptide_iigi_auc = apply(dt_imputed_iigi, 2, function(x){
  MESS::auc(timepoints, x)
})

## Incretin Effect Calculations
ie = 100*(cpeptide_ogtt_auc - cpeptide_iigi_auc)/cpeptide_ogtt_auc
write.csv(ie, file = "ie_cpeptide-based_calculations_09102021.csv")
ie_div_iigi = 100*(cpeptide_ogtt_auc - cpeptide_iigi_auc)/cpeptide_iigi_auc
ie_no_div = (cpeptide_ogtt_auc - cpeptide_iigi_auc)
inc = data.frame(ie_div_ogtt = ie, ie_div_iigi = ie_div_iigi, ie_no_div = ie_no_div)
write.csv(inc, file = "ie_variants_calculations_09102021.csv")

### Visualize c-peptide OGTT/IIG
cpeptide_ogtt_melted = reshape2::melt(as.matrix(dt_imputed_ogtt))
cpeptide_ogtt_melted$test = "OGTT"
cpeptide_iigi_melted = reshape2::melt(as.matrix(dt_imputed_iigi))
cpeptide_iigi_melted$test = "IIGI" 
cpeptide_combined = rbind(cpeptide_ogtt_melted, cpeptide_iigi_melted)
colnames(cpeptide_combined) = c("time", "SubjectID", "c_pepdite", "test")

#### Rank Subject to visualize their IE based on A1C order
dem_a1c = demographics[, c("SubjectID", "a1c")]
a1c_myList <- as.list(dem_a1c$a1c)
names(a1c_myList) <- dem_a1c$SubjectID
a1c_myList = unlist(a1c_myList)
a1c_order = names(sort(a1c_myList, decreasing = TRUE))
cpeptide_combined_33 = cpeptide_combined[which(cpeptide_combined$SubjectID %in% dem_a1c$SubjectID),]
cpeptide_combined_33$SubjectID = factor(cpeptide_combined_33$SubjectID, levels=a1c_order)

gp = ggplot(cpeptide_combined_33, aes(time, c_pepdite, col = test)) + 
  geom_line(aes(group=test), alpha=0.8, size = 1.5) + 
  scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999")) + 
  facet_wrap(~SubjectID, nrow = 4) +
  theme_bw() +
  theme(legend.position="top") + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=45, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  labs(y = "C-Peptide concentration (ng/mL)", x = "Time (mins)") +
  ggtitle("C-Peptide concentration (OGTT vs IIGI)")
gp
ggsave('cpeptide_OGTT_IIGI_ordered_09102021.pdf', h=6, w=10)


#### Visualize GLP-1 OGTT/IIGI
ogtt_glp1 = ogtt_analytes[,c("SubjectID", "Assay", "TimePoint","GLP1")]
ogtt_glp1 = na.omit(ogtt_glp1)
iigi_glp1 = iigi_analytes[,c("SubjectID", "Assay", "TimePoint","GLP1")]
iigi_glp1 = na.omit(iigi_glp1)
ogtt_iigi_glp1_combined = rbind(ogtt_glp1, iigi_glp1)
ogtt_iigi_glp1_combined$SubjectID = factor(ogtt_iigi_glp1_combined$SubjectID, levels=a1c_order)
gp = ggplot(ogtt_iigi_glp1_combined, aes(TimePoint, GLP1, col = Assay)) + 
  geom_line(aes(group = Assay), alpha=0.8, size = 1.5) + 
  geom_point(aes(group = Assay), alpha=0.8, size = 1.5)  + 
  facet_wrap(~SubjectID, nrow=4) +
  theme_bw() +
  theme(legend.position="top") + 
  scale_color_brewer(palette="Pastel1") + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=45, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  labs(y = "GLP1 (pmol/L)" , x = "Time (mins)") +
  ggtitle("GLP-1 (OGTT vs IIGI)")
gp
ggsave('GLP1_OGTT_IIGI_09102021.pdf', h=6, w=10) 


#### Visualize GIP OGTT/IIGI
ogtt_gip = ogtt_analytes[,c("SubjectID", "Assay", "TimePoint","GIP")]
ogtt_gip = na.omit(ogtt_gip)
iigi_gip = iigi_analytes[,c("SubjectID", "Assay", "TimePoint","GIP")]
iigi_gip = na.omit(iigi_gip)
ogtt_iigi_gip_combined = rbind(ogtt_gip, iigi_gip)
ogtt_iigi_gip_combined$SubjectID = factor(ogtt_iigi_gip_combined$SubjectID, levels=a1c_order)
gp = ggplot(ogtt_iigi_gip_combined, aes(TimePoint, GIP, col = Assay)) + 
  geom_line(aes(group = Assay), alpha=0.8, size = 1.5) + 
  geom_point(aes(group = Assay), alpha=0.8, size = 1.5)  + 
  facet_wrap(~SubjectID, nrow = 4) +
  theme_bw() +
  theme(legend.position="top") + 
  scale_color_brewer(palette="Set2") + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=45, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  labs(y = "GIP (pg/mL)" , x = "Time (mins)") +
  ggtitle("GIP (OGTT vs IIGI)")
gp
ggsave('GIP_OGTT_IIGI_09102021.pdf', h=6, w=10) 


############### Incretins features
## AUC GLP1
glp1_dcast =  reshape2::dcast(ogtt_glp1, SubjectID  ~ TimePoint, value.var = "GLP1")
rownames(glp1_dcast) = glp1_dcast[,1]
glp1_dcast = glp1_dcast[,-1]
glp1_dcast_imputed = apply(glp1_dcast, 1, na_interpolation)
glp1_auc = apply(glp1_dcast_imputed, 2, function(x){
  MESS::auc(as.numeric(rownames(glp1_dcast_imputed)), x,  type = "linear")
})
glp1_auc_df = as.data.frame(glp1_auc)
glp1_auc_df$log_glp1_auc = log(glp1_auc_df$glp1_auc)

## AUC GIP
gip_dcast =  reshape2::dcast(ogtt_gip, SubjectID  ~ TimePoint, value.var = "GIP")
rownames(gip_dcast) = gip_dcast[,1]
gip_dcast = gip_dcast[,-1]
gip_dcast_imputed = apply(gip_dcast, 1, na_interpolation)
gip_auc = apply(gip_dcast_imputed, 2, function(x){
  MESS::auc(as.numeric(rownames(gip_dcast_imputed)), x,  type = "linear")
})
gip_auc_df = as.data.frame(gip_auc)
gip_auc_df$log_gip_auc = log(gip_auc_df$gip_auc)

## retrieve max GLP1 and GIP
ogtt_glp1_max = ddply(ogtt_glp1, ~SubjectID, summarise, GLP1_max = max(GLP1))  
ogtt_gip_max = ddply(ogtt_gip, ~SubjectID, summarise, GIP_max = max(GIP)) 

## retrieve GLP-1 and GIP at 30 min
ogtt_glp1_30 = ogtt_glp1[which(ogtt_glp1$TimePoint == 30), c("SubjectID", "GLP1")]
colnames(ogtt_glp1_30) = c("SubjectID", "GLP1_30min")
ogtt_gip_30 = ogtt_gip[which(ogtt_gip$TimePoint == 30), c("SubjectID", "GIP")]
colnames(ogtt_gip_30) = c("SubjectID", "GIP_30min")

## retrieve GLP-1 and GIP at 120 min
ogtt_glp1_120 = ogtt_glp1[which(ogtt_glp1$TimePoint == 120), c("SubjectID", "GLP1")]
colnames(ogtt_glp1_120) = c("SubjectID", "GLP1_120min")
ogtt_gip_120 = ogtt_gip[which(ogtt_gip$TimePoint == 120), c("SubjectID", "GIP")]
colnames(ogtt_gip_120) = c("SubjectID", "GIP_120min")

## Combine ie, glp1_auc, gip_auc
ie_df = as.data.frame(ie)
ie_df = merge(ie_df, glp1_auc_df, by = "row.names", all.x = TRUE)
ie_df = merge(ie_df, gip_auc_df, by.x = "Row.names", by.y = "row.names", all.x = TRUE)
colnames(ie_df)[which(colnames(ie_df) == "Row.names")] = "SubjectID"
ie_df = merge(ie_df, ogtt_glp1_max, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, ogtt_gip_max, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, ogtt_glp1_30, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, ogtt_gip_30, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, ogtt_glp1_120, by = "SubjectID", all.x = TRUE)
ie_df = merge(ie_df, ogtt_gip_120, by = "SubjectID", all.x = TRUE)
ie_df$log_GLP1_GIP_max_comb = log(ie_df$GLP1_max) + log(ie_df$GIP_max)
dim(ie_df)

##### Thresholding IE
summary(ie_df$ie)
ie_q1 = unname(floor(summary(ie_df$ie)[2]))
ie_q2 = unname(floor(summary(ie_df$ie)[5]))
ie_df$ie_3_classes = "Unknown"
ie_df$ie_3_classes[which(ie_df$ie>ie_q2)] = "Normal"
ie_df$ie_3_classes[which(ie_df$ie<ie_q1)] = "Dysfunction"
ie_df$ie_3_classes[which(ie_df$ie>=ie_q1 & ie_df$ie<=ie_q2)] = "Intermediate"
write.csv(ie_df, "ie_df_glp1_gip_features_09102021.csv")

## Merge IE with demographics. Onluy include samples that have GLP1 and GIP
ie_df_demographics =  merge(ie_df, demographics, by = "SubjectID")
dim(ie_df_demographics)
ie_df_demographics = na.omit(ie_df_demographics)
write.csv(ie_df_demographics, "ie_df_glp1_gip_features_demographics_09102021.csv")

### LM: GIP and GLP1 vs IE
ie_lm = lm(ie ~ GLP1_120min +  bmi + age, data = ie_df_demographics)
summary(ie_lm)
summ(ie_lm, scale = TRUE)


group.color = c(Normal = "#66c2a5", Dysfunction = "#fc8d62", Intermediate = "#984ea3") 
ie_df_demographics$ie_3_classes <- factor(ie_df_demographics$ie_3_classes, 
                                            levels = c("Normal", "Intermediate", "Dysfunction"))

### Visualize IE vs GLP1_120
gp = ggplot(ie_df_demographics, aes(x = GLP1_120min, y = ie, color = ie_3_classes, label = SubjectID)) + 
  geom_point() +
  geom_text_repel() +
  stat_smooth(method = "lm", color= "red") +
  scale_fill_manual(values = group.color) + 
  scale_color_manual(values = group.color) + 
  annotate("text", x = 50, y = 22, label = "beta == 0.54",
           parse = TRUE) +
  annotate("text", x = 50, y = 17, label = "p == 0.08",
           parse = TRUE) +
  annotate("text", x = 50, y = 12, label = "R ^ 2 == 0.15",
           parse = TRUE) +
  ggtitle("Incretin Effect (%) vs GLP-1 (pmol/L)") + 
  labs(x = bquote(''*~GLP1[120] (pmol/L)*''),
       y = "IE (%)") + 
  theme_bw() +
  ylim(c(0,100)) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5,
                                   vjust=0, face="plain"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5,
                                   vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5,
                                    vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5,
                                    vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="top",
        strip.text = element_text(size=12, face = "bold"))
gp
ggsave('ie_glp1_120_12202021.pdf', h=5, w=5);


ie_lm = lm(ie ~ GIP_120min +  bmi + age, data = ie_df_demographics)
summary(ie_lm)
summ(ie_lm, scale = TRUE)


### Visualize IE vs GIP_120
gp = ggplot(ie_df_demographics, aes(x = GIP_120min, y = ie, color = ie_3_classes, label = SubjectID)) + 
  geom_text_repel() +
  geom_point() +
  stat_smooth(method = "lm", color= "red") +
  scale_fill_manual(values = group.color) + 
  scale_color_manual(values = group.color) + 
  annotate("text", x = 400, y = 22, label = "beta == 0.09",
           parse = TRUE) +
  annotate("text", x = 400, y = 17, label = "p == 0.001",
           parse = TRUE) +
  annotate("text", x = 400, y = 12, label = "R ^ 2 == 0.37",
           parse = TRUE) +
  
  ggtitle("Incretin Effect (%) vs GIP (pg/mL)") +  
  labs(x = bquote(''*~GIP[120]  (pg/mL)*''),
       y = "IE (%)") + 
  theme_bw() +
  ylim(c(0,100)) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5,
                                   vjust=0, face="plain"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5,
                                   vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5,
                                    vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5,
                                    vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="top",
        strip.text = element_text(size=12, face = "bold"))
gp
ggsave('ie_gip_120_12202021.pdf', h=5, w=5);


### Visualize IE vs GLP1_120
gp = ggplot(ie_df_demographics, aes(x = GLP1_120min, y = ie, color = ie_3_classes, label = SubjectID)) + 
  geom_point() +
  geom_text_repel() +
  facet_wrap(~SSPG_class_2, scale = "free")+
  stat_smooth(method = "lm", color= "red") +
  scale_fill_manual(values = group.color) + 
  scale_color_manual(values = group.color) + 
  ggtitle("Incretin Effect (%) vs GLP-1 (pmol/L)") + 
  labs(x = bquote(''*~GLP1[120] (pmol/L)*''),
       y = "IE (%)") + 
  theme_bw() +
  ylim(c(0,100)) + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5,
                                   vjust=0, face="plain"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5,
                                   vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5,
                                    vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5,
                                    vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position="top",
        strip.text = element_text(size=12, face = "bold"))
gp
ggsave('ie_glp1_120__perIR09102021.pdf', h=5, w=8);


############################
####### HOMA (IR/B)
############################
### Load baseline glucose/insulin from OGTT
ogtt_glucose_baseline = subset(ogtt, ogtt_min == 0)[,c("SubjectID", "ogtt_min", "ogtt_value")]
colnames(ogtt_glucose_baseline) = c("SubjectID", "TimePoint", "Glucose")
ogtt_insulin_baseline = subset(ogtt_analytes, TimePoint == 0 & Assay == "OGTT")[,c("SubjectID", "TimePoint", "Insulin")]
ogtt_insulin_baseline = na.omit(ogtt_insulin_baseline)

### HOMA calculations
homa_df = merge(ogtt_glucose_baseline, ogtt_insulin_baseline)
homa_df$Glucose_mmol_L = round(homa_df$Glucose/18, 2) ## Divide by 18 to convert glucose units from mg/dL to mmol/L
homa_df$HOMA_B = (20*homa_df$Insulin)/(homa_df$Glucose_mmol_L - 3.5)
homa_df$HOMA_B = round(homa_df$HOMA_B, digits = 2)
homa_df$HOMA_IR = (homa_df$Insulin * homa_df$Glucose_mmol_L)/22.5
homa_df$HOMA_IR = round(homa_df$HOMA_IR, digits = 2)
homa_df$HOMA_S = 1/homa_df$HOMA_IR
homa_df$HOMA_S  = round(homa_df$HOMA_S, digits = 2)
homa_df$HOMA_IR_group = NA
homa_df$HOMA_IR_group[homa_df$HOMA_IR<=2.5] = "IS"
homa_df$HOMA_IR_group[homa_df$HOMA_IR>2.5] = "IR"
write.csv(homa_df, "homa_df_09102021.csv")



#### Plot HOMA-IR and HOMA-B
homa_df_orderd = merge(homa_df, demographics, by = "SubjectID")
homa_df_orderd = homa_df_orderd[order(homa_df_orderd$a1c, decreasing = TRUE),]
homa_df_orderd$SubjectID <- factor(homa_df_orderd$SubjectID, 
                                                levels = homa_df_orderd$SubjectID)

## Plot HOMA-IR
gp = ggplot(homa_df_orderd, aes(x=SubjectID, y=HOMA_IR)) +
  geom_segment( aes(x=SubjectID, xend=SubjectID, y=0, yend=HOMA_IR), color = "pink")+
  geom_point(size=2) +
  coord_flip() +
  scale_fill_manual(values = group.color) + 
  scale_color_manual(values = group.color) + 
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="top",
    legend.title = element_blank()
  ) +
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=12, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("HOMA-IR") +
  xlab("") +
  ylab("HOMA-IR")
gp
ggsave('HOMA_IR_09102021.pdf', h=6, w=3); 

## Plot HOMA-B
gp = ggplot(homa_df_orderd, aes(x=SubjectID, y=HOMA_B)) +
  geom_segment( aes(x=SubjectID, xend=SubjectID, y=0, yend=HOMA_B), color = "green")+
  geom_point(size=2) +
  coord_flip() +
  scale_fill_manual(values = group.color) + 
  scale_color_manual(values = group.color) + 
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="top",
    legend.title = element_blank()
  ) +
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=12, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("HOMA-B") +
  xlab("") +
  ylab("HOMA-B")
gp
ggsave('HOMA_B_09102021.pdf', h=6, w=3); 


##################################
####### Matsuda Index  
##################################
### Mean Glucose
mean_glucose_df = subset(ogtt, ogtt_min >=0 & ogtt_min<=120) 
ogtt_glucose_mean = ddply(mean_glucose_df, ~SubjectID, summarise, ogtt_glucose_mean = mean(ogtt_value, na.rm = TRUE))
ogtt_glucose_mean$ogtt_glucose_mean_mmol_L = ogtt_glucose_mean$ogtt_glucose_mean/18

### Mean Insulin
mean_insulin_df = subset(ogtt_analytes,  TimePoint >=0 & TimePoint<=120)[, c("SubjectID", "TimePoint", "Insulin")]
mean_insulin_df = na.omit(mean_insulin_df)
ogtt_insulin_mean = ddply(mean_insulin_df, ~SubjectID, summarise, ogtt_insulin_mean = mean(Insulin, na.rm = TRUE))

### Matsuda calculations
metsuda_df = merge(ogtt_glucose_mean, ogtt_insulin_mean)
metsuda_df = merge(metsuda_df, homa_df)
metsuda_df$matsuda_index = 10000/sqrt(metsuda_df$ogtt_glucose_mean_mmol_L*metsuda_df$ogtt_insulin_mean*metsuda_df$Glucose_mmol_L*metsuda_df$Insulin)
write.csv(metsuda_df, file = "matsuda_index_09102021.csv", row.names = FALSE)


#### Plot Matusda index
matsuda_orderd = merge(metsuda_df, demographics, by = "SubjectID")
matsuda_orderd = matsuda_orderd[order(matsuda_orderd$a1c, 
                                      decreasing = TRUE),]
matsuda_orderd$SubjectID <- factor(matsuda_orderd$SubjectID, 
                                   levels = matsuda_orderd$SubjectID)

gp = ggplot(matsuda_orderd, aes(x=SubjectID, y=matsuda_index)) +
  geom_segment( aes(x=SubjectID, xend=SubjectID, y=0, yend=matsuda_index), color = "brown")+
  geom_point(size=2) +
  coord_flip() +
  scale_fill_manual(values = group.color) + 
  scale_color_manual(values = group.color) + 
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="top",
    legend.title = element_blank()
  ) +
  theme(axis.text.x = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=12, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Matsuda-Index") +
  xlab("") +
  ylab("Matsuda-Index")
gp
ggsave('matsuda_index_09102021.pdf', h=6, w=3); 



################################
###### Disposition Index
################################
###### Analysis of Insulin Secretion curves from OGTT ########
ogtt_insulin_secretion_files = list.files("data/isec_results/ogtt/", recursive=TRUE, pattern = "RESULTS.TXT$")
insulin_rate = data.frame()
s_name = vector()
for(f in ogtt_insulin_secretion_files){
  xx = read.table(paste("data/isec_results/ogtt/", f, sep = ""), fill = TRUE)
  s =  strsplit(f, "/")[[1]][1]
  s_name = c(s_name, s)
  insulin_rate_tmp = c(as.numeric(as.character(xx[11,2])), as.numeric(as.character(xx[12:23,3])))
  insulin_rate = rbind(insulin_rate, insulin_rate_tmp)
}
rownames(insulin_rate) = paste("43883-", s_name, sep = "")
colnames(insulin_rate) = c("before_0", "Interval_0_15","Interval_15_30", "Interval_30_45", "Interval_45_60", "Interval_60_75", "Interval_75_90", 
                           "Interval_90_105","Interval_105_120", "Interval_120_135", "Interval_135_150", "Interval_150_165", 
                           "Interval_165_180") 
insulin_rate_ogtt = insulin_rate 
dim(insulin_rate_ogtt)
insulin_rate_ogtt = subset(insulin_rate_ogtt, rownames(insulin_rate_ogtt) != "43883-018")
dim(insulin_rate_ogtt)
rownames(insulin_rate_ogtt) = gsub("43883-0","S", rownames(insulin_rate_ogtt))


###### Analysis of Insulin Secretion curves from IIGI ########
iigi_insulin_secretion_files = list.files("data/isec_results/iigi/", recursive=TRUE, pattern = "RESULTS.TXT$")
insulin_rate = data.frame()
s_name = vector()
for(f in iigi_insulin_secretion_files){
  xx = read.table(paste("data/isec_results/iigi/", f, sep = ""), fill = TRUE)
  s =  strsplit(f, "/")[[1]][1]
  s_name = c(s_name, s)
  insulin_rate_tmp = c(as.numeric(as.character(xx[11,2])), as.numeric(as.character(xx[12:23,3])))
  insulin_rate = rbind(insulin_rate, insulin_rate_tmp)
}
rownames(insulin_rate) = paste("43883-", s_name, sep = "")
colnames(insulin_rate) = c("before_0", "Interval_0_15","Interval_15_30", "Interval_30_45", "Interval_45_60", "Interval_60_75", "Interval_75_90", 
                           "Interval_90_105","Interval_105_120", "Interval_120_135", "Interval_135_150", "Interval_150_165", 
                           "Interval_165_180")
insulin_rate_iigi = insulin_rate
dim(insulin_rate_iigi)
insulin_rate_iigi = subset(insulin_rate_iigi, rownames(insulin_rate_iigi) != "43883-018")
dim(insulin_rate_iigi)
rownames(insulin_rate_iigi) = gsub("43883-0","S", rownames(insulin_rate_iigi))


### Visualize OGTT insulin secretion
ogtt_insulin_rate_clinical_demographics_annotation = merge(insulin_rate_iigi, demographics, by.x = "row.names", by.y = "SubjectID", all.x = TRUE)
rownames(ogtt_insulin_rate_clinical_demographics_annotation) = ogtt_insulin_rate_clinical_demographics_annotation$Row.names
ogtt_insulin_rate_clinical_demographics_annotation = ogtt_insulin_rate_clinical_demographics_annotation[,-1]
write.csv(ogtt_insulin_rate_clinical_demographics_annotation, 
          file = "insulin_rate_clinical_demographics_annotation_OGTT_09102021.csv")
my_sample_row <- data.frame(SSPG = ogtt_insulin_rate_clinical_demographics_annotation$sspg,
                            MuscleIR = ogtt_insulin_rate_clinical_demographics_annotation$SSPG_class_2,
                            A1C = ogtt_insulin_rate_clinical_demographics_annotation$a1c_t2d_status,
                            OGTT_2h = ogtt_insulin_rate_clinical_demographics_annotation$ogtt_t2d_status,
                            FPG = ogtt_insulin_rate_clinical_demographics_annotation$fpg_t2d_status)
rownames(my_sample_row) <- rownames(ogtt_insulin_rate_clinical_demographics_annotation)
dev.off()
pdf("heatmap_insulinSecretion_OGTT_10272023.pdf", w = 10, h = 10)
pheatmap(ogtt_insulin_rate_clinical_demographics_annotation[,2:13], 
         color = viridis(25), 
         cluster_cols = FALSE, 
         angle_col = 45,
         annotation_row = my_sample_row, 
         main = "Insulin Secretion based on C-peptide from OGTT")
dev.off()


######## Compare Insulin Secretion from IIGI and OGTT
insulin_rate_ogtt_melted = reshape2::melt(as.matrix(insulin_rate_ogtt))
insulin_rate_ogtt_melted$test = "OGTT"
insulin_rate_iigi_melted = reshape2::melt(as.matrix(insulin_rate_iigi))
insulin_rate_iigi_melted$test = "IIGI" 
insulin_rate_combined = rbind(insulin_rate_ogtt_melted, insulin_rate_iigi_melted)
colnames(insulin_rate_combined) = c("SubjectID", "interval", "insulin_secretion", "test")
insulin_rate_combined$interval = as.character(insulin_rate_combined$interval)
insulin_rate_combined$interval[which(insulin_rate_combined$interval == "before_0")] = "before_before_0"
xx = strsplit(as.character(insulin_rate_combined$interval), "_")
insulin_rate_combined$interval = as.numeric(unlist(lapply(xx, '[[', 3)))

insulin_rate_combined$SubjectID = factor(insulin_rate_combined$SubjectID, levels=a1c_order)
gp = ggplot(insulin_rate_combined, aes(interval, insulin_secretion, col = test)) + 
  geom_line(aes(group=test), alpha=0.8, size = 1.5) + 
  facet_wrap(~SubjectID, nrow = 4) +
  theme_bw() +
  theme(legend.position="top") + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=45, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  labs(y = "Insulin Secretion (pmol/kg/min)", x = "Time (mins)") +
  ggtitle("Insulin Secretion (OGTT vs IIGI)")
gp
ggsave('Insulin_Secretion_OGTT_IIGI_09102021.pdf', h=6, w=10) 


##### Calculate InsulinSecreation/SSPG
ogtt_insulinsecretion_over_ssp = insulin_rate_ogtt 
for(sub in rownames(ogtt_insulinsecretion_over_ssp)){
  ogtt_insulinsecretion_over_ssp[sub,] = ogtt_insulinsecretion_over_ssp[sub,]/demographics$sspg[demographics$SubjectID == sub]
}

### Calculate Disposition Index (DI)
### TODO: Review Units of Disposition Index
timepoints = seq(from = 0, to = 30, by = 15)
disposition_index = apply(ogtt_insulinsecretion_over_ssp[, c("before_0", "Interval_0_15","Interval_15_30")], 1, function(x){
  MESS::auc(timepoints, x)
})
disposition_index = as.data.frame(disposition_index)
colnames(disposition_index) = "DI"
disposition_index$SubjectID = rownames(disposition_index)


## Thresholding Beta-Cell
summary(disposition_index$DI)
disposition_index$di_3classes = "Unknown"
disposition_index$di_3classes[which(disposition_index$DI>2.2)] = "Normal"
disposition_index$di_3classes[which(disposition_index$DI<1.2)] = "Dysfunction"
disposition_index$di_3classes[which(disposition_index$DI>=1.2 & disposition_index$DI <=2.2)] = "Intermediate"

## Save to file
write.csv(disposition_index, file = "disposition_index_09102021.csv")


### Visualize InsulinSecreation/SSPG
disposition_index_clinical_demographics_annotation = merge(ogtt_insulinsecretion_over_ssp, demographics, by.x = "row.names", by.y = "SubjectID", all.x = TRUE)
dim(disposition_index_clinical_demographics_annotation)
disposition_index_clinical_demographics_annotation$SubjectID = disposition_index_clinical_demographics_annotation$Row.names

rownames(disposition_index_clinical_demographics_annotation) = disposition_index_clinical_demographics_annotation$SubjectID
disposition_index_clinical_demographics_annotation = disposition_index_clinical_demographics_annotation[,-1]
disposition_index_clinical_demographics_annotation = merge(disposition_index_clinical_demographics_annotation, disposition_index, by = "row.names", all.x = TRUE)
rownames(disposition_index_clinical_demographics_annotation) = disposition_index_clinical_demographics_annotation$Row.names
disposition_index_clinical_demographics_annotation = disposition_index_clinical_demographics_annotation[,-1]
write.csv(disposition_index_clinical_demographics_annotation, 
          file = "disposition_index_clinical_demographics_annotation_OGTT_09102021.csv")
my_sample_row <- data.frame(DI = disposition_index_clinical_demographics_annotation$DI,
                            MuscleIR = disposition_index_clinical_demographics_annotation$SSPG_class_2,
                            BetaCell = disposition_index_clinical_demographics_annotation$di_3classes,
                            A1C = disposition_index_clinical_demographics_annotation$a1c_t2d_status)
rownames(my_sample_row) <- rownames(disposition_index_clinical_demographics_annotation)
dev.off()
pdf("InsulinSecretionRate_over_SSPG_OGTT_09102021.pdf", w = 10, h = 10)
pheatmap(disposition_index_clinical_demographics_annotation[,2:13], color = viridis(10), 
         cluster_cols = FALSE, annotation_row = my_sample_row, main = "InsulinSecretionRate/SSPG from OGTT", angle_col = 45)
dev.off()



### Visualize 8 selected subjects InsulinSecreation/SSPG
disposition_index_clinical_demographics_annotation = merge(ogtt_insulinsecretion_over_ssp, demographics, by.x = "row.names", by.y = "SubjectID", all.x = TRUE)
disposition_index_clinical_demographics_annotation$SubjectID = disposition_index_clinical_demographics_annotation$Row.names
selected_subjects = c("S49", "S14", "S27", "S42", "S58", "S53", "S28","S57","S45")
disposition_index_clinical_demographics_annotation = disposition_index_clinical_demographics_annotation[which(disposition_index_clinical_demographics_annotation$SubjectID %in% selected_subjects),]
disposition_index_clinical_demographics_annotation$SubjectID = factor(disposition_index_clinical_demographics_annotation$SubjectID, levels=a1c_order)

rownames(disposition_index_clinical_demographics_annotation) = disposition_index_clinical_demographics_annotation$Row.names
disposition_index_clinical_demographics_annotation = disposition_index_clinical_demographics_annotation[,-1]
disposition_index_clinical_demographics_annotation = merge(disposition_index_clinical_demographics_annotation, disposition_index, by = "row.names", all.x = TRUE)
rownames(disposition_index_clinical_demographics_annotation) = disposition_index_clinical_demographics_annotation$Row.names
disposition_index_clinical_demographics_annotation = disposition_index_clinical_demographics_annotation[,-1]

my_sample_row <- data.frame(#mod_DI = disposition_index_clinical_demographics_annotation$DI,
                            Muscle_IR = disposition_index_clinical_demographics_annotation$SSPG_class,
                            BetaCell = disposition_index_clinical_demographics_annotation$di_3classes,
                            A1C = disposition_index_clinical_demographics_annotation$a1c_t2d_status)
rownames(my_sample_row) <- rownames(disposition_index_clinical_demographics_annotation)
dev.off()
pdf("InsulinSecretionRate_over_SSPG_OGTT_selected_09302021.pdf", w = 8, h = 4)
pheatmap(disposition_index_clinical_demographics_annotation[,2:13], 
         color = viridis(10), 
         cluster_cols = FALSE, 
         annotation_row = my_sample_row, 
         fontsize = 12,
         main = "", 
         angle_col = 45)
dev.off()


#############################
##### Hepatic IR     ########
#############################
#### Visualize insulin from OGTT and IIGI
ogtt_insulin = ogtt_analytes[,c("SubjectID", "TimePoint","Insulin")]
dim(ogtt_insulin)
ogtt_insulin = na.omit(ogtt_insulin)
dim(ogtt_insulin)
length(unique(ogtt_insulin$SubjectID))
ogtt_insulin$test = "OGTT"
iigi_insulin = iigi_analytes[,c("SubjectID", "TimePoint","Insulin")]
dim(iigi_insulin)
iigi_insulin = na.omit(iigi_insulin)
dim(iigi_insulin)
length(unique(iigi_insulin$SubjectID))
iigi_insulin$test = "IIGI"
ogtt_iigi_insulin_combined = rbind(ogtt_insulin, iigi_insulin)

ogtt_iigi_insulin_combined$SubjectID = factor(ogtt_iigi_insulin_combined$SubjectID, levels=a1c_order)
gp = ggplot(ogtt_iigi_insulin_combined, aes(TimePoint, Insulin, col = test)) + 
  geom_line(aes(group=test), alpha=0.8, size = 1.5) + 
  geom_point(aes(group=test), alpha=0.8, size = 1.5)  + 
  facet_wrap(~SubjectID, nrow = 4) +
  theme_bw() +
  theme(legend.position="top") + 
  scale_color_brewer(palette="Dark2") + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=45, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  labs(y = "Insulin (microIU/mL)" , x = "Time (mins)") +
  ggtitle("Insulin (OGTT vs IIGI)")
gp
ggsave('figures/Insulin_OGTT_IIGI_02222023.pdf', h=6, w=10) 


### Calculate AUC for Insulin OGTT 0-120
ogtt_insulin_dcast = reshape2::dcast(ogtt_insulin, SubjectID ~ TimePoint, value.var = "Insulin")
rownames(ogtt_insulin_dcast) = ogtt_insulin_dcast$SubjectID
ogtt_insulin_dcast = ogtt_insulin_dcast[,-1]
ogtt_insulin_dcast_imputed = t(apply(ogtt_insulin_dcast, 1, na_interpolation))
write.csv(ogtt_insulin_dcast_imputed, file = "ogtt_insulin_dcast_imputed_02282021.csv")
timepoints = c(0, 15, 30, 60, 90, 120)
insulin_ogtt_auc_0_120 = apply(ogtt_insulin_dcast_imputed[,c("0", "15", "30", "60", "90", "120")], 1, function(x){
  MESS::auc(timepoints, x)
})
insulin_ogtt_auc_0_120 = as.data.frame(insulin_ogtt_auc_0_120)


### Calculate body fat percentage (BFP%)
demographics$sex_01 = NA
demographics$sex_01[demographics$sex == "M"] = 1
demographics$sex_01[demographics$sex == "F"] = 0
demographics$bf_percentage = 1.2*demographics$bmi + 0.23*demographics$age - 10.8*demographics$sex_01 - 5.4
demographics_insulinAUC = merge(demographics, insulin_ogtt_auc_0_120, by.x = "SubjectID", by.y = "row.names", all.x  = TRUE)

### metabolic panel (hdl)
hdl = subset(metabolic_df, blood == "hdl")
hdl = hdl[,c("SubjectID", "value")]
colnames(hdl) = c("SubjectID", "hdl")
hdl$hdl = as.numeric(as.character(hdl$hdl))
hdl = ddply(hdl, ~SubjectID, summarise, hdl = mean(hdl))



### Calculate Hepatic IR
hepatic_ir_df = merge(demographics_insulinAUC, hdl, by = "SubjectID", all.x  = TRUE)
hepatic_ir_df$hepatic_IR = -0.091 + 0.4*log(hepatic_ir_df$insulin_ogtt_auc_0_120) +
  0.346*log(hepatic_ir_df$bf_percentage) - 0.408*log(hepatic_ir_df$hdl) +
  0.435*log(hepatic_ir_df$bmi)
summary(hepatic_ir_df$hepatic_IR)
dim(hepatic_ir_df)
hepatic_ir_df = hepatic_ir_df[which(!is.na(hepatic_ir_df$hepatic_IR)),]
dim(hepatic_ir_df)
write.csv(hepatic_ir_df, file = "demographics_hepatic_IR_09102021.csv")

### Thresholding Hepatic-IR
hepatic_ir_df$hepatic_ir_3classes = "Unknown"
hepatic_ir_df$hepatic_ir_3classes[which(hepatic_ir_df$hepatic_IR<3.95)] = "IS"
hepatic_ir_df$hepatic_ir_3classes[which(hepatic_ir_df$hepatic_IR>4.8)] = "IR"
hepatic_ir_df$hepatic_ir_3classes[which(hepatic_ir_df$hepatic_IR>=3.95 & hepatic_ir_df$hepatic_IR<=4.8)] = "Intermediate"
table(hepatic_ir_df$hepatic_ir_3classes)


#############################################
######## Aggregate Metabolic indicator ######
#############################################
### reformat each metabolic indicator df to be able to merge them later
hepatic_ir_df = hepatic_ir_df
ie_df = ie_df
disposition_index_df = disposition_index
homa_matsuda_df = metsuda_df[, c("SubjectID", "HOMA_B", "HOMA_IR", "HOMA_S", "HOMA_IR_group" , "matsuda_index")]

## Aggregate metabolic indicators
aggregated_metabolic_indicators = merge(hepatic_ir_df, ie_df, by = "SubjectID")
aggregated_metabolic_indicators = merge(aggregated_metabolic_indicators, disposition_index_df, by = "SubjectID")
aggregated_metabolic_indicators = merge(aggregated_metabolic_indicators, homa_matsuda_df, by = "SubjectID")
dim(aggregated_metabolic_indicators)

#### Save aggregated metabolic indicators
write.csv(aggregated_metabolic_indicators, "aggregated_metabolic_indicators_09102021.csv", row.names = FALSE)



########################
#######  Figure 2: Metabolic Subphenotyping levels
########################
aggregated_metabolic_indicators = read.csv("../aggregated_metabolic_indicators__4_predictions_w_prs_01042022.csv")

metabolic_indicators_orderd = aggregated_metabolic_indicators[order(aggregated_metabolic_indicators$a1c, 
                                                         decreasing = TRUE),]
metabolic_indicators_orderd$SubjectID <- factor(metabolic_indicators_orderd$SubjectID, 
                                                levels = metabolic_indicators_orderd$SubjectID) ## This to lock the order

#### Plot A1C levels
# group.color = c(Normal = "#00AFBB", PreDM = "#E7B800", T2D ="#FC4E07")
group.color = c(Normoglycemic = "#00AFBB", PreDM = "#E7B800", T2D ="#FC4E07")
gp = ggplot(metabolic_indicators_orderd, aes(x=SubjectID, y=a1c, color  = a1c_t2d_status)) +
  geom_segment( aes(x=SubjectID, xend=SubjectID, y=0, yend=a1c, color  = a1c_t2d_status))+
  geom_point(size=3) +
  coord_flip() +
  scale_fill_manual(values = group.color) + 
  scale_color_manual(values = group.color) + 
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="top",
    legend.title = element_blank()
  ) +
  theme(axis.text.x = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("HbA1c") +
  xlab("") +
  ylab("HbA1c %")
gp
ggsave('metabolicindicator_a1c_10272023.pdf', h=7, w=3.5); 


### Plot SSPG levels
group.color = c(IR = "#333BFF", IS = "#CC6600")
metabolic_indicators_orderd$SSPG_class_2 <- factor(metabolic_indicators_orderd$SSPG_class_2, levels = c("IS", "IR"))
gp = ggplot(metabolic_indicators_orderd, aes(x=SubjectID, y=sspg, color  = SSPG_class_2)) +
  geom_segment( aes(x=SubjectID, xend=SubjectID, y=0, yend=sspg, color  = SSPG_class_2))+
  geom_point(size=3) +
  coord_flip() +
  scale_fill_manual(values = group.color) + 
  scale_color_manual(values = group.color) + 
  theme_light() +
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="top",
    legend.title = element_blank()
  ) +
  theme(axis.text.x = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Muscle IR") +
  xlab("") +
  ylab("SSPG mg/dL")
gp
ggsave('metabolicindicator_sspg_09102021.pdf', h=7, w=3.5); 

## Plot DI levels
group.color = c(Normal = "green", Dysfunction = "red", Intermediate = "blue", Unknown ="grey86")
metabolic_indicators_orderd$di_3classes <- factor(metabolic_indicators_orderd$di_3classes, 
                                                  levels = c("Normal", "Intermediate", "Dysfunction"))
gp = ggplot(metabolic_indicators_orderd, aes(x=SubjectID, y=DI, color  = di_3classes)) +
  geom_segment( aes(x=SubjectID, xend=SubjectID, y=0, yend=DI, color  = di_3classes))+
  geom_point(size=3) +
  coord_flip() +
  scale_fill_manual(values = group.color) + 
  scale_color_manual(values = group.color) + 
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="top",
    legend.title = element_blank()
  ) +
  theme(axis.text.x = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Beta-cell Function") +
  xlab("") +
  ylab("Disposition Index")
  # ylab("Disposition Index (pmol*dL)/(kg*ml)")
gp
ggsave('metabolicindicator_di_09102021.pdf', h=7, w=3.5)



## Plot IE levels
group.color = c(Normal = "#66c2a5", Dysfunction = "#fc8d62", Intermediate = "#984ea3",
                Unknown ="grey86") 
metabolic_indicators_orderd$ie_3_classes <- factor(metabolic_indicators_orderd$ie_3_classes, 
                                                   levels = c("Normal", "Intermediate", "Dysfunction"))
gp = ggplot(metabolic_indicators_orderd, aes(x=SubjectID, y=ie, color  = ie_3_classes)) +
  geom_segment( aes(x=SubjectID, xend=SubjectID, y=0, yend=ie, color  = ie_3_classes))+
  geom_point(size=3) +
  coord_flip() +
  scale_fill_manual(values = group.color) + 
  scale_color_manual(values = group.color) + 
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="top",
    legend.title = element_blank()
  ) +
  theme(axis.text.x = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Incretin Effect") +
  xlab("") +
  ylab("Incretin Effect %")
gp
ggsave('metabolicindicator_ie_09102021.pdf', h=7, w=3.5); 



## Plot Hepatic IR index levels
group.color = c(IS = "#3288BD", Intermediate = "gray", IR = "#D53E4F")
# group.color = c(IS = "#00AFBB", Intermediate = "#E7B800", IR = "#FC4E07")
metabolic_indicators_orderd$hepatic_ir_3classes <- factor(metabolic_indicators_orderd$hepatic_ir_3classes, 
                                                          levels = c("IS", "Intermediate", "IR"))
gp = ggplot(metabolic_indicators_orderd, aes(x=SubjectID, y=hepatic_IR, color  = hepatic_ir_3classes)) +
  geom_segment( aes(x=SubjectID, xend=SubjectID, y=3, yend=hepatic_IR, color  = hepatic_ir_3classes))+
  geom_point(size=3) +
  coord_flip() +
  scale_fill_manual(values = group.color) + 
  scale_color_manual(values = group.color) + 
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="top",
    legend.title = element_blank()
  ) +
  
  # scale_x_continuous(limits=c(2, 5.5)) + 
  theme(axis.text.x = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
  ggtitle("Hepatic IR") +
  ylim(3,5.25) + 
  xlab("") +
  ylab("Hepatic IR Index")
gp
ggsave("/Users/ahmedm/Library/CloudStorage/Box-Box/Ahmed Metwally's Files/Stanford/cgm/43883/Manuscripts/Manuscript_1_MetabolicSubphenotyping/metabolicindicator_hepaticIR_07102022.pdf", h=7, w=3.5)



########################
#### Classes Histogram
########################
### MuscleIR
table(metabolic_indicators_orderd$SSPG_class_2)
group.color = c(IR = "#333BFF", IS = "#CC6600")
p = ggplot(data = metabolic_indicators_orderd) +
  geom_histogram(aes(x=sspg, fill=SSPG_class_2), color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = group.color) + 
  theme_bw() +
  theme(legend.position="top") +
  ggtitle("Muscle IR") +
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  xlab("SSPG (mg/dL)") +
  ylab("Count") +
  labs(fill="")
p
ggsave('metabolicindicator_histogram_sspg_09102021.pdf', h=3, w=5)

## DI
table(metabolic_indicators_orderd$di_3classes)
metabolic_indicators_orderd$di_3classes <- factor(metabolic_indicators_orderd$di_3classes, 
                                                  levels = c("Normal", "Intermediate", "Dysfunction"))
group.color = c(Normal = "green", Dysfunction = "red", Intermediate = "blue", Unknown ="grey86")
p = ggplot(data = metabolic_indicators_orderd) +
  geom_histogram(aes(x=DI, fill=di_3classes), color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = group.color) + 
  theme_bw() +
  theme(legend.position="top") +
  ggtitle("Beta-cell Dysfunction") +
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  xlab("Disposition Index (pmol*dL)/(kg*ml)") +
  ylab("Count") +
  labs(fill="")
p
ggsave('metabolicindicator_histogram_modi_12202021.pdf', h=3, w=5)



### IE
table(metabolic_indicators_orderd$ie_3_classes)
metabolic_indicators_orderd$ie_3_classes <- factor(metabolic_indicators_orderd$ie_3_classes, 
                                                   levels = c("Normal", "Intermediate", "Dysfunction"))
group.color = c(Normal = "#66c2a5", Dysfunction = "#fc8d62", Intermediate = "#984ea3",
                Unknown ="grey86")   
p = ggplot(data = metabolic_indicators_orderd) +
  geom_histogram(aes(x=ie, fill=ie_3_classes), color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = group.color) + 
  theme_bw() +
  theme(legend.position="top") +
  ggtitle("Impaired Incretin Effect") +
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  xlab("Incretin Effect %") +
  ylab("Count") +
  labs(fill="")
p
ggsave('metabolicindicator_histogram_ie_12202021.pdf', h=3, w=5)


#### HepaticIR
table(metabolic_indicators_orderd$hepatic_ir_3classes)
metabolic_indicators_orderd$hepatic_ir_3classes <- factor(metabolic_indicators_orderd$hepatic_ir_3classes, 
                                                          levels = c("IS", "Intermediate", "IR"))

group.color = c(IS = "#3288BD", Intermediate = "gray", IR = "#D53E4F")
p = ggplot(data = metabolic_indicators_orderd,) +
  geom_histogram(aes(x=hepatic_IR, fill= hepatic_ir_3classes_v2), color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values = group.color) + 
  theme_bw() +
  theme(legend.position="top") +
  ggtitle("Hepatic IR") +
  theme(axis.text.x = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  xlab("Hepatic IR Index") +
  ylab("Count") +
  labs(fill="")
p
ggsave('metabolicindicator_histogram_hepatic_12202021.pdf', h=3, w=5)




############################################################################# 
############## Fig 2: Determining Dominant Metabolic phenotype ##############
#############################################################################
## Calculate deviance from average cohort
metabolic_indicators_orderd$MuscleIR_Deviance = (metabolic_indicators_orderd$sspg - mean(metabolic_indicators_orderd$sspg))/ sd(metabolic_indicators_orderd$sspg)
metabolic_indicators_orderd$BetaCell_Deviance = -(metabolic_indicators_orderd$DI - mean(metabolic_indicators_orderd$DI))/ sd(metabolic_indicators_orderd$DI)
metabolic_indicators_orderd$Incretin_Deviance = -(metabolic_indicators_orderd$ie - mean(metabolic_indicators_orderd$ie))/ sd(metabolic_indicators_orderd$ie)
metabolic_indicators_orderd$HepaticIR_Deviance = (metabolic_indicators_orderd$hepatic_IR - mean(metabolic_indicators_orderd$hepatic_IR))/ sd(metabolic_indicators_orderd$hepatic_IR)
metabolic_indicators_orderd$Total_MSP_Deviance = rowSums(metabolic_indicators_orderd[, c("MuscleIR_Deviance", "BetaCell_Deviance", 
                                                                  "Incretin_Deviance", "HepaticIR_Deviance")], na.rm = TRUE)
rownames(metabolic_indicators_orderd) = metabolic_indicators_orderd$SubjectID
write.csv(metabolic_indicators_orderd, file = "metabolic_indicators_risk_score_10272023.csv")


#### Heatmap based on deviance of metabolic phenotypes
dev.off()
pdf("DominantMetabolicPhenotype_with_Total_MSP_Deviance_10272023.pdf", width = 5, height = 7)
pheatmap(metabolic_indicators_orderd[, c("MuscleIR_Deviance", "BetaCell_Deviance", 
                                  "Incretin_Deviance", "HepaticIR_Deviance", "Total_MSP_Deviance")],
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         color = RColorBrewer::brewer.pal(100, "YlOrRd"), #RColorBrewer::brewer.pal(4, "Purples"), #myColor,
         fontsize_row = 10,
         fontsize_col = 12,
         main = "",
         angle_col = 45
)
dev.off()


###### summarize metabolic subphenotype
apply(metabolic_indicators_orderd[, c("sspg", "ie", "DI", "hepatic_IR", "a1c", "fpg", "ogtt_2h")], 2, summary)

### Correlation between the deviance of metabolic subphenotypes
metabolic_indicators_tstat = metabolic_indicators_orderd[, c("MuscleIR_Deviance", "BetaCell_Deviance", 
                                  "Incretin_Deviance", "HepaticIR_Deviance")]
pdf("DominantMetabolicPhenotype_correlation_09102021.pdf", width = 6, height = 6)
metabolic_indicators_tstat %>% correlate() %>% network_plot(min_cor=0.1, 
                                                           repel = TRUE,
                                                           legend = TRUE,
                                                           curved = TRUE,
                                                           colors = c("#00AFBB", "white", "#FC4E07"))#c("blue", "white", "red"))
dev.off()


pdf("DominantMetabolicPhenotype_correlation_pvalue_09102021.pdf", width = 7, height = 7)
chart.Correlation(metabolic_indicators_tstat, histogram=TRUE, pch=19)
dev.off()

rcorr(as.matrix(metabolic_indicators_tstat))


################################################
############# Preprocess OGTT timeseries
################################################
dt = as.data.frame(ogtt)
dt = dt[,c("SubjectID", "ogtt_min", "ogtt_value")]

######## Impute OGTT missing timepoints 
dt = reshape2::dcast(dt, ogtt_min ~ SubjectID)
rownames(dt) = dt$ogtt_min
dt = dt[,-1]
ogtt_dt_imputed = apply(dt, 2, na_interpolation)


### Calculate number of missing data per subjects
sum(is.na(dt$S21))



na_count <-sapply(dt, function(y) sum(length(which(is.na(y)))))
na_count
barplot(na_count, main="Car Distribution",
        xlab="Number of Gears")
gp

names(na_count)
ee = as.data.frame(na_count)
ee$subject = names(na_count)

gp = ggplot(ee, aes(x=subject, y=na_count)) +
  geom_bar(stat = "identity") + 
  # scale_fill_manual(values=cbPalette) + 
  theme_bw() +
  coord_cartesian(ylim = c(0, 3)) + 
  theme(legend.position="top") + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=90, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"), #legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 14, color = "black", face = "bold")) +
  labs(y = "# of missing timepoints", x = "") +
  ggtitle("Missing glucose timepoints during OGTT")

gp
ggsave(sprintf('/Users/ahmedm/Dropbox/Diabetes_v4_resubmission/Metabolic_Subphenotype_Predictor_07102022/figures/Missing_data.pdf'), h=3, w=8)

######## Normalize imputed OGTT timeseries 
znorm <- function(ts){
  ts.mean <- mean(ts)
  ts.dev <- sd(ts)
  (ts - ts.mean)/ts.dev
}
ogtt_imputed = as.data.frame(ogtt_dt_imputed)
ogtt_imputed_normalized = apply(ogtt_imputed, 2, znorm)
write.csv(ogtt_imputed_normalized, "ogtt_imputed_normalized_09102021.csv")
ogtt_imputed_normalized_melted = reshape2::melt(as.matrix(ogtt_imputed_normalized))
names(ogtt_imputed_normalized_melted) = c("ogtt_min", "SubjectID", "ogtt_value")


######## Smooth normalized/imputed OGTT timeseries
ogtt_imputed_normalized_smoothed = apply(ogtt_imputed_normalized, 2, function(x){
  smoothingSpline = smooth.spline(x, spar=0.35)
  smoothingSpline$y
})
rownames(ogtt_imputed_normalized_smoothed) = rownames(ogtt_imputed_normalized)
ogtt_imputed_normalized_smoothed_melted = reshape2::melt(as.matrix(ogtt_imputed_normalized_smoothed))
names(ogtt_imputed_normalized_smoothed_melted) = c("ogtt_min", "SubjectID", "ogtt_value")


#### Visualized smoothed vs nonsmoothes curves
ogtt_imputed_normalized_melted
ogtt_imputed_normalized_smoothed_melted


tmp_normalized = ogtt_imputed_normalized_melted
tmp_normalized$type = "non_smoothed"

tmp_normalized_smoothed = ogtt_imputed_normalized_smoothed_melted
tmp_normalized_smoothed$type = "smoothed"

tmp_smoothed_nonsmoothed = rbind(tmp_normalized, tmp_normalized_smoothed)

gp = ggplot(tmp_smoothed_nonsmoothed, aes(ogtt_min, ogtt_value, col = type)) + 
  geom_line(aes(group=type), alpha=0.8, size = 1.5) + 
  scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999")) + 
  facet_wrap(~SubjectID, nrow = 4) +
  theme_bw() +
  theme(legend.position="top") + 
  theme(axis.text.x = element_text(colour="black", size=10, angle=45, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=10, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  labs(y = "Z-normalized glucose concentration", x = "Time (mins)") +
  ggtitle("Z-normalized smoothed vs non-smoothed OGTT glucose timeseries")
gp
ggsave('figures/smoothed_nonsmoothed.pdf', h=6, w=10)



#### Order OGTT curves based on a1c and choose the 32 subjects
dem_a1c = demographics[, c("SubjectID", "a1c")]
a1c_myList <- as.list(dem_a1c$a1c)
names(a1c_myList) <- dem_a1c$SubjectID
a1c_myList = unlist(a1c_myList)
ogtt_order = names(sort(a1c_myList, decreasing = FALSE))


#### Plot imputed OGTT curves Ordered based on a1c 
ogtt_imputed = reshape2::melt(as.matrix(ogtt_dt_imputed))
names(ogtt_imputed) = c("ogtt_min", "SubjectID", "ogtt_value")
ogtt_imputed_33 = ogtt_imputed[which(ogtt_imputed$SubjectID %in% dem_a1c$SubjectID),]
ogtt_imputed_33$SubjectID = factor(ogtt_imputed_33$SubjectID, levels=ogtt_order)
ogtt_imputed_33_a1c_group = merge(ogtt_imputed_33, demographics, by = "SubjectID")

# group.color = c(Normal = "#00AFBB", PreDM = "#E7B800", T2D ="#FC4E07")
group.color = c(Normoglycemic = "#00AFBB", PreDM = "#E7B800", T2D ="#FC4E07")
gp = ggplot(ogtt_imputed_33_a1c_group, aes(ogtt_min, ogtt_value)) + 
  geom_line(aes(group=SubjectID, color = a1c_t2d_status), alpha=0.8, size = 1.5) + 
  geom_hline(yintercept=50, linetype="dashed", color = "yellow", alpha = 0.6) +
  geom_hline(yintercept=140, linetype="dashed", color = "orange", alpha = 0.6) +  
  geom_hline(yintercept=200, linetype="dashed", color = "red", alpha = 0.6) +
  geom_vline(xintercept=120, linetype="dashed", alpha = 0.6) +  
  facet_wrap(~SubjectID, nrow = 4) +
  scale_fill_manual(values = group.color) + 
  scale_color_manual(values = group.color) + 
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
ggsave('ogtt_imputed_ordered_colored_10272023.pdf', h=6, w=10); 



#### Plot imputed OGTT curves Ordered based on a1c 
ogtt_imputed = reshape2::melt(as.matrix(ogtt_dt_imputed))
names(ogtt_imputed) = c("ogtt_min", "SubjectID", "ogtt_value")
ogtt_imputed_33 = ogtt_imputed[which(ogtt_imputed$SubjectID %in% dem_a1c$SubjectID),]
ogtt_imputed_33$SubjectID = factor(ogtt_imputed_33$SubjectID, levels=ogtt_order)
gp = ggplot(ogtt_imputed_33, aes(ogtt_min, ogtt_value)) + 
  geom_line(aes(group=SubjectID), alpha=0.8, size = 1.5) + 
  geom_hline(yintercept=50, linetype="dashed", color = "yellow", alpha = 0.6) +
  geom_hline(yintercept=140, linetype="dashed", color = "orange", alpha = 0.6) +  
  geom_hline(yintercept=200, linetype="dashed", color = "red", alpha = 0.6) +
  geom_vline(xintercept=120, linetype="dashed", alpha = 0.6) +  
  facet_wrap(~SubjectID, nrow = 4) +
  theme_bw() +
  theme(legend.position="top") + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=45, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  labs(y = "Glucose (mg/dL)", x = "Time relative to the start of drinking Glucola (mins)") +
  ggtitle("Glucose Concentration during Extended OGTT (ordered by A1C)")
gp
ggsave('ogtt_imputed_ordered_09102020.pdf', h=6, w=10); 


#### Plot normalized and imputed OGTT curves Ordered based on a1c and choose the 32 subjects
ogtt_imputed_normalized_melted_33 = ogtt_imputed_normalized_melted[which(ogtt_imputed_normalized_melted$SubjectID %in% dem_a1c$SubjectID),]
ogtt_imputed_normalized_melted_33$SubjectID = factor(ogtt_imputed_normalized_melted_33$SubjectID, levels=ogtt_order)
gp = ggplot(ogtt_imputed_normalized_melted_33, aes(ogtt_min, ogtt_value)) + 
  geom_line(aes(group=SubjectID), alpha=0.8, size = 1.5) + 
  facet_wrap(~SubjectID, nrow = 4) +
  theme_bw() +
  theme(legend.position="top") + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=45, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  labs(y = "Z-normalized Glucose", x = "Time relative to the start of drinking Glucola (mins)") +
  ggtitle("Normalized and Imputed Glucose Concentration from Extended OGTT (ordered by A1C)")
gp
ggsave('ogtt_imputed_normalized_ordered_09102020.pdf', h=6, w=10); 


### Plot imputed, normalized and smoothed data and ordered OGTT curves 
ogtt_imputed_normalized_smoothed_melted$SubjectID = factor(ogtt_imputed_normalized_smoothed_melted$SubjectID, levels=ogtt_order)
gp = ggplot(ogtt_imputed_normalized_smoothed_melted, aes(ogtt_min, ogtt_value)) + 
  geom_line(aes(group=SubjectID), alpha=0.8, size = 1.5) + 
  facet_wrap(~SubjectID, nrow = 4) +
  theme_bw() +
  theme(legend.position="top") + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=45, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  labs(y = "Z-normalized Glucose", x = "Time relative to the start of drinking Glucola (mins)") +
  ggtitle("Normalized, Imputed, and Smoothed Extended OGTT (ordered by A1C)")
gp
ggsave('ogtt_imputed_normalized_smoothed_ordered_09102020.pdf', h=6, w=10); 


### IIGI QC: Remove low quality curves (remove subject 18 because of missing data beyond 75 min) & select 32 subjects
dt = as.data.frame(iigi)
dt = dt[,c("SubjectID", "iigi_min", "iigi_value")]
dt = reshape2::dcast(dt, iigi_min ~ SubjectID)
rownames(dt) = dt$iigi_min
dt = dt[,-1]
iigi_dt_imputed = apply(dt, 2, na_interpolation)
iigi_imputed = reshape2::melt(as.matrix(iigi_dt_imputed))
names(iigi_imputed) = c("iigi_min", "SubjectID", "iigi_value")

iigi_imputed$SubjectID = factor(iigi_imputed$SubjectID, levels=ogtt_order)
gp = ggplot(iigi_imputed, aes(iigi_min, iigi_value)) + 
  geom_line(aes(group=SubjectID), alpha=0.8, size = 1.5, color ="blue") + 
  geom_hline(yintercept=50, linetype="dashed", color = "yellow", alpha = 0.6) +
  geom_hline(yintercept=140, linetype="dashed", color = "orange", alpha = 0.6) +  
  geom_hline(yintercept=200, linetype="dashed", color = "red", alpha = 0.6) +
  geom_vline(xintercept=120, linetype="dashed", alpha = 0.6) +  
  facet_wrap(~SubjectID, nrow = 4) +
  theme_bw() +
  theme(legend.position="top") + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=45, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  labs(y = "Glucose (mg/dL)", x = "Time relative to the start of IIGI (mins)") +
  ggtitle("Extended IIGI")
gp
ggsave('iigi_imputed_09102020.pdf', h=6, w=10); 


######## Plot glucose comparison between OGTT and IIGI
colnames(ogtt_imputed) = c("time", "SubjectID", "value")
colnames(iigi_imputed) = c("time", "SubjectID", "value")
ogtt_imputed$test = "OGTT" 
iigi_imputed$test = "IIGI"
glucose_ogtt_iigi = rbind(ogtt_imputed, iigi_imputed)
dim(glucose_ogtt_iigi)
comon_subj = intersect(unique(ogtt_imputed$SubjectID), unique(iigi_imputed$SubjectID))
glucose_ogtt_iigi = subset(glucose_ogtt_iigi, SubjectID %in% comon_subj)
dim(glucose_ogtt_iigi)

glucose_ogtt_iigi$SubjectID = factor(glucose_ogtt_iigi$SubjectID, levels=ogtt_order)
gp = ggplot(glucose_ogtt_iigi, aes(time, value, col = test)) + 
  geom_point(aes(group=test), alpha=0.8, size = 3) +  
  geom_line(aes(group=test), alpha=0.4, size = 1.5) + 
  facet_wrap(~SubjectID, nrow = 4) +
  theme_bw() +
  theme(legend.position="top") + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=45, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=12, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=15, face="plain"), legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
  labs(y = "Glucose (mg/dL)", x = "Time relative to the start of OGTT or IIGI (mins)") +
  ggtitle("Glucose from OGTT vs IIGI")
gp
ggsave('ogtt_iigi_imputed_12202021.pdf', h=6, w=10); 




######## Correlation between OGTT and IIGT
glucose_ogtt_iigi_dcast = reshape2::dcast(glucose_ogtt_iigi, SubjectID + time ~ test)
subjects = as.character(unique(glucose_ogtt_iigi_dcast$SubjectID))
ogtt_iigi_cor = vector()
ogtt_iigi_dist = vector()
ogtt_iigi_rmse = vector()
ogtt_iigi_diff = vector()
for(subj in subjects){
  # print(subj)
  tmp = subset(glucose_ogtt_iigi_dcast, SubjectID == subj)
  tmp$IIGI = na_interpolation(tmp$IIGI)
  tmp$OGTT = na_interpolation(tmp$OGTT)
  ogtt_iigi_cor = c(ogtt_iigi_cor, cor(tmp$OGTT, tmp$IIGI))
  ogtt_iigi_dist = c(ogtt_iigi_dist, dist(rbind(tmp$OGTT, tmp$IIGI)))
  ogtt_iigi_rmse = c(ogtt_iigi_rmse, rmse(tmp$OGTT, tmp$IIGI))
  ogtt_iigi_diff = c(ogtt_iigi_diff, sum(tmp$OGTT -  tmp$IIGI))
  glucose_ogtt_iigi_dcast[which(glucose_ogtt_iigi_dcast$SubjectID == subj), ] = tmp
}
names(ogtt_iigi_cor) = subjects
names(ogtt_iigi_dist) = subjects
names(ogtt_iigi_rmse) = subjects
names(ogtt_iigi_diff) = subjects
ogtt_iigi_evaluation = cbind(ogtt_iigi_cor, ogtt_iigi_dist, ogtt_iigi_rmse, ogtt_iigi_diff)
write.csv(ogtt_iigi_evaluation, file = "ogtt_iigi_evaluation_09102021.csv")
ogtt_iigi_evaluation_summary = apply(ogtt_iigi_evaluation, 2, summary)
write.csv(ogtt_iigi_evaluation_summary, file = "ogtt_iigi_evaluation_summary_09102021.csv")


## Stat
stat_tmp = as.data.frame(ogtt_iigi_evaluation)

mean(stat_tmp$ogtt_iigi_cor)
sd(stat_tmp$ogtt_iigi_cor)

# ######################
# ### PCA of OGTT TS
# ######################
df_pca = prcomp(t(ogtt_imputed_normalized))
write.csv(df_pca$x, file = "ogtt_df_pca_pcs_09102021.csv")
df_pca_annotation = merge(df_pca$x, metabolic_indicators_orderd, by.x = "row.names", by.y = "SubjectID", all.x = TRUE)
colnames(df_pca_annotation)[which(colnames(df_pca_annotation) == "Row.names")] = "SubjectID"
write.csv(df_pca_annotation, "ogtt_df_pca_annotation_09102021.csv", row.names = FALSE)


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
  ggsave(paste(text, '_12202021.pdf', sep = ""), h=4, w=6)
}


### PCA colored based on Muscle IR
df_pca_annotation$SSPG_class_2 <- factor(df_pca_annotation$SSPG_class_2, levels = c("IS", "IR"))
plot_pca_metabolic_subphenotype(df = df_pca_annotation, subphenotype = "SSPG_class_2", 
                                group_color = c(IR = "#333BFF", IS = "#CC6600"), 
                                text = "Reduced representation of glucose-series colored by Muscle IR")

### PCA colored based on Beta Cell
df_pca_annotation$di_3classes_v2 <- factor(df_pca_annotation$di_3classes_v2, 
                                        levels = c("Normal", "Intermediate", "Dysfunction"))
plot_pca_metabolic_subphenotype(df = df_pca_annotation, subphenotype = "di_3classes_v2",
                                text = "Reduced representation of glucose-series colored by Beta-cell Dysfunction",
                                group_color = c(Normal = "green", Dysfunction = "red", Intermediate = "blue"))

### PCA colored based on IE
df_pca_annotation$ie_3_classes <- factor(df_pca_annotation$ie_3_classes, 
                                         levels = c("Normal", "Intermediate", "Dysfunction"))
plot_pca_metabolic_subphenotype(df = df_pca_annotation, subphenotype = "ie_3_classes", 
                                group_color = c(Normal = "#66c2a5", Dysfunction = "#fc8d62", Intermediate = "#984ea3"),
                                text = "Reduced representation of glucose-series colored by Incretin Effect")

### PCA colored based on Hepatic IR
df_pca_annotation$hepatic_ir_3classes_v2 <- factor(df_pca_annotation$hepatic_ir_3classes_v2, 
                                                levels = c("IS", "Intermediate", "IR", "Unknown"))
plot_pca_metabolic_subphenotype(df = df_pca_annotation, subphenotype = "hepatic_ir_3classes_v2",
                                text = "Reduced representation of glucose-series colored by Hepatic IR",
                                group_color = c(IS = "#3288BD", Intermediate = "gray", IR = "#D53E4F"))


####################################################
### Feature Extraction from the OGTT curves
####################################################
## remove timepoint -10 from imputed ogtt
ogtt_imputed = ogtt_dt_imputed[-1, ]  

###### Total AUC (tAUC)
ogtt_auc = apply(ogtt_imputed, 2, function(x){
  MESS::auc(as.numeric(rownames(ogtt_imputed)), x, from = 0, to = 180,  type = "linear")
})


###### Incremental AUC (iAUC). 
ogtt_iauc = ogtt_auc - 180*ogtt_imputed["0", ]


###### pAUC: positive AUC
qq = apply(ogtt_imputed, 2, function(x){
  x-x["0"]
})
qq[which(qq<0)]=0
ogtt_pauc = apply(qq, 2, function(x){
  MESS::auc(as.numeric(rownames(qq)), x, from = 0, to = 180,  type = "linear")
})


###### nAUC: negative AUC
ogtt_nauc = ogtt_pauc - ogtt_iauc


###### OGTT Peak
ogtt_max = apply(ogtt_imputed, 2, function(x){
  max(x)
})


###### FPG (Glucose at 0)
ogtt_fpg = ogtt_imputed["0", ]


###### Glucose at 60
ogtt_60 = ogtt_imputed["60", ]


###### Glucose at 120
ogtt_120 = ogtt_imputed["120", ]


###### Glucose at 180
ogtt_180 = ogtt_imputed["180", ]


###### CV
ogtt_cv = apply(ogtt_imputed, 2, cv)


###### Curve size
ogtt_curve_size = apply(diff(ogtt_imputed), 2, function(x){
  sum(abs(x))
})


###### Time to peak
time = as.numeric(rownames(ogtt_imputed))
ogtt_time_baseline_peak = apply(ogtt_imputed, 2, function(x){
  mx = max(x)
  idx = which(x == mx)
  t = time[idx[1]]
})
ogtt_time_baseline_peak = unlist(ogtt_time_baseline_peak)


###### Time it takes to arrive to first sample below baseline after peak
time = as.numeric(rownames(ogtt_imputed))
ogtt_time_peak_baseline = apply(ogtt_imputed, 2, function(x){
  mx = max(x)
  idx_max = which(x == mx)
  idx_less_baseline = which(x<x[1])
  idx_less_baseline_more_max = which(idx_less_baseline > idx_max[1])
  idx = idx_less_baseline[idx_less_baseline_more_max[1]]
  t_peak_baseline = time[idx[1]] - time[idx_max[1]]
})
ogtt_time_peak_baseline = unlist(ogtt_time_peak_baseline)


###### Slope from baseline to peak
time = as.numeric(rownames(ogtt_imputed))
ogtt_slope_baseline_peak = apply(ogtt_imputed, 2, function(x){
  mx = max(x)
  idx = which(x == mx)
  slope = as.numeric((x[idx[1]]-x[1])/(time[idx[1]]-time[1]))
})
ogtt_slope_baseline_peak = unlist(ogtt_slope_baseline_peak)


###### Slope from peak to last time point
time = as.numeric(rownames(ogtt_imputed))
ogtt_slope_peak_last = apply(ogtt_imputed, 2, function(x){
  mx = max(x)
  idx = which(x == mx)
  slope = as.numeric(x[length(time)] - x[idx[1]])/(time[length(time)] - time[idx[1]])
})
ogtt_slope_peak_last = unlist(ogtt_slope_peak_last)


###### Mark subjects that went below basline
ogtt_time_below_basline = apply(ogtt_imputed, 2, function(x){
  mx = max(x)
  idx_max = which(x == mx)
  idx_less_baseline = which(x<x[1])
  idx_less_baseline_more_max = which(idx_less_baseline > idx_max[1])
  length(idx_less_baseline_more_max)>0
})


#####################################################
######## Visualize features in annotation bar  ######
#####################################################
ogtt_features = data.frame(ogtt_fpg = ogtt_fpg, ogtt_60 = as.numeric(ogtt_60), ogtt_120 = as.numeric(ogtt_120), 
                           ogtt_180 = as.numeric(ogtt_180), ogtt_auc = ogtt_auc, ogtt_iauc = as.numeric(ogtt_iauc),
                           ogtt_pauc = as.numeric(ogtt_pauc), ogtt_nauc = as.numeric(ogtt_nauc), ogtt_max = ogtt_max,  
                           ogtt_curve_size = ogtt_curve_size, ogtt_cv = ogtt_cv,
                           ogtt_time_baseline_peak = ogtt_time_baseline_peak, ogtt_time_peak_baseline = ogtt_time_peak_baseline, 
                           ogtt_slope_baseline_peak = ogtt_slope_baseline_peak, ogtt_slope_peak_last = ogtt_slope_peak_last, 
                           ogtt_time_below_basline = ogtt_time_below_basline)
write.csv(ogtt_features, file = "ogtt_features_09102021.csv")
ogtt_features_demographics_clinical = merge(demographics, ogtt_features, by.x = "SubjectID", by.y = "row.names")
write.csv(ogtt_features_demographics_clinical, file = "ogtt_features_demographics_clinical_09102021.csv")



#############
###  Correlations between PC1 and curve features
##########

### Variance explained
summary(df_pca)[2,]



df_2_pcs = df_pca$x[,c("PC1", "PC2")]

correlation_PCs_OGTT_extracted_features = t(cor(df_2_pcs, ogtt_features))
# cor.test(df_pca$x[,c("PC1")], ogtt_features$ogtt_fpg)

write.csv(correlation_PCs_OGTT_extracted_features, "data/correlation_PCs_OGTT_extracted_features.csv")




ogtt_features_mod = ogtt_features[,-which(colnames(ogtt_features) == "ogtt_time_below_basline")] 
cor.matrix = rcorr(x = as.matrix(df_2_pcs),
                   y = as.matrix(ogtt_features_mod), type = "pearson")
R <- data.frame(cor.matrix$r)
P <- data.frame(cor.matrix$P)

# Data frame with r values
R.df <- R[,-c(1:ncol(df_2_pcs))]
R.df <- R.df[-c((ncol(df_2_pcs))+1:ncol(R)),]
R.df

# Data frame with p values
P.df <- P[,-c(1:ncol(df_2_pcs))]
P.df <- P.df[-c((ncol(df_2_pcs))+1:ncol(P)),]
P.df


##### Visualize Correlation  ###
test_labels = as.matrix(R.df)
test_labels[,] = ""
# seen = vector()
for(i in 1:nrow(P.df)){
  for(j in 1:ncol(P.df)){
    if(P.df[i,j] < 0.05){
      test_labels[i,j] = "x"
    } 
  }
}

## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))
paletteLength = 15
brewer.pal(n = 8, name = "Spectral")
myColor = colorRampPalette(c("#3288BD", "#f7f7f7", "#D53E4F"))(paletteLength)
myColor = colorRampPalette(c("Purple", "#f7f7f7", "Green"))(paletteLength)
myBreaks <- c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(1/paletteLength, 1, length.out=floor(paletteLength/2)))




dev.off()
pdf("figures/correlation_PCs_Features.pdf", w = 11, h = 3)
pheatmap(R.df,
         #color = viridis(100),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         display_numbers = test_labels,
         color = myColor,
         breaks = myBreaks,
         fontsize = 12,
         fontsize_number = 13,
         angle_col = 45,
         # main = "",
         main = "Correlation between OGTT curve features and reduced representation (top 2 PC)"
)
dev.off()






# ####################################################################
# ##### Correlation of curve features with metabolic phenotypes   ####
# ####################################################################
metabolic_indicators_4_corr = metabolic_indicators_orderd[, c("SubjectID", "sspg", "DI", "ie", "hepatic_IR")]
colnames(metabolic_indicators_4_corr)[which(colnames(metabolic_indicators_4_corr) == "sspg")] = "Muscle IR"
colnames(metabolic_indicators_4_corr)[which(colnames(metabolic_indicators_4_corr) == "DI")] = "Beta-cell Function"
colnames(metabolic_indicators_4_corr)[which(colnames(metabolic_indicators_4_corr) == "ie")] = "Incretin Effect"
colnames(metabolic_indicators_4_corr)[which(colnames(metabolic_indicators_4_corr) == "hepatic_IR")] = "Hepatic IR"
rownames(metabolic_indicators_4_corr) = metabolic_indicators_4_corr$SubjectID
metabolic_indicators_4_corr = metabolic_indicators_4_corr[,-1]

idx_rm_features = which(colnames(ogtt_features) %in% c("ogtt_time_peak_baseline", "ogtt_time_below_basline", "phasic_annotation"))
ogtt_features_4_corr = ogtt_features[, -idx_rm_features]
dim(ogtt_features_4_corr)
ogtt_features_4_corr = na.omit(ogtt_features_4_corr)
dim(ogtt_features_4_corr)
ogtt_features_4_corr = ogtt_features_4_corr[rownames(metabolic_indicators_4_corr),]
dim(ogtt_features_4_corr)

cor.matrix = rcorr(x = as.matrix(metabolic_indicators_4_corr),
                   y = as.matrix(ogtt_features_4_corr), type = "pearson")
R <- data.frame(cor.matrix$r)
P <- data.frame(cor.matrix$P)

# Data frame with r values
R.df <- R[,-c(1:ncol(metabolic_indicators_4_corr))]
R.df <- R.df[-c((ncol(metabolic_indicators_4_corr))+1:ncol(R)),]
R.df

# Data frame with p values
P.df <- P[,-c(1:ncol(metabolic_indicators_4_corr))]
P.df <- P.df[-c((ncol(metabolic_indicators_4_corr))+1:ncol(P)),]
P.df


##### Visualize Correlation  ###
test_labels = as.matrix(R.df)
test_labels[,] = ""
# seen = vector()
for(i in 1:nrow(P.df)){
  for(j in 1:ncol(P.df)){
    if(P.df[i,j] < 0.05){
      test_labels[i,j] = "x"
    } 
  }
}

## Edit body of pheatmap:::draw_colnames, customizing it to your liking
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))
paletteLength = 15
brewer.pal(n = 8, name = "Spectral")
myColor = colorRampPalette(c("#3288BD", "#f7f7f7", "#D53E4F"))(paletteLength)
myBreaks <- c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1),
              seq(1/paletteLength, 1, length.out=floor(paletteLength/2)))


### Rename columns
names(R.df)
colnames(R.df)[which(colnames(R.df) == "ogtt_fpg")] = "G0"
colnames(R.df)[which(colnames(R.df) == "ogtt_60")] = "G60"
colnames(R.df)[which(colnames(R.df) == "ogtt_120")] = "G120"
colnames(R.df)[which(colnames(R.df) == "ogtt_180")] = "G180"
colnames(R.df)[which(colnames(R.df) == "ogtt_auc")] = "AUC"
colnames(R.df)[which(colnames(R.df) == "ogtt_iauc")] = "iAUC"
colnames(R.df)[which(colnames(R.df) == "ogtt_pauc")] = "pAUC"
colnames(R.df)[which(colnames(R.df) == "ogtt_nauc")] = "nAUC"
colnames(R.df)[which(colnames(R.df) == "ogtt_max")] = "G_Peak"
colnames(R.df)[which(colnames(R.df) == "ogtt_curve_size")] = "CurveSize"
colnames(R.df)[which(colnames(R.df) == "ogtt_cv")] = "CV"
colnames(R.df)[which(colnames(R.df) == "ogtt_time_baseline_peak")] = "T_baseline2peak"
colnames(R.df)[which(colnames(R.df) == "ogtt_slope_baseline_peak")] = "S_baseline2peak"
colnames(R.df)[which(colnames(R.df) == "ogtt_slope_peak_last")] = "S_peak2end"


dev.off()
pdf("ogtt_features_metabolicpanels_corr_12202021.pdf", w = 8, h = 3)
pheatmap(R.df,
         #color = viridis(100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         display_numbers = test_labels,
         color = myColor,
         breaks = myBreaks,
         fontsize = 12,
         fontsize_number = 15,
         angle_col = 45,
         main = ""
         # main = "Correlation between OGTT curve features and Metabolic indicators"
)
dev.off()



####################################################################################
##### Association between metabolic indicator vs metabolic subphenotype
####################################################################################
aggregated_metabolic_indicators_prs = merge(aggregated_metabolic_indicators, prs, by = "SubjectID")
metabolic_indicators_4_lm = aggregated_metabolic_indicators_prs[,c("a1c", "fpg", "ogtt_2h", "sspg", "DI", "ie", "hepatic_IR")]
colnames(metabolic_indicators_4_lm) = c("HbA1c", "FPG", "OGTT_2h", "Muscle IR", "Beta-cell Function", "Incretin Effect", "Hepatic IR")

mps_melted = reshape2::melt(metabolic_indicators_4_lm,  
                                    id.vars= c("HbA1c", "FPG", "OGTT_2h"))
colnames(mps_melted)[which(colnames(mps_melted) == "variable")] = "msp"
colnames(mps_melted)[which(colnames(mps_melted) == "value")] = "msp_value"
mps_melted_indicators = reshape2::melt(mps_melted,  id.vars= c("msp", "msp_value"))
colnames(mps_melted_indicators)[which(colnames(mps_melted_indicators) == "variable")] = "glucoseindicator"
colnames(mps_melted_indicators)[which(colnames(mps_melted_indicators) == "value")] = "glucoseindicator_value"


gm = ggplot(mps_melted_indicators, aes(x=msp_value, y=glucoseindicator_value))+
  geom_point(aes(color = msp, fill = msp))+
  geom_smooth(method="lm", aes(color = msp, fill = msp)) + 
  facet_grid(glucoseindicator~ msp, scales = "free") + 
  theme_bw() +
  labs(title = "", x="", y="") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.spacing = unit(1, "lines")) + 
  theme(axis.text.x = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=12, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=12, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text = element_text(size=15, face="plain"),
        legend.title = element_blank(),
        legend.position="None",
        strip.text = element_text(size = 13))
gm
ggsave('lm_msp_glucoseindicator_12202021.pdf', h=5, w=8);


####### A1c ~ MSP
lm_mod = lm(a1c~ bmi + age + sex + T2D_PRS + sspg, aggregated_metabolic_indicators_prs)
summary(lm_mod)
summ(lm_mod, scale = TRUE)

lm_mod = lm(a1c~ bmi + age + sex + T2D_PRS + DI, aggregated_metabolic_indicators_prs)
summary(lm_mod)
summ(lm_mod, scale = TRUE)

lm_mod = lm(a1c~ bmi + age + sex + T2D_PRS + ie, aggregated_metabolic_indicators_prs)
summary(lm_mod)
summ(lm_mod, scale = TRUE)

lm_mod = lm(a1c~ bmi + age + sex + T2D_PRS + hepatic_IR, aggregated_metabolic_indicators_prs)
summary(lm_mod)
summ(lm_mod, scale = TRUE)


####### FPG ~ MSP
lm_mod = lm(fpg~ bmi + age + sex + T2D_PRS + sspg, aggregated_metabolic_indicators_prs)
summary(lm_mod)
summ(lm_mod, scale = TRUE)

lm_mod = lm(fpg~ bmi + age + sex + T2D_PRS + DI, aggregated_metabolic_indicators_prs)
summary(lm_mod)
summ(lm_mod, scale = TRUE)

lm_mod = lm(fpg~ bmi + age + sex + T2D_PRS + ie, aggregated_metabolic_indicators_prs)
summary(lm_mod)
summ(lm_mod, scale = TRUE)

lm_mod = lm(fpg~ bmi + age + sex + T2D_PRS + hepatic_IR, aggregated_metabolic_indicators_prs)
summary(lm_mod)
summ(lm_mod, scale = TRUE)


####### ogtt_2h ~ MSP
lm_mod = lm(ogtt_2h~ bmi + age + sex + T2D_PRS + sspg, aggregated_metabolic_indicators_prs)
summary(lm_mod)
summ(lm_mod, scale = TRUE)

lm_mod = lm(ogtt_2h~ bmi + age + sex + T2D_PRS + DI, aggregated_metabolic_indicators_prs)
summary(lm_mod)
summ(lm_mod, scale = TRUE)

lm_mod = lm(ogtt_2h~ bmi + age + sex + T2D_PRS + ie, aggregated_metabolic_indicators_prs)
summary(lm_mod)
summ(lm_mod, scale = TRUE)

lm_mod = lm(ogtt_2h~ bmi + age + sex + T2D_PRS + hepatic_IR, aggregated_metabolic_indicators_prs)
summary(lm_mod)
summ(lm_mod, scale = TRUE)


############################
###### MSP prediction
###############################
ml_performance = read.csv("2021-09-16_22-21-21_ML_Performance.csv")
ml_performance[which(ml_performance$Subphenotype == "muscle_ir_status"), "Subphenotype"] = "Muscle IR"
ml_performance[which(ml_performance$Subphenotype == "di_status"), "Subphenotype"] = "Beta-cell Dysfunction"
ml_performance[which(ml_performance$Subphenotype == "ie_status"), "Subphenotype"] = "Impaired Incretin Effect"
ml_performance[which(ml_performance$Subphenotype == "hepatic_ir_status"), "Subphenotype"] = "Hepatic IR"


ml_performance[which(ml_performance$Features == "demographics_lab"), "Features"] = "Demographics+Lab"
ml_performance[which(ml_performance$Features == "demographics"), "Features"] = "Demographics"
ml_performance[which(ml_performance$Features == "lab"), "Features"] = "Lab"
ml_performance[which(ml_performance$Features == "homair"), "Features"] = "HOMA_IR"
ml_performance[which(ml_performance$Features == "homab"), "Features"] = "HOMA_B"
ml_performance[which(ml_performance$Features == "matsuda"), "Features"] = "Matsuda"
ml_performance[which(ml_performance$Features == "incretins"), "Features"] = "Incretins"
ml_performance[which(ml_performance$Features == "ogtt_timeseries"), "Features"] = "OGTT_G_ReducedRep"
ml_performance[which(ml_performance$Features == "ogtt_features"), "Features"] = "OGTT_G_Features"
ml_performance[which(ml_performance$Features == "combine_ogtt_feature_timeseries"), "Features"] = "OGTT_G_all"


table(ml_performance$Features)

###### Draw auROC for best performers
cdata <- ddply(ml_performance, c("Features", "Subphenotype", "Classifier"), summarise,
               N    = length(auROC),
               mean = mean(auROC),
               sd   = sd(auROC),
               se   = sd / sqrt(N)
)
cdata = subset(cdata, Classifier != "logistic")
cdata = subset(cdata, Features != "OGTT_G_all")
cdata = subset(cdata, Features != "Lab")

cdata_seletced = cdata %>% group_by(Features, Subphenotype) %>% slice(which.max(mean))
cdata_seletced$Features = factor(cdata_seletced$Features, levels = c("PRS", "Demographics", "Lab", "Demographics+Lab",
                                                                     "HOMA_B", "HOMA_IR", "Matsuda", "Incretins",
                                                                     "OGTT_G_Features",  "OGTT_G_ReducedRep")) 

cdata_seletced = subset(cdata_seletced, Subphenotype %in% c("Muscle IR", "Beta-cell Dysfunction", 
                                                        "Impaired Incretin Effect", "Hepatic IR"))
cdata_seletced$Subphenotype = factor(cdata_seletced$Subphenotype, levels = c("Muscle IR", "Beta-cell Dysfunction", 
                                                                             "Impaired Incretin Effect", "Hepatic IR"))

cbPalette <- c("#756bb1", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#c51b8a")
gp = ggplot(cdata_seletced, aes_string(x="Subphenotype", y="mean", fill = "Features")) + 
  geom_bar(stat = "identity", position="dodge", color = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position=position_dodge(.9), width=0.2) + 
  scale_fill_manual(values=cbPalette) + 
  theme_bw() +
  coord_cartesian(ylim = c(0.4, 1.23)) + 
  theme(legend.position="top") + 
  theme(axis.text.x = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="bold"),
        axis.text.y = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"), #legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 14, color = "black", face = "bold")) +
  labs(y = "auROC", x = "") +
  ggtitle("Prediction of Metabolic Subphenotypes")
gp
ggsave(sprintf('ML_performance_auc_best_all_demographics_01042022_1.pdf'), h=10, w=11)


################################
#### Draw best classifier for each metabolic subphenotype
################################
names(cdata_seletced)
cdata_seletced_4_best = cdata_seletced[, c("Features", "Subphenotype", "Classifier")]

gp = ggplot(cdata_seletced_4_best, aes(x = Subphenotype, y = Features)) + 
  geom_raster(aes(fill=Classifier, color = "grey50")) + 
  scale_fill_manual(values=c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")) + 
  # geom_tile(aes(fill = Classifier, width=0.7, height=0.7), size=2) + 
  theme_bw() +
  theme(legend.position="top") + 
  theme(axis.text.x = element_text(colour="black", size=14, angle=45, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"), #legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  theme(strip.text.x = element_text(size = 14, color = "black", face = "bold")) +
  labs(y = "Feature set", x = "") +
  ggtitle("Best Performing Classifier")
gp
ggsave(sprintf('ML_performance_auc_best_performing_clasifier_01042022_1.pdf'), h=8, w=10)
ggsave(sprintf('ML_performance_auc_best_performing_clasifier_01042022_1.jpeg'), h=8, w=10)




### Wilcox Rank Test
## MuscleIR
x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Muscle IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "Demographics" & Subphenotype == "Muscle IR" & Classifier == "Linear_SVC")$auROC
wilcox.test(x, y) 

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Muscle IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "Demographics+Lab" & Subphenotype == "Muscle IR" & Classifier == "RBF_SVC")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Muscle IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "HOMA_B" & Subphenotype == "Muscle IR" & Classifier == "Linear_SVC")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Muscle IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "HOMA_IR" & Subphenotype == "Muscle IR" & Classifier == "L1_logistic")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Muscle IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "Matsuda" & Subphenotype == "Muscle IR" & Classifier == "L1_logistic")$auROC
wilcox.test(x, y) 

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Muscle IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "Incretins" & Subphenotype == "Muscle IR" & Classifier == "GPC")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Muscle IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "OGTT_G_Features" & Subphenotype == "Muscle IR" & Classifier == "RBF_SVC")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Muscle IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "PRS" & Subphenotype == "Muscle IR" & Classifier == "L2_Logistic")$auROC
wilcox.test(x, y)



## BetaCell
x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "L2_Logistic")$auROC
y = subset(ml_performance, Features == "Demographics" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "GPC")$auROC
wilcox.test(x, y) 

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "L2_Logistic")$auROC
y = subset(ml_performance, Features == "Demographics+Lab" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "RBF_SVC")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "L2_Logistic")$auROC
y = subset(ml_performance, Features == "HOMA_B" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "RBF_SVC")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "L2_Logistic")$auROC
y = subset(ml_performance, Features == "HOMA_IR" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "Linear_SVC")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "L2_Logistic")$auROC
y = subset(ml_performance, Features == "Matsuda" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "Linear_SVC")$auROC
wilcox.test(x, y) 

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "L2_Logistic")$auROC
y = subset(ml_performance, Features == "Incretins" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "GPC")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "L2_Logistic")$auROC
y = subset(ml_performance, Features == "OGTT_G_Features" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "L1_logistic")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "L2_Logistic")$auROC
y = subset(ml_performance, Features == "PRS" & Subphenotype == "Beta-cell Dysfunction" & Classifier == "GPC")$auROC
wilcox.test(x, y)



## IE
x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Impaired Incretin Effect" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "Demographics" & Subphenotype == "Impaired Incretin Effect" & Classifier == "GPC")$auROC
wilcox.test(x, y) 

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Impaired Incretin Effect" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "Demographics+Lab" & Subphenotype == "Impaired Incretin Effect" & Classifier == "GPC")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Impaired Incretin Effect" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "HOMA_B" & Subphenotype == "Impaired Incretin Effect" & Classifier == "L1_logistic")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Impaired Incretin Effect" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "HOMA_IR" & Subphenotype == "Impaired Incretin Effect" & Classifier == "L1_logistic")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Impaired Incretin Effect" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "Matsuda" & Subphenotype == "Impaired Incretin Effect" & Classifier == "GPC")$auROC
wilcox.test(x, y) 

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Impaired Incretin Effect" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "Incretins" & Subphenotype == "Impaired Incretin Effect" & Classifier == "GPC")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Impaired Incretin Effect" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "OGTT_G_Features" & Subphenotype == "Impaired Incretin Effect" & Classifier == "L1_logistic")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Impaired Incretin Effect" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "PRS" & Subphenotype == "Impaired Incretin Effect" & Classifier == "RBF_SVC")$auROC
wilcox.test(x, y)



## HepaticIR
x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Hepatic IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "Demographics" & Subphenotype == "Hepatic IR" & Classifier == "Linear_SVC")$auROC
wilcox.test(x, y) 

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Hepatic IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "Demographics+Lab" & Subphenotype == "Hepatic IR" & Classifier == "GPC")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Hepatic IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "HOMA_B" & Subphenotype == "Hepatic IR" & Classifier == "Linear_SVC")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Hepatic IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "HOMA_IR" & Subphenotype == "Hepatic IR" & Classifier == "L2_Logistic")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Hepatic IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "Matsuda" & Subphenotype == "Hepatic IR" & Classifier == "L1_logistic")$auROC
wilcox.test(x, y) 

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Hepatic IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "Incretins" & Subphenotype == "Hepatic IR" & Classifier == "GPC")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Hepatic IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "OGTT_G_Features" & Subphenotype == "Impaired Incretin Effect" & Classifier == "Linear_SVC")$auROC
wilcox.test(x, y)

x = subset(ml_performance, Features == "OGTT_G_ReducedRep" & Subphenotype == "Hepatic IR" & Classifier == "GPC")$auROC
y = subset(ml_performance, Features == "PRS" & Subphenotype == "Impaired Incretin Effect" & Classifier == "GPC")$auROC
wilcox.test(x, y)



################
### T2D PRS
################
genotyping_df = read.csv("data/genotyping_OR_12202021.csv")
genotyping_df$subject_id = gsub("43883-0", "S", genotyping_df$subject_id)
genotyping_df = subset(genotyping_df, subject_id !="S18")

## Stat
genotyping_df$OR
genotyping_df$a1c_avg
cor(genotyping_df$OR, genotyping_df$a1c_avg)
cor.test(genotyping_df$OR, genotyping_df$a1c_avg, method="pearson")

### Read Sean's Genotyping
genotyping_df_sean = read.csv("data/CGM_genotyping_updated_08222023.csv")
colnames(genotyping_df_sean)[1] = "subject_id"
genotyping_df_sean = genotyping_df_sean[,c("subject_id", "T2D.GRS_norm")]
colnames(genotyping_df_sean)
dim(genotyping_df_sean)

# Stats
cor(indictaos_metabolic_prs_4_prediction$T2D.GRS_norm, indictaos_metabolic_prs_4_prediction$a1c)
cor.test(indictaos_metabolic_prs_4_prediction$T2D.GRS_norm, indictaos_metabolic_prs_4_prediction$a1c, method="pearson")
colnames(indictaos_metabolic_prs_4_prediction)

###################################################
##### Prepare aggregate file for classification ###
######################################################
# aggregated_metabolic_indicators_4_prediction = read.csv("aggregated_metabolic_indicators_09102021.csv")
names(genotyping_df)
genotyping_df_sub = genotyping_df[,c("subject_id", "OR")]
# indictaos_metabolic_prs_4_prediction = merge(aggregated_metabolic_indicators_4_prediction, genotyping_df_sub, by.x = "SubjectID", by.y="subject_id")
indictaos_metabolic_prs_4_prediction = merge(metabolic_indicators_orderd, genotyping_df_sub, by.x = "SubjectID", by.y="subject_id")
write.csv(indictaos_metabolic_prs_4_prediction, "aggregated_metabolic_indicators__4_predictions_w_prs_01042022.csv")


## Do it for Sean's genotyping
indictaos_metabolic_prs_4_prediction = merge(metabolic_indicators_orderd, genotyping_df_sean, by.x = "SubjectID", by.y="subject_id")


### Visualize T2D ordered by A1C
# group.color = c(Normal = "#00AFBB", PreDM = "#E7B800", T2D ="#FC4E07")
group.color = c(Normoglycemic = "#00AFBB", PreDM = "#E7B800", T2D ="#FC4E07")
indictaos_metabolic_prs_4_prediction$a1c_t2d_status <- factor(indictaos_metabolic_prs_4_prediction$a1c_t2d_status, 
                                       levels = c("Normoglycemic", "PreDM", "T2D"))

# gp = ggplot(indictaos_metabolic_prs_4_prediction, aes(x =  reorder(SubjectID, -OR), y = OR)) + 
gp = ggplot(indictaos_metabolic_prs_4_prediction, aes(x =  reorder(SubjectID, -T2D.GRS_norm), y = T2D.GRS_norm)) + 
  geom_bar(aes(fill = a1c_t2d_status), stat = "identity") + 
  scale_fill_manual(values=group.color) +
  theme_bw() +
  annotate("text", x = 27, y = 0.8, label = "r(PRS, HbA1c)== 0.48",
           parse = TRUE) +
  annotate("text", x = 27, y = 0.7, label = "p == 0.005",
           parse = TRUE) +
  theme(legend.position="top") + 
  theme(axis.text.x = element_text(colour="black", size=14, angle=90, hjust=0.5, vjust=0.5, face="plain"),
        axis.text.y = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"),# legend.title =  element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(fill=guide_legend(title="Glycemic status")) +
  theme(strip.text.x = element_text(size = 14, color = "black", face = "bold")) +
  labs(y = "T2D Polygenic Risk Score", x = "") +
  ggtitle("")
gp
ggsave(sprintf('genotyping_a1c_10272023_corrected_.pdf'), h=4, w=7)


### Correlations pvalues
pgrs_p = read.csv("genotyping/pGSR_01042021/corr_pval.csv", row.names = 1)
pgrs_r = read.csv("genotyping/pGSR_01042021/Correlation.csv", row.names = 1)

##### Visualize Correlation  ###
R.df <- data.frame(pgrs_r)
P.df <- data.frame(pgrs_p)

### Rename columns
# names(R.df)
colnames(R.df)[which(colnames(R.df) == "PRS")] = "PRS_T2D"
colnames(R.df)[which(colnames(R.df) == "Combination.of.Insulin.Secretion.2.and.Insulin.Action")] = "pGRS_Insulin_Secretion+Action"
colnames(R.df)[which(colnames(R.df) == "Insulin.Action")] = "pGRS_InsulinAction"
colnames(R.df)[which(colnames(R.df) == "Insulin.Secretion.2")] = "pGRS_InsulinSecretion_2"
colnames(R.df)[which(colnames(R.df) == "Adiposity")] = "pGRS_Adiposity"
colnames(R.df)[which(colnames(R.df) == "Impaired.Lipid.Metabolism")] = "pGRS_LipidMetabolism"
colnames(R.df)[which(colnames(R.df) == "Insulin.Secretion.1")] = "pGRS_InsulinSecretion_1"


rownames(R.df)[which(rownames(R.df) == "Total cholesterol")] = "Cholesterol"
rownames(R.df)[which(rownames(R.df) == "IE 180min")] = "IE %"
rownames(R.df)[which(rownames(R.df) == "Modified DI")] = "DI"
rownames(R.df)[which(rownames(R.df) == "HOMA B")] = "HOMA_B"
rownames(R.df)[which(rownames(R.df) == "HOMA S")] = "HOMA_S"
rownames(R.df)[which(rownames(R.df) == "HOMA IR")] = "HOMA_IR"


R.df_2 = R.df
R.df_2$indicator = rownames(R.df_2)
pgrs_r_melted = melt(R.df_2, variables="indicator")
colnames(pgrs_r_melted)[3] = "R"


P.df_2 = P.df
P.df_2$indicator = rownames(P.df_2)
pgrs_p_melted = melt(P.df_2, variables="indicator")
colnames(pgrs_p_melted)[3] = "P"


df_corr_plot = pgrs_r_melted
df_corr_plot$P = pgrs_p_melted$P
df_corr_plot$P_log = -log10(df_corr_plot$P)

### PRS vs indicators
df_corr_plot_subset = subset(df_corr_plot, variable == "PRS_T2D")
gp = ggplot(df_corr_plot_subset, aes(x=variable, y=indicator, color=R)) +
# gp = ggplot(df_corr_plot_subset, aes(x=indicator, y=variable, color=R)) +
  geom_point(aes(size=P_log)) +
  scale_color_distiller(palette="RdGy") +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="right"#,
   # legend.title = element_blank()
  ) +
  theme(axis.text.x = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=1, face="plain"),
        axis.text.y = element_text(colour="black", size=14, angle=0, hjust=1, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"), #legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  guides(size=guide_legend(title="p-value(-log10)")) +
  ggtitle("") +
  xlab("") +
  ylab("")
gp
ggsave(sprintf('prs_indicators_01032021_1.pdf'), h=7, w=4)



################
### PRS vs A1C
##################

aggregated_metabolic_indicators
ggplot(aggregated_metabolic_indicators)
gp = ggplot(aggregated_metabolic_indicators, aes(x=a1c, y=OR)) +
  geom_point(aes(color=a1c_t2d_status)) +
  geom_smooth(method='lm') +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position="right"#,
    # legend.title = element_blank()
  ) +
  theme(axis.text.x = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=1, face="plain"),
        axis.text.y = element_text(colour="black", size=14, angle=0, hjust=1, vjust=0.5, face="plain"),
        axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),
        axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),
        legend.text=element_text(size=14, face="plain"), #legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("") +
  xlab("HbA1c") +
  ylab("PRS")
gp

ggsave("/Users/ahmedm/Library/CloudStorage/Box-Box/Ahmed Metwally's Files/Stanford/cgm/43883/Manuscripts/Manuscript_1_MetabolicSubphenotyping/PRS_A1C.pdf", h=5, w=5)

model_prs_a1c = lm(a1c~OR, aggregated_metabolic_indicators)
summary(model_prs_a1c)


