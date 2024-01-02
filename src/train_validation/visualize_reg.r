rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(tidyr)
require(plyr)
require(dplyr)
require(reshape2)
require(ComplexHeatmap)
# 
pardir="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_dataset/"
datadir=paste0(pardir,"data_Metabolic_subphenotype_and_CGM_analysis/")
resdir=paste0(pardir,"res/cgm_ogtt/")
setwd(resdir)
# 
restab=read.table(paste0("2023-01-27_10-38-25_iterations_20_ML_regression_Performance.csv"),sep=",",header=TRUE)
# name replace
replist=list("timeseries_regular"="PCA","timeseries_fpca"="FPCA","extrfeatures_NAN"="Manual defined")
vec=restab[,"Extraction"]
for(prename in names(replist)){
    vec[vec==prename]=replist[[prename]]
}
restab[,"Extraction"]=vec
# vis performance
perflist=c("RRMSE_validation","RRMSE","Correlation","Accuracy","F1","Precision","Recall")
summ_list=list()
for(perfterm in perflist){
    restab_auc=restab[!is.na(restab[,perfterm]),]
    cdata<-ddply(restab_auc, c("Features","Extraction","Regresser"),summarise,
               N=length(get(perfterm)),
               mean=mean(get(perfterm)),
               sd=sd(get(perfterm)),
               se=sd/sqrt(N)
    )
    summ_list[[perfterm]]=cdata
}
cdata_seletced_valid<-summ_list[["RRMSE_validation"]] %>% group_by(Features,Extraction) %>% slice(which.min(mean))
matchcol=c("Features","Extraction","Regresser")
for(perfterm in perflist){
    cdata_seletced=merge(cdata_seletced_valid[,matchcol],summ_list[[perfterm]],by=matchcol,all.x=TRUE)
    # merge
    cdata_seletced$Features=factor(cdata_seletced$Features,levels=unique(cdata_seletced$Features)) 
    cdata_seletced$Extraction=factor(cdata_seletced$Extraction,levels=unique(cdata_seletced$Extraction))
    p<-ggplot(cdata_seletced,aes(x=Extraction,y=mean,fill=Features))+
    geom_bar(stat="identity",position="dodge",color="black")+
    geom_errorbar(aes(ymin=mean-2*se,ymax=mean+2*se),position=position_dodge(.9),width=0.2)+ 
    theme_bw()+
    # coord_cartesian(ylim=c(0.2,1.23))+ 
    theme(legend.position="top")+ 
    theme(axis.text.x = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="bold"),axis.text.y = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),legend.text=element_text(size=14, face="plain"),plot.title = element_text(hjust = 0.5, face = "bold"))+
    theme(strip.text.x=element_text(size = 14, color = "black", face = "bold")) +
    labs(y=perfterm,x="") +
    ggtitle("Prediction of SSPG")
    ggsave(paste0("ML_performance_auc_best_all_regress",perfterm,".pdf"),h=10,w=11,plot=p)
}
# regression RRMSE targetted vis
perfterm="RRMSE"
cdata_seletced=merge(cdata_seletced_valid[,matchcol],summ_list[[perfterm]],by=matchcol,all.x=TRUE)
cdata_seletced=cdata_seletced[cdata_seletced[,"Features"] %in% c("CTRU_Venous","HOME_CGM_MEAN") & cdata_seletced[,"Extraction"] %in% c("Manual defined","PCA"),]
# merge
cdata_seletced$Features=factor(cdata_seletced$Features,levels=unique(cdata_seletced$Features)) 
cdata_seletced$Extraction=factor(cdata_seletced$Extraction,levels=unique(cdata_seletced$Extraction))
p<-ggplot(cdata_seletced,aes(x=Extraction,y=mean,fill=Features))+
geom_bar(stat="identity",position="dodge",color="black")+
geom_errorbar(aes(ymin=mean-2*se,ymax=mean+2*se),position=position_dodge(.9),width=0.2)+ 
theme_bw()+
theme(legend.position="top")+ 
theme(axis.text.x = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="bold"),axis.text.y = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),legend.text=element_text(size=14, face="plain"),plot.title = element_text(hjust = 0.5, face = "bold"))+
theme(strip.text.x=element_text(size = 14, color = "black", face = "bold")) +
labs(y="Relative root mean squared error",x="") +
ggtitle("Prediction of SSPG")
ggsave(paste0("ML_performance_auc_best_all_regress.rrmse.targetvis.pdf"),h=10,w=11,plot=p)
# best methods heatmap
cdata_seletced_valid_df=as.data.frame(cdata_seletced_valid)
featvec=unique(cdata_seletced_valid_df[,"Features"])
extrvec=unique(cdata_seletced_valid_df[,"Extraction"])
platmat=matrix(NA,nrow=length(featvec),ncol=length(extrvec))
rownames(platmat)=featvec
colnames(platmat)=extrvec
for(rowi in seq(dim(cdata_seletced_valid_df)[1])){
    record=cdata_seletced_valid_df[rowi,]
    platmat[record[,"Features"],record[,"Extraction"]]=record[,"Regresser"]
}
# 
methlist=unique(platmat)
colors=structure(seq(length(methlist)),names=methlist)
pdf("method.best.regressor.pdf")
Heatmap(platmat,col=colors)
dev.off()