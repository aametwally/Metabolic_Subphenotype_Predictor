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
restab=read.table(paste0("2023-01-29_09-51-37_iterations_20_ML_Performance.csv"),sep=",",header=TRUE)
# name replace
replist=list("timeseries_regular"="PCA","timeseries_fpca"="FPCA","extrfeatures_NAN"="Manual defined")
vec=restab[,"Extraction"]
for(prename in names(replist)){
    vec[vec==prename]=replist[[prename]]
}
restab[,"Extraction"]=vec
# vis performance
perflist=c("auROC","Accuracy","F1","Precision","Recall")
for(perfterm in perflist){
    restab_auc=restab[!is.na(restab[,perfterm]),]
    cdata<-ddply(restab_auc, c("Features","Extraction","Classifier"),summarise,
                N=length(get(perfterm)),
                mean=mean(get(perfterm)),
                sd=sd(get(perfterm)),
                se=sd/sqrt(N)
    )
    cdata_seletced<-cdata %>% group_by(Features,Extraction) %>% slice(which.max(mean))
    cdata_seletced$Features=factor(cdata_seletced$Features,levels=unique(cdata_seletced$Features)) 
    cdata_seletced$Extraction=factor(cdata_seletced$Extraction,levels=unique(cdata_seletced$Extraction))
    p<-ggplot(cdata_seletced,aes(x=Extraction,y=mean,fill=Features))+
    geom_bar(stat="identity",position="dodge",color="black")+
    geom_errorbar(aes(ymin=mean-2*se,ymax=mean+2*se),position=position_dodge(.9),width=0.2)+ 
    theme_bw()+
    coord_cartesian(ylim=c(0.2,1.23))+ 
    theme(legend.position="top")+ 
    theme(axis.text.x = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="bold"),axis.text.y = element_text(colour="black", size=14, angle=0, hjust=0.5, vjust=0.5, face="plain"),axis.title.x = element_text(colour="black", size=14, angle=0, hjust=.5, vjust=0.5, face="bold"),axis.title.y = element_text(colour="black", size=14, angle=90, hjust=.5, vjust=.5, face="bold"),legend.text=element_text(size=14, face="plain"),plot.title = element_text(hjust = 0.5, face = "bold"))+
    theme(strip.text.x=element_text(size = 14, color = "black", face = "bold")) +
    labs(y=perfterm,x="") +
    ggtitle("Prediction of Metabolic Subphenotypes")
    ggsave(paste0("ML_performance_auc_best_all_validation_class",perfterm,".pdf"),h=10,w=11,plot=p)
}
# best methods heatmap
perfterm="auROC"
restab_auc=restab[!is.na(restab[,perfterm]),]
cdata<-ddply(restab_auc, c("Features","Extraction","Classifier"),summarise,
            N=length(get(perfterm)),
            mean=mean(get(perfterm)),
            sd=sd(get(perfterm)),
            se=sd/sqrt(N)
)
cdata_seletced_valid<- cdata %>% group_by(Features,Extraction) %>% slice(which.max(mean))
cdata_seletced_valid_df=as.data.frame(cdata_seletced_valid)
featvec=unique(cdata_seletced_valid_df[,"Features"])
extrvec=unique(cdata_seletced_valid_df[,"Extraction"])
platmat=matrix(NA,nrow=length(featvec),ncol=length(extrvec))
rownames(platmat)=featvec
colnames(platmat)=extrvec
for(rowi in seq(dim(cdata_seletced_valid_df)[1])){
    record=cdata_seletced_valid_df[rowi,]
    platmat[record[,"Features"],record[,"Extraction"]]=record[,"Classifier"]
}
# 
methlist=unique(platmat)
colors=structure(seq(length(methlist)),names=methlist)
pdf("method.best.classifier.validation.pdf")
Heatmap(platmat,col=colors)
dev.off()