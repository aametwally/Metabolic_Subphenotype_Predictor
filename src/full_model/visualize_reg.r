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
pardir="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_dataset/ogtt_cgm/"
resdir=paste0(pardir,"res/")
setwd(resdir)
# 
restab=read.table(paste0("2023-09-27_10-18-04_iterations_100_ML_regression_Performance.csv"),sep=",",header=TRUE)
# name replace
replist=list("extrfeatures_NAN"="Manual defined")
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
    cdata<-ddply(restab_auc, c("Features","Extraction","Regressor","hyperparameter"),summarise,
               N=length(get(perfterm)),
               meanval=mean(get(perfterm)),
               sdval=sd(get(perfterm)),
               seval=sdval/sqrt(N)
    )
    summ_list[[perfterm]]=cdata
}
cdata_seletced_valid<-summ_list[["RRMSE_validation"]]%>%group_by(Extraction,Regressor,hyperparameter)%>%summarize(avgmean=mean(meanval),.group="drop")%>%as.data.frame()%>%slice(which.min(avgmean))
matchcol=c("Extraction","Regressor","hyperparameter")
grcolors=c(venous_exp_with_matching="red",CGM_ALL="blue",ctru_cgm="green")
for(perfterm in perflist){
    cdata_seletced=merge(cdata_seletced_valid[,matchcol],summ_list[[perfterm]],by=matchcol,all.x=TRUE)
    cdata_seletced$Features=factor(cdata_seletced$Features,levels=sort(unique(cdata_seletced$Features)))
    cdata_seletced$Extraction=factor(cdata_seletced$Extraction,levels=unique(cdata_seletced$Extraction))
    p<-ggplot(cdata_seletced[order(cdata_seletced$Features),],aes(x=Extraction,y=meanval,fill=Features))+geom_bar(stat="identity",position="dodge",color="black")+geom_errorbar(aes(ymin=meanval-2*seval,ymax=meanval+2*seval),position=position_dodge(.9),width=0.2)+theme_bw()+coord_cartesian(ylim=c(0.2,1.23))+scale_fill_manual(values=grcolors)+theme(legend.position="top")+theme(axis.text.x=element_text(colour="black",size=14,angle=0,hjust=0.5,vjust=0.5,face="bold"),axis.text.y=element_text(colour="black",size=14,angle=0,hjust=0.5,vjust=0.5,face="plain"),axis.title.x=element_text(colour="black",size=14,angle=0,hjust=.5,vjust=0.5,face="bold"),axis.title.y=element_text(colour="black",size=14,angle=90,hjust=.5,vjust=.5,face="bold"),legend.text=element_text(size=14,face="plain"),plot.title=element_text(hjust=0.5,face="bold"))+theme(strip.text.x=element_text(size=14,color="black",face="bold"))+labs(y=perfterm,x="")+ggtitle("Prediction of SSPG")
    ggsave(paste0("ML_performance_auc_best_all_regress",perfterm,".pdf"),h=10,w=11,plot=p)
}
