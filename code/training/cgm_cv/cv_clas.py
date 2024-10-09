import os
import numpy as np
import pandas as pd
import datetime
import random
from math import log10
# 
import sklearn as sk
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
import skfda
from skfda.representation.basis import BSplineBasis
from skfda.preprocessing.smoothing import BasisSmoother
from skfda.representation.grid import FDataGrid
from skfda.preprocessing.dim_reduction import FPCA
from skfda.preprocessing.smoothing.validation import SmoothingParameterSearch
from skfda.misc.regularization import L2Regularization
from skfda.misc.operators import LinearDifferentialOperator
# 
import xgboost as xgb
from sklearn.preprocessing import LabelEncoder
# 
import matplotlib.pyplot as plt
## Set project Path and get datetimestamp
pardir="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_dataset/ogtt_cgm/"
datadir=pardir+"data/"
resdir=pardir+"res/"
os.chdir(resdir)
extr_list=["extrfeatures","timeseries"]
setnames=["train","validation"]
features_list=["ctru_cgm","home_cgm_1_matching","home_cgm_2_matching"]
model_list=["CGM_ALL","HOME_CGM_MEAN"]
# Read metabolic indicators 
metabolic_indicators=pd.read_csv(datadir+'metabolic_subphenotyping_all_cohort_metadata_09102023.csv',na_values=['--'])
setsep={"train": np.unique(np.array(metabolic_indicators['subject_id'][metabolic_indicators['exp_type']=="venous_without_matching_cgm_and_without_planned_athome_cgm"])),
"validation": np.unique(np.array(metabolic_indicators['subject_id'][metabolic_indicators['exp_type']=="venous_with_matching_cgm_and_with_planned_athome_cgm"]))}
# setsep['validation']=np.setdiff1d(setsep['validation'],setsep['train'])
setsep['validation']=np.setdiff1d(setsep['validation'],np.array(['S112','S115','S122','S105','S108']))
# 
metabolic_indicators=metabolic_indicators[metabolic_indicators['exp_type']=='venous_with_matching_cgm_and_with_planned_athome_cgm']
# 
demogrlist=["age","bmi","sex","ethnicity"]
metabolic_indicators=metabolic_indicators[["subject_id","sspg_2_classes"]+demogrlist]
for cate in ["sex","ethnicity"]:
    metabolic_indicators=pd.get_dummies(metabolic_indicators,prefix=[cate],columns=[cate], drop_first=True)
#### load timeseries, and features
# cv data
combtab_list_cv=dict()
combtab_list_cv["timeseries"]=dict()
combtab_list_cv["extrfeatures"]=dict()
# load data and formulate table
for feature in features_list:
    features=pd.read_csv(datadir+"metabolic_subphenotyping_ogtt_glucose_"+feature+"_features_09102023.csv", index_col=0)
    features['ogtt_time_peak_baseline']=features['ogtt_time_peak_baseline'].fillna(180)
    features["ogtt_time_below_basline"]=features["ogtt_time_below_basline"].astype(int)
    ### Combine extracted features info with metadata
    combtab=features.merge(metabolic_indicators,left_on="subject_id",right_on="subject_id")
    combtab=combtab[combtab['subject_id'].isin(setsep['validation'])]
    combtab_list_cv["extrfeatures"][feature]=combtab
# combined dataset
combtab=combtab_list_cv["extrfeatures"]["home_cgm_1_matching"]
combtab=pd.concat([combtab,combtab_list_cv["extrfeatures"]["home_cgm_2_matching"]])
combtab=pd.concat([combtab,combtab_list_cv["extrfeatures"]["ctru_cgm"]])
combtab_list_cv["extrfeatures"]["CGM_ALL"]=combtab
# 
locdf_conc=pd.concat([combtab_list_cv["extrfeatures"]["home_cgm_1_matching"],combtab_list_cv["extrfeatures"]["home_cgm_2_matching"]])
locdf_mean=locdf_conc.drop(['sspg_2_classes'],axis=1).groupby("subject_id",as_index=False).mean()
locdf_mean=locdf_mean.merge(locdf_conc[['subject_id','sspg_2_classes']],left_on="subject_id",right_on="subject_id")
combtab_list_cv["extrfeatures"]["HOME_CGM_MEAN"]=locdf_mean
## Combine timeseries with metadata
timeseries=pd.read_csv(datadir+"metabolic_subphenotyping_all_cohort_ogtt_glucose_09102023.csv",index_col=False)
timeseries_sub=timeseries[timeseries["sample_location_extraction_method"].str.contains("CGM")]
timeseries_sub=timeseries_sub[['subject_id','sample_location_extraction_method','timepoint','glucose']]
timeseries_sub_wid=timeseries_sub.pivot(index=["subject_id","sample_location_extraction_method"],columns="timepoint",values="glucose")
timeseries_sub_wid=timeseries_sub_wid.fillna(method='ffill',axis=1)
combtab=timeseries_sub_wid.merge(metabolic_indicators,left_on="subject_id",right_on='subject_id')
combtab.index=combtab["subject_id"]
combtab=combtab[combtab['subject_id'].isin(setsep['validation'])]
combtab_list_cv["timeseries"]["CGM_ALL"]=combtab
timeseq=np.array(timeseries_sub_wid.columns)
timeind=list(range(len(timeseq)))
# 
random.seed(10)
num_iterations=100
num_folds=5
C=10
classifiers={
    'L1_logistic': LogisticRegression(C=C,penalty='l1',solver='saga',max_iter=10000,random_state=0),
    'L2_Logistic': LogisticRegression(C=C,penalty='l2',solver='saga',max_iter=10000,random_state=0),
    # 'Linear_SVC': SVC(kernel='linear',C=C,probability=True,random_state=0),
    'RBF_SVC': SVC(kernel='rbf',C=C,probability=True,random_state=0),
    'RFR': RandomForestClassifier(random_state=23234),
    'XGB': xgb.XGBClassifier(learning_rate=1.0)
}
parameters_grid={
    'C': [0.01, 0.1, 1, 10, 100],
    'n_estimators': [20, 50, 100, 150, 200],
    'learning_rate': [0.001, 0.01, 0.1, 0.5, 1.0]
}
Cset={'L1_logistic': [0.01, 0.1, 1, 10, 100],
      'L2_Logistic': [0.01, 0.1, 1, 10, 100],
    #   'Linear_SVC': [0.001,0.005,0.01,0.05,0.1],
      'RBF_SVC': [0.01, 0.1, 1, 10, 100]}
paraind=[0,0,0,1,2]
tune_para=list(parameters_grid.keys())
extr_list_redim=["timeseries","timeseries","extrfeatures"]
pcmethods=["regular","fpca","NAN"]
loglambda_rag=[-15,5]
lambdavec=np.power(10,np.linspace(start=loglambda_rag[0],stop=loglambda_rag[1],num=int((loglambda_rag[1]-loglambda_rag[0])/0.25+1.0)))
basisaddon=3#1st derivative
### Prepare files for logging evaluation results:
datetime_stmp=str(datetime.datetime.now())
datetime_stmp=datetime_stmp.replace(" ", "_")
datetime_stmp=datetime_stmp.replace(":", "-")
datetime_stmp=datetime_stmp.split('.')[0]
logfile=resdir+datetime_stmp+"_iterations_"+str(num_iterations)+"_ML_classification_Performance_cv.csv"
f=open(logfile,"a")
f.write("Features,Extraction,Iteration,Fold,Classifier,hyperparameter,auROC,Accuracy,Recall,Precision,F1"+"\n")
f.close()
# 
ycol=["sspg_2_classes"]
unxcol=["subject_id","sspg_2_classes"]
# 
for modelset in model_list:
    for idx,extrmeth in enumerate(extr_list_redim):
        pcm=pcmethods[idx]
        if extrmeth=='timeseries' and modelset=='HOME_CGM_MEAN':
            continue
        # training set
        datatemp=combtab_list_cv[extrmeth][modelset]
        xcol=datatemp.columns.difference(unxcol)
        X_train=np.array(datatemp[xcol])
        y_train=np.array(datatemp[ycol])
        le=LabelEncoder()
        groups=le.fit_transform(datatemp['subject_id'])
        for iteri in range(num_iterations):
            skf=StratifiedKFold(n_splits=num_folds,shuffle=True,random_state=1)
            [group_u,ind_p]=np.unique(groups,return_index=True)
            fold=0
            # random order of hyperparameters
            hyp_index=list(range(num_folds))
            random.shuffle(hyp_index)
            for train_index,validation_index in skf.split(X_train[ind_p,],y_train[ind_p]):
                validation_index_new=np.where(np.isin(groups,group_u[validation_index]))[0]
                train_index_new=np.where(np.isin(groups,group_u[train_index]))[0]
                x_trans=dict()
                grplist={"train": le.inverse_transform(groups[train_index_new]),"validation": le.inverse_transform(groups[validation_index_new])}
                if extrmeth=="timeseries":
                    x_train_tmp,x_validation_tmp=X_train[train_index_new],X_train[validation_index_new]
                    x_tmp={'train': x_train_tmp[:,timeind],'validation': x_validation_tmp[:,timeind]}
                    if pcm=="regular":
                        r_pca=sk.decomposition.PCA()
                        principalComponents=dict()
                        principalComponents["train"]=r_pca.fit_transform(x_tmp["train"])
                        principalComponents["validation"]=r_pca.transform(x_tmp["validation"])
                        for setname in setnames:
                            principalDf=pd.DataFrame(data=principalComponents[setname])
                            principalDf_2pcs=principalDf.loc[:,[0,1]]
                            principalDf_2pcs.columns=["PC1","PC2"]
                            x_trans[setname]=principalDf_2pcs
                    elif pcm=="fpca":
                        score_dic=dict()
                        nxvec=x_tmp["train"].shape[1]
                        grid_points=np.linspace(0,1,nxvec)
                        fd_tr=FDataGrid(x_tmp["train"],grid_points)
                        nbasis=nxvec+basisaddon
                        basis_fd=BSplineBasis(n_basis=nbasis)
                        grid=SmoothingParameterSearch(BasisSmoother(basis_fd,regularization=L2Regularization(LinearDifferentialOperator(3))),lambdavec,param_name='smoothing_parameter')
                        _=grid.fit(fd_tr)
                        lambda_sele=grid.best_params_['smoothing_parameter']
                        fpca=FPCA(n_components=2)
                        # 
                        fd=dict()
                        fd["train"]=fd_tr
                        fd["validation"]=FDataGrid(x_tmp["validation"],grid_points)
                        # fpca.components_.plot()
                        for setname in setnames:
                            fd_s=grid.transform(fd[setname])
                            if setname=="train":
                                score_dic=fpca.fit_transform(fd_s)
                            else:
                                score_dic=fpca.transform(fd_s)
                            principalDf=pd.DataFrame(data=score_dic)
                            principalDf_2pcs=principalDf.loc[:,[0,1]]
                            principalDf_2pcs.columns=["PC1","PC2"]
                            x_trans[setname]=principalDf_2pcs
                    for setele in ["train","validation"]:
                        temptab=pd.DataFrame(data=x_trans[setele])
                        temptab['subject_id']=grplist[setele]
                        temptab=temptab.set_index('subject_id')
                        locdf=temptab.merge(metabolic_indicators,left_on=temptab.index,right_on='subject_id')
                        locdf=locdf.drop(['subject_id','sspg_2_classes'],axis=1)
                        x_trans[setele]=np.array(locdf)
                else:
                    x_trans["train"],x_trans["validation"]=X_train[train_index_new],X_train[validation_index_new]
                y=dict()
                y["train"],y["validation"]=y_train[train_index_new].ravel(),y_train[validation_index_new].ravel()
                for index,(name,classifier) in enumerate(classifiers.items()):
                    print("### Features=",modelset,"  Extraction=",extrmeth,"  pca =",pcm,"   Iteration=",iteri,"   Fold=",fold,"   classifier=",name)
                    tun_parameter=tune_para[paraind[index]]
                    if tun_parameter=='C':
                        value_para=Cset[name][hyp_index[fold]]
                    else:
                        value_para=parameters_grid[tun_parameter][hyp_index[fold]]
                    kwargs={tun_parameter: value_para}
                    hyperpara_str=tun_parameter+'='+str(value_para)
                    classifier.set_params(**kwargs)
                    le2=LabelEncoder()
                    classifier.fit(x_trans["train"],le2.fit_transform(y["train"]))
                    y_pred=dict()
                    y_pred["train"]=le2.inverse_transform(classifier.predict(x_trans["train"]))
                    y_pred["validation"]=le2.inverse_transform(classifier.predict(x_trans["validation"]))
                    ### Accuracy, AUC, Recall, Precision, F1
                    accuracy=metrics.accuracy_score(y["validation"],y_pred["validation"])
                    recallflag=len(np.unique(y["validation"]))==1
                    precflag=len(np.unique(y_pred["validation"]))==1
                    probs=classifier.predict_proba(x_trans["validation"])
                    if recallflag:
                        recall=float("nan")
                        classifier_auc=float("nan")
                    else:
                        if probs.shape[1]==2:
                            classifier_auc=metrics.roc_auc_score(y["validation"],probs[:,1],multi_class='ovr')
                            recall=metrics.recall_score(y["validation"],y_pred["validation"],pos_label="IR")
                        else:
                            classifier_auc=metrics.roc_auc_score(y["validation"],probs,multi_class='ovr')
                            recall=metrics.recall_score(y["validation"],y_pred["validation"],average="micro")
                    if precflag:
                        precision=float("nan")
                    else:
                        if probs.shape[1]==2:
                            precision=metrics.precision_score(y["validation"],y_pred["validation"],pos_label="IR")
                        else:
                            precision=metrics.precision_score(y["validation"],y_pred["validation"],average="micro")
                    if recallflag or precflag:
                        f1=float("nan")
                    else:
                        if probs.shape[1]==2:
                            f1=metrics.f1_score(y["validation"],y_pred["validation"],pos_label="IR")
                        else:
                            f1=metrics.f1_score(y["validation"],y_pred["validation"],average="micro")
                    ## Logging evaluation results
                    f = open(logfile,"a")
                    f.write(modelset + "," + extrmeth + "_" + pcm + "," + str(iteri) + "," + str(fold) + "," + name + "," + hyperpara_str + "," + str(round(classifier_auc,2)) + ","+ str(round(accuracy,2)) +","+ str(round(recall,2)) +","+ str(round(precision,2)) +","+ str(round(f1,2)) +"\n")
                    f.close()
                #
                fold += 1

