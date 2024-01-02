import os
import numpy as np
import pandas as pd
import datetime
import random
# 
import sklearn as sk
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import normalize
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
extr_list=["extrfeatures"]#
setnames=["train","validation","test"]
features_list_test=["venous_exp_with_matching","ctru_cgm","home_cgm_1_matching","home_cgm_2_matching"]
model_list=["venous_exp_with_matching","ctru_cgm","CGM_ALL","HOME_CGM_MEAN","HOMECGM1","HOMECGM2","baseline_homa_ir","baseline_demo","baseline_demo_lab"]# Different test set for the model
# Read metabolic indicators 
metabolic_indicators=pd.read_csv(datadir+'metabolic_subphenotyping_all_cohort_metadata_09102023.csv',na_values=['--'])
# refine train-validaiton and test sets (no overlap)
setsep={"train": np.unique(np.array(metabolic_indicators['subject_id'][metabolic_indicators['exp_type']=="venous_without_matching_cgm_and_without_planned_athome_cgm"])),
"test": np.unique(np.array(metabolic_indicators['subject_id'][metabolic_indicators['exp_type']=="venous_with_matching_cgm_and_with_planned_athome_cgm"]))}
setsep['test']=np.setdiff1d(setsep['test'],setsep['train'])
setsep['test']=np.setdiff1d(setsep['test'],np.array(['S112','S115','S122','S105','S108']))
# 
demogrlist=["age","bmi","sex","ethnicity"]
cliniclist=["homa_ir","a1c","fpg"]
metabolic_indicators=metabolic_indicators[["subject_id","sspg_2_classes"]+demogrlist+cliniclist]
for cate in ["sex","ethnicity"]:
    metabolic_indicators=pd.get_dummies(metabolic_indicators,prefix=[cate],columns=[cate], drop_first=True)
addmodelfeat={"baseline_homa_ir": ["homa_ir"],"baseline_demo": ["sex_M","ethnicity_Caucasian","ethnicity_Hispanic"],"baseline_demo_lab": ["sex_M","ethnicity_Caucasian","ethnicity_Hispanic","a1c","fpg"]}
# 
# idremv=['S112','S115']
# metabolic_indicators=metabolic_indicators[~metabolic_indicators['subject_id'].isin(idremv)]
#### load CGM features
# training data
combtab_list_train=dict()
features=pd.read_csv(datadir+"metabolic_subphenotyping_ogtt_glucose_venous_exp_without_matching_features_09102023.csv",index_col=0)
features['ogtt_time_peak_baseline']=features['ogtt_time_peak_baseline'].fillna(180)
# features=features.drop(columns=['ogtt_time_peak_baseline'])
features["ogtt_time_below_basline"]=features["ogtt_time_below_basline"].astype(int)
### Combine extracted features info with metadata
combtab=features.merge(metabolic_indicators,left_on="subject_id",right_on="subject_id")
combtab=combtab[combtab['subject_id'].isin(setsep['train'])]
combtab_list_train["extrfeatures"]=combtab
# test data
metabolic_indicators_noncl=metabolic_indicators.drop(cliniclist,axis=1)
combtab_list_test=dict()
combtab_list_test["extrfeatures"]=dict()
# load data and formulate table
for feature in features_list_test:
    features=pd.read_csv(datadir+"metabolic_subphenotyping_ogtt_glucose_"+feature+"_features_09102023.csv", index_col=0)
    features['ogtt_time_peak_baseline']=features['ogtt_time_peak_baseline'].fillna(180)
    # features=features.drop(columns=['ogtt_time_peak_baseline'])
    features["ogtt_time_below_basline"]=features["ogtt_time_below_basline"].astype(int)
    ### Combine extracted features info with metadata
    combtab=features.merge(metabolic_indicators_noncl,left_on="subject_id",right_on="subject_id")
    combtab=combtab[combtab['subject_id'].isin(setsep['test'])]
    combtab_list_test["extrfeatures"][feature]=combtab
# combined dataset
combtab=combtab_list_test["extrfeatures"]["home_cgm_1_matching"]
combtab=pd.concat([combtab,combtab_list_test["extrfeatures"]["home_cgm_2_matching"]])
combtab=pd.concat([combtab,combtab_list_test["extrfeatures"]["ctru_cgm"]])
combtab_list_test["extrfeatures"]["CGM_ALL"]=combtab
# 
combtab_list_test["extrfeatures"]["HOMECGM1"]=combtab_list_test["extrfeatures"]["home_cgm_1_matching"]
combtab_list_test["extrfeatures"]["HOMECGM2"]=combtab_list_test["extrfeatures"]["home_cgm_2_matching"]
# 
locdf_conc=pd.concat([combtab_list_test["extrfeatures"]["home_cgm_1_matching"],combtab_list_test["extrfeatures"]["home_cgm_2_matching"]])
locdf_mean=locdf_conc.drop(['sspg_2_classes'],axis=1).groupby("subject_id",as_index=False).mean()
locdf_mean=locdf_mean.merge(locdf_conc[['subject_id','sspg_2_classes']],left_on="subject_id",right_on="subject_id")
combtab_list_test["extrfeatures"]["HOME_CGM_MEAN"]=locdf_mean
# 
for key in addmodelfeat:
    subtab=metabolic_indicators[metabolic_indicators["subject_id"].isin(setsep['test'])]
    selelist=addmodelfeat[key]
    combtab_list_test["extrfeatures"][key]=subtab[selelist+["sspg_2_classes"]]
# 
random.seed(10)
num_iterations=100
num_folds=5
C=10
classifiers={
    'L1_logistic': LogisticRegression(C=C,penalty='l1',solver='saga',max_iter=10000,random_state=0),
    'L2_Logistic': LogisticRegression(C=C,penalty='l2',solver='saga',max_iter=10000,random_state=0),
    'Linear_SVC': SVC(kernel='linear',C=C,probability=True,random_state=0),
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
      'Linear_SVC': [0.001,0.005,0.01,0.05,0.1],
      'RBF_SVC': [0.01, 0.1, 1, 10, 100]}
paraind=[0,0,0,0,1,2]
tune_para=list(parameters_grid.keys())
extr_list_redim=["extrfeatures"]
pcmethods=["NAN"]
### Prepare files for logging evaluation results:
datetime_stmp=str(datetime.datetime.now())
datetime_stmp=datetime_stmp.replace(" ", "_")
datetime_stmp=datetime_stmp.replace(":", "-")
datetime_stmp=datetime_stmp.split('.')[0]
logfile=resdir+datetime_stmp+"_iterations_"+str(num_iterations)+"_ML_classification_Performance.csv"
f=open(logfile,"a")
f.write("Features,Extraction,Iteration,Fold,Classifier,hyperparameter,auROC_validation,auROC,Accuracy,Recall,Precision,F1"+"\n")
f.close()
# 
ycol=["sspg_2_classes"]
unxcol=["subject_id","sspg_2_classes"]
# 
for modelset in model_list:
    for idx,extrmeth in enumerate(extr_list_redim):
        pcm=pcmethods[idx]
        # training set
        datatemp=combtab_list_train[extrmeth]
        if modelset in addmodelfeat.keys():
            xcol=addmodelfeat[modelset]
        else:
            xcol=datatemp.columns.difference(unxcol+cliniclist)
        # X_train=normalize(np.array(datatemp[xcol]),axis=0)
        X_train=np.array(datatemp[xcol])
        y_train=np.array(datatemp[ycol])
        le=LabelEncoder()
        groups=le.fit_transform(datatemp['subject_id'])
        # test set for evaluation
        datatemp=combtab_list_test[extrmeth][modelset]
        # X_test=normalize(np.array(datatemp[xcol]),axis=0)
        X_test=np.array(datatemp[xcol])
        y_test=np.array(datatemp[ycol])        
        for iteri in range(num_iterations):
            # skf=StratifiedGroupKFold(n_splits=num_folds,shuffle=True,random_state=None)
            skf=StratifiedKFold(n_splits=num_folds,shuffle=True,random_state=None)
            # skf=KFold(n_splits=num_folds,shuffle=True,random_state=None)
            [group_u,ind_p]=np.unique(groups,return_index=True)
            fold=0
            # random order of hyperparameters
            hyp_index=list(range(num_folds))
            random.shuffle(hyp_index)
            for train_index,validaiton_index in skf.split(X_train[ind_p,],y_train[ind_p]):
                validation_index_new=np.where(np.isin(groups,group_u[validaiton_index]))
                train_index_new=np.where(np.isin(groups,group_u[train_index]))
                validation_index=validation_index_new
                train_index=train_index_new
                x_trans=dict()
                x_trans["train"],x_trans["validation"],x_trans["test"]=X_train[train_index],X_train[validation_index],X_test
                y=dict()
                y["train"],y["validation"],y["test"]=y_train[train_index].ravel(),y_train[validation_index].ravel(),y_test.ravel()
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
                    y_pred["test"]=le2.inverse_transform(classifier.predict(x_trans["test"]))
                    ### Accuracy, AUC, Recall, Precision, F1
                    accuracy=metrics.accuracy_score(y["test"],y_pred["test"])
                    recallflag=len(np.unique(y["test"]))==1
                    precflag=len(np.unique(y_pred["test"]))==1
                    probs_test=classifier.predict_proba(x_trans["test"])
                    if recallflag:
                        recall=float("nan")
                        classifier_auc=float("nan")
                    else:
                        if probs_test.shape[1]==2:
                            classifier_auc=metrics.roc_auc_score(y["test"],probs_test[:,1],multi_class='ovr')
                            recall=metrics.recall_score(y["test"],y_pred["test"],pos_label="IR")
                        else:
                            classifier_auc=metrics.roc_auc_score(y["test"],probs_test,multi_class='ovr')
                            recall=metrics.recall_score(y["test"],y_pred["test"],average="micro")
                    if precflag:
                        precision=float("nan")
                    else:
                        if probs_test.shape[1]==2:
                            precision=metrics.precision_score(y["test"],y_pred["test"],pos_label="IR")
                        else:
                            precision=metrics.precision_score(y["test"],y_pred["test"],average="micro")
                    if recallflag or precflag:
                        f1=float("nan")
                    else:
                        if probs_test.shape[1]==2:
                            f1=metrics.f1_score(y["test"],y_pred["test"],pos_label="IR")
                        else:
                            f1=metrics.f1_score(y["test"],y_pred["test"],average="micro")
                    if len(np.unique(y["validation"]))==1:
                        classifier_auc_validaiton=float("nan")
                    else:
                        probs_test=classifier.predict_proba(x_trans["validation"])
                        if probs_test.shape[1]==2:
                            classifier_auc_validaiton=metrics.roc_auc_score(y["validation"],probs_test[:,1],multi_class='ovr')
                        else:
                            classifier_auc_validaiton=metrics.roc_auc_score(y["validation"],probs_test,multi_class='ovr')
                    ## Logging evaluation results
                    f = open(logfile,"a")
                    f.write(modelset + "," + extrmeth + "_" + pcm + "," + str(iteri) + "," + str(fold) + "," + name + "," + hyperpara_str + ","+ str(round(classifier_auc_validaiton,2)) + "," + str(round(classifier_auc,2)) + ","+ str(round(accuracy,2)) +","+ str(round(recall,2)) +","+ str(round(precision,2)) +","+ str(round(f1,2)) +"\n")
                    f.close()
                #
                fold += 1

