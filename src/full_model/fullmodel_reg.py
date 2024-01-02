import os
import numpy as np
import pandas as pd
import datetime
import random
# 
import sklearn as sk
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import Lasso, Ridge
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor
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
groupsep=80
IR_thres=120
extr_list=["extrfeatures"]
setnames=["train","validation","test"]
features_list_test=["venous_exp_with_matching","ctru_cgm","home_cgm_1_matching","home_cgm_2_matching"]
model_list=["venous_exp_with_matching","ctru_cgm","CGM_ALL"]# Different test set for the model
# Read metabolic indicators 
metabolic_indicators=pd.read_csv(datadir+'metabolic_subphenotyping_all_cohort_metadata_09102023.csv',na_values=['--'])
# refine train-validaiton and test sets (no overlap)
setsep={"train": np.unique(np.array(metabolic_indicators['subject_id'][metabolic_indicators['exp_type']=="venous_without_matching_cgm_and_without_planned_athome_cgm"])),
"test": np.unique(np.array(metabolic_indicators['subject_id'][metabolic_indicators['exp_type']=="venous_with_matching_cgm_and_with_planned_athome_cgm"]))}
setsep['test']=np.setdiff1d(setsep['test'],setsep['train'])
# 
demogrlist=["age","bmi","sex","ethnicity"]
metabolic_indicators=metabolic_indicators[["subject_id","sspg"]+demogrlist]
for cate in ["sex","ethnicity"]:
    metabolic_indicators=pd.get_dummies(metabolic_indicators,prefix=[cate],columns=[cate], drop_first=True)
metabolic_indicators["muscle_ir_status"]=0
metabolic_indicators.loc[metabolic_indicators.sspg>groupsep,"muscle_ir_status"]=1
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
combtab_list_test=dict()
combtab_list_test["extrfeatures"]=dict()
# load data and formulate table
for feature in features_list_test:
    features=pd.read_csv(datadir+"metabolic_subphenotyping_ogtt_glucose_"+feature+"_features_09102023.csv", index_col=0)
    features['ogtt_time_peak_baseline']=features['ogtt_time_peak_baseline'].fillna(180)
    # features=features.drop(columns=['ogtt_time_peak_baseline'])
    features["ogtt_time_below_basline"]=features["ogtt_time_below_basline"].astype(int)
    ### Combine extracted features info with metadata
    combtab=features.merge(metabolic_indicators,left_on="subject_id",right_on="subject_id")
    combtab=combtab[combtab['subject_id'].isin(setsep['test'])]
    combtab_list_test["extrfeatures"][feature]=combtab
# combined dataset
combtab=combtab_list_test["extrfeatures"]["home_cgm_1_matching"]
combtab=pd.concat([combtab,combtab_list_test["extrfeatures"]["home_cgm_2_matching"]])
combtab=pd.concat([combtab,combtab_list_test["extrfeatures"]["ctru_cgm"]])
combtab_list_test["extrfeatures"]["CGM_ALL"]=combtab
# 
random.seed(10)
num_iterations=100
num_folds=5
C=10
alpha=1/C
regressors={
    'L1_linear': Lasso(alpha=alpha,max_iter=10000,random_state=0),
    'L2_linear': Ridge(alpha=alpha,solver='saga',max_iter=10000,random_state=0),
    'Linear_SVR': SVR(kernel='linear',C=C),
    'RBF_SVR': SVR(kernel='rbf',C=C),
    'RFR': RandomForestRegressor(random_state=23234),
    'XGB': xgb.XGBRegressor(learning_rate=1.0)
}
parameters_grid={
    'alpha': [0.01, 0.1, 1, 10, 100],
    'C': [0.01, 0.1, 1, 10, 100],
    'n_estimators': [20, 50, 100, 150, 200],
    'learning_rate': [0.001, 0.01, 0.1, 0.5, 1.0]
}
Cset={'Linear_SVR': [0.0005,0.001,0.005,0.01,0.05],
      'RBF_SVR': [0.01, 0.1, 1, 10, 100]}
paraind=[0,0,1,1,2,3]#[0,0,1,1,2,3]
tune_para=list(parameters_grid.keys())
extr_list_redim=["extrfeatures"]
pcmethods=["NAN"]
### Prepare files for logging evaluation results:
datetime_stmp=str(datetime.datetime.now())
datetime_stmp=datetime_stmp.replace(" ", "_")
datetime_stmp=datetime_stmp.replace(":", "-")
datetime_stmp=datetime_stmp.split('.')[0]
logfile=resdir+datetime_stmp+"_iterations_"+str(num_iterations)+"_ML_regression_Performance.csv"
f=open(logfile,"a")
f.write("Features,Extraction,Iteration,Fold,Regressor,hyperparameter,RRMSE_validation,RRMSE,Correlation,Accuracy,Precision,Recall,F1"+"\n")
f.close()
# 
ycol=["sspg"]
ygroupcol=["muscle_ir_status"]
unxcol=["subject_id","sspg","muscle_ir_status"]
# 
for modelset in model_list:
    for idx,extrmeth in enumerate(extr_list_redim):
        pcm=pcmethods[idx]
        # training set
        datatemp=combtab_list_train[extrmeth]
        xcol=datatemp.columns.difference(unxcol)
        # X_train=normalize(np.array(datatemp[xcol]),axis=0)
        X_train=np.array(datatemp[xcol])
        y_train=np.array(datatemp[ycol])
        ygroups=np.array(datatemp[ygroupcol])
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
            for train_index,validaiton_index in skf.split(X_train[ind_p,],ygroups[ind_p]):
                validation_index_new=np.where(np.isin(groups,group_u[validaiton_index]))
                train_index_new=np.where(np.isin(groups,group_u[train_index]))
                validation_index=validation_index_new
                train_index=train_index_new
                x_trans=dict()
                x_trans["train"],x_trans["validation"],x_trans["test"]=X_train[train_index],X_train[validation_index],X_test
                y=dict()
                y["train"],y["validation"],y["test"]=y_train[train_index].ravel(),y_train[validation_index].ravel(),y_test.ravel()
                for index,(name,regressor) in enumerate(regressors.items()):
                    print("### Features=",modelset,"  Extraction=",extrmeth,"  pca =",pcm,"   Iteration=",iteri,"   Fold=",fold,"   regressor=",name)
                    tun_parameter=tune_para[paraind[index]]
                    if tun_parameter=='C':
                        value_para=Cset[name][hyp_index[fold]]
                    else:
                        value_para=parameters_grid[tun_parameter][hyp_index[fold]]
                    kwargs={tun_parameter: value_para}
                    hyperpara_str=tun_parameter+'='+str(value_para)
                    regressor.set_params(**kwargs)
                    regressor.fit(x_trans["train"],y["train"])
                    y_pred=dict()
                    y_pred["train"]=regressor.predict(x_trans["train"])
                    y_pred["validation"]=regressor.predict(x_trans["validation"])
                    y_pred["test"]=regressor.predict(x_trans["test"])
                    ### RMSE and R2
                    corr=np.corrcoef(y["test"],y_pred["test"])[0,1]
                    maxscal=np.concatenate((y["test"],y_pred["test"]),axis=None).max()
                    rrmse=metrics.mean_squared_error(y["test"]/maxscal,y_pred["test"]/maxscal)**0.5
                    maxscal=np.concatenate((y["validation"],y_pred["validation"]),axis=None).max()
                    rrmse_validation=metrics.mean_squared_error(y["validation"]/maxscal,y_pred["validation"]/maxscal)**0.5
                    accuracy=metrics.accuracy_score(y["test"]>IR_thres,y_pred["test"]>IR_thres)
                    y_test_clas=y["test"]>IR_thres
                    y_pred_test_clas=y_pred["test"]>IR_thres
                    recallflag=len(np.unique(y_test_clas))==1
                    precflag=len(np.unique(y_pred_test_clas))==1
                    if recallflag:
                        recall=float("nan")
                        classifier_auc=float("nan")
                    else:
                        recall=metrics.recall_score(y_test_clas,y_pred_test_clas)
                    if precflag:
                        precision=float("nan")
                    else:
                        precision=metrics.precision_score(y_test_clas,y_pred_test_clas)
                    if recallflag or precflag:
                        f1=float("nan")
                    else:
                        f1=metrics.f1_score(y_test_clas,y_pred_test_clas)
                    ## Logging evaluation results
                    f = open(logfile,"a")
                    f.write(modelset + "," + extrmeth + "_" + pcm + "," + str(iteri) + "," + str(fold) + "," + name + "," + hyperpara_str + ","+ str(round(rrmse_validation,2)) + "," + str(round(rrmse,2)) + ","+ str(round(corr,2)) +","+str(round(accuracy,2)) +","+ str(round(precision,2)) +","+ str(round(recall,2)) +","+ str(round(f1,2)) +"\n")
                    f.close()
                #
                fold += 1
