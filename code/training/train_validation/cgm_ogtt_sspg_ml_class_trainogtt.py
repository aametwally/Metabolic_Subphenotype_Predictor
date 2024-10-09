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
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn import metrics
from sklearn.preprocessing import LabelEncoder
import skfda
from skfda.representation.basis import BSplineBasis
from skfda.preprocessing.smoothing import BasisSmoother
from skfda.representation.grid import FDataGrid
from skfda.preprocessing.dim_reduction import FPCA
from skfda.preprocessing.smoothing.validation import SmoothingParameterSearch
from skfda.misc.regularization import L2Regularization
from skfda.misc.operators import LinearDifferentialOperator
# 
import matplotlib.pyplot as plt
## Set project Path and get datetimestamp
pardir="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_dataset/"
datadir=pardir+"data_Metabolic_subphenotype_and_CGM_analysis/dataupdate/"
resdir=pardir+"res/cgm_ogtt/"
os.chdir(resdir)
IR_thres=120
filelist={"training": "aggregated_metabolic_indicators__4_predictions_w_prs_01042022.csv", "test": "validation_cohort_processed_sspg_01192023.csv"}
extr_list=["timeseries","extrfeatures"]
setnames=["train","validation","test"]
features_list_test=["CTRU_Venous","CTRU_CGM","Home_CGM_1","Home_CGM_2"]
model_list=[ele for ele in features_list_test]
model_list.extend(["HOME_CGM_ALL","CGM_ALL","HOME_CGM_MEAN"])
# Read metabolic indicators 
metabolic_indicators_dict=dict()
for dataset in filelist.keys():
    metabolic_indicators=pd.read_csv(datadir+filelist[dataset],na_values=['--'])# Prepare annotation columns for the regressers
    metabolic_indicators=metabolic_indicators[metabolic_indicators["SSPG"].notna()]
    metabolic_indicators["SSPG"]=metabolic_indicators["SSPG"].astype(float)
    metabolic_indicators=metabolic_indicators[["Subject","SSPG"]]
    metabolic_indicators["muscle_ir_status"]=0
    metabolic_indicators.loc[metabolic_indicators.SSPG>IR_thres,"muscle_ir_status"]=1
    metabolic_indicators_dict[dataset]=metabolic_indicators
#### load timeseries, and features
# training data
combtab_list_train=dict()
timeseries=pd.read_csv(resdir+"timeseries.CTRU_Venous.train.csv",index_col=0)
timeseries=timeseries.transpose()
features=pd.read_csv(resdir+"feature.CTRU_Venous.train.csv", index_col=0)
## Combine timeseries with metadata
combtab=timeseries.merge(metabolic_indicators_dict["training"],left_on=timeseries.index,right_on='Subject')
combtab_list_train["timeseries"]=combtab
### Combine extracted features info with metadata
combtab=features.merge(metabolic_indicators_dict["training"],left_on=features.index,right_on="Subject")
combtab_list_train["extrfeatures"]=combtab
# test data
combtab_list_test=dict()
combtab_list_test["timeseries"]=dict()
combtab_list_test["extrfeatures"]=dict()
# load data and formulate table
for feature in features_list_test:
    timeseries=pd.read_csv(resdir+"timeseries."+feature+".csv",index_col=0)
    timeseries=timeseries.transpose()
    features=pd.read_csv(resdir+"feature."+feature+".csv", index_col=0)
    ## Combine timeseries with metadata
    combtab=timeseries.merge(metabolic_indicators_dict["test"],left_on=timeseries.index,right_on='Subject')
    # combtab.index=combtab["Subject"]
    combtab_list_test["timeseries"][feature]=combtab
    ### Combine extracted features info with metadata
    combtab=features.merge(metabolic_indicators_dict["test"],left_on=features.index,right_on="Subject")
    combtab_list_test["extrfeatures"][feature]=combtab
# combined dataset
for extr in extr_list:
    combtab=combtab_list_test[extr]["Home_CGM_1"]
    combtab=pd.concat([combtab,combtab_list_test[extr]["Home_CGM_2"]])
    combtab_list_test[extr]["HOME_CGM_ALL"]=combtab
    combtab=pd.concat([combtab,combtab_list_test[extr]["CTRU_CGM"]])
    combtab_list_test[extr]["CGM_ALL"]=combtab
    locdf1=combtab_list_test[extr]["Home_CGM_1"]
    locdf2=combtab_list_test[extr]["Home_CGM_2"]
    locdf_conc=pd.concat([locdf1,locdf2])
    locdf_mean=locdf_conc.groupby("Subject",as_index=False).mean()
    combtab_list_test[extr]["HOME_CGM_MEAN"]=locdf_mean

random.seed(10)
num_iterations=20
num_folds=4
C=10
kernel=1.0*RBF(1)
classifiers={
    'L1_logistic': LogisticRegression(C=C,penalty='l1',solver='saga',max_iter=10000),
    'L2_Logistic': LogisticRegression(C=C,penalty='l2',solver='saga',max_iter=10000),
    'Linear_SVC': SVC(kernel='linear', C=C,probability=True,random_state=0),
    'RBF_SVC': SVC(kernel='rbf',C=C,probability=True,random_state=0),
    'GPC': GaussianProcessClassifier(kernel,warm_start=True,random_state=23234)
}
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
logfile=resdir+datetime_stmp+"_iterations_"+str(num_iterations)+"_ML_classification_Performance.csv"
f=open(logfile,"a")
f.write("Features,Extraction,Iteration,Fold,Classifier,auROC_validation,auROC,Accuracy,Recall,Precision,F1"+"\n")
f.close()
# 
ycol=["muscle_ir_status"]
unxcol=["Subject","SSPG","muscle_ir_status"]
# 
for modelset in model_list:
    for idx,extrmeth in enumerate(extr_list_redim):
        pcm=pcmethods[idx]
        # training set
        datatemp=combtab_list_train[extrmeth]
        xcol=datatemp.columns.difference(unxcol)
        X_train=np.array(datatemp[xcol])
        y_train=np.array(datatemp[ycol])
        le=LabelEncoder()
        groups=le.fit_transform(datatemp['Subject'])
        # test set
        datatemp=combtab_list_test[extrmeth][modelset]
        X_test=np.array(datatemp[xcol])
        y_test=np.array(datatemp[ycol])        
        for iteri in range(num_iterations):
            # skf=StratifiedGroupKFold(n_splits=num_folds,shuffle=True,random_state=None)
            skf=StratifiedKFold(n_splits=num_folds,shuffle=True,random_state=None)
            # skf=KFold(n_splits=num_folds,shuffle=True,random_state=None)
            [group_u,ind_p]=np.unique(groups,return_index=True)
            fold=0
            for train_index,validaiton_index in skf.split(X_train[ind_p,],y_train[ind_p]):
                validation_index_new=np.where(np.isin(groups,group_u[validaiton_index]))
                train_index_new=np.where(np.isin(groups,group_u[train_index]))
                validation_index=validation_index_new
                train_index=train_index_new
                fold += 1
                if extrmeth=="timeseries":
                    x_train_tmp,x_validation_tmp,x_test_tmp=X_train[train_index],X_train[validation_index],X_test
                    x_tmp={'train': x_train_tmp,'validation': x_validation_tmp,'test': x_test_tmp}
                    x_trans=dict()
                    if pcm=="regular":
                        r_pca=sk.decomposition.PCA()
                        principalComponents=dict()
                        principalComponents["train"]=r_pca.fit_transform(x_tmp["train"])
                        principalComponents["validation"]=r_pca.transform(x_tmp["validation"])
                        principalComponents["test"]=r_pca.transform(x_tmp["test"])
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
                        fd["test"]=FDataGrid(x_tmp["test"],grid_points)
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
                else:
                    x_trans["train"],x_trans["validation"],x_trans["test"]=X_train[train_index],X_train[validation_index],X_test
                y=dict()
                y["train"],y["validation"],y["test"]=y_train[train_index].ravel(),y_train[validation_index].ravel(),y_test.ravel()
                # 
                for index,(name,classifier) in enumerate(classifiers.items()):
                    print("### Features=",modelset,"  Extraction=",extrmeth,"  pca =",pcm,"   Iteration=",iteri,"   Fold=",fold,"   classifier=",name)
                    # x_trans
                    classifier.fit(x_trans["train"],y["train"])
                    y_pred=dict()
                    y_pred["train"]=classifier.predict(x_trans["train"])
                    y_pred["validation"]=classifier.predict(x_trans["validation"])
                    y_pred["test"]=classifier.predict(x_trans["test"])
                    ### Accuracy, AUC, Recall, Precision, F1
                    accuracy=metrics.accuracy_score(y["test"],y_pred["test"])
                    recallflag=len(np.unique(y["test"]))==1
                    precflag=len(np.unique(y_pred["test"]))==1
                    if recallflag:
                        recall=float("nan")
                        classifier_auc=float("nan")
                    else:
                        recall=metrics.recall_score(y["test"],y_pred["test"])
                        # roc curve and auc
                        probs_test=classifier.predict_proba(x_trans["test"])
                        classifier_auc=metrics.roc_auc_score(y["test"],probs_test[:, 1],multi_class='ovr')
                    if precflag:
                        precision=float("nan")
                    else:
                        precision=metrics.precision_score(y["test"],y_pred["test"])
                    if recallflag or precflag:
                        f1=float("nan")
                    else:
                        f1=metrics.f1_score(y["test"],y_pred["test"])
                    if len(np.unique(y["validation"]))==1:
                        classifier_auc_validaiton=float("nan")
                    else:
                        probs_test=classifier.predict_proba(x_trans["validation"])
                        classifier_auc_validaiton=metrics.roc_auc_score(y["validation"],probs_test[:, 1],multi_class='ovr')
                    ## Logging evaluation results
                    f = open(logfile,"a")
                    f.write(modelset + "," + extrmeth + "_" + pcm + "," + str(iteri) + "," + str(fold) + "," + name + "," + str(round(classifier_auc_validaiton,2)) + "," + str(round(classifier_auc,2)) + ","+ str(round(accuracy,2)) +","+ str(round(recall,2)) +","+ str(round(precision,2)) +","+ str(round(f1,2)) +"\n")
                    f.close()

