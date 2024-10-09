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
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
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
datadir=pardir+"data_Metabolic_subphenotype_and_CGM_analysis/"
resdir=pardir+"res/cgm_ogtt/"
os.chdir(resdir)
#### Read metabolic indicators, timeseries, and features
metabolic_indicators=pd.read_csv(datadir+"dataupdate/validation_cohort_processed_sspg_01192023.csv",na_values=['--'])# Prepare annotation columns for the regressers
metabolic_indicators=metabolic_indicators[metabolic_indicators["SSPG"].notna()]
metabolic_indicators["muscle_ir_status"]=0
metabolic_indicators.loc[metabolic_indicators.SSPG>80,"muscle_ir_status"]=1
# 
features_list=["CTRU_Venous","CTRU_CGM","Home_CGM_1","Home_CGM_2"]
model_list=[ele for ele in features_list]
model_list.extend(["HOME_CGM_ALL","CGM_ALL","HOME_CGM_MEAN"])
extr_list=["timeseries","extrfeatures"]
combtab_list=dict()
combtab_list["timeseries"]=dict()
combtab_list["extrfeatures"]=dict()
# load data and formulate table
for feature in features_list:
    timeseries=pd.read_csv(resdir+"timeseries."+feature+".csv",index_col=0)
    timeseries=timeseries.transpose()
    features=pd.read_csv(resdir+"feature."+feature+".csv", index_col=0)
    ## Combine timeseries with metadata
    combtab=timeseries.merge(metabolic_indicators,left_on=timeseries.index,right_on='Subject')
    # combtab.index=combtab["Subject"]
    combtab_list["timeseries"][feature]=combtab
    ### Combine extracted features info with metadata
    combtab=features.merge(metabolic_indicators,left_on=features.index,right_on="Subject")
    combtab_list["extrfeatures"][feature]=combtab
# combined dataset
for extr in extr_list:
    combtab=combtab_list[extr]["Home_CGM_1"]
    combtab=pd.concat([combtab,combtab_list[extr]["Home_CGM_2"]])
    combtab_list[extr]["HOME_CGM_ALL"]=combtab
    combtab=pd.concat([combtab,combtab_list[extr]["CTRU_CGM"]])
    combtab_list[extr]["CGM_ALL"]=combtab
    locdf1=combtab_list[extr]["Home_CGM_1"]
    locdf2=combtab_list[extr]["Home_CGM_2"]
    locdf_conc=pd.concat([locdf1,locdf2])
    locdf_mean=locdf_conc.groupby("Subject",as_index=False).mean()
    combtab_list[extr]["HOME_CGM_MEAN"]=locdf_mean

random.seed(10)
num_iterations=20
num_folds=4
C=10
alpha=1/C
kernel=1.0*RBF(1)
regressers={
    'L1_linear': Lasso(alpha=alpha,max_iter=10000),
    'L2_linear': Ridge(alpha=alpha,solver='saga',max_iter=10000),
    'Linear_SVR': SVR(kernel='linear',C=C),
    'RBF_SVR': SVR(kernel='rbf',C=C),
    'GPR': GaussianProcessRegressor(kernel,random_state=23234),
    'RFR': RandomForestRegressor(random_state=23234),
    'gboostR': GradientBoostingRegressor(random_state=23234)
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
logfile=resdir+datetime_stmp+"_iterations_"+str(num_iterations)+"_ML_regression_Performance.csv"
f=open(logfile,"a")
f.write("Features,Extraction,Iteration,Fold,Regresser,RRMSE,Correlation,Accuracy,Precision,Recall,F1"+"\n")
f.close()
# 
ycol=["SSPG"]
ygroupcol=["muscle_ir_status"]
unxcol=["Subject","SSPG","muscle_ir_status"]
IR_thres=120
# 
for modelset in model_list:
    for idx,extrmeth in enumerate(extr_list_redim):
        pcm=pcmethods[idx]
        datatemp=combtab_list[extrmeth][modelset]
        xcol=datatemp.columns.difference(unxcol)
        X=np.array(datatemp[xcol])
        y=np.array(datatemp[ycol])
        ygroups=np.array(datatemp[ygroupcol])
        le=LabelEncoder()
        groups=le.fit_transform(datatemp['Subject'])
        for iteri in range(num_iterations):
            # skf=StratifiedGroupKFold(n_splits=num_folds,shuffle=True,random_state=None)
            skf=StratifiedKFold(n_splits=num_folds,shuffle=True,random_state=None)
            # skf=KFold(n_splits=num_folds,shuffle=True,random_state=None)
            [group_u,ind_p]=np.unique(groups,return_index=True)
            fold=0
            for train_index,test_index in skf.split(X[ind_p,],ygroups[ind_p]):
                test_index_new=np.where(np.isin(groups,group_u[test_index]))
                train_index_new=np.where(np.isin(groups,group_u[train_index]))
                test_index=test_index_new
                train_index=train_index_new
                fold += 1
                if extrmeth=="timeseries":
                    x_train_tmp,x_test_tmp=X[train_index],X[test_index]
                    if pcm=="regular":
                        r_pca=sk.decomposition.PCA()
                        principalComponents=r_pca.fit_transform(x_train_tmp)
                        principalDf=pd.DataFrame(data=principalComponents)
                        principalDf_2pcs=principalDf.loc[:,[0,1]]
                        principalDf_2pcs.columns =["PC1","PC2"]
                        x_train=principalDf_2pcs
                        test_principalComponents=r_pca.transform(x_test_tmp)
                        test_principalDf=pd.DataFrame(data=test_principalComponents)
                        test_principalDf_2pcs=test_principalDf.loc[:,[0,1]]
                        test_principalDf_2pcs.columns=["PC1","PC2"]
                        x_test=test_principalDf_2pcs
                    elif pcm=="fpca":
                        score_dic=dict()
                        nxvec=x_train_tmp.shape[1]
                        grid_points=np.linspace(0,1,nxvec)
                        fd_tr=FDataGrid(x_train_tmp,grid_points)
                        nbasis=nxvec+basisaddon
                        basis_fd=BSplineBasis(n_basis=nbasis)
                        grid=SmoothingParameterSearch(BasisSmoother(basis_fd,regularization=L2Regularization(LinearDifferentialOperator(3))),lambdavec,param_name='smoothing_parameter')
                        _=grid.fit(fd_tr)
                        lambda_sele=grid.best_params_['smoothing_parameter']
                        fd_tr_s=grid.transform(fd_tr)
                        fpca=FPCA(n_components=2)
                        score_dic["train"]=fpca.fit_transform(fd_tr_s)
                        # fpca.components_.plot()
                        # test data
                        fd_test=FDataGrid(x_test_tmp,grid_points)
                        basis_fd_s_test=grid.transform(fd_test)
                        score_dic["test"]=fpca.transform(basis_fd_s_test)
                        for setsep in ["train","test"]:
                            principalDf=pd.DataFrame(data=score_dic[setsep])
                            principalDf_2pcs=principalDf.loc[:,[0,1]]
                            principalDf_2pcs.columns=["PC1","PC2"]
                            globals()['x_'+setsep]=principalDf_2pcs
                else:
                    x_train, x_test = X[train_index],X[test_index]
                y_train,y_test = y[train_index].ravel(),y[test_index].ravel()
                # 
                for index,(name,regresser) in enumerate(regressers.items()):
                    print("### Features=",modelset,"  Extraction=",extrmeth,"  pca =",pcm,"   Iteration=",iteri,"   Fold=",fold,"   regresser=",name)
                    regresser.fit(x_train,y_train)
                    y_pred_train=regresser.predict(x_train)
                    y_pred_test=regresser.predict(x_test)
                    ### RMSE and R2
                    # r2=metrics.r2_score(y_test,y_pred_test)
                    corr=np.corrcoef(y_test,y_pred_test)[0,1]
                    maxscal=np.concatenate((y_test,y_pred_test),axis=None).max()
                    rrmse=metrics.mean_squared_error(y_test/maxscal,y_pred_test/maxscal)**0.5
                    accuracy=metrics.accuracy_score(y_test>IR_thres,y_pred_test>IR_thres)
                    y_test_clas=y_test>IR_thres
                    y_pred_test_clas=y_pred_test>IR_thres
                    recallflag=len(np.unique(y_test_clas))==1
                    precflag=len(np.unique(y_pred_test_clas))==1
                    if recallflag:
                        recall=float("nan")
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
                    f.write(modelset + "," + extrmeth + "_" + pcm + "," + str(iteri) + "," + str(fold) + "," + name + "," + str(round(rrmse,2)) + ","+ str(round(corr,2)) + ","+ str(round(accuracy,2)) + ","+ str(round(precision,2)) + ","+ str(round(recall,2)) + ","+ str(round(f1,2)) +"\n")
                    f.close()

