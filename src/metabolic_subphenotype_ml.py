import os
import numpy as np
import pandas as pd
import datetime
import random

import sklearn as sk
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn import metrics


## Set project Path and get datetimestamp
project_path = "/Users/ahmedm/Library/CloudStorage/Box-Box/Ahmed Metwally's Files/Stanford/cgm/43883/Manuscripts/Manuscript_1_MetabolicSubphenotyping/"
os.chdir(project_path)

#### Read metabolic indicators, incretin, ogtt features, and ogtt timeseries
metabolic_indicators = pd.read_csv("aggregated_metabolic_indicators__4_predictions_w_prs_01042022.csv")

ogtt_timeseries = pd.read_csv("ogtt_imputed_normalized_09102021.csv", index_col = 0)
ogtt_timeseries = ogtt_timeseries.transpose()
ogtt_features = pd.read_csv("ogtt_features_09102021.csv", index_col=0)

## Select the 31 subjects in our study 
ogtt_timeseries = ogtt_timeseries[ogtt_timeseries.index.isin(metabolic_indicators.SubjectID)]

## Prepare annotation columns for the classifiers
metabolic_indicators["ethnicity_numeric"] = 0
metabolic_indicators.ethnicity_numeric[metabolic_indicators.ethnicity=="Asian"] = 1
metabolic_indicators.ethnicity_numeric[metabolic_indicators.ethnicity=="Hispanic"] = 2

metabolic_indicators["muscle_ir_status"] = 0
metabolic_indicators.muscle_ir_status[metabolic_indicators.SSPG_class_2=="IR"] = 1

metabolic_indicators["ie_status"] = 0
metabolic_indicators.ie_status[metabolic_indicators.ie_3_classes=="Intermediate"] = 1
metabolic_indicators.ie_status[metabolic_indicators.ie_3_classes=="Dysfunction"] = 2

metabolic_indicators["di_status"] = 0
metabolic_indicators.di_status[metabolic_indicators.di_3classes=="Intermediate"] = 1
metabolic_indicators.di_status[metabolic_indicators.di_3classes=="Dysfunction"] = 2

metabolic_indicators["di_status_v2"] = 0
metabolic_indicators.di_status_v2[metabolic_indicators.di_3classes_v2=="Intermediate"] = 1
metabolic_indicators.di_status_v2[metabolic_indicators.di_3classes_v2=="Dysfunction"] = 2


metabolic_indicators["hepatic_ir_status"] = 0
metabolic_indicators.hepatic_ir_status[metabolic_indicators.hepatic_ir_3classes=="Intermediate"]=1
metabolic_indicators.hepatic_ir_status[metabolic_indicators.hepatic_ir_3classes=="IR"]=2

metabolic_indicators["hepatic_ir_status_v2"] = 0
metabolic_indicators.hepatic_ir_status_v2[metabolic_indicators.hepatic_ir_3classes_v2=="Intermediate"]=1
metabolic_indicators.hepatic_ir_status_v2[metabolic_indicators.hepatic_ir_3classes_v2=="IR"]=2


## Combine OGGT glucose timeseries with metadata
metabolic_indicators_ogtt_timeseries = ogtt_timeseries.merge(metabolic_indicators, left_on=ogtt_timeseries.index, right_on='SubjectID')
metabolic_indicators_ogtt_timeseries.index = metabolic_indicators_ogtt_timeseries["SubjectID"]
metabolic_indicators_ogtt_timeseries.head()

### Combine ogtt features info with metadata
metabolic_indicators_ogtt_features = metabolic_indicators.merge(ogtt_features, left_on='SubjectID', right_on=ogtt_features.index)

### Combine ogtt features info wih timeseries with metadata
metabolic_indicators_ogtt_features_timeseries = metabolic_indicators_ogtt_features.merge(ogtt_timeseries, left_on='SubjectID', right_on=ogtt_timeseries.index)



## Prediction of metabolic Subphenotype (w/ gold standard comparison)
subphenotypes = ["muscle_ir_status", "di_status", "ie_status", "hepatic_ir_status", 'hepatic_ir_status_v2']
features_list = ["demographics_lab", "demographics", "lab", "homair", "homab", "matsuda", "incretins", "ogtt_timeseries", "ogtt_features", "PRS"]

random.seed(10)
num_iterations = 20
num_folds = 4
C = 10
kernel = 1.0 * RBF(1)
classifiers = {
    'L1_logistic': LogisticRegression(C=C, penalty='l1',
                                      solver='saga',
                                      #multi_class='multinomial',
                                      max_iter=10000),
    'L2_Logistic': LogisticRegression(C=C, penalty='l2',
                                                    solver='saga',
                                                    #multi_class='multinomial',
                                                    max_iter=10000),
    'Linear_SVC': SVC(kernel='linear', C=C, probability=True,
                      random_state=0),
    'RBF_SVC': SVC(kernel='rbf', C=C, probability=True,
                      random_state=0),#,
    'GPC': GaussianProcessClassifier(kernel, warm_start = True, random_state = 23234)
}


### Prepare files for logging evaluation results:
datetime_stmp = str(datetime.datetime.now())
datetime_stmp = datetime_stmp.replace(" ", "_")
datetime_stmp = datetime_stmp.replace(":", "-")
datetime_stmp = datetime_stmp.split('.')[0]


f = open("/Users/ahmedm/Documents/" + datetime_stmp + "_iterations_"+ str(num_iterations) + "_ML_Performance.csv", "a")
f.write("Features,Subphenotype,Iteration,Fold,Classifier,auROC,Accuracy,Recall_macro,Precision_macro,F1_macro,Recall_micro,Precision_micro,F1_micro" + "\n")
f.close()
for features in features_list:
    for subphenotype in subphenotypes:
        if features == "demographics_lab": 
            tmp = metabolic_indicators[["age", "sex_01", "bmi", "PreDM_T2D_history", "ethnicity_numeric", "a1c", "fpg", subphenotype]]
            tmp = tmp.dropna(axis=0)
            X = tmp[["age", "sex_01", "bmi", "PreDM_T2D_history", "ethnicity_numeric", "a1c", "fpg"]]
            y = tmp[subphenotype] 
        elif features == "demographics": 
            tmp = metabolic_indicators[["age", "sex_01", "bmi", "PreDM_T2D_history", "ethnicity_numeric", subphenotype]]
            tmp = tmp.dropna(axis=0)
            X = tmp[["age", "sex_01", "bmi", "PreDM_T2D_history", "ethnicity_numeric"]]
            y = tmp[subphenotype]
        elif features == "lab": 
            tmp = metabolic_indicators[["a1c", "fpg", subphenotype]]
            tmp = tmp.dropna(axis=0)
            X = tmp[["a1c", "fpg"]]
            y = tmp[subphenotype] 
        elif features == "homair":
            X = metabolic_indicators[["HOMA_IR"]]
            y = metabolic_indicators[subphenotype]
        elif features == "homab":
            X = metabolic_indicators[["HOMA_B"]]
            y = metabolic_indicators[subphenotype]
        elif features == "matsuda":
            X = metabolic_indicators[["matsuda_index"]]
            y = metabolic_indicators[subphenotype]
            
        elif features == "PRS":
            print("IN PRS")
            X = metabolic_indicators[["OR"]]
            y = metabolic_indicators[subphenotype]   
            
        elif features == "incretins":
            print("in incretins")        
            tmp = metabolic_indicators[['GLP1_120min', 'GIP_120min', subphenotype]]
            tmp = tmp.dropna(axis=0)
            X = tmp[['GLP1_120min', 'GIP_120min']]
            y = tmp[[subphenotype]]
        elif features == "ogtt_features":
            X = metabolic_indicators_ogtt_features[['ogtt_fpg', 'ogtt_60', 'ogtt_120',
                     'ogtt_180', 'ogtt_auc', 'ogtt_iauc', 'ogtt_pauc', 'ogtt_nauc',
                     'ogtt_max', 'ogtt_curve_size', 'ogtt_cv', 'ogtt_time_baseline_peak',
                     'ogtt_slope_baseline_peak','ogtt_slope_peak_last']]
            y = metabolic_indicators_ogtt_features[subphenotype]
        elif features == "ogtt_timeseries":
            X = metabolic_indicators_ogtt_timeseries.iloc[:,0:16]
            y = metabolic_indicators_ogtt_timeseries[subphenotype] 
        else:
            print("Wrong Featureset")
            exit(0)

        X = np.array(X)
        y = np.array(y)
        for j in range(num_iterations):
            skf = StratifiedKFold(n_splits=num_folds, shuffle = True, random_state= None)
            fold=0
            print("X.shape", X.shape)
            print("y.shape", y.shape)
            
            for train_index, test_index in skf.split(X, y):
                fold += 1
                if features == "ogtt":
                    x_train_tmp, x_test_tmp = X[train_index], X[test_index]
                    ogtt_pca = sk.decomposition.PCA()
                    principalComponents = ogtt_pca.fit_transform(x_train_tmp)
                    principalDf = pd.DataFrame(data = principalComponents)
                    principalDf
                    principalDf_2pcs = principalDf.loc[:,[0,1]]
                    principalDf_2pcs.columns = ["PC1", "PC2"]
                    x_train = principalDf_2pcs
                    test_principalComponents = ogtt_pca.transform(x_test_tmp)
                    test_principalDf = pd.DataFrame(data = test_principalComponents)
                    test_principalDf_2pcs = test_principalDf.loc[:,[0,1]]
                    test_principalDf_2pcs.columns = ["PC1", "PC2"]
                    x_test = test_principalDf_2pcs
                    y_train, y_test = y[train_index], y[test_index]
                else:
                    x_train, x_test = X[train_index], X[test_index]
                    y_train, y_test = y[train_index], y[test_index]

                for index, (name, classifier) in enumerate(classifiers.items()):
                    print("### Features=", features, "   Subphenotye=", subphenotype, "   Iteration = ", j, "   Fold = ", fold, "   classifier=", name)
                    classifier.fit(x_train, y_train)
                    y_pred_train = classifier.predict(x_train)
                    y_pred_test = classifier.predict(x_test)


                    ### Accuracy, AUC, Recall, Precision, F1
                    accuracy = metrics.accuracy_score(y_test, y_pred_test)
                    f1_macro = metrics.f1_score(y_test, y_pred_test, average = "macro")
                    f1_micro = metrics.f1_score(y_test, y_pred_test, average = "micro")
                    prec_macro = metrics.precision_score(y_test, y_pred_test, average = "macro")
                    prec_micro = metrics.precision_score(y_test, y_pred_test, average = "micro")
                    recall_macro = metrics.recall_score(y_test, y_pred_test, average = "macro")
                    recall_micro = metrics.recall_score(y_test, y_pred_test, average = "micro")

                    ## Confusion matrix
                    cm = metrics.confusion_matrix(y_test, y_pred_test)

                    # roc curve and auc
                    probs_test = classifier.predict_proba(x_test)
                    if subphenotype == "muscle_ir_status":
                        classifier_auc = metrics.roc_auc_score(y_test, probs_test[:, 1], multi_class='ovr')
                    else:
                        classifier_auc = metrics.roc_auc_score(y_test, probs_test, multi_class='ovr')


                    ## Logging evaluation results
                    f = open("/Users/ahmedm/Documents/" + datetime_stmp + "_iterations_"+ str(num_iterations) + "_ML_Performance.csv", "a")
                    f.write(features + "," + subphenotype + "," + str(j) + "," + str(fold) + "," + name + "," + str(round(classifier_auc,2)) +
                            ","+ str(round(accuracy,2)) +","+ str(round(recall_macro,2)) +","+ 
                            str(round(prec_macro,2)) +","+ str(round(f1_macro,2))+","+ str(round(recall_micro,2)) +
                            ","+ str(round(prec_micro,2)) +","+ str(round(f1_micro,2)) +"\n")
                    f.close()