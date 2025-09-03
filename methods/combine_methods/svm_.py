#!/usr/bin/env Python
# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
import matplotlib
from tqdm import tqdm
tqdm.pandas(ascii=True)
# from sklearn.externals import joblib
# import sklearn.external.joblib as extjoblib
import joblib

from sklearn.model_selection import train_test_split,StratifiedKFold
from sklearn import preprocessing,metrics
from sklearn.pipeline import make_pipeline
from sklearn.metrics import recall_score,accuracy_score,matthews_corrcoef
from sklearn.metrics import make_scorer
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GridSearchCV
import os,sys
import time

from thundersvm import *
from thundersvm import SVC



random_state = np.random.RandomState(0)

def plot_roc(label, predict_score, label1, predict_score1, path,title):
    fpr, tpr, _ = roc_curve(label, predict_score)
    fpr1, tpr1, _ = roc_curve(label1, predict_score1)
    roc_auc = auc(fpr, tpr)
    roc_auc1 = auc(fpr1, tpr1)
    plt.figure(figsize=(10, 10))
    # plt.rcParams['figure.figsize'] = (10, 10)  # ͼ�δ�С
    plt.rcParams['savefig.dpi'] = 300  # ͼƬ����
    plt.rcParams['figure.dpi'] = 300  # �ֱ���
    plt.rcParams.update({'font.size': 15})

    # ����ͼ���ߴ�ϸ
    bwith = 2.0  # �߿��������Ϊ2
    TK = plt.gca()  # ��ȡ�߿�
    TK.spines['bottom'].set_linewidth(bwith)  # ͼ���±�
    TK.spines['left'].set_linewidth(bwith)  # ͼ�����
    TK.spines['top'].set_linewidth(bwith)  # ͼ���ϱ�
    TK.spines['right'].set_linewidth(bwith)  # ͼ���ұ�
    lw = 2
    # plt.figure(figsize=(10,10))
    plt.plot(fpr, tpr, color='darkorange', marker='o', markersize=4, lw=lw, linestyle='dashed',
             label='Train ROC Curve (area = %0.3f)' % roc_auc)
    plt.plot(fpr1, tpr1, color='green', marker='o', lw=lw, markersize=4, linestyle='dashed',
             label='Test ROC Curve (area = %0.3f)' % roc_auc1)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel('False Positive Rate', fontsize=25)
    plt.ylabel('True Positive Rate', fontsize=25)
    # plt.title('Receiver Operating Characteristic', fontsize=20, pad=10)
    plt.title('ROC Curve of Support Vector Machine', fontsize=25, pad=10)
    # plt.legend(loc="lower right")
    plt.legend(frameon=False, loc="lower right", fontsize='large')  # ����ͼ���ޱ߿򣬽�ͼ���������Ͻ�

    plt.savefig(path + title+'_svm_auc.png')
    plt.close()

def train_svm(xp, yp, path, title,kernel_choose = 'linear'):
    random_state = np.random.RandomState(0)
    # xp = np.load(xp)
    xp = xp

    if len(xp.shape) == 2:

        X = xp
        y = np.load(yp).ravel()
        # shuffle and split training and test sets
    if len(xp.shape) == 3:

        X = xp.reshape(xp.shape[0], -1)
        y = np.load(yp).ravel()
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0,stratify = y)

    # Learn to predict each class against the other
    model = svm.SVC(kernel=kernel_choose, probability=True, random_state=random_state).fit(X_train, y_train)
    y_score_train = model.decision_function(X_train)
    y_score = model.decision_function(X_test)
    plot_roc(y_train, y_score_train, y_test, y_score, path,title)


def mean_fun(onelist):
    count = 0
    for i in onelist:
        count += i
    return float(count/len(onelist))

#评估模型的分类效果
def calc_metrics(y_label, y_proba,predict_):
    con_matrix = metrics.confusion_matrix(y_label, predict_)
    # con_matrix = confusion_matrix(y_label, [1 if x >= 0.5 else 0 for x in y_proba])
    TN = con_matrix[0][0]
    FP = con_matrix[0][1]
    FN = con_matrix[1][0]
    TP = con_matrix[1][1]
    P = TP + FN
    N = TN + FP
    Sn = TP / P if P > 0 else 0
    Sp = TN / N if N > 0 else 0

    ACC = metrics.accuracy_score(y_label, predict_)
    MCC = metrics.matthews_corrcoef(y_label, predict_)
    precision = metrics.precision_score(y_label, predict_)
    f1 = metrics.f1_score(y_label, predict_)
    recall = metrics.recall_score(y_label, predict_)
    

    auc_ = roc_auc_score(y_label,y_proba)

    return ACC, Sn, Sp, precision, MCC, f1,auc_,recall

def twoclass_svm(ComResultAB_np,label,path,c_v = 1):

    # #############################################################################
    # Classification and ROC analysis
    # Run classifier with cross-validation and plot ROC curves
    if c_v != 1:
        mcc = -1

        cv = StratifiedKFold(n_splits=c_v,random_state=random_state,shuffle=True)
        # classifier = svm.SVC(decision_function_shape='ovo',random_state=random_state,probability=True)
        # classifier = svm.SVC(kernel='rbf', probability=True, random_state=random_state)
        ACCs = []
        MCCs = []
        precisions = []
        F1s = []
        recalls = []
        AUCs = []
        Sns = []
        Sps = []

        # for epoch in tqdm(range(1), ascii=True):
        for i, (train, test) in enumerate(cv.split(ComResultAB_np, label)):
            time_start = time.time()
            print('ComResultAB_np[train].shape')
            print(ComResultAB_np[train].shape)
            #数据标准化
            scaler = preprocessing.StandardScaler().fit(ComResultAB_np[train])
            X_train_transformed = scaler.transform(ComResultAB_np[train])

            #网格搜索最优参数
            # parameters = {'C':[1, 5, 1000],'gamma':[0.1, 0.01]}
            
            parameters = {'C':[1, 10, 100,1024], 'gamma':[0.125, 0.5 ,1, 4]}
            # parameters = {'C':[1024], 'gamma':[0.5]}
            svc_greid = svc = SVC(kernel='rbf', probability=True, random_state=random_state)
            # svc = svm.SVC(kernel='rbf', probability=True, random_state=random_state)
            classifier_grid = GridSearchCV(svc_greid, parameters)
            classifier_grid.fit(X_train_transformed, label[train])
            # clf = GridSearchCV(svc, parameters, cv=5, n_jobs=1)
            best_para = classifier_grid.best_params_
            print(best_para)
            svc = SVC(kernel='rbf', C = best_para['C'], gamma = best_para['gamma'] , probability=True, random_state=random_state)
            classifier = svc
            #模型训练
            classifier.fit(X_train_transformed, label[train])
            #model validation
            X_test_transformed = scaler.transform(ComResultAB_np[test])
            predict_ = classifier.predict(X_test_transformed)
            predict_score = classifier.decision_function(X_test_transformed)
            y_prob = classifier.predict_proba(X_test_transformed)


            # Modeal Evaluation
            print('label and predict')
            print(label[test])
            print(predict_)

            ACC, Sn, Sp, precision, MCC, f1,auc_,recall = calc_metrics(label[test], predict_score,predict_)

            # auc_ = metrics.roc_auc_score(label[test], y_prob, multi_class="ovo",
            #                                      average="weighted")
            ACCs.append(ACC)
            MCCs.append(MCC)
            precisions.append(precision)
            F1s.append(f1)
            recalls.append(recall)
            AUCs.append(auc_)
            Sns.append(Sn)
            Sps.append(Sp)

            time_end = time.time()
            time_sum = time_end - time_start
            print('The %s fold cost %s minutes ' % (i, time_sum/60))

            if mcc < MCC or mcc == MCC:
                # joblib.dump(classifier, path)
                # print(classifier.best_params_)
                print(path)
                classifier.save_to_file(path)
                print('save the model ending')
                mcc = MCC



        MCC_mean = mean_fun(MCCs)
        print(MCC_mean)



        

        mcc = MCC_mean
        print(mcc)
        ACC_mean = mean_fun(ACCs)
        precision_mean = mean_fun(precisions)
        F1s_mean = mean_fun(F1s)
        recalls_mean = mean_fun(recalls)
        AUCs_mean = mean_fun(AUCs)
        Sns_mean = mean_fun(Sns)
        Sps_mean = mean_fun(Sps)
    else:
        time_start = time.time()
        print('traing data shape is {shape}'.format(shape = ComResultAB_np.shape))
        #数据标准化
        scaler = preprocessing.StandardScaler().fit(ComResultAB_np)
        X_train_transformed = scaler.transform(ComResultAB_np)
        # 网格搜索
        # parameters = {'C':[1, 10, 100], 'gamma':[0.01, 0.125, 1 ]}
        parameters = {'C':[1, 10, 100], 'gamma':[0.001, 0.01,0.125]} 
        svc_greid = SVC(kernel='rbf', probability=True, random_state=random_state)
        classifier_grid = GridSearchCV(svc_greid, parameters, n_jobs=10)
        classifier_grid.fit(X_train_transformed, label)
        best_para = classifier_grid.best_params_
        print('The best traing parameters is {para}'.format(para = best_para))
        #模型训练
        thsvc = SVC(kernel='rbf', C = best_para['C'], gamma = best_para['gamma'] , probability=True, random_state=random_state)
        thsvc.fit(X_train_transformed, label)
        #save the model
        thsvc.save_to_file(path)

        mcc,ACC_mean,precision_mean,F1s_mean,recalls_mean,AUCs_mean,Sns_mean,Sps_mean = [0 for x in range(8)]
        print('The traing svm model on %s human sample cost %s minutes ' % (ComResultAB_np.shape[0],(time.time() - time_start)/60))

    return mcc,ACC_mean,precision_mean,F1s_mean,recalls_mean,AUCs_mean,Sns_mean,Sps_mean,best_para

def two_svm_predict(modelpath,X_test,y_test):
    import joblib
    print('predict starting')
    clf = SVC()
    clf.load_from_file(modelpath)
    # clf.set_params(**tfidf_params)
    
    scaler = preprocessing.StandardScaler().fit(X_test)
    X_test_transformed = scaler.transform(X_test)
    predict_ = clf.predict(X_test_transformed)
    predict_score = clf.decision_function(X_test_transformed)
    y_prob = clf.predict_proba(X_test_transformed)
    print('prediction ending')
    # Modeal Evaluation
    ACC, Sn, Sp, precision, MCC, f1,auc_,recall = calc_metrics(y_test, predict_score,predict_)

    return ACC, Sn, Sp, precision, MCC, f1,auc_,recall
    
def multiclass_svm(ComResultAB_np,label,path,Epochs = 1):

    # #############################################################################
    # Classification and ROC analysis
    # Run classifier with cross-validation and plot ROC curves
    mcc = -1

    for epoch in tqdm(range(Epochs), ascii=True):
        random_state = np.random.RandomState(epoch)
        cv = StratifiedKFold(n_splits=5,random_state=random_state)
        classifier = svm.SVC(decision_function_shape='ovo',random_state=random_state,probability=True)

        ACCs = []
        MCCs = []
        precisions = []
        F1s = []
        recalls = []
        AUCs = []



        for i, (train, test) in enumerate(cv.split(ComResultAB_np, label)):
            scaler = preprocessing.StandardScaler().fit(ComResultAB_np[train])
            X_train_transformed = scaler.transform(ComResultAB_np[train])

            classifier.fit(X_train_transformed, label[train])
            X_test_transformed = scaler.transform(ComResultAB_np[test])
            predict_ = classifier.predict(X_test_transformed)
            predict_score = classifier.decision_function(X_test_transformed)
            y_prob = classifier.predict_proba(X_test_transformed)

            # Modeal Evaluation
            print('label and predict')
            print(label[test])
            print(predict_)
            ACC = metrics.accuracy_score(label[test], predict_)
            MCC = metrics.matthews_corrcoef(label[test], predict_)
            precision = metrics.precision_score(label[test], predict_, average='macro')
            f1 = metrics.f1_score(label[test], predict_, average='weighted')
            recall = metrics.recall_score(label[test], predict_, average='micro')
            auc_ = metrics.roc_auc_score(label[test], y_prob, multi_class="ovo",
                                                 average="weighted")
            ACCs.append(ACC)
            MCCs.append(MCC)
            precisions.append(precision)
            F1s.append(f1)
            recalls.append(recall)
            AUCs.append(auc_)



        MCC_mean = mean_fun(MCCs)
        print(MCC_mean)



        if mcc < MCC_mean or mcc == MCC_mean:
            # joblib.dump(classifier, path)
            classifier.save_to_file(path)

            mcc = MCC_mean
            print(mcc)
            ACC_mean = mean_fun(ACCs)
            precision_mean = mean_fun(precisions)
            F1s_mean = mean_fun(F1s)
            recalls_mean = mean_fun(recalls)
            AUCs_mean = mean_fun(AUCs)

    return mcc,ACC_mean,precision_mean,F1s_mean,recalls_mean,AUCs_mean


def mul_svm_predict(modelpath,X_test,y_test):
    import joblib

    clf = joblib.load(modelpath)
    scaler = preprocessing.StandardScaler().fit(X_test)
    X_test_transformed = scaler.transform(X_test)
    predict_ = clf.predict(X_test_transformed)
    y_prob = clf.predict_proba(X_test_transformed)

    # Modeal Evaluation

    ACC = metrics.accuracy_score(y_test, predict_)
    MCC = metrics.matthews_corrcoef(y_test, predict_)
    precision = metrics.precision_score(y_test, predict_, average='macro')
    f1 = metrics.f1_score(y_test, predict_, average='weighted')
    recall = metrics.recall_score(y_test, predict_, average='micro')
    fbeta = metrics.fbeta_score(y_test, predict_, average='macro', beta=0.5)
    auc_ = metrics.roc_auc_score(y_test, y_prob, multi_class="ovo",
                                 average="weighted")

    return ACC,MCC,precision,f1,recall,auc_,fbeta