#!/usr/bin/env Python
# coding=utf-8

from argparse import ArgumentParser
from collections import Counter
from functools import reduce
from mrmr import mrmr_classif
from sklearn import metrics
from sklearn import metrics,preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import auc as calculate_auc
from sklearn.metrics import f1_score
from sklearn.metrics import mean_squared_error
from sklearn.metrics import roc_auc_score, precision_recall_curve
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.model_selection import train_test_split, PredefinedSplit
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from tensorflow.keras import Input
from tensorflow.keras.layers import Dense, Dropout,BatchNormalization,concatenate,Conv1D,Activation,MaxPooling1D,Flatten
from tensorflow.keras.models import Sequential,Model
from tensorflow.keras.utils import to_categorical
from tqdm import tqdm
tqdm.pandas(ascii=True)

import glob
import itertools
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import random
import tensorflow as tf
import time
import xgboost as xgb

def make_comcod(com_n = None):
    combination_methods =  [['Open reading frame (1D)', 'Entropy density of transcript (1D)', 'Global descriptor (1D)', 'K-mer (1D)', 'Codon related (1D)', 'Pseudo protein related (1D)', 'Guanine-cytosine related (1D)', 'Nucleotide related (1D)', 'Secondary structure (1D)', 'EIIP based spectrum (1D)','Solubility lipoaffinity (1D)','Partition coefficient (1D)','Polarizability refractivity (1D)','Hydrogen bond related (1D)','Topological indice (1D)','Molecular fingerprint (1D)']]

    return combination_methods
    
def standardization(X):
    scaler = StandardScaler()
    scaler.fit(X)
    X = scaler.transform(X)
    return X, scaler

def mknpy_RNAonly(combin2,datapath):
    # generate the npy data according to the combination lists
    for numc in range(len(combin2)):
        # print(combin2[numc])
        file_tempA = np.load(os.path.join(datapath, combin2[numc] + '.npy'))
        if numc == 0:
            combnpysA = file_tempA
        else:
            combnpysA = np.concatenate((combnpysA, file_tempA), axis=1)
    # print('Feature A shape is {shapeA}'.format(shapeA = combnpysA.shape))
    return combnpysA

def mknpy_RNA_RNA(combin2,datapath):
    # generate the npy data according to the combination lists
    for numc in range(len(combin2)):
#         print(combin2[numc])
        file_tempA = np.load(os.path.join(datapath, combin2[numc] + '_A.npy'))
        file_tempB = np.load(os.path.join(datapath, combin2[numc] + '_B.npy'))
        if numc == 0:
            combnpysA = file_tempA
            combnpysB = file_tempB
        else:
            combnpysA = np.concatenate((combnpysA, file_tempA), axis=1)
            combnpysB = np.concatenate((combnpysB, file_tempB), axis=1)
    # print('Feature A and B file shape is {shapeA},{shapeB}'.format(shapeA = combnpysA.shape,shapeB = combnpysB.shape))
    # combineFeat = np.concatenate((combnpysA,combnpysB),axis = 1)
    return combnpysA,combnpysB

def mknpy_RNA_pro(combin2,datapath):
    # generate the npy data according to the combination lists
    combnpysB = np.load(os.path.join(datapath, 'protein_B.npy'))
    for numc in range(len(combin2)):
        file_tempA = np.load(os.path.join(datapath, combin2[numc] + '_A.npy'))        
        if numc == 0:
            combnpysA = file_tempA
        else:
            combnpysA = np.concatenate((combnpysA, file_tempA), axis=1)
    # combineFeat = np.concatenate((combnpysA,combnpysB),axis = 1)
    return combnpysA,combnpysB

def mknpy_RNA_compound(combin2,datapath):
    # generate the npy data according to the combination lists
    combnpysB = np.load(os.path.join(datapath, 'compound_B.npy'))
    for numc in range(len(combin2)):
        file_tempA = np.load(os.path.join(datapath, combin2[numc] + '_A.npy'))        
        if numc == 0:
            combnpysA = file_tempA
        else:
            combnpysA = np.concatenate((combnpysA, file_tempA), axis=1)
    print('combnpysA')
    print(combnpysA.shape)
    print('combnpysB')
    print(combnpysB.shape)
    # combineFeat = np.concatenate((combnpysA,combnpysB),axis = 1)
    return combnpysA,combnpysB

def svm_two(train_x, valid_x, train_y, valid_y):
    # 合并训练集和验证集
    train_val_features = np.concatenate((train_x,valid_x),axis=0)
    train_val_labels = np.concatenate((train_y,valid_y),axis=0)
    test_fold = np.zeros(train_val_features.shape[0])   # 将所有index初始化为0,0表示第一轮的验证集
    test_fold[:train_x.shape[0]] = -1            # 将训练集对应的index设为-1，表示永远不划分到验证集中
    ps = PredefinedSplit(test_fold=test_fold)
    
    random_state = 43
    # parameters = {'C':[0.1, 1, 10, 100,1024], 'gamma':[0.001, 0.01,0.125, 0.5 ,1, 4]}
    # 网格搜索最优参数    
    parameters = {'C':[10], 'gamma':[0.001]}    
    svc_greid = SVC(kernel='rbf', probability=True, random_state=random_state)  
    
    grid_search_params = {'estimator': svc_greid,             # 目标分类器
                      'param_grid': parameters,  # 前面定义的我们想要优化的参数
                      'cv': ps,                     # 使用前面自定义的split验证策略
                      'n_jobs': 20}
    classifier_grid = GridSearchCV(**grid_search_params)
    classifier_grid.fit(train_val_features, train_val_labels)

    best_para = classifier_grid.best_params_
    # print('SVM the best parameters are {}'.format(best_para))
    best_estimator = classifier_grid.best_estimator_

    # print("thundersvm Training 20 group Time: %s minutes" % ((time.time() - time0)/60))
    return best_estimator,best_para

def RF(train_x, valid_x, train_y, valid_y):
    
    # 合并训练集和验证集
    train_val_features = np.concatenate((train_x,valid_x),axis=0)
    train_val_labels = np.concatenate((train_y,valid_y),axis=0)
    test_fold = np.zeros(train_val_features.shape[0])   # 将所有index初始化为0,0表示第一轮的验证集
    test_fold[:train_x.shape[0]] = -1            # 将训练集对应的index设为-1，表示永远不划分到验证集中
    ps = PredefinedSplit(test_fold=test_fold)
    
    random_state = 43
    # 网格搜索最优参数
    rfc = RandomForestClassifier(random_state=random_state)    
    tuned_parameters = [{'n_estimators':[70,170]}] #[50,70,100,130,170,200]
    grid_search_params = {'estimator': rfc,             # 目标分类器
                      'param_grid': tuned_parameters,  # 前面定义的我们想要优化的参数
                      'cv': ps,                     # 使用前面自定义的split验证策略
                      'n_jobs': 10}                  # 并行运行的任务数，-1表示使用所有               # 输出信息，数字越大输出信息越多
    clf = GridSearchCV(**grid_search_params)
    clf.fit(train_val_features, train_val_labels)
    best_para = clf.best_params_    
    # print('RF the best parameters are {}'.format(best_para))
    # 根据最优参数得到最好训练模型
    best_estimator = clf.best_estimator_  
    # rfc = RandomForestClassifier(n_estimators = best_para['n_estimators'], random_state=random_state)
    return best_estimator,best_para

def xgboost_model(traindata,trainlabel,X_test,y_test):
    # Convert input data from numpy to XGBoost format
    random_state = 43   

    dtrain = xgb.DMatrix(traindata, label=trainlabel)
    dtest = xgb.DMatrix(X_test, label=y_test)
    
    # Specify sufficient boosting iterations to reach a minimum
    num_rounds = [100]#[50, 100, 200, 300, 500]
    max_depths = [3, 6]#[3, 6, 10]
    learning_rates = [0.05]#[0.05, 0.1, 0.15]

    parameters_all = [[x,y,z] for x in num_rounds for y in max_depths for z in learning_rates] 
    parameters = random.sample(parameters_all, 2) 
    res_dict = {}

    for index, parameter in enumerate(parameters):
        # print('Current parameter is {}/{} dropout, learning rate, batchsize: {}'.format(index,len(parameters),parameter))
        num_round = parameter[0]
        max_depth = parameter[1]
        learning_rate = parameter[2]   

        # Leave most parameters as default
        param = {
            'booster': 'gbtree',
            'objective': 'binary:logistic',  # 多分类的问题multi:softmax
            # 'num_class': 2,               # 类别数，与 multisoftmax 并用
            'gamma': 0.1,                  # 用于控制是否后剪枝的参数,越大越保守，一般0.1、0.2这样子。
            'max_depth': max_depth,               # 构建树的深度，越大越容易过拟合
            'lambda': 2,                   # 控制模型复杂度的权重值的L2正则化项参数，参数越大，模型越不容易过拟合。
            'subsample': 1,              # 随机采样训练样本
            'colsample_bytree': 0.7,       # 生成树时进行的列采样
            'min_child_weight': 3,
            'verbosity': 0,                   # 设置成1则没有运行信息输出，最好是设置为0.
            'eta': learning_rate,                  # 如同学习率
            'seed': 1000,
            'nthread': 4,                  # cpu 线程数
            'tree_method': 'gpu_hist' # Use GPU accelerated algorithm
        }

        # Repeat for CPU algorithm
        tmp = time.time()
        param['tree_method'] = 'hist'
        cpu_res = {}
        bst  = xgb.train(param, dtrain, num_round,  evals_result=cpu_res)
        # print("CPU Training Time: %s seconds" % (str(time.time() - tmp)))
        
        score_valid = bst.predict(dtest)
        group_valid = (score_valid >= 0.5)*1
        Acc_valid, Sn_valid, Sp_valid, Pre_valid, MCC_valid, AUC_valid = calc_metrics(y_test,score_valid,group_valid)

        score = bst.predict(dtest)
        # print(ans)
        group = (score >= 0.5)*1
        Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(y_test,score,group)

        res_dict[index] = [parameter, Acc, Sn, Sp, Pre, MCC, AUC, MCC_valid]
        # print('The parameter is dropout, learning rate, batchsize: {}, and test MCC: {}'.format(parameter,MCC))
        # print('Each group parameters selection use time: {} min'.format((time.time() - tmp)/60))
    df = pd.DataFrame(res_dict)
    df_list = df.iloc[7,:].tolist()
    maxindex = df_list.index(max(df_list))
    dfres = df.iloc[:,maxindex].tolist()
    # print('The best parameter is dropout, learning rate, batchsize: {}'.format(dfres[0]))
    best_para = dfres[0]
    Acc, Sn, Sp, Pre, MCC, AUC = dfres[1],dfres[2],dfres[3],dfres[4],dfres[5],dfres[6]
   
    return Acc, Sn, Sp, Pre, MCC, AUC,best_para


#evaluation the performence
def calc_metrics(y_label, y_proba,y_predict):

    con_matrix = metrics.confusion_matrix(y_label, y_predict)
    # print(con_matrix)
    TN = con_matrix[0][0]
    FP = con_matrix[0][1]
    FN = con_matrix[1][0]
    TP = con_matrix[1][1]
    P = TP + FN
    N = TN + FP
    Sn = TP / P if P > 0 else 0
    Sp = TN / N if N > 0 else 0
    Acc = (TP + TN) / (P + N) if (P + N) > 0 else 0
    Pre = (TP) / (TP + FP) if (TP+FP) > 0 else 0
    MCC = 0
    tmp = math.sqrt((TP + FP) * (TP + FN)) * math.sqrt((TN + FP) * (TN + FN))
    if tmp != 0:
        MCC = (TP * TN - FP * FN) / tmp
    fpr, tpr, thresholds = metrics.roc_curve(y_label, y_proba)
    AUC = metrics.auc(fpr, tpr)
    return Acc, Sn, Sp, Pre, MCC, AUC

class DNNModel(tf.keras.models.Model):
    def __init__(self, classes, shape1, shape2, dropout):
        super(DNNModel, self).__init__()

        seed = 1234

        self.d1 = tf.keras.layers.Dense(128, activation='relu')
        # self.normal1 = tf.keras.layers.BatchNormalization()
        self.drop1 = tf.keras.layers.Dropout(dropout, seed = seed)
        self.d2 = tf.keras.layers.Dense(112, activation='relu')
        self.drop8 = tf.keras.layers.Dropout(dropout, seed = seed)
        self.d9 = tf.keras.layers.Dense(classes, activation='softmax')

    def call(self, inputs,training=None):
        # print('inputs.shape')
        # print(inputs.shape)

        x1 = self.d1(inputs)
        x1 = self.drop1(x1,training = training)
        x1 = self.d2(x1)
        x2 = self.drop8(x1,training = training)
        
        out = self.d9(x2)


        return out, x2

class CNN1DModel(tf.keras.models.Model):
    def __init__(self, classes, shape1, shape2, dropout):
        super(CNN1DModel, self).__init__()

        self.conv1 = tf.keras.layers.Conv1D(
            filters=32,  # ��������Ԫ�������ˣ���Ŀ
            kernel_size=3,  # ����Ұ��С
            padding="same",  # padding���ԣ�vaild �� same��
            activation=tf.nn.relu  # �����
        )
        self.pool1 = tf.keras.layers.MaxPool1D(pool_size=2, strides=2, padding='valid')
        self.normal1 = tf.keras.layers.BatchNormalization()
        self.drop1 = tf.keras.layers.Dropout(dropout)

        self.conv2 = tf.keras.layers.Conv1D(
            filters=32,
            kernel_size=3,
            padding='same',
            activation=tf.nn.relu
        )
        self.pool2 = tf.keras.layers.MaxPool1D(pool_size=2, strides=2, padding='valid')
        self.normal2 = tf.keras.layers.BatchNormalization()
        self.drop2 = tf.keras.layers.Dropout(dropout)

        self.conv3 = tf.keras.layers.Conv1D(
            filters=64,
            kernel_size=3,
            padding='same',
            activation=tf.nn.relu
        )
        self.pool3 = tf.keras.layers.MaxPool1D(pool_size=2, strides=2, padding='valid')
        self.normal3 = tf.keras.layers.BatchNormalization()
        self.drop3 = tf.keras.layers.Dropout(dropout)

        # self.flatten = tf.keras.layers.Reshape(target_shape=(6 * 1 * 64,)) # 3mers
        #         self.flatten = tf.keras.layers.Reshape(target_shape=(30 * 1 * 64,))# 4mers
        self.flatten = tf.keras.layers.Reshape(target_shape=(shape1 * 1 * 64,))  # Transcript
        # self.dense1 = tf.keras.layers.Dense(units=shape1, activation=tf.nn.relu)
        # self.drop4 = tf.keras.layers.Dropout(dropout)
        self.dense1 = tf.keras.layers.Dense(units=shape2, activation=tf.nn.relu)
        self.dense2 = tf.keras.layers.Dense(classes, activation='softmax')

    def call(self, inputs, training):
        inputs = tf.expand_dims(inputs, axis=2)
        x = self.conv1(inputs)
        x = self.pool1(x)
        x = self.normal1(x, training=training)
        x = self.drop1(x)  # (15, 58, 32)

        x = self.conv2(x)
        x = self.pool2(x)
        x = self.normal2(x, training=training)
        x = self.drop2(x)  # (15, 28, 32)

        x = self.conv3(x)
        x = self.pool3(x)
        x = self.normal3(x, training=training)
        x = self.drop3(x)  # (15, 13, 64)
        # print('The shape before reshape')
        # print(x.shape)
        x = self.flatten(x)
        x2 = self.dense1(x)
        output = self.dense2(x2)

        return output, x2

class Dataset_make(tf.keras.utils.Sequence):

    def __init__(self, x_set, y_set, batch_size):
        self.x, self.y = x_set, y_set
        self.batch_size = batch_size

    def __len__(self):
        return math.ceil(len(self.x) / self.batch_size)

    def __getitem__(self, idx):
        batch_x = self.x[idx * self.batch_size:(idx + 1) *
                                               self.batch_size]
        batch_y = self.y[idx * self.batch_size:(idx + 1) *
                                               self.batch_size]
        c = Counter(list(batch_y))
        # print(dict(c))

        return batch_x, batch_y

def train_model(X_train,y_train,X_test,y_test,modelnm = 'DNN'):
    
    # parameters
    classes = len(np.unique(y_train))
    shape1 = 128
    shape2 = 32
    Epochs = 50

    # dropouts = [0.2,0.5,0.7]
    # learnrantes = [0.001,0.01,0.1]
    # batch_sizes = [32,64,128,256]

    # 超参数选择
    dropouts = [0.1,0.2]
    learnrantes = [0.001]
    batch_sizes = [64]
    parameters_all = [[x,y,z] for x in dropouts for y in learnrantes for z in batch_sizes] 
    parameters = random.sample(parameters_all, 2) 
    
    #traing, validation 8:2 分层抽样
    train_x, test_x, train_y, test_y = train_test_split(X_train, y_train, test_size=0.2, random_state=42,stratify=y_train)
    
    res_dict = {}
    for index, parameter in enumerate(parameters):
        # print('Current parameter is {}/{} dropout, learning rate, batchsize: {}'.format(index,len(parameters),parameter))
        dropout = parameter[0]
        learnrante = parameter[1]
        batch_size = parameter[2]

        time0 = time.time()        
    
        # make the model
        if modelnm == 'DNN':
            # DNN
            model = DNNModel(classes, shape1, shape2, dropout)
        elif modelnm == 'CNN':
            # CNN
            shape1 = int((int((int((train_x.shape[1]) / 2)) / 2)) / 2)
            model = CNN1DModel(classes, shape1, shape2, dropout)

        # define loss and optimizer
        loss_obj = tf.keras.losses.SparseCategoricalCrossentropy()
        optimizer = tf.keras.optimizers.Adam(lr=learnrante)

        # training model
        def train_step(data, labels):
            with tf.GradientTape() as tape:
                N = len(labels)
                logits, x2 = model(data, training=True)
                

                loss = loss_obj(labels, logits)

                loss = tf.reduce_mean(loss)
                loss = loss / N

                # print(loss)
            grads = tape.gradient(loss, model.trainable_variables)
            optimizer.apply_gradients(zip(grads, model.trainable_variables))

            group_t = np.argmax(logits.numpy(), axis=1)
            # train_acc(labels, logits)
            return loss,logits,group_t
        
        # testing model
        def test_step(data, labels):
            N = len(labels)
            logits, x2 = model(data, training=False)
            # logits,x2 = model.predict(data, batch_size=32, verbose=1)
            test_loss = loss_obj(labels, logits)
            test_loss = tf.reduce_mean(test_loss)
            test_loss = test_loss / N
            logits_ = np.argmax(logits.numpy(), axis=1)
            # test_a = test_acc(labels, logits)
            return test_loss, logits_, logits
        
        # data preparing
        train_dataset = Dataset_make(train_x, train_y, batch_size)
        test_dataset = Dataset_make(test_x, test_y, batch_size)

        mcc = -1
        for epoch in range(Epochs):
            # train_acc.reset_states()  
            # test_acc.reset_states()

            group_all = np.array([], dtype=int)
            labels_all = np.array([], dtype=int)
            score_all = np.zeros(shape=(1,2))

            group_train = np.array([], dtype=int)
            labels_train = np.array([], dtype=int)
            score_trains = np.zeros(shape=(1,2))

            for images, labels in train_dataset:
                t_loss,score_train,group_t = train_step(images, labels)
                group_train = np.append(group_train, group_t)
                labels_train = np.append(labels_train, labels)
                score_trains = np.concatenate((score_trains, score_train),axis = 0)


            for images, labels in test_dataset:
                loss, group, score = test_step(images, labels)
                # print('images')
                # print(images.shape)
                group_all = np.append(group_all, group)
                labels_all = np.append(labels_all, labels)

                score_all = np.concatenate((score_all, score),axis = 0)
            score_all = np.delete(score_all, 0, 0)
            score_trains = np.delete(score_trains, 0, 0)

            Acc_train, Sn_train, Sp_train, Pre_train, MCC_train, AUC_train = calc_metrics(labels_train,score_trains[:,1],group_train)

            MCC_v = metrics.matthews_corrcoef(labels_all, group_all)
            score, feature_out = model(X_test, training=False)
            # score = model.predict(X_test, batch_size=32, verbose=1)
            group = np.argmax(score.numpy(), axis=1)
            # Acc_e, Sn_e, Sp_e, Pre_e, MCC_e, AUC_e = calc_metrics(y_test,score,group)

            if mcc < MCC_v or mcc == MCC_v:
                mcc = MCC_v
                # Modeal Evaluation for epoches
                # print('score_all')
                # print(score_all)
                
                Acc_v, Sn_v, Sp_v, Pre_v, MCC_v, AUC_v = calc_metrics(labels_all,score_all[:,1],group_all)

                last_improve = epoch
                # print('last_improve: %s' % last_improve)           
                Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(y_test,score[:,1],group)



            if epoch - last_improve >= 10:
                break

            tmp = 'Epoch {:.3f}, traing Acc {:.3f}, validation Acc {:.3f}, Test Acc {:.3f}, train_loss{:.3f}, validation loss{:.3f}, test auc{:.3f}, test mcc{:.3f}'
            # print(tmp.format(epoch + 1, Acc_train,
            #                 Acc_v,
            #                 Acc, t_loss, loss, AUC,MCC))
            # print(report)      
        # print(model.summary())

        res_dict[index] = [parameter, Acc, Sn, Sp, Pre, MCC, AUC]
        # print('The parameter is dropout, learning rate, batchsize: {}, and test MCC: {}'.format(parameter,MCC))
        # print('Each group parameters selection use time: {} min'.format((time.time() - time0)/60))
    df = pd.DataFrame(res_dict)
    df_list = df.iloc[5,:].tolist()
    maxindex = df_list.index(max(df_list))
    dfres = df.iloc[:,maxindex].tolist()
    # print('The best parameter is dropout, learning rate, batchsize: {}'.format(dfres[0]))
    Acc, Sn, Sp, Pre, MCC, AUC = dfres[1],dfres[2],dfres[3],dfres[4],dfres[5],dfres[6]
    return Acc, Sn, Sp, Pre, MCC, AUC,dfres[0]

def rna_encoder_1(data_train_r_npy,data_vaild_r_npy,lr_r,epoch_r,batch_size_r,save_path):
    # lr = 1e-4
    # epoch = 100
    # batch_size = 64

    input_shape = data_train_r_npy.shape[1]
    x_input_r = Input(input_shape)
    x_input01 = tf.expand_dims(x_input_r, axis=2)
    x_input01 = tf.expand_dims(x_input01, axis=3)
    print('x_input: {}'.format(x_input01.shape))

    # ENCODER
    conv_1 = tf.keras.layers.Conv2D(filters=45,  kernel_size=(3, 1), padding="valid",  activation=tf.nn.relu)(x_input01)
    batch_norm_1 = tf.keras.layers.BatchNormalization()(conv_1)
    pooling_1 = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=2, padding='valid')(batch_norm_1)
    print('pooling_1: {}'.format(pooling_1.shape))

    conv_2 = tf.keras.layers.Conv2D(filters=32,  kernel_size=(3, 1), padding="valid",  activation=tf.nn.relu)(pooling_1)
    # print(conv_2.shape)
    batch_norm_2 = tf.keras.layers.BatchNormalization()(conv_2)
    pooling_2 = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=4, padding='valid')(batch_norm_2)
    print('pooling_2: {}'.format(pooling_2.shape))

    conv_3 = tf.keras.layers.Conv2D(filters=16,  kernel_size=(3, 1), padding="valid",  activation=tf.nn.relu)(pooling_2)
    # print(conv_3.shape)
    batch_norm_3 = tf.keras.layers.BatchNormalization()(conv_3)
    pooling_3 = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=2, padding='valid')(batch_norm_3)
    print('pooling_3: {}'.format(pooling_3.shape))

    Flatten1 = tf.keras.layers.Flatten()(pooling_3)
    print('Flatten1: {}'.format(Flatten1.shape))

    Dense_encoder1 = tf.keras.layers.Dense(units=100, activation=tf.nn.relu)(Flatten1)
    print('Dense_encoder1: {}'.format(Dense_encoder1.shape))

    encode_r = Dense_encoder1
    print('encode_r: {}'.format(encode_r.shape))

    Dense_decoder1 = tf.keras.layers.Dense(units=Flatten1.shape[1], activation=tf.nn.relu)(encode_r)
    print('Dense_decoder1: {}'.format(Dense_decoder1.shape))
    print('Dense_decoder1 dtype: {}'.format(Dense_decoder1.dtype))
    
    Dense_decoder1_npy_reshape = tf.keras.layers.Reshape((41, 1, 16))(Dense_decoder1)
    print('Dense_decoder1_npy_reshape: {}'.format(Dense_decoder1_npy_reshape.shape))

    # DECODER
    conv_trans_1 = tf.keras.layers.Conv2DTranspose(16, kernel_size=(3, 1), strides=2, activation=tf.nn.relu)(Dense_decoder1_npy_reshape)  # =>(26, 26, 128)
    batch_norm_trans_1 = tf.keras.layers.BatchNormalization()(conv_trans_1)
    print('batch_norm_trans_1: {}'.format(batch_norm_trans_1.shape))

    conv_trans_2 = tf.keras.layers.Conv2DTranspose(32, kernel_size=(3, 1), strides=2, activation=tf.nn.relu)(batch_norm_trans_1)  # =>(26, 26, 128)
    batch_norm_trans_2 = tf.keras.layers.BatchNormalization()(conv_trans_2)
    print('batch_norm_trans_2: {}'.format(batch_norm_trans_2.shape))

    conv_trans_3 = tf.keras.layers.Conv2DTranspose(45, kernel_size=(3, 1),strides=1,activation=tf.nn.relu)(batch_norm_trans_2)  # =>(26, 26, 128)
    batch_norm_trans_3 = tf.keras.layers.BatchNormalization()(conv_trans_3)
    print('batch_norm_trans_3: {}'.format(batch_norm_trans_3.shape))

    output00_1 = tf.keras.layers.Conv2DTranspose(1, kernel_size=(2, 1), padding='same',strides=1, activation=tf.nn.relu)(batch_norm_trans_3)  # =>(26, 26, 128)
    print('output00_1: {}'.format(output00_1.shape))
    output00_1_1 = tf.keras.layers.Reshape((output00_1.shape[1]*output00_1.shape[2], 1, 1))(output00_1)
    print('output00_1_1: {}'.format(output00_1_1.shape))
    output01 = tf.keras.layers.Conv2DTranspose(1, kernel_size=(2, 1), strides=1, activation='sigmoid')(output00_1_1)
    print('output01: {}'.format(output01.shape))
    # output layer
    output02 = tf.squeeze(output01,  axis=2)
    output = tf.squeeze(output02, axis=2)
    print('output: {}'.format(output.shape))
    autoencoder = tf.keras.Model(inputs=x_input_r, outputs=output)
    # define autoencoder model

    weight_decay = 1e-4
    opt = tf.keras.optimizers.Adam(lr= lr_r, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay= weight_decay) #
    autoencoder.compile(optimizer = opt,loss='mse')

    # fit the autoencoder model to reconstruct input
    history = autoencoder.fit(data_train_r_npy, data_train_r_npy, epochs=epoch_r, batch_size=batch_size_r, verbose=2, validation_data=(data_vaild_r_npy,data_vaild_r_npy))

    encoder_rna = tf.keras.Model(inputs=x_input_r, outputs=encode_r)
    encoder_rna.save( save_path + '/encoder_rna_1.h5')
    return encoder_rna

def rna_encoder_2(data_train_r_npy,data_vaild_r_npy,lr_r,epoch_r,batch_size_r,save_path):
    # lr = 1e-4
    # epoch = 100
    # batch_size = 64

    input_shape = data_train_r_npy.shape[1]
    x_input_r = Input(input_shape)
    x_input01 = tf.expand_dims(x_input_r, axis=2)
    x_input01 = tf.expand_dims(x_input01, axis=3)
    print('x_input: {}'.format(x_input01.shape))

    # ENCODER
    conv_1 = tf.keras.layers.Conv2D(filters=45,  kernel_size=(3, 1), padding="valid",  activation=tf.nn.relu)(x_input01)
    batch_norm_1 = tf.keras.layers.BatchNormalization()(conv_1)
    pooling_1 = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=2, padding='valid')(batch_norm_1)
    print('pooling_1: {}'.format(pooling_1.shape))

    conv_2 = tf.keras.layers.Conv2D(filters=32,  kernel_size=(3, 1), padding="valid",  activation=tf.nn.relu)(pooling_1)
    # print(conv_2.shape)
    batch_norm_2 = tf.keras.layers.BatchNormalization()(conv_2)
    pooling_2 = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=4, padding='valid')(batch_norm_2)
    print('pooling_2: {}'.format(pooling_2.shape))

    conv_3 = tf.keras.layers.Conv2D(filters=16,  kernel_size=(3, 1), padding="valid",  activation=tf.nn.relu)(pooling_2)
    # print(conv_3.shape)
    batch_norm_3 = tf.keras.layers.BatchNormalization()(conv_3)
    pooling_3 = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=2, padding='valid')(batch_norm_3)
    print('pooling_3: {}'.format(pooling_3.shape))

    Flatten1 = tf.keras.layers.Flatten()(pooling_3)
    print('Flatten1: {}'.format(Flatten1.shape))

    Dense_encoder1 = tf.keras.layers.Dense(units=100, activation=tf.nn.relu)(Flatten1)
    print('Dense_encoder1: {}'.format(Dense_encoder1.shape))

    encode_r = Dense_encoder1
    print('encode_r: {}'.format(encode_r.shape))

    Dense_decoder1 = tf.keras.layers.Dense(units=Flatten1.shape[1], activation=tf.nn.relu)(encode_r)
    print('Dense_decoder1: {}'.format(Dense_decoder1.shape))
    print('Dense_decoder1 dtype: {}'.format(Dense_decoder1.dtype))
    
    Dense_decoder1_npy_reshape = tf.keras.layers.Reshape((41, 1, 16))(Dense_decoder1)
    print('Dense_decoder1_npy_reshape: {}'.format(Dense_decoder1_npy_reshape.shape))

    # DECODER
    conv_trans_1 = tf.keras.layers.Conv2DTranspose(16, kernel_size=(3, 1), strides=2, activation=tf.nn.relu)(Dense_decoder1_npy_reshape)  # =>(26, 26, 128)
    batch_norm_trans_1 = tf.keras.layers.BatchNormalization()(conv_trans_1)
    print('batch_norm_trans_1: {}'.format(batch_norm_trans_1.shape))

    conv_trans_2 = tf.keras.layers.Conv2DTranspose(32, kernel_size=(3, 1), strides=2, activation=tf.nn.relu)(batch_norm_trans_1)  # =>(26, 26, 128)
    batch_norm_trans_2 = tf.keras.layers.BatchNormalization()(conv_trans_2)
    print('batch_norm_trans_2: {}'.format(batch_norm_trans_2.shape))

    conv_trans_3 = tf.keras.layers.Conv2DTranspose(45, kernel_size=(3, 1),strides=1,activation=tf.nn.relu)(batch_norm_trans_2)  # =>(26, 26, 128)
    batch_norm_trans_3 = tf.keras.layers.BatchNormalization()(conv_trans_3)
    print('batch_norm_trans_3: {}'.format(batch_norm_trans_3.shape))

    output00_1 = tf.keras.layers.Conv2DTranspose(1, kernel_size=(2, 1), padding='same',strides=1, activation=tf.nn.relu)(batch_norm_trans_3)  # =>(26, 26, 128)
    print('output00_1: {}'.format(output00_1.shape))
    output00_1_1 = tf.keras.layers.Reshape((output00_1.shape[1]*output00_1.shape[2], 1, 1))(output00_1)
    print('output00_1_1: {}'.format(output00_1_1.shape))
    output01 = tf.keras.layers.Conv2DTranspose(1, kernel_size=(2, 1), strides=1, activation='sigmoid')(output00_1_1)
    print('output01: {}'.format(output01.shape))
    # output layer
    output02 = tf.squeeze(output01,  axis=2)
    output = tf.squeeze(output02, axis=2)
    print('output: {}'.format(output.shape))
    autoencoder = tf.keras.Model(inputs=x_input_r, outputs=output)
    # define autoencoder model

    weight_decay = 1e-4
    opt = tf.keras.optimizers.Adam(lr= lr_r, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay= weight_decay) #
    autoencoder.compile(optimizer = opt,loss='mse')

    # fit the autoencoder model to reconstruct input
    history = autoencoder.fit(data_train_r_npy, data_train_r_npy, epochs=epoch_r, batch_size=batch_size_r, verbose=2, validation_data=(data_vaild_r_npy,data_vaild_r_npy))

    encoder_rna = tf.keras.Model(inputs=x_input_r, outputs=encode_r)
    encoder_rna.save( save_path + '/encoder_rna_2.h5')
    return encoder_rna

def pro_encoder(data_train_p_npy,data_vaild_p_npy,lr_p,epoch_p,batch_size_p,save_path):
    # lr_p = 1e-4
    # epoch_p = 100
    # batch_size_p = 64

    input_shape = data_train_p_npy.shape[1]
    x_input_p = Input(input_shape)
    x_input01 = tf.expand_dims(x_input_p, axis=2)
    x_input01 = tf.expand_dims(x_input01, axis=3)
    print('x_input01.shape: {}'.format(x_input01.shape))

    # ENCODER
    conv_1 = tf.keras.layers.Conv2D(filters=45,  kernel_size=(3, 1), padding="valid",  activation=tf.nn.relu)(x_input01)
    batch_norm_1 = tf.keras.layers.BatchNormalization()(conv_1)
    pooling_1 = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=2, padding='valid')(batch_norm_1)
    print('pooling_1: {}'.format(pooling_1.shape))

    conv_2 = tf.keras.layers.Conv2D(filters=32,  kernel_size=(3, 1), padding="valid",  activation=tf.nn.relu)(pooling_1)
    # print(conv_2.shape)
    batch_norm_2 = tf.keras.layers.BatchNormalization()(conv_2)
    pooling_2 = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=4, padding='valid')(batch_norm_2)
    print('pooling_2: {}'.format(pooling_2.shape))

    conv_3 = tf.keras.layers.Conv2D(filters=26,  kernel_size=(3, 1), padding="valid",  activation=tf.nn.relu)(pooling_2)
    # print(conv_3.shape)
    batch_norm_3 = tf.keras.layers.BatchNormalization()(conv_3)
    pooling_3 = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=2, padding='valid')(batch_norm_3)
    print('pooling_3: {}'.format(pooling_3.shape))

    Flatten1 = tf.keras.layers.Flatten()(pooling_3)
    print('Flatten1: {}'.format(Flatten1.shape))

    Dense_encoder1 = tf.keras.layers.Dense(units=100, activation=tf.nn.relu)(Flatten1)
    print('Dense_encoder1: {}'.format(Dense_encoder1.shape))

    encode_r = Dense_encoder1
    print('encode_p: {}'.format(encode_r.shape))

    Dense_decoder1 = tf.keras.layers.Dense(units=Flatten1.shape[1], activation=tf.nn.relu)(encode_r)
    print('Dense_decoder1: {}'.format(Dense_decoder1.shape))
    print('Dense_decoder1 dtype: {}'.format(Dense_decoder1.dtype))
    
    Dense_decoder1_npy_reshape = tf.keras.layers.Reshape((10, 1, 26))(Dense_decoder1)
    print('Dense_decoder1_npy_reshape: {}'.format(Dense_decoder1_npy_reshape.shape))

    conv_trans_2 = tf.keras.layers.Conv2DTranspose(32, kernel_size=(4, 1), strides=2, activation=tf.nn.relu)(Dense_decoder1_npy_reshape)  # =>(26, 26, 128)
    batch_norm_trans_2 = tf.keras.layers.BatchNormalization()(conv_trans_2)
    print('batch_norm_trans_2: {}'.format(batch_norm_trans_2.shape))

    conv_trans_3 = tf.keras.layers.Conv2DTranspose(45, kernel_size=(3, 1),strides=2,activation=tf.nn.relu)(batch_norm_trans_2)  # =>(26, 26, 128)
    batch_norm_trans_3 = tf.keras.layers.BatchNormalization()(conv_trans_3)
    print('batch_norm_trans_3: {}'.format(batch_norm_trans_3.shape))


    output00_1 = tf.keras.layers.Conv2DTranspose(1, kernel_size=(3, 1), strides=1, activation='sigmoid')(batch_norm_trans_3)  # =>(26, 26, 128)
    print('output00_1: {}'.format(output00_1.shape))
    output00_1_1 = tf.keras.layers.Reshape((output00_1.shape[1]*output00_1.shape[2], 1, 1))(output00_1)
    print('output00_1_1: {}'.format(output00_1_1.shape))
    # output01 = tf.keras.layers.Conv2DTranspose(1, kernel_size=(3, 1), strides=1, activation='sigmoid')(output00_1_1)
    # print('output01: {}'.format(output01.shape))
    # output layer
    output02 = tf.squeeze(output00_1_1,  axis=2)
    output_p = tf.squeeze(output02, axis=2)
    print('output_p: {}'.format(output_p.shape))

    autoencoder_p = tf.keras.Model(inputs=x_input_p, outputs=output_p)
    # define autoencoder model
    weight_decay = 1e-4
    opt = tf.keras.optimizers.Adam(lr= lr_p, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay= weight_decay) #
    autoencoder_p.compile(optimizer = opt,loss='mse')

    # fit the autoencoder model to reconstruct input
    history = autoencoder_p.fit(data_train_p_npy, data_train_p_npy, epochs=epoch_p, batch_size=batch_size_p, verbose=2, validation_data=(data_vaild_p_npy,data_vaild_p_npy))

    encoder_pro = tf.keras.Model(inputs=x_input_p, outputs=encode_r)
    encoder_pro.save( save_path + '/encoder_pro.h5')
    
    return encoder_pro


def compound_encoder(data_train_r_npy,data_vaild_r_npy,lr_r,epoch_r,batch_size_r,save_path):
    # lr = 1e-4
    # epoch = 100
    # batch_size = 64

    input_shape = data_train_r_npy.shape[1]
    x_input_r = Input(input_shape)
    x_input01 = tf.expand_dims(x_input_r, axis=2)
    x_input01 = tf.expand_dims(x_input01, axis=3)
    print('x_input: {}'.format(x_input01.shape))

    # ENCODER
    conv_1 = tf.keras.layers.Conv2D(filters=64,  kernel_size=(3, 1), padding="valid",  activation=tf.nn.relu)(x_input01)
    batch_norm_1 = tf.keras.layers.BatchNormalization()(conv_1)
    pooling_1 = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=2, padding='valid')(batch_norm_1)
    print('pooling_1: {}'.format(pooling_1.shape))

    conv_2 = tf.keras.layers.Conv2D(filters=48,  kernel_size=(3, 1), padding="valid",  activation=tf.nn.relu)(pooling_1)
    # print(conv_2.shape)
    batch_norm_2 = tf.keras.layers.BatchNormalization()(conv_2)
    pooling_2 = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=4, padding='valid')(batch_norm_2)
    print('pooling_2: {}'.format(pooling_2.shape))

    conv_3 = tf.keras.layers.Conv2D(filters=32,  kernel_size=(3, 1), padding="valid",  activation=tf.nn.relu)(pooling_2)
    # print(conv_3.shape)
    batch_norm_3 = tf.keras.layers.BatchNormalization()(conv_3)
    pooling_3 = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=2, padding='valid')(batch_norm_3)
    print('pooling_3: {}'.format(pooling_3.shape))

    conv_4 = tf.keras.layers.Conv2D(filters=16,  kernel_size=(3, 1), padding="valid",  activation=tf.nn.relu)(pooling_3)
    # print(conv_3.shape)
    batch_norm_4 = tf.keras.layers.BatchNormalization()(conv_4)
    pooling_4 = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=2, padding='valid')(batch_norm_4)
    print('pooling_4: {}'.format(pooling_4.shape))

    conv_5 = tf.keras.layers.Conv2D(filters=8,  kernel_size=(3, 1), padding="valid",  activation=tf.nn.relu)(pooling_4)
    # print(conv_3.shape)
    batch_norm_5 = tf.keras.layers.BatchNormalization()(conv_5)
    pooling_5 = tf.keras.layers.MaxPool2D(pool_size=(2, 1), strides=2, padding='valid')(batch_norm_5)
    print('pooling_5: {}'.format(pooling_5.shape))

    Flatten1 = tf.keras.layers.Flatten()(pooling_5)
    print('Flatten1: {}'.format(Flatten1.shape))

    Dense_encoder1 = tf.keras.layers.Dense(units=100, activation=tf.nn.relu)(Flatten1)
    print('Dense_encoder1: {}'.format(Dense_encoder1.shape))

    encode_r = Dense_encoder1
    print('encode_r: {}'.format(encode_r.shape))

    Dense_decoder1 = tf.keras.layers.Dense(units=Flatten1.shape[1], activation=tf.nn.relu)(encode_r)
    print('Dense_decoder1: {}'.format(Dense_decoder1.shape))
    print('Dense_decoder1 dtype: {}'.format(Dense_decoder1.dtype))
    
    Dense_decoder1_npy_reshape = tf.keras.layers.Reshape((34, 1, 8))(Dense_decoder1)
    print('Dense_decoder1_npy_reshape: {}'.format(Dense_decoder1_npy_reshape.shape))

    # DECODER
    conv_trans_1 = tf.keras.layers.Conv2DTranspose(16, kernel_size=(5, 1), strides=2, activation=tf.nn.relu)(Dense_decoder1_npy_reshape)  # =>(26, 26, 128)
    batch_norm_trans_1 = tf.keras.layers.BatchNormalization()(conv_trans_1)
    print('batch_norm_trans_1: {}'.format(batch_norm_trans_1.shape))

    conv_trans_2 = tf.keras.layers.Conv2DTranspose(32, kernel_size=(3, 1), strides=2, activation=tf.nn.relu)(batch_norm_trans_1)  # =>(26, 26, 128)
    batch_norm_trans_2 = tf.keras.layers.BatchNormalization()(conv_trans_2)
    print('batch_norm_trans_2: {}'.format(batch_norm_trans_2.shape))

    conv_trans_3 = tf.keras.layers.Conv2DTranspose(45, kernel_size=(3, 1),strides=2,activation=tf.nn.relu)(batch_norm_trans_2)  # =>(26, 26, 128)
    batch_norm_trans_3 = tf.keras.layers.BatchNormalization()(conv_trans_3)
    print('batch_norm_trans_3: {}'.format(batch_norm_trans_3.shape))

    output00_1 = tf.keras.layers.Conv2DTranspose(1, kernel_size=(4, 1), strides=1, activation=tf.nn.relu)(batch_norm_trans_3)  # =>(26, 26, 128)
    print('output00_1: {}'.format(output00_1.shape))
    output00_1_1 = tf.keras.layers.Reshape((output00_1.shape[1]*output00_1.shape[2], 1, 1))(output00_1)
    print('output00_1_1: {}'.format(output00_1_1.shape))
    output01 = tf.keras.layers.Conv2DTranspose(1, kernel_size=(5, 1), strides=1, activation='sigmoid')(output00_1_1)
    print('output01: {}'.format(output01.shape))
    # output layer
    output02 = tf.squeeze(output01,  axis=2)
    output = tf.squeeze(output02, axis=2)
    print('output: {}'.format(output.shape))
    autoencoder = tf.keras.Model(inputs=x_input_r, outputs=output)
    # define autoencoder model

    weight_decay = 1e-4
    opt = tf.keras.optimizers.Adam(lr= lr_r, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay= weight_decay) #
    autoencoder.compile(optimizer = opt,loss='mse')

    # fit the autoencoder model to reconstruct input
    history = autoencoder.fit(data_train_r_npy, data_train_r_npy, epochs=epoch_r, batch_size=batch_size_r, verbose=2, validation_data=(data_vaild_r_npy,data_vaild_r_npy))

    encoder_rna = tf.keras.Model(inputs=x_input_r, outputs=encode_r)
    encoder_rna.save( save_path + '/encoder_compound.h5')
    return encoder_rna

def train_auto_encoder(X_train, layers, lr_r, drop_rate,batch_size=100, nb_epoch=100, activation='sigmoid'):
    trained_encoders = []
    trained_decoders = []
    X_train_tmp = np.copy(X_train)
    # X_test_tmp = np.copy(X_test)
  

    for n_in, n_out in zip(layers[:-1], layers[1:]):
        print('Pre-training the layer: Input {} -> Output {}'.format(n_in, n_out))
        ae = Sequential(
            [Dense(n_out, input_dim=X_train_tmp.shape[1], activation=activation, ),
             Dense(n_in, activation=activation),
             Dropout(drop_rate)]
        )

        weight_decay = 1e-4
        opt = tf.keras.optimizers.Adam(lr= lr_r, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay= weight_decay) #

        ae.compile(loss='mean_squared_error', optimizer=opt)
        ae.fit(X_train_tmp, X_train_tmp, batch_size=batch_size, epochs=nb_epoch, verbose=0, shuffle=True)
        # store trained encoder
        trained_encoders.append(ae.layers[0])
        trained_decoders.append(ae.layers[1])
        # update training data
        encoder = Sequential([ae.layers[0]])
        # encoder.evaluate(X_train_tmp, X_train_tmp, batch_size=batch_size)
        X_train_tmp = encoder.predict(X_train_tmp)
        # X_test_tmp = encoder.predict(X_test_tmp)
    
    return trained_encoders

def conv_block(x, growth_rate):
    x1 = Conv1D(growth_rate, 3, padding='same')(x)
    x1 = BatchNormalization()(x1)
    x1 = Activation('relu')(x1)
    x = concatenate(([x, x1]), axis=-1)
    return x

def dense_block(x, blocks, size):
    if blocks == 0:
        x = x
    elif blocks > 0:
        for i in range(blocks):
            x = conv_block(x, size)
    return x

def conjoint_sae(encoders_rna, encoders_protein, rna_coding_length, pro_coding_length, drop_rate):

    Blocks = 3
    Block_size = 18
    Filter_num = 24
    Kernel_sizes = 6

    # NN for RNA feature analysis
    xr_in_conjoint = Input(shape=(rna_coding_length,))
    xr_encoded = encoders_rna[0](xr_in_conjoint)
    xr_encoded = Dropout(drop_rate)(xr_encoded)
    xr_encoded = encoders_rna[1](xr_encoded)
    xr_encoded = Dropout(drop_rate)(xr_encoded)
    xr_encoded = encoders_rna[2](xr_encoded)
    xr_encoder = Dropout(drop_rate)(xr_encoded)
    xr_encoder = BatchNormalization()(xr_encoder)
    # xr_encoder = PReLU()(xr_encoder)
    xr_encoder = Dropout(drop_rate)(xr_encoder)

    # NN for protein feature analysis
    xp_in_conjoint = Input(shape=(pro_coding_length,))
    xp_encoded = encoders_protein[0](xp_in_conjoint)
    xp_encoded = Dropout(drop_rate)(xp_encoded)
    xp_encoded = encoders_protein[1](xp_encoded)
    xp_encoded = Dropout(drop_rate)(xp_encoded)
    xp_encoder = encoders_protein[2](xp_encoded)
    xp_encoder = Dropout(drop_rate)(xp_encoder)
    xp_encoder = BatchNormalization()(xp_encoder)
    # xp_encoder = PReLU()(xp_encoder)
    xp_encoder = Dropout(drop_rate)(xp_encoder)

    x_out_conjoint_combine = concatenate([xp_encoder, xr_encoder], name="concatenate001")
    print('x_out_conjoint_combine: {}'.format(x_out_conjoint_combine.shape))
    x_input01 = tf.expand_dims(x_out_conjoint_combine, axis=2)
    # x_input01 = tf.expand_dims(x_input01, axis=3)
    print('x_input01: {}'.format(x_input01.shape))

    seq_in_conjoint = Conv1D(filters=Filter_num + Blocks * Block_size, kernel_size=1)(x_input01)
    seq_in_conjoint = BatchNormalization()(seq_in_conjoint)
    seq_in_conjoint = Activation('relu')(seq_in_conjoint)
    print('seq_in_conjoint: {}'.format(seq_in_conjoint.shape))

    seq_cnn = Conv1D(filters=Filter_num, kernel_size=Kernel_sizes, padding='same')(seq_in_conjoint)
    seq_cnn = BatchNormalization()(seq_cnn)
    seq_cnn = Activation('relu')(seq_cnn)
    print('seq_cnn: {}'.format(seq_cnn.shape))

    seq_cnn = MaxPooling1D(pool_size=Kernel_sizes)(seq_cnn)
    seq_cnn = dense_block(seq_cnn, Blocks, Block_size)
    seq_cnn = BatchNormalization()(seq_cnn)

    print('seq_cnn.shape:{}'.format(seq_cnn.shape))
    x_out_conjoint = Flatten()(seq_cnn)
    print('x_out_conjoint.shape:{}'.format(x_out_conjoint.shape))

    x_out_conjoint = Dense(192, activation='relu', use_bias=False)(
        x_out_conjoint)
    x_out_conjoint = Dropout(drop_rate)(x_out_conjoint)
    x_out_conjoint = Dense(64, activation='relu', use_bias=False)(
        x_out_conjoint)
    x_out_conjoint = Dropout(drop_rate)(x_out_conjoint)
    y_conjoint = Dense(2, activation='softmax')(x_out_conjoint)

    model_conjoint = Model(inputs=[xr_in_conjoint, xp_in_conjoint], outputs=y_conjoint)

    return model_conjoint,x_out_conjoint_combine

class CLA_EarlyStoppingAndPerformance(tf.keras.callbacks.Callback):

    def __init__(self, train_data, valid_data, MASK = -1, patience=5, criteria = 'val_loss', metric = 'ROC', last_avf = None, verbose = 0):
        super(CLA_EarlyStoppingAndPerformance, self).__init__()
        
        sp = ['val_loss', 'val_auc']
        assert criteria in sp, 'not support %s ! only %s' % (criteria, sp)
        self.x, self.y  = train_data
        self.x_val, self.y_val = valid_data
        self.last_avf = last_avf
        
        self.history = {'loss':[],
                        'val_loss':[],
                        'auc':[],
                        'val_auc':[],
                        
                        'epoch':[]}
        self.MASK = MASK
        self.patience = patience
        # best_weights to store the weights at which the minimum loss occurs.
        self.best_weights = None
        self.criteria = criteria
        self.metric = metric
        self.best_epoch = 0
        self.verbose = verbose
        
    def sigmoid(self, x):
        s = 1/(1+np.exp(-x))
        return s

    
    def roc_auc(self, y_true, y_pred):
        if self.last_avf == None:
            y_pred_logits = self.sigmoid(y_pred)
        else:
            y_pred_logits = y_pred
            
        N_classes = y_pred_logits.shape[1]

        aucs = []
        for i in range(N_classes):
            y_pred_one_class = y_pred_logits[:,i]
            y_true_one_class = y_true[:, i]
            mask = ~(y_true_one_class == self.MASK)
            try:
                if self.metric == 'ROC':
                    auc = roc_auc_score(y_true_one_class[mask], y_pred_one_class[mask]) #ROC_AUC
                elif self.metric == 'PRC': 
                    auc = prc_auc_score(y_true_one_class[mask], y_pred_one_class[mask]) #PRC_AUC
                elif self.metric == 'ACC':
                    auc = accuracy_score(y_true_one_class[mask], np.round(y_pred_one_class[mask])) #ACC
            except:
                auc = np.nan
            aucs.append(auc)
        return aucs  
    
        
        
    def on_train_begin(self, logs=None):
        # The number of epoch it has waited when loss is no longer minimum.
        self.wait = 0
        # The epoch the training stops at.
        self.stopped_epoch = 0
        # Initialize the best as infinity.
        if self.criteria == 'val_loss':
            self.best = np.Inf  
        else:
            self.best = -np.Inf
               
    def on_epoch_end(self, epoch, logs={}):
        
        y_pred = self.model.predict(self.x)
        roc_list = self.roc_auc(self.y, y_pred)
        roc_mean = np.nanmean(roc_list)
        
        y_pred_val = self.model.predict(self.x_val)
        roc_val_list = self.roc_auc(self.y_val, y_pred_val)        
        roc_val_mean = np.nanmean(roc_val_list)
        
        self.history['loss'].append(logs.get('loss'))
        self.history['val_loss'].append(logs.get('val_loss'))
        self.history['auc'].append(roc_mean)
        self.history['val_auc'].append(roc_val_mean)
        self.history['epoch'].append(epoch)
        
        
        eph = str(epoch+1).zfill(4)        
        loss = '{0:.4f}'.format((logs.get('loss')))
        val_loss = '{0:.4f}'.format((logs.get('val_loss')))
        auc = '{0:.4f}'.format(roc_mean)
        auc_val = '{0:.4f}'.format(roc_val_mean)    
        
        if self.verbose:
            if self.metric == 'ACC':
                print('\repoch: %s, loss: %s - val_loss: %s; acc: %s - val_acc: %s' % (eph,
                                                                                   loss, 
                                                                                   val_loss, 
                                                                                   auc,
                                                                                   auc_val), end=100*' '+'\n')

            else:
                print('\repoch: %s, loss: %s - val_loss: %s; auc: %s - val_auc: %s' % (eph,
                                                                                   loss, 
                                                                                   val_loss, 
                                                                                   auc,
                                                                                   auc_val), end=100*' '+'\n')


        if self.criteria == 'val_loss':
            current = logs.get(self.criteria)
            if current <= self.best:
                self.best = current
                self.wait = 0
                # Record the best weights if current results is better (less).
                self.best_weights = self.model.get_weights()
                self.best_epoch = epoch

            else:
                self.wait += 1
                if self.wait >= self.patience:
                    self.stopped_epoch = epoch
                    self.model.stop_training = True
                    print('\nRestoring model weights from the end of the best epoch.')
                    self.model.set_weights(self.best_weights)    
                    
        else:
            current = roc_val_mean
            if current >= self.best:
                self.best = current
                self.wait = 0
                # Record the best weights if current results is better (less).
                self.best_weights = self.model.get_weights()
                self.best_epoch = epoch

            else:
                self.wait += 1
                if self.wait >= self.patience:
                    self.stopped_epoch = epoch
                    self.model.stop_training = True
                    print('\nRestoring model weights from the end of the best epoch.')
                    self.model.set_weights(self.best_weights)              
    
    def on_train_end(self, logs=None):
        self.model.set_weights(self.best_weights)
        if self.stopped_epoch > 0:
            print('\nEpoch %05d: early stopping' % (self.stopped_epoch + 1))

        
    def evaluate(self, testX, testY):
        
        y_pred = self.model.predict(testX)
        roc_list = self.roc_auc(testY, y_pred)
        return roc_list 

def split_training_validation_test(classes00, size=[0.75,0.1,0.15], shuffle=True, SEED = 888):
    """split sampels based on balnace classes"""
    random.seed(SEED)
    num_samples = len(classes00)

    classes01 = np.array(classes00)    
    classes = np.squeeze(classes01)

    classes_unique = np.unique(classes)

    indices = np.arange(num_samples)

    training_indice = []
    validation_indice = []
    test_indice = []

    for cl in classes_unique:
        indices_cl = indices[classes == cl]

        num_samples_cl = len(indices_cl)
       
        
        if shuffle:
            random.shuffle(indices_cl)  # in-place shuffle

        # module and residual
        num_samples_training = int(num_samples_cl * size[0])
        num_samples_validation = int(num_samples_cl * size[1])
        num_samples_test = int(num_samples_cl * size[2])

        test_indice = test_indice + [test for test in indices_cl[0:num_samples_test]]
        validation_indice = validation_indice + [valid for valid in indices_cl[num_samples_test:(num_samples_test + num_samples_validation)]]
        training_indice = training_indice + [train for train in indices_cl[(num_samples_test + num_samples_validation):]]

    random.shuffle(training_indice)
    random.shuffle(validation_indice)
    random.shuffle(test_indice)

    print('training_indice length: {},validation_indice length: {},test_indice length: {}, '.format(len(training_indice),len(validation_indice), len(test_indice)))
    return training_indice, validation_indice, test_indice

# =========================================================================================================
# start the classification

def evaluation_method(datapath,label_path,resultpath,type = 'RNAonly',com_num = 16,modelnm = 'RF'):
    # make result file path
    resultpath = resultpath + '/classification_result'
    os.makedirs(resultpath, exist_ok= True)
    # get all combination methods
    combination_methods = make_comcod(com_num)
    # get label
    labeldata = pd.read_csv(label_path)
    if type =='RNAonly':
        label = np.array(labeldata.iloc[:,1])
    else:
        label = np.array(labeldata.iloc[:,2])
    # check the label classes
    if len(set(label))==1:
        print('The sample class is only one, so can not classificate!')
        exit()

    res_eval = {}
    combin_num = 0
    combine_X = mknpy_RNAonly(combination_methods[combin_num],datapath)

    print('The {}/{} combination of coding methods is {method}'.format(combin_num+1,len(combination_methods),method = combination_methods[combin_num]))
    print('The data shape:{}'.format(combine_X.shape))
    # split the training and test data
    train_x, valid_x, train_y, valid_y = train_test_split(combine_X, label, test_size=0.2, random_state=42,stratify=label)
    # normalization data
    scaler = preprocessing.StandardScaler().fit(train_x)
    train_x = scaler.transform(train_x)
    valid_x = scaler.transform(valid_x)

    # embedded feature extracted by CNN-AE==========================================================
    # Specify sufficient boosting iterations to reach a minimum
    batch_size_rs = [32,64,96,128] # 8,16,48,80,16,32,64,96,128
    lr_rs = [0.0001,0.001,0.01]
    epoch_rs = [100]

    parameters_all = [[x,y,z] for x in batch_size_rs for y in lr_rs for z in epoch_rs] 
    # parameters = random.sample(parameters_all, 10)
    parameters = parameters_all

    for index_emb, parameter_emb in enumerate(parameters):
        print('Current parameter_emb is {}/{} dropout, learning rate, batchsize: {}'.format(index_emb,len(parameter_emb),parameter_emb))
        
        save_path = resultpath + '/Para_' + '_'.join([str(i) for i in parameter_emb])
        os.makedirs(save_path, exist_ok = True)

        batch_size_r = parameter_emb[0]
        lr_r = parameter_emb[1]
        epoch_r = parameter_emb[2]
    
        encoder_r = rna_encoder_1(train_x,valid_x,lr_r,epoch_r,batch_size_r,save_path)

        train_x_encode = encoder_r.predict(train_x)
        valid_x_encode = encoder_r.predict(valid_x)               

        print('train_x_encode.shape: {}'.format(train_x_encode.shape))
        print('valid_x_encode.shape: {}'.format(valid_x_encode.shape))

        np.save(save_path + '/' + 'train_x_encode.npy', train_x_encode)
        np.save(save_path + '/' + 'valid_x_encode.npy', valid_x_encode)

        # starting traing       
        print('The classification model is: {}'.format(modelnm))
        if modelnm == 'svm':
            model,best_para = svm_two(train_x_encode, valid_x_encode, train_y, valid_y)
            group = model.predict(valid_x_encode) # test
            score = model.predict_proba(valid_x_encode) # get the confidence probability
            Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(valid_y,score[:, 1],group)        
        elif modelnm == 'RF':
            model,best_para = RF(train_x_encode, valid_x_encode, train_y, valid_y)
            group = model.predict(valid_x_encode) # test
            score = model.predict_proba(valid_x_encode) # get the confidence probability
            Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(valid_y,score[:, 1],group)
        elif modelnm == 'xgboost':
            Acc, Sn, Sp, Pre, MCC, AUC,best_para = xgboost_model(train_x_encode,train_y,valid_x_encode,valid_y)
        elif modelnm == 'DNN' or modelnm == 'CNN':
            Acc, Sn, Sp, Pre, MCC, AUC,best_para = train_model(train_x_encode,train_y,valid_x_encode,valid_y,modelnm = modelnm)

        res_eval[index_emb+1] = [modelnm,parameter_emb,best_para, Acc, Sn, Sp, Pre, MCC, AUC]
    indexs = ['model name','embedding parameters','best parameters','Acc','Sn','Sp','Pre','MCC','AUC']
    evaluation_result = pd.DataFrame(res_eval,index = indexs).T
    evaluation_result = evaluation_result.sort_values(by="MCC" , ascending=False)
    evaluation_result.to_csv(resultpath + '/Evaluation_result.csv')
    print(evaluation_result)

# 
def evaluation_interaction_cnn_ae(datapath,label_path,resultpath,type = 'RNAonly',com_num = 16,modelnm = 'RF'):
    # make result file path
    resultpath = resultpath + '/classification_result'
    os.makedirs(resultpath, exist_ok= True)
    # get all combination methods
    combination_methods = make_comcod(com_num)
    # get label
    labeldata = pd.read_csv(label_path)
    label = np.array(labeldata.iloc[:,2])
    # check the label classes
    if len(set(label))==1:
        print('The sample class is only one, so can not classificate!')
        exit()

    res_eval = {}
    combin_num = 0

    if type =='RNA-RNA':
        combine_X_A, combine_X_B = mknpy_RNA_RNA(combination_methods[combin_num],datapath)
    elif type =='RNA-pro':
        combine_X_A, combine_X_B = mknpy_RNA_pro(combination_methods[combin_num],datapath)
    elif type =='RNA-compound':
        combine_X_A, combine_X_B = mknpy_RNA_compound(combination_methods[combin_num],datapath)
    print('The {}/{} combination of coding methods is {method}'.format(combin_num+1,len(combination_methods),method = combination_methods[combin_num]))
    print('The data A shape:{}, The data B shape:{} '.format(combine_X_A.shape, combine_X_B.shape))

    train_label = labeldata.iloc[:,2].tolist()
    train_idx, valid_idx, test_idx = split_training_validation_test(train_label, size=[0.8,0.2,0])    
    
    train_X_A = combine_X_A[train_idx,:]
    valid_X_A = combine_X_A[valid_idx,:]

    train_X_B = combine_X_B[train_idx,:]
    valid_X_B = combine_X_B[valid_idx,:]

    train_y = np.array([int(train_label[i]) for i in train_idx])
    valid_y = np.array([int(train_label[i]) for i in valid_idx])    

    # normalization
    scaler_r = MinMaxScaler().fit(train_X_A)
    train_feature_RNA = scaler_r.transform(train_X_A)
    test_feature_RNA = scaler_r.transform(valid_X_A)

    scaler_p = MinMaxScaler().fit(train_X_B)
    train_feature_pro = scaler_p.transform(train_X_B)
    test_feature_pro = scaler_p.transform(valid_X_B)

    # Specify sufficient boosting iterations to reach a minimum
    batch_size_rs = [32,64,96,128] # 8,16,48,80,16,32,64,96,128
    lr_rs = [0.0001,0.001,0.01]
    epoch_rs = [100]

    parameters_all = [[x,y,z] for x in batch_size_rs for y in lr_rs for z in epoch_rs] 
    # parameters = random.sample(parameters_all, 10)
    parameters = parameters_all

    for index_emb, parameter_emb in enumerate(parameters):
        print('Current parameter_emb is {}/{} dropout, learning rate, batchsize: {}'.format(index_emb,len(parameter_emb),parameter_emb))
        
        save_path = resultpath + '/Para_' + '_'.join([str(i) for i in parameter_emb])
        os.makedirs(save_path, exist_ok = True)

        batch_size_r = parameter_emb[0]
        lr_r = parameter_emb[1]
        epoch_r = parameter_emb[2]

        if type =='RNA-RNA':
            encoder_r = rna_encoder_1(train_feature_RNA,test_feature_RNA,lr_r,epoch_r,batch_size_r,save_path)
            encoder_p = rna_encoder_2(train_feature_pro,test_feature_pro,lr_r,epoch_r,batch_size_r,save_path)
        elif type =='RNA-pro':            
            encoder_r = rna_encoder_1(train_feature_RNA,test_feature_RNA,lr_r,epoch_r,batch_size_r,save_path)
            encoder_p = pro_encoder(train_feature_pro,test_feature_pro,lr_r,epoch_r,batch_size_r,save_path)
            
        elif type =='RNA-compound':            
            encoder_p = compound_encoder(train_feature_pro,test_feature_pro,lr_r,epoch_r,batch_size_r,save_path)
            encoder_r = rna_encoder_1(train_feature_RNA,test_feature_RNA,lr_r,epoch_r,batch_size_r,save_path)           

        train_feature_RNA_encode = encoder_r.predict(train_feature_RNA)
        test_feature_RNA_encode = encoder_r.predict(test_feature_RNA)
        train_feature_pro_encode = encoder_p.predict(train_feature_pro)
        test_feature_pro_encode = encoder_p.predict(test_feature_pro)                    

        print('train_feature_RNA_encode.shape: {}'.format(train_feature_RNA_encode.shape))
        print('test_feature_RNA_encode.shape: {}'.format(test_feature_RNA_encode.shape))
        print('train_feature_B_encode.shape: {}'.format(train_feature_pro_encode.shape))
        print('test_feature_B_encode.shape: {}'.format(test_feature_pro_encode.shape))


        np.save(save_path + '/' + 'train_feature_RNA_encode.npy', train_feature_RNA_encode)
        np.save(save_path + '/' + 'test_feature_RNA_encode.npy', test_feature_RNA_encode)
        np.save(save_path + '/' + 'train_feature_B_encode.npy', train_feature_pro_encode)
        np.save(save_path + '/' + 'test_feature_B_encode.npy', test_feature_pro_encode)

        train_x = np.concatenate((train_feature_RNA_encode, train_feature_pro_encode), axis = 1)
        valid_x = np.concatenate((test_feature_RNA_encode, test_feature_pro_encode), axis = 1)

        # starting traing       
        print('The classification model is: {}'.format(modelnm))
        if modelnm == 'svm':
            model,best_para = svm_two(train_x, valid_x, train_y, valid_y)
            group = model.predict(valid_x) # test
            score = model.predict_proba(valid_x) # get the confidence probability
            Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(valid_y,score[:, 1],group)        
        elif modelnm == 'RF':
            model,best_para = RF(train_x, valid_x, train_y, valid_y)
            group = model.predict(valid_x) # test
            score = model.predict_proba(valid_x) # get the confidence probability
            Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(valid_y,score[:, 1],group)
        elif modelnm == 'xgboost':
            Acc, Sn, Sp, Pre, MCC, AUC,best_para = xgboost_model(train_x,train_y,valid_x,valid_y)
        elif modelnm == 'DNN' or modelnm == 'CNN':
            Acc, Sn, Sp, Pre, MCC, AUC,best_para = train_model(train_x,train_y,valid_x,valid_y,modelnm = modelnm)

        res_eval[index_emb+1] = [modelnm,parameter_emb,best_para, Acc, Sn, Sp, Pre, MCC, AUC]
    indexs = ['model name','embedding parameters','best parameters','Acc','Sn','Sp','Pre','MCC','AUC']
    evaluation_result = pd.DataFrame(res_eval,index = indexs).T
    evaluation_result = evaluation_result.sort_values(by="MCC" , ascending=False)
    evaluation_result.to_csv(resultpath + '/Evaluation_result.csv')
    print(evaluation_result)

def evaluation_interaction(datapath,label_path,resultpath,type = 'RNAonly',com_num = 16,modelnm = 'RF'):
    # make result file path
    resultpath = resultpath + '/classification_result'
    os.makedirs(resultpath, exist_ok= True)
    # get all combination methods
    combination_methods = make_comcod(com_num)
    # get label
    labeldata = pd.read_csv(label_path)
    label = np.array(labeldata.iloc[:,2])
    # check the label classes
    if len(set(label))==1:
        print('The sample class is only one, so can not classificate!')
        exit()

    res_eval = {}
    combin_num = 0

    if type =='RNA-RNA':
        combine_X_A, combine_X_B = mknpy_RNA_RNA(combination_methods[combin_num],datapath)
    elif type =='RNA-pro':
        combine_X_A, combine_X_B = mknpy_RNA_pro(combination_methods[combin_num],datapath)
    elif type =='RNA-compound':
        combine_X_A, combine_X_B = mknpy_RNA_compound(combination_methods[combin_num],datapath)
    print('The {}/{} combination of coding methods is {method}'.format(combin_num+1,len(combination_methods),method = combination_methods[combin_num]))
    print('The data A shape:{}, The data B shape:{} '.format(combine_X_A.shape, combine_X_B.shape))

    train_label = labeldata.iloc[:,2].tolist()
    train_idx, valid_idx, test_idx = split_training_validation_test(train_label, size=[0.8,0.2,0])    
    
    train_X_A = combine_X_A[train_idx,:]
    valid_X_A = combine_X_A[valid_idx,:]

    train_X_B = combine_X_B[train_idx,:]
    valid_X_B = combine_X_B[valid_idx,:]

    train_y = np.array([int(train_label[i]) for i in train_idx])
    valid_y = np.array([int(train_label[i]) for i in valid_idx])    

    # normalization
    scaler_r = MinMaxScaler().fit(train_X_A)
    train_feature_RNA = scaler_r.transform(train_X_A)
    test_feature_RNA = scaler_r.transform(valid_X_A)

    scaler_p = MinMaxScaler().fit(train_X_B)
    train_feature_pro = scaler_p.transform(train_X_B)
    test_feature_pro = scaler_p.transform(valid_X_B)

    # Specify sufficient boosting iterations to reach a minimum
    batch_size_rs = [32,64,96,128] # 8,16,48,80,16,32,64,96,128
    lr_rs = [0.0001,0.001,0.01]
    epoch_rs = [100]

    parameters_all = [[x,y,z] for x in batch_size_rs for y in lr_rs for z in epoch_rs] 
    # parameters = random.sample(parameters_all, 10)
    parameters = parameters_all

    for index_emb, parameter_emb in enumerate(parameters):
        print('Current parameter_emb is {}/{} dropout, learning rate, batchsize: {}'.format(index_emb,len(parameter_emb),parameter_emb))
        
        save_path = resultpath + '/Para_' + '_'.join([str(i) for i in parameter_emb])
        os.makedirs(save_path, exist_ok = True)

        batch_size_r = parameter_emb[0]
        lr_r = parameter_emb[1]
        epoch_r = parameter_emb[2]

        # if type =='RNA-RNA':

        # =====================================================================
        # start SAE-based embedded feature extraction

        X_train_RNA_sub = train_feature_RNA
        X_train_pro = train_feature_pro
        test_RNA_sub = test_feature_RNA

        trainX = (X_train_RNA_sub, X_train_pro) 
        validX = (test_RNA_sub, test_feature_pro)



        batch_size = parameter_emb[0]
        lr_r = parameter_emb[1]
        drop_rate = parameter_emb[2]
        epoch_num = 100

        layers_RNA = [X_train_RNA_sub.shape[1], 512, 256, 128]
        encoders_rna =  train_auto_encoder(X_train=X_train_RNA_sub,layers=layers_RNA, lr_r = lr_r, drop_rate = drop_rate,batch_size=batch_size, nb_epoch=epoch_num)

        layers_pro = [X_train_pro.shape[1], 512, 256, 128] #2048, 512, 64
        encoders_pro =  train_auto_encoder(X_train=X_train_pro,layers=layers_pro, lr_r = lr_r, drop_rate = drop_rate,batch_size=batch_size, nb_epoch=epoch_num)




        # model parameters
        batch_size_clf = parameter_emb[0]
        lr_clf = parameter_emb[1]
        drop_rate_clf = 0.2



        epochs = 200 #300
        patience = 20 #30 early stopping                
        weight_decay = 1e-4
        monitor = 'val_auc'
        metric = 'ROC'


        pro_coding_length = X_train_pro.shape[1]
        rna_coding_length = X_train_RNA_sub.shape[1] 
        # create model
        clf, x_out_conjoint = conjoint_sae(encoders_rna,encoders_pro,  rna_coding_length, pro_coding_length,drop_rate_clf)


        opt = tf.keras.optimizers.Adam(lr= lr_clf, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay= weight_decay) #
        clf.compile(optimizer = opt, loss = 'categorical_crossentropy', metrics = ['accuracy'])
        performance = CLA_EarlyStoppingAndPerformance((trainX, train_y), 
                                                                (validX, valid_y), 
                                                                patience = patience,
                                                                metric = metric,
                                                                criteria = monitor)
        # fit your model
        clf.fit(trainX, train_y, batch_size= batch_size_clf, 
                                epochs= epochs, verbose= 1, shuffle = True, 
                                validation_data = (validX, valid_y), 
                                callbacks=[performance]) 

        SAE_model = Model( inputs = clf.input, outputs = clf.get_layer('concatenate001').output )
        train_feature_RNA_encode = SAE_model.predict(trainX)
        test_feature_RNA_encode = SAE_model.predict(validX)
        # end SAE-based embedded feature extraction
        # =====================================================================
                     

        print('train_feature_encode.shape: {}'.format(train_feature_RNA_encode.shape))
        print('test_feature_RNA_encode.shape: {}'.format(test_feature_RNA_encode.shape))

        np.save(save_path + '/' + 'train_feature_encode.npy', train_feature_RNA_encode)
        np.save(save_path + '/' + 'test_feature_encode.npy', test_feature_RNA_encode)

        train_x = train_feature_RNA_encode
        valid_x = test_feature_RNA_encode

        # starting traing       
        print('The classification model is: {}'.format(modelnm))
        if modelnm == 'svm':
            model,best_para = svm_two(train_x, valid_x, train_y, valid_y)
            group = model.predict(valid_x) # test
            score = model.predict_proba(valid_x) # get the confidence probability
            Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(valid_y,score[:, 1],group)        
        elif modelnm == 'RF':
            model,best_para = RF(train_x, valid_x, train_y, valid_y)
            group = model.predict(valid_x) # test
            score = model.predict_proba(valid_x) # get the confidence probability
            Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(valid_y,score[:, 1],group)
        elif modelnm == 'xgboost':
            Acc, Sn, Sp, Pre, MCC, AUC,best_para = xgboost_model(train_x,train_y,valid_x,valid_y)
        elif modelnm == 'DNN' or modelnm == 'CNN':
            Acc, Sn, Sp, Pre, MCC, AUC,best_para = train_model(train_x,train_y,valid_x,valid_y,modelnm = modelnm)

        res_eval[index_emb+1] = [modelnm,parameter_emb,best_para, Acc, Sn, Sp, Pre, MCC, AUC]
    indexs = ['model name','embedding parameters','best parameters','Acc','Sn','Sp','Pre','MCC','AUC']
    evaluation_result = pd.DataFrame(res_eval,index = indexs).T
    evaluation_result = evaluation_result.sort_values(by="MCC" , ascending=False)
    evaluation_result.to_csv(resultpath + '/Evaluation_result.csv')
    print(evaluation_result)

# ========================================================================================

# datapath = '../out/RNA-pro/encoding_features'
# label_path = '../demo/RNA-Protein/RNA-Protein-Interacting.csv'
# resultpath = '../out/RNA-pro/classification_result'

# datapath = '../out/RNA-compound/encoding_features'
# label_path = '../demo/RNA-compound/RNA-small-molecule-Interacting.csv'
# resultpath = '../out/RNA-compound/classification_result'

# # # modelnm = 'RF' # 'RF','svm','xgboost','DNN','CNN'
# # # com_num = 16  # number of combination methods
# # # type = 'RNA-pro'  # RNAonly, RNA-RNA, RNA-pro, RNA-compound
# # # evaluation_method(datapath,label_path,resultpath,type = type,com_num = com_num,modelnm = modelnm)
# evaluation_interaction(datapath,label_path,resultpath,type = 'RNA-compound',com_num = 16,modelnm = 'svm')
