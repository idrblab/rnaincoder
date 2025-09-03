# coding=utf-8
import numpy as np
import time
from sklearn.manifold import TSNE
import tensorflow as tf
from tensorflow.keras import layers
# from tensorflow import keras
from sklearn import metrics
from sklearn.model_selection import train_test_split
import math
# import matplotlib.pyplot as plt
# from sklearn.metrics import roc_curve, auc
# from random import randint
# import seaborn as sns
import pandas as pd
# from sklearn.model_selection import StratifiedShuffleSplit
# from scipy import interp
# import os
# import keras
# from tensorflow.keras.layers.core import *
# from tensorflow.keras.models import *
from tensorflow.keras import optimizers
# from sklearn import preprocessing
from sklearn.model_selection import GridSearchCV,StratifiedKFold
from collections import Counter
from tensorflow.keras.layers import Dense, Conv1D, BatchNormalization, MaxPooling1D,Bidirectional,LSTM,GRU
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from tensorflow.keras.layers import Dropout, Flatten, Input
from tensorflow.keras.models import Model
# from keras.regularizers import l2
# from keras.optimizers import Adam, SGD
from tensorflow.keras.layers import concatenate
import random
import xgboost as xgb


class DNNModel(tf.keras.models.Model):
    def __init__(self, classes, shape1, shape2, dropout):
        super(DNNModel, self).__init__()

        seed = 1234
        # fc_layers = 
        # for fc_layer in fc_layers:
        #     model_t = Dense(units=fc_layer,#, input_dim=input_dim,
        #                     **params_dic)(model_t)
        #     model_t = BatchNormalization()(model_t)
        #     model_t = Activation(activation)(model_t)
        #     # model_t = Dropout(dropout)(model_t)
        #     input_dim = fc_layer

        self.d1 = tf.keras.layers.Dense(128, activation='relu')
        # self.normal1 = tf.keras.layers.BatchNormalization()
        self.drop1 = tf.keras.layers.Dropout(dropout, seed = seed)
        self.d2 = tf.keras.layers.Dense(112, activation='relu')
        # self.drop2 = tf.keras.layers.Dropout(dropout, seed = seed)
        # self.d3 = tf.keras.layers.Dense(96, activation='relu')
        # self.drop3 = tf.keras.layers.Dropout(dropout, seed = seed)
        # self.d4 = tf.keras.layers.Dense(80, activation='relu')
        # self.drop4 = tf.keras.layers.Dropout(dropout, seed = seed)
        # self.d5 = tf.keras.layers.Dense(64, activation='relu')
        # self.drop5 = tf.keras.layers.Dropout(dropout, seed = seed)
        # self.d6 = tf.keras.layers.Dense(48, activation='relu')
        # self.drop6 = tf.keras.layers.Dropout(dropout, seed = seed)
        # self.d7 = tf.keras.layers.Dense(32, activation='relu')
        # self.drop7 = tf.keras.layers.Dropout(dropout, seed = seed)
        # self.d8 = tf.keras.layers.Dense(16, activation='relu')
        self.drop8 = tf.keras.layers.Dropout(dropout, seed = seed)
        self.d9 = tf.keras.layers.Dense(classes, activation='softmax')

    def call(self, inputs,training=None):
        # print('inputs.shape')
        # print(inputs.shape)

        x1 = self.d1(inputs)
        x1 = self.drop1(x1,training = training)
        x1 = self.d2(x1)
        # x1 = self.drop2(x1,training = training)
        # x1 = self.d3(x1)
        # x1 = self.drop3(x1,training = training)
        # x1 = self.d4(x1)
        # x1 = self.drop4(x1,training = training)
        # x1 = self.d5(x1)
        # x1 = self.drop5(x1,training = training)
        # x1 = self.d6(x1)
        # x1 = self.drop6(x1,training = training)
        # x1 = self.d7(x1)
        # x1 = self.drop7(x1,training = training)
        # x1 = self.d8(x1)
        x2 = self.drop8(x1,training = training)
        
        out = self.d9(x2)


        return out, x2

class CNN1DModel(tf.keras.models.Model):
    def __init__(self, classes, shape1, shape2, dropout):
        super(CNN1DModel, self).__init__()

        self.conv1 = tf.keras.layers.Conv1D(
            filters=32,  # ��������Ԫ�������ˣ���Ŀ
            kernel_size=3,  # ����Ұ��С
            padding="valid",  # padding���ԣ�vaild �� same��
            activation=tf.nn.relu  # �����
        )
        self.pool1 = tf.keras.layers.MaxPool1D(pool_size=2, strides=2, padding='valid')
        self.normal1 = tf.keras.layers.BatchNormalization()
        self.drop1 = tf.keras.layers.Dropout(dropout)

        self.conv2 = tf.keras.layers.Conv1D(
            filters=32,
            kernel_size=3,
            padding='valid',
            activation=tf.nn.relu
        )
        self.pool2 = tf.keras.layers.MaxPool1D(pool_size=2, strides=2, padding='valid')
        self.normal2 = tf.keras.layers.BatchNormalization()
        self.drop2 = tf.keras.layers.Dropout(dropout)

        self.conv3 = tf.keras.layers.Conv1D(
            filters=64,
            kernel_size=3,
            padding='valid',
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

#评估模型的分类效果
def calc_metrics(y_label, y_proba,y_predict):
    # print('y_label')
    # print(y_label)
    # print('y_predict')
    # print(y_predict)
    con_matrix = metrics.confusion_matrix(y_label, y_predict)
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
    if len(y_proba.shape) == 1:
        AUC = metrics.roc_auc_score(y_label,y_proba)
    elif len(y_proba.shape)==2:
        fpr, tpr, thresholds = metrics.roc_curve(y_label, y_proba[:, 1])
        AUC = metrics.auc(fpr, tpr)
    return Acc, Sn, Sp, Pre, MCC, AUC



def train_model(X_train,y_train,X_test,y_test,modelnm = 'DNN'):
    
    # parameters
    classes = len(np.unique(y_train))
    shape1 = 128
    shape2 = 32
    Epochs = 50

    dropouts = [0.2,0.5,0.7]
    learnrantes = [0.001,0.01,0.1]
    batch_sizes = [32,64,128,256]

    # 超参数选择
    dropouts = [0.2,0.5]
    learnrantes = [0.001, 0.01, 0.1]
    batch_sizes = [32,64,128]
    parameters_all = [[x,y,z] for x in dropouts for y in learnrantes for z in batch_sizes] 
    parameters = random.sample(parameters_all, 10) 
    
    #traing, validation 8:2 分层抽样
    train_x, test_x, train_y, test_y = train_test_split(X_train, y_train, test_size=0.2, random_state=42,stratify=y_train)
    
    res_dict = {}
    for index, parameter in enumerate(parameters):
        print('Current parameter is {}/{} dropout, learning rate, batchsize: {}'.format(index,len(parameters),parameter))
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
            shape1 = int((int((int((train_x.shape[1] - 2) / 2) - 2) / 2) - 2) / 2)
            model = CNN1DModel(classes, shape1, shape2, dropout)


        # define loss and optimizer
        loss_obj = tf.keras.losses.SparseCategoricalCrossentropy()
        optimizer = tf.keras.optimizers.Adam(lr=learnrante)

        # train_acc = tf.keras.metrics.SparseCategoricalAccuracy()
        # test_acc = tf.keras.metrics.SparseCategoricalAccuracy()

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

            Acc_train, Sn_train, Sp_train, Pre_train, MCC_train, AUC_train = calc_metrics(labels_train,score_trains,group_train)

            MCC_v = metrics.matthews_corrcoef(labels_all, group_all)
            # accuracy = metrics.accuracy_score(labels_all, group_all)
            # classesnms = [str(i) for i in range(0, classes)]
            # report = metrics.classification_report(labels_all, group_all, target_names=classesnms)
            print(X_test.shape)
            # Modeal Evaluation on testing dataset
            score, feature_out = model(X_test, training=False)
            # score = model.predict(X_test, batch_size=32, verbose=1)
            group = np.argmax(score.numpy(), axis=1)
            # Acc_e, Sn_e, Sp_e, Pre_e, MCC_e, AUC_e = calc_metrics(y_test,score,group)

            if mcc < MCC_v or mcc == MCC_v:
                mcc = MCC_v
                # Modeal Evaluation for epoches
                # print('score_all')
                # print(score_all)
                
                Acc_v, Sn_v, Sp_v, Pre_v, MCC_v, AUC_v = calc_metrics(labels_all,score_all,group_all)

                last_improve = epoch
                print('last_improve: %s' % last_improve)           
                Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(y_test,score,group)



            if epoch - last_improve >= 10:
                break

            tmp = 'Epoch {:.3f}, traing Acc {:.3f}, validation Acc {:.3f}, Test Acc {:.3f}, train_loss{:.3f}, validation loss{:.3f}, test auc{:.3f}, test mcc{:.3f}'
            print(tmp.format(epoch + 1, Acc_train,
                            Acc_v,
                            Acc, t_loss, loss, AUC,MCC))
            # print(report)      
        print(model.summary())

        res_dict[index] = [parameter, Acc, Sn, Sp, Pre, MCC, AUC]
        print('The parameter is dropout, learning rate, batchsize: {}, and test MCC: {}'.format(parameter,MCC))
        print('Each group parameters selection use time: {} min'.format((time.time() - time0)/60))
    df = pd.DataFrame(res_dict)
    df_list = df.iloc[5,:].tolist()
    maxindex = df_list.index(max(df_list))
    dfres = df.iloc[:,maxindex].tolist()
    print('The best parameter is dropout, learning rate, batchsize: {}'.format(dfres[0]))
    Acc, Sn, Sp, Pre, MCC, AUC = dfres[1],dfres[2],dfres[3],dfres[4],dfres[5],dfres[6]
    return Acc, Sn, Sp, Pre, MCC, AUC,dfres[0]

def xgboost_model(traindata,trainlabel,X_test,y_test):
    # Convert input data from numpy to XGBoost format
    #traing, validation 8:2 分层抽样
    # train_x, valid_x, train_y, valid_y = train_test_split(X_train, y_train, test_size=0.2, random_state=42,stratify=y_train)
    random_state = 43
    c_v = 5
    cv = StratifiedKFold(n_splits=c_v,random_state=random_state,shuffle=True)
    mcc = -1 
    res_eval = {}
    for fold, (train, valid) in enumerate(cv.split(traindata, trainlabel)):    
        print('the %s fold is starting' % fold)
        print('traindata[train].shape')
        print(traindata[train].shape)
        #traing, validation 8:2 分层抽样
        # train_x, valid_x, train_y, valid_y = train_test_split(traindata, trainlabel, test_size=0.2, random_state=42,stratify=trainlabel)
        train_x = traindata[train]
        valid_x = traindata[valid]
        train_y = trainlabel[train]
        valid_y = trainlabel[valid]

        dtrain = xgb.DMatrix(train_x, label=train_y)
        dvalid = xgb.DMatrix(valid_x, label=valid_y)

        dtest = xgb.DMatrix(X_test, label=y_test)
        
        # Specify sufficient boosting iterations to reach a minimum
        num_rounds = [50, 100, 200, 300, 500]
        max_depths = [3, 6, 10]
        learning_rates = [0.05, 0.1, 0.15]

        parameters_all = [[x,y,z] for x in num_rounds for y in max_depths for z in learning_rates] 
        parameters = random.sample(parameters_all, 10) 
        res_dict = {}

        for index, parameter in enumerate(parameters):
            print('Current parameter is {}/{} dropout, learning rate, batchsize: {}'.format(index,len(parameters),parameter))
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
            # plst = param.items()

            # gpu_res = {} # Store accuracy result
            # tmp = time.time()
            # # Train model
            # param['tree_method'] = 'gpu_hist'
            # bst  = xgb.train(param, dtrain, num_round, evals_result=gpu_res)
            # print("GPU Training Time: %s seconds" % (str(time.time() - tmp)))


            # Repeat for CPU algorithm
            tmp = time.time()
            param['tree_method'] = 'hist'
            cpu_res = {}
            bst  = xgb.train(param, dtrain, num_round,  evals_result=cpu_res)
            print("CPU Training Time: %s seconds" % (str(time.time() - tmp)))
            
            score_valid = bst.predict(dvalid)
            group_valid = (score_valid >= 0.5)*1
            Acc_valid, Sn_valid, Sp_valid, Pre_valid, MCC_valid, AUC_valid = calc_metrics(valid_y,score_valid,group_valid)

            score = bst.predict(dtest)
            # print(ans)
            group = (score >= 0.5)*1
            Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(y_test,score,group)

            res_dict[index] = [parameter, Acc, Sn, Sp, Pre, MCC, AUC, MCC_valid]
            print('The parameter is dropout, learning rate, batchsize: {}, and test MCC: {}'.format(parameter,MCC))
            print('Each group parameters selection use time: {} min'.format((time.time() - tmp)/60))
        df = pd.DataFrame(res_dict)
        df_list = df.iloc[7,:].tolist()
        maxindex = df_list.index(max(df_list))
        dfres = df.iloc[:,maxindex].tolist()
        print('The best parameter is dropout, learning rate, batchsize: {}'.format(dfres[0]))
        best_para = dfres[0]
        Acc, Sn, Sp, Pre, MCC, AUC = dfres[1],dfres[2],dfres[3],dfres[4],dfres[5],dfres[6]
        res_eval[fold] = [best_para, Acc, Sn, Sp, Pre, MCC, AUC]
    df_eval = pd.DataFrame(res_eval)
    df_list_all = df_eval.iloc[5,:].tolist()
    maxindex_all = df_list_all.index(max(df_list_all))
    dfpara = df_eval.iloc[:,maxindex_all].tolist()
    print('The best parameter is dropout, learning rate, batchsize: {}'.format(dfpara[0]))
    best_para_all = dfpara[0]
    dfres_all = df_eval.iloc[1:,:].mean(axis=1).tolist()
    print('df_eval')
    print(df_eval)
    print(dfres_all)
    Acc, Sn, Sp, Pre, MCC, AUC = dfres_all[0],dfres_all[1],dfres_all[2],dfres_all[3],dfres_all[4],dfres_all[5]
    
    return Acc, Sn, Sp, Pre, MCC, AUC,best_para_all
    


def train_model_5fold(X_train,y_train,X_test,y_test):
    # parameters
    classes = len(np.unique(y_train))
    shape1 = 64
    shape2 = 32
    dropout = 0.2
    learnrante = 0.01
    batch_size = 64
    Epochs = 50
      
    
    
    Accs = []
    Sns = []
    Sps = []
    Pres = []
    MCCs = []
    AUCs = []
    
    kf = StratifiedKFold(n_splits=5,shuffle=True,random_state=43)
    fold = 1
    for train_index, test_index in kf.split(X_train,y_train):
        
        # make the model
        model = DNNModel(classes, shape1, shape2, dropout)


        # define loss and optimizer
        loss_obj = tf.keras.losses.SparseCategoricalCrossentropy()
        optimizer = tf.keras.optimizers.Adam(lr=learnrante)

        train_acc = tf.keras.metrics.SparseCategoricalAccuracy()
        test_acc = tf.keras.metrics.SparseCategoricalAccuracy()

        # training model
        def train_step(data, labels):
            with tf.GradientTape() as tape:
                N = len(labels)
                logits, x2 = model(data)

                loss = loss_obj(labels, logits)

                loss = tf.reduce_mean(loss)
                loss = loss / N

                # print(loss)
            grads = tape.gradient(loss, model.trainable_variables)
            optimizer.apply_gradients(zip(grads, model.trainable_variables))

            train_acc(labels, logits)
            return loss
        
        # testing model
        def test_step(data, labels):
            N = len(labels)
            logits, x2 = model(data,training = False)
            test_loss = loss_obj(labels, logits)
            test_loss = tf.reduce_mean(test_loss)
            test_loss = test_loss / N
            logits_ = np.argmax(logits.numpy(), axis=1)
            test_a = test_acc(labels, logits)
            return test_loss, logits_, logits,test_a
        
        # data preparing
        print('train_index:', train_index)
        print('test_index', test_index)
        train_x, test_x, train_y, test_y = X_train[train_index],X_train[test_index],y_train[train_index],y_train[test_index]
        train_dataset = Dataset_make(train_x, train_y, batch_size)
        test_dataset = Dataset_make(test_x, test_y, batch_size)

        mcc = -1
        for epoch in range(Epochs):
            train_acc.reset_states()  
            test_acc.reset_states()

            group_all = np.array([], dtype=int)
            labels_all = np.array([], dtype=int)
            score_all = np.zeros(shape=(1,2))

            for images, labels in train_dataset:
                t_loss = train_step(images, labels)

            for images, labels in test_dataset:
                loss, group, score, test_a = test_step(images, labels)
                group_all = np.append(group_all, group)
                labels_all = np.append(labels_all, labels)
                # print('score.shape')
                # print(score.shape)
                # print(score)
                score_all = np.concatenate((score_all, score),axis = 0)
            score_all = np.delete(score_all, 0, 0)

            MCC_v = metrics.matthews_corrcoef(labels_all, group_all)
            accuracy = metrics.accuracy_score(labels_all, group_all)
            classesnms = [str(i) for i in range(0, classes)]
            report = metrics.classification_report(labels_all, group_all, target_names=classesnms)

            if mcc < MCC_v or mcc == MCC_v:
                mcc = MCC_v
               
                Acc_v, Sn_v, Sp_v, Pre_v, MCC_v, AUC_v = calc_metrics(labels_all,score_all,group_all)

                last_improve = epoch
                print('last_improve: %s' % last_improve)           

                # Modeal Evaluation for cross folds
                score, feature_out = model(X_test)
                group = np.argmax(score.numpy(), axis=1)
                Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(y_test,score,group)


            if epoch - last_improve >= 10:
                break

            tmp = 'Epoch {:.3f}, Acc {:.3f}, Test Acc {:.3f}, acc{:.3f}, t_loss{:.3f}, loss{:.3f}, mcc{:.3f}'
            print(tmp.format(epoch + 1,
                            train_acc.result(),
                            test_a, accuracy, t_loss, loss, MCC))
            # print(report)
        fold_info = 'fold {:.3f}, Acc {:.3f}, Presion {:.3f}, AUC{:.3f}, MCC{:.3f}'
        print(fold_info.format(fold, Acc, Pre, AUC, MCC))
        fold += 1
        Accs.append(Acc)
        Sns.append(Sn)
        Sps.append(Sp)
        Pres.append(Pre)
        MCCs.append(MCC)
        AUCs.append(AUC)
    Acc_mean, Sn_mean, Sp_mean, Pre_mean, MCC_mean, AUC_mean = np.mean(Accs), np.mean(Sns), np.mean(Sps), np.mean(Pres), np.mean(MCCs), np.mean(AUCs)

        

    return Acc_mean, Sn_mean, Sp_mean, Pre_mean, MCC_mean, AUC_mean
def train_model_5fold_alone(X_train,y_train):
    # parameters
    classes = len(np.unique(y_train))
    shape1 = 64
    shape2 = 32
    dropout = 0.2
    learnrante = 0.01
    batch_size = 64
    Epochs = 50
      
    
    
    # train_x, test_x, train_y, test_y = train_test_split(X_train, y_train, test_size=0.2, random_state=42,stratify=y_train)
    Accs = []
    Sns = []
    Sps = []
    Pres = []
    MCCs = []
    AUCs = []
    
    kf = StratifiedKFold(n_splits=5,shuffle=True,random_state=43)
    fold = 1
    for train_index, test_index in kf.split(X_train,y_train):
        
        # make the model
        model = DNNModel(classes, shape1, shape2, dropout)


        # define loss and optimizer
        loss_obj = tf.keras.losses.SparseCategoricalCrossentropy()
        optimizer = tf.keras.optimizers.Adam(lr=learnrante)

        train_acc = tf.keras.metrics.SparseCategoricalAccuracy()
        test_acc = tf.keras.metrics.SparseCategoricalAccuracy()

        # training model
        def train_step(data, labels):
            with tf.GradientTape() as tape:
                N = len(labels)
                logits, x2 = model(data)

                loss = loss_obj(labels, logits)

                loss = tf.reduce_mean(loss)
                loss = loss / N

                # print(loss)
            grads = tape.gradient(loss, model.trainable_variables)
            optimizer.apply_gradients(zip(grads, model.trainable_variables))

            train_acc(labels, logits)
            return loss
        
        # testing model
        def test_step(data, labels):
            N = len(labels)
            logits, x2 = model(data,training = False)
            test_loss = loss_obj(labels, logits)
            test_loss = tf.reduce_mean(test_loss)
            test_loss = test_loss / N
            logits_ = np.argmax(logits.numpy(), axis=1)
            test_a = test_acc(labels, logits)
            return test_loss, logits_, logits,test_a
        
        # data preparing
        print('train_index:', train_index)
        print('test_index', test_index)
        train_x, test_x, train_y, test_y = X_train[train_index],X_train[test_index],y_train[train_index],y_train[test_index]
        train_dataset = Dataset_make(train_x, train_y, batch_size)
        test_dataset = Dataset_make(test_x, test_y, batch_size)

        mcc = -1
        for epoch in range(Epochs):
            train_acc.reset_states()  
            test_acc.reset_states()

            group_all = np.array([], dtype=int)
            labels_all = np.array([], dtype=int)
            score_all = np.zeros(shape=(1,2))

            for images, labels in train_dataset:
                t_loss = train_step(images, labels)

            for images, labels in test_dataset:
                loss, group, score, test_a = test_step(images, labels)
                group_all = np.append(group_all, group)
                labels_all = np.append(labels_all, labels)
                # print('score.shape')
                # print(score.shape)
                # print(score)
                score_all = np.concatenate((score_all, score),axis = 0)
            score_all = np.delete(score_all, 0, 0)

            MCC_v = metrics.matthews_corrcoef(labels_all, group_all)
            accuracy = metrics.accuracy_score(labels_all, group_all)
            classesnms = [str(i) for i in range(0, classes)]
            report = metrics.classification_report(labels_all, group_all, target_names=classesnms)

            if mcc < MCC_v or mcc == MCC_v:
                mcc = MCC_v
                # Modeal Evaluation for epoches
                # print('score_all')
                # print(score_all)
                
                Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(labels_all,score_all,group_all)

                last_improve = epoch
                print('last_improve: %s' % last_improve)           

                # Modeal Evaluation for cross folds
                # score, feature_out = model(X_test)
                # group = np.argmax(score.numpy(), axis=1)
                # Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(y_test,score,group)


            if epoch - last_improve >= 10:
                break

            tmp = 'Epoch {:.3f}, Acc {:.3f}, Test Acc {:.3f}, acc{:.3f}, t_loss{:.3f}, loss{:.3f}, mcc{:.3f}'
            print(tmp.format(epoch + 1,
                            train_acc.result(),
                            test_a, accuracy, t_loss, loss, MCC))
            # print(report)
        fold_info = 'fold {:.3f}, Acc {:.3f}, Presion {:.3f}, AUC{:.3f}, MCC{:.3f}'
        print(fold_info.format(fold, Acc, Pre, AUC, MCC))
        fold += 1
        Accs.append(Acc)
        Sns.append(Sn)
        Sps.append(Sp)
        Pres.append(Pre)
        MCCs.append(MCC)
        AUCs.append(AUC)
    Acc_mean, Sn_mean, Sp_mean, Pre_mean, MCC_mean, AUC_mean = np.mean(Accs), np.mean(Sns), np.mean(Sps), np.mean(Pres), np.mean(MCCs), np.mean(AUCs)

        

    return Acc_mean, Sn_mean, Sp_mean, Pre_mean, MCC_mean, AUC_mean

def CNN_RNN(X_t,y_t,X_test,y_test, model_dir):


    X_train, X_valid, y_train, y_valid = train_test_split(X_t, y_t, test_size=0.2, random_state=42,stratify=y_t)
    X_train = tf.expand_dims(X_train, axis=2)
    X_valid = tf.expand_dims(X_valid, axis=2)
    X_test = tf.expand_dims(X_test, axis=2)

    model = Sequential()

    model.add(Conv1D(filters=45, kernel_size=6, strides=1, activation='relu'))
    model.add(MaxPooling1D(pool_size=2))
    model.add(BatchNormalization())
    model.add(Dropout(0.2))

    model.add(Conv1D(filters=64, kernel_size=5, strides=1, activation='relu'))
    model.add(MaxPooling1D(pool_size=2))
    model.add(BatchNormalization())
    model.add(Dropout(0.2))

    # model.add(Bidirectional(LSTM(100, consume_less='gpu')))
    model.add(Bidirectional(GRU(units=64, activation='tanh', recurrent_activation='hard_sigmoid', use_bias=True, kernel_initializer='glorot_uniform',
                  recurrent_initializer='orthogonal', bias_initializer='zeros', kernel_regularizer=None, recurrent_regularizer=None,
                  bias_regularizer=None, activity_regularizer=None, kernel_constraint=None, recurrent_constraint=None, bias_constraint=None,
                  dropout=0, recurrent_dropout=0, implementation=1, return_sequences=False, return_state=False, go_backwards=False,
                  stateful=False, unroll=False, reset_after=False)))
    model.add(Dropout(0.2))

    model.add(Flatten())
    model.add(Dense(128, activation='relu'))
    model.add(Dropout(0.1))
    model.add(Dense(64, activation='relu'))
    model.add(Dropout(0.1))
    model.add(Dense(2, activation='softmax', name='myfeatures'))

    sgd = optimizers.SGD(lr=0.02, decay=1e-6, momentum=0.9, nesterov=True)

    model.compile(loss='mean_squared_error', optimizer=sgd, metrics=['accuracy'])
    

    # model.compile(loss='binary_crossentropy', optimizer='rmsprop', metrics=['accuracy'])
    # print(model.summary())


    checkpointer = ModelCheckpoint(
        filepath=model_dir+'/bestmodel_ACNN_BLSTM.hdf5',
        verbose=1,
        save_best_only=True)
    earlystopper = EarlyStopping(monitor='val_loss', patience=6, verbose=1)

    print('Training model...')
    history = model.fit(X_train, y_train, epochs=20, batch_size=128, shuffle=True,
                        validation_data=(X_valid, y_valid),
                        callbacks=[checkpointer, earlystopper],
                        verbose=1)
    print(model.summary())
    tresults = model.evaluate(X_test, y_test)
    print(tresults)

    y_pred = model.predict(X_test, batch_size=32, verbose=1)
    pred_group = np.argmax(y_pred, axis=1)
    y = y_test
    print ('Calculating AUC...')
    print('y')
    print(y)
    print('y_pred')
    print(y_pred)
    print('pred_group')
    print(pred_group)

    auroc = metrics.roc_auc_score(y, y_pred[:, 1])
    auprc = metrics.average_precision_score(y, y_pred[:, 1])
    print(auroc, auprc)
    Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(y,y_pred,pred_group)
    # intermediate_layer_model = Model(inputs=model.input,outputs=model.get_layer('myfeatures').output)

    return Acc, Sn, Sp, Pre, MCC, AUC

def dense_factory(x, nb_filter, dropout_rate=None):

    x = layers.Dense(nb_filter, activation='relu')(x)
    x = Dropout(dropout_rate)(x)

    return x

def denseblock(x, nb_layers, nb_filter, growth_rate, dropout_rate=None):

    list_feat = [x]
    concat_axis = -1

    for i in range(nb_layers):
        x = dense_factory(x, nb_filter, dropout_rate)
        list_feat.append(x)
        # x = merge._Merge(list_feat, mode='concat', concat_axis=concat_axis)
        x = concatenate(list_feat, axis = concat_axis)
        
        nb_filter += growth_rate
    return x

def DNN_stack(nb_classes, nb_layers,img_dim1, init_form,growth_rate, nb_filter, dropout_rate,dropout_dense,weight_decay):
    """ Build the DenseNet model

    :param nb_classes: int -- number of classes
    :param img_dim: tuple -- (channels, rows, columns)
    :param depth: int -- how many layers
    :param nb_dense_block: int -- number of dense blocks to add to end
    :param growth_rate: int -- number of filters to add
    :param nb_filter: int -- number of filters
    :param dropout_rate: float -- dropout rate
    :param weight_decay: float -- weight decay
    :param nb_layers:int --numbers of layers in a dense block
    :param filter_size_ori: int -- filter size of first conv1d
    :param dropout_dense: float---drop out rate of dense

    :returns: keras model with nb_layers of conv_factory appended
    :rtype: keras model

    """
    # # first input of 33 seq #
    main_input = Input(shape=img_dim1)


    # Add dense blocks
    # for block_idx in range(nb_dense_block - 1):
    x = denseblock(main_input, nb_layers, nb_filter, growth_rate, dropout_rate=dropout_rate)

  
    x = Dense(32,name ='Dense_2', activation='relu')(x)
    
    x = Dropout(dropout_dense)(x)
    
    x = Dense(nb_classes,name = 'Dense_softmax', activation='softmax')(x)

    phosidnseq = Model(input=[main_input], output=[x], name="PhosIDNSeq")
    #feauture_dense = Model(input=[main_input, input2, input3], output=[x], name="multi-DenseNet")

    return phosidnseq


def model_denseblock(X_t, y_t,X_test, y_test,model_dir,nb_layers = 5):

    X_train, X_valid, y_train, y_valid = train_test_split(X_t, y_t, test_size=0.2, random_state=42,stratify=y_t)
    X_train = tf.expand_dims(X_train, axis=2)
    X_valid = tf.expand_dims(X_valid, axis=2)

    nb_classes = 2
    img_dim1 = X_train.shape[1:]


    ##########parameters#########

    nb_batch_size = 512

    init_form = 'RandomUniform'   #'glorot_normal' #'
    learning_rate = 0.0003
    nb_dense_block = 1
    # nb_layers = 5
    nb_filter = 32
    growth_rate = 24
   
    filter_size_block = 7
    filter_size_ori=1
    dense_number = 32
    self_number = 128
    dropout_rate =0.5
    dropout_dense = 0.3
    weight_decay=0.0001
    nb_epoch = 20



    ###################
    # Construct model #
    ###################
    model = DNN_stack(nb_classes, nb_layers,img_dim1, init_form,growth_rate, nb_filter, dropout_rate,dropout_dense,weight_decay)


    # choose optimazation
    opt = optimizers.Adam(lr=learning_rate, beta_1=0.9, beta_2=0.999, epsilon=1e-08)

    # model compile
    model.compile(loss='binary_crossentropy',
                    optimizer=opt,
                    metrics=['accuracy'])

    # # load weights#

    checkpointer = ModelCheckpoint(
        filepath=model_dir+'/bestmodel_DNN.hdf5',
        verbose=1,
        save_best_only=True)
    earlystopper = EarlyStopping(monitor='val_loss', patience=6, verbose=1)
    # if nb_epoch > 0 :
    model.fit(X_train, y_train, batch_size=nb_batch_size,
                        # validation_data=([X_val1, X_val2, X_val3, y_val),
                        # validation_split=0.1,
                        epochs= nb_epoch, shuffle=True, verbose=1)
    print('Training model...')


    history = model.fit(X_train, y_train, epochs=20, batch_size=128, shuffle=True,
                        validation_data=(X_valid, y_valid),
                        callbacks=[checkpointer, earlystopper],
                        verbose=1)
    print(model.summary())
    tresults = model.evaluate(X_test, y_test)
    print(tresults)

    y_pred = model.predict(X_test, batch_size=32, verbose=1)
    y = y_test
    print ('Calculating AUC...')
    auroc = metrics.roc_auc_score(y, y_pred)
    auprc = metrics.average_precision_score(y, y_pred)
    print(auroc, auprc)
    # intermediate_layer_model = Model(inputs=model.input,outputs=model.get_layer('Dense_softmax').output)

    return auprc, auroc, auprc, auroc, auprc, auroc


def DNN_stackblock(X_t,y_t,X_test,y_test, model_dir,nb_layers):
    #parameters
    nb_filter = 128
    growth_rate = 16
    dropout_rate = 0.5
    nb_epoch = 50
    learnrante = 0.0001


    X_train, X_valid, y_train, y_valid = train_test_split(X_t, y_t, test_size=0.2, random_state=42,stratify=y_t)
    # X_train = tf.expand_dims(X_train, axis=2)
    # X_valid = tf.expand_dims(X_valid, axis=2)
    # X_test = tf.expand_dims(X_test, axis=2)

    model = Sequential()

    nb_filters = [128,64,42,32,16]
    for i in range(nb_layers):
        model.add(Dense(nb_filters[i], activation='relu'))
        # model.add(BatchNormalization())
        model.add(Dropout(dropout_rate))
        nb_filter -= growth_rate
    
    model.add(Dense(2, activation='softmax', name='myfeatures'))

    # sgd = optimizers.SGD(lr=0.02, decay=1e-4, momentum=0.9, nesterov=True)
    optimizer = tf.keras.optimizers.Adam(lr=learnrante)

    # model.compile(loss='mean_squared_error', optimizer=optimizer, metrics=['accuracy'])
    loss_obj = tf.keras.losses.SparseCategoricalCrossentropy()

    # model.compile(loss='binary_crossentropy', optimizer=optimizer, metrics=['accuracy'])
    model.compile(loss=loss_obj, optimizer=optimizer, metrics=['accuracy'])
    # print(model.summary())


    checkpointer = ModelCheckpoint(
        filepath=model_dir+'/bestmodel_ACNN_BLSTM.hdf5',
        verbose=1,
        save_best_only=True)
    earlystopper = EarlyStopping(monitor='val_loss', patience=6, verbose=1)

    print('Training model...')
    history = model.fit(X_train, y_train, epochs=nb_epoch, batch_size=128, shuffle=True,
                        validation_data=(X_valid, y_valid),
                        callbacks=[checkpointer, earlystopper],
                        verbose=1)
    print(model.summary())
    tresults = model.evaluate(X_test, y_test)
    print(tresults)

    y_pred = model.predict(X_test, batch_size=32, verbose=1)
    pred_group = np.argmax(y_pred, axis=1)
    y = y_test
    print ('Calculating AUC...')
    print('y')
    print(y)
    print('y_pred')
    print(y_pred)
    print('pred_group')
    print(pred_group)

    auroc = metrics.roc_auc_score(y, y_pred[:, 1])
    auprc = metrics.average_precision_score(y, y_pred[:, 1])
    print(auroc, auprc)
    Acc, Sn, Sp, Pre, MCC, AUC = calc_metrics(y,y_pred,pred_group)

    return Acc, Sn, Sp, Pre, MCC, AUC