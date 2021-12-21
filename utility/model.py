#!/usr/bin/env Python
# coding=utf-8

from collections import Counter
from sklearn import metrics,preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.model_selection import train_test_split, PredefinedSplit
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from tqdm import tqdm
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
tqdm.pandas(ascii=True)


def make_comcod(com_n = None):
    # make the combination list of coding methods
    methods_1Ds  = ['Open reading frame (1D)', 'Entropy density of transcript (1D)', 'Global descriptor (1D)', 'K-mer (1D)', 'Codon related (1D)', 'Pseudo protein related (1D)', 'Guanine-cytosine related (1D)', 'Nucleotide related (1D)', 'Secondary structure (1D)', 'EIIP based spectrum (1D)']
    if com_n:
        combination_methods = list(itertools.combinations(methods_1Ds, com_n))
    else:
        for num in range(len(methods_1Ds)):
            com_temp = list(itertools.combinations(methods_1Ds, num+1))
            if num == 0:
                combination_methods = com_temp
            else:
                combination_methods = combination_methods + com_temp
                     
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
    combineFeat = np.concatenate((combnpysA,combnpysB),axis = 1)
    return combineFeat

def mknpy_RNA_pro(combin2,datapath):
    # generate the npy data according to the combination lists
    combnpysB = np.load(os.path.join(datapath, 'protein_B.npy'))
    for numc in range(len(combin2)):
        file_tempA = np.load(os.path.join(datapath, combin2[numc] + '_A.npy'))        
        if numc == 0:
            combnpysA = file_tempA
        else:
            combnpysA = np.concatenate((combnpysA, file_tempA), axis=1)
    combineFeat = np.concatenate((combnpysA,combnpysB),axis = 1)
    return combineFeat

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
    combineFeat = np.concatenate((combnpysA,combnpysB),axis = 1)
    return combineFeat

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


#评估模型的分类效果
def calc_metrics(y_label, y_proba,y_predict):
    # print('y_label')
    # print(y_label)
    # print('y_proba')
    # print(y_proba)
    # print('y_predict')
    # print(y_predict)
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

            # print('labels_train.shape:{}'.format(labels_train.shape))
            # print('score_trains.shape:{}'.format(score_trains.shape))
            # print('group_train.shape:{}'.format(group_train.shape))

            Acc_train, Sn_train, Sp_train, Pre_train, MCC_train, AUC_train = calc_metrics(labels_train,score_trains[:,1],group_train)

            MCC_v = metrics.matthews_corrcoef(labels_all, group_all)
            # accuracy = metrics.accuracy_score(labels_all, group_all)
            # classesnms = [str(i) for i in range(0, classes)]
            # report = metrics.classification_report(labels_all, group_all, target_names=classesnms)
            # print(X_test.shape)
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

def get_abbre(com_lists):
    dict_abre = {
        'Open reading frame (1D)': 'ORF',
        'Entropy density of transcript (1D)': 'EDT',
        'Global descriptor (1D)': 'GBD',
        'K-mer (1D)': 'KME',
        'Codon related (1D)': 'CDR',
        'Pseudo protein related (1D)': 'PPR',
        'Guanine-cytosine related (1D)': 'GCR',
        'Nucleotide related (1D)': 'NTR',
        'Secondary structure (1D)': 'SST',
        'EIIP based spectrum (1D)': 'EBS'
    }
    abres = [dict_abre[methnm] for methnm in com_lists]
    return '+'.join(abres)

def evaluation_method(datapath,label_path,resultpath,type = 'RNAonly',com_num = 2,modelnm = 'RF'):
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
    for combin_num in tqdm(range(0,len(combination_methods),1)):
        if type =='RNAonly':
            combine_X = mknpy_RNAonly(combination_methods[combin_num],datapath)
        elif type =='RNA-RNA':
            combine_X = mknpy_RNA_RNA(combination_methods[combin_num],datapath)
        elif type =='RNA-pro':
            combine_X = mknpy_RNA_pro(combination_methods[combin_num],datapath)
        elif type =='RNA-compound':
            combine_X = mknpy_RNA_compound(combination_methods[combin_num],datapath)
        print('The {}/{} combination of coding methods is {method}'.format(combin_num+1,len(combination_methods),method = combination_methods[combin_num]))
        print('The data shape:{}'.format(combine_X.shape))
        # split the training and test data
        train_x, valid_x, train_y, valid_y = train_test_split(combine_X, label, test_size=0.2, random_state=42,stratify=label)
        # normalization data
        scaler = preprocessing.StandardScaler().fit(train_x)
        train_x = scaler.transform(train_x)
        valid_x = scaler.transform(valid_x)

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

        res_eval[combin_num+1] = [combination_methods[combin_num],get_abbre(combination_methods[combin_num]),best_para, Acc, Sn, Sp, Pre, MCC, AUC]
    indexs = ['combination methods','abbreviation','best parameters','Acc','Sn','Sp','Pre','MCC','AUC']
    evaluation_result = pd.DataFrame(res_eval,index = indexs).T
    evaluation_result = evaluation_result.sort_values(by="MCC" , ascending=False)
    evaluation_result.to_csv(resultpath + '/Evaluation_result.csv')
    print(evaluation_result)
    # make plot
    evaluation_result.iloc[:10,:].plot('abbreviation', 'MCC', kind='barh', legend=False, title ="Evaluation", figsize=(15, 10))
    plt.xlabel("MCC", fontsize=12)
    plt.ylabel("Abbreviation", fontsize=12)
    plt.savefig(resultpath + '/Evaluation_result.png', bbox_inches='tight')


# ========================================================================================

# datapath = '../output/RNA_only/encoding_features'
# label_path = '../demo/RNA-only/Homo38_small.csv'
# resultpath = '../output/RNA_only/classification_result'


# modelnm = 'svm' # 'RF','svm','xgboost','DNN','CNN'
# com_num = 2  # number of combination methods
# type = 'RNAonly'  # RNAonly, RNA-RNA, RNA-pro, RNA-compound
# evaluation_method(datapath,label_path,resultpath,type = type,com_num = com_num,modelnm = modelnm)

