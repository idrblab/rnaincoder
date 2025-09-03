# coding=utf-8
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
import tensorflow as tf
from tensorflow.keras import layers
from tensorflow import keras
from sklearn import metrics
from sklearn.model_selection import train_test_split
import math
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, roc_auc_score
from random import randint
import seaborn as sns
from pandas import DataFrame
from utils5 import kmeans_visual
from utils4 import rnaheatmap2
from sklearn.model_selection import StratifiedShuffleSplit
from scipy import interp
import os
import keras
from keras.layers import Embedding, Bidirectional
from keras.layers.core import *
from keras.models import *
from keras.layers.recurrent import LSTM
from sklearn.model_selection import train_test_split,StratifiedKFold
from collections import Counter


# from tornado.concurrent import run_on_executor

plt.switch_backend('agg')
tf.keras.backend.set_floatx('float64')
gpus = tf.config.experimental.list_physical_devices(device_type='GPU')
tf.config.experimental.set_visible_devices(devices=gpus[1:6], device_type='GPU')
from sklearn import preprocessing


# from multilabel_classification import plot_confusion_matrix


def plot_confusion_matrix(cm, classes, path, n_name,
                          normalize=False,
                          title=None,
                          cmap=plt.cm.Blues):
    if not title:
        if normalize:
            title = 'Normalized confusion matrix'
        else:
            title = 'Confusion Matrix'

    # Compute confusion matrix
    classes = ['class' + i for i in classes]
    # Only use the labels that appear in the data
    #     classes = classes[unique_labels(y_true, y_pred)]
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        # print("Normalized confusion matrix")
    else:
        print('Confusion Matrix')

    # print(cm)
    plt.rcParams.update({"font.size": 10})

    fig, ax = plt.subplots()
    fig.set_figheight(10)
    fig.set_figwidth(10)
    ax.set_title(title, fontsize=25)
    ax.set_xlabel('Predicted label', fontsize=20)
    ax.set_ylabel('True label', fontsize=20)
    # ax.set_xticklabels(classes, fontsize=15)
    # ax.set_xticklabels(classes, fontsize=15)

    im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
    cb = ax.figure.colorbar(im, ax=ax)
    cb.ax.tick_params(labelsize=15)
    # We want to show all ticks...
    ax.set(xticks=np.arange(cm.shape[1]),
           yticks=np.arange(cm.shape[0]), xticklabels=classes, yticklabels=classes)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    fmt = '.2f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(j, i, format(cm[i, j], fmt),
                    ha="center", va="center",
                    color="white" if cm[i, j] > thresh else "black", fontsize=15)
    fig.tight_layout()
    plt.savefig(path + '/' + str(n_name) + '_Confusion Matrix.png')
    plt.close()


def to_categorical(y, num_classes=None, dtype='float64'):
    y = np.array(y, dtype='int')

    input_shape = y.shape
    if input_shape and input_shape[-1] == 1 and len(input_shape) > 1:
        input_shape = tuple(input_shape[:-1])

    y = y.ravel()
    if not num_classes:
        num_classes = np.max(y) + 1
    n = y.shape[0]

    categorical = np.zeros((n, num_classes), dtype=dtype)

    categorical[np.arange(n), y] = 1

    output_shape = input_shape + (num_classes,)
    categorical = np.reshape(categorical, output_shape)
    return categorical


def z_score(x):
    # mu = np.average(x, axis=0)
    # sigma = np.std(x, axis=0)
    # x = (x-mu)/sigma
    min_max_scaler = preprocessing.MinMaxScaler()
    x_nomal = min_max_scaler.fit_transform(x)
    return x_nomal


def rnaheatmap1(rnadfdata, Label, png_path, sename='before', title = 'before',method='ward', metric='euclidean'):
    sns.set(font_scale=2)
    sns.set_style('white')

    Label = Label.reshape(len(Label), 1)
    df = DataFrame(np.hstack((rnadfdata, Label)), index=[('Pair' + str(i)) for i in range(0, rnadfdata.shape[0])],
                   columns=[('Feature' + str(i)) for i in range(0, rnadfdata.shape[1])] + ['Label'])

    Colors = ['#49759c', '#a2cffe', '#448ee4', '#8ab8fe', '#CEFFCE', '#28FF28', '#007500', '#FFFF93', '#8C8C00',
              '#FFB5B5', '#FF0000', '#CE0000', '#750000']
    len_labelj = len(df['Label'].unique())
    row_c = dict(zip(df['Label'].unique(), list(Colors[0:len_labelj])))
    # sns.set(font_scale=3.0)
    # plt.rc('font', family='Times New Roman', size=50)
    cm = sns.clustermap(df.drop(columns=['Label']), method=method, metric=metric,
                        cmap=plt.get_cmap('Blues'), row_colors=df['Label'].map(row_c), figsize=(10, 10),
                        tree_kws={'linewidths': 2})
    cm.fig.suptitle('Cluster Map ' + title + ' Model Training', x=0.5, y=1.02, fontsize=25)
    cm.savefig(png_path + '/' + sename + '_heatmap.png', dpi=300)
    plt.close()


def plot_clustering_matplotlib(z_run, labels, path, folder_name, title1):
    import random
    random.seed(0)
    # labels = labels[:z_run.shape[0]] # because of weird batch_size

    # hex_colors = []
    # for _ in np.unique(labels):
    #     hex_colors.append('#%06X' % randint(0, 0xFFFFFF))
    if len(np.unique(labels)) <= 11:

        hex_colors = ['#437A8B', '#C23147', '#5F86CC', '#F09150', '#AA65C7', '#14BA5E', '#8490BC', '#474EE2', '#904D0C',
                      '#478CC2', '#BEDFB8']
    else:
        hex_colors = []
        for _ in np.unique(labels):
            hex_colors.append('#%06X' % randint(0, 0xFFFFFF))
        print(hex_colors)

    # colors = [hex_colors[int(i)] for i in labels]

    z_run_tsne = TSNE(perplexity=80, min_grad_norm=1E-12, n_iter=250).fit_transform(z_run)
    all_x = [[] for i in range(0, len(np.unique(labels)))]
    all_y = [[] for i in range(0, len(np.unique(labels)))]
    for i in range(z_run_tsne.shape[0]):
        for j in range(0, len(np.unique(labels))):
            if labels[i] == j:
                all_x[j].append(z_run_tsne[i][0])
                all_y[j].append(z_run_tsne[i][1])

    plt.figure(figsize=(10, 10))
    plt.rcParams['savefig.dpi'] = 300  # Í¼Æ¬ÏñËØ
    plt.rcParams['figure.dpi'] = 300  # ·Ö±æÂÊ
    plt.rcParams.update({'font.size': 10})
    bwith = 2.0
    TK = plt.gca()
    TK.spines['bottom'].set_linewidth(bwith)
    TK.spines['left'].set_linewidth(bwith)
    TK.spines['top'].set_linewidth(bwith)
    TK.spines['right'].set_linewidth(bwith)
    types_all = []

    for i in range(0, len(np.unique(labels))):
        types = plt.scatter(all_x[i], all_y[i], c=hex_colors[i], s=20, marker='.')
        types_all.append(types)

    # plt.title('tSNE')

    re_labels = [str(i) for i in range(0, len(types_all
                                              ))]

    plt.legend(types_all, re_labels, loc="upper left")
    # plt.legend(frameon=False,loc="lower right",fontsize='large') #ÉèÖÃÍ¼ÀýÎÞ±ß¿ò£¬½«Í¼Àý·ÅÔÚ×óÉÏ½Ç

    plt.xticks()
    plt.yticks()

    plt.title('TSNE' + title1)

    plt.savefig(path + '/' + folder_name + "tsne.png")
    plt.close()


####################»­Í¼

# return ax
def calculate_multilabel_auc(true_label_, pre_label, logits_, class_list, path, n_name,
                             title='ROC Curve of Deep Neural Network'):
    colors = ['#437A8B', '#C23147', '#5F86CC', '#F09150', '#AA65C7', '#E68223', '#D52685', '#EF7670', '#00A4C5',
              '#9184C1', '#FF9900', '#BEDFB8', '#60C1BD', '#00704A', '#CEFFCE', '#28FF28', '#007500', '#FFFF93',
              '#8C8C00', '#FFB5B5']
    all_true = [[] for i in range(0, len(class_list))]
    all_pre = [[] for i in range(0, len(class_list))]
    true_label = to_categorical(true_label_)
    pre_label = np.array(pre_label)

    plt.figure(figsize=(10, 10))
    # plt.rcParams['figure.figsize'] = (10, 10)  # Í¼ÐÎ´óÐ¡
    plt.rcParams['savefig.dpi'] = 300  # Í¼Æ¬ÏñËØ
    plt.rcParams['figure.dpi'] = 300  # ·Ö±æÂÊ
    plt.rcParams.update({'font.size': 10})

    # ÉèÖÃÍ¼¿òÏß´ÖÏ¸
    bwith = 2.0  # ±ß¿ò¿í¶ÈÉèÖÃÎª2
    TK = plt.gca()  # »ñÈ¡±ß¿ò
    TK.spines['bottom'].set_linewidth(bwith)  # Í¼¿òÏÂ±ß
    TK.spines['left'].set_linewidth(bwith)  # Í¼¿ò×ó±ß
    TK.spines['top'].set_linewidth(bwith)  # Í¼¿òÉÏ±ß
    TK.spines['right'].set_linewidth(bwith)  # Í¼¿òÓÒ±ß
    lw = 1

    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    n_classes = len(class_list)
    for i in range(0, len(class_list)):
        fpr[i], tpr[i], _ = roc_curve(true_label[:, i], pre_label[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])
        plt.plot(fpr[i], tpr[i], color=colors[i], marker='.', markersize=2, lw=lw, linestyle='dashed',
                 # plt.plot(fpr[i], tpr[i], color=colors[i], marker='.', markersize=2, lw=lw,
                 label=class_list[i] + ' AUC = %0.3f' % (roc_auc[i]))

    # # micro����������
    fpr["micro"], tpr["micro"], _ = roc_curve(true_label.ravel(), pre_label.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
    plt.plot(fpr["micro"], tpr["micro"], color='Red', marker='.', markersize=2, lw=lw, linestyle=':',
             label='micro AUC = %0.3f' % (roc_auc["micro"]))
    print('roc_auc["micro"]')
    print(roc_auc["micro"])
    # macro������һ��
    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))
    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(n_classes):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])
    # Finally average it and compute AUC
    mean_tpr /= n_classes
    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])
    plt.plot(fpr["macro"], tpr["macro"], color='green', marker='.', markersize=2, lw=lw, linestyle=':',
             label='macro AUC = %0.3f' % (roc_auc["macro"]))

    print('roc_auc["macro"]')
    print(roc_auc["macro"])

    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    plt.xlabel('False Positive Rate', fontsize=25)
    plt.ylabel('True Positive Rate', fontsize=25)
    # plt.title('Receiver Operating Characteristic', fontsize=25, pad=10)
    plt.title(title, fontsize=25, pad=10)
    # plt.legend(loc="lower right")
    plt.legend(frameon=False, loc="lower right", fontsize='large')

    plt.savefig(path + '/' + str(n_name) + '.png')
    plt.close()

    # Modeal Evaluation
    y_test = list(true_label_)
    predict_ = list(logits_)
    ACC = metrics.accuracy_score(y_test, predict_)
    MCC = metrics.matthews_corrcoef(y_test, predict_)
    precision = metrics.precision_score(y_test, predict_, average='macro')
    f1 = metrics.f1_score(y_test, predict_, average='weighted')
    recall = metrics.recall_score(y_test, predict_, average='micro')
    fbeta = metrics.fbeta_score(y_test, predict_, average='macro', beta=0.5)

    return roc_auc["micro"], ACC, MCC, precision, f1, recall, fbeta

def evaluate_Binarycase(labels_all,predict_all):
    accuracy = metrics.accuracy_score(labels_all, predict_all)
    MCC = metrics.matthews_corrcoef(labels_all, predict_all)
    precision = metrics.precision_score(labels_all, predict_all)
    f1 = metrics.f1_score(labels_all, predict_all)
    recall = metrics.recall_score(labels_all, predict_all)
    fbeta = metrics.fbeta_score(labels_all, predict_all,beta=0.5)
    return accuracy, MCC, precision, f1, recall, fbeta

def calculate_twolabel_auc(true_label_, pre_label, logits_, class_list, path, n_name,
                             title='ROC Curve of Deep Neural Network'):
    colors = ['#437A8B', '#C23147', '#5F86CC', '#F09150', '#AA65C7', '#E68223', '#D52685', '#EF7670', '#00A4C5',
              '#9184C1', '#FF9900', '#BEDFB8', '#60C1BD', '#00704A', '#CEFFCE', '#28FF28', '#007500', '#FFFF93',
              '#8C8C00', '#FFB5B5']
    all_true = [[] for i in range(0, len(class_list))]
    all_pre = [[] for i in range(0, len(class_list))]
    true_label = to_categorical(true_label_)
    pre_label = np.array(pre_label)

    plt.figure(figsize=(10, 10))
    # plt.rcParams['figure.figsize'] = (10, 10)  # Í¼ÐÎ´óÐ¡
    plt.rcParams['savefig.dpi'] = 300  # Í¼Æ¬ÏñËØ
    plt.rcParams['figure.dpi'] = 300  # ·Ö±æÂÊ
    plt.rcParams.update({'font.size': 10})

    # ÉèÖÃÍ¼¿òÏß´ÖÏ¸
    bwith = 2.0  # ±ß¿ò¿í¶ÈÉèÖÃÎª2
    TK = plt.gca()  # »ñÈ¡±ß¿ò
    TK.spines['bottom'].set_linewidth(bwith)  # Í¼¿òÏÂ±ß
    TK.spines['left'].set_linewidth(bwith)  # Í¼¿ò×ó±ß
    TK.spines['top'].set_linewidth(bwith)  # Í¼¿òÉÏ±ß
    TK.spines['right'].set_linewidth(bwith)  # Í¼¿òÓÒ±ß


    auc = roc_auc_score(true_label,pre_label )
    print('two class auc')
    print(auc)

    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    plt.xlabel('False Positive Rate', fontsize=25)
    plt.ylabel('True Positive Rate', fontsize=25)
    # plt.title('Receiver Operating Characteristic', fontsize=25, pad=10)
    plt.title(title, fontsize=25, pad=10)
    # plt.legend(loc="lower right")
    plt.legend(frameon=False, loc="lower right", fontsize='large')

    plt.savefig(path + '/' + str(n_name) + '.png')
    plt.close()

    # Modeal Evaluation
    y_test = list(true_label_)
    predict_ = list(logits_)
    
    ACC, MCC, precision, f1, recall, fbeta = evaluate_Binarycase(y_test,predict_)

    return auc, ACC, MCC, precision, f1, recall, fbeta




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

        return batch_x, batch_y


class DNNModel(tf.keras.models.Model):
    def __init__(self, classes, shape1, shape2, dropout):
        super(DNNModel, self).__init__()

        self.d1 = tf.keras.layers.Dense(shape1, activation='relu')
        self.drop = tf.keras.layers.Dropout(dropout)
        self.d2 = tf.keras.layers.Dense(shape2, activation='relu', name='dense_out')
        self.d3 = tf.keras.layers.Dense(classes, activation='softmax')

    def call(self, inputs):
        # print('inputs.shape')
        # print(inputs.shape)

        x1 = self.d1(inputs)

        # print('x1.shape')
        # print(x1.shape)

        x1 = self.drop(x1)

        # print('x1.shape')
        # print(x1.shape)
        x2 = self.d2(x1)
        # print('x2.shape')
        # print(x2.shape)

        out = self.d3(x2)
        # print('out.shape')
        # print(out.shape)

        return out, x2


class fourierDNN(tf.keras.models.Model):
    def __init__(self, classes, shape1, shape2, dropout):
        super(fourierDNN, self).__init__()

        self.d1 = tf.keras.layers.Dense(shape1, activation='relu')
        self.drop = tf.keras.layers.Dropout(dropout)
        self.d2 = tf.keras.layers.Dense(shape2, activation='relu', name='dense_out')
        self.d3 = tf.keras.layers.Dense(classes, activation='softmax')

    def call(self, inputs):
        # print('inputs.shape')
        # print(inputs.shape)

        x1 = self.d1(inputs)

        out = self.d3(x1)
        # print('out.shape')
        # print(out.shape)

        return out, inputs


# train_CNNModel(train_dataset, test_dataset, Epochs = 30, learnrante = 0.001)
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


class RNNModel(tf.keras.models.Model):

    def __init__(self, classes, units, total_words, embedding_len, max_review_len):
        super(RNNModel, self).__init__()

        self.rnn = keras.Sequential([
            layers.SimpleRNN(units, dropout=0.5, return_sequences=True, unroll=False),
            layers.SimpleRNN(units, dropout=0.5, unroll=False)
        ])

        self.Dense1 = Dense(64, activation='relu')
        self.Dropout1 = Dropout(0.4)
        # fc, [b, 80, 100] => [b, 64] => [b, 1]
        self.outlayer = Dense(classes, activation='softmax')
        # self.outlayer = layers.Dense(1)

    def call(self, inputs, training=True):
        """
        net(x) net(x, training=True) :train mode
        net(x, training=False): test
        :param inputs: [b, 80]
        :param training:
        :return:
        """
        # [b, 80]
        x = tf.expand_dims(inputs, axis=2)
        
        x = self.rnn(x)
        # x = self.lstm(x)
        # print('rnn.shape')
        # print(x.shape)
        x = self.Dense1(x)
        # print('Dense1.shape')
        # print(x.shape)
        x2 = self.Dropout1(x)
        # print('Dropout1.shape')
        # print(x2.shape)

        # out: [b, 64] => [b, 1]
        output = self.outlayer(x2)
        # print('output.shape')
        # print(output.shape)
        # p(y is pos|x)
        # prob = tf.sigmoid(x)

        return output, x2


def train_auto_encoder(X_train, X_test, layers, batch_size=100, nb_epoch=100, activation='sigmoid', optimizer='adam'):
    trained_encoders = []
    trained_decoders = []
    X_train_tmp = np.copy(X_train)
    X_test_tmp = np.copy(X_test)
    for n_in, n_out in zip(layers[:-1], layers[1:]):
        print('Pre-training the layer: Input {} -> Output {}'.format(n_in, n_out))
        ae = Sequential(
            [Dense(n_out, input_dim=X_train_tmp.shape[1], activation=activation, ),
             Dense(n_in, activation=activation),
             Dropout(0.2)]
        )
        ae.compile(loss='mean_squared_error', optimizer=optimizer)
        ae.fit(X_train_tmp, X_train_tmp, batch_size=batch_size, epochs=nb_epoch, verbose=0, shuffle=True)
        # store trained encoder
        trained_encoders.append(ae.layers[0])
        trained_decoders.append(ae.layers[1])
        # update training data
        encoder = Sequential([ae.layers[0]])
        # encoder.evaluate(X_train_tmp, X_train_tmp, batch_size=batch_size)
        X_train_tmp = encoder.predict(X_train_tmp)
        X_test_tmp = encoder.predict(X_test_tmp)

    return trained_encoders, trained_decoders, X_train_tmp, X_test_tmp


def get_auto_encoders(X_train, X_test, batch_size=32):
    encoders_protein, decoders_protein, train_tmp_p, test_tmp_p = train_auto_encoder(
        X_train=X_train,
        X_test=X_test,
        layers=[X_train.shape[1], 256, 128, 64], batch_size=batch_size)
    return encoders_protein


class autoEncoderModel(tf.keras.models.Model):

    def __init__(self, classes, encoders_protein):
        super(autoEncoderModel, self).__init__()
        self.xp_encoded0 = encoders_protein[0]
        self.xp_encoded1 = encoders_protein[1]
        self.xp_encoded2 = encoders_protein[2]

        self.Dropout = Dropout(0.2)
        self.Normal = tf.keras.layers.BatchNormalization()

        self.D1 = Dense(64, kernel_initializer='random_uniform', activation='relu')
        self.D2 = Dense(64, kernel_initializer='random_uniform', activation='relu')
        self.D3 = Dense(classes, activation='softmax')

    def call(self, inputs, training=True):
        inputAs = inputs

        # xp_in_conjoint = Input(shape=(pro_coding_length,))
        xp_encoded = self.xp_encoded0(inputAs)
        xp_encoded = self.Dropout(xp_encoded)
        xp_encoded = self.xp_encoded1(xp_encoded)
        xp_encoded = self.Dropout(xp_encoded)
        xp_encoder = self.xp_encoded2(xp_encoded)
        xp_encoder = self.Dropout(xp_encoder)
        xp_encoder = self.Normal(xp_encoder)
        # xp_encoder = PReLU()(xp_encoder)
        xp_encoder = self.Dropout(xp_encoder)

        # x_out_conjoint = concatenate([xp_encoder, xr_encoder])
        x_out_conjoint = self.D1(xp_encoder)
        x_out_conjoint = self.Normal(x_out_conjoint, training=training)
        x_out_conjoint = self.D2(x_out_conjoint)
        y_conjoint = self.D3(x_out_conjoint)

        return y_conjoint, x_out_conjoint


################
# @run_on_executor
def train_DNNModel00(targetpath, labelpath, Result_path, modelnm, Epochs=50, learnrante=0.001, batch_size=32, shape1=13,
                   shape2=64, dropout=0.5):
    batch_size = batch_size

    # targets = np.load(targetpath).astype(np.float64)
    targets = targetpath

    if len(targets.shape) == 2:
        targets = z_score(targets)
    else:
        targets = targets.reshape(targets.shape[0], -1).astype(np.float64)
        targets = targets[:, :10000]

    # print('targets.shape')
    # print(targets.shape)

    # labels = np.load(labelpath).astype(np.float64)
    labels_original = labelpath
    classes = len(np.unique(labels_original))

    classes_ = [str(i) for i in range(0, classes)]


    ACC_dict = {}
    MCC_dict = {}
    precision_dict = {}
    f1_dict = {}
    recall_dict = {}
    fbeta_dict = {}

    last_improves = {}
    
    cv = StratifiedKFold(n_splits=5,shuffle=True, random_state=42)
    for i, (train, test) in enumerate(cv.split(targets, labels_original)):
        print('cross is %s' % i)
        print('ComResultAB_np[train].shape')
        print(targets[train].shape)
        X_train = targets[train]
        print('X_train.head()')
        print(X_train[1:5,])
        print('labels_original length is %s' % len(labels_original))
        y_train = labels_original[train]
        X_test = targets[test]
        y_test = labels_original[test]
        print('TR0 dataset is %s' % X_train.shape[0] )
        print(X_train.shape)
        print('VAO dataset is %s' % X_test.shape[0] )
        print(X_test.shape)
        print('Counter(y_train)')
        print(Counter(y_train))
        print('Counter(y_test)')
        print(Counter(y_test))

        train_dataset = Dataset_make(X_train, y_train, batch_size)
        test_dataset = Dataset_make(X_test, y_test, batch_size)

        if modelnm == 'DNN':
            model = DNNModel(classes, shape1, shape2, dropout)
        elif modelnm == 'CNN':
            model = CNN1DModel(classes, shape1, shape2, dropout)
        elif modelnm == 'RNN':
            # the most frequest words
            total_words = 10000
            max_review_len = 80
            embedding_len = 1000
            units = 64

            model = RNNModel(classes, units, total_words, embedding_len, max_review_len)
        elif modelnm == 'autoEncoder':
            # generate train and test data
            X_train_conjoint = X_train
            X_test_conjoint = X_test

            # create model
            encoders_pro = get_auto_encoders(X_train_conjoint, X_test_conjoint, batch_size=batch_size)

            model = autoEncoderModel(classes, encoders_pro)
        elif modelnm == 'fourierTransform':
            model = fourierDNN(classes, shape1, shape2, dropout)

        loss_obj = tf.keras.losses.SparseCategoricalCrossentropy()
        optimizer = tf.keras.optimizers.Adam(lr=learnrante)

        train_acc = tf.keras.metrics.SparseCategoricalAccuracy()
        test_acc = tf.keras.metrics.SparseCategoricalAccuracy()

        # @tf.function
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

        # @tf.function
        def test_step(data, labels):
            N = len(labels)
            logits, x2 = model(data)
            test_loss = loss_obj(labels, logits)
            test_loss = tf.reduce_mean(test_loss)
            test_loss = test_loss / N
            logits_ = np.argmax(logits.numpy(), axis=1)
            test_a = test_acc(labels, logits)
            return test_loss, logits_, test_a

        Epochs = Epochs
        acc = 0
        is_early_stoping = 0
        last_improve = 0

        for epoch in range(Epochs):
            train_acc.reset_states() 
            test_acc.reset_states()

            predict_all = np.array([], dtype=int)
            labels_all = np.array([], dtype=int)

            for images, labels in train_dataset:
                t_loss = train_step(images, labels)

            for images, labels_test in test_dataset:
                loss, logits_2, test_a = test_step(images, labels_test)
                predict_all = np.append(predict_all, logits_2)
                labels_all = np.append(labels_all, labels_test)

            accuracy = metrics.accuracy_score(labels_all, predict_all)
            precision = metrics.precision_score(labels_all, predict_all, average='macro')
            recall = metrics.recall_score(labels_all, predict_all, average='macro')

            report = metrics.classification_report(labels_all, predict_all, target_names=classes_)
            confusion = metrics.confusion_matrix(labels_all, predict_all)

            if acc < accuracy or acc == accuracy:
                # tf.saved_model.save(model, os.path.join(Result_path, 'multilabel_DNNmodel.tf'))
                # model.save(os.path.join(Result_path, 'multilabel_bestmodel'), save_format='tf')
                model.save(os.path.join(Result_path, 'multilabel_bestmodel'), save_format='tf')
                # tf.saved_model.save(model, os.path.join(Result_path, 'multilabel_bestmodel'))

                # model.save_weights(os.path.join(Result_path, 'multilabel_DNNmodel.h5'))
                acc = accuracy

                # Modeal Evaluation
                if classes == 2:
                    ACC, MCC, precision, f1, recall, fbeta = evaluate_Binarycase(labels_all, predict_all)
                else:

                    ACC = accuracy
                    MCC = metrics.matthews_corrcoef(labels_all, predict_all)
                    precision = metrics.precision_score(labels_all, predict_all, average='macro')
                    f1 = metrics.f1_score(labels_all, predict_all, average='weighted')
                    recall = metrics.recall_score(labels_all, predict_all, average='micro')
                    fbeta = metrics.fbeta_score(labels_all, predict_all, average='macro', beta=0.5)

                ACC_dict[i] = ACC
                MCC_dict[i] = MCC
                precision_dict[i] = precision
                f1_dict[i] = f1
                recall_dict[i] = recall
                fbeta_dict[i] = fbeta

                last_improve = epoch
                last_improves[i] = last_improve


            if epoch - last_improve >= 10:
                break

            tmp = 'Epoch {:.3f}, Acc {:.3f}, Test Acc {:.3f}, acc{:.3f}, t_loss{:.3f}, loss{:.3f}, precision{:.3f},  recall{:.3f}'
            print(tmp.format(epoch + 1,
                             train_acc.result(),
                             test_a, accuracy, t_loss, loss, precision, recall))
            # print(report)
    print('last_improve')
    print(last_improves)

    ACC = sum(ACC_dict.values())/len(ACC_dict.values())
    MCC = sum(MCC_dict.values())/len(MCC_dict.values())
    precision = sum(precision_dict.values())/len(precision_dict.values())
    f1 = sum(f1_dict.values())/len(f1_dict.values())
    recall = sum(recall_dict.values())/len(recall_dict.values())
    fbeta = sum(fbeta_dict.values())/len(fbeta_dict.values())
    
    return report, ACC, MCC, precision, f1, recall, fbeta




################
# @run_on_executor
def train_DNNModel(targetpath, labelpath, Result_path, modelnm, Epochs=50, learnrante=0.001, batch_size=32, shape1=13,
                   shape2=64, dropout=0.5):
    batch_size = batch_size

    # targets = np.load(targetpath).astype(np.float64)
    targets = targetpath

    if len(targets.shape) == 2:
        targets = z_score(targets)
    else:
        targets = targets.reshape(targets.shape[0], -1).astype(np.float64)
        targets = targets[:, :10000]

    labels = labelpath
    classes = len(np.unique(labels))

    classes_ = [str(i) for i in range(0, classes)]

    if classes > 2:
        split = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
        # print(split)
        # print('labels')
        # print(labels)
        for train_index, test_index in split.split(targets, labels):
            # print("TRAIN:", train_index, "TEST:", test_index)
            X_train, X_test = targets[train_index], targets[test_index]
            y_train, y_test = labels[train_index], labels[test_index]
    else:
        X_train, X_test, y_train, y_test = train_test_split(targets, labels, test_size=0.11, random_state=42,
                                                            stratify=labels)
        # split = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
    print('The traing number is %s' % X_train.shape[0])    
    print('Counter(traing number)')
    print(Counter(y_train))
    print('The validating number is %s' % X_test.shape[0])
    print('Counter(validating number)')
    print(Counter(y_test))
    train_dataset = Dataset_make(X_train, y_train, batch_size)
    test_dataset = Dataset_make(X_test, y_test, batch_size)

    if modelnm == 'DNN':
        model = DNNModel(classes, shape1, shape2, dropout)
    elif modelnm == 'CNN':
        model = CNN1DModel(classes, shape1, shape2, dropout)
    elif modelnm == 'RNN':
        # the most frequest words
        total_words = 10000
        max_review_len = 80
        embedding_len = 1000
        units = 64

        model = RNNModel(classes, units, total_words, embedding_len, max_review_len)
    elif modelnm == 'autoEncoder':
        # generate train and test data
        X_train_conjoint = X_train
        X_test_conjoint = X_test

        # create model
        encoders_pro = get_auto_encoders(X_train_conjoint, X_test_conjoint, batch_size=batch_size)
        model = autoEncoderModel(classes, encoders_pro)
    elif modelnm == 'fourierTransform':
        model = fourierDNN(classes, shape1, shape2, dropout)

    loss_obj = tf.keras.losses.SparseCategoricalCrossentropy()  ########y¡ªtrueÎª±êÁ¿£¬ y_preÎªone hot
    optimizer = tf.keras.optimizers.Adam(lr=learnrante)

    train_acc = tf.keras.metrics.SparseCategoricalAccuracy()
    test_acc = tf.keras.metrics.SparseCategoricalAccuracy()

    # @tf.function
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

    # @tf.function
    def test_step(data, labels):
        N = len(labels)
        logits, x2 = model(data)
        test_loss = loss_obj(labels, logits)
        test_loss = tf.reduce_mean(test_loss)
        test_loss = test_loss / N
        logits_ = np.argmax(logits.numpy(), axis=1)
        test_a = test_acc(labels, logits)
        return test_loss, logits_, test_a

    Epochs = Epochs
    acc = 0
    is_early_stoping = 0
    last_improve = 0

    for epoch in range(Epochs):
        train_acc.reset_states()  ######ÇåÁã
        test_acc.reset_states()

        predict_all = np.array([], dtype=int)
        labels_all = np.array([], dtype=int)

        for images, labels in train_dataset:
            t_loss = train_step(images, labels)

        for images, labels in test_dataset:
            loss, logits_2, test_a = test_step(images, labels)
            predict_all = np.append(predict_all, logits_2)
            labels_all = np.append(labels_all, labels)

        accuracy = metrics.accuracy_score(labels_all, predict_all)
        precision = metrics.precision_score(labels_all, predict_all, average='macro')
        recall = metrics.recall_score(labels_all, predict_all, average='macro')

        report = metrics.classification_report(labels_all, predict_all, target_names=classes_)
        confusion = metrics.confusion_matrix(labels_all, predict_all)

        if acc < accuracy or acc == accuracy:
            # tf.saved_model.save(model, os.path.join(Result_path, 'multilabel_DNNmodel.tf'))
            # model.save(os.path.join(Result_path, 'multilabel_bestmodel'), save_format='tf')
            model.save(os.path.join(Result_path, 'multilabel_bestmodel'), save_format='tf')
            # tf.saved_model.save(model, os.path.join(Result_path, 'multilabel_bestmodel'))

            # model.save_weights(os.path.join(Result_path, 'multilabel_DNNmodel.h5'))
            acc = accuracy

            # Modeal Evaluation
            ACC = accuracy
            MCC = metrics.matthews_corrcoef(labels_all, predict_all)
            precision = metrics.precision_score(labels_all, predict_all, average='macro')
            f1 = metrics.f1_score(labels_all, predict_all, average='weighted')
            recall = metrics.recall_score(labels_all, predict_all, average='micro')
            fbeta = metrics.fbeta_score(labels_all, predict_all, average='macro', beta=0.5)

            last_improve = epoch
            print('last_improve')
            print(last_improve)
        if epoch - last_improve >= 10:
            break

        tmp = 'Epoch {:.3f}, Acc {:.3f}, Test Acc {:.3f}, acc{:.3f}, t_loss{:.3f}, loss{:.3f}, precision{:.3f},  recall{:.3f}'
        print(tmp.format(epoch + 1,
                         train_acc.result(),
                         test_a, accuracy, t_loss, loss, precision, recall))
        # print(report)
    return report, ACC, MCC, precision, f1, recall, fbeta


################
# @run_on_executor
def train_DNNModel_corss(targetpath, labelpath, X_test_cross,y_test_cross,Result_path, modelnm, Epochs=50, learnrante=0.001, batch_size=32, shape1=13,
                   shape2=64, dropout=0.5):
    batch_size = batch_size

    # targets = np.load(targetpath).astype(np.float64)
    targets = targetpath

    if len(targets.shape) == 2:
        targets = z_score(targets)
        X_test_cross = z_score(X_test_cross)

    else:
        targets = targets.reshape(targets.shape[0], -1).astype(np.float64)
        targets = targets[:, :10000]

    # print('targets.shape')
    # print(targets.shape)

    # labels = np.load(labelpath).astype(np.float64)
    labels = labelpath
    classes = len(np.unique(labels))

    classes_ = [str(i) for i in range(0, classes)]

    X_train, X_test, y_train, y_test = targets,X_test_cross,labels,y_test_cross


    print('The traing number is %s' % X_train.shape[0])    
    print('Counter(traing number)')
    print(Counter(y_train))
    print('The validating number is %s' % X_test.shape[0])
    print('Counter(validating number)')
    print(Counter(y_test))
    train_dataset = Dataset_make(X_train, y_train, batch_size)
    test_dataset = Dataset_make(X_test, y_test, batch_size)

    if modelnm == 'DNN':
        model = DNNModel(classes, shape1, shape2, dropout)
    elif modelnm == 'CNN':
        model = CNN1DModel(classes, shape1, shape2, dropout)
    elif modelnm == 'RNN':
        # the most frequest words
        total_words = 10000
        max_review_len = 80
        embedding_len = 1000
        units = 64

        model = RNNModel(classes, units, total_words, embedding_len, max_review_len)
    elif modelnm == 'autoEncoder':
        # generate train and test data
        X_train_conjoint = X_train
        X_test_conjoint = X_test

        # create model
        encoders_pro = get_auto_encoders(X_train_conjoint, X_test_conjoint, batch_size=batch_size)
        model = autoEncoderModel(classes, encoders_pro)
    elif modelnm == 'fourierTransform':
        model = fourierDNN(classes, shape1, shape2, dropout)

    loss_obj = tf.keras.losses.SparseCategoricalCrossentropy()  ########y¡ªtrueÎª±êÁ¿£¬ y_preÎªone hot
    optimizer = tf.keras.optimizers.Adam(lr=learnrante)

    train_acc = tf.keras.metrics.SparseCategoricalAccuracy()
    test_acc = tf.keras.metrics.SparseCategoricalAccuracy()

    # @tf.function
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

    # @tf.function
    def test_step(data, labels):
        N = len(labels)
        logits, x2 = model(data)
        test_loss = loss_obj(labels, logits)
        test_loss = tf.reduce_mean(test_loss)
        test_loss = test_loss / N
        logits_ = np.argmax(logits.numpy(), axis=1)
        test_a = test_acc(labels, logits)
        return test_loss, logits_, test_a

    Epochs = Epochs
    acc = 0
    is_early_stoping = 0
    last_improve = 0

    for epoch in range(Epochs):
        train_acc.reset_states()  ######ÇåÁã
        test_acc.reset_states()

        predict_all = np.array([], dtype=int)
        labels_all = np.array([], dtype=int)

        for images, labels in train_dataset:
            # # print('images.shape')
            # print(images.shape)
            # print(images)
            t_loss = train_step(images, labels)

        for images, labels in test_dataset:
            # print('images.shape')
            # print(images.shape)
            loss, logits_2, test_a = test_step(images, labels)
            predict_all = np.append(predict_all, logits_2)
            labels_all = np.append(labels_all, labels)

        accuracy = metrics.accuracy_score(labels_all, predict_all)
        precision = metrics.precision_score(labels_all, predict_all, average='macro')
        recall = metrics.recall_score(labels_all, predict_all, average='macro')

        report = metrics.classification_report(labels_all, predict_all, target_names=classes_)
        con_matrix = metrics.confusion_matrix(labels_all, predict_all)

        if acc < accuracy or acc == accuracy:
            # tf.saved_model.save(model, os.path.join(Result_path, 'multilabel_DNNmodel.tf'))
            # model.save(os.path.join(Result_path, 'multilabel_bestmodel'), save_format='tf')
            model.save(os.path.join(Result_path, 'multilabel_bestmodel'), save_format='tf')
            # tf.saved_model.save(model, os.path.join(Result_path, 'multilabel_bestmodel'))

            # model.save_weights(os.path.join(Result_path, 'multilabel_DNNmodel.h5'))
            acc = accuracy

            # Modeal Evaluation
            if classes == 2:
                ACC, MCC, precision, f1, recall, fbeta = evaluate_Binarycase(labels_all, predict_all)

                TN = con_matrix[0][0]
                FP = con_matrix[0][1]
                FN = con_matrix[1][0]
                TP = con_matrix[1][1]
                P = TP + FN
                N = TN + FP
                Sn = TP / P if P > 0 else 0
            else:

                ACC = accuracy
                MCC = metrics.matthews_corrcoef(labels_all, predict_all)
                precision = metrics.precision_score(labels_all, predict_all, average='macro')
                f1 = metrics.f1_score(labels_all, predict_all, average='weighted')
                recall = metrics.recall_score(labels_all, predict_all, average='micro')
                fbeta = metrics.fbeta_score(labels_all, predict_all, average='macro', beta=0.5)

            last_improve = epoch
            print('last_improve')
            print(last_improve)
        if epoch - last_improve >= 10:
            break

        tmp = 'Epoch {:.3f}, Acc {:.3f}, Test Acc {:.3f}, acc{:.3f}, t_loss{:.3f}, loss{:.3f}, precision{:.3f},  recall{:.3f}'
        print(tmp.format(epoch + 1,
                         train_acc.result(),
                         test_a, accuracy, t_loss, loss, precision, recall))
        # print(report)
    return report, ACC, MCC, precision, f1, recall, fbeta,Sn

def train_TransferModel(targetpath, labelpath, X_test_cross,y_test_cross,Result_path, model_path, Epochs=50, learnrante=0.001, batch_size=32):
    batch_size = batch_size

    # targets = np.load(targetpath).astype(np.float64)
    targets = targetpath

    if len(targets.shape) == 2:
        targets = z_score(targets)
        X_test_cross = z_score(X_test_cross)
    else:
        targets = targets.reshape(targets.shape[0], -1).astype(np.float64)
        targets = targets[:, :10000]

    # print('targets.shape')
    # print(targets.shape)

    # labels = np.load(labelpath).astype(np.float64)
    labels = labelpath
    classes = len(np.unique(labels))

    classes_ = [str(i) for i in range(0, classes)]

   
    # 5 folds cross validation 
    X_train, X_test, y_train, y_test = targets,X_test_cross,labels,y_test_cross

    print('The traing number is %s' % X_train.shape[0])    
    print('Counter(traing number)')
    print(Counter(y_train))
    print('The validating number is %s' % X_test.shape[0])
    print('Counter(validating number)')
    print(Counter(y_test))

    train_dataset = Dataset_make(X_train, y_train, batch_size)
    test_dataset = Dataset_make(X_test, y_test, batch_size)

    # Load the model for transfer learning
    model = tf.keras.models.load_model(os.path.join(model_path, 'multilabel_bestmodel'))
    model.summary()
    print("Number of layers in the my_model: ", len(model.layers))
    model.trainable = True

    loss_obj = tf.keras.losses.SparseCategoricalCrossentropy()  ########y¡ªtrueÎª±êÁ¿£¬ y_preÎªone hot
    optimizer = tf.keras.optimizers.Adam(lr=learnrante)

    train_acc = tf.keras.metrics.SparseCategoricalAccuracy()
    test_acc = tf.keras.metrics.SparseCategoricalAccuracy()

    # @tf.function
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

    # @tf.function
    def test_step(data, labels):
        N = len(labels)
        logits, x2 = model(data)
        test_loss = loss_obj(labels, logits)
        test_loss = tf.reduce_mean(test_loss)
        test_loss = test_loss / N
        logits_ = np.argmax(logits.numpy(), axis=1)
        test_a = test_acc(labels, logits)
        return test_loss, logits_, test_a

    Epochs = Epochs
    acc = 0
    is_early_stoping = 0
    last_improve = 0

    for epoch in range(Epochs):
        train_acc.reset_states()  ######ÇåÁã
        test_acc.reset_states()

        predict_all = np.array([], dtype=int)
        labels_all = np.array([], dtype=int)

        for images, labels in train_dataset:
            t_loss = train_step(images, labels)

        for images, labels in test_dataset:
            loss, logits_2, test_a = test_step(images, labels)
            predict_all = np.append(predict_all, logits_2)
            labels_all = np.append(labels_all, labels)

        accuracy = metrics.accuracy_score(labels_all, predict_all)
        # precision = metrics.precision_score(labels_all, predict_all, average='macro')
        # recall = metrics.recall_score(labels_all, predict_all, average='macro')

        report = metrics.classification_report(labels_all, predict_all, target_names=classes_)
        con_matrix = metrics.confusion_matrix(labels_all, predict_all)

        #### evaluation -------------
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
        # print('evaluation by myself Acc, Pre and MCC is %s,%s,%s' % (Acc,Pre,MCC))


        if acc < accuracy or acc == accuracy:
            # tf.saved_model.save(model, os.path.join(Result_path, 'multilabel_DNNmodel.tf'))
            # model.save(os.path.join(Result_path, 'multilabel_bestmodel'), save_format='tf')
            model.save(os.path.join(Result_path, 'multilabel_bestmodel'), save_format='tf')
            # tf.saved_model.save(model, os.path.join(Result_path, 'multilabel_bestmodel'))

            # model.save_weights(os.path.join(Result_path, 'multilabel_DNNmodel.h5'))
            acc = accuracy

            # Modeal Evaluation
            if classes == 2:
                ACC, MCC, precision, f1, recall, fbeta = evaluate_Binarycase(labels_all, predict_all)

                TN = con_matrix[0][0]
                FP = con_matrix[0][1]
                FN = con_matrix[1][0]
                TP = con_matrix[1][1]
                P = TP + FN
                N = TN + FP
                Sn = TP / P if P > 0 else 0
                # print('evaluation by sklearn code Acc, Pre and MCC is %s,%s,%s' % (ACC,precision,MCC))
            else:

                ACC = accuracy
                MCC = metrics.matthews_corrcoef(labels_all, predict_all)
                precision = metrics.precision_score(labels_all, predict_all, average='macro')
                f1 = metrics.f1_score(labels_all, predict_all, average='weighted')
                recall = metrics.recall_score(labels_all, predict_all, average='micro')
                fbeta = metrics.fbeta_score(labels_all, predict_all, average='macro', beta=0.5)

            

            last_improve = epoch
            print('last_improve')
            print(last_improve)
        if epoch - last_improve >= 10:
            break

        tmp = 'Epoch {:.3f}, Acc {:.3f}, Test Acc {:.3f}, acc{:.3f}, t_loss{:.3f}, loss{:.3f}, precision{:.3f},  recall{:.3f}'
        print(tmp.format(epoch + 1,
                         train_acc.result(),
                         test_a, accuracy, t_loss, loss, precision, recall))
        # print(report)
    return report, ACC, MCC, precision, f1, recall, fbeta,Sn

def evaluate_DNN(alldata, datapath, png_path, n_name, labelpath=None, label=None, shape1=128, shape2=32, dropout=0.2):
    # data = np.load(datapath).astype(np.float64)
    # print('rnadfdata01.shape')
    # print(rnadfdata.shape)
    if label == 1:
        data = datapath
        if len(data.shape) == 2:
            k = 5
            data = z_score(data)
            alldata01 = z_score(alldata)
        else:
            k = 3
            # ȥ��ȫΪ0����

            data = data.reshape(data.shape[0], -1).astype(np.float64)
            alldata01 = alldata.reshape(alldata.shape[0], -1).astype(np.float64)

            mask = (data == 0).all(0)
            column_indices = np.where(mask)[0]
            data = data[:, ~mask]

            data = data[:, :10000]
            alldata01 = alldata01[:, :10000]
        # label = np.load(labelpath).astype(np.float64)
        label = labelpath
        class_list = len(np.unique(label))
        class_ = [str(i) for i in range(0, class_list)]

        data_feature = data.shape[1]
        # ��������ģ�ͺ����¼���ʱֱ�Ӽ�������ģ��
        my_model = tf.keras.models.load_model(os.path.join(png_path, 'multilabel_bestmodel'))

        print('data.shape')
        print(data.shape)
        l_all, neuralfeatures = my_model(alldata01)
        l, m_f = my_model(data)

        logits_ = np.argmax(l.numpy(), axis=1)
        label = np.squeeze(label)

        con_matrix = metrics.confusion_matrix(label, logits_)
        n_nameauc = n_name
        n_name = n_name[:-4]

        print('lllllll')
        print(l)
        print('logits_')
        print(logits_)

        title = 'ROC Curve of Neural Networks'
        if class_list == 2:
            aucvalue, ACC, MCC, precision, f1, recall, fbeta = calculate_twolabel_auc(label, l, logits_, class_,
                                                                                        png_path,
                                                                                        n_nameauc, title)
            TN = con_matrix[0][0]
            FP = con_matrix[0][1]
            FN = con_matrix[1][0]
            TP = con_matrix[1][1]
            P = TP + FN
            N = TN + FP
            Sn = TP / P if P > 0 else 0

        else:
            aucvalue, ACC, MCC, precision, f1, recall, fbeta = calculate_multilabel_auc(label, l, logits_, class_,
                                                                                        png_path,
                                                                                        n_nameauc, title)


        return neuralfeatures, aucvalue, ACC, MCC, precision, f1, recall, fbeta, Sn
    else:
        if len(alldata.shape) == 2:
            k = 5
            alldata01 = z_score(alldata)
        else:
            k = 3
            # ȥ��ȫΪ0����

            alldata01 = alldata.reshape(alldata.shape[0], -1).astype(np.float64)

            mask = (alldata01 == 0).all(0)
            column_indices = np.where(mask)[0]
            alldata01 = alldata01[:, ~mask]

            alldata01 = alldata01[:, :10000]
        kmeans_visual(alldata01, png_path, k=k, title='/multilabel_label')
        rnaheatmap2(alldata01, png_path, method='ward', metric='euclidean', title=n_name)
