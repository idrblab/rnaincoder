from tensorflow.keras import Model, Sequential
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.layers import Dense, Input, Dropout
import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from matplotlib import pyplot as plt
from sklearn import manifold
from random import randint
from sklearn.metrics import roc_curve, auc
from pandas import DataFrame
import matplotlib

import seaborn as sns

import matplotlib.pyplot as plt

plt.switch_backend('agg')

from utils1 import plot_clustering_matplotlib
from utils2 import plot_metrics
from utils3 import plot_roc
from utils4 import rnaheatmap1, rnaheatmap2
from utils5 import kmeans_visual


def z_score(x):
    # mu = np.average(x, axis=0)
    # sigma = np.std(x, axis=0)
    # x = (x-mu)/sigma
    min_max_scaler = preprocessing.MinMaxScaler()
    x_nomal = min_max_scaler.fit_transform(x)
    return x_nomal

METRICS = [
    keras.metrics.TruePositives(name='tp'),
    keras.metrics.FalsePositives(name='fp'),
    keras.metrics.TrueNegatives(name='tn'),
    keras.metrics.FalseNegatives(name='fn'),
    keras.metrics.BinaryAccuracy(name='accuracy'),
    keras.metrics.Precision(name='precision'),
    keras.metrics.Recall(name='recall'),
    keras.metrics.AUC(name='auc'),
]

def training(targets_path, save_path, labels_path=None, dropout=0.2, shape1=128, shape2=32, lr=0.001):
    targets = np.load(targets_path)
    # targets = targets_path
    print('type(targets_path)')
    print(type(targets))
    print(targets.shape)
    targets = z_score(targets)
    if labels_path:

        labels = np.load(labels_path).astype(np.float64)
        print(type(labels))
        print(labels.shape)
        rnaheatmap1(targets, labels, save_path, sename='before', method='ward', metric='euclidean')
        plot_clustering_matplotlib(targets, labels, save_path, '/rna_rna_before', ' before Model Learning')
        X_train, X_test, y_train, y_test = train_test_split(targets, labels, random_state=42)
        # dataset = tf.data.Dataset.from_tensor_slices((X_train, y_train)).shuffle(1000).batch(16)
        feature_ori_dim = targets.shape[1]

        import random
        xname = random.sample('zyxwvutsrqponmlkjihgfedcba', 1)[0]
        xname01 = random.sample('zyxwvutsrqponmlkjihgfedcba', 1)[0]
        xname02 = random.sample('zyxwvutsrqponmlkjihgfedcba', 1)[0]
        xnma = xname + xname01 + xname02

        input_ = Input(shape=feature_ori_dim)
        x = Dense(units=shape1, name='dense_1001', activation='relu')(input_)
        x = Dropout(dropout)(x)
        x = Dense(units=shape2, name=xnma, activation='relu')(x)

        output = Dense(1, activation='sigmoid')(x)
        model = Model(input_, output)
        model.summary()

        checkpoint_path = save_path + "/rna_rna_oned.h5"
        callback = tf.keras.callbacks.ModelCheckpoint(checkpoint_path, monitor='val_loss', mode='min',
                                                      save_best_only=True,
                                                      save_weights_only=True,
                                                      verbose=1)
        callback2 = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10)
        model.compile(optimizer=tf.keras.optimizers.Adam(lr=lr), loss='binary_crossentropy', metrics=METRICS)

        history = model.fit(X_train, y_train, epochs=15, validation_data=(X_test, y_test), verbose=1,
                            callbacks=[callback, callback2])

        plot_metrics(history, save_path)

        train_predictions = model.predict(X_train)
        test_predictions = model.predict(X_test)

        plot_roc(y_train, train_predictions, y_test, test_predictions, save_path)

        model.load_weights(checkpoint_path)

        dense1_layer = Model(inputs=input_, outputs=model.get_layer(xnma).output)
        # dense1_layer = Model(inputs=input_, outputs=model.layers[1].output)
        dense1_output = dense1_layer.predict(targets)
        plot_clustering_matplotlib(dense1_output, labels, save_path, '/rna_rna_after', ' after Model Leaning')
        rnaheatmap1(dense1_output, labels, save_path, sename='after', method='ward', metric='euclidean')
    else:
        kmeans_visual(targets, save_path, k=3, title='/rna_rna')
        rnaheatmap2(targets, save_path, method='ward', metric='euclidean', title='/dd_')










