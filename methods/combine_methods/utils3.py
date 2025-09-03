#coding:gbk
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from sklearn import manifold
from sklearn.metrics import roc_curve, auc
import os
def plot_roc(label, predict_score, label1, predict_score1, path):
    fpr, tpr, _ = roc_curve(label, predict_score)
    fpr1, tpr1, _ = roc_curve(label1, predict_score1)
    roc_auc = auc(fpr, tpr)
    roc_auc1 = auc(fpr1, tpr1)
    plt.figure(figsize=(10, 10))
    # plt.rcParams['figure.figsize'] = (10, 10)  # 图形大小
    plt.rcParams['savefig.dpi'] = 300  # 图片像素
    plt.rcParams['figure.dpi'] = 300  # 分辨率
    plt.rcParams.update({'font.size': 20})

    # 设置图框线粗细
    bwith = 2.0  # 边框宽度设置为2
    TK = plt.gca()  # 获取边框
    TK.spines['bottom'].set_linewidth(bwith)  # 图框下边
    TK.spines['left'].set_linewidth(bwith)  # 图框左边
    TK.spines['top'].set_linewidth(bwith)  # 图框上边
    TK.spines['right'].set_linewidth(bwith)  # 图框右边
    lw = 2
    # plt.figure(figsize=(10,10))
    plt.plot(fpr, tpr, color='darkorange', marker='o', markersize=5, lw=lw, linestyle='dashed',
             label='Train ROC Curve (area = %0.3f)' % roc_auc)
    plt.plot(fpr1, tpr1, color='green', marker='o', lw=lw, markersize=5, linestyle='dashed',
             label='Test ROC Curve (area = %0.3f)' % roc_auc1)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)


    plt.xlabel('False Positive Rate', fontsize=25)
    plt.ylabel('True Positive Rate', fontsize=25)
    plt.title('ROC Curve of Deep Neural Network', fontsize=25, pad=10)
    # plt.legend(loc="lower right")
    plt.legend(frameon=False, loc="lower right", fontsize='large')  # 设置图例无边框，将图例放在左上角

    plt.savefig(path + '/rna_rna_AUC.png')
    plt.close()