import seaborn as sns
import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as plt
import os
import sys
sys.setrecursionlimit(10000)
plt.switch_backend('agg')
def rnaheatmap1(rnadfdata,Label,png_path, sename = 'before',method ='ward',metric='euclidean'):

    sns.set(font_scale=2)
    sns.set_style('white')

    Label = Label.reshape(len(Label), 1)
    df = DataFrame(np.hstack((rnadfdata, Label)),index=[('Pair' + str(i)) for i in range(0, rnadfdata.shape[0])], columns=[('Feature' + str(i)) for i in range(0, rnadfdata.shape[1])] + ['Label'])

    Colors=['#49759c','#a2cffe','#448ee4','#8ab8fe','#CEFFCE','#28FF28','#007500','#FFFF93','#8C8C00','#FFB5B5','#FF0000','#CE0000','#750000']
    len_labelj = len(df['Label'].unique())
    row_c = dict(zip(df['Label'].unique(), list(Colors[0:len_labelj])))
    # sns.set(font_scale=3.0)
    # plt.rc('font', family='Times New Roman', size=50)
    cm = sns.clustermap(df.drop(columns=['Label']),method = method,metric= metric,
                  cmap=plt.get_cmap('Blues'), row_colors=df['Label'].map(row_c),figsize=(10, 10), tree_kws={'linewidths':2})
    cm.fig.suptitle('Cluster Map ' + sename + ' Model Training', x=0.5,y=1.02, fontsize=25)
    cm.savefig(png_path + '/' + sename + '_rna_rna_heatmap.png', dpi=300)
    plt.close()
    
def rnaheatmap2(rnadfdata,png_path, title,method ='ward',metric='euclidean'):


    print('rnadfdata.shape')
    print(rnadfdata.shape)
    sns.set(font_scale=2)
    sns.set_style('white')
    #去除全为0的列
    mask = (rnadfdata == 0).all(0)
    column_indices = np.where(mask)[0]
    rnadfdata = rnadfdata[:, ~mask]

    print('rnadfdata01.shape')
    print(rnadfdata.shape)
    if rnadfdata.shape[1] >10000:
        rnadfdata = rnadfdata[:,:10000]
        # rnadfdata = rnadfdata

    df = DataFrame(rnadfdata,index=[('Pair' + str(i)) for i in range(0, rnadfdata.shape[0])], columns=[('Feature' + str(i)) for i in range(0, rnadfdata.shape[1])])

    Colors=['#49759c','#a2cffe','#448ee4','#8ab8fe','#CEFFCE','#28FF28','#007500','#FFFF93','#8C8C00','#FFB5B5','#FF0000','#CE0000','#750000']

    cm = sns.clustermap(df,method = method, metric= metric,
                  cmap=plt.get_cmap('Blues'),figsize=(10, 10), tree_kws={'linewidths':2})

    cm.fig.suptitle('Cluster Map', x=0.5,y=1.02, fontsize=25)
    cm.savefig(png_path + title + '_heatmap.png', dpi=300)
    plt.close()