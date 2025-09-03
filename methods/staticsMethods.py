# -*- coding:utf-8 -*-

from sklearn.manifold import TSNE
from Bio import SeqIO
import Bio.SeqIO as Seq
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
# import pylab as pl
# import shutil
import numpy as np
# from pandas import Series,DataFrame
import seaborn as sns
# import palettable#python颜色库
from sklearn import preprocessing

import csv
from sklearn.decomposition import PCA
import pygal
# from pygal.style import Style
import re
import zipfile
###############方法函数————————————————————————————————————————————————————————————
def fastalen(infasta1,infasta2):
    seqname = []
    for seq in Seq.parse(infasta1, 'fasta'):
        seqid = seq.id
        seqname.append(seqid)
    for seq in Seq.parse(infasta2, 'fasta'):
        seqid = seq.id
        seqname.append(seqid)
    leng = len(seqname)
    return leng

###画序列长度分布图，将两个文件的序列长度合在一起画图----------------------------------------------
def plotfaA(infasta1, infasta2,png_path, num_lis = 20):
    pngfile = os.path.join(png_path, '_fatsaA.svg')
    if os.path.exists(pngfile):  # 如果文件存在
        # 删除文件，可使用以下两种方法。
        os.remove(pngfile)
    seqname = []
    seqlen = []
    for seq in Seq.parse(infasta1, 'fasta'):
        seqid = seq.id
        seqname.append(seqid)
        RNA_size = len(seq.seq)
        seqlen.append(RNA_size)

    for seq in Seq.parse(infasta2, 'fasta'):
        seqid = seq.id
        seqname.append(seqid)
        RNA_size = len(seq.seq)
        seqlen.append(RNA_size)

    qujian = (max(seqlen) - min(seqlen)) / int(num_lis)
    seqlen.sort()
    count = []
    nameall = []
    n = 1
    x = 0
    xn = min(seqlen)
    for i, value in enumerate(seqlen):
        # print(i, value)
        if (value > qujian * n):
            h = i - x
            if x == 0:
                count.append(h + 1)
            else:
                count.append(h)

            name = (str(xn) + '-' + str(value))
            nameall.append(name)
            xn = value

            x = i
            n += 1
            # print(h+1)
    import pandas as pd
    df = pd.DataFrame({
        'variable': nameall,
        'value': count,
    })
    print(df.shape)
    bar_chart = pygal.Bar(height=250,show_legend=False)

    # bar_chart.title = 'test'  # 设置标题
    bar_chart.x_labels = df['variable']
    bar_chart.x_title = "Sequence Length"
    bar_chart.y_title = "Counts"
    bar_chart.add('Distribution', df['value'])

    bar_chart.render_to_file(os.path.join(png_path, '_fatsaA.svg'))
    imgpath = os.path.join(png_path, '_fatsaA.svg')
    return imgpath
###画序列长度分布图，将两个文件的序列长度合在一起画图----------------------------------------------
def gen_dict(fiel):
    name = []
    seq_all = []
    for seq in SeqIO.parse(fiel, 'fasta'):
        seq_name = seq.id
        seq_name = '>' + seq_name
        sequence = str(seq.seq)
        name.append(seq_name)
        seq_all.append(sequence)
    name_seq = dict(zip(name, seq_all))
    return name,name_seq
def plotfa_interact(infasta1, infasta2, png_path, num_lis = 20):
    #删除已有的图片
    pngfile = os.path.join(png_path, '_fatsaB.svg')
    if os.path.exists(pngfile):  # 如果文件存在
        # 删除文件，可使用以下两种方法。
        os.remove(pngfile)
    #将两个fasta文件生成字典
    sename01,sequence01 = gen_dict(infasta1)
    sename02,sequence02 = gen_dict(infasta2)
    #按照顺序一对一将相互作用的两个文件将名字和序列长度记录
    seqlen = []
    for n in range(0,len(sename01),1):
        intername = sename01[n]+sename01[n]
        interlen = len(sequence01[sename01[n]])+len(sequence02[sename02[n]])
        seqlen.append(interlen)


    qujian = (max(seqlen) - min(seqlen)) / int(num_lis)
    seqlen.sort()
    count = []
    nameall = []
    n = 1
    x = 0
    xn = min(seqlen)
    for i, value in enumerate(seqlen):
        # print(i, value)
        if (value > qujian * n):
            h = i - x
            if x == 0:
                count.append(h + 1)
            else:
                count.append(h)

            name = (str(xn) + '-' + str(value))
            nameall.append(name)
            xn = value

            x = i
            n += 1
            # print(h+1)
    import pandas as pd
    df = pd.DataFrame({
        'variable': nameall,
        'value': count,
    })
    # custom_style = Style(colors='#ff0808')
    # red = Style(colors=('#C23147',))
    print(df.shape)
    bar_chart = pygal.Bar(height=250,show_legend=False)

    # bar_chart.title = 'test'  # 设置标题
    bar_chart.x_labels = df['variable']
    bar_chart.x_title = "Sequence Length"
    bar_chart.y_title = "Counts"
    bar_chart.add('Distribution', df['value'])

    bar_chart.render_to_file(os.path.join(png_path, '_fatsaB.svg'))
    imgpath = os.path.join(png_path, '_fatsaB.svg')
    return imgpath

def plotfaB(infasta,numlist =15):
    ################SVG-----------------------------------------------
    import os
    png_path = os.path.join(os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."),
                            'statics/images')
    pngfile = os.path.join(png_path, '_fatsaB.svg')
    if os.path.exists(pngfile):  # 如果文件存在
        # 删除文件，可使用以下两种方法。
        os.remove(pngfile)

    seqname = []
    seqlen = []
    import Bio.SeqIO as Seq
    import os

    for seq in Seq.parse(infasta, 'fasta'):
        seqid = seq.id
        seqname.append(seqid)
        RNA_size = len(seq.seq)
        seqlen.append(RNA_size)

    qujian = (max(seqlen) - min(seqlen)) / int(numlist)
    seqlen.sort()
    count = []
    nameall = []
    n = 1
    x = 0
    xn = min(seqlen)
    for i, value in enumerate(seqlen):
        # print(i, value)
        if (value > qujian * n):
            h = i - x
            if x == 0:
                count.append(h + 1)
            else:
                count.append(h)

            name = (str(xn) + '-' + str(value))
            nameall.append(name)
            xn = value

            x = i
            n += 1
            # print(h+1)
    import pandas as pd
    df = pd.DataFrame({
        'variable': nameall,
        'value': count,
    })
    # red = Style(colors=('#C23147',))
    bar_chart = pygal.Bar(height=250,show_legend=False)
    # bar_chart.title = 'test'  # 设置标题
    bar_chart.x_labels = df['variable']
    bar_chart.x_title = "Sequence Length"
    bar_chart.y_title = "Counts"
    bar_chart.add('Distribution', df['value'])


    png_path = os.path.join(os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."),
                                 'statics/images')
    bar_chart.render_to_file(os.path.join(png_path, '_fatsaB.svg'))
    # plt.savefig(os.path.join(png_path, '_fatsa.jpg'))
    imgpath = os.path.join(png_path, '_fatsaB.svg')
    return imgpath


def plotfaRNA(infasta):

    ################SVG-----------------------------------------------
    import os
    png_path = os.path.join(os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."),
                            'statics/images')
    pngfile = os.path.join(png_path, '_RNAfatsaRNA.svg')
    if os.path.exists(pngfile):  # 如果文件存在
        # 删除文件，可使用以下两种方法。
        os.remove(pngfile)

    seqname = []
    seqlen = []
    import Bio.SeqIO as Seq
    import os

    for seq in Seq.parse(infasta, 'fasta'):
        seqid = seq.id
        seqname.append(seqid)
        RNA_size = len(seq.seq)
        seqlen.append(RNA_size)

    qujian = (max(seqlen) - min(seqlen)) / 10
    seqlen.sort()
    count = []
    nameall = []
    n = 1
    x = 0
    xn = min(seqlen)
    for i, value in enumerate(seqlen):
        # print(i, value)
        if (value > qujian * n):
            h = i - x
            if x == 0:
                count.append(h + 1)
            else:
                count.append(h)

            name = (str(xn) + '-' + str(value))
            nameall.append(name)
            xn = value

            x = i
            n += 1
            # print(h+1)
    import pandas as pd
    df = pd.DataFrame({
        'variable': nameall,
        'value': count,
    })
    # red = Style(colors=('#C23147',))
    bar_chart = pygal.Bar(height=450,show_legend=False)
    # bar_chart.title = 'test'  # 设置标题
    bar_chart.x_labels = df['variable']
    bar_chart.x_title = "Sequence Length"
    bar_chart.y_title = "Counts"
    bar_chart.add('Distribution', df['value'])



    png_path = os.path.join(os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."),
                                 'statics/images')

    bar_chart.render_to_file(os.path.join(png_path, '_RNAfatsaRNA.svg'))
    # plt.savefig(os.path.join(png_path, '_fatsa.jpg'))
    imgpath = os.path.join(png_path, '_RNAfatsaRNA.svg')
    return imgpath

def plotsmall(datadf01,numlist):
    ################SVG-----------------------------------------------
    import os
    png_path = os.path.join(os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."),
                            'statics/images')
    pngfile = os.path.join(png_path, '_small.svg')
    if os.path.exists(pngfile):  # 如果文件存在
        # 删除文件，可使用以下两种方法。
        os.remove(pngfile)

    # datadf = pd.read_csv(pddata)
    # datadf01 = datadf.dropna(axis=1, how='any', thresh=None, subset=None, inplace=False)
    seqlen = datadf01['MW'].tolist()

    qujian = (max(seqlen) - min(seqlen)) / int(numlist)
    seqlen.sort()
    count = []
    nameall = []
    n = 1
    x = 0
    xn = int(min(seqlen))
    for i, value in enumerate(seqlen):
        #     print(i, value)
        value = int(value)
        if (value > qujian * n + min(seqlen)):
            h = i - x
            if x == 0:
                count.append(h + 1)
            else:
                count.append(h)

            name = (str(xn) + '-' + str(value))
            nameall.append(name)
            xn = value

            x = i
            n += 1
    df = pd.DataFrame({
        'variable': nameall,
        'value': count,
    })

    bar_chart = pygal.Bar(height=250,show_legend=False)
    # bar_chart.title = 'test'  # 设置标题
    bar_chart.x_labels = df['variable']
    bar_chart.x_title = "Molecular Weight"
    bar_chart.y_title = "Counts"
    bar_chart.add('Distribution', df['value'])


    png_path = os.path.join(os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."),
                                 'statics/images')

    bar_chart.render_to_file(os.path.join(png_path, '_small.svg'))
    imgpath = os.path.join(png_path, '_small.svg')
    return imgpath

def plotsmall_singel(datadf01, fastafile, png_path,numlist):
    ################SVG-----------------------------------------------
    import os
    
    pngfile = os.path.join(png_path, '_fatsaA.svg')
    if os.path.exists(pngfile):  # 如果文件存在
        # 删除文件，可使用以下两种方法。
        os.remove(pngfile)
    #小分子的质量列表
    seqlen = datadf01['MW'].tolist()
    # RNA的序列长度列表
    seqlenRNA = []
    import Bio.SeqIO as Seq
    import os

    for seq in Seq.parse(fastafile, 'fasta'):
        RNA_size = len(seq.seq)
        seqlenRNA.append(RNA_size)
    seqlen.extend(seqlenRNA)
    qujian = (max(seqlen) - min(seqlen)) / int(numlist)
    seqlen.sort()
    count = []
    nameall = []
    n = 1
    x = 0
    xn = int(min(seqlen))
    for i, value in enumerate(seqlen):
        #     print(i, value)
        value = int(value)
        if (value > qujian * n + min(seqlen)):
            h = i - x
            if x == 0:
                count.append(h + 1)
            else:
                count.append(h)

            name = (str(xn) + '-' + str(value))
            nameall.append(name)
            xn = value

            x = i
            n += 1
    df = pd.DataFrame({
        'variable': nameall,
        'value': count,
    })

    bar_chart = pygal.Bar(height=250,show_legend=False)
    # bar_chart.title = 'test'  # 设置标题
    bar_chart.x_labels = df['variable']
    bar_chart.x_title = "Molecular Weight"
    bar_chart.y_title = "Counts"
    bar_chart.add('Distribution', df['value'])

    bar_chart.render_to_file(os.path.join(png_path, '_fatsaA.svg'))
    imgpath = os.path.join(png_path, '_fatsaA.svg')
    return imgpath


def plotsmall_inter(datadf01, fastafile, filepath3,png_path, numlist):
    ################SVG-----------------------------------------------
    import os

    pngfile = os.path.join(png_path, '_small.svg')
    if os.path.exists(pngfile):  # 如果文件存在
        # 删除文件，可使用以下两种方法。
        os.remove(pngfile)
    # 小分子的名字和分子质量字典
    seqlensma = datadf01['MW'].tolist()
    rowname = datadf01['Name'].tolist()
    dict_small = dict(zip(rowname,seqlensma))
    # print(dict_small)
    # RNA的名字和序列长度字典
    seqlenRNA = []
    seqname = []
    import Bio.SeqIO as Seq
    import os
    for seq in Seq.parse(fastafile, 'fasta'):
        seqname.append(seq.id)
        RNA_size = len(seq.seq)
        seqlenRNA.append(RNA_size)
    dict_RNA = dict(zip(seqname,seqlenRNA))
    # print(dict_RNA)
    # label的序列order
    A_RNAname = []
    B_SMAname = []
    csvfile01 = open(filepath3, 'r')
    for line in csvfile01:
        if '>' in line:
            line = line.strip('\n').split(',')
            A_RNAname.append(line[0][1:])
            B_SMAname.append(line[1])


    #按照顺序一对一将相互作用的两个文件将名字和序列长度记录
    seqlen = []
    for n in range(0,len(A_RNAname),1):
        interlen = dict_RNA[A_RNAname[n]]+dict_small[B_SMAname[n]]
        seqlen.append(interlen)
    #将相互作用的一对一的关系长度用于计算
    qujian = (max(seqlen) - min(seqlen)) / int(numlist)
    seqlen.sort()
    count = []
    nameall = []
    n = 1
    x = 0
    xn = int(min(seqlen))
    for i, value in enumerate(seqlen):
        #     print(i, value)
        value = int(value)
        if (value > qujian * n + min(seqlen)):
            h = i - x
            if x == 0:
                count.append(h + 1)
            else:
                count.append(h)

            name = (str(xn) + '-' + str(value))
            nameall.append(name)
            xn = value

            x = i
            n += 1
    df = pd.DataFrame({
        'variable': nameall,
        'value': count,
    })

    bar_chart = pygal.Bar(height=250, show_legend=False)
    # bar_chart.title = 'test'  # 设置标题
    bar_chart.x_labels = df['variable']
    bar_chart.x_title = "Molecular Weight"
    bar_chart.y_title = "Counts"
    bar_chart.add('Distribution', df['value'])

    bar_chart.render_to_file(os.path.join(png_path, '_small.svg'))
    # plt.savefig(os.path.join(png_path, '_fatsa.jpg'))
    imgpath = os.path.join(png_path, '_small.svg')
    return imgpath
###将用户自定义上传的fasta文件根据label文件分开------------------------------
class rnaProcess_File(object):
    def __init__(self, filepath1, filepath3):
        self.filepath1 = filepath1
        self.filepath3 = filepath3
    def gen_dict(self, fiel):
        name = []
        seq_all = []
        for seq in SeqIO.parse(fiel, 'fasta'):
            seq_name = seq.id
            seq_name = '>' + seq_name
            sequence = str(seq.seq)
            name.append(seq_name)
            seq_all.append(sequence)
        name_seq = dict(zip(name, seq_all))
        return name_seq

    def get_interaction_pair(self):
        A_dict = self.gen_dict(self.filepath1)
        A_save = open(self.filepath1, 'w')

        csvfile01 = open(self.filepath3, 'r')
        for line in csvfile01:
            if '>' in line:
                line = line.strip('\n').split(',')
                A = line[0]
                a_seq = A_dict[A]
                A_save.write(A + '\n' + a_seq + '\n')
        A_save.close()

#####将用户自定义上传的fasta文件根据是否存在除了ACGT的其他字符，将U替换为T然后删除这些非法序列 ------------------------------
class remove_uncorrect(object):
    # 判断一个RNAfasta文件中是否有除了ACGT的其他字符，将U替换为T然后删除这些非法序列
    def __init__(self, filepath1, filepath2):
        self.filepath1 = filepath1
        self.filepath2 = filepath2

    def gen_dict(self, file):
        name = []
        seq_all = []
        for seq in Seq.parse(file, 'fasta'):
            seq_name = seq.id
            seq_name = '>' + seq_name
            sequence = str(seq.seq)
            name.append(seq_name)
            seq_all.append(sequence)
        name_seq = dict(zip(name, seq_all))
        return name_seq

    def get_file(self):
        seqnames = []
        sequences = []

        for seq in Seq.parse(self.filepath1, 'fasta'):
            seqid = seq.id
            seqid = '>' + seqid
            n = 0
            sequence_o = str(seq.seq)
            sequence_r = sequence_o.replace("U", "T")

            all_leters = 'bdefhijklmnopqrsuvwxyz'
            all_leter = []
            for leter in all_leters:
                all_leter.append(leter.upper())
            for all_lete in all_leter:
                print(sequence_r)
                sequence_r = sequence_r.replace(all_lete, "")
            for element in set(sequence_r):
                if element not in ['A', 'C', 'G', 'T']:
                    n = + 1
            if n == 0:
                seqnames.append(seqid)
                sequences.append(sequence_r)
        name_seq = dict(zip(seqnames, sequences))
        # 重新写入新的删除非法字符的fasta文件
        A_save = open(self.filepath1, 'w')
        for seqname in seqnames:
            a_seq = name_seq[seqname]
            A_save.write(str(seqname) + '\n' + str(a_seq) + '\n')
        A_save.close()

        # 根据非法字符所在的序列删除label文件中的序列名字
        csvfile01 = open(self.filepath2, 'r')
        labelnames = []
        labels = []

        for line in csvfile01:
            if '>' in line:
                line = line.strip('\n').split(',')
                for seqname in seqnames:
                    if seqname == line[0]:
                        labelnames.append(line[0])
                        labels.append(line[1])
        csvfile01.close()
        name_label = {
            'Seqname': labelnames,
            'Label': labels
        }
        df_label = pd.DataFrame(name_label)
        df_label.to_csv(self.filepath2, index=False, header=True)
###画有label的序列长度分布图----------------------------------------------
def plotRNA(infasta,csvpath):
    png_path = os.path.join(os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."),
                            'files/download')
    pngfile = os.path.join(png_path, '_RNAfatsa.jpg')
    if os.path.exists(pngfile):  # 如果文件存在
        # 删除文件，可使用以下两种方法。
        os.remove(pngfile)
    csvfile = open(csvpath,'r')
    reader = csv.DictReader(csvfile)
    labelcolumn = [row['Label'] for row in reader]
    labelj = list(set(labelcolumn))

    seqname = []
    seqlen_a = []
    dictdata = {}
    #每一对匹配关系存进字典里
    line1 = []
    line0 = []
    csvfile01 = open(csvpath,'r')
    for line in csvfile01:
        if '>' in line:
            line = line.strip('\n').split(',')
            line1.append(line[1])
            line0.append(line[0])
    dictcsv = dict(zip(line0, line1))
    # print(dictcsv)
    #把序列按照组分为num组
    for num in list(range(0,len(labelj),1)):
        seqlen = []
        for seq in Seq.parse(infasta, 'fasta'):
            seqid = seq.id
            seqidnew = ">" + seqid
            seqname.append(seqid)

            if dictcsv[seqidnew] == labelj[num]:
                RNA_size = len(seq.seq)
                seqlen.append(RNA_size)
                seqlen_a.append(RNA_size)
                dictname = 'group' +str(num)
                dictdata[dictname] =seqlen
#     print(seqlen_a)
    max_size = max(seqlen)
    min_size = min(seqlen)
    max_size1 = max_size + 1
    import seaborn as sns

    plt.figure(figsize=(8,3.5))
    # plt.subplots_adjust(top=0.1, bottom=0.05, right=0.1, left=0.05, hspace=0, wspace=0)
    # plt.margins(1, 1)
    for num in list(range(0,len(labelj),1)):
        dictname = 'group' +str(num)
#         print(dictname)
        dataplot = np.array(dictdata[dictname])
        plt.xlabel('Sequence Length')
        plt.ylabel('Density')
        plt.xticks(np.arange(0, max_size1, 1000))
        sns.kdeplot(dataplot, shade=True,label='Label' +str(num))
        plt.title(r'Sequence Length Distribution')
    png_path = os.path.join(os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."),
                            'files/download')
    plt.savefig(os.path.join(png_path, '_RNAfatsa.jpg'),bbox_inches='tight',dpi=100,pad_inches=0.05)
    plt.close()
    imgpath = os.path.join(png_path, '_RNAfatsa.jpg')
    return imgpath
########画聚类热图---------------------------------------------
def makeHeatmap(npdata):
    tsne_png_path = os.path.join(os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."),
                                 'statics/images')
    pngfile = os.path.join(tsne_png_path, '_Resheatmap.jpg')
    if os.path.exists(pngfile):  # 如果文件存在
        # 删除文件，可使用以下两种方法。
        os.remove(pngfile)
    #将三维矩阵转为二维然后像二维矩阵一样处理
    if len(npdata.shape) == 3:
        s1 = npdata.shape[0]
        s2 = npdata.shape[1]
        s3 = npdata.shape[2]

        data1 = npdata.reshape([s1,s2*s3])
        npdata = data1
    # 画heatmap聚类图##################
    m_scaler = preprocessing.MinMaxScaler()
    heatdata = m_scaler.fit_transform(npdata)

    df = pd.DataFrame(heatdata,
                      index=[('Pair' + str(i)) for i in range(1, npdata.shape[0]+1)],#DataFrame的行标签设置为大写字母
                      columns=[('Feature' + str(i)) for i in range(1, npdata.shape[1]+1)])#设置DataFrame的列标签
    plt.clf()
    plt.figure(dpi=300)
    sns.clustermap(df,method ='ward',metric='euclidean',
                  cmap=plt.get_cmap('Blues'),#matplotlib中的颜色盘'Greens'
                  )
    plt.title('Heatmap')
    png_path = os.path.join(os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."),
                            'statics/images')

    plt.savefig(os.path.join(png_path, '_Resheatmap.jpg'))
    plt.close()


#画tsne无监督聚类图##################
def makeTsne(npdata):
    tsne_png_path = os.path.join(os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."),
                                 'statics/images')
    pngfile = os.path.join(tsne_png_path, '_ResTsnemap.jpg')
    if os.path.exists(pngfile):  # 如果文件存在
        # 删除文件，可使用以下两种方法。
        os.remove(pngfile)
    #将三维矩阵转为二维然后像二维矩阵一样处理
    if len(npdata.shape) == 3:
        s1 = npdata.shape[0]
        s2 = npdata.shape[1]
        s3 = npdata.shape[2]

        data1 = npdata.reshape([s1,s2*s3])
        npdata = data1
        # 画tsne无监督聚类图##################
        # 数据归一化
    sc = preprocessing.StandardScaler()
    data_std = sc.fit_transform(npdata)

    tsne = TSNE(n_components=2, learning_rate=200)
    data_tsne = tsne.fit_transform(data_std)
    plt.clf()
    plt.figure(dpi=300)
    plt.scatter(data_tsne[:, 0], data_tsne[:, 1])
    plt.title('t-SNE Map')
    tsne_png_path = os.path.join(os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."),
                                 'statics/images')
    plt.savefig(os.path.join(tsne_png_path, '_ResTsnemap.jpg'))
    plt.close()
###将两个上传的fasta文件根据第三个文件排序------------------------------
class Process_File(object):
    def __init__(self, filepath1, filepath2, filepath3):
        self.filepath1 = filepath1
        self.filepath2 = filepath2
        self.filepath3 = filepath3
        # self.savepath_a = savepath_a
        # self.savepath_b = savepath_b
    def gen_dict(self, fiel):
        name = []
        seq_all = []
        # print('fiel')
        # print(fiel)
        for seq in Seq.parse(fiel, 'fasta'):
            seq_name = seq.id
            seq_name = '>' + str(seq_name)
            sequence = str(seq.seq)
            name.append(seq_name)
            seq_all.append(sequence)
        # print(name)
        # print(seq_all)
        name_seq = dict(zip(name, seq_all))
        return name_seq

    def get_interaction_pair(self):
        A_dict = self.gen_dict(self.filepath1)
        B_dict = self.gen_dict(self.filepath2)
        # A_save = open(self.savepath_a + '/' + 'new_a.txt', 'w')
        # B_save = open(self.savepath_b + '/' + 'new_b.txt', 'w')
        A_save = open(self.filepath1, 'w')
        B_save = open(self.filepath2, 'w')
        # print(A_dict)
        with open(self.filepath3, 'r') as f:
            for line in f:
                if '>' in line:
                    line = line.strip('\n').split(',')
                    # print('line')
                    # print(line)
                    A = line[0]
                    B = line[1]
                    a_seq = A_dict[A]
                    b_seq = B_dict[B]
                    A_save.write(A + '\n' + a_seq + '\n')
                    B_save.write(B + '\n' + b_seq + '\n')

        A_save.close()
        B_save.close()

###将单个上传的fasta文件根据第二个cav文件排序------------------------------
class order_File(object):
    def __init__(self, filepath1, filepath3):
        self.filepath1 = filepath1
        # self.filepath2 = filepath2
        self.filepath3 = filepath3
        # self.savepath_a = savepath_a
        # self.savepath_b = savepath_b
    def gen_dict(self, fiel):
        name = []
        seq_all = []
        for seq in SeqIO.parse(fiel, 'fasta'):
            seq_name = seq.id
            seq_name = '>' + seq_name
            sequence = str(seq.seq)
            name.append(seq_name)
            seq_all.append(sequence)
        name_seq = dict(zip(name, seq_all))
        return name_seq

    def get_interaction_pair(self):
        A_dict = self.gen_dict(self.filepath1)
        # B_dict = self.gen_dict(self.filepath2)
        # A_save = open(self.savepath_a + '/' + 'new_a.txt', 'w')
        # B_save = open(self.savepath_b + '/' + 'new_b.txt', 'w')
        A_save = open(self.filepath1, 'w')
        # B_save = open(self.filepath2, 'w')
        # print(A_dict)
        with open(self.filepath3, 'r') as f:
            for line in f:
                if '>' in line:
                    line = line.strip('\n').split(',')
                    A = line[0]
                    # B = line[1]
                    a_seq = A_dict[A]
                    # b_seq = B_dict[B]
                    A_save.write(A + '\n' + a_seq + '\n')
                    # B_save.write(B + '\n' + b_seq + '\n')

        A_save.close()
        # B_save.close()
###将单个上传的df文件根据第二个cav文件排序------------------------------
class order_DfFile(object):
    def __init__(self, filepath1, filepath3):
        self.filepath1 = filepath1
        # self.filepath2 = filepath2
        self.filepath3 = filepath3
        # self.savepath_a = savepath_a
        # self.savepath_b = savepath_b
    def gen_dict(self, file):
        dffile = pd.read_csv(file)
        indexcol = list(dffile.iloc[:,1])
        dict = {}
        for i, element in enumerate(indexcol):
            dict[element] = i
        return dict

    def gen_index(self, file):
        dffile = pd.read_csv(file)
        indexrow = list(dffile['Name'])

        dict = {}
        for i, element in enumerate(indexrow):
            dict[element] = i
        return dict

    def get_interaction_pair(self):
        rowindex_dict = self.gen_index(self.filepath1)
        # print('rowindex_dict')
        # print(rowindex_dict)

        labelfile = pd.read_csv(self.filepath3)
        samllfile = pd.read_csv(self.filepath1)
        # print('samllfile')
        # print(samllfile)
        internames = list(labelfile.iloc[:,1])
        # print('internames')
        # print(internames)
        arrange_index = []
        for intername in internames:

            arrange_index.append(rowindex_dict[intername])


        select_row = samllfile.loc[arrange_index, :]
        select_row.index = range(len(arrange_index))


        select_row.to_csv(os.path.join(self.filepath1), index=False, header=True)
########根据用户上传的label文件对画聚类热图---------------------------------------------

def rnaheatmap(rnadfdata,tsne_png_path,Labelfilepath = None,sename = None,method ='ward',metric='euclidean'):
    tsne_png_path = tsne_png_path
    pngfile = os.path.join(tsne_png_path, '_rnaheatmap.jpg')
    if os.path.exists(pngfile):  # 如果文件存在
        # 删除文件，可使用以下两种方法。
        os.remove(pngfile)
    # 将输入数据整理为画图所用数据##################
    m_scaler = preprocessing.MinMaxScaler()
    heatdata = m_scaler.fit_transform(rnadfdata)
    if sename:
        df = pd.DataFrame(heatdata,
                          index=sename,  # DataFrame的行标签设置为大写字母
                          columns=[('Feature' + str(i)) for i in range(1, heatdata.shape[1]+1)])  # 设置DataFrame的列标签
    else:
        df = pd.DataFrame(heatdata,
                          index=rnadfdata.index.tolist(),  # DataFrame的行标签设置为大写字母
                          columns=rnadfdata.columns.tolist())  # 设置DataFrame的列标签


    # 画heatmap聚类图##################
    if Labelfilepath:

        labelcol = []
        labelc = []
        Labelfile = open(Labelfilepath, 'r')
        for line in Labelfile:
            if '>' in line:
                # line = line.strip('\n').split(',')
                line = line.strip('\n').split('\t')
                lab = line[2]
                lab01 = int(lab)
                labelc.append(lab01)

        df['Label'] = labelc
        plt.clf()
        plt.figure(dpi=300)
        Colors=['#49759c','#a2cffe','#448ee4','#8ab8fe','#CEFFCE','#28FF28','#007500','#FFFF93','#8C8C00','#FFB5B5','#FF0000','#CE0000','#750000']
        len_labelj = len(df['Label'].unique())
        row_c = dict(zip(df['Label'].unique(), list(Colors[0:len_labelj])))
        sns.clustermap(df.drop(columns=['Label']),method ='ward',metric='euclidean',
                      cmap=plt.get_cmap('Blues'),#matplotlib中的颜色盘'Greens'
                      row_colors=df['Label'].map(row_c), #行方向聚类用颜色区分不同类
                      )
    else:
        plt.clf()
        plt.figure(dpi=300)
        sns.clustermap(df, method='ward', metric='euclidean',
                       cmap=plt.get_cmap('Blues'),  # matplotlib中的颜色盘'Greens'
                       )
    plt.title('Heatmap')
    png_path = tsne_png_path
    # shutil.rmtree(png_path)
    # os.mkdir(png_path)
    plt.savefig(os.path.join(png_path, '_rnaheatmap.jpg'))
    plt.close()

########根据用户上传的RNA-RNA 相互作用 label文件对画聚类热图---------------------------------------------

def rnabefheatmap(rnadfdata,tsne_png_path,Labelfilepath = None,sename = None,method ='ward',metric='euclidean'):
    tsne_png_path = tsne_png_path
    pngfile = os.path.join(tsne_png_path, '_rnaheatmap.jpg')
    if os.path.exists(pngfile):  # 如果文件存在
        # 删除文件，可使用以下两种方法。
        os.remove(pngfile)
    # 将输入数据整理为画图所用数据##################
    m_scaler = preprocessing.MinMaxScaler()
    heatdata = m_scaler.fit_transform(rnadfdata)
    if sename:
        df = pd.DataFrame(heatdata,
                          index=sename,  # DataFrame的行标签设置为大写字母
                          columns=[('Feature' + str(i)) for i in range(1, heatdata.shape[1]+1)])  # 设置DataFrame的列标签
    else:
        df = pd.DataFrame(heatdata,
                          index=rnadfdata.index.tolist(),  # DataFrame的行标签设置为大写字母
                          columns=rnadfdata.columns.tolist())  # 设置DataFrame的列标签


    # 画heatmap聚类图##################
    if Labelfilepath:

        labelcol = []
        labelc = []
        Labelfile = open(Labelfilepath, 'r')
        for line in Labelfile:
            if '>' in line:
                # line = line.strip('\n').split(',')
                line = line.strip('\n').split('\t')
                lab = line[2]
                lab01 = int(lab)
                labelc.append(lab01)

        df['Label'] = labelc
        plt.clf()
        plt.figure(dpi=300)
        Colors=['#49759c','#a2cffe','#448ee4','#8ab8fe','#CEFFCE','#28FF28','#007500','#FFFF93','#8C8C00','#FFB5B5','#FF0000','#CE0000','#750000']
        len_labelj = len(df['Label'].unique())
        row_c = dict(zip(df['Label'].unique(), list(Colors[0:len_labelj])))
        sns.clustermap(df.drop(columns=['Label']),method ='ward',metric='euclidean',
                      cmap=plt.get_cmap('Blues'),#matplotlib中的颜色盘'Greens'
                      row_colors=df['Label'].map(row_c), #行方向聚类用颜色区分不同类
                      )
        plt.title('Sequence Similarity ClusterMap before Model Training')
    else:
        plt.clf()
        plt.figure(dpi=300)
        sns.clustermap(df, method='ward', metric='euclidean',
                       cmap=plt.get_cmap('Blues'),  # matplotlib中的颜色盘'Greens'
                       )
        plt.title('Sequence Similarity ClusterMap before Model Training')

    png_path = tsne_png_path
    # shutil.rmtree(png_path)
    # os.mkdir(png_path)
    plt.savefig(os.path.join(png_path, '_rnaheatmap.jpg'))
    plt.close()
#画RNA PCA聚类图##################
def makePCA(rnadfdata,Labelfilepath = None,n_components=2):
    tsne_png_path = os.path.join(os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."),
                                 'statics/images')
    pngfile = os.path.join(tsne_png_path, '_PCAmap.jpg')
    if os.path.exists(pngfile):  # 如果文件存在
        # 删除文件，可使用以下两种方法。
        os.remove(pngfile)

    #将三维矩阵转为二维然后像二维矩阵一样处理
    if len(rnadfdata.shape) == 3:
        s1 = rnadfdata.shape[0]
        s2 = rnadfdata.shape[1]
        s3 = rnadfdata.shape[2]

        data1 = rnadfdata.reshape([s1,s2*s3])
        heatdata = data1
    else:
        heatdata = np.array(rnadfdata)
    pca = PCA(n_components=2)
    data_pca = pca.fit_transform(heatdata)
    if Labelfilepath:
        labelcol = []
        labelc = []
        Labelfile = open(Labelfilepath, 'r')
        for line in Labelfile:
            if '>' in line:
                line = line.strip('\n').split(',')
                lab = line[1]
                lab01 = int(lab)
                labelc.append(lab01)
        labelc01 = np.array(labelc)
        x_all = {}
        y_all = {}
        for x in range(len(set(labelc))):
            xlistname = 'listx_' + str(x)
            x_all.update({xlistname: []})
            ylistname = 'listy_' + str(x)
            y_all.update({ylistname: []})
        for i in range(len(data_pca)):
            for n in range(len(list(set(labelc)))):
                if labelc01[i] == list(set(labelc))[n]:
                    xlistname = 'listx_' + str(n)
                    ylistname = 'listy_' + str(n)
                    x_all[xlistname].append(data_pca[i][0])
                    y_all[ylistname].append(data_pca[i][1])
        Colors = ['#49759c', '#a2cffe', '#448ee4', '#8ab8fe', '#CEFFCE', '#28FF28', '#007500', '#FFFF93', '#8C8C00',
                  '#FFB5B5', '#FF0000', '#CE0000', '#750000']
        markers = ['x', '.', 'D']
        plt.clf()
        plt.figure(dpi=300)
        for n in range(len(list(set(labelc)))):
            xlistname = 'listx_' + str(n)
            ylistname = 'listy_' + str(n)
            plt.scatter(x_all[xlistname], y_all[ylistname], c=Colors[n], marker=markers[n])
    else:
        plt.clf()
        plt.figure(dpi=300)
        plt.scatter(data_pca[:, 0], data_pca[:, 1], c='#49759c')

    plt.title('PCA Map')
    tsne_png_path = os.path.join(os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + os.path.sep + ".."),
                                 'statics/images')
    plt.savefig(os.path.join(tsne_png_path, '_PCAmap.jpg'))
    plt.close()
##########根据特征值对每一个特征画图————————————————————————————————————————————

# @tornado.concurrent.run_on_executor
def displaypng(dfdata,filepath,labellist = None):
    colname = dfdata.columns.values.tolist()
    col_len = dfdata.shape[1]

    # if labellist:
    #     dfdata['label'] = labellist
    for num in list(range(0,col_len,1)):
        pngname = re.sub(r'\s+', '_', colname[num])
        #Bar 图
        plt.clf()
        plt.figure(dpi=300)
        plt.hist(dfdata.iloc[:,num], 50, rwidth=0.8,facecolor='#427A8B', align='left')
        pngname0 =  pngname+'0.jpg'
        plt.savefig(os.path.join(filepath,pngname0))
        plt.close()
        #箱线图
        if labellist:
            plt.clf()
            plt.figure(dpi=300)
            sns.boxplot(x=labellist, y = dfdata.iloc[:,num])
            pngname1 =  pngname+'1.jpg'
            plt.savefig(os.path.join(filepath,pngname1))
            plt.close()
        else:
            plt.clf()
            plt.figure(dpi=300)
            sns.boxplot(y = dfdata.iloc[:,num])
            pngname1 =  pngname+'1.jpg'
            plt.savefig(os.path.join(filepath,pngname1))
            plt.close()
         #密度分布图
        plt.clf()
        plt.figure(dpi=300)
        # sns.kdeplot(x = dfdata.iloc[:,num],shade=True)
        sns.kdeplot(dfdata.iloc[:, num], shade=True)
        pngname2 =  pngname+'2.jpg'
        plt.savefig(os.path.join(filepath,pngname2))
        plt.close()
        #Warm图
        if labellist:
            plt.clf()
            plt.figure(dpi=300)
            sns.swarmplot(x=labellist,y = dfdata.iloc[:,num])
            pngname3 =  pngname+'3.jpg'
            plt.savefig(os.path.join(filepath,pngname3))
            plt.close()
        else:
            plt.clf()
            plt.figure(dpi=300)
            sns.swarmplot(y = dfdata.iloc[:,num])
            pngname3 =  pngname+'3.jpg'
            plt.savefig(os.path.join(filepath,pngname3))
            plt.close()

        # 小提琴图
        plt.clf()
        plt.figure(dpi=300)
        sns.violinplot(y = dfdata.iloc[:,num])
        pngname4 =  pngname+'4.jpg'
        plt.savefig(os.path.join(filepath,pngname4))
        plt.close()

##二维的数据结果展示图
def displayplot2d(npdata,B_2Dmeths,savepath,labellist =None):
    leng = npdata.shape[2] + 1
    colnames = []
    for i in range(1, leng):
        strs = ('Feature%i' % i)
        colnames.append(strs)
    #确定一级特征名字
    dict_2Dfeature = {
        3:['Physicochemical features'],
        4:['Sequence intrinsic features'],
        7:['Sequence intrinsic features','Physicochemical features'],
        10:['Sequence intrinsic features','Structural features'],
        11:['Sequence intrinsic features','Structural features'],
        14:['Sequence intrinsic features','Physicochemical features','Structural features']
    }
    if B_2Dmeths == 'one1':
        Feat01s = ['Structural features']
    else:
        Feat01s = dict_2Dfeature[npdata.shape[2]]
    Feat01s = sorted(set(Feat01s), key=Feat01s.index)
    for num_Feat01 in list(range(0, len(Feat01s), 1)):
        for num in list(range(0, len(colnames), 1)):
            # pngname = re.sub(r'\s+', '_', Feat01s[num_Feat01])+colnames[num])
            pngname = re.sub(r'\s+', '_', colnames[num])
            # print(pngname)
            #Bar 图
            plt.clf()
            plt.figure(dpi=300)
            plt.hist(np.mean(npdata[:,:,num], axis=1), 50, rwidth=0.8,facecolor='#427A8B', align='left')
            pngname0 =  pngname+'0.jpg'
            plt.savefig(os.path.join(savepath,pngname0))
            plt.close()
            #箱线图
            if labellist:
                plt.clf()
                plt.figure(dpi=300)
                sns.boxplot(x=labellist, y = np.mean(npdata[:,:,num], axis=1))
                pngname1 =  pngname+'1.jpg'
                plt.savefig(os.path.join(savepath,pngname1))
                plt.close()
            else:
                plt.clf()
                plt.figure(dpi=300)
                sns.boxplot(y = np.mean(npdata[:,:,num], axis=1))
                pngname1 =  pngname+'1.jpg'
                plt.savefig(os.path.join(savepath,pngname1))
                plt.close()
            #密度分布图
            plt.clf()
            plt.figure(dpi=300)
            sns.kdeplot(data = np.mean(npdata[:,:,num], axis=1),shade=True)
            pngname2 =  pngname+'2.jpg'
            plt.savefig(os.path.join(savepath,pngname2))
            plt.close()

#只对一列数字画图
def displayOnepng(dfdata,num,filepath,labellist = None):
    colname = dfdata.columns.values.tolist()
    pngname = re.sub(r'\s+', '_', colname[num])
    print('colname length')
    print(dfdata.iloc[:,num])
    print('label')
    print(labellist)
    #Bar 图
    plt.clf()
    plt.figure(dpi=300)
    plt.hist(dfdata.iloc[:,num], 50, rwidth=0.8,facecolor='#427A8B', align='left')
    pngname0 =  pngname+'0.jpg'
    plt.savefig(os.path.join(filepath,pngname0))
    plt.close()
    #箱线图
    if labellist:
        plt.clf()
        plt.figure(dpi=300)
        sns.boxplot(x=labellist, y = dfdata.iloc[:,num])
        pngname1 =  pngname+'1.jpg'
        plt.savefig(os.path.join(filepath,pngname1))
        plt.close()
    else:
        plt.clf()
        plt.figure(dpi=300)
        sns.boxplot(y = dfdata.iloc[:,num])
        pngname1 =  pngname+'1.jpg'
        plt.savefig(os.path.join(filepath,pngname1))
        plt.close()
    #密度分布图
    if len(set(dfdata.iloc[:,num])) == 1:
        pass
    else:
        plt.clf()
        plt.figure(dpi=300)
        sns.kdeplot(x = dfdata.iloc[:,num],shade=True)
        pngname2 =  pngname+'2.jpg'
        plt.savefig(os.path.join(filepath,pngname2))
        plt.close()



#二维单列数据展示---------------------------------------------------------

def displayOneplot2d(npdata,savepath,labellist,colnames,num):
    # pngname = re.sub(r'\s+', '_', Feat01s[num_Feat01])+colnames[num])
    pngname = re.sub(r'\s+', '_', colnames[num])
    # print(pngname)
    #Bar 图
    plt.clf()
    plt.figure(dpi=300)
    plt.hist(np.mean(npdata[:,:,num], axis=1), 50, rwidth=0.8,facecolor='#427A8B', align='left')
    pngname0 =  pngname+'0.jpg'
    plt.savefig(os.path.join(savepath,pngname0))
    plt.close()
    #箱线图
    if labellist:
        plt.clf()
        plt.figure(dpi=300)
        sns.boxplot(x=labellist, y = np.mean(npdata[:,:,num], axis=1))
        pngname1 =  pngname+'1.jpg'
        plt.savefig(os.path.join(savepath,pngname1))
        plt.close()
    else:
        plt.clf()
        plt.figure(dpi=300)
        sns.boxplot(y = np.mean(npdata[:,:,num], axis=1))
        pngname1 =  pngname+'1.jpg'
        plt.savefig(os.path.join(savepath,pngname1))
        plt.close()
    #密度分布图
    plt.clf()
    plt.figure(dpi=300)
    sns.kdeplot(data = np.mean(npdata[:,:,num], axis=1),shade=True)
    pngname2 =  pngname+'2.jpg'
    plt.savefig(os.path.join(savepath,pngname2))
    plt.close()

#### 给所有的名字中间空格加上下划线-------------------------
def exchange(x):
    featurenames = str(x).split(": ")
    newname = str(featurenames[1]) + " (" + str(featurenames[0]) + ")"
    newname = re.sub(r'\s+', '_', newname)
    print(newname)
    return newname
##一维的数据结果按照特征数量展示图
def displaypng01(dfdata, filepath, requireFeature, figuretype,labellist=None):
    colname = dfdata.columns.values.tolist()
    col_len = dfdata.shape[1]

    for i, j in enumerate(colname):
        featurenames = j.split(": ")
        newname = str(featurenames[1]) + " (" + str(featurenames[0]) + ")"
        newname = re.sub(r'\s+', '_', newname)
        # j = re.sub(r'\s+', '_', j)
        # print(newname)
        # print(type(newname))


        if newname == str(requireFeature):
            IDindex = i
            print('IDindex')
            print(IDindex)

    pngname = colname[IDindex]

    if figuretype == '0':
        # 小提琴图
        if labellist:
            plt.clf()
            plt.figure(dpi=300)
            #给画的图按照label分配不同的颜色
            labellist01s = sorted(set(labellist), key=labellist.index)
            my_pal = {}
            colors = ['#437A8B','#C23147','#5F86CC','#F09150','#AA65C7']
            for order,labellist01 in enumerate(labellist01s):
                my_pal[labellist01] = colors[order]
            # print(my_pal)
            sns.violinplot(x=labellist, y=dfdata.iloc[:, IDindex],palette=my_pal)
            pngname4 = pngname + '0.jpg'
            pngpath = os.path.join(filepath, pngname4)
            plt.savefig(os.path.join(filepath, pngname4))
            plt.close()
        else:
            plt.clf()
            plt.figure(dpi=300,figsize=(8, 5))
            # 给画的图按照label分配不同的颜色
            labellist01s = sorted(set(labellist), key=labellist.index)
            my_pal = {}
            colors = ['#437A8B', '#C23147', '#5F86CC', '#F09150', '#AA65C7']
            for order, labellist01 in enumerate(labellist01s):
                my_pal[labellist01] = colors[order]

            sns.violinplot(y=dfdata.iloc[:, IDindex],palette=my_pal)
            pngname4 = pngname + '0.jpg'
            pngpath = os.path.join(filepath, pngname4)
            plt.savefig(os.path.join(filepath, pngname4))
            plt.close()


    elif figuretype == '1':
        # 密度分布图
        plt.clf()
        plt.figure(dpi=300,figsize=(8, 5))
        # sns.kdeplot(x = dfdata.iloc[:,num],shade=True)
        sns.kdeplot(dfdata.iloc[:, IDindex], shade=True,color= "#5F86CC")
        pngname2 = pngname + '1.jpg'
        pngpath = os.path.join(filepath, pngname2)
        plt.savefig(pngpath)
        plt.close()
    elif figuretype == '2':

        # Bar 图
        plt.clf()
        plt.figure(dpi=300,figsize=(8, 5))
        plt.hist(dfdata.iloc[:, IDindex], 50, rwidth=0.8, facecolor='#F09150', align='left')
        pngname0 = pngname + '2.jpg'
        pngpath = os.path.join(filepath, pngname0)
        plt.savefig(pngpath)
        plt.close()

    elif figuretype == '3':
        plt.clf()
        plt.figure(dpi=300,figsize=(8, 5))
        sns.violinplot(y=dfdata.iloc[:, IDindex],color = '#AA65C7')
        pngname4 = pngname + '3.jpg'
        pngpath = os.path.join(filepath, pngname4)
        plt.savefig(os.path.join(filepath, pngname4))
        plt.close()



    return pngpath
##二维的数据结果按照特征数量展示图
def displayplot2d01(npdata,feature,figuretype,savepath,labellist =None):
    leng = npdata.shape[2] + 1
    # print('npdata.shape')
    # print(npdata.shape)
    if feature == 'Sequence-based_one-hot_encoding_(SeqOH)':
        npdata = npdata[:,:,0:3]
    elif feature == 'Physico-chemical_sparse_encoding_(SECod)':
        if leng == 4:
            npdata = npdata
            # print('npdata.shape')
            # print(npdata.shape)
            # print(npdata.head)
        else:
            npdata = npdata[:,:,3:7]
    elif feature == 'Structural_features':
        npdata = npdata[:, :, -7:]
    npdata_me = np.mean(npdata, axis=2)
    # pngname = re.sub(r'\s+', '_', Feat01s[num_Feat01])+colnames[num])
    pngname = re.sub(r'_', ' ', feature)
    # print(pngname)
    if figuretype == '0':
        #小提琴图
        if labellist:
            plt.clf()
            plt.figure(dpi=300,figsize=(8, 5))

            #给画的图按照label分配不同的颜色
            labellist01s = sorted(set(labellist), key=labellist.index)
            my_pal = {}
            colors = ['#437A8B','#C23147','#5F86CC','#F09150','#AA65C7']
            for order,labellist01 in enumerate(labellist01s):
                my_pal[labellist01] = colors[order]
            # print(my_pal)
            # my_pal = {"versicolor": "g", "setosa": "b", "virginica": "m"}

            sns.violinplot(x=labellist, y=np.mean(npdata_me, axis=1),palette=my_pal)
            pngname4 =  pngname+'0.jpg'
            pngpath = os.path.join(savepath, pngname4)
            plt.savefig(os.path.join(savepath,pngname4))
            plt.close()
        else:
            plt.clf()
            plt.figure(dpi=300,figsize=(8, 5))
            #给画的图按照label分配不同的颜色
            labellist01s = sorted(set(labellist), key=labellist.index)
            my_pal = {}
            colors = ['#437A8B','#C23147','#5F86CC','#F09150','#AA65C7']
            for order,labellist01 in enumerate(labellist01s):
                my_pal[labellist01] = colors[order]
            # print(my_pal)

            sns.violinplot(y=np.mean(npdata_me, axis=1),palette=my_pal)
            pngname4 =  pngname+'0.jpg'
            pngpath = os.path.join(savepath, pngname4)
            plt.savefig(os.path.join(savepath,pngname4))
            plt.close()


    elif figuretype == '1':
        # 密度分布图
        plt.clf()
        plt.figure(dpi=300, figsize=(8, 5))
        sns.kdeplot(data=np.mean(npdata_me, axis=1), shade=True,color= "#5F86CC")
        # 设置刻度字体大小
        # plt.xticks(fontsize=20)
        # plt.yticks(fontsize=20)
        pngname2 = pngname + '1.jpg'
        pngpath = os.path.join(savepath, pngname2)
        plt.savefig(os.path.join(savepath, pngname2))
        plt.close()

    elif figuretype == '2':
        #Bar 图
        plt.clf()
        plt.figure(dpi=300,figsize=(8, 5))
        plt.hist(np.mean(npdata_me, axis=1), 50, rwidth=0.8,facecolor='#F09150', align='left')
        # 设置刻度字体大小
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        pngname0 =  pngname+'2.jpg'
        pngpath = os.path.join(savepath, pngname0)
        plt.savefig(os.path.join(savepath,pngname0))
        plt.close()
    elif figuretype == '3':
        #小提琴图
        plt.clf()
        plt.figure(dpi=300,figsize=(8, 5))
        sns.violinplot(y = np.mean(npdata_me, axis=1),color = '#AA65C7')
        pngname4 =  pngname+'3.jpg'
        pngpath = os.path.join(savepath, pngname4)
        plt.savefig(os.path.join(savepath,pngname4))
        plt.close()


    return pngpath
#############生成合并后的一维下载数据————————————————————————————————————————————————————————
def makeDownload_1D(Aresult_F,Bresult_F,Result_path):
    A_res = np.array(Aresult_F)
    B_res = np.array(Bresult_F)
    # 生成序列名字
    Afileseqname = Aresult_F._stat_axis.values.tolist()
    Bfileseqname = Bresult_F._stat_axis.values.tolist()
    combineNameAB = []
    for i in range(0, len(Afileseqname)):
        combineNameAB.append('%s_%s' % (Afileseqname[i], Bfileseqname[i]))
    combineNameBA = []
    for i in range(0, len(Afileseqname)):
        combineNameBA.append('%s_%s' % (Bfileseqname[i], Afileseqname[i]))
    # 生成特征名字
    Afilefeaname = Aresult_F.columns.values.tolist()
    Bfilefeaname = Bresult_F.columns.values.tolist()
    ABfeaname = []
    ABfeaname.extend(Afilefeaname)
    ABfeaname.extend(Bfilefeaname)
    BAfeaname = []
    BAfeaname.extend(Bfilefeaname)
    BAfeaname.extend(Afilefeaname)
    # 按照AB方式合并
    ComResultAB_np = np.concatenate((A_res, B_res), axis=1)
    ComResultAB = pd.DataFrame(data=ComResultAB_np, index=combineNameAB, columns=ABfeaname)
    print('shape(ComResultAB)')
    print(ComResultAB.shape)

    makeLonelyDownload_1D(ComResultAB, Result_path, 'Result—AB', npfile=ComResultAB_np, demension=2)
    # ComResultAB.to_json(os.path.join(Result_path, 'Result—AB.json'), orient='split')
    # 按照BA方式合并
    ComResultBA_np = np.concatenate((B_res, A_res), axis=1)
    ComResultBA = pd.DataFrame(data=ComResultBA_np, index=combineNameBA, columns=BAfeaname)
    print('shape(ComResultBA)')
    print(ComResultBA.shape)

    makeLonelyDownload_1D(ComResultBA, Result_path, 'Result—BA', npfile=ComResultBA_np, demension=2)
    return ComResultAB_np,ComResultAB

#############生成合并后的二维下载数据————————————————————————————————————————————————————————
def makeDownload_2D(Aresult_F,Bresult_F,Result_path):
    print('Aresult_F.shape')
    print(Aresult_F.shape)
    print('Bresult_F.shape')
    print(Bresult_F.shape)
    # 按照AB方式合并
    ComResultAB_np = np.concatenate((Aresult_F, Bresult_F), axis=2)
    np.save(os.path.join(Result_path, 'Result—AB.npy'), ComResultAB_np)
    # 按照BA方式合并
    ComResultBA_np = np.concatenate((Bresult_F, Aresult_F), axis=2)
    np.save(os.path.join(Result_path, 'Result—BA.npy'), ComResultBA_np)
    return ComResultAB_np
#############单个一维数据生成下载数据————————————————————————————————————————————————————————
def makeLonelyDownload_1D(csvfiledf, Result_path,title,npfile = None,demension = 1):
    if demension == 1:
        csvfiledf.to_csv(os.path.join(Result_path, title + '.csv'), index=True, header=True)
        csvfiledf.to_csv(os.path.join(Result_path, title + '.tsv'), sep="\t", index=True, header=True)
        np.save(os.path.join(Result_path, title + '.npy'), np.array(csvfiledf))
    else:
        csvfiledf.to_csv(os.path.join(Result_path, title + '.csv'), index=True, header=True)
        csvfiledf.to_csv(os.path.join(Result_path, title + '.tsv'), sep="\t", index=True, header=True)
        np.save(os.path.join(Result_path, title + '.npy'), npfile)
    # csvfiledf.to_json(os.path.join(Result_path, title + '.npy'), orient="records", force_ascii=False)

def namemoded(origanlname):
    featurenames = str(origanlname).split(": ")
    newname = str(featurenames[1]) + " (" + str(featurenames[0]) + ")"
    return newname



#根据样本数确定batchsize
def make_batsz(samplenum):
    batchsize = samplenum / 15

    for index in range(4, 20, 1):
        if 2 ** index > batchsize:
            break
    bias0 = batchsize - 2 ** (index - 1)
    bias1 = 2 ** index - batchsize
    if bias0 < bias1:
        batchsize01 = 2 ** (index - 1)
    else:
        batchsize01 = 2 ** index
    return batchsize01

# make a zip file from a dir

#打包目录为zip文件（未压缩）
def make_zip(source_dir, output_filename):
    zipf = zipfile.ZipFile(output_filename, 'w')
    pre_len = len(os.path.dirname(source_dir))
    for parent, _, filenames in os.walk(source_dir):
        for filename in filenames:
            pathfile = os.path.join(parent, filename)
            arcname = pathfile[pre_len:].strip(os.path.sep)     #相对路径
            zipf.write(pathfile, arcname)
    zipf.close()
#make independent data
def makeindedata(ComResultAB_np,Interlabel):
    label = np.array(Interlabel['Label'].tolist())
    print('ComResultAB_np.shape')
    print(ComResultAB_np.shape)
    print('len(label)')
    print(len(label))
    # judge whether The number of sequence in coding reault is match the label files length or not
    try:
        if ComResultAB_np.shape[0] != len(label):
            raise ValueError(
                "The number of sequence in coding reault is not match the label files!")
    except ValueError as e:
        print(e)

    noindeindex = [i for i, x in enumerate(label) if x != 'I']
    indeindex = [i for i, x in enumerate(label) if x == 'I']

    indedata = ComResultAB_np[indeindex, :]
    ComResultAB_np = ComResultAB_np[noindeindex, :]

    label = Interlabel.iloc[noindeindex]['Label'].tolist()
    print('Make independent data successfully!')
    label = [int(i) for i in label]
    x = np.array(label[1:])
    x0 = np.array(label[0:1])
    label = np.concatenate([x0, x])

    return ComResultAB_np,label


############### end 方法函数 end—————————————————————————————————————————————————————————