import matplotlib.pyplot as plt
plt.switch_backend('agg')
from sklearn import manifold

import os


def plot_clustering_matplotlib(z_run, labels, folder_name, title, title1):

    pngfile = 'tsne.png.png'
    if os.path.exists(pngfile):

        os.remove(pngfile)

    hex_colors = ['#E68223', '#518C2D']

    z_run_tsne = manifold.TSNE(perplexity=80, min_grad_norm=1E-12, n_iter=1000).fit_transform(z_run)
    plt.figure(figsize=(10, 10))
    plt.rcParams['savefig.dpi'] = 300  # Í¼Æ¬ÏñËØ
    plt.rcParams['figure.dpi'] = 300  # ·Ö±æÂÊ
    plt.rcParams.update({'font.size': 25})
    bwith = 2.0
    TK = plt.gca()
    TK.spines['bottom'].set_linewidth(bwith)
    TK.spines['left'].set_linewidth(bwith)
    TK.spines['top'].set_linewidth(bwith)
    TK.spines['right'].set_linewidth(bwith)
    # ax = plt.subplot(111)
    type1_x = []
    type1_y = []
    type2_x = []
    type2_y = []

    for i in range(z_run_tsne.shape[0]):
        if labels[i] == 0:
            type1_x.append(z_run_tsne[i][0])
            type1_y.append(z_run_tsne[i][1])
        if labels[i] == 1:
            type2_x.append(z_run_tsne[i][0])
            type2_y.append(z_run_tsne[i][1])

    type1 = plt.scatter(type1_x, type1_y, s=150, c=hex_colors[0], marker='o')
    type2 = plt.scatter(type2_x, type2_y, s=150, c=hex_colors[1], marker='o')

    plt.legend((type1, type2), ('0', '1'), loc="upper left")
    # plt.legend(frameon=False,loc="lower right",fontsize='large') #ÉèÖÃÍ¼ÀýÎÞ±ß¿ò£¬½«Í¼Àý·ÅÔÚ×óÉÏ½Ç

    plt.xticks()
    plt.yticks()

    plt.title('TSNE' + title1)
    plt.savefig(folder_name + title + "tsne.png")
    plt.close()