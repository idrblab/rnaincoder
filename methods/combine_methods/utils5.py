

from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from random import randint
import numpy as np
plt.switch_backend('agg')
def plot_clustering1(z_run_tsne, labels, folder_name, title):



    hex_colors = []
    for _ in np.unique(labels):
        hex_colors.append('#%06X' % randint(0, 0xFFFFFF))

    colors = [hex_colors[int(i)] for i in labels]
    # 更新属性设置
    params = {
        'font.family': 'Calibri',
        'font.size': 12
    }
    plt.rcParams.update(params)
    bwith = 1.0
    TK = plt.gca()
    TK.spines['bottom'].set_linewidth(bwith)
    TK.spines['left'].set_linewidth(bwith)
    TK.spines['top'].set_linewidth(bwith)
    TK.spines['right'].set_linewidth(bwith)

    plt.scatter(z_run_tsne[:, 0], z_run_tsne[:, 1], c=colors, s=100, marker='*')
    print("tsne done")

    plt.rcParams['figure.figsize'] = (10, 10)  
    plt.rcParams['savefig.dpi'] = 300  
    plt.rcParams['figure.dpi'] = 300  
    plt.title('TSNE after KMeans')

    plt.savefig(folder_name + title + "tsne.png")
    plt.close()

from sklearn.manifold import TSNE

def kmeans_visual(data, save_path, title,k=5):
    ######
    # data = np.load(data)
    
    if len(data.shape) == 2:
        model = KMeans(n_clusters=k)
        print('np.any(np.isnan(data))')
        np.any(np.isnan(data))
        print('np.all(np.isfinite(data))')
        np.all(np.isfinite(data))
        model.fit(data)
        label = model.labels_
        plot_data = model.cluster_centers_
        z_tsne = TSNE(perplexity=80, min_grad_norm=1E-12, n_iter=1000).fit_transform(data)
        plot_clustering1(z_tsne, label, save_path,
                                   title+ '_kmeans')

    if len(data.shape) == 3:
        data = data.reshape(data.shape[0], -1)
        model = KMeans(n_clusters=k)
        model.fit(data)
        label = model.labels_

        z_tsne = TSNE(perplexity=80, min_grad_norm=1E-12, n_iter=1000).fit_transform(data)
        plot_clustering1(z_tsne, label, save_path,
                                   title + '_kmeans')
        print('2d tsne done')