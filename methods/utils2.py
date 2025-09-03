import matplotlib.pyplot as plt
plt.switch_backend('agg')
from sklearn import manifold
import os
def plot_metrics(history, path):
    metrics = ['loss', 'accuracy', 'precision', 'recall']
    bwith = 2
    plt.figure(figsize=(30, 7))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    for n, metric in enumerate(metrics):
        name = metric
        plt.subplot(1, 4, n + 1)
        # plt.suptitle('Various Measurements of Model Performance', fontsize=25, x=0.5, y=1.001)
        TK = plt.gca()
        TK.spines['bottom'].set_linewidth(bwith)
        TK.spines['left'].set_linewidth(bwith)
        TK.spines['top'].set_linewidth(bwith)
        TK.spines['right'].set_linewidth(bwith)


        plt.plot(history.epoch, history.history[metric], color='green', lw=2, label='Train')
        plt.plot(history.epoch, history.history['val_' + metric], lw=2,
                 color='red', linestyle="--", label='Val')
        plt.xlabel('Epochs')
        plt.ylabel(name)
        if metric == 'loss':
            plt.ylim([0, plt.ylim()[1]])
        # plt.legend(fontsize=20, loc='upper right')
        elif metric == 'precision':
            plt.ylim([0.5, 1])

        elif metric == 'recall':
            plt.ylim([0.5, 1])
        else:
            plt.ylim([0, 1])  ###plt.ylimÉèÖÃ×ø±êÖáµÄ·¶Î§

        plt.legend(fontsize=20)
        # plt.tight_layout()
    # plt.figure(figsize=(10, 10))
    # plt.suptitle('Various Measurements of Model', x=0.5, y=1.01, fontsize=25)
    plt.suptitle('Various Measurements of Model Performance', fontsize=25, x=0.5, y=1.001)
    plt.savefig(path + '/rna_rna_metric.png', dpi=500)
    plt.close()