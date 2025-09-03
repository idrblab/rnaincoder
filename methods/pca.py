#coding=gbk
import numpy as np
from sklearn.decomposition import PCA
from Bio import SeqIO

def reduce_pca(file):

    # encoded = np.load(file)
    encoded = file

    n = encoded.shape[0]

    encoded = encoded.reshape(n, -1)
    if n >= 256:
        n_components = 256
    if n < 256:
        n_components = n
    pca = PCA(n_components)
    pca.fit(encoded)
    newX = pca.fit_transform(encoded)
    return newX


