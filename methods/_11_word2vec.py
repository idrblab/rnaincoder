from gensim.models import Word2Vec
import pandas as pd
import numpy as np
from Bio import SeqIO

def seq_to_kmers(seq, k):
    N = len(seq)
    return [seq[i:i + k] for i in range(N - k + 1)]

def get_protein_embedding(model, protein):
    """get protein embedding,infer a list of 3-mers to (num_word,100) matrix"""
    num_seq = len(protein)

    vec = np.zeros((len(protein), 100))
    i = 0
    for word in protein:
        vec[i,] = model.wv[word]  ########model.wv找到这个词对应的词向量
        i += 1
    return vec

def process_fasta(filepath, modelpath, N):
    model = Word2Vec.load(modelpath)

    all_rna_array = []
    seq_seq = []
    seqname = []
    for seq in SeqIO.parse(filepath,'fasta'):
        seqid = seq.id
        seqname.append(seqid)
        seq_seq.append(seq.seq)
    for lincrna in seq_seq:

        result = get_protein_embedding(model, seq_to_kmers(str(lincrna), 4))
        pad_array = np.zeros((N, 100))
        if len(result) >= N:
            pad_array = result[0:N, :]
        else:
            pad_array[0:len(result), :] = result
        all_rna_array.append(pad_array)
    all_rna_array01 = np.array(all_rna_array)
    return seqname,all_rna_array01
#if __name__ == '__main__':
    #all_rna_array = process_fasta('lincrna_try.txt', 'word2vec_4gram_100_20.model', 2000)