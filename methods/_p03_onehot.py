"""
theme: one_hot
author: xiaweiqi
data: 2020.11.03
"""

from numpy import argmax
from Bio import SeqIO
AA = ["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
AA_to_num = {k:v for v, k in enumerate(AA)}
num_to_AA = {k:v for v, k in AA_to_num.items()}
import numpy as np
def encode(protein_sequence):
    integer_encoded = []
    for i in protein_sequence:
        if i in AA:
            integer_encoded.append(AA_to_num[i])
        else:
            integer_encoded.append(int(90))

    onehot_encoded = list()
    for value in integer_encoded:
        if value != 90:

            letter = [0 for _ in range(len(num_to_AA))]
            letter[value] = 1
            onehot_encoded.append(letter)
        else:
            letter = [0 for _ in range(len(num_to_AA))]
            onehot_encoded.append(letter)

    return onehot_encoded
def process_file(path, N):
    seq_seq = []
    seqname = []
    all_protein_array = []
    for seq in SeqIO.parse(path, 'fasta'):
        seq_seq.append(seq.seq)
        seqname.append(seq.id)
    for protein in seq_seq:
        protein = str(protein)

        if len(protein) > N:
            protein = protein[-N:]
        else:
            for i in range(len(protein), N):

                protein = protein + 'X'

        seq = encode(protein)
        all_protein_array.append(seq)
    all_protein_array = np.array(all_protein_array)
    return seqname,all_protein_array



