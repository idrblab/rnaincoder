# -*- coding: utf-8 -*-
"""
theme: one_hot
author: xiaweiqi
data: 2020.11.03
"""

import re
import numpy as np

from Bio import SeqIO
from pandas import DataFrame
AALetter=["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
#############################################################################################
def CalculateAAComposition(ProteinSequence):

	LengthSequence=len(ProteinSequence)
	Result={}
	for i in AALetter:
		Result[i]=round(float(ProteinSequence.count(i))/LengthSequence*100,3)
	return Result
def convert_to_array(protein_dict, N):
    protein_array = np.zeros(N, dtype = np.float)
    for i, (k, v) in enumerate(protein_dict.items()):
        protein_array[i] = v
    return protein_array
def process_fasta(filepath):
    AALetter=['Amino acid A', 'Amino acid R', 'Amino acid N', 'Amino acid D', 'Amino acid C', 'Amino acid E', 'Amino acid Q', 'Amino acid G', 'Amino acid H', 'Amino acid I', 'Amino acid L', 'Amino acid K', 'Amino acid M', 'Amino acid F', 'Amino acid P', 'Amino acid S', 'Amino acid T', 'Amino acid W', 'Amino acid Y', 'Amino acid V']

    seqname = []
    seq_seq = []
    all_protein_array = []
    for seq in SeqIO.parse(filepath,'fasta'):
        seqid = seq.id
        seqname.append(seqid)
        seq_seq.append(seq.seq)
    for protein in seq_seq:
        result = CalculateAAComposition(protein)
        N = len(result)
        protein_array = convert_to_array(result, N)
        all_protein_array.append(protein_array)

    df_cod = DataFrame(data=all_protein_array, index=seqname, columns=AALetter)

    return df_cod


