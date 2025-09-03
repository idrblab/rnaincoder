"""
theme: psssa, pss, psa
data:2020.11.03
author:xiaweiqi
"""

import pickle as pkl
import numpy as np
from Bio import SeqIO
_SecondaryStr={'1':'EALMQKRH','2':'VIYCWFT','3':'GNPSD'}
#'1'stand for Helix; '2'stand for Strand, '3' stand for coil
_SolventAccessibility={'1':'ALFCGIVW','2':'RKQEND'}
#'1'stand for Buried; '2'stand for Exposed, '3' stand for Intermediate

HSC = {}
for k, v in _SecondaryStr.items():
	for j in list(v):
		letter = [0, 0, 0]
		k = int(k)
		letter[k-1] = 1
		HSC[j] = letter


BE = {}
for k, v in _SolventAccessibility.items():
	for j in list(v):
		letter = [0, 0]
		k = int(k)
		letter[k-1] = 1
		BE[j] = letter
k = 'MPSTHY'
for i in k:
	BE[i] = [0, 0]


def encode(protein_seq):
	PSSSA = []
	PSS = []
	PSA = []
	for s in protein_seq:
		if s in BE:
			be = BE[s]
		else:
			be = [0, 0]
		if s in HSC:
			hse = HSC[s]
		else:
			hse = [0, 0, 0]
		hse_be = hse + be

		PSSSA.append(hse_be)
		PSA.append(be)
		PSS.append(hse)
	return PSSSA, PSA, PSS

def process_psssa(path, N):
	PSSSA_all, PSA_all, PSS_all = [], [], []
	seq_seq = []
	seqname = []
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

		PSSSA, PSA, PSS = encode(protein)
		PSSSA_all.append(PSSSA)
		PSA_all.append(PSA)
		PSS_all.append(PSS)
	PSSSA_all = np.array(PSSSA_all)
	return seqname,PSSSA_all
