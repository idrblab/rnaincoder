import time
import os
import numpy as np
from Bio import SeqIO

import methods._p01_CTD2 as CTD2
import methods._p02_AAC as AAC
import methods._p03_onehot as one_hot
import methods._p04_PSSSA as PSSSA
import methods._p051_mk_PSSM as PSSM
import methods._p052_simplePSSM as sPSSM

dictMe = {
    'Amino acid composition (1D)': 'p1',
    'Position specific scoring (2D)': 'p4',
    'One-hot encoding (2D)':'p3',

    "Electric charge based (1D)": "p2_01",
    "Hydrophobicity based (1D)": "p2_02",
    "Polarity based (1D)": "p2_03",
    "Polarizability based (1D)": "p2_04",
    "Solvent accessibility based (1D)": "p2_05",
    "Surface tension based (1D)": "p2_06",
    "Van der waals volume (1D)": "p2_07",

    "Structural CTD based (1D)": "p2_08",
    'Solvent accessibility (2D)':'p5'
}


def switch_prometh(fun, textPath):
    if fun == 'p1':
        # return AAC.process_fasta(textPath)
        return CTD2.process_fasta(textPath, 9)
    elif fun == 'p2_01':
        return CTD2.process_fasta(textPath,1)
    elif fun == 'p2_02':
        return CTD2.process_fasta(textPath,2)
    elif fun == 'p2_03':
        return CTD2.process_fasta(textPath,3)
    elif fun == 'p2_04':
        return CTD2.process_fasta(textPath,4)
    elif fun == 'p2_05':
        return CTD2.process_fasta(textPath,5)
    elif fun == 'p2_06':
        return CTD2.process_fasta(textPath,6)
    elif fun == 'p2_07':
        return CTD2.process_fasta(textPath,7)
    elif fun == 'p2_08':
        return CTD2.process_fasta(textPath,8)
    elif fun == 'p3':
        seqname, narrydata = one_hot.process_file(textPath, 1000)
        return seqname, narrydata
    elif fun == 'p4':
        #####暂时模板
        PSSM_sample = np.load(os.path.join(os.path.dirname(__file__), 'Data/Protein_pssm.npy'))
        seq_nameall = []
        for seq in SeqIO.parse(textPath, 'fasta'):
            seq_name = seq.id
            seq_nameall.append(seq_name)
        length = len(seq_nameall)
        PSSM_sample01 = PSSM_sample[:length, :,:]

        return seq_nameall,PSSM_sample01


    elif fun == 'p5':
        seqname, narrydata = PSSSA.process_psssa(textPath, 1000)
        return seqname, narrydata

    else:
        return None