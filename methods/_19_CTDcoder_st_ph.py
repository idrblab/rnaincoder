import sys
import re
import numpy as np
from pandas import DataFrame
from Bio import SeqIO  # https://blog.csdn.net/weixin_30600197/article/details/97106545
import Methods_all_16_methods as Methods_all_16_methods


class CTDcoder:
    """Generates CTD counts for a fasta file"""

    def __init__(self, infasta,ACGT_encode):

        # input fasta_file
        self.infasta = infasta

        # input feature_encode
        self.ACGT_encode = ACGT_encode
        # self.ACGT_encode = ['0111', '0100', '0110', '1000', '1010', '0100',]

    def CTD(self, fas, encode='0101'):

        # one-hot encoding
        seq = str(fas.seq)
        seq = seq.replace('A', encode[0]).replace('C', encode[1]).replace('G', encode[2]).replace('T', encode[3])

        # X stands for 0, Y stands for 1.
        num_X = seq.count('0')
        num_Y = seq.count('1')
        XY_trans = seq.count('01') + seq.count('10')

        n = len(seq)

        # distributions
        x,y=0,0
        X0_dis,X1_dis,X2_dis,X3_dis,X4_dis=0.0,0.0,0.0,0.0,0.0
        Y0_dis,Y1_dis,Y2_dis,Y3_dis,Y4_dis=0.0,0.0,0.0,0.0,0.0


        for i in range(len(seq)):
            if seq[i]=="0":
                x=x+1
                if x == 1:
                    X0_dis=((i*1.0)+1)/n
                if x == int(round(num_X/4.0)):
                    X1_dis=((i*1.0)+1)/n
                if x == int(round(num_X/2.0)):
                    X2_dis=((i*1.0)+1)/n
                if x == int(round((num_X*3/4.0))):
                    X3_dis=((i*1.0)+1)/n
                if x == num_X:
                    X4_dis=((i*1.0)+1)/n
            if seq[i]=="1":
                y=y+1
                if y == 1:
                    Y0_dis=((i*1.0)+1)/n
                if y == int(round(num_Y/4.0)):
                    Y1_dis=((i*1.0)+1)/n
                if y == int(round((num_Y/2.0))):
                    Y2_dis=((i*1.0)+1)/n
                if y == int(round((num_Y*3/4.0))):
                    Y3_dis=((i*1.0)+1)/n
                if y == num_Y:
                    Y4_dis=((i*1.0)+1)/n

        return (list(map(float, [num_X / n, num_Y / n, XY_trans / (n - 1),
                               X0_dis, X1_dis, X2_dis, X3_dis, X4_dis,
                               Y0_dis, Y1_dis, Y2_dis, Y3_dis, Y4_dis])))

    def get_ctd(self):

        # feaname = feature_num * ["X", "X", "XY", "X0", "X1", "X2", "X3", "X4", "Y0", "Y1", "Y2", "Y3", "Y4"]
        # feaname = [str(i + 1) + '_' + str(j + 1) for i in range(len(self.ACGT_encode)) for j in range(13)]
        dictname = {
            '0010': 'Strong H-Bond donors',
            '0110': 'Linear free energy',
            '0101': 'Molar refractivity',
            '1000': 'Lipoaffinity index',
            '0100': 'Gas-hexadecane PC', 
            '0011': 'NH- count',
            '0001': 'Primary or secondary nitrogens'}
        dictname_sh = {'0010': 'SHD', '0110': 'MLF',
                    '0101': 'MRe', 
                    '1000': 'LFI', '0100': 'HPC',
                    '0011': 'CNH',
                '0001': 'PSN'}

        # SEQ = ['01','02','03','04','05','06','07','08','09','10','11','12','13']
        SEQ_sh = ['CA','CB','AB','A0','A1','A2','A3','A4','B0','B1','B2','B3','B4']
        SEQ = {'CA': ' composition of A',
                'CB': ' composition of B',
                'AB': ' transition between A and B',
                'A0': ' distribution of 0.00A',
                'A1': ' distribution of 0.25A',
                'A2': ' distribution of 0.50A',
                'A3': ' distribution of 0.75A',
                'A4': ' distribution of 1.00A',
                'B0': ' distribution of 0.00B',
                'B1': ' distribution of 0.25B',
                'B2': ' distribution of 0.50B',
                'B3': ' distribution of 0.75B',
                'B4': ' distribution of 1.00B'}
        feaname = [dictname_sh[encode]+ j + ": " + dictname[encode] + SEQ[j] for encode in self.ACGT_encode for j in SEQ_sh]
        seqname = []
        feature = []
        for fas in SeqIO.parse(self.infasta, 'fasta'):
            seqid = fas.id
            seqname.append(seqid)
            for encode in self.ACGT_encode:
                feature += self.CTD(fas=fas, encode=encode)

        np.set_printoptions(precision=6)
        data = np.array(feature).reshape((len(seqname), 13 * len(self.ACGT_encode)))

        df = DataFrame(data=data, index=seqname, columns=feaname)
        new_feaname = [i for i in feaname if i in  Methods_all_16_methods.list_new_featrures]
        df_new = df[new_feaname]
        # df = round(df.iloc[:, :], 6)
        return df_new
        # df = round(df.iloc[:, :], 6)
        # return df

class CTDcoder_3class:
    """Generates CTD counts for a fasta file"""

    def __init__(self, infasta,ACGT_encode):

        # input fasta_file
        self.infasta = infasta

        # input feature_encode
        self.ACGT_encode = ACGT_encode
        # self.ACGT_encode = ['0111', '0100', '0110', '1000', '1010', '0100',]

    def CTD(self, fas, encode='0021'):

        # one-hot encoding
        seq = str(fas.seq)
        seq = seq.replace('A', encode[0]).replace('C', encode[1]).replace('G', encode[2]).replace('T', encode[3])

        # X stands for 0, Y stands for 1.
        num_X = seq.count('0')
        num_Y = seq.count('1')
        num_Z = seq.count('2')
        XY_trans = seq.count('01') + seq.count('10')
        XZ_trans = seq.count('02') + seq.count('20')
        YZ_trans = seq.count('12') + seq.count('21')

        n = len(seq)

        # distributions
        x,y,z=0,0,0
        X0_dis,X1_dis,X2_dis,X3_dis,X4_dis=0.0,0.0,0.0,0.0,0.0
        Y0_dis,Y1_dis,Y2_dis,Y3_dis,Y4_dis=0.0,0.0,0.0,0.0,0.0
        Z0_dis,Z1_dis,Z2_dis,Z3_dis,Z4_dis=0.0,0.0,0.0,0.0,0.0


        for i in range(len(seq)):
            if seq[i]=="0":
                x=x+1
                if x == 1:
                    X0_dis=((i*1.0)+1)/n
                if x == int(round(num_X/4.0)):
                    X1_dis=((i*1.0)+1)/n
                if x == int(round(num_X/2.0)):
                    X2_dis=((i*1.0)+1)/n
                if x == int(round((num_X*3/4.0))):
                    X3_dis=((i*1.0)+1)/n
                if x == num_X:
                    X4_dis=((i*1.0)+1)/n
            if seq[i]=="1":
                y=y+1
                if y == 1:
                    Y0_dis=((i*1.0)+1)/n
                if y == int(round(num_Y/4.0)):
                    Y1_dis=((i*1.0)+1)/n
                if y == int(round((num_Y/2.0))):
                    Y2_dis=((i*1.0)+1)/n
                if y == int(round((num_Y*3/4.0))):
                    Y3_dis=((i*1.0)+1)/n
                if y == num_Y:
                    Y4_dis=((i*1.0)+1)/n
            if seq[i]=="2":
                z=z+1
                if z == 1:
                    Z0_dis=((i*1.0)+1)/n
                if z == int(round(num_Z/4.0)):
                    Z1_dis=((i*1.0)+1)/n
                if z == int(round((num_Z/2.0))):
                    Z2_dis=((i*1.0)+1)/n
                if z == int(round((num_Z*3/4.0))):
                    Z3_dis=((i*1.0)+1)/n
                if z == num_Z:
                    Z4_dis=((i*1.0)+1)/n    


        return (list(map(float, [num_X / n, num_Y / n, num_Z / n, 
                XY_trans / (n - 1),XZ_trans / (n - 1),YZ_trans / (n - 1),
                               X0_dis, X1_dis, X2_dis, X3_dis, X4_dis,
                               Y0_dis, Y1_dis, Y2_dis, Y3_dis, Y4_dis,
                               Z0_dis, Z1_dis, Z2_dis, Z3_dis, Z4_dis])))

    def get_ctd(self):

        # feaname = feature_num * ["X", "X", "XY", "X0", "X1", "X2", "X3", "X4", "Y0", "Y1", "Y2", "Y3", "Y4"]
        # feaname = [str(i + 1) + '_' + str(j + 1) for i in range(len(self.ACGT_encode)) for j in range(13)]
        # dictname = {'0010': 'Strong H-Bond donors','0110': 'Linear free energy','0111': 'Molar refractivity','1011': 'Atomic polarizability','0100': 'Solubility','1000': 'Lipoaffinity index','1010': 'Gas-hexadecane PC','0101': 'Quadruple bonds','0011': 'NH- count','0001': 'Hydrogen atoms','1100': 'Secondary carbons','1110': 'Primary or secondary nitrogens'}
        dictname = {'1020': 'Lipoaffinity index_3','0102': 'Gas-hexadecane PC_3', '1200':'Strong H-Bond acceptors_3','0120':'Potential Hydrogen Bonds_3','1002':'Sum of path lengths starting from oxygens_3','0021':'Topological polar surface area_3'}

        dictname_sh = {'1020': 'LFI','0102': 'HPC','1200':'SHA','0120':'PHB','1002':'SLF','0021':'TPS'}


        # SEQ = ['01','02','03','04','05','06','07','08','09','10','11','12','13']
        SEQ_sh = ['CA','CB','CC','AB','AC','BC','A0','A1','A2','A3','A4','B0','B1','B2','B3','B4','C0','C1','C2','C3','C4']
        SEQ = {'CA': ' composition of A',
                'CB': ' composition of B',
                'CC': ' composition of C',
                'AB': ' transition between A and B',
                'AC': ' transition between A and C',
                'BC': ' transition between B and C',
                'A0': ' distribution of 0.00A',
                'A1': ' distribution of 0.25A',
                'A2': ' distribution of 0.50A',
                'A3': ' distribution of 0.75A',
                'A4': ' distribution of 1.00A',
                'B0': ' distribution of 0.00B',
                'B1': ' distribution of 0.25B',
                'B2': ' distribution of 0.50B',
                'B3': ' distribution of 0.75B',
                'B4': ' distribution of 1.00B',
                'C0': ' distriCution of 0.00C',
                'C1': ' distriCution of 0.25C',
                'C2': ' distriCution of 0.50C',
                'C3': ' distriCution of 0.75C',
                'C4': ' distriCution of 1.00C'}
        feaname = [dictname_sh[encode]+ j + ": " + dictname[encode] + SEQ[j] for encode in self.ACGT_encode for j in SEQ_sh]
        seqname = []
        feature = []
        for fas in SeqIO.parse(self.infasta, 'fasta'):
            seqid = fas.id
            seqname.append(seqid)
            for encode in self.ACGT_encode:
                feature += self.CTD(fas=fas, encode=encode)

        np.set_printoptions(precision=6)
        data = np.array(feature).reshape((len(seqname), 21 * len(self.ACGT_encode)))

        df = DataFrame(data=data, index=seqname, columns=feaname)
        new_feaname = [i for i in feaname if i in  Methods_all_16_methods.list_new_featrures]
        df_new = df[new_feaname]
        # df = round(df.iloc[:, :], 6)
        return df_new

# textPath = 'RPI1807_rna_seq.fa'
# # ACGT_encode = ['0001']
# # Partition = CTDcoder(textPath,ACGT_encode).get_ctd()


        
# ACGT_encode_3 = ['0102']
# Partition_3 = CTDcoder_3class(textPath,ACGT_encode_3).get_ctd()
# Partition_3.to_csv('Partition.csv')