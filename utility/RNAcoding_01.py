#!/usr/bin/env Python
# coding=utf-8
import os
import os.path
import sys
import pandas as pd
import numpy as np
import itertools
import Bio.SeqIO as Seq
import argparse as agp

sys.path.append("..")

import methods.Methods_all_16_methods as Methods_all
import methods.proteinMeth as proteinMeth
import methods.MoleculerMeths as MoleculerMeths

############## 常用函数 ################################
def press_oneRNA(methodsAs, Afilepath):
    # allmethod = len(methodsAs) * 2
    if len(methodsAs) == 1:
        for method in methodsAs:
            Aresult = Methods_all.switch_meth(Methods_all.dictMe[method], Afilepath)
            Aresult = round(Aresult.iloc[:, :], 6)

    else:
        methodsA1 = methodsAs[0]
        Aresult = Methods_all.switch_meth(Methods_all.dictMe[methodsA1], Afilepath)
        Aresult = round(Aresult.iloc[:, :], 6)
        for i in range(1, len(methodsAs)):
            result_n = Methods_all.switch_meth(Methods_all.dictMe[methodsAs[i]], Afilepath)
            Aresult = pd.concat([Aresult, result_n], axis=1, join='inner')
            Aresult = round(Aresult.iloc[:, :], 6)
    return Aresult

def press_twoRNA(methodsAs, Afilepath):
    # allmethod = len(methodsAs) * 2
    if len(methodsAs) == 1:
        for method in methodsAs:
            seqname2D, Aresult = Methods_all.switch_meth(Methods_all.dictMe[method], Afilepath)
            Aresult = np.around(Aresult, decimals=6)
    else:
        methodsA1 = methodsAs[0]
        seqname2D, Aresult = Methods_all.switch_meth(Methods_all.dictMe[methodsA1], Afilepath)
        Aresult = np.around(Aresult, decimals=6)
        # print(Aresult.shape)
        for i in range(1, len(methodsAs)):
            seqname2D, result_n = Methods_all.switch_meth(Methods_all.dictMe[methodsAs[i]], Afilepath)
            # print(result_n.shape)
            Aresult = np.concatenate((Aresult, result_n), axis=2)
            Aresult = np.around(Aresult, decimals=6)
    return seqname2D, Aresult

def press_oneProtein(methodsAs, Afilepath):
    if len(methodsAs) == 1:
        for method in methodsAs:
            Aresult = proteinMeth.switch_prometh(proteinMeth.dictMe[method], Afilepath)
            Aresult = round(Aresult.iloc[:, :], 6)

    else:
        methodsA1 = methodsAs[0]
        Aresult = proteinMeth.switch_prometh(proteinMeth.dictMe[methodsA1], Afilepath)
        Aresult = round(Aresult.iloc[:, :], 6)
        for i in range(1, len(methodsAs)):
            result_n = proteinMeth.switch_prometh(proteinMeth.dictMe[methodsAs[i]], Afilepath)
            # print(result_n)
            # result_n.to_csv('/public/home/wangyx/LncRNA/smallRNA/corain_standlone/output/result_n.csv')
            # print(Aresult)
            # Aresult.to_csv('/public/home/wangyx/LncRNA/smallRNA/corain_standlone/output/Aresult.csv')
            Aresult = pd.concat([Aresult, result_n], axis=1, join='inner')
            Aresult = round(Aresult.iloc[:, :], 6)
    return Aresult

def press_twoProtein(methodsAs, Afilepath):
    allmethod = len(methodsAs) * 2
    if len(methodsAs) == 1:
        seqname2D, Aresult = proteinMeth.switch_prometh(proteinMeth.dictMe[methodsAs[0]], Afilepath)
        Aresult = np.around(Aresult, decimals=6)
        # print(methodsAs[0])
        # print(Aresult.shape)
    else:
        methodsA1 = methodsAs[0]
        # print(methodsA1)
        seqname2D, Aresult = proteinMeth.switch_prometh(proteinMeth.dictMe[methodsA1], Afilepath)
        Aresult = np.around(Aresult, decimals=6)
        # print(Aresult.shape)
        for i in range(1, len(methodsAs)):
            # print(methodsAs[i])
            seqname2D, result_n = proteinMeth.switch_prometh(proteinMeth.dictMe[methodsAs[i]], Afilepath)
            # print(result_n.shape)
            Aresult = np.concatenate((Aresult, result_n), axis=2)
            # print(Aresult.shape)
            Aresult = np.around(Aresult, decimals=6)
    return seqname2D, Aresult

def remove_uncorrectchacter(infasta):
    seqnames = []
    sequences = []
    dict_data = {}
    for seq in Seq.parse(infasta, 'fasta'):        
        seqid = seq.id
        # print(len(seqid))
        if len(seqid) == 0:
            continue
        else:
            seqid = '>' + seqid
            n = 0
            sequence_o = str(seq.seq)
            sequence_o01 = sequence_o.upper()
            sequence_r = sequence_o01.replace("U", "T")
            charterlist = ['B',  'D', 'E', 'F',  'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S',  'U', 'V', 'W', 'X', 'Y', 'Z']
            for unchar in charterlist:
                # print(unchar)
                # print(type(unchar))
                sequence_r = sequence_r.replace(unchar, "")
            seqnames.append(seqid)
            dict_data[seqid] = sequence_r
            sequences.append(sequence_r)

    
    name_seq = dict(zip(seqnames, sequences))
    
    A_save = open(infasta, 'w')
    for seqname in list(dict_data.keys()):
        a_seq = dict_data[seqname]
        if a_seq != '':
            A_save.write(str(seqname) + '\n' + str(a_seq) + '\n')
    A_save.close()

def remove_uncorrectchacter_protein(infasta):
    seqnames = []
    sequences = []
    dict_data = {}
    for seq in Seq.parse(infasta, 'fasta'):        
        seqid = seq.id
        # print(len(seqid))
        if len(seqid) == 0:
            continue
        else:
            seqid = '>' + seqid
            n = 0
            sequence_o = str(seq.seq)
            seqnames.append(seqid)
            dict_data[seqid] = sequence_o
            sequences.append(sequence_o)

    A_save = open(infasta, 'w')
    for seqname in list(dict_data.keys()):
        a_seq = dict_data[seqname]
        if a_seq != '':
            A_save.write(str(seqname) + '\n' + str(a_seq) + '\n')
    A_save.close()

def mk_numpy(Aresult_F_2D,seqname2D_A,nameA):
    dictdata = dict(zip(seqname2D_A,Aresult_F_2D.tolist()))
    datalist = []
    for indexname in nameA.tolist():
        datalist.append(dictdata[indexname])
    A_res = np.array(datalist).reshape(-1, 1000, len(datalist[0][0]))
    return A_res


def make_comcod(com_n):
    # make the combination list of coding methods
    methods_1Ds  = ['Open reading frame (1D)', 'Entropy density of transcript (1D)', 'Global descriptor (1D)', 'K-mer (1D)', 'Codon related (1D)', 'Pseudo protein related (1D)', 'Guanine-cytosine related (1D)', 'Nucleotide related (1D)', 'Secondary structure (1D)', 'EIIP based spectrum (1D)']
    if com_n:
        combination_methods = list(itertools.combinations(methods_1Ds, com_n))
    else:
        for num in range(len(methods_1Ds)):
            com_temp = list(itertools.combinations(methods_1Ds, num+1))
            if num == 0:
                combination_methods = com_temp
            else:
                combination_methods = combination_methods + com_temp

    return combination_methods






# RNA coding analysis---------------------------------------------------------------------------
def RNA_coding(Afastapath, Interfilepath, Resultpath, dimension = '1', savetype = 'npycsv',n_select = None):
    '''
    --n_select:
    '1', 'Open reading frame (1D)'
    '2', 'Entropy density of transcript (1D)'
    '3', 'Global descriptor (1D)'
    '4', 'K-mer (1D)'
    '5', 'Codon related (1D)'
    '6', 'Pseudo protein related (1D)'
    '7', 'Guanine-cytosine related (1D)'
    '8', 'Nucleotide related (1D)'
    '9', 'Secondary structure (1D)'
    '10', 'EIIP based spectrum (1D)'
    '''
    ######### get the methods of file A ------------------------------------------------------------------------------
    #all methods
    methodsAs00s_all = list(Methods_all.dictMe.keys())

    methodsAs_1D = []
    methodsAs_2D = []

    #RNA
    for methodsAs00 in methodsAs00s_all:
        if methodsAs00[-3:-1] == '1D':
            methodsAs_1D.append(methodsAs00)
        elif methodsAs00[-3:-1] == '2D':
            methodsAs_2D.append(methodsAs00)

    methodsAs_1D  = ['Open reading frame (1D)', 'Entropy density of transcript (1D)', 'Global descriptor (1D)', 'K-mer (1D)', 'Codon related (1D)', 'Pseudo protein related (1D)', 'Guanine-cytosine related (1D)', 'Nucleotide related (1D)', 'Secondary structure (1D)', 'EIIP based spectrum (1D)','Solubility lipoaffinity (1D)','Partition coefficient (1D)','Polarizability refractivity (1D)','Hydrogen bond related (1D)','Topological indice (1D)','Molecular fingerprint (1D)']
    methodsAs_2D = []
    # print('methods: ')
    # print(methodsAs_2D)

    Resultpath = Resultpath + '/encoding_features'
    os.makedirs(Resultpath, exist_ok= True)

    remove_uncorrectchacter(Afastapath)
    print('remove uncorrectchacter is right')

    # RNA
    if dimension == '1':
        for n in range(1, 2, 1):
            methodsAs_1Ds = list(itertools.combinations(methodsAs_1D, n))
            
            if n_select:
                n_select = int(n_select)
                n_len_be = n_select-1
                n_len_af = n_select
            else:
                n_len_be = 0
                n_len_af = len(methodsAs_1Ds)
            # print('n_len_af')
            # print(n_len_af)
            for num in range(n_len_be, n_len_af, 1):
                methodsAs_1D01 = methodsAs_1Ds[num]
                print('The %s method is starting!' % num)
                print(list(methodsAs_1D01))
                methodsAs_1D02 = list(methodsAs_1D01)
                print('dimension %s !' % dimension)
            
                #### one demensional methods
                # print('The method is starting!')
                Aresult_F_1D = press_oneRNA(methodsAs_1D02, Afastapath)
                
                # print('Aresult_F_1D.head()')
                # print(Aresult_F_1D.head())
                trainval_seq_data = pd.read_csv(Interfilepath)
                trainval_seq_data01  = trainval_seq_data.replace('>', '',regex =True)
                tranval_A = trainval_seq_data01.iloc[:,0]
                # print('tranval_A')
                # print(tranval_A)

                A_fea = pd.merge(tranval_A, Aresult_F_1D, left_on='Seqname', right_index=True, how='left', sort=False)
                # replace NA with 0
                A_fea = A_fea.fillna(0)
                # print('A_fea.head()')
                # print(A_fea.head())
                # combine the file
                A_res = np.array(A_fea.iloc[:,1:], np.float64)
                # A_res = A_res.astype(np.float64)
                # print('A_res.shape')
                # print(A_res.shape)
                # print('A_res.dtype')
                # print(A_res.dtype)

                methodsAs_1D01_01 = str(methodsAs_1D01).split("'")[1]
                # print(savetype)
                if 'npy' in savetype:
                    FilepathA = os.path.join(Resultpath, str(methodsAs_1D01_01) + '.npy')
                    np.save(FilepathA, A_res)
                if 'csv' in savetype:
                    A_fea.to_csv(os.path.join(Resultpath, str(methodsAs_1D01_01) + '.csv'))
                print('The %s method is ending!' % num)
    elif dimension == '2':
        # two demensional methods
        for n in range(1, 2, 1):
            methodsAs_2Ds = list(itertools.combinations(methodsAs_2D, n))
            for num in range(0, len(methodsAs_2Ds), 1):
                methodsAs_2D01 = methodsAs_2Ds[num]
                print('The %s method is starting!' % num)
                print(list(methodsAs_2D01))
                methodsAs_2D02 = list(methodsAs_2D01)

                seqname2D_A, Aresult_F_2D = press_twoRNA(methodsAs_2D02, Afastapath)
                A_res = Aresult_F_2D
        
                methodsAs_2D01_01 = str(methodsAs_2D01).split("'")[1]
                FilepathA = os.path.join(Resultpath, str(methodsAs_2D01_01) + '.npy')
                np.save(FilepathA, A_res)
                print('The %s method is ending!' % num)
# RNA RNA coding analysis---------------------------------------------------------------------------
def RNA_RNA_coding(Afastapath,Bfastapath,Interfilepath,Resultpath,dimension,savetype,n_select =None):
    #########------------------------------------------------------------------------------

    methodsAs00s_all = list(Methods_all.dictMe.keys())

    methodsAs_1D = []
    methodsAs_2D = []
    
    #RNA
    for methodsAs00 in methodsAs00s_all:
        if methodsAs00[-3:-1] == '1D':
            methodsAs_1D.append(methodsAs00)
        elif methodsAs00[-3:-1] == '2D':
            methodsAs_2D.append(methodsAs00)
    # print('protein methods: ')
    # print(methodsAs_1D)
    # print(methodsAs_2D)
    methodsAs_1D  = ['Open reading frame (1D)', 'Entropy density of transcript (1D)', 'Global descriptor (1D)', 'K-mer (1D)', 'Codon related (1D)', 'Pseudo protein related (1D)', 'Guanine-cytosine related (1D)', 'Nucleotide related (1D)', 'Secondary structure (1D)', 'EIIP based spectrum (1D)','Solubility lipoaffinity (1D)','Partition coefficient (1D)','Polarizability refractivity (1D)','Hydrogen bond related (1D)','Topological indice (1D)','Molecular fingerprint (1D)']

    Resultpath = Resultpath + '/encoding_features'
    os.makedirs(Resultpath, exist_ok= True)

    remove_uncorrectchacter(Afastapath)
    remove_uncorrectchacter(Bfastapath)
    print('remove_uncorrectchacter is right')
    # get label file
    trainval_seq_data = pd.read_csv(Interfilepath)
    trainval_seq_data01  = trainval_seq_data.replace('>', '',regex =True)
    tranval_A = trainval_seq_data01['A']
    tranval_B = trainval_seq_data01['B']
    # print('tranval_A')
    # print(tranval_A)

    # RNA
    for n in range(1, 2, 1):
        methodsAs_1Ds = list(itertools.combinations(methodsAs_1D, n))
        
        if n_select:
            n_select = int(n_select)
            n_len_be = n_select-1
            n_len_af = n_select
        else:
            n_len_be = 0
            n_len_af = len(methodsAs_1Ds)
        print('n_len_af')
        print(n_len_af)
        for num in range(n_len_be, n_len_af, 1):
            methodsAs_1D01 = methodsAs_1Ds[num]
            print('The %s method is starting!' % num)
            print(list(methodsAs_1D01))
            methodsAs_1D02 = list(methodsAs_1D01)
            if dimension == '1':
                ####
                Aresult_F_1D = press_oneRNA(methodsAs_1D02, Afastapath)
                Bresult_F_1D = press_oneRNA(methodsAs_1D02, Bfastapath)

                A_fea = pd.merge(tranval_A, Aresult_F_1D, left_on='A', right_index=True, how='left', sort=False)
                B_fea = pd.merge(tranval_B, Bresult_F_1D, left_on='B', right_index=True, how='left', sort=False)
                # replace NA with 0
                try:
                    A_fea = A_fea.fillna(0)
                    B_fea = B_fea.fillna(0)
                except Exception as e:
                    continue
                # 
                A_res = np.array(A_fea.iloc[:,1:], np.float64)
                # A_res = A_res.astype(np.float64)
                B_res = np.array(B_fea.iloc[:,1:], np.float64)

                methodsAs_1D01_01 = str(methodsAs_1D01).split("'")[1]
                if 'npy' in savetype:
                    FilepathA = os.path.join(Resultpath, str(methodsAs_1D01_01) + '_A.npy')
                    np.save(FilepathA, A_res)

                    FilepathB = os.path.join(Resultpath, str(methodsAs_1D01_01) + '_B.npy')
                    np.save(FilepathB, B_res)
                if 'csv' in savetype:
                    A_fea.to_csv(os.path.join(Resultpath, str(methodsAs_1D01_01) + '_A.csv'))
                    B_fea.to_csv(os.path.join(Resultpath, str(methodsAs_1D01_01) + '_B.csv'))
                print('The %s method is ending!' % num)

    if dimension == '2':
        # 
        nameA =np.array(tranval_A)
        nameB =np.array(tranval_B)
        for n in range(1, 2, 1):
            methodsAs_2Ds = list(itertools.combinations(methodsAs_2D, n))
            for num in range(0, len(methodsAs_2Ds), 1):
                methodsAs_2D01 = methodsAs_2Ds[num]
                print('The %s method is starting!' % num)
                print(list(methodsAs_2D01))
                methodsAs_2D02 = list(methodsAs_2D01)
                ####
                seqname2D_A, Aresult_F_2D = press_twoRNA(methodsAs_2D02, Afastapath)
                seqname2D_B, Bresult_F_2D = press_twoRNA(methodsAs_2D02, Bfastapath)
                print(nameA)
                print(seqname2D_A)
                A_res = mk_numpy(Aresult_F_2D,seqname2D_A,nameA)
                B_res = mk_numpy(Bresult_F_2D,seqname2D_B,nameB)           
       
                methodsAs_2D01_01 = str(methodsAs_2D01).split("'")[1]
                FilepathA = os.path.join(Resultpath, str(methodsAs_2D01_01) + '_A.npy')
                np.save(FilepathA, A_res)
                FilepathB = os.path.join(Resultpath, str(methodsAs_2D01_01) + '_B.npy')
                np.save(FilepathB, B_res)
                print('The %s method is ending!' % num)
# RNA protein coding analysis---------------------------------------------------------------------------
def RNA_protein_coding(Afastapath,Bfastapath,Interfilepath,Resultpath,dimension,savetype,n_select =None):
    ######### ------------------------------------------------------------------------------
    #
    methodsAs00s_all = list(Methods_all.dictMe.keys())
    methodsBs00s_all = list(proteinMeth.dictMe.keys())


    methodsAs_1D = []
    methodsAs_2D = []
    #RNA
    for methodsAs00 in methodsAs00s_all:
        if methodsAs00[-3:-1] == '1D':
            methodsAs_1D.append(methodsAs00)
        elif methodsAs00[-3:-1] == '2D':
            methodsAs_2D.append(methodsAs00)
    print('methods: ')
    print(methodsAs_1D)
    print(methodsAs_2D)
    # 
    methodsBs_1D = []
    methodsBs_2D = []
    #protein
    for methodsBs00 in methodsBs00s_all:
        if methodsBs00[-3:-1] == '1D':
            methodsBs_1D.append(methodsBs00)
        elif methodsBs00[-3:-1] == '2D':
            methodsBs_2D.append(methodsBs00)
    print('protein methods: ')
    print(methodsBs_1D)
    print(methodsBs_2D)    

    # 
    Resultpath = Resultpath + '/encoding_features'
    os.makedirs(Resultpath, exist_ok= True)

    remove_uncorrectchacter(Afastapath)
    remove_uncorrectchacter_protein(Bfastapath)
    print('remove_uncorrectchacter is right')
    # get label file
    trainval_seq_data = pd.read_csv(Interfilepath)
    trainval_seq_data01  = trainval_seq_data.replace('>', '',regex =True)
    tranval_A = trainval_seq_data01['A']
    tranval_B = trainval_seq_data01['B']
    # print('tranval_A')
    # print(tranval_A)

    # RNA
    if dimension == '1':
        # encoding protein
        Bresult_F_1D = press_oneProtein(methodsBs_1D, Bfastapath)
        B_fea = pd.merge(tranval_B, Bresult_F_1D, left_on='B', right_index=True, how='left', sort=False)
        B_res = np.array(B_fea.iloc[:,1:], np.float64)
        # print('B_fea.head()')
        # print(B_fea.head())                        
        # print('B_res.shape')
        # print(B_res.shape)

        if 'npy' in savetype:
            FilepathB = os.path.join(Resultpath, 'protein_B.npy')
            np.save(FilepathB, B_res)
        if 'csv' in savetype:
            B_fea.to_csv(os.path.join(Resultpath, 'protein_B.csv'))
        # encoding RNA
        for n in range(1, 2, 1):
            methodsAs_1Ds = list(itertools.combinations(methodsAs_1D, n))
            
            if n_select:
                n_select = int(n_select)
                n_len_be = n_select-1
                n_len_af = n_select
            else:
                n_len_be = 0
                n_len_af = len(methodsAs_1Ds)
            print('n_len_af')
            print(n_len_af)

            for num in range(n_len_be, n_len_af, 1):
                methodsAs_1D01 = methodsAs_1Ds[num]
                print('The %s method is starting!' % num)
                print(list(methodsAs_1D01))
                methodsAs_1D02 = list(methodsAs_1D01)
                
                ###
                Aresult_F_1D = press_oneRNA(methodsAs_1D02, Afastapath)               
                
                # print('Aresult_F_1D.head()')
                # print(Aresult_F_1D.head())

                A_fea = pd.merge(tranval_A, Aresult_F_1D, left_on='A', right_index=True, how='left', sort=False)
                # replace NA with 0
                try:
                    A_fea = A_fea.fillna(0)
                except Exception as e:
                    continue
                # print('A_fea.head()')
                # print(A_fea.head())

                A_res = np.array(A_fea.iloc[:,1:], np.float64)
                # A_res = A_res.astype(np.float64)
                # print('A_res.shape')
                # print(A_res.shape)

                methodsAs_1D01_01 = str(methodsAs_1D01).split("'")[1]
                if 'npy' in savetype:
                    FilepathA = os.path.join(Resultpath, str(methodsAs_1D01_01) + '_A.npy')
                    np.save(FilepathA, A_res)
                if 'csv' in savetype:
                    A_fea.to_csv(os.path.join(Resultpath, str(methodsAs_1D01_01) + '_A.csv'))
                print('The %s method is ending!' % num)

    if dimension == '2':
        # 
        nameB =np.array(tranval_B)
        seqname2D_B, Bresult_F_2D = press_twoProtein(methodsBs_2D, Bfastapath)
        B_res = mk_numpy(Bresult_F_2D,seqname2D_B,nameB)
        FilepathB = os.path.join(Resultpath, 'protein_B.npy')
        np.save(FilepathB, B_res)

        # 
        nameA =np.array(tranval_A)        
        for n in range(1, 2, 1):
            methodsAs_2Ds = list(itertools.combinations(methodsAs_2D, n))
            for num in range(0, len(methodsAs_2Ds), 1):
                methodsAs_2D01 = methodsAs_2Ds[num]
                print('The %s method is starting!' % num)
                print(list(methodsAs_2D01))
                methodsAs_2D02 = list(methodsAs_2D01)
                ###
                seqname2D_A, Aresult_F_2D = press_twoRNA(methodsAs_2D02, Afastapath)                
                print(nameA)
                print(seqname2D_A)
                A_res = mk_numpy(Aresult_F_2D,seqname2D_A,nameA)                    
                methodsAs_2D01_01 = str(methodsAs_2D01).split("'")[1]
                FilepathA = os.path.join(Resultpath, str(methodsAs_2D01_01) + '_A.npy')
                np.save(FilepathA, A_res)
                print('The %s method is ending!' % num)
# RNA compound coding analysis---------------------------------------------------------------------------
def RNA_compound_coding(Afastapath,Bfastapath,Interfilepath,Resultpath,dimension,savetype,n_select =None):
    #########  ------------------------------------------------------------------------------
    #
    methodsAs00s_all = list(Methods_all.dictMe.keys())
    methodsBs00s_all = list(MoleculerMeths.dictMe.keys())


    methodsAs_1D = []
    methodsAs_2D = []
    #RNA
    for methodsAs00 in methodsAs00s_all:
        if methodsAs00[-3:-1] == '1D':
            methodsAs_1D.append(methodsAs00)
        elif methodsAs00[-3:-1] == '2D':
            methodsAs_2D.append(methodsAs00)
    # print('methods: ')
    # print(methodsAs_1D)
    # print(methodsAs_2D)

    Resultpath = Resultpath + '/encoding_features'
    os.makedirs(Resultpath, exist_ok= True)

    remove_uncorrectchacter(Afastapath)
    print('remove_uncorrectchacter is right')
    # get label file
    trainval_seq_data = pd.read_csv(Interfilepath)
    trainval_seq_data01  = trainval_seq_data.replace('>', '',regex =True)
    tranval_A = trainval_seq_data01['A']
    tranval_B = trainval_seq_data01['B']
    # print('tranval_A')
    # print(tranval_A)

    # encoding compound
    ResfilePath = os.path.join(Resultpath, 'MoleculerResult.csv')
    Bresult_F_1D = MoleculerMeths.Moleculer_methods(Bfastapath,ResfilePath,Interfilepath,methodsBs00s_all)
    # print('tranval_B')
    # print(tranval_B)
    # print('Bresult_F_1D')
    # print(Bresult_F_1D)
    # B_fea = pd.merge(tranval_B, Bresult_F_1D, left_on='B', right_index=True, how='left', sort=False)
    B_fea = Bresult_F_1D
    B_res = np.array(B_fea.iloc[:,1:], np.float64)
    # print('B_fea.head()')
    # print(B_fea.head())                        
    # print('B_res.shape')
    # print(B_res.shape)


    # RNA
    if dimension == '1':
        if 'npy' in savetype:
            FilepathB = os.path.join(Resultpath, 'compound_B.npy')
            np.save(FilepathB, B_res)
        if 'csv' in savetype:
            B_fea.to_csv(os.path.join(Resultpath, 'compound_B.csv'))
        # encoding RNA
        for n in range(1, 2, 1):
            methodsAs_1Ds = list(itertools.combinations(methodsAs_1D, n))
            
            if n_select:
                n_select = int(n_select)
                n_len_be = n_select-1
                n_len_af = n_select
            else:
                n_len_be = 0
                n_len_af = len(methodsAs_1Ds)
            print('n_len_af')
            print(n_len_af)

            for num in range(n_len_be, n_len_af, 1):
                methodsAs_1D01 = methodsAs_1Ds[num]
                print('The %s method is starting!' % num)
                print(list(methodsAs_1D01))
                methodsAs_1D02 = list(methodsAs_1D01)
                
                ####
                Aresult_F_1D = press_oneRNA(methodsAs_1D02, Afastapath)               
                
                # print('Aresult_F_1D.head()')
                # print(Aresult_F_1D.head())

                A_fea = pd.merge(tranval_A, Aresult_F_1D, left_on='A', right_index=True, how='left', sort=False)
                # replace NA with 0
                try:
                    A_fea = A_fea.fillna(0)
                except Exception as e:
                    continue
                # print('A_fea.head()')
                # print(A_fea.head())

                A_res = np.array(A_fea.iloc[:,1:], np.float64)

                # print('A_res.shape')
                # print(A_res.shape)

                methodsAs_1D01_01 = str(methodsAs_1D01).split("'")[1]
                if 'npy' in savetype:
                    FilepathA = os.path.join(Resultpath, str(methodsAs_1D01_01) + '_A.npy')
                    np.save(FilepathA, A_res)
                if 'csv' in savetype:
                    A_fea.to_csv(os.path.join(Resultpath, str(methodsAs_1D01_01) + '_A.csv'))
                print('The %s method is ending!' % num)

    if dimension == '2':
        # 
        nameA =np.array(tranval_A)        
        for n in range(1, 2, 1):
            methodsAs_2Ds = list(itertools.combinations(methodsAs_2D, n))
            for num in range(0, len(methodsAs_2Ds), 1):
                methodsAs_2D01 = methodsAs_2Ds[num]
                print('The %s method is starting!' % num)
                print(list(methodsAs_2D01))
                methodsAs_2D02 = list(methodsAs_2D01)
                ####
                seqname2D_A, Aresult_F_2D = press_twoRNA(methodsAs_2D02, Afastapath)                
                print(nameA)
                print(seqname2D_A)
                A_res = mk_numpy(Aresult_F_2D,seqname2D_A,nameA)                    
                methodsAs_2D01_01 = str(methodsAs_2D01).split("'")[1]
                FilepathA = os.path.join(Resultpath, str(methodsAs_2D01_01) + '_A.npy')
                np.save(FilepathA, A_res)
                print('The %s method is ending!' % num)

