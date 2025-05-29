import numpy as np

import methods._01_ORF_code as ORF_code
import methods._02_CTDcode as CTDcode
import methods._03_Fickettcode as Fickettcode
import methods._04_kmer_counts_molmap as kmer_counts
import methods._05_Hexamercode as Hexamercode
import methods._06_proparcoder as proparcoder
import methods._07_GCcounts as GCcounts
import methods._08_edpfeature as edpfeature
import methods._09_StopCodon as StopCodon
import methods._10_onehot as onehot
import methods._11_word2vec as word2vec
import methods._15_SparseEncoding as SparseEncod
import methods._16_RNAFOLD as RNAFold
import methods._18_SStructure as SStructure
import methods._19_CTDcoder_st_ph as CTDcoderstph

import sys
import os
ac_path = os.path.join(os.path.dirname(__file__), 'repDNA')
sys.path.append(ac_path)
import ac as RNA_ac
import psenac as RNA_psenac


import pandas as pd
from Bio import SeqIO

dictMe = {
    'Open reading frame (1D)': '1',
    'Transcript related (1D)': '2',
    'EDP (1D)': '2_1',
    'CTD (1D)': '2_2',
    'Kmers (1D)': '2_3',
    'Codon related (1D)': '3',
    'Pseudo protein related (1D)': '6',
    'Guanine-cytosine related (1D)': '7',
    'Sequence-intrinsic Features:One-hot encoding (2D)': '13',
    # 'Word2vec embeddings': '14',
    # 'Sequence intrinsic features': '13',
    'Sparse encoding (2D)': '15',
    # 'Physicochemical features': '15',
    'Structure-based Features:One-hot encoding (2D)':'16',
    # 'Structural features':'16',
    'Nucleotide related (1D)':'17',
    'Secondary structure (1D)':'18',
    'EIIP based spectrum (1D)':'18_1',
    # 'Solubility and partition coefficient':'19_1',
    'Solubility lipoaffinity (1D)':'19_1',
    'Partition coefficient (1D)':'19_101',
    'Polarizability refractivity (1D)':'19_2',
    'Hydrogen bond related (1D)':'19_3',
    'Topological indice (1D)':'19_4',
    'Molecular fingerprint (1D)':'19_5'
}



def switch_meth(fun, textPath):
    if fun == '1':

        
        ORF = ORF_code.ORF_count(textPath).get_ORF()
        ORFedp = edpfeature.EDPcoder(textPath).getEDP_orf()
        # hexmerORF = SStructure.makeORFEucDist(textPath)

        print(ORF.shape)
        print(ORFedp.shape)
        print(ORF)
        print(ORFedp)
        # ORF.to_csv('/public/home/wangyx/LncRNA/smallRNA/code/ORF.csv')
        # ORFedp.to_csv('/public/home/wangyx/LncRNA/smallRNA/code/ORFedp.csv')

        T1 = pd.concat([ORF, ORFedp], axis=1, join='inner')

        # T1 = pd.concat([T1, hexmerORF], axis=1, join='inner')
        return T1
    elif fun == '2':
        coding_file = os.path.join(os.path.dirname(__file__), 'Data/gencode.v34.pc_transcripts_test.fa')
        # coding_file = '/home/wangyunxia/RNACode/tornadoBulid/LncRNAcoder/methods/Data/gencode.v34.pc_transcripts_test.fa'
        noncoding_file = os.path.join(os.path.dirname(__file__), 'Data/gencode.v34.lncRNA_transcripts_test.fa')
        # Hexamer = Hexamercode.Hexamercoder(textPath, 6, 3, coding_file, noncoding_file).get_hexamer()
        tran_len = edpfeature.EDPcoder(textPath).get_tran_len()
        UTR_len = edpfeature.EDPcoder(textPath).getUTR_len()
        UTR_cov = edpfeature.EDPcoder(textPath).getUTR_cov()
        EDP = edpfeature.EDPcoder(textPath).getEDP()
        CTD = CTDcode.CTDcoder(textPath).get_ctd()
        # Kmer1 = kmer_counts.BasicCounter(textPath, int(1)).get_counts()
        # Kmer2 = kmer_counts.BasicCounter(textPath, int(2)).get_counts()
        # Kmer3 = kmer_counts.BasicCounter(textPath, int(3)).get_counts()    
        # EucLogDist = SStructure.makeEucDist(textPath)

        T1 = pd.concat([UTR_len, tran_len], axis=1, join='inner')
        # T1 = pd.concat([T1, UTR_len], axis=1, join='inner')
        T1 = pd.concat([T1, UTR_cov], axis=1, join='inner')
        T1 = pd.concat([T1, EDP], axis=1, join='inner')
        T1 = pd.concat([T1, CTD], axis=1, join='inner')
        # T1 = pd.concat([T1, Kmer1], axis=1, join='inner')
        # T1 = pd.concat([T1, Kmer2], axis=1, join='inner')
        # T1 = pd.concat([T1, Kmer3], axis=1, join='inner')
        # T1 = pd.concat([T1, EucLogDist], axis=1, join='inner')
        return T1
    elif fun == '2_1':
        # coding_file = os.path.join(os.path.dirname(__file__), 'Data/gencode.v34.pc_transcripts_test.fa')
        # # coding_file = '/home/wangyunxia/RNACode/tornadoBulid/LncRNAcoder/methods/Data/gencode.v34.pc_transcripts_test.fa'
        # noncoding_file = os.path.join(os.path.dirname(__file__), 'Data/gencode.v34.lncRNA_transcripts_test.fa')
        # Hexamer = Hexamercode.Hexamercoder(textPath, 6, 3, coding_file, noncoding_file).get_hexamer()
        # tran_len = edpfeature.EDPcoder(textPath).get_tran_len()
        # # UTR_len = edpfeature.EDPcoder(textPath).getUTR_len()
        # UTR_cov = edpfeature.EDPcoder(textPath).getUTR_cov()
        EDP = edpfeature.EDPcoder(textPath).getEDP()
        
        # CTD = CTDcode.CTDcoder(textPath).get_ctd()
        # Kmer = kmer_counts.BasicCounter(textPath, int(3)).get_counts()
        # EucLogDist = SStructure.makeEucDist(textPath)

        # T1 = pd.concat([Hexamer, tran_len], axis=1, join='inner')
        # # T1 = pd.concat([T1, UTR_len], axis=1, join='inner')
        # T1 = pd.concat([T1, UTR_cov], axis=1, join='inner')
        # T1 = pd.concat([T1, EDP], axis=1, join='inner')
        # T1 = pd.concat([T1, CTD], axis=1, join='inner')
        # T1 = pd.concat([T1, Kmer], axis=1, join='inner')
        # T1 = pd.concat([T1, EucLogDist], axis=1, join='inner')
        return EDP
    elif fun == '2_2':
        # coding_file = os.path.join(os.path.dirname(__file__), 'Data/gencode.v34.pc_transcripts_test.fa')
        # # coding_file = '/home/wangyunxia/RNACode/tornadoBulid/LncRNAcoder/methods/Data/gencode.v34.pc_transcripts_test.fa'
        # noncoding_file = os.path.join(os.path.dirname(__file__), 'Data/gencode.v34.lncRNA_transcripts_test.fa')
        # Hexamer = Hexamercode.Hexamercoder(textPath, 6, 3, coding_file, noncoding_file).get_hexamer()
        # tran_len = edpfeature.EDPcoder(textPath).get_tran_len()
        # # UTR_len = edpfeature.EDPcoder(textPath).getUTR_len()
        # UTR_cov = edpfeature.EDPcoder(textPath).getUTR_cov()
        # EDP = edpfeature.EDPcoder(textPath).getEDP()
        CTD = CTDcode.CTDcoder(textPath).get_ctd()
        # Kmer = kmer_counts.BasicCounter(textPath, int(3)).get_counts()
        # # EucLogDist = SStructure.makeEucDist(textPath)
        #
        # T1 = pd.concat([Hexamer, tran_len], axis=1, join='inner')
        # # T1 = pd.concat([T1, UTR_len], axis=1, join='inner')
        # T1 = pd.concat([T1, UTR_cov], axis=1, join='inner')
        # T1 = pd.concat([T1, EDP], axis=1, join='inner')
        # T1 = pd.concat([T1, CTD], axis=1, join='inner')
        # T1 = pd.concat([T1, Kmer], axis=1, join='inner')
        # T1 = pd.concat([T1, EucLogDist], axis=1, join='inner')
        return CTD
    elif fun == '2_3':
        # coding_file = os.path.join(os.path.dirname(__file__), 'Data/gencode.v34.pc_transcripts_test.fa')
        # # coding_file = '/home/wangyunxia/RNACode/tornadoBulid/LncRNAcoder/methods/Data/gencode.v34.pc_transcripts_test.fa'
        # noncoding_file = os.path.join(os.path.dirname(__file__), 'Data/gencode.v34.lncRNA_transcripts_test.fa')
        # Hexamer = Hexamercode.Hexamercoder(textPath, 6, 3, coding_file, noncoding_file).get_hexamer()
        # tran_len = edpfeature.EDPcoder(textPath).get_tran_len()
        # # UTR_len = edpfeature.EDPcoder(textPath).getUTR_len()
        # UTR_cov = edpfeature.EDPcoder(textPath).getUTR_cov()
        # EDP = edpfeature.EDPcoder(textPath).getEDP()
        # CTD = CTDcode.CTDcoder(textPath).get_ctd()
        # Kmer = kmer_counts.BasicCounter(textPath, int(3)).get_counts()

        Kmer1 = kmer_counts.BasicCounter(textPath, int(1)).get_counts()
        Kmer2 = kmer_counts.BasicCounter(textPath, int(2)).get_counts()
        Kmer3 = kmer_counts.BasicCounter(textPath, int(3)).get_counts()
        # Kmer4 = kmer_counts.BasicCounter(textPath, int(4)).get_counts()

        T1 = pd.concat([Kmer1, Kmer2], axis=1, join='inner')
        T1 = pd.concat([T1, Kmer3], axis=1, join='inner')
        # T1 = pd.concat([T1, Kmer4], axis=1, join='inner')  
        # EucLogDist = SStructure.makeEucDist(textPath)

        # T1 = pd.concat([Hexamer, tran_len], axis=1, join='inner')
        # # T1 = pd.concat([T1, UTR_len], axis=1, join='inner')
        # T1 = pd.concat([T1, UTR_cov], axis=1, join='inner')
        # T1 = pd.concat([T1, EDP], axis=1, join='inner')
        # T1 = pd.concat([T1, CTD], axis=1, join='inner')
        # T1 = pd.concat([T1, Kmer], axis=1, join='inner')
        # T1 = pd.concat([T1, EucLogDist], axis=1, join='inner')
        return T1
    elif fun == '3':
        Fickett = Fickettcode.Fickettcoder(textPath).get_fickett()
        StopCod = StopCodon.get_stop(textPath)
        return pd.concat([Fickett, StopCod], axis=1, join='inner')
    elif fun == '6':
        return proparcoder.ProtPar(textPath).get_protper()
    elif fun == '7':
        return GCcounts.GCconder(textPath).get_gc()
    elif fun == '13':
        return onehot.Onehot(textPath).get_onehot(1000) #2D 方法
    # elif fun == '14':
    #     return word2vec.process_fasta(textPath, '/home/wangyunxia/RNACode/tornadoBulid/LncRNAcoder/methods/_1101_word2vec_4gram_100_20.model', 2000)#2D 方法
    elif fun == '15':
        return SparseEncod.get_encoding(textPath, 1000)
    elif fun == '16':
        # path2 = '/home/wangyunxia/RNACode/tornadoBulid/LncRNAcoder/methods/RNAFoldOut'
        # RNAFold.setDir(path2)
        # RNAFold.rnaflod(textPath, path2)
        # RNAFold.process(path2, path2)
        # RNAFold.second_file(path2)
        # seqname, rnaresult = RNAFold.encode_feature(path2, 1000)
        seqname1, rnaresult1 = onehot.Onehot(textPath).get_onehot(1000)
        seqname2, rnaresult2 = SparseEncod.get_encoding(textPath, 1000)
        rnaresult3 = np.concatenate((rnaresult1, rnaresult2), axis=2)
        # print("rnaresult3")
        # print(type(rnaresult3))
        # print(seqname1)
        return seqname1, rnaresult3
    elif fun == '17':
        # print('This step is wrong')
        dac = RNA_ac.rna_dac(textPath,2)
        # print('dac is wrong')        
        dcc = RNA_ac.rna_dcc(textPath,2)
        # print('dcc is wrong')
        psenac = RNA_psenac.rna_pc_psednc(textPath)
        # print('psenac is wrong')
        SCPseDNC = RNA_psenac.rna_SCPseDNC(textPath)


        T1 = pd.concat([dac, dcc], axis=1, join='inner')
        T1 = pd.concat([T1, psenac], axis=1, join='inner')
        T1 = pd.concat([T1, SCPseDNC], axis=1, join='inner')

        return T1
    elif fun == '18':
        #正常算法------
        base_pairs = SStructure.cal_base_pairs(textPath)
        # print(base_pairs)

        # base_pairs.to_csv('/public/home/wangyx/LncRNA/smallRNA/code/base_pairs.csv')
        # index_nms = base_pairs.index
        # SStruc = SStructure.extract_SSfeatures(textPath)
        # SStruc.index = index_nms        
        # T1 = pd.concat([SStruc, base_pairs], axis=1, join='inner')
        # print(SStruc.head())

        # SStruc.to_csv('/public/home/wangyx/LncRNA/smallRNA/code/SStruc.csv')

        return base_pairs
        # 正常算法------end
        # SStruc.to_csv(os.path.join(os.path.dirname(__file__), 'Data/Secondary_structureResult.csv'), index=True, header=True)
        #替代算法--------
        # SStruc = pd.read_csv(os.path.join(os.path.dirname(__file__), 'Data/Secondary_structureResult.csv'), sep=',', encoding='utf-8-sig')
        # seq_nameall = []
        # for seq in SeqIO.parse(textPath, 'fasta'):
        #     seq_name = seq.id
        #     seq_nameall.append(seq_name)
        # length = len(seq_nameall)
        # SStruc01 = SStruc.iloc[:length,1:]
        # SStruc01.index = seq_nameall
        # return SStruc01
        # 替代算法--------end
    elif fun == '18_1':
        SStruc = SStructure.makeEIIP(textPath)
        return SStruc
    elif fun == '19_1':
        ACGT_encode = ['1000']
        Solubility = CTDcoderstph.CTDcoder(textPath,ACGT_encode).get_ctd()
        
        
        ACGT_encode_3 = ['1020']
        Solubility_3 = CTDcoderstph.CTDcoder_3class(textPath,ACGT_encode_3).get_ctd()
        T1 = pd.concat([Solubility, Solubility_3], axis=1, join='inner')
        # Kmers
        for i in range(1,4,1):
            Kmer_Solubility = kmer_counts.BasicCounter_2(textPath, int(i), encode = ACGT_encode).get_counts() 
            Kmer_Solubility_3 = kmer_counts.BasicCounter_3(textPath, int(i), encode = ACGT_encode_3).get_counts() 
            
            T1 = pd.concat([T1, Kmer_Solubility], axis=1, join='inner')
            T1 = pd.concat([T1, Kmer_Solubility_3], axis=1, join='inner')
        return T1
    elif fun == '19_101':
        ACGT_encode = ['0100']
        Partition = CTDcoderstph.CTDcoder(textPath,ACGT_encode).get_ctd()
        
        ACGT_encode_3 = ['0102']
        Partition_3 = CTDcoderstph.CTDcoder_3class(textPath,ACGT_encode_3).get_ctd()
        T1 = pd.concat([Partition, Partition_3], axis=1, join='inner')
        # Kmers
        for i in range(1,4,1):
            Kmer_Partition = kmer_counts.BasicCounter_2(textPath, int(i), encode = ACGT_encode).get_counts() 
            Kmer_Partition_3 = kmer_counts.BasicCounter_3(textPath, int(i), encode = ACGT_encode_3).get_counts() 
       
            T1 = pd.concat([T1, Kmer_Partition], axis=1, join='inner')
            T1 = pd.concat([T1, Kmer_Partition_3], axis=1, join='inner')
        
        return T1
    elif fun == '19_2':
        ACGT_encode = ['0101']
        T1 = CTDcoderstph.CTDcoder(textPath,ACGT_encode).get_ctd()

        # Kmers
        for i in range(1,4,1):
            Kmer_Polarizability = kmer_counts.BasicCounter_2(textPath, int(i), encode = ACGT_encode).get_counts() 
            T1 = pd.concat([T1, Kmer_Polarizability], axis=1, join='inner')
      
        return T1
    elif fun == '19_3':
        ACGT_encode = ['0010', '0110']
        Hydrogen = CTDcoderstph.CTDcoder(textPath,ACGT_encode).get_ctd()

        ACGT_encode_3 = ['1200', '0120']
        Hydrogen_3 = CTDcoderstph.CTDcoder_3class(textPath,ACGT_encode_3).get_ctd()
        T1 = pd.concat([Hydrogen, Hydrogen_3], axis=1, join='inner')

        ACGT_encode_00 = ['0010']
        ACGT_encode_01 = ['0110']
        # Kmers
        for i in range(1,4,1):
            Kmer_Hydrogen = kmer_counts.BasicCounter_2(textPath, int(i), encode = ACGT_encode_00).get_counts() 
            Kmer_Hydrogen_01 = kmer_counts.BasicCounter_2(textPath, int(i), encode = ACGT_encode_01).get_counts() 

            Kmer_Hydrogen_3 = kmer_counts.BasicCounter_3(textPath, int(i), encode = ACGT_encode_3).get_counts() 
       
            T1 = pd.concat([T1, Kmer_Hydrogen], axis=1, join='inner')
            T1 = pd.concat([T1, Kmer_Hydrogen_01], axis=1, join='inner')
            T1 = pd.concat([T1, Kmer_Hydrogen_3], axis=1, join='inner')
        
        return T1
    elif fun == '19_4':
        ACGT_encode = ['0001']
        Topological = CTDcoderstph.CTDcoder(textPath,ACGT_encode).get_ctd()

        ACGT_encode_3 = ['1002', '0021']
        Topological_3 = CTDcoderstph.CTDcoder_3class(textPath,ACGT_encode_3).get_ctd()
        T1 = pd.concat([Topological, Topological_3], axis=1, join='inner')
        # Kmers
        for i in range(1,4,1):
            Kmer_Topological = kmer_counts.BasicCounter_2(textPath, int(i), encode = ACGT_encode).get_counts() 
            Kmer_Topological_3 = kmer_counts.BasicCounter_3(textPath, int(i), encode = ACGT_encode_3).get_counts() 
           
            T1 = pd.concat([T1, Kmer_Topological], axis=1, join='inner')
            T1 = pd.concat([T1, Kmer_Topological_3], axis=1, join='inner')
        
        return T1
    elif fun == '19_5':
        ACGT_encode = ['0011']
        T1 = CTDcoderstph.CTDcoder(textPath,ACGT_encode).get_ctd()

        # Kmers
        for i in range(1,4,1):
            Kmer_Molecular = kmer_counts.BasicCounter_2(textPath, int(i), encode = ACGT_encode).get_counts() 

            T1 = pd.concat([T1, Kmer_Molecular], axis=1, join='inner')
        return T1
    else:
        return None

