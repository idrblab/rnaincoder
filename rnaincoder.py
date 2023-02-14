#!/usr/bin/env Python
# coding=utf-8
from utility import RNAcoding_01,model
import argparse as agp


def corain_process(type,Afastapath,Interfilepath,Resultpath,dimension,savetype,n_select,com_num,modelnm,Bfastapath=None):
    # calcuation ===========================================
    resultpath = Resultpath + '/' + type
    # encoding part
    if type =='RNAonly':
        RNAcoding_01.RNA_coding(Afastapath, Interfilepath, resultpath, dimension = dimension, savetype = savetype,n_select = n_select)

    elif type =='RNA-RNA':
        RNAcoding_01.RNA_RNA_coding(Afastapath,Bfastapath,Interfilepath,resultpath,dimension = dimension, savetype = savetype,n_select =n_select)

    elif type =='RNA-pro':
        RNAcoding_01.RNA_protein_coding(Afastapath,Bfastapath,Interfilepath,resultpath,dimension = dimension, savetype = savetype,n_select =n_select)

    elif type =='RNA-compound':
        RNAcoding_01.RNA_compound_coding(Afastapath,Bfastapath,Interfilepath,resultpath,dimension = dimension, savetype = savetype,n_select =n_select)
    # evaluation part
    datapath = Resultpath + '/' + type + '/' +  'encoding_features'
    if type =='RNAonly':
        model.evaluation_method(datapath,Interfilepath,resultpath,type = type,com_num = com_num,modelnm = modelnm)
    else:
        model.evaluation_interaction(datapath,Interfilepath,resultpath,type = type,com_num = com_num,modelnm = modelnm)

def main():
    # parameters ===========================================
    # encoding part
    Resultpath = './output'
    type = 'RNAonly'  # RNAonly, RNA-RNA, RNA-pro, RNA-compound
    dimension = '2'  # only type is RNAonly, RNA-RNA, RNA-pro, it is 2 
    savetype = 'npycsv'
    n_select = None

    # evaluation part
    modelnm = 'DNN' # 'RF','svm','xgboost','DNN','CNN'
    com_num = 2  # number of combination methods
    # data part
    Afastapath = './demo/RNA-only/Homo38_small.fasta'
    Bfastapath = None
    Interfilepath = './demo/RNA-only/Homo38_small.csv'

    # Afastapath = './demo/RNA-RNA/SampleData-lncRNA-A.fasta'
    # Bfastapath = './demo/RNA-RNA/SampleData-miRNA-B.fasta'
    # Interfilepath = './demo/RNA-RNA/RNA-RNA-Interacting.csv'

    # Afastapath = './demo/RNA-Protein/SampleData-RNA-A.fasta'
    # Bfastapath = './demo/RNA-Protein/SampleData-Protein-B.fasta'
    # Interfilepath = './demo/RNA-Protein/RNA-Protein-Interacting.csv'

    # Afastapath = './demo/RNA-SamllMoleculer/SampleData-RNA-A.fasta'
    # Bfastapath = './demo/RNA-SamllMoleculer/Small-Moleculer.smi'
    # Interfilepath = './demo/RNA-SamllMoleculer/RNA-small-molecule-Interacting.csv'

    # corain_process(type,Afastapath,Interfilepath,Resultpath,dimension,savetype,n_select,com_num,modelnm,Bfastapath = Bfastapath)
    # 外部输入参数计算----------------------------
    parser = agp.ArgumentParser()
    parser.add_argument('-t','--type',help="The encoding file type: RNAonly, RNA-RNA, RNA-pro, RNA-compound. It is the coding task of a single RNA fasta file, or RNA interaction task.")
    parser.add_argument('-a','--Afastapath',help="The RNA fasta A file path. RNA fasta format is necessary and standard fasta format can be seen https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp.")
    parser.add_argument('-b','--Bfastapath',default = None,help="The encoding fasta B file. If the type is RNAonly, this file is not necessary.")
    parser.add_argument('-l','--Interfilepath',help="The label file. This file format is .csv format.")
    parser.add_argument('-o','--Resultpath',help="The encoding result file path.")
    parser.add_argument('-d','--dimension',type = str,default = '1', help="The feature dimension of encoding feature: 1 or 2.")
    parser.add_argument('-s','--savetype', type = str,default = 'csvnpy',help="The encoding result file type: csvnpy or csv or npy, defatult is csv and npy.")
    parser.add_argument('-n','--n_select', type = int,default = None,help="The encoding feature number: 1-10. Default is all encoding features.")
    parser.add_argument('-c','--com_num',  type = int,default = None,help="Number of combination features: 1-10. Default is all encoding features combinations (1023).")
    parser.add_argument('-m','--modelnm', type = str,default = 'RF',help="Classification model of evaluation feature combination:'RF','svm','xgboost','DNN','CNN'.")
    args = parser.parse_args()
    corain_process(args.type,args.Afastapath,args.Interfilepath,args.Resultpath,args.dimension,args.savetype,args.n_select,args.com_num,args.modelnm,args.Bfastapath)

if __name__ == '__main__':
    main()
