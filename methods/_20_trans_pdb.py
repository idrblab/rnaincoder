import sys
import numpy as np
from Bio.PDB.PDBParser import PDBParser
from pdb_process.PixelateResidue import *
from pdb_process.ModifyName import *
import pandas as pd

from Bio.PDB import DSSP, PDBParser
import os
import re

def preprocess_input(input_pixels):
    """ preprocess input of CNN """

    input_pixels[:, 1] /= 31.0
    input_pixels[:, 2] += 1.0
    input_pixels[:, 2] /= 2.2
    return input_pixels

def process_pdb(datafile,datapath):
    # f = open(datafile, 'r')
    # pdblist = []
    # for line in f.readlines():
    #     line = line.strip('\n')
    #     pdblist.append(line)
    pdbdata = pd.read_csv(datafile)
    pdblist = pdbdata.iloc[:,0].tolist()
    print(pdblist)

    pdbnames = []
    for index, RNA in enumerate(pdblist):
        
        pdbnames.append(RNA)
        
        p = PDBParser(QUIET=True)
        print(RNA)
        pdbfile = datapath + '/' + str(RNA) + '.pdb'

        s = p.get_structure(RNA, pdbfile)
        """ only consider the first model in PDB file """
        model = s[0]


        residues = list(model.get_residues())
        length =  len(residues)
        
        pixels = np.zeros((length, 3, NBINS, NBINS, NBINS))

        try:
            modify_residue_atom_name(residues)
        except:
            print("Not a canonical residue name")
            pass
  

        try:        
            pixels = pixelate_atoms_in_box(model, pixels)
        except:
            print("Not a canonical residue name")
            pass

        pixels = preprocess_input(pixels)
        # print('RNA {} shape: {}'.format(RNA,pixels.shape))
        
        if pixels.shape[0] <= 100:
            pixels_pad = np.pad(pixels,((0,100-(pixels.shape[0])),(0,0),(0,0),(0,0),(0,0)),'constant',constant_values = (0,0))
        else:
            pixels_pad = pixels[0:100,:,:,:,:]        

        pixels_npy = np.expand_dims(pixels_pad, axis = 0)
        # print('pixels_npy {} shape: {}'.format(RNA,pixels_npy.shape))
        if index == 0:
            encode_npy = pixels_npy
        else:
            encode_npy = np.concatenate((encode_npy,pixels_npy), axis = 0)
        # name = Result_path + '/encode_RNA.npy'
        # print('encode_npy {} shape: {}'.format(RNA,encode_npy.shape))
        # np.save(name, encode_npy)
    # print('encode_npy.shape: {}'.format(encode_npy.shape))
    # encode_npy = np.load()
    np.save('/home/wangyunxia/RNACode/tornadoBulid/LncRNAcoder/handlers/files/uploaddata/encode_npy.npy',encode_npy)
    return pdbnames,encode_npy

# def pdbtoclash(datafile,datapath,result_path):

#     pdbdata = pd.read_csv(datafile)
#     pdblist = pdbdata.iloc[:,0].tolist()
#     print(pdblist)

    
#     for index, RNA in enumerate(pdblist):
#         p = PDBParser(QUIET=True)
#         print(RNA)
#         pdbfile = datapath + '/' + str(RNA) + '.pdb'

#         # get the clash score---------------------------------------------------------------------------------------------------------------- 
#         s = p.get_structure(RNA, pdbfile)
#         """ only consider the first model in PDB file """
#         model = s[0]

#         command = 'java -jar /home/wangyunxia/RNACode/tornadoBulid/LncRNAcoder/methods/process_pdb/rnaqua-master/target/rnaqua.jar -c CLASH-SCORE -s ' + pdbfile + ' -o '+ result_path +"/"+ RNA +'.xml'
#         os.system(command)

#         f = open(result_path +"/"+ RNA +'.xml','r')
#         data = f.read()

#         pattern = re.compile(r'score\>(.*?)\<\/score')  # 创建个正着表达式对象，查找数字
#         result1 = pattern.findall(data)
#         if len(result1)<10:
#             buchong = [0 for i in range(10-len(result1))]
#             result1.extend(buchong)
#             clash_score = result1
#         else:
#             clash_score = result1[0:10]  

#         result_npy = np.array(result2)
#         pixels_npy = np.expand_dims(result_npy, axis = 0)

#         if index == 0:
#             encode_npy = pixels_npy
#         else:
#             encode_npy = np.concatenate((encode_npy,pixels_npy), axis = 0)
        
#         f.close()
#     # print('encode_npy.shape: {}'.format(encode_npy.shape))
#     return encode_npy

def pdbto_clash_ssr(datapath,datafile,result_path):


    pdbdata = pd.read_csv(datafile)
    pdblist = pdbdata.iloc[:,0].tolist()
    print(pdblist)

    
    for index, RNA in enumerate(pdblist):
        p = PDBParser(QUIET=True)
        print(RNA)
        pdbfile = datapath + '/' + str(RNA) + '.pdb'


        # get the clash score---------------------------------------------------------------------------------------------------------------- 
        s = p.get_structure(RNA, pdbfile)
        """ only consider the first model in PDB file """
        model = s[0]

        command = 'java -jar /public/home/wangyx/LncRNA/resubmission_code/Case5_standlone_code/corain_standlone/methods/pdb_process/rnaqua-master/target/rnaqua.jar -c CLASH-SCORE -s ' + pdbfile + ' -o '+ result_path +"/"+ RNA +'.xml'
        os.system(command)

        f = open(result_path +"/"+ RNA +'.xml','r')
        data = f.read()

        pattern = re.compile(r'score\>(.*?)\<\/score')  # 创建个正着表达式对象，查找数字
        result1 = pattern.findall(data)
        if len(result1)<10:
            buchong = [0 for i in range(10-len(result1))]
            result1.extend(buchong)
            clash_score = result1
        else:
            clash_score = result1[0:10]

        f.close()
        clash_score_names = ['CSCO' + str(i)+ ': Clash score ' + str(i) for i in range(10)]
        # get the hexlical parameters ------------------------------------------------------------------------------------------------------

        command = 'find_pair ' + pdbfile + ' | analyze'
        os.system(command)
        outputpath = '/public/home/wangyx/LncRNA/resubmission_code/Case5_standlone_code/corain_standlone/' + str(RNA) + '.out'
        # print(outputpath)


        f = open(outputpath,'r')
        data_lines = f.readlines()

        xi = 0
        for data_line in data_lines:
            
            if xi == 2 and 'ave' in data_line:
                # print(data_line)

                datas = data_line.split('\n')[0]
                helical_para = re.split(r"[ ]+", datas)[2:]
                # print(datas)                
                xi = 3

            if xi == 1:
                # print(data_line)
                # names = data_line.split('\n')[0]
                # names = re.split(r"[ ]+", names)[2:]
                xi = 2
                # print(names)

            if 'Local base-pair helical parameters' in data_line:
                # print(data_line)
                xi = 1
        #  concat ----------------------------------------------------------------------------------
        names = ['XDISP: X-disp', 'YDISP: Y-disp', 'HRISE: h-Rise', 'INCLP: Incl parameter', 'TIPPA: Tip parameter', 'HTIST: h-Twist']
        helical_para.extend(clash_score)
        names.extend(clash_score_names)
        df_data = pd.DataFrame(helical_para, index = names)
        df_data.columns = [RNA]
        
        df_data_T =df_data.T
        # print(df_data_T)


        if index == 0:
            df_data_all = df_data_T
        else:
            df_data_all = pd.concat([df_data_all,df_data_T], axis = 0)
        # print(df_data_all)
        f.close()
    print(df_data_all)
    df_data_all.to_csv('/public/home/wangyx/LncRNA/resubmission_code/Case5_standlone_code/corain_standlone/demo/RNA-pdb/df_data_all.csv')
    print('df_data_all.shape: {}'.format(df_data_all.shape))
    df_data_all = pd.read_csv('/public/home/wangyx/LncRNA/resubmission_code/Case5_standlone_code/corain_standlone/demo/RNA-pdb/df_data_all.csv',index_col = 0) 
    print(df_data_all)
    average = df_data_all.mean(axis = 0).tolist()
    print('average')
    print(average)
    return df_data_all


# input_file = '/public/home/wangyx/LncRNA/resubmission_code/Case5_standlone_code/corain_standlone/demo/RNA-pdb/RNA-Label.csv'
# datapath = '/public/home/wangyx/LncRNA/resubmission_code/Case5_standlone_code/corain_standlone/demo/RNA-pdb'
# result_path = '/public/home/wangyx/LncRNA/resubmission_code/Case5_standlone_code/corain_standlone/out/RNApdb'
# # # # datares = process_pdb(input_file,datapath)      #(17, 100, 3, 32, 32, 32)
# # # # datares = pdbtoclash(input_file,datapath,result_path)   #(17, 10)
# datares = pdbto_clash_ssr(datapath,input_file,result_path)    #(17, 6)
# # print('datares.shape: {}'.format(datares.shape))   