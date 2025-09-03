import numpy as np 
import os

def formateachline(eachline):
    str_list = eachline.split()
    col = str_list[0]
    for i in range(1,22):
        col += '\t' + str_list[i]
    col += '\n'
    return col
def pssm(eachline):
    str_list = eachline.split()
    
    str_list = str_list[2:22]
    str_list = [int(i) for i in str_list]
    col = str_list
    
    return col
def pad_trunc(pro_array, N):
    protein = np.zeros((N, 20))
    length = len(pro_array)
    if length <= N:
        protein[0:length, :] = pro_array
    else:
        protein[0:N, :] = pro_array[0:N, :]
    return protein
def simplifypssm(pssmdir, N):
    listfile = os.listdir(pssmdir)
    seqnamelist = []
    for x in listfile:
        x01 = x.split('.')[0]
        seqnamelist.append(x01)

    big_array = []
    for eachfile in listfile:
        with open(pssmdir + '/' + eachfile,'r') as inputpssm:
           
            count = 0
            pssm_array = []
            for eachline in inputpssm:
                count +=1
                if count <= 3:
                    continue
                if not len(eachline.strip()):
                    break
                # oneline = formateachline(eachline)
                one_array = pssm(eachline)
                pssm_array.append(one_array)
                   
                    # outfile.write(''.join(oneline))
            pssm_array = np.array(pssm_array)
            pssm_array = pad_trunc(pssm_array, N)
            big_array.append(pssm_array)
    
    big_array = np.array(big_array)
    return seqnamelist,big_array

