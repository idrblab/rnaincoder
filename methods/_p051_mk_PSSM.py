import sys
import os
import numpy as np
import shutil
def command_pssm(content, output_file,pssm_file):
     os.system("./ncbi-blast-2.9.0+/bin/psiblast -query %s -db ./swissprot/swissprot -num_iterations 3 -out %s -out_ascii_pssm %s &" %(content, output_file, pssm_file))
    
def openFile(path, fileName):
    url = path + fileName
    if os.path.exists(url):
        shutil.rmtree(url)  #删除目录，包括目录下的所有文件
    os.mkdir(url)
    return url
def PSSM(proseq):
    outdir1 = openFile(os.path.join(os.path.dirname(__file__), 'pssm/'), 'PSSM')
    outdir0 = openFile(os.path.join(os.path.dirname(__file__), 'pssm/'),'S')
    outdir2 = os.path.join(os.path.dirname(__file__), 'pssm/new')
    inputfile = open(proseq,'r')
    content = ''
    input_file = ''
    output_file = ''
    pssm_file = ''
    #chain_name = []
    for eachline in inputfile:
        if '>' in eachline:
            if len(content):
                #temp_file = open(outdir + '/fasta/' + chain_name,'w')
                temp_file = open(outdir0 + '/' + chain_name,'w')
                temp_file.write(content)
                #input_file = outdir + '/fasta/' + chain_name
                input_file = outdir0 + '/' + chain_name
                output_file = outdir2 + '/' + chain_name + '.out'
                pssm_file = outdir1 + '/' + chain_name + '.pssm'                
                command_pssm(input_file, output_file, pssm_file)
                temp_file.close
            content = ''
            chain_name = eachline[1:(len(eachline)-1)]
            print(type(chain_name))
        else:
            content += ''.join(eachline)
        #print content
        #print chain_name
    if len(content):
        temp_file = open(outdir0 + '/' +chain_name,'w')
        temp_file.write(content)
        input_file = outdir0 + '/' +chain_name
        output_file = outdir2 + '/' + chain_name + '.out'
        pssm_file = outdir1 + '/' + chain_name + '.pssm'
        command_pssm(input_file, output_file,pssm_file)
        temp_file.close
    inputfile.close()

