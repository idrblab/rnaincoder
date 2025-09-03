#!/usr/bin/env python
'''the python script is downloaded from https://sourceforge.net/projects/rna-cpat/files/v1.2.2/'''
'''deal with Kmer. DNA sequence should only A, C, G, T. python2.7 or newer'''

import os,sys
import math
import string
from optparse import OptionParser
import warnings
#from cpmodule.FrameKmer import kmer_freq_file
import Bio.SeqIO as Seq
import numpy as np
from pandas import DataFrame

####################################
# !/usr/bin/env python
'''deal with Kmer. DNA sequence should only A, C, G, T. python2.7 or newer'''

# import built-in modules
import os, sys
import numpy
import math
from collections import Counter
import re
import itertools

#from cpmodule import ireader


def word_generator(seq, word_size, step_size, frame=0):
    '''generate DNA word from sequence using word_size and step_size. Frame is 0, 1 or2'''
    for i in range(frame, len(seq), step_size):
        word = seq[i:i + word_size]
        if len(word) == word_size:
            yield word


def seq_generator(fastafile):
    '''DNA sequence only contains A,C,G,T,N. sequence with other characters will be removed'''
    tmpseq = ''
    name = ''
    DNA_pat = re.compile(r'^[ACGTN]+$')
    for line in reader(fastafile):
        line = line.strip().upper()
        line = str(line)
        if line.startswith(('#', ' ', '\n')): continue
        if line.startswith(('>', '@')):
            if tmpseq:
                yield [name, tmpseq]
                tmpseq = ''
            name = line.split()[0][1:]
        elif DNA_pat.match(line):
            tmpseq += line
    yield [name, tmpseq]


def all_possible_kmer(l):
    '''return all possible combinations of A,C,G,T,N. only support A,C,G,T,N. l is length of kmer'''
    for i in itertools.product(['A', 'C', 'G', 'T', 'N'], repeat=l):
        yield ''.join(i)


def kmer_freq_file(fastafile, word_size, step_size=1, frame=0, min_count=0):
    '''Calculate kmer frequency from fasta file'''
    seq_num = 0
    ret_dict = {}
    for n, s in seq_generator(fastafile):
        seq_num += 1
        if seq_num == 1:
            count_table = Counter(word_generator(s, word_size=word_size, step_size=step_size, frame=frame))
        else:
            count_table.update(word_generator(s, word_size=word_size, step_size=step_size, frame=frame))

    # return count_table
    for kmer in all_possible_kmer(word_size):
        if not count_table.__contains__(kmer): count_table[kmer] = 0
        if count_table[kmer] >= min_count:
            # print kmer + '\t' + str(count_table[kmer])
            if 'N' in kmer: continue
            ret_dict[kmer] = count_table[kmer]
    return ret_dict


def kmer_freq_seq(seq, word_size, step_size=1, frame=0, min_count=0):
    '''Calculate kmer frequency from DNA sequence. coding. genome is hexamer table calculated
    from coding region and whole genome (as background control)
    '''
    count_table = Counter(word_generator(seq, word_size=word_size, step_size=step_size, frame=frame))
    for kmer in all_possible_kmer(word_size):
        if not count_table.__contains__(kmer): count_table[kmer] = 0
        if count_table[kmer] >= min_count:
            print
            kmer + '\t' + str(count_table[kmer])


def kmer_ratio(seq, word_size, step_size, coding, noncoding):
    if len(seq) < word_size:
        return 0

    sum_of_log_ratio_0 = 0.0
    sum_of_log_ratio_1 = 0.0
    sum_of_log_ratio_2 = 0.0
    frame0_count = 0.0
    frame1_count = 0.0
    frame2_count = 0.0
    for k in word_generator(seq=seq, word_size=word_size, step_size=step_size, frame=0):
        if (not coding.__contains__(k)) or (not noncoding.__contains__(k)):
            continue
        if coding[k] > 0 and noncoding[k] > 0:
            sum_of_log_ratio_0 += math.log(coding[k] / noncoding[k])
        elif coding[k] > 0 and noncoding[k] == 0:
            sum_of_log_ratio_0 += 1
        elif coding[k] == 0 and noncoding[k] == 0:
            continue
        elif coding[k] == 0 and noncoding[k] > 0:
            sum_of_log_ratio_0 -= 1
        else:
            continue
        frame0_count += 1
    '''	
    for k in word_generator(seq=seq, word_size = word_size, step_size=step_size,frame=1):
        if (not coding.has_key(k)) or (not noncoding.has_key(k)):
            continue
        if coding[k]>0 and noncoding[k] >0:
            sum_of_log_ratio_1  +=  math.log( coding[k] / noncoding[k])
        elif coding[k]>0 and noncoding[k] == 0:
            sum_of_log_ratio_1 += 1
        elif coding[k] == 0 and noncoding[k] == 0:
            continue
        elif coding[k] == 0 and noncoding[k] >0 :
            sum_of_log_ratio_1 -= 1
        else:
            continue
        frame1_count += 1

    for k in word_generator(seq=seq, word_size = word_size, step_size=step_size,frame=2):
        if (not coding.has_key(k)) or (not noncoding.has_key(k)):
            continue
        if coding[k]>0 and noncoding[k] >0:
            sum_of_log_ratio_2  +=  math.log( coding[k] / noncoding[k])
        elif coding[k]>0 and noncoding[k] == 0:
            sum_of_log_ratio_2 += 1
        elif coding[k] == 0 and noncoding[k] == 0:
            continue
        elif coding[k] == 0 and noncoding[k] >0 :
            sum_of_log_ratio_2 -= 1
        else:
            continue
        frame2_count += 1
    return max(sum_of_log_ratio_0/frame0_count, sum_of_log_ratio_1/frame1_count,sum_of_log_ratio_2/frame2_count)	
    '''
    try:
        return sum_of_log_ratio_0 / frame0_count
    except:
        return -1


##############################################

"""
read compressed (.gz .bz) files
"""
# !/usr/bin/env python
# encoding: utf-8

import bz2
import gzip
import urllib


def nopen(f, mode="rb"):
    if not isinstance(f, str):
        return f
    if f.startswith("|"):
        p = Popen(f[1:], stdout=PIPE, stdin=PIPE, shell=True)
        if mode[0] == "r": return p.stdout
        return p
    return {"r": sys.stdin, "w": sys.stdout}[mode[0]] if f == "-" \
        else gzip.open(f, mode) if f.endswith((".gz", ".Z", ".z")) \
        else bz2.BZ2File(f, mode) if f.endswith((".bz", ".bz2", ".bzip2")) \
        else urllib.urlopen(f) if f.startswith(("http://", "https://", "ftp://")) \
        else open(f, mode)


def reader(fname):
    for l in nopen(fname):
        yield l.strip()
#########################################################

# !/usr/bin/env python
'''the python script is downloaded from https://sourceforge.net/projects/rna-cpat/files/v1.2.2/'''
'''deal with Kmer. DNA sequence should only A, C, G, T. python2.7 or newer'''

import os, sys
import math
import string
from optparse import OptionParser
import warnings
#from cpmodule.FrameKmer import kmer_freq_file
import Bio.SeqIO as Seq
import numpy as np
from pandas import DataFrame


class Hexamercoder():
    def __init__(self, infasta=None, word_size=6, step_size=3, coding_file=None, noncoding_file=None):
        self.infasta = infasta
        self.coding_file = coding_file
        self.noncoding_file = noncoding_file
        self.word_size = word_size
        self.step_size = step_size

    def word_generator(self, seq, frame=0):
        '''generate DNA word from sequence using word_size and step_size. Frame is 0, 1 or2'''
        for i in range(frame, len(seq), self.step_size):
            word = seq[i:i + self.word_size]
            if len(word) == self.word_size:
                yield word

    def kmer_ratio(self, seq, coding, noncoding):
        if len(seq) < self.word_size:
            return 0

        sum_of_log_ratio_0 = 0.0
        sum_of_log_ratio_1 = 0.0
        sum_of_log_ratio_2 = 0.0
        frame0_count = 0.0
        frame1_count = 0.0
        frame2_count = 0.0
        for k in self.word_generator(seq=seq, frame=0):
            if (not coding.__contains__(k)) or (not noncoding.__contains__(k)):
                continue
            if coding[k] > 0 and noncoding[k] > 0:
                sum_of_log_ratio_0 += math.log(coding[k] / noncoding[k])
            elif coding[k] > 0 and noncoding[k] == 0:
                sum_of_log_ratio_0 += 1
            elif coding[k] == 0 and noncoding[k] == 0:
                continue
            elif coding[k] == 0 and noncoding[k] > 0:
                sum_of_log_ratio_0 -= 1
            else:
                continue
            frame0_count += 1

        try:
            return sum_of_log_ratio_0 / frame0_count
        except:
            return -1

    def coding_nocoding_potential(self):
        ## Make Hexamer frequency table
        cod = kmer_freq_file(fastafile=self.coding_file, word_size=6, step_size=3, frame=0)
        noncod = kmer_freq_file(fastafile=self.noncoding_file, word_size=6, step_size=1, frame=0)

        cod_sum = 0.000001
        cod_sum += sum(cod.values())
        noncod_sum = 0.000001
        noncod_sum += sum(noncod.values())

        colname = (["hexamer", 'coding', 'noncoding'])

        hexamer_all = []
        cod_all = []
        noncod_all = []
        for kmer in cod:
            if 'N' in kmer:
                continue
            hexamer_all.append(kmer)
            cod_l = str(float(cod[kmer] / cod_sum))
            noncod_l = str(float(noncod[kmer] / noncod_sum))

            cod_all.append(cod_l)
            noncod_all.append(noncod_l)

        cod_all = np.expand_dims(np.array(cod_all), axis=1)
        noncod_all = np.expand_dims(np.array(noncod_all), axis=1)
        hexamer_all = np.expand_dims(np.array(hexamer_all), axis=1)
        all_ = np.concatenate((hexamer_all, cod_all, noncod_all), axis=1)

        df_cod = DataFrame(data=all_, columns=colname)

        coding = {}
        noncoding = {}
        for i in range(np.array(df_cod).shape[0]):
            line = np.array(df_cod)[i]
            coding[line[0]] = float(line[1])
            noncoding[line[0]] = float(line[2])
        return coding, noncoding

    def get_hexamer(self):

        # feaname = (["Hexamer"])
        # feaname = (["HS"])
        feaname = (["TraHS: Hexamer score on transcript"])
        seqname = []
        hexamer_all = []
        self.coding, self.noncoding = self.coding_nocoding_potential()

        for seq in Seq.parse(self.infasta, 'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            hexamer_fe = self.kmer_ratio(seq.seq, self.coding, self.noncoding)
            hexamer_all.append(hexamer_fe)

        hexamer_all = np.array(hexamer_all)
        df = DataFrame(data=hexamer_all, index=seqname, columns=feaname)
        return df