#!/usr/bin/env python
import os
import sys
current_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(current_dir)
sys.path.append("..")

import sys
import Bio.SeqIO as Seq
import numpy as np
from pandas import DataFrame

class ORF_count():

    def __init__(self, infasta=None):
        self.infasta = infasta
        self.start_codons = 'ATG'
        self.stop_codons = 'TAG,TAA,TGA'
        self.Coverage = 0

    def extract_feature_from_seq(self, seq, stt, stp):
        '''extract features of sequence from fasta entry'''

        stt_coden = stt.strip().split(',')
        stp_coden = stp.strip().split(',')
        transtab = str.maketrans("ACGTNX", "TGCANX")
        mRNA_seq = seq.upper()
        mRNA_size = len(seq)
        tmp = ExtractORF(mRNA_seq)
        (CDS_size1, CDS_integrity, CDS_seq1) = tmp.longest_ORF(start=stt_coden, stop=stp_coden)
        return (mRNA_size, CDS_size1, CDS_integrity)

    def len_cov(self, seq):
        # print('seq')
        # print(seq)
        (mRNA_size, CDS_size, CDS_integrity) = self.extract_feature_from_seq(seq=seq, stt=self.start_codons,
                                                                             stp=self.stop_codons)
        mRNA_len = mRNA_size
        CDS_len = CDS_size
        self.Coverage = float(CDS_len) / mRNA_len
        Integrity = CDS_integrity
        return (CDS_len, self.Coverage, Integrity)

    def get_orf(self, seq, stt, stp):
        stt_coden = stt.strip().split(',')
        stp_coden = stp.strip().split(',')
        transtab = str.maketrans("ACGTNX", "TGCANX")
        mRNA_seq = seq.upper()
        tmp = ExtractORF(mRNA_seq)
        (CDS_size1, CDS_integrity, CDS_seq1) = tmp.longest_ORF(start=stt_coden, stop=stp_coden)
        return (CDS_seq1, CDS_integrity)

    def get_orf_frame_score(self,seq):
        ORF_length_in_frame1, _ = self.get_orf(seq, stt=self.start_codons,stp=self.stop_codons)
        ORF_length_in_frame2, _ = self.get_orf(seq[1:], stt=self.start_codons,stp=self.stop_codons)
        ORF_length_in_frame3, _ = self.get_orf(seq[2:], stt=self.start_codons,stp=self.stop_codons)

        ORF_length_in_frame1 = len(ORF_length_in_frame1)
        ORF_length_in_frame2 = len(ORF_length_in_frame2)
        ORF_length_in_frame3 = len(ORF_length_in_frame3)

        ORF_len = [ORF_length_in_frame1, ORF_length_in_frame2, ORF_length_in_frame3]
        ORF_frame = ((ORF_len[0] - ORF_len[1]) ** 2 + (ORF_len[0] - ORF_len[2]) ** 2 + (
                    ORF_len[1] - ORF_len[2]) ** 2) / 2
        return ORF_frame

    def get_ORFcov(self):
        Len_all = []
        Cov_all = []
        inte_fe_all = []
        orf_frame_score_all = []
        seqname = []

        #colna = (["ORF-length", 'ORF-coverage', 'ORF-integrity', 'orf_frame_score'])
        colna = (['ORFCo: ORF coverage'])

        for seq in Seq.parse(self.infasta, 'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            sequence = seq.seq
            sequence = sequence.replace('U', 'T')
            Len, Cov, inte_fe = self.len_cov(sequence)
            Len_all.append(Len)
            Cov_all.append(Cov)
            inte_fe_all.append(inte_fe)
            orf_frame_score = self.get_orf_frame_score(sequence)
            orf_frame_score_all.append(orf_frame_score)

        Len_all = np.expand_dims(np.array(Len_all), axis=1)
        Cov_all = np.expand_dims(np.array(Cov_all), axis=1)
        inte_fe_all = np.expand_dims(np.array(inte_fe_all), axis=1)
        orf_frame_score_all = np.expand_dims(np.array(orf_frame_score_all), axis=1)

        all_ = np.concatenate((Len_all, Cov_all, inte_fe_all, orf_frame_score_all), axis=1)
        df_cod = DataFrame(data=Cov_all, index=seqname, columns=colna)
        return df_cod

    def get_ORFleng(self):
        Len_all = []
        Cov_all = []
        inte_fe_all = []
        orf_frame_score_all = []
        seqname = []

        colna = (['ORF-coverage',"ORF-length"])

        for seq in Seq.parse(self.infasta, 'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            sequence = seq.seq
            sequence = sequence.replace('U', 'T')
            Len, Cov, inte_fe = self.len_cov(sequence)
            Len_all.append(Len)
            Cov_all.append(Cov)
            inte_fe_all.append(inte_fe)
            orf_frame_score = self.get_orf_frame_score(sequence)
            orf_frame_score_all.append(orf_frame_score)

        Len_all = np.expand_dims(np.array(Len_all), axis=1)
        Cov_all = np.expand_dims(np.array(Cov_all), axis=1)
        inte_fe_all = np.expand_dims(np.array(inte_fe_all), axis=1)
        orf_frame_score_all = np.expand_dims(np.array(orf_frame_score_all), axis=1)

        all_ = np.concatenate((Len_all, Cov_all, inte_fe_all, orf_frame_score_all), axis=1)
        df_cod = DataFrame(data=Len_all, index=seqname, columns=colna)
        return df_cod

    def get_ORFinte(self):
        Len_all = []
        Cov_all = []
        inte_fe_all = []
        orf_frame_score_all = []
        seqname = []

        colna = (['ORF-coverage',"ORF-length",'ORF-integrity'])

        for seq in Seq.parse(self.infasta, 'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            sequence = seq.seq
            sequence = sequence.replace('U', 'T')
            Len, Cov, inte_fe = self.len_cov(sequence)
            Len_all.append(Len)
            Cov_all.append(Cov)
            inte_fe_all.append(inte_fe)
            orf_frame_score = self.get_orf_frame_score(sequence)
            orf_frame_score_all.append(orf_frame_score)

        Len_all = np.expand_dims(np.array(Len_all), axis=1)
        Cov_all = np.expand_dims(np.array(Cov_all), axis=1)
        inte_fe_all = np.expand_dims(np.array(inte_fe_all), axis=1)
        orf_frame_score_all = np.expand_dims(np.array(orf_frame_score_all), axis=1)

        all_ = np.concatenate((Len_all, Cov_all, inte_fe_all, orf_frame_score_all), axis=1)
        df_cod = DataFrame(data=inte_fe_all, index=seqname, columns=colna)
        return df_cod

    def get_ORF(self):
        Len_all = []
        Cov_all = []
        inte_fe_all = []
        orf_frame_score_all = []
        seqname = []

        # colna = (['ORF-length',"ORF-coverage",'ORF-integrity','orf_frame_score'])
        # colna = (['OL', "OC", 'OI', 'OR'])
        colna = (['LORFL: Longest ORF length', "ORFCo: ORF coverage", 'ORFIn: ORF integrity', 'ORFFS: ORF frame score'])
        print(self.infasta)
        # inde = 0
        for seq in Seq.parse(self.infasta, 'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            # sequence = seq.seq
            # sequence = sequence.replace('U', 'T')
            Len, Cov, inte_fe = self.len_cov(seq.seq)
            Len_all.append(Len)
            Cov_all.append(Cov)
            inte_fe_all.append(inte_fe)

            orf_frame_score = self.get_orf_frame_score(seq.seq)

            orf_frame_score_all.append(orf_frame_score)
            # print('aaaaaaaaaaa')
        #     inde =inde + 1
        #     print(inde)
        # print('bbbbbbbbbbbbbb')
        print(orf_frame_score_all)    
        Len_all = np.expand_dims(np.array(Len_all), axis=1)
        Cov_all = np.expand_dims(np.array(Cov_all), axis=1)
        inte_fe_all = np.expand_dims(np.array(inte_fe_all), axis=1)
        orf_frame_score_all = np.expand_dims(np.array(orf_frame_score_all), axis=1)
        all_ = np.concatenate((Len_all, Cov_all, inte_fe_all, orf_frame_score_all), axis=1)
        df_cod = DataFrame(data=all_, index=seqname, columns=colna)
        return df_cod


# !/usr/bin/env python

'''
Extract the most probable ORF in a given sequence 
The most probable ORF is the longest open reading frame found in the sequence
When having same length, the upstream ORF is selected
modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/v1.2.1/
'''


class ExtractORF:
    def __init__(self, seq):
        self.seq = seq
        self.result = (0, 0, 0, 0)
        self.longest = 0

    def codons(self, frame):
        start_coord = frame
        while start_coord + 3 <= len(self.seq):
            yield (self.seq[start_coord:start_coord + 3], start_coord)
            start_coord += 3

    def longest_orf_in_seq(self, frame_number, start_codon, stop_codon):
        codon_posi = self.codons(frame_number)
        start_codons = start_codon
        stop_codons = stop_codon
        while True:
            try:
                codon, index = next(codon_posi)
            except StopIteration:
                break
            if codon in start_codons and codon not in stop_codons:
                ORF_start = index
                end = False
                while True:
                    try:
                        codon, index = next(codon_posi)
                    except StopIteration:
                        end = True
                        integrity = -1
                    if codon in stop_codons:
                        integrity = 1
                        end = True
                    if end:
                        ORF_end = index + 3
                        ORF_Length = (ORF_end - ORF_start)
                        if ORF_Length > self.longest:
                            self.longest = ORF_Length
                            self.result = (integrity, ORF_start, ORF_end, ORF_Length)
                        if ORF_Length == self.longest and ORF_start < self.result[1]:
                            self.result = (integrity, ORF_start, ORF_end, ORF_Length)
                        break

    def longest_ORF(self, start=['ATG'], stop=['TAA', 'TAG', 'TGA']):
        orf_seq = ""
        for frame in range(3):
            self.longest_orf_in_seq(frame, start, stop)
        orf_seq = self.seq[self.result[1]:self.result[2]]
        ORF_integrity = self.result[0]
        ORF_length = self.result[3]
        return ORF_length, ORF_integrity, orf_seq
