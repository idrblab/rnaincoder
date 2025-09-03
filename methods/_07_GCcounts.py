import numpy as np
from pandas import DataFrame
from Bio import Seq
import Bio.SeqIO as SeqO
import pandas as pd


class GCconder():
    def __init__(self, infasta=None):
        self.infasta = infasta




    def GetGC_Content(self,seq):
        '''calculate GC content of sequence'''
        A, C, G, T = 1e-9, 1e-9, 1e-9, 1e-9
        for i in range(0, len(seq)):
            if seq[i] == 'A':
                A += 1
            elif seq[i] == 'C':
                C += 1
            elif seq[i] == 'G':
                G += 1
            elif seq[i] == 'T':
                T += 1

        GC = (G + C) / (A + C + G + T)

        return GC


    def get_GC(self):

        with open(self.infasta) as infasta:
            data = [l.strip() for l in infasta]


        new_data = []
        seq = ''
        for i, line in enumerate(data[0:]):
            if line[0] == '>':
                if seq:
                    new_data.append(seq.upper())
                    seq = ''
                else:
                    assert i == 0, 'There may be a header without a sequence at line {}.'.format(i)
                new_data.append(line)
            else:
                seq += line
        new_data.append(seq.upper())
        data = new_data

        clean_data = data
        seqs = clean_data[1::2]
        headers = clean_data[::2]

        counts = np.zeros([len(seqs), 1], dtype=np.float32)

        for i, seq in enumerate(seqs):
            counts[i] = self.GetGC_Content(seq)
        names = headers
        # colname = ['GC content']
        colname = (["GCCoW: GC content of whole sequence"])

        seqname = []
        for seq in SeqO.parse(self.infasta,'fasta'):
            seqid = seq.id
            seqname.append(seqid)

        df = DataFrame(data=counts, index=seqname, columns= colname,dtype='double')

        return df


    def GC1(self,mRNA):
        if len(mRNA) < 3:
            numGC = 0
            mRNA = 'ATG'
        else:
            numGC = mRNA[0::3].count("C") + mRNA[0::3].count("G")
        return numGC * 1.0 / len(mRNA) * 3

    def GC2(self,mRNA):
        if len(mRNA) < 3:
            numGC = 0
            mRNA = 'ATG'
        else:
            numGC = mRNA[1::3].count("C") + mRNA[1::3].count("G")
        return numGC * 1.0 / len(mRNA) * 3

    def GC3(self,mRNA):
        if len(mRNA) < 3:
            numGC = 0
            mRNA = 'ATG'
        else:
            numGC = mRNA[2::3].count("C") + mRNA[2::3].count("G")
        return numGC * 1.0 / len(mRNA) * 3

    def pair_base(self,mRNA):
        if len(mRNA) < 3:
            numpd = 0
        else:
            numpd = mRNA.count("GC") + mRNA.count("CG") + mRNA.count("AT") + mRNA.count("TA")
        return numpd

    def gc1_frame_score(self,seq):
        GC1_in_frame1 = self.GC1(seq)
        GC1_in_frame2 = self.GC1(seq[1:])
        GC1_in_frame3 = self.GC1(seq[2:])
        GC1_all = [GC1_in_frame1, GC1_in_frame2, GC1_in_frame3]
        GC1_frame = ((GC1_all[0] - GC1_all[1]) ** 2 + (GC1_all[0] - GC1_all[2]) ** 2 + (GC1_all[1] - GC1_all[2]) ** 2) / 2
        return GC1_frame

    def gc2_frame_score(self,seq):
        GC2_in_frame1 = self.GC2(seq)
        GC2_in_frame2 = self.GC2(seq[1:])
        GC2_in_frame3 = self.GC2(seq[2:])
        GC2_all = [GC2_in_frame1, GC2_in_frame2, GC2_in_frame3]
        GC2_frame = ((GC2_all[0] - GC2_all[1]) ** 2 + (GC2_all[0] - GC2_all[2]) ** 2 + (GC2_all[1] - GC2_all[2]) ** 2) / 2
        return GC2_frame

    def gc3_frame_score(self,seq):
        GC3_in_frame1 = self.GC3(seq)
        GC3_in_frame2 = self.GC3(seq[1:])
        GC3_in_frame3 = self.GC3(seq[2:])
        GC3_all = [GC3_in_frame1, GC3_in_frame2, GC3_in_frame3]
        GC3_frame = ((GC3_all[0] - GC3_all[1]) ** 2 + (GC3_all[0] - GC3_all[2]) ** 2 + (GC3_all[1] - GC3_all[2]) ** 2) / 2
        return GC3_frame

    def get_GC1(self):

        feaname = (["GCCo1: GC content of 1st position of codons"])
        seqname = []
        GC1_all = []

        for seq in SeqO.parse(self.infasta,'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            GC1_ = self.GC1(seq.seq)
            GC1_all.append(GC1_)

        GC1_all = np.array(GC1_all)
        df = DataFrame(data=GC1_all, index=seqname, columns=feaname,dtype='double')
        return df

    def get_GC2(self):

        feaname = (["GCCo2: GC content of 2nd position of codons"])
        seqname = []
        GC2_all = []

        for seq in SeqO.parse(self.infasta,'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            GC2_ = self.GC2(seq.seq)
            GC2_all.append(GC2_)

        GC2_all = np.array(GC2_all)
        df = DataFrame(data=GC2_all, index=seqname, columns=feaname,dtype='double')
        return df

    def get_GC3(self):

        feaname = (["GCCo3: GC content of 3rd position of codons"])
        seqname = []
        GC3_all = []

        for seq in SeqO.parse(self.infasta,'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            GC3_ = self.GC3(seq.seq)
            GC3_all.append(GC3_)

        GC3_all = np.array(GC3_all)
        df = DataFrame(data=GC3_all, index=seqname, columns=feaname,dtype='double')
        return df

    def get_gc1_frame_score(self):

        # feaname = (["gc1_frame_score"])
        feaname = (["GCC1V: GCCo1 variance frame score"])
        seqname = []
        gc1_frame_score_all = []

        for seq in SeqO.parse(self.infasta,'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            gc1_frame_score_ = self.gc1_frame_score(seq.seq)
            gc1_frame_score_all.append(gc1_frame_score_)

        gc1_frame_score_all = np.array(gc1_frame_score_all)
        df = DataFrame(data=gc1_frame_score_all, index=seqname, columns=feaname,dtype='double')
        return df

    def get_gc2_frame_score(self):

        # feaname = (["gc2_frame_score"])
        feaname = (["GCC2V: GCCo2 variance frame score"])
        seqname = []
        gc2_frame_score_all = []

        for seq in SeqO.parse(self.infasta,'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            gc2_frame_score_ = self.gc2_frame_score(seq.seq)
            gc2_frame_score_all.append(gc2_frame_score_)

        gc2_frame_score_all = np.array(gc2_frame_score_all)
        df = DataFrame(data=gc2_frame_score_all, index=seqname, columns=feaname,dtype='double')
        return df

    def get_gc3_frame_score(self):

        # feaname = (["gc3_frame_score"])
        feaname = (["GCC3V: GCCo3 variance frame score"])
        seqname = []
        gc3_frame_score_all = []

        for seq in SeqO.parse(self.infasta,'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            gc3_frame_score_ = self.gc3_frame_score(seq.seq)
            gc3_frame_score_all.append(gc3_frame_score_)

        gc3_frame_score_all = np.array(gc3_frame_score_all)
        df = DataFrame(data=gc3_frame_score_all, index=seqname, columns=feaname,dtype='double')
        return df

    def get_pair_base(self):

        feaname = (["PB: pair base"])
        seqname = []
        pbcounts_all = []

        for seq in SeqO.parse(self.infasta,'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            pbcounts_ = self.pair_base(seq.seq)
            pbcounts_all.append(pbcounts_)

        pbcount_all = np.array(pbcounts_all)
        df = DataFrame(data=pbcount_all, index=seqname, columns=feaname,dtype='double')
        return df

    def get_gc(self):
        gc = self.get_GC()
        gc1 = self.get_GC1()
        gc2 = self.get_GC2()
        gc3 = self.get_GC3()

        gc1_f = self.get_gc1_frame_score()
        gc2_f = self.get_gc2_frame_score()
        gc3_f = self.get_gc3_frame_score()
        # gpb_f = self.get_pair_base()

        GC_a = pd.concat([gc, gc1,gc2,gc3,gc1_f,gc2_f,gc3_f], axis=1, join='inner')

        return GC_a







