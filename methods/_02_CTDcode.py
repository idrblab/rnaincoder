import sys
import Bio.SeqIO as Seq
import numpy as np
from pandas import DataFrame

class CTDcoder:
    """Generates CTD counts for a fasta file
    """
    def __init__(self, infasta):
        self.infasta = infasta



    def CTD(self,seq):
        # n=len(seq)-1
        n=len(seq)
        n=float(n)
        num_A,num_T,num_G,num_C=0,0,0,0
        AT_trans,AG_trans,AC_trans,TG_trans,TC_trans,GC_trans=0,0,0,0,0,0
        for i in range(len(seq)):
            if seq[i]=="A":
                num_A=num_A+1
            if seq[i]=="T":
                num_T=num_T+1
            if seq[i]=="G":
                num_G=num_G+1
            if seq[i]=="C":
                num_C=num_C+1
        for i in range(len(seq)-1):
            if (seq[i]=="A" and seq[i+1]=="T") or (seq[i]=="T" and seq[i+1]=="A"):
                AT_trans=AT_trans+1
            if (seq[i]=="A" and seq[i+1]=="G") or (seq[i]=="G" and seq[i+1]=="A"):
                AG_trans=AG_trans+1
            if (seq[i]=="A" and seq[i+1]=="C") or (seq[i]=="C" and seq[i+1]=="A"):
                AC_trans=AC_trans+1
            if (seq[i]=="T" and seq[i+1]=="G") or (seq[i]=="G" and seq[i+1]=="T"):
                TG_trans=TG_trans+1
            if (seq[i]=="T" and seq[i+1]=="C") or (seq[i]=="C" and seq[i+1]=="T"):
                TC_trans=TC_trans+1
            if (seq[i]=="G" and seq[i+1]=="C") or (seq[i]=="C" and seq[i+1]=="G"):
                GC_trans=GC_trans+1

        a,t,g,c=0,0,0,0
        A0_dis,A1_dis,A2_dis,A3_dis,A4_dis=0.0,0.0,0.0,0.0,0.0
        T0_dis,T1_dis,T2_dis,T3_dis,T4_dis=0.0,0.0,0.0,0.0,0.0
        G0_dis,G1_dis,G2_dis,G3_dis,G4_dis=0.0,0.0,0.0,0.0,0.0
        C0_dis,C1_dis,C2_dis,C3_dis,C4_dis=0.0,0.0,0.0,0.0,0.0
        for i in range(len(seq)):
            if seq[i]=="A":
                a=a+1
                if a == 1:
                    A0_dis=((i*1.0)+1)/n
                if a == int(round(num_A/4.0)):
                    A1_dis=((i*1.0)+1)/n
                if a == int(round(num_A/2.0)):
                    A2_dis=((i*1.0)+1)/n
                if a == int(round((num_A*3/4.0))):
                    A3_dis=((i*1.0)+1)/n
                if a == num_A:
                    A4_dis=((i*1.0)+1)/n
            if seq[i]=="T":
                t=t+1
                if t == 1:
                    T0_dis=((i*1.0)+1)/n
                if t == int(round(num_T/4.0)):
                    T1_dis=((i*1.0)+1)/n
                if t == int(round((num_T/2.0))):
                    T2_dis=((i*1.0)+1)/n
                if t == int(round((num_T*3/4.0))):
                    T3_dis=((i*1.0)+1)/n
                if t == num_T:
                    T4_dis=((i*1.0)+1)/n
            if seq[i]=="G":
                g=g+1
                if g == 1:
                    G0_dis=((i*1.0)+1)/n
                if g == int(round(num_G/4.0)):
                    G1_dis=((i*1.0)+1)/n
                if g == int(round(num_G/2.0)):
                    G2_dis=((i*1.0)+1)/n
                if g == int(round(num_G*3/4.0)):
                    G3_dis=((i*1.0)+1)/n
                if g == num_G:
                    G4_dis=((i*1.0)+1)/n
            if seq[i]=="C":
                c=c+1
                if c == 1:
                    C0_dis=((i*1.0)+1)/n
                if c == int(round(num_C/4.0)):
                    C1_dis=((i*1.0)+1)/n
                if c == int(round(num_C/2.0)):
                    C2_dis=((i*1.0)+1)/n
                if c == int(round(num_C*3/4.0)):
                    C3_dis=((i*1.0)+1)/n
                if c == num_C:
                    C4_dis=((i*1.0)+1)/n
        return(str(num_A/n),str(num_T/n),str(num_G/n),str(num_C/n),str(AT_trans/(n-1)),str(AG_trans/(n-1)),str(AC_trans/(n-1)),str(TG_trans/(n-1)),str(TC_trans/(n-1)),str(GC_trans/(n-1)),str(A0_dis),str(A1_dis),str(A2_dis),str(A3_dis),str(A4_dis),str(T0_dis),str(T1_dis),str(T2_dis),str(T3_dis),str(T4_dis),str(G0_dis),str(G1_dis),str(G2_dis),str(G3_dis),str(G4_dis),str(C0_dis),str(C1_dis),str(C2_dis),str(C3_dis),str(C4_dis))
    def get_ctd(self):
        # feaname = (["A", "T", "G", "C", "AT", "AG", "AC", "TG", "TC", "GC", "A0", "A1", "A2", "A3", "A4", "T0", "T1", "T2", "T3", "T4", "G0", "G1", "G2", "G3", "G4", "C0", "C1", "C2", "C3", "C4"])
        # feaname = (["CA", "CT", "CG", "CC", "AT", "AG", "AC", "TG", "TC", "GC", "A0", "A1", "A2", "A3", "A4", "T0", "T1", "T2", "T3", "T4", "G0", "G1", "G2", "G3", "G4", "C0", "C1", "C2", "C3", "C4"])
        feaname = (
        ['ComRA: Composition of A','ComRT: Composition of T','ComRG: Composition of G','ComRC: Composition of C','TraAT: Transition between A and T','TraAG: Transition between A and G','TraAC: Transition between A and C','TraTG: Transition between T and G','TraTC: Transition between T and C','TraGC: Transition between G and C','DPRA0: Distribution of 0.00A','DPRA1: Distribution of 0.25A','DPRA2: Distribution of 0.50A','DPRA3: Distribution of 0.75A','DPRA4: Distribution of 1.00A','DPRT0: Distribution of 0.00T','DPRT1: Distribution of 0.25T','DPRT2: Distribution of 0.50T','DPRT3: Distribution of 0.75T','DPRT4: Distribution of 1.00T','DPRG0: Distribution of 0.00G','DPRG1: Distribution of 0.25G','DPRG2: Distribution of 0.50G','DPRG3: Distribution of 0.75G','DPRG4: Distribution of 1.00G','DPRC0: Distribution of 0.00C','DPRC1: Distribution of 0.25C','DPRC2: Distribution of 0.50C','DPRC3: Distribution of 0.75C','DPRC4: Distribution of 1.00C'])

        seqname = []
        num=0
        fea = None
        for seq in Seq.parse(self.infasta,'fasta'):
            num += 1
            seqid = seq.id
            seqname.append(seqid)
            # sequence = seq.seq
            # sequence = sequence.replace('U', 'T')
            if num == 1:
                A,T,G,C,AT,AG,AC,TG,TC,GC,A0,A1,A2,A3,A4,T0,T1,T2,T3,T4,G0,G1,G2,G3,G4,C0,C1,C2,C3,C4 = self.CTD(seq.seq)
                counts = np.array([A,T,G,C,AT,AG,AC,TG,TC,GC,A0,A1,A2,A3,A4,T0,T1,T2,T3,T4,G0,G1,G2,G3,G4,C0,C1,C2,C3,C4]).reshape(1,30)
            else:
                A,T,G,C,AT,AG,AC,TG,TC,GC,A0,A1,A2,A3,A4,T0,T1,T2,T3,T4,G0,G1,G2,G3,G4,C0,C1,C2,C3,C4 = self.CTD(seq.seq)
                fea = np.array([A,T,G,C,AT,AG,AC,TG,TC,GC,A0,A1,A2,A3,A4,T0,T1,T2,T3,T4,G0,G1,G2,G3,G4,C0,C1,C2,C3,C4]).reshape(1,30)
                counts = np.append(counts,fea,axis = 0) 

        df = DataFrame(data=counts, index=seqname, columns=feaname,dtype='double')
        # df = round(df.iloc[:, :], 6)
        #print(df)
        return df

# textPath = 'RPI1807_rna_seq.fa'
# Partition = CTDcoder(textPath).get_ctd()
# Partition.to_csv('Partition_ctd.csv')
