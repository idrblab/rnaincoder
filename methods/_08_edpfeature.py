import os,sys
import numpy as np
import Bio.SeqIO as Seq
from pandas import DataFrame

class EDPcoder:
    def __init__(self, infasta=None):
        self.infasta = infasta

        ## AA list 
        self._AA_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

        ## Di-codon list
        self._Di_Codon_list = []
        for aa1 in self._AA_list:
            for aa2 in self._AA_list:
                self._Di_Codon_list.append(aa1+aa2)

        ## 2mer list
        self._DNA = ['A', 'C', 'G', 'T']

        ## 3mer list
        self._3mer_list = []
        for dna1 in self._DNA:
            for dna2 in self._DNA:
                for dna3 in self._DNA:
                    self._3mer_list.append(dna1+dna2+dna3)

        ## 6mer list
        self._6mer_list = []
        for mer1 in self._3mer_list:
            for mer2 in self._3mer_list:
                self._6mer_list.append(mer1+mer2)



    def IUPAC_2mer(self,seq):
        '''Return a list of all possible 2mers of the sequence'''

        ## IUPAC code
        _IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'M': 'AC', 'K': 'GT', 'S': 'CG',
                  'W': 'AT', 'H': 'ACT', \
                  'B': 'CGT', 'V': 'ACG', 'D': 'AGT', 'N': 'ACGT'}

        kmer_list = []
        for dna1 in _IUPAC[seq[0]]:
            for dna2 in _IUPAC[seq[1]]:
                kmer_list.append(dna1 + dna2)
        return kmer_list



    def Codon2AA2(self,codon):
        '''convert codon to aa'''
        if codon == "TTT" or codon == "TTC":
            return 'F'
        elif codon == 'TTA' or codon == 'TTG' or codon == 'CTT' or codon == 'CTA' or codon == 'CTC' or codon == 'CTG':
            return 'L'
        elif codon == 'ATT' or codon == 'ATC' or codon == 'ATA':
            return 'I'
        elif codon == 'ATG':
            return 'M'
        elif codon == 'GTA' or codon == 'GTC' or codon == 'GTG' or codon == 'GTT':
            return 'V'
        elif codon == 'GAT' or codon == 'GAC':
            return 'D'
        elif codon == 'GAA' or codon == 'GAG':
            return 'E'
        elif codon == 'TCA' or codon == 'TCC' or codon == 'TCG' or codon == 'TCT':
            return 'S'
        elif codon == 'CCA' or codon == 'CCC' or codon == 'CCG' or codon == 'CCT':
            return 'P'
        elif codon == 'ACA' or codon == 'ACG' or codon == 'ACT' or codon == 'ACC':
            return 'T'
        elif codon == 'GCA' or codon == 'GCC' or codon == 'GCG' or codon == 'GCT':
            return 'A'
        elif codon == 'TAT' or codon == 'TAC':
            return 'Y'
        elif codon == 'CAT' or codon == 'CAC':
            return 'H'
        elif codon == 'CAA' or codon == 'CAG':
            return 'Q'
        elif codon == 'AAT' or codon == 'AAC':
            return 'N'
        elif codon == 'AAA' or codon == 'AAG':
            return 'K'
        elif codon == 'TGT' or codon == 'TGC':
            return 'C'
        elif codon == 'TGG':
            return 'W'
        elif codon == 'CGA' or codon == 'CGC' or codon == 'CGG' or codon == 'CGT':
            return 'R'
        elif codon == 'AGT' or codon == 'AGC':
            return 'S'
        elif codon == 'AGA' or codon == 'AGG':
            return 'R'
        elif codon == 'GGA' or codon == 'GGC' or codon == 'GGG' or codon == 'GGT':
            return 'G'
        # stop codon
        elif codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
            return 'J'
        else:
            return 'Z'     ## IUPAC Ambiguity Codes

    def IUPAC_3mer(self, seq):
        '''Return a list of all possible 3mers of the sequence'''
        ## IUPAC code
        _IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'M': 'AC', 'K': 'GT', 'S': 'CG',
                  'W': 'AT', 'H': 'ACT', \
                  'B': 'CGT', 'V': 'ACG', 'D': 'AGT', 'N': 'ACGT'}

        kmer_list = []
        for dna1 in _IUPAC[seq[0]]:
            for dna2 in _IUPAC[seq[1]]:
                for dna3 in _IUPAC[seq[2]]:
                    if self.Codon2AA2(dna1 + dna2 + dna3) != "J":
                        kmer_list.append(dna1 + dna2 + dna3)
        return kmer_list

    def GetORF_UTR(self,seq):
        '''Get ORF and UTR from sequence'''
        STP = {}
        STP[0] = []
        STP[1] = []
        STP[2] = []
        STP[0].append(0)
        STP[1].append(1)
        STP[2].append(2)
        AAnum = int(len(seq) / 3)

        for i in range(0, 3):
            for j in range(0, AAnum):
                tmp = seq[(i+3*j):(i+3*(j+1))]
                if tmp == 'TAG' or tmp == 'TAA' or tmp == 'TGA':
                    STP[i].append(i+3*j)

        ORF = {}

        for i in range(0,3):
            for j in range(1, len(STP[i])):
                tmpN = int((STP[i][j] - STP[i][j-1])/3)
                for k in range(0, tmpN):
                    tmpS = seq[ (STP[i][j-1] + 3*k):(STP[i][j-1] + 3*(k+1)) ] 
                    if tmpS == 'ATG':
                        ORF[3*k + STP[i][j-1]] = STP[i][j] + 3
                        break

            codonNum = int((len(seq) - STP[i][-1]) / 3)
            for k in range(codonNum):
                if seq[ (STP[i][-1] + 3*k): (STP[i][-1] + 3*(k+1)) ] == "ATG":
                    ORF[ STP[i][-1] + 3*k ] = len(seq)
                    break

        # longest ORF
        ORFseq = []
        ORFlen = []
        ORFstart = []
        ORFend = []
        for (k,v) in ORF.items():
            ORFseq.append(seq[k:(v)])
            ORFlen.append(v - k)
            ORFstart.append(k)
            ORFend.append(v)

        if ORF:
            idx = np.argmax(ORFlen)
            ORF_l = ORFseq[idx]
            UTR5 = ''
            UTR3 = ''
            if len(seq[ 0:ORFstart[idx] ]) > 0:
                UTR5 = seq[ 0:ORFstart[idx] ]
            if len(seq[ ORFend[idx]: ]) > 0:
                UTR3 = seq[ ORFend[idx]: ]
            return [ ORF_l, UTR5, UTR3, ORFstart[idx], ORFend[idx] ] 
        else:
            return ['', '', '', 0, 0]

    def GetEDP_noORF(self):

        Codon = {}
        for aa in self._AA_list:
            Codon[aa] = 1e-9 

        sum_codon = 1e-9 * 20 

        H = 0.0
        for (k,v) in Codon.items():
            Codon[k] /= sum_codon
            Codon[k] = -Codon[k] * np.log2(Codon[k])
            H += Codon[k]

        EDP = {}
        value = []
        for (k,v) in Codon.items():
            EDP[k] = Codon[k] / H

        outline = ''
        value = []
        for i in range(20):
            outline += str(0) + "\t"
            value.append(str(0))

        return value

    def GetEDP(self, seq, transcript_len):
        '''get features including: ORF length, ORF ratio, ORF EDP of codon'''

        # entropy density
        Codon = {}
        for aa in self._AA_list:
            Codon[aa] = 1e-9 

        sum_codon = 1e-9 * 20 

        if(len(seq) > 3):
            num = int(len(seq) / 3)
            # print('num')
            # print(num)
            for i in range(0,num) :
                if self.Codon2AA2( seq[i*3:(i+1)*3] ) == "J":
                    continue
                ## consider the IUPAC codon
                elif self.Codon2AA2( seq[i*3:(i+1)*3] ) == "Z":
                    tmp_kmer_list = self.IUPAC_3mer(seq[i*3:(i+1)*3])
                    for tmp_kmer in tmp_kmer_list:
                        Codon[ self.Codon2AA2(tmp_kmer) ] += 1.0 / len(tmp_kmer_list)
                    sum_codon += 1.0
                else:
                    Codon[ self.Codon2AA2( seq[i*3:(i+1)*3] ) ] += 1.0
                    sum_codon += 1.0

            H = 0.0
            for (k,v) in Codon.items():
                Codon[k] /= sum_codon
                Codon[k] = -Codon[k] * np.log2(Codon[k])
                H += Codon[k]

            EDP = {}
            for (k,v) in Codon.items():
                EDP[k] = Codon[k] / H

            value = []
            for (k,v) in EDP.items():
                #outline += str(v) + "\t"
                value.append(v)
            
            return value
        
    def GetKmerEDP(self, seq):

        Kmer = {}
        for aa in _Kmer_list:
            Kmer[aa] = 1e-9

        sum_Kmer = 1e-9 * 16 

        if(len(seq) > 3):
            for i in range(0,len(seq)-1) :
                ## consider IUPAC kmer
                if seq[ i:(i+2) ] not in _Kmer_list:
                    tmp_kmer_list = self.IUPAC_2mer(seq[i:(i+2)])
                    for tmp_kmer in tmp_kmer_list:
                        Kmer[ tmp_kmer ] += 1.0 / len(tmp_kmer_list)
                else:
                    Kmer[ seq[ i:(i+2)] ] += 1.0 
                sum_Kmer += 1.0

            H = 0.0
            for (k,v) in Kmer.items():
                Kmer[k] /= sum_Kmer
                Kmer[k] = -Kmer[k] * np.log2(Kmer[k])
                H += Kmer[k]

            EDP = {}
            for (k,v) in Kmer.items():
                EDP[k] = Kmer[k] / H

            outline = ''
            for (k,v) in EDP.items():
                outline += str(v) + "\t"

            return outline

        else:
            return GetKmerEDP_Default()
#############################################
    def getEDP_orf(self):

        seqname = []
        EDP_all = []

        for seq in Seq.parse(self.infasta,'fasta'):
            #print(seq)
            seqid = seq.id
            seqname.append(seqid)
            ORF, UTR5, UTR3, start, end = self.GetORF_UTR(seq.seq)
            #print("this is %s" % ORF)
            transcript_len = len(seq.seq)

            ## EDP feature for longest ORF
            tmp_seq = ORF
            if len(tmp_seq) < 6:
                EDP_fea = self.GetEDP_noORF()
                EDP_all.append(EDP_fea)
            else:
                EDP_fea = self.GetEDP(tmp_seq, transcript_len)
                EDP_all.append(EDP_fea)
                #print(len(EDP_fea))
            #Kmer_EDP_fea = GetKmerEDP(tmp_seq)

        EDP_all = np.array(EDP_all)
        from pandas import DataFrame
        # coname = (['Aorf', 'Corf', 'Dorf', 'Eorf', 'Forf', 'Gorf', 'Horf', 'Iorf', 'Korf', 'Lorf', 'Morf', 'Norf', 'Porf', 'Qorf', 'Rorf', 'Sorf', 'Torf', 'Vorf', 'Worf', 'Yorf'])

        coname = (['ORFEA: Entropy density A on ORF','ORFEC: Entropy density C on ORF','ORFED: Entropy density D on ORF','ORFEE: Entropy density E on ORF','ORFEF: Entropy density F on ORF','ORFEG: Entropy density G on ORF','ORFEH: Entropy density H on ORF','ORFEI: Entropy density I on ORF','ORFEK: Entropy density K on ORF','ORFEL: Entropy density L on ORF','ORFEM: Entropy density M on ORF','ORFEN: Entropy density N on ORF','ORFEP: Entropy density P on ORF','ORFEQ: Entropy density Q on ORF','ORFER: Entropy density R on ORF','ORFES: Entropy density S on ORF','ORFET: Entropy density T on ORF','ORFEV: Entropy density V on ORF','ORFEW: Entropy density W on ORF','ORFEY: Entropy density Y on ORF'])
        feaname = seqname
        df = DataFrame(data=EDP_all, index=feaname, columns=coname)

        return df

    def getEDP(self):

        seqname = []
        EDP_all = []

        for seq in Seq.parse(self.infasta,'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            transcript_len = len(seq.seq)

            ## EDP feature for longest ORF
            tmp_seq = seq.seq
            if len(tmp_seq) < 6:
                EDP_fea = self.GetEDP_noORF()
                EDP_all.append(EDP_fea)
            else:
                EDP_fea = self.GetEDP(tmp_seq, transcript_len)
                EDP_all.append(EDP_fea)
            #Kmer_EDP_fea = GetKmerEDP(tmp_seq)

        EDP_all = np.array(EDP_all)
        from pandas import DataFrame
        # coname = (['PA', 'PC', 'D', 'E', 'F', 'PG', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'PT', 'V', 'W', 'Y'])
        coname = (['TraEA: Entropy density A on transcript','TraEC: Entropy density C on transcript','TraED: Entropy density D on transcript','TraEE: Entropy density E on transcript','TraEF: Entropy density F on transcript','TraEG: Entropy density G on transcript','TraEH: Entropy density H on transcript','TraEI: Entropy density I on transcript','TraEK: Entropy density K on transcript','TraEL: Entropy density L on transcript','TraEM: Entropy density M on transcript','TraEN: Entropy density N on transcript','TraEP: Entropy density P on transcript','TraEQ: Entropy density Q on transcript','TraER: Entropy density R on transcript','TraES: Entropy density S on transcript','TraET: Entropy density T on transcript','TraEV: Entropy density V on transcript','TraEW: Entropy density W on transcript','TraEY: Entropy density Y on transcript'])
        feaname = seqname

        df = DataFrame(data=EDP_all, index=feaname, columns=coname)

        return df

    def get_tran_len(self):

        seqname = []
        EDP_all = []

        for seq in Seq.parse(self.infasta,'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            transcript_len = len(seq.seq)

            EDP_all.append(transcript_len)

        EDP_all = np.array(EDP_all)
        
        coname = (['TrLen: Transcript length'])
        # coname = (['TL'])
        feaname = seqname
        df = DataFrame(data=EDP_all, index=feaname, columns=coname)

        return df

    def getUTR_cov(self):

        seqname = []
        UTR5_all = []
        UTR3_all = []

        for seq in Seq.parse(self.infasta,'fasta'):
            #print(seq)
            seqid = seq.id
            seqname.append(seqid)
            ORF, UTR5, UTR3, start, end = self.GetORF_UTR(seq.seq)
            #print("this is %s" % ORF)
            UTR5_len = len(UTR5)
            UTR3_len = len(UTR3)
            transcript_len = len(seq.seq)

            UTR5_cov = len(UTR5) / transcript_len
            UTR3_cov = len(UTR3) / transcript_len

            UTR5_all.append(UTR5_cov)
            UTR3_all.append(UTR3_cov)

        UTR5_all = np.expand_dims(np.array(UTR5_all),axis = 1)
        UTR3_all = np.expand_dims(np.array(UTR3_all),axis = 1)
        UTR = np.concatenate((UTR5_all,UTR3_all),axis=1)

        coname = (["C5UTR: Coverage of 5 untranslated region", "C3UTR: Coverage of 3 untranslated region"])
        # coname = (['U5C', 'U3C'])
        feaname = seqname
        df = DataFrame(data=UTR, index=feaname, columns=coname)

        return df

    def getUTR_len(self):

        seqname = []
        UTR5_all = []
        UTR3_all = []
        UTR5_covall = []
        UTR3_covall = []

        for seq in Seq.parse(self.infasta,'fasta'):
            #print(seq)
            seqid = seq.id
            seqname.append(seqid)
            ORF, UTR5, UTR3, start, end = self.GetORF_UTR(seq.seq)
            #print("this is %s" % ORF)
            UTR5_len = len(UTR5)
            UTR3_len = len(UTR3)
            transcript_len = len(seq.seq)

            UTR5_cov = len(UTR5) / transcript_len
            UTR3_cov = len(UTR3) / transcript_len

            UTR5_all.append(UTR5_len)
            UTR3_all.append(UTR3_len)
            UTR5_covall.append(UTR5_cov)
            UTR3_covall.append(UTR3_cov)

        UTR5_all = np.expand_dims(np.array(UTR5_all),axis = 1)
        UTR3_all = np.expand_dims(np.array(UTR3_all),axis = 1)
        UTR5_covall = np.expand_dims(np.array(UTR5_covall),axis = 1)
        UTR3_covall = np.expand_dims(np.array(UTR3_covall),axis = 1)
        UTR = np.concatenate((UTR5_all,UTR3_all),axis=1)

        coname = (['L5UTR: Length of 5 untranslated region', 'L3UTR: Length of 3 untranslated region'])
        # coname = (['U5L', 'U3L'])
        feaname = seqname
        df = DataFrame(data=UTR, index=feaname, columns=coname)

        return df