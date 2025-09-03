import sys
import re
from Bio.Seq import Seq
#from ORF_exact import ExtractORF
from Bio.SeqUtils import ProtParam
import numpy as np
from pandas import DataFrame
from Bio import SeqIO

class ProtPar():

	def __init__(self, infasta=None):
		self.infasta = infasta


	def mRNA_translate(self,mRNA):
		return Seq(mRNA).translate()

	def protein_param(self,putative_seqprot):
		return (putative_seqprot.instability_index(),putative_seqprot.isoelectric_point(),putative_seqprot.gravy(),putative_seqprot.molecular_weight())

	def param(self,seq):
		strinfoAmbiguous = re.compile("X|B|Z|J|U",re.I)
		ptU = re.compile("U",re.I)
		seqRNA = ptU.sub("T",str(seq).strip())
		seqRNA = seqRNA.upper()
		CDS_size1,CDS_integrity,seqCDS= ExtractORF(seqRNA).longest_ORF(start=['ATG'],stop=['TAA','TAG','TGA'])
		seqprot = self.mRNA_translate(seqCDS)
		if '*' in seqprot:
			nPos = seqprot.index('*')
			if nPos != len(seqprot)-1:
				seqprot = seqprot.split('*')[0]
				seqprot = seqprot+'*'
				print('seqprot')
				print(seqprot)
		pep_len = len(seqprot.strip("*"))
		newseqprot = strinfoAmbiguous.sub("",str(seqprot))
		protparam_obj = ProtParam.ProteinAnalysis(str(newseqprot.strip("*")))
		if pep_len > 0:
			Instability_index,PI,Gravy,Mw = self.protein_param(protparam_obj)
			pI_Mw = np.log10((float(Mw)/PI) + 1)
		else:
			Instability_index = 0.0
			PI=0.0
			Gravy=0.0
			Mw = 0.0
			pI_Mw = 0.0
		return(Instability_index,PI,Gravy,Mw,pI_Mw)


	def get_Instab(self):

		insta_fe_all = []
		PI_fe_all = []
		gra_fe_all =[]
		Mw_all = []
		pI_Mw_all = []
		seqname = []
		for seq in SeqIO.parse(self.infasta,'fasta'):
			seqid = seq.id
			seqname.append(seqid)

			insta_fe,PI_fe,gra_fe,Mw,pI_Mw = self.param(seq.seq)
			insta_fe_all.append(insta_fe)
			PI_fe_all.append(PI_fe)
			gra_fe_all.append(gra_fe)
			Mw_all.append(Mw)
			pI_Mw_all.append(pI_Mw)


		insta_fe_all = np.expand_dims(np.array(insta_fe_all),axis = 1)
		PI_fe_all = np.expand_dims(np.array(PI_fe_all),axis = 1)
		gra_fe_all = np.expand_dims(np.array(gra_fe_all),axis = 1)
		Mw_all = np.expand_dims(np.array(Mw_all),axis = 1)
		pI_Mw_all = np.expand_dims(np.array(pI_Mw_all),axis = 1)

		colna = (["Instability"])
		all_ = np.concatenate((insta_fe_all,PI_fe_all,gra_fe_all,Mw_all,pI_Mw_all),axis=1)
		df_cod = DataFrame(data=insta_fe_all, index=seqname, columns=colna)

		return df_cod
	def get_PI(self):

		insta_fe_all = []
		PI_fe_all = []
		gra_fe_all =[]
		Mw_all = []
		pI_Mw_all = []
		seqname = []

		for seq in SeqIO.parse(self.infasta,'fasta'):
			seqid = seq.id
			seqname.append(seqid)

			insta_fe,PI_fe,gra_fe,Mw,pI_Mw = self.param(seq.seq)
			insta_fe_all.append(insta_fe)
			PI_fe_all.append(PI_fe)
			gra_fe_all.append(gra_fe)
			Mw_all.append(Mw)
			pI_Mw_all.append(pI_Mw)


		insta_fe_all = np.expand_dims(np.array(insta_fe_all),axis = 1)
		PI_fe_all = np.expand_dims(np.array(PI_fe_all),axis = 1)
		gra_fe_all = np.expand_dims(np.array(gra_fe_all),axis = 1)
		Mw_all = np.expand_dims(np.array(Mw_all),axis = 1)
		pI_Mw_all = np.expand_dims(np.array(pI_Mw_all),axis = 1)

		colna = (["Instability",'isoelectric_point'])
		all_ = np.concatenate((insta_fe_all,PI_fe_all,gra_fe_all,Mw_all,pI_Mw_all),axis=1)
		df_cod = DataFrame(data=PI_fe_all, index=seqname, columns=colna)

		return df_cod
	def get_Grav(self):

		insta_fe_all = []
		PI_fe_all = []
		gra_fe_all =[]
		Mw_all = []
		pI_Mw_all = []
		seqname = []

		for seq in SeqIO.parse(self.infasta,'fasta'):
			seqid = seq.id
			seqname.append(seqid)

			insta_fe,PI_fe,gra_fe,Mw,pI_Mw = self.param(seq.seq)
			insta_fe_all.append(insta_fe)
			PI_fe_all.append(PI_fe)
			gra_fe_all.append(gra_fe)
			Mw_all.append(Mw)
			pI_Mw_all.append(pI_Mw)


		insta_fe_all = np.expand_dims(np.array(insta_fe_all),axis = 1)
		PI_fe_all = np.expand_dims(np.array(PI_fe_all),axis = 1)
		gra_fe_all = np.expand_dims(np.array(gra_fe_all),axis = 1)
		Mw_all = np.expand_dims(np.array(Mw_all),axis = 1)
		pI_Mw_all = np.expand_dims(np.array(pI_Mw_all),axis = 1)

		colna = (["Instability",'isoelectric_point','Gravy'])
		all_ = np.concatenate((insta_fe_all,PI_fe_all,gra_fe_all,Mw_all,pI_Mw_all),axis=1)
		df_cod = DataFrame(data=gra_fe_all, index=seqname, columns=colna)

		return df_cod
	def get_MW(self):

		insta_fe_all = []
		PI_fe_all = []
		gra_fe_all =[]
		Mw_all = []
		pI_Mw_all = []
		seqname = []

		for seq in SeqIO.parse(self.infasta,'fasta'):
			seqid = seq.id
			seqname.append(seqid)

			insta_fe,PI_fe,gra_fe,Mw,pI_Mw = self.param(seq.seq)
			insta_fe_all.append(insta_fe)
			PI_fe_all.append(PI_fe)
			gra_fe_all.append(gra_fe)
			Mw_all.append(Mw)
			pI_Mw_all.append(pI_Mw)


		insta_fe_all = np.expand_dims(np.array(insta_fe_all),axis = 1)
		PI_fe_all = np.expand_dims(np.array(PI_fe_all),axis = 1)
		gra_fe_all = np.expand_dims(np.array(gra_fe_all),axis = 1)
		Mw_all = np.expand_dims(np.array(Mw_all),axis = 1)
		pI_Mw_all = np.expand_dims(np.array(pI_Mw_all),axis = 1)

		colna = (["Instability",'isoelectric_point','Gravy','Molecular_weight'])
		all_ = np.concatenate((insta_fe_all,PI_fe_all,gra_fe_all,Mw_all,pI_Mw_all),axis=1)
		df_cod = DataFrame(data=Mw_all, index=seqname, columns=colna)

		return df_cod
	def get_protper(self):

		insta_fe_all = []
		PI_fe_all = []
		gra_fe_all =[]
		Mw_all = []
		pI_Mw_all = []
		seqname = []
		# n = 0
		for seq in SeqIO.parse(self.infasta,'fasta'):
			seqid = seq.id
			seqname.append(seqid)
			# n += 1
			# print('n')
			# print(n)
			insta_fe,PI_fe,gra_fe,Mw,pI_Mw = self.param(seq.seq)
			insta_fe_all.append(insta_fe)
			PI_fe_all.append(PI_fe)
			gra_fe_all.append(gra_fe)
			Mw_all.append(Mw)
			pI_Mw_all.append(pI_Mw)


		insta_fe_all = np.expand_dims(np.array(insta_fe_all),axis = 1)
		PI_fe_all = np.expand_dims(np.array(PI_fe_all),axis = 1)
		gra_fe_all = np.expand_dims(np.array(gra_fe_all),axis = 1)
		Mw_all = np.expand_dims(np.array(Mw_all),axis = 1)
		pI_Mw_all = np.expand_dims(np.array(pI_Mw_all),axis = 1)
		# colna = (["Instability",'isoelectric_point','Gravy','Molecular_weight','pI/Mw_frame_score'])
		# colna = (["II", 'PI', 'GA', 'MW', 'IM'])
		colna = (["ProII: Pseudo protein instability index", 'ProPI: Pseudo protein isoelectric point', 'ProAH: Pseudo protein average hydropathy', 'ProMW: Pseudo protein molecular weight', 'PPMFS: Pseudo protein PI-MW frame score'])

		all_ = np.concatenate((insta_fe_all,PI_fe_all,gra_fe_all,Mw_all,pI_Mw_all),axis=1)
		df_cod = DataFrame(data=all_, index=seqname, columns=colna)

		return df_cod


#########################################################################
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

