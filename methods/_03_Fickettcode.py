#!/usr/bin/env python
'''the python script is downloaded from https://sourceforge.net/projects/rna-cpat/files/v1.2.2/'''
'''calculate coding potential'''

# Fickett TESTCODE data
# NAR 10(17) 5303-531
import Bio.SeqIO as Seq
import numpy as np
from pandas import DataFrame

class Fickettcoder:
	def __init__(self, infasta=None):
		self.infasta = infasta

		self.position_prob ={
		'A':[0.94,0.68,0.84,0.93,0.58,0.68,0.45,0.34,0.20,0.22],
		'C':[0.80,0.70,0.70,0.81,0.66,0.48,0.51,0.33,0.30,0.23],
		'G':[0.90,0.88,0.74,0.64,0.53,0.48,0.27,0.16,0.08,0.08],
		'T':[0.97,0.97,0.91,0.68,0.69,0.44,0.54,0.20,0.09,0.09]
		}
		self.position_weight={'A':0.26,'C':0.18,'G':0.31,'T':0.33}
		self.position_para  =[1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,0.0]

		self.content_prob={
		'A':[0.28,0.49,0.44,0.55,0.62,0.49,0.67,0.65,0.81,0.21],
		'C':[0.82,0.64,0.51,0.64,0.59,0.59,0.43,0.44,0.39,0.31],
		'G':[0.40,0.54,0.47,0.64,0.64,0.73,0.41,0.41,0.33,0.29],
		'T':[0.28,0.24,0.39,0.40,0.55,0.75,0.56,0.69,0.51,0.58]
		}
		self.content_weight={'A':0.11,'C':0.12,'G':0.15,'T':0.14}
		self.content_para  =[0.33,0.31,0.29,0.27,0.25,0.23,0.21,0.17,0]

	def look_up_position_prob(self, value, base):
		'''look up positional probability by base and value'''
		if float(value)<0:
			return None
		for idx,val in enumerate (self.position_para):
			if (float(value) >= val):
				return float(self.position_prob[base][idx]) * float(self.position_weight[base])

	def look_up_content_prob(self, value, base):
		'''look up content probability by base and value'''
		if float(value)<0:
			return None
		for idx,val in enumerate (self.content_para):
			if (float(value) >= val):
				return float(self.content_prob[base][idx]) * float(self.content_weight[base])

	def fickett_value(self, dna):
		'''calculate Fickett value. Input is DNA sequence'''
		if len(dna)<2:
			return 0
		fickett_score=0
		dna=dna.upper()
		total_base = len(dna)
		A_content = float(dna.count('A'))/total_base
		C_content = float(dna.count('C'))/total_base
		G_content = float(dna.count('G'))/total_base
		T_content = float(dna.count('T'))/total_base
		#print "A content\t" + str(A_content)
		#print "C content\t" + str(C_content)
		#print "G content\t" + str(G_content)
		#print "T content\t" + str(T_content)
		
		phase_0 = [dna[i] for i in range(0,len(dna)) if i % 3==0]
		phase_1 = [dna[i] for i in range(0,len(dna)) if i % 3==1]
		phase_2 = [dna[i] for i in range(0,len(dna)) if i % 3==2]
		
		A_position=max(phase_0.count('A'),phase_1.count('A'),phase_2.count('A'))/(min(phase_0.count('A'),phase_1.count('A'),phase_2.count('A')) +1.0)
		C_position=max(phase_0.count('C'),phase_1.count('C'),phase_2.count('C'))/(min(phase_0.count('C'),phase_1.count('C'),phase_2.count('C')) +1.0)
		G_position=max(phase_0.count('G'),phase_1.count('G'),phase_2.count('G'))/(min(phase_0.count('G'),phase_1.count('G'),phase_2.count('G')) +1.0)
		T_position=max(phase_0.count('T'),phase_1.count('T'),phase_2.count('T'))/(min(phase_0.count('T'),phase_1.count('T'),phase_2.count('T')) +1.0)
		#print "A position\t" + str(A_position)
		#print "C position\t" + str(C_position)
		#print "G position\t" + str(G_position)
		#print "T position\t" + str(T_position)

		
		#for i (A_content,C_content,G_content,T_content):
		fickett_score += self.look_up_content_prob(A_content,'A')
		fickett_score += self.look_up_content_prob(C_content,'C')
		fickett_score += self.look_up_content_prob(G_content,'G')
		fickett_score += self.look_up_content_prob(T_content,'T')
		
		fickett_score += self.look_up_position_prob(A_position,'A')
		fickett_score += self.look_up_position_prob(C_position,'C')
		fickett_score += self.look_up_position_prob(G_position,'G')
		fickett_score += self.look_up_position_prob(T_position,'T')

		return fickett_score

	# def fickett_value(self, dna):
	# 	'''calculate Fickett value. Input is DNA sequence'''
	# 	if len(dna)<2:
	# 		return [0,0,0,0,0,0,0,0]
	# 	fickett_score=0
	# 	dna=dna.upper()
	# 	total_base = len(dna)
	# 	A_content = float(dna.count('A'))/total_base
	# 	C_content = float(dna.count('C'))/total_base
	# 	G_content = float(dna.count('G'))/total_base
	# 	T_content = float(dna.count('T'))/total_base
	# 	#print "A content\t" + str(A_content)
	# 	#print "C content\t" + str(C_content)
	# 	#print "G content\t" + str(G_content)
	# 	#print "T content\t" + str(T_content)
		
	# 	phase_0 = [dna[i] for i in range(0,len(dna)) if i % 3==0]
	# 	phase_1 = [dna[i] for i in range(0,len(dna)) if i % 3==1]
	# 	phase_2 = [dna[i] for i in range(0,len(dna)) if i % 3==2]
		
	# 	A_position=max(phase_0.count('A'),phase_1.count('A'),phase_2.count('A'))/(min(phase_0.count('A'),phase_1.count('A'),phase_2.count('A')) +1.0)
	# 	C_position=max(phase_0.count('C'),phase_1.count('C'),phase_2.count('C'))/(min(phase_0.count('C'),phase_1.count('C'),phase_2.count('C')) +1.0)
	# 	G_position=max(phase_0.count('G'),phase_1.count('G'),phase_2.count('G'))/(min(phase_0.count('G'),phase_1.count('G'),phase_2.count('G')) +1.0)
	# 	T_position=max(phase_0.count('T'),phase_1.count('T'),phase_2.count('T'))/(min(phase_0.count('T'),phase_1.count('T'),phase_2.count('T')) +1.0)
	# 	#print "A position\t" + str(A_position)
	# 	#print "C position\t" + str(C_position)
	# 	#print "G position\t" + str(G_position)
	# 	#print "T position\t" + str(T_position)

		
	# 	#for i (A_content,C_content,G_content,T_content):
	# 	content_A_prob = self.look_up_content_prob(A_content,'A')
	# 	content_C_prob = self.look_up_content_prob(C_content,'C')
	# 	content_G_prob = self.look_up_content_prob(G_content,'G')
	# 	content_T_prob = self.look_up_content_prob(T_content,'T')
		
	# 	position_A_prob = self.look_up_position_prob(A_position,'A')
	# 	position_C_prob = self.look_up_position_prob(C_position,'C')
	# 	position_G_prob = self.look_up_position_prob(G_position,'G')
	# 	position_T_prob = self.look_up_position_prob(T_position,'T')

	# 	return [content_A_prob,content_C_prob,content_G_prob,content_T_prob,position_A_prob,position_C_prob,position_G_prob,position_T_prob]


	def get_fickett(self):


		#feaname = (["Fickett"])
		# feaname = (["FS"])
		# feaname = (["FickS: Fickett score"])
		feaname = (["FickS: Fickett score"])
		# seqname = (["Aconp: A Content Probability","Cconp: C Content Probability","Gconp: G Content Probability","Tconp: T Content Probability","Aposp: A Position Probability","Cposp: C Position Probability","Gposp: G Position Probability","Tposp: T Position Probability"])
		fickett_all = []
		seqname = []
		for seq in Seq.parse(self.infasta,'fasta'):
			seqid = seq.id
			seqname.append(seqid)
			fickett_fe = self.fickett_value(seq.seq)
			fickett_all.append(fickett_fe)


		fickett_all = np.array(fickett_all)
		df = DataFrame(data=fickett_all, index=seqname, columns=feaname)
		return df