
import numpy as np 
from Bio import SeqIO
group_dict = {
	'A':[1, 1, 1],
	'T':[0, 0, 1],
	'C':[1, 0, 0],
	'G':[0, 1, 0]
}

def encoding(sequence):
	seq_array = []
	for i in sequence:
		map_seq = group_dict[i]
		seq_array.append(map_seq)
	seq_array = np.array(seq_array)
	return seq_array


def get_encoding(filepath, N):
	seq_seq = []
	seqname = []
	for seq in SeqIO.parse(filepath, 'fasta'):
		seq_seq.append(seq.seq)
		seqname.append(seq.id)
	encoding_array = np.zeros((len(seq_seq), N, 3))

	for i, seq in enumerate(seq_seq):
		seq = str(seq)
		seq_array = encoding(seq)
		if len(seq) >= N:
			encoding_array[i, :, :] = seq_array[0:N]
		if len(seq) < N:
			encoding_array[i, 0:len(seq), :] = seq_array
	return seqname,encoding_array

