from Bio import Seq
import Bio.SeqIO as SeqO
import numpy as np
from pandas import DataFrame

def _stop_codon_num(seq):
    # print(Seq)
    translate_prot = Seq.translate(seq)
    stop_num = translate_prot.count("*")
    return stop_num

def get_length(seq):
    return np.log(len(seq)+1)

def _stop_codon_frequency(seq):
    stop_num = _stop_codon_num(seq)
    transript_length = get_length(seq)
    stop_freq = float(stop_num) / transript_length
    return stop_freq

def _stop_num_frame_score(seq):
    stop_num_in_frame1 = _stop_codon_num(seq)
    stop_num_in_frame2 = _stop_codon_num(seq[1:])
    stop_num_in_frame3 = _stop_codon_num(seq[2:])
    stop_num_all = [stop_num_in_frame1, stop_num_in_frame2, stop_num_in_frame3]
    stop_num_frame = ((stop_num_all[0] - stop_num_all[1]) ** 2 + (stop_num_all[0] - stop_num_all[2]) ** 2 + (
            stop_num_all[1] - stop_num_all[2]) ** 2) / 2
    return stop_num_frame

def _stop_frequency_frame_score(seq):
    stop_num_in_frame1 = _stop_codon_frequency(seq)
    stop_num_in_frame2 = _stop_codon_frequency(seq[1:])
    stop_num_in_frame3 = _stop_codon_frequency(seq[2:])
    stop_num_all = [stop_num_in_frame1, stop_num_in_frame2, stop_num_in_frame3]
    stop_num_frame = ((stop_num_all[0] - stop_num_all[1]) ** 2 + (stop_num_all[0] - stop_num_all[2]) ** 2 + (
            stop_num_all[1] - stop_num_all[2]) ** 2) / 2
    return stop_num_frame

def get_stop_codon_num(infasta):

    feaname = (["SCCou: Stop codon count"])
    seqname = []
    stop_codon_num_all = []

    for seq in SeqO.parse(infasta,'fasta'):
        seqid = seq.id
        seqname.append(seqid)
        stop_codon_num = _stop_codon_num(seq.seq)
        stop_codon_num_all.append(stop_codon_num)



    stop_codon_num_all = np.array(stop_codon_num_all)
    df = DataFrame(data=stop_codon_num_all, index=seqname, columns=feaname)
    return df

def get_stop_codon_frequency(infasta):

    # feaname = (["stop_codon_frequency"])
    feaname = (["SCFre: Stop codon frequency"])

    seqname = []
    stop_codon_frequency_all = []

    for seq in SeqO.parse(infasta,'fasta'):
        seqid = seq.id
        seqname.append(seqid)
        stop_codon_frequency = _stop_codon_frequency(seq.seq)
        stop_codon_frequency_all.append(stop_codon_frequency)


    stop_codon_frequency_all = np.array(stop_codon_frequency_all)
    df = DataFrame(data=stop_codon_frequency_all, index=seqname, columns=feaname)
    return df

def get_stop_frequency_frame_score(infasta):

    feaname = (["SCFFS: Stop codon frequency frame score"])
    seqname = []
    stop_codon_frequency_all = []

    for seq in SeqO.parse(infasta,'fasta'):
        seqid = seq.id
        seqname.append(seqid)
        stop_codon_frequency = _stop_frequency_frame_score(seq.seq)
        stop_codon_frequency_all.append(stop_codon_frequency)


    stop_codon_frequency_all = np.array(stop_codon_frequency_all)
    df = DataFrame(data=stop_codon_frequency_all, index=seqname, columns=feaname)
    return df

def get_stop_num_frame_score(infasta):

    # feaname = (["stop_codon_num_score"])
    feaname = (["SCCFS: Stop codon count frame score"])
    seqname = []
    stop_codon_frequency_all = []

    for seq in SeqO.parse(infasta,'fasta'):
        seqid = seq.id
        seqname.append(seqid)
        stop_codon_frequency = _stop_num_frame_score(seq.seq)
        stop_codon_frequency_all.append(stop_codon_frequency)


    stop_codon_frequency_all = np.array(stop_codon_frequency_all)
    df = DataFrame(data=stop_codon_frequency_all, index=seqname, columns=feaname)
    return df

def get_stop(infasta):


    seqname = []
    stop_codon_num_all = []
    stop_codon_frequency_all = []
    stop_codon_numfram_all = []
    stop_codon_frefram_all = []

    for seq in SeqO.parse(infasta,'fasta'):
        seqid = seq.id
        seqname.append(seqid)
        stop_codon_num = _stop_codon_num(seq.seq)
        stop_codon_num_all.append(stop_codon_num)
        stop_codon_frequency = _stop_codon_frequency(seq.seq)
        stop_codon_frequency_all.append(stop_codon_frequency)

        stop_codon_numfram = _stop_num_frame_score(seq.seq)
        stop_codon_numfram_all.append(stop_codon_numfram)
        stop_codon_frequencyfram = _stop_frequency_frame_score(seq.seq)
        stop_codon_frefram_all.append(stop_codon_frequencyfram)


    stop_codon_num_all = np.expand_dims(np.array(stop_codon_num_all),axis = 1)
    stop_codon_frequency_all = np.expand_dims(np.array(stop_codon_frequency_all),axis = 1)
    stop_codon_numfram_all = np.expand_dims(np.array(stop_codon_numfram_all),axis = 1)
    stop_codon_frefram_all = np.expand_dims(np.array(stop_codon_frefram_all),axis = 1)



    stop_codon = np.concatenate((stop_codon_num_all, stop_codon_frequency_all,stop_codon_numfram_all,stop_codon_frefram_all), axis=1)

    # coname = (["stop_codon_num", "stop_codon_frequency", "stop_codon_num_score", "stop_codon_frequency_score"])
    coname = (["SCCou: Stop codon count", "SCFre: Stop codon frequency", "SCCFS: Stop codon count frame score", "SCFFS: Stop codon frequency frame score"])

    df = DataFrame(data=stop_codon, index=seqname, columns=coname)
    return df



