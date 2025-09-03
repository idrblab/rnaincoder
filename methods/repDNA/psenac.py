__author__ = 'Fule Liu'

from util import get_data
from psenacutil import extend_phyche_index
from Bio import SeqIO
from pandas import DataFrame

def check_psenac(lamada, w, k):
    """Check the validation of parameter lamada, w and k.
    """
    try:
        if not isinstance(lamada, int) or lamada <= 0:
            raise ValueError("Error, parameter lamada must be an int type and larger than and equal to 0.")
        elif w > 1 or w < 0:
            raise ValueError("Error, parameter w must be ranged from 0 to 1.")
        elif not isinstance(k, int) or k <= 0:
            raise ValueError("Error, parameter k must be an int type and larger than 0.")
    except ValueError:
        raise


def get_sequence_list_and_phyche_value_psednc(input_data, extra_phyche_index=None):
    """For PseDNC, PseKNC, make sequence_list and phyche_value.

    :param input_data: file type or handle.
    :param extra_phyche_index: dict, the key is the dinucleotide (string),
                                     the value is its physicochemical property value (list).
                               It means the user-defined physicochemical indices.
    """
    if extra_phyche_index is None:
        extra_phyche_index = {}

    original_phyche_value = {'AA': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11],
                             'AC': [1.50, 0.50, 0.80, 0.13, 1.29, 1.04],
                             'AG': [0.78, 0.36, 0.09, 0.68, -0.24, -0.62],
                             'AT': [1.07, 0.22, 0.62, -1.02, 2.51, 1.17],
                             'CA': [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25],
                             'CC': [0.06, 1.08, 0.09, 0.56, -0.82, 0.24],
                             'CG': [-1.66, -1.22, -0.44, -0.82, -0.29, -1.39],
                             'CT': [0.78, 0.36, 0.09, 0.68, -0.24, -0.62],
                             'GA': [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71],
                             'GC': [-0.08, 0.22, 1.33, -0.35, 0.65, 1.59],
                             'GG': [0.06, 1.08, 0.09, 0.56, -0.82, 0.24],
                             'GT': [1.50, 0.50, 0.80, 0.13, 1.29, 1.04],
                             'TA': [-1.23, -2.37, -0.44, -2.24, -1.51, -1.39],
                             'TC': [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71],
                             'TG': [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25],
                             'TT': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11]}

    sequence_list = get_data(input_data)
    phyche_value = extend_phyche_index(original_phyche_value, extra_phyche_index)

    return sequence_list, phyche_value


def get_sequence_list_and_phyche_value_pseknc(input_data, extra_phyche_index=None):
    """For PseDNC, PseKNC, make sequence_list and phyche_value.

    :param input_data: file type or handle.
    :param extra_phyche_index: dict, the key is the dinucleotide (string),
                                     the value is its physicochemical property value (list).
                               It means the user-defined physicochemical indices.
    """
    if extra_phyche_index is None:
        extra_phyche_index = {}

    original_phyche_value = {
        'AA': [0.06, 0.5, 0.09, 1.59, 0.11, -0.11],
        'AC': [1.5, 0.5, 1.19, 0.13, 1.29, 1.04],
        'GT': [1.5, 0.5, 1.19, 0.13, 1.29, 1.04],
        'AG': [0.78, 0.36, -0.28, 0.68, -0.24, -0.62],
        'CC': [0.06, 1.08, -0.28, 0.56, -0.82, 0.24],
        'CA': [-1.38, -1.36, -1.01, -0.86, -0.62, -1.25],
        'CG': [-1.66, -1.22, -1.38, -0.82, -0.29, -1.39],
        'TT': [0.06, 0.5, 0.09, 1.59, 0.11, -0.11],
        'GG': [0.06, 1.08, -0.28, 0.56, -0.82, 0.24],
        'GC': [-0.08, 0.22, 2.3, -0.35, 0.65, 1.59],
        'AT': [1.07, 0.22, 0.83, -1.02, 2.51, 1.17],
        'GA': [-0.08, 0.5, 0.09, 0.13, -0.39, 0.71],
        'TG': [-1.38, -1.36, -1.01, -0.86, -0.62, -1.25],
        'TA': [-1.23, -2.37, -1.38, -2.24, -1.51, -1.39],
        'TC': [-0.08, 0.5, 0.09, 0.13, -0.39, 0.71],
        'CT': [0.78, 0.36, -0.28, 0.68, -0.24, -0.62]}

    sequence_list = get_data(input_data)
    phyche_value = extend_phyche_index(original_phyche_value, extra_phyche_index)

    return sequence_list, phyche_value


def get_sequence_list_and_phyche_value(input_data, k, phyche_index, extra_phyche_index, all_property):
    """For PseKNC-general make sequence_list and phyche_value.

    :param input_data: file type or handle.
    :param k: int, the value of k-tuple.
    :param k: physicochemical properties list.
    :param extra_phyche_index: dict, the key is the dinucleotide (string),
                                     the value is its physicochemical property value (list).
                               It means the user-defined physicochemical indices.
    :param all_property: bool, choose all physicochemical properties or not.
    """
    if phyche_index is None:
        phyche_index = []
    if extra_phyche_index is None:
        extra_phyche_index = {}

    diphyche_list = ['Base stacking', 'Protein induced deformability', 'B-DNA twist', 'Dinucleotide GC Content',
                     'A-philicity', 'Propeller twist', 'Duplex stability:(freeenergy)',
                     'Duplex tability(disruptenergy)', 'DNA denaturation', 'Bending stiffness', 'Protein DNA twist',
                     'Stabilising energy of Z-DNA', 'Aida_BA_transition', 'Breslauer_dG', 'Breslauer_dH',
                     'Breslauer_dS', 'Electron_interaction', 'Hartman_trans_free_energy', 'Helix-Coil_transition',
                     'Ivanov_BA_transition', 'Lisser_BZ_transition', 'Polar_interaction', 'SantaLucia_dG',
                     'SantaLucia_dH', 'SantaLucia_dS', 'Sarai_flexibility', 'Stability', 'Stacking_energy',
                     'Sugimoto_dG', 'Sugimoto_dH', 'Sugimoto_dS', 'Watson-Crick_interaction', 'Twist', 'Tilt', 'Roll',
                     'Shift', 'Slide', 'Rise']
    triphyche_list = ['Dnase I', 'Bendability (DNAse)', 'Bendability (consensus)', 'Trinucleotide GC Content',
                      'Nucleosome positioning', 'Consensus_roll', 'Consensus-Rigid', 'Dnase I-Rigid', 'MW-Daltons',
                      'MW-kg', 'Nucleosome', 'Nucleosome-Rigid']

    # Set and check physicochemical properties.
    phyche_list = []
    if k == 2:
        phyche_list = diphyche_list
    elif k == 3:
        phyche_list = triphyche_list

    try:
        if all_property is True:
            phyche_index = phyche_list
        else:
            for e in phyche_index:
                if e not in phyche_list:
                    error_info = 'Sorry, the physicochemical properties ' + e + ' is not exit.'
                    raise NameError(error_info)
    except NameError:
        raise

    # Generate phyche_value and sequence_list.
    from psenacutil import get_phyche_index

    phyche_value = extend_phyche_index(get_phyche_index(k, phyche_index), extra_phyche_index)
    sequence_list = get_data(input_data)

    return sequence_list, phyche_value


class PseDNC():
    def __init__(self, lamada=3, w=0.05):
        self.lamada = lamada
        self.w = w
        self.k = 2
        check_psenac(self.lamada, self.w, self.k)

    def make_psednc_vec(self, input_data, extra_phyche_index=None):
        """Make PseDNC vector.

        :param input_data: file type or handle.
        :param extra_phyche_index: dict, the key is the dinucleotide (string),
                                         the value is its physicochemical property value (list).
                                   It means the user-defined physicochemical indices.
        """
        sequence_list, phyche_value = get_sequence_list_and_phyche_value_psednc(input_data, extra_phyche_index)
        from psenacutil import make_pseknc_vector

        vector = make_pseknc_vector(sequence_list, self.lamada, self.w, self.k, phyche_value, theta_type=1)

        return vector


class PseKNC():
    """This class should be used to make PseKNC vector."""

    def __init__(self, k=3, lamada=1, w=0.5):
        """
        :param k: k-tuple.
        """
        self.k = k
        self.lamada = lamada
        self.w = w
        check_psenac(self.lamada, self.w, self.k)

    def make_pseknc_vec(self, input_data, extra_phyche_index=None):
        """Make PseKNC vector.

        :param input_data: file type or handle.
        :param extra_phyche_index: dict, the key is the dinucleotide (string),
                                         the value is its physicochemical property value (list).
                                   It means the user-defined physicochemical indices.
        """
        sequence_list, phyche_value = get_sequence_list_and_phyche_value_pseknc(input_data, extra_phyche_index)
        from psenacutil import make_old_pseknc_vector

        return make_old_pseknc_vector(sequence_list, self.lamada, self.w, self.k, phyche_value, theta_type=1)


class PCPseDNC():
    def __init__(self, lamada=1, w=0.05):
        self.lamada = lamada
        self.w = w
        self.k = 2
        check_psenac(self.lamada, self.w, self.k)

    def make_pcpsednc_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        """Make a PCPseDNC vector.

        :param input_data: file object or sequence list.
        :param phyche_index: physicochemical properties list.
        :param all_property: choose all physicochemical properties or not.
        :param extra_phyche_index: dict, the key is the dinucleotide (string),
                                         the value is its physicochemical property value (list).
                                   It means the user-defined physicochemical indices.
        """
        # Make vector.
        sequence_list, phyche_value = get_sequence_list_and_phyche_value(input_data, self.k, phyche_index,
                                                                         extra_phyche_index, all_property)
        from psenacutil import make_pseknc_vector

        vector = make_pseknc_vector(sequence_list, self.lamada, self.w, self.k, phyche_value, theta_type=1)

        return vector


class PCPseTNC():
    def __init__(self, lamada=1, w=0.05):
        self.lamada = lamada
        self.w = w
        self.k = 3
        check_psenac(self.lamada, self.w, self.k)

    def make_pcpsetnc_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        """Make a PCPseDNC vector.

        :param input_data: file object or sequence list.
        :param phyche_index: physicochemical properties list.
        :param all_property: choose all physicochemical properties or not.
        :param extra_phyche_index: dict, the key is the dinucleotide (string),
                                         the value is its physicochemical property value (list).
                                   It means the user-defined physicochemical indices.
        """
        sequence_list, phyche_value = get_sequence_list_and_phyche_value(input_data, self.k, phyche_index,
                                                                         extra_phyche_index, all_property)
        # Make vector.
        from psenacutil import make_pseknc_vector

        vector = make_pseknc_vector(sequence_list, self.lamada, self.w, self.k, phyche_value, theta_type=1)

        return vector


class SCPseDNC():
    def __init__(self, lamada=1, w=0.05):
        self.lamada = lamada
        self.w = w
        self.k = 2
        check_psenac(self.lamada, self.w, self.k)

    def make_scpsednc_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        """Make a SCPseDNC vector.

        :param input_data: file object or sequence list.
        :param phyche_index: physicochemical properties list.
        :param all_property: choose all physicochemical properties or not.
        :param extra_phyche_index: dict, the key is the dinucleotide (string),
                                         the value is its physicochemical property value (list).
                                   It means the user-defined physicochemical indices.
        """
        sequence_list, phyche_value = get_sequence_list_and_phyche_value(input_data, self.k, phyche_index,
                                                                         extra_phyche_index, all_property)
        # Make vector.
        from psenacutil import make_pseknc_vector

        vector = make_pseknc_vector(sequence_list, self.lamada, self.w, self.k, phyche_value, theta_type=2)

        return vector


class SCPseTNC():
    def __init__(self, lamada=1, w=0.05):
        self.lamada = lamada
        self.w = w
        self.k = 3
        check_psenac(self.lamada, self.w, self.k)

    def make_scpsetnc_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        """Make a SCPseTNC vector.

        :param input_data: file object or sequence list.
        :param phyche_index: physicochemical properties list.
        :param all_property: choose all physicochemical properties or not.
        :param extra_phyche_index: dict, the key is the dinucleotide (string),
                                         the value is its physicochemical property value (list).
                                   It means the user-defined physicochemical indices.
        """
        sequence_list, phyche_value = get_sequence_list_and_phyche_value(input_data, self.k, phyche_index,
                                                                         extra_phyche_index, all_property)
        # Make vector.
        from psenacutil import make_pseknc_vector

        vector = make_pseknc_vector(sequence_list, self.lamada, self.w, self.k, phyche_value, theta_type=2)

        return vector
import numpy as np
def rna_pc_psednc(filepath):
    seq_seq = []
    seqname = []
    for seq in SeqIO.parse(filepath, 'fasta'):
        seq_seq.append(seq.seq)
        seqid = seq.id
        seqname.append(seqid)

    pc_psednc = PCPseDNC(lamada=2, w=0.05)
    
    vec_all = []

    for seq in seq_seq:
        seq = str(seq)

        vec = pc_psednc.make_pcpsednc_vec([seq], all_property=True)[0]
        vec_all.append(vec)

    vec_all = np.array(vec_all)
    colname = ['PCDAA: Parallel correlation PseDN composition of AA','PCDAT: Parallel correlation PseDN composition of AT','PCDAC: Parallel correlation PseDN composition of AC','PCDAG: Parallel correlation PseDN composition of AG','PCDTA: Parallel correlation PseDN composition of TA','PCDTT: Parallel correlation PseDN composition of TT','PCDTC: Parallel correlation PseDN composition of TC','PCDTG: Parallel correlation PseDN composition of TG','PCDCA: Parallel correlation PseDN composition of CA','PCDCT: Parallel correlation PseDN composition of CT','PCDCC: Parallel correlation PseDN composition of CC','PCDCG: Parallel correlation PseDN composition of CG','PCDGA: Parallel correlation PseDN composition of GA','PCDGT: Parallel correlation PseDN composition of GT','PCDGC: Parallel correlation PseDN composition of GC','PCDGG: Parallel correlation PseDN composition of GG','PCDl1: Parallel correlation PseDN composition of lamda1','PCDl2: Parallel correlation PseDN composition of lamda2']
    df = DataFrame(data=vec_all, index=seqname,columns=colname)
    return df

def rna_SCPseDNC(filepath):
    seq_seq = []
    seqname = []
    for seq in SeqIO.parse(filepath, 'fasta'):
        seq_seq.append(seq.seq)
        seqid = seq.id
        seqname.append(seqid)
    sc_psednc = SCPseDNC(lamada=2, w=0.05)
    vec_all = []
    for seq in seq_seq:
        seq = str(seq)

        vec = sc_psednc.make_scpsednc_vec([seq], phyche_index=['Protein induced deformability', 'Bending stiffness','Electron_interaction','Stacking_energy','Watson-Crick_interaction','Slide'])[0]
        vec_all.append(vec)
    vec_all = np.array(vec_all)
    colname = ['SCDAA: Series correlation PseDN composition of AA','SCDAT: Series correlation PseDN composition of AT','SCDAC: Series correlation PseDN composition of AC','SCDAG: Series correlation PseDN composition of AG','SCDTA: Series correlation PseDN composition of TA','SCDTT: Series correlation PseDN composition of TT','SCDTC: Series correlation PseDN composition of TC','SCDTG: Series correlation PseDN composition of TG','SCDCA: Series correlation PseDN composition of CA','SCDCT: Series correlation PseDN composition of CT','SCDCC: Series correlation PseDN composition of CC','SCDCG: Series correlation PseDN composition of CG','SCDGA: Series correlation PseDN composition of GA','SCDGT: Series correlation PseDN composition of GT','SCDGC: Series correlation PseDN composition of GC','SCDGG: Series correlation PseDN composition of GG','L1PID: Series correlation PseDN composition of lamda1-PID','L1BDS: Series correlation PseDN composition of lamda1-BDS','L1ETI: Series correlation PseDN composition of lamda1-ETI','L1SKE: Series correlation PseDN composition of lamda1-SKE','L1WCI: Series correlation PseDN composition of lamda1-WCI','L1SLI: Series correlation PseDN composition of lamda1-SLI','L2PID: Series correlation PseDN composition of lamda2-PID','L2BDS: Series correlation PseDN composition of lamda2-BDS','L2ETI: Series correlation PseDN composition of lamda2-ETI','L2SKE: Series correlation PseDN composition of lamda2-SKE','L2WCI: Series correlation PseDN composition of lamda2-WCI','L2SLI: Series correlation PseDN composition of lamda2-SLI']
    df = DataFrame(data=vec_all, index=seqname,columns=colname)
    return df

