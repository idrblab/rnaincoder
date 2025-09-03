__author__ = 'Fule Liu'

from util import get_data, generate_phyche_value
from functools import reduce
from pandas import DataFrame

def check_acc(lag, k):
    """Check ACC parameter validation.
    """
    try:
        if not isinstance(lag, int) or lag <= 0:
            raise ValueError("Error, parameter lag must be an int type and larger than 0.")
        elif not isinstance(k, int) or lag <= 0:
            raise ValueError("Error, parameter k must be an int type and larger than 0.")
    except ValueError:
        raise


def ready_acc(input_data, k, phyche_index=None, all_property=False, extra_phyche_index=None):
    """Public function for get sequence_list and phyche_value.
    """
    sequence_list = get_data(input_data)
    if phyche_index is None:
        phyche_index = []
    if extra_phyche_index is None:
        extra_phyche_index = {}
    phyche_value = generate_phyche_value(k, phyche_index, all_property, extra_phyche_index)

    return sequence_list, phyche_value


class DAC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 2
        check_acc(self.lag, self.k)

    def make_dac_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        """Make DAC vector.

        :param input_data: file object or sequence list.
        :param phyche_index: physicochemical properties list.
        :param all_property: bool, choose all physicochemical properties or not.
        :param extra_phyche_index: dict, the key is the dinucleotide (string), and its corresponding value is a list.
                                   It means user-defined phyche_index.
        """
        sequence_list, phyche_value = ready_acc(input_data, self.k, phyche_index, all_property, extra_phyche_index)
        from acutil import make_ac_vector
        return make_ac_vector(sequence_list, self.lag, phyche_value, self.k)


class DCC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 2
        check_acc(self.lag, self.k)

    def make_dcc_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        """Make DCC vector.

        :param input_data: file object or sequence list.
        :param phyche_index: physicochemical properties list.
        :param all_property: bool, choose all physicochemical properties or not.
        :param extra_phyche_index: dict, the key is the dinucleotide (string), and its corresponding value is a list.
                                   It means user-defined phyche_index.
        """
        sequence_list, phyche_value = ready_acc(input_data, self.k, phyche_index, all_property, extra_phyche_index)
        from acutil import make_cc_vector
        return make_cc_vector(sequence_list, self.lag, phyche_value, self.k)


class DACC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 2
        check_acc(self.lag, self.k)

    def make_dacc_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        """Make DACC vector.

        :param input_data: file object or sequence list.
        :param phyche_index: physicochemical properties list.
        :param all_property: bool, choose all physicochemical properties or not.
        :param extra_phyche_index: dict, the key is the dinucleotide (string), and its corresponding value is a list.
                                   It means user-defined phyche_index.
        """
        sequence_list, phyche_value = ready_acc(input_data, self.k, phyche_index, all_property, extra_phyche_index)
        from acutil import make_ac_vector, make_cc_vector
        zipped = list(zip(make_ac_vector(sequence_list, self.lag, phyche_value, self.k),
                     make_cc_vector(sequence_list, self.lag, phyche_value, self.k)))
        vector = [reduce(lambda x, y: x + y, e) for e in zipped]

        return vector


class TAC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 3
        check_acc(self.lag, self.k)

    def make_tac_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        """Make TAC vector.

        :param input_data: file object or sequence list.
        :param phyche_index: physicochemical properties list.
        :param all_property: bool, choose all physicochemical properties or not.
        :param extra_phyche_index: dict, the key is the dinucleotide (string), and its corresponding value is a list.
                                   It means user-defined phyche_index.
        """
        sequence_list, phyche_value = ready_acc(input_data, self.k, phyche_index, all_property, extra_phyche_index)
        from acutil import make_ac_vector
        return make_ac_vector(sequence_list, self.lag, phyche_value, self.k)


class TCC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 3
        check_acc(self.lag, self.k)

    def make_tcc_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        """Make DAC vector.

        :param input_data: file object or sequence list.
        :param phyche_index: physicochemical properties list.
        :param all_property: bool, choose all physicochemical properties or not.
        :param extra_phyche_index: dict, the key is the dinucleotide (string), and its corresponding value is a list.
                                   It means user-defined phyche_index.
        """
        sequence_list, phyche_value = ready_acc(input_data, self.k, phyche_index, all_property, extra_phyche_index)
        from acutil import make_cc_vector
        return make_cc_vector(sequence_list, self.lag, phyche_value, self.k)


class TACC():
    def __init__(self, lag):
        self.lag = lag
        self.k = 3
        check_acc(self.lag, self.k)

    def make_tacc_vec(self, input_data, phyche_index=None, all_property=False, extra_phyche_index=None):
        """Make DAC vector.

        :param input_data: file object or sequence list.
        :param phyche_index: physicochemical properties list.
        :param all_property: bool, choose all physicochemical properties or not.
        :param extra_phyche_index: dict, the key is the dinucleotide (string), and its corresponding value is a list.
                                   It means user-defined phyche_index.
        """
        sequence_list, phyche_value = ready_acc(input_data, self.k, phyche_index, all_property, extra_phyche_index)

        from acutil import make_ac_vector, make_cc_vector

        zipped = list(zip(make_ac_vector(sequence_list, self.lag, phyche_value, self.k),
                     make_cc_vector(sequence_list, self.lag, phyche_value, self.k)))
        vector = [reduce(lambda x, y: x + y, e) for e in zipped]

        return vector
import numpy as np
from Bio import SeqIO
def rna_dac(filepath,lag):
    seq_seq = []
    seqname = []
    # print('This step is ran dac 001')
    for seq in SeqIO.parse(filepath, 'fasta'):
        seq_seq.append(seq.seq)
        seqid = seq.id
        seqname.append(seqid)
    # print('This step is ran dac 002')
    dac = DAC(lag)
    vec_all = []
    # print('This step is ran dac 003')
    # print(len(seq_seq))
    for seq in seq_seq:
        seq = str(seq)
        # print(seq)
        vec = dac.make_dac_vec([seq], all_property=True)[0]
        
        # print('This step is ran dac 004')
        # print(vec)
        vec_all.append(vec)
    vec_all = np.array(vec_all)
    # print('This step is ran dac 005')
    # print(vec_all.shape)
    colname = ['ABST1: Dinucleotide auto covariance of lag1 on BST',
               'APID1: Dinucleotide auto covariance of lag1 on PID',
               'ABDT1: Dinucleotide auto covariance of lag1 on BDT',
               'ADGC1: Dinucleotide auto covariance of lag1 on DGC',
               'AAPH1: Dinucleotide auto covariance of lag1 on APH',
               'APTW1: Dinucleotide auto covariance of lag1 on PTW',
               'ADST1: Dinucleotide auto covariance of lag1 on DST',
               'ADTB1: Dinucleotide auto covariance of lag1 on DTB',
               'ADDT1: Dinucleotide auto covariance of lag1 on DDT',
               'ABDS1: Dinucleotide auto covariance of lag1 on BDS',
               'APDT1: Dinucleotide auto covariance of lag1 on PDT',
               'ASEZ1: Dinucleotide auto covariance of lag1 on SEZ',
               'AABA1: Dinucleotide auto covariance of lag1 on ABA',
               'ABDG1: Dinucleotide auto covariance of lag1 on BDG',
               'ABDH1: Dinucleotide auto covariance of lag1 on BDH',
               'ABLS1: Dinucleotide auto covariance of lag1 on BLS','AETI1: Dinucleotide auto covariance of lag1 on ETI','AHTF1: Dinucleotide auto covariance of lag1 on HTF','AHCT1: Dinucleotide auto covariance of lag1 on HCT','AIBA1: Dinucleotide auto covariance of lag1 on IBA','ALBZ1: Dinucleotide auto covariance of lag1 on LBZ','APIA1: Dinucleotide auto covariance of lag1 on PIA','ASDG1: Dinucleotide auto covariance of lag1 on SDG','ASDH1: Dinucleotide auto covariance of lag1 on SDH','ASDS1: Dinucleotide auto covariance of lag1 on SDS','ASFB1: Dinucleotide auto covariance of lag1 on SFB','ASTB1: Dinucleotide auto covariance of lag1 on STB','ASKE1: Dinucleotide auto covariance of lag1 on SKE','ASMG1: Dinucleotide auto covariance of lag1 on SMG','ASMH1: Dinucleotide auto covariance of lag1 on SMH','ASMS1: Dinucleotide auto covariance of lag1 on SMS','AWCI1: Dinucleotide auto covariance of lag1 on WCI','ATWI1: Dinucleotide auto covariance of lag1 on TWI','ATIL1: Dinucleotide auto covariance of lag1 on TIL','AROL1: Dinucleotide auto covariance of lag1 on ROL','ASHI1: Dinucleotide auto covariance of lag1 on SHI','ASLI1: Dinucleotide auto covariance of lag1 on SLI','ARIS1: Dinucleotide auto covariance of lag1 on RIS','ABST2: Dinucleotide auto covariance of lag2 on BST','APID2: Dinucleotide auto covariance of lag2 on PID','ABDT2: Dinucleotide auto covariance of lag2 on BDT','ADGC2: Dinucleotide auto covariance of lag2 on DGC','AAPH2: Dinucleotide auto covariance of lag2 on APH','APTW2: Dinucleotide auto covariance of lag2 on PTW','ADST2: Dinucleotide auto covariance of lag2 on DST','ADTB2: Dinucleotide auto covariance of lag2 on DTB','ADDT2: Dinucleotide auto covariance of lag2 on DDT','ABDS2: Dinucleotide auto covariance of lag2 on BDS','APDT2: Dinucleotide auto covariance of lag2 on PDT','ASEZ2: Dinucleotide auto covariance of lag2 on SEZ','AABA2: Dinucleotide auto covariance of lag2 on ABA','ABDG2: Dinucleotide auto covariance of lag2 on BDG','ABDH2: Dinucleotide auto covariance of lag2 on BDH','ABLS2: Dinucleotide auto covariance of lag2 on BLS','AETI2: Dinucleotide auto covariance of lag2 on ETI','AHTF2: Dinucleotide auto covariance of lag2 on HTF','AHCT2: Dinucleotide auto covariance of lag2 on HCT','AIBA2: Dinucleotide auto covariance of lag2 on IBA','ALBZ2: Dinucleotide auto covariance of lag2 on LBZ','APIA2: Dinucleotide auto covariance of lag2 on PIA','ASDG2: Dinucleotide auto covariance of lag2 on SDG','ASDH2: Dinucleotide auto covariance of lag2 on SDH','ASDS2: Dinucleotide auto covariance of lag2 on SDS','ASFB2: Dinucleotide auto covariance of lag2 on SFB','ASTB2: Dinucleotide auto covariance of lag2 on STB','ASKE2: Dinucleotide auto covariance of lag2 on SKE','ASMG2: Dinucleotide auto covariance of lag2 on SMG','ASMH2: Dinucleotide auto covariance of lag2 on SMH','ASMS2: Dinucleotide auto covariance of lag2 on SMS','AWCI2: Dinucleotide auto covariance of lag2 on WCI','ATWI2: Dinucleotide auto covariance of lag2 on TWI','ATIL2: Dinucleotide auto covariance of lag2 on TIL','AROL2: Dinucleotide auto covariance of lag2 on ROL','ASHI2: Dinucleotide auto covariance of lag2 on SHI','ASLI2: Dinucleotide auto covariance of lag2 on SLI','ARIS2: Dinucleotide auto covariance of lag2 on RIS']

    df_cod = DataFrame(data=vec_all, index=seqname,columns=colname)
    return df_cod
    # return vec_all

def rna_dcc(filepath,lag):
    seq_seq = []
    seq_id = []
    for seq in SeqIO.parse(filepath, 'fasta'):
        seq_seq.append(seq.seq)
        seq_id.append(seq.id)
    dcc = DCC(lag)
    vec_all = []
    for seq in seq_seq:
        seq = str(seq)
        vec = dcc.make_dcc_vec([seq], phyche_index=['Protein induced deformability', 'Bending stiffness','Electron_interaction','Stacking_energy','Watson-Crick_interaction','Slide'])[0]


        # vec = dcc.make_dcc_vec([seq], all_property=True)[0]
        vec_all.append(vec)
    vec_all = np.array(vec_all)
    colname = ['CC101: Cross covariance of lag1 on PID-BDS','CC102: Cross covariance of lag1 on PID-ETI','CC103: Cross covariance of lag1 on PID-SKE','CC104: Cross covariance of lag1 on PID-WCI','CC105: Cross covariance of lag1 on PID-SLI','CC106: Cross covariance of lag1 on BDS-PID','CC107: Cross covariance of lag1 on BDS-ETI','CC108: Cross covariance of lag1 on BDS-SKE','CC109: Cross covariance of lag1 on BDS-WCI','CC110: Cross covariance of lag1 on BDS-SLI','CC111: Cross covariance of lag1 on ETI-PID','CC112: Cross covariance of lag1 on ETI-BDS','CC113: Cross covariance of lag1 on ETI-SKE','CC114: Cross covariance of lag1 on ETI-WCI','CC115: Cross covariance of lag1 on ETI-SLI','CC116: Cross covariance of lag1 on SKE-PID','CC117: Cross covariance of lag1 on SKE-BDS','CC118: Cross covariance of lag1 on SKE-ETI','CC119: Cross covariance of lag1 on SKE-WCI','CC120: Cross covariance of lag1 on SKE-SLI','CC121: Cross covariance of lag1 on WCI-PID','CC122: Cross covariance of lag1 on WCI-BDS','CC123: Cross covariance of lag1 on WCI-ETI','CC124: Cross covariance of lag1 on WCI-SKE','CC125: Cross covariance of lag1 on WCI-SLI','CC126: Cross covariance of lag1 on SLI-PID','CC127: Cross covariance of lag1 on SLI-BDS','CC128: Cross covariance of lag1 on SLI-ETI','CC129: Cross covariance of lag1 on SLI-SKE','CC130: Cross covariance of lag1 on SLI-WCI','CC201: Cross covariance of lag2 on PID-BDS','CC202: Cross covariance of lag2 on PID-ETI','CC203: Cross covariance of lag2 on PID-SKE','CC204: Cross covariance of lag2 on PID-WCI','CC205: Cross covariance of lag2 on PID-SLI','CC206: Cross covariance of lag2 on BDS-PID','CC207: Cross covariance of lag2 on BDS-ETI','CC208: Cross covariance of lag2 on BDS-SKE','CC209: Cross covariance of lag2 on BDS-WCI','CC210: Cross covariance of lag2 on BDS-SLI','CC211: Cross covariance of lag2 on ETI-PID','CC212: Cross covariance of lag2 on ETI-BDS','CC213: Cross covariance of lag2 on ETI-SKE','CC214: Cross covariance of lag2 on ETI-WCI','CC215: Cross covariance of lag2 on ETI-SLI','CC216: Cross covariance of lag2 on SKE-PID','CC217: Cross covariance of lag2 on SKE-BDS','CC218: Cross covariance of lag2 on SKE-ETI','CC219: Cross covariance of lag2 on SKE-WCI','CC220: Cross covariance of lag2 on SKE-SLI','CC221: Cross covariance of lag2 on WCI-PID','CC222: Cross covariance of lag2 on WCI-BDS','CC223: Cross covariance of lag2 on WCI-ETI','CC224: Cross covariance of lag2 on WCI-SKE','CC225: Cross covariance of lag2 on WCI-SLI','CC226: Cross covariance of lag2 on SLI-PID','CC227: Cross covariance of lag2 on SLI-BDS','CC228: Cross covariance of lag2 on SLI-ETI','CC229: Cross covariance of lag2 on SLI-SKE','CC230: Cross covariance of lag2 on SLI-WCI']
    df_cod = DataFrame(data=vec_all, index=seq_id,columns=colname)
    return df_cod
