import Bio.SeqIO as Seq
import numpy as np
from pandas import DataFrame
#from One_Hot_Encoder import One_Hot_Encoder


class Onehot:
    def __init__(self, infasta=None):
        self.infasta = infasta
        self.alphabet = "ACGT"

    def get_max(self,infasta):
        seq_len_a = []
        for seq in Seq.parse(infasta,'fasta'):
            sequence = seq.seq
            seq_len = len(sequence)
            seq_len_a.append(seq_len)
        return max(seq_len_a)

    def get_onehot(self, N):   

        seqname = []
        data = []
        max_len = N 

        for seq in Seq.parse(self.infasta,'fasta'):
            seqid = seq.id
            seqname.append(seqid)
            sequence = seq.seq
            one_hot_matrix = One_Hot_Encoder(self.alphabet).encode(sequence)
            seq_len = len(sequence)
            if  seq_len < max_len:
                equal_ele = np.array([float(0)] * 4)
                embed_len = max_len - seq_len
                tmp = np.tile(equal_ele, (embed_len, 1))
                one_hot_matrix = np.concatenate((one_hot_matrix, tmp), axis=0)
            else:
                one_hot_matrix = one_hot_matrix[:max_len]
            data.append(one_hot_matrix)
        
        data01 = np.array(data)
        return seqname,data01


#######################################################################
class One_Hot_Encoder:
    """
    The One_Hot_Encoder class provides functions to encode a string over a
    given alphabet into an integer matrix of shape (len(string), len(alphabet))
    where each row represents a position in the string and each column
    represents a character from the alphabet. Each row has exactly one 1 at the
    matching alphabet character and consists of 0s otherwise.
    """

    def __init__(self, alphabet):
        """ Initialize the object with an alphabet.

        Parameters
        ----------
        alphabet : str
            The alphabet that will be used for encoding/decoding (e.g. "ACGT").
        """
        self.alphabet = alphabet
        self.table = {symbol: i for i, symbol in enumerate(alphabet)}
        self.table_rev = {v: k for k, v in self.table.items()}

    def encode(self, sequence):
        """ Encode a sequence into a one-hot integer matrix.

        The sequence should only contain characters from the alphabet provided to __init__.

        Parameters
        ----------
        sequence : str
            The sequence that should be encoded.

        Returns
        -------
        one_hot: numpy.ndarray
            A numpy array with shape (len(sequence), len(alphabet)).
        """
        one_hot = np.zeros((len(sequence), len(self.table)), np.uint8)
        for i, x in enumerate(sequence):
            if x in self.table.keys():
                one_hot[i, self.table[x]] = 1

        # one_hot[np.arange(len(sequence)), [self.table[x] for x in sequence]] = 1
        return one_hot

    def decode(self, one_hot):
        """ Decode a one-hot integer matrix into the original sequence.

        Parameters
        ----------
        one_hot : numpy.ndarray
            A one-hot matrix (e.g. as created by the encode function).

        Returns
        -------
        sequence: str
            The sequence that is represented by the one-hot matrix.
        """
        return ''.join(self.table_rev[x] for x in np.argmax(one_hot, axis=1))



