#!/usr/bin/env python
import sys
import argparse
import numpy as np

from collections import defaultdict
from itertools import product
from pandas import DataFrame
import pandas as pd

# from my_tqdm import my_tqdm
# from fasta_reader import Reader
from tqdm import tnrange, trange, tqdm, tqdm_notebook



list_new_featrures = ['SHDCA: Strong H-Bond donors composition of A','SHDAB: Strong H-Bond donors transition between A and B','SHDA0: Strong H-Bond donors distribution of 0.00A','SHDA1: Strong H-Bond donors distribution of 0.25A','SHDA2: Strong H-Bond donors distribution of 0.50A','SHDA3: Strong H-Bond donors distribution of 0.75A','SHDA4: Strong H-Bond donors distribution of 1.00A','MLFCA: Linear free energy composition of A','MLFCB: Linear free energy composition of B','MLFAB: Linear free energy transition between A and B','MLFA0: Linear free energy distribution of 0.00A','MLFA1: Linear free energy distribution of 0.25A','MLFA2: Linear free energy distribution of 0.50A','MLFA3: Linear free energy distribution of 0.75A','MLFA4: Linear free energy distribution of 1.00A','MLFB0: Linear free energy distribution of 0.00B','MLFB1: Linear free energy distribution of 0.25B','MLFB2: Linear free energy distribution of 0.50B','MLFB3: Linear free energy distribution of 0.75B','MLFB4: Linear free energy distribution of 1.00B','SHAAB: Strong H-Bond acceptors_3 transition between A and B','SHAAC: Strong H-Bond acceptors_3 transition between A and C','PHBAB: Potential Hydrogen Bonds_3 transition between A and B','PHBAC: Potential Hydrogen Bonds_3 transition between A and C','SHD: Strong H-Bond donors0','MLF: Linear free energy0','MLF: Linear free energy1','SHD: Strong H-Bond donors00','SHD: Strong H-Bond donors01','SHD: Strong H-Bond donors10','MLF: Linear free energy00','MLF: Linear free energy01','MLF: Linear free energy10','MLF: Linear free energy11','SHA: Strong H-Bond acceptors_301','SHA: Strong H-Bond acceptors_302','SHA: Strong H-Bond acceptors_310','SHA: Strong H-Bond acceptors_320','PHB: Potential Hydrogen Bonds_301','PHB: Potential Hydrogen Bonds_302','PHB: Potential Hydrogen Bonds_310','PHB: Potential Hydrogen Bonds_320','SHD: Strong H-Bond donors000','SHD: Strong H-Bond donors001','SHD: Strong H-Bond donors010','SHD: Strong H-Bond donors011','SHD: Strong H-Bond donors100','SHD: Strong H-Bond donors101','SHD: Strong H-Bond donors110','MLF: Linear free energy000','MLF: Linear free energy001','MLF: Linear free energy010','MLF: Linear free energy011','MLF: Linear free energy100','MLF: Linear free energy101','MLF: Linear free energy110','MLF: Linear free energy111','SHA: Strong H-Bond acceptors_3001','SHA: Strong H-Bond acceptors_3002','SHA: Strong H-Bond acceptors_3010','SHA: Strong H-Bond acceptors_3011','SHA: Strong H-Bond acceptors_3012','SHA: Strong H-Bond acceptors_3020','SHA: Strong H-Bond acceptors_3021','SHA: Strong H-Bond acceptors_3022','SHA: Strong H-Bond acceptors_3100','SHA: Strong H-Bond acceptors_3101','SHA: Strong H-Bond acceptors_3102','SHA: Strong H-Bond acceptors_3110','SHA: Strong H-Bond acceptors_3120','SHA: Strong H-Bond acceptors_3200','SHA: Strong H-Bond acceptors_3201','SHA: Strong H-Bond acceptors_3202','SHA: Strong H-Bond acceptors_3210','SHA: Strong H-Bond acceptors_3220','PHB: Potential Hydrogen Bonds_3001','PHB: Potential Hydrogen Bonds_3002','PHB: Potential Hydrogen Bonds_3010','PHB: Potential Hydrogen Bonds_3011','PHB: Potential Hydrogen Bonds_3012','PHB: Potential Hydrogen Bonds_3020','PHB: Potential Hydrogen Bonds_3021','PHB: Potential Hydrogen Bonds_3022','PHB: Potential Hydrogen Bonds_3100','PHB: Potential Hydrogen Bonds_3101','PHB: Potential Hydrogen Bonds_3102','PHB: Potential Hydrogen Bonds_3110','PHB: Potential Hydrogen Bonds_3120','PHB: Potential Hydrogen Bonds_3200','PHB: Potential Hydrogen Bonds_3201','PHB: Potential Hydrogen Bonds_3202','PHB: Potential Hydrogen Bonds_3210','PHB: Potential Hydrogen Bonds_3220','CNHCA: NH- count composition of A','CNHCB: NH- count composition of B','CNHAB: NH- count transition between A and B','CNHA0: NH- count distribution of 0.00A','CNHA1: NH- count distribution of 0.25A','CNHA2: NH- count distribution of 0.50A','CNHA3: NH- count distribution of 0.75A','CNHA4: NH- count distribution of 1.00A','CNHB0: NH- count distribution of 0.00B','CNHB1: NH- count distribution of 0.25B','CNHB2: NH- count distribution of 0.50B','CNHB3: NH- count distribution of 0.75B','CNHB4: NH- count distribution of 1.00B','CNH: NH- count0','CNH: NH- count1','CNH: NH- count00','CNH: NH- count01','CNH: NH- count10','CNH: NH- count11','CNH: NH- count000','CNH: NH- count001','CNH: NH- count010','CNH: NH- count011','CNH: NH- count100','CNH: NH- count101','CNH: NH- count110','CNH: NH- count111','MReCA: Molar refractivity composition of A','MReCB: Molar refractivity composition of B','MReAB: Molar refractivity transition between A and B','MReA0: Molar refractivity distribution of 0.00A','MReA1: Molar refractivity distribution of 0.25A','MReA2: Molar refractivity distribution of 0.50A','MReA3: Molar refractivity distribution of 0.75A','MReA4: Molar refractivity distribution of 1.00A','MReB0: Molar refractivity distribution of 0.00B','MReB1: Molar refractivity distribution of 0.25B','MReB2: Molar refractivity distribution of 0.50B','MReB3: Molar refractivity distribution of 0.75B','MReB4: Molar refractivity distribution of 1.00B','MRe: Molar refractivity0','MRe: Molar refractivity1','MRe: Molar refractivity00','MRe: Molar refractivity01','MRe: Molar refractivity10','MRe: Molar refractivity11','MRe: Molar refractivity000','MRe: Molar refractivity001','MRe: Molar refractivity010','MRe: Molar refractivity011','MRe: Molar refractivity100','MRe: Molar refractivity101','MRe: Molar refractivity110','MRe: Molar refractivity111','PSNCA: Primary or secondary nitrogens composition of A','PSNAB: Primary or secondary nitrogens transition between A and B','PSNA0: Primary or secondary nitrogens distribution of 0.00A','PSNA1: Primary or secondary nitrogens distribution of 0.25A','PSNA2: Primary or secondary nitrogens distribution of 0.50A','PSNA3: Primary or secondary nitrogens distribution of 0.75A','PSNA4: Primary or secondary nitrogens distribution of 1.00A','SLFAB: Sum of path lengths starting from oxygens_3 transition between A and B','SLFAC: Sum of path lengths starting from oxygens_3 transition between A and C','TPSAB: Topological polar surface area_3 transition between A and B','TPSAC: Topological polar surface area_3 transition between A and C','PSN: Primary or secondary nitrogens0','PSN: Primary or secondary nitrogens00','PSN: Primary or secondary nitrogens01','PSN: Primary or secondary nitrogens10','SLF: Sum of path lengths starting from oxygens_301','SLF: Sum of path lengths starting from oxygens_302','SLF: Sum of path lengths starting from oxygens_310','SLF: Sum of path lengths starting from oxygens_320','TPS: Topological polar surface area_301','TPS: Topological polar surface area_302','TPS: Topological polar surface area_310','TPS: Topological polar surface area_320','PSN: Primary or secondary nitrogens000','PSN: Primary or secondary nitrogens001','PSN: Primary or secondary nitrogens010','PSN: Primary or secondary nitrogens011','PSN: Primary or secondary nitrogens100','PSN: Primary or secondary nitrogens101','PSN: Primary or secondary nitrogens110','SLF: Sum of path lengths starting from oxygens_3001','SLF: Sum of path lengths starting from oxygens_3002','SLF: Sum of path lengths starting from oxygens_3010','SLF: Sum of path lengths starting from oxygens_3011','SLF: Sum of path lengths starting from oxygens_3012','SLF: Sum of path lengths starting from oxygens_3020','SLF: Sum of path lengths starting from oxygens_3021','SLF: Sum of path lengths starting from oxygens_3022','SLF: Sum of path lengths starting from oxygens_3100','SLF: Sum of path lengths starting from oxygens_3101','SLF: Sum of path lengths starting from oxygens_3102','SLF: Sum of path lengths starting from oxygens_3110','SLF: Sum of path lengths starting from oxygens_3120','SLF: Sum of path lengths starting from oxygens_3200','SLF: Sum of path lengths starting from oxygens_3201','SLF: Sum of path lengths starting from oxygens_3202','SLF: Sum of path lengths starting from oxygens_3210','SLF: Sum of path lengths starting from oxygens_3220','TPS: Topological polar surface area_3001','TPS: Topological polar surface area_3002','TPS: Topological polar surface area_3010','TPS: Topological polar surface area_3011','TPS: Topological polar surface area_3012','TPS: Topological polar surface area_3020','TPS: Topological polar surface area_3021','TPS: Topological polar surface area_3022','TPS: Topological polar surface area_3100','TPS: Topological polar surface area_3101','TPS: Topological polar surface area_3102','TPS: Topological polar surface area_3110','TPS: Topological polar surface area_3120','TPS: Topological polar surface area_3200','TPS: Topological polar surface area_3201','TPS: Topological polar surface area_3202','TPS: Topological polar surface area_3210','TPS: Topological polar surface area_3220','HPCCA: Gas-hexadecane PC composition of A','HPCAB: Gas-hexadecane PC transition between A and B','HPCA0: Gas-hexadecane PC distribution of 0.00A','HPCA1: Gas-hexadecane PC distribution of 0.25A','HPCA2: Gas-hexadecane PC distribution of 0.50A','HPCA3: Gas-hexadecane PC distribution of 0.75A','HPCA4: Gas-hexadecane PC distribution of 1.00A','HPCAB: Gas-hexadecane PC_3 transition between A and B','HPCAC: Gas-hexadecane PC_3 transition between A and C','HPC: Gas-hexadecane PC0','HPC: Gas-hexadecane PC00','HPC: Gas-hexadecane PC01','HPC: Gas-hexadecane PC10','HPC: Gas-hexadecane PC_301','HPC: Gas-hexadecane PC_302','HPC: Gas-hexadecane PC_310','HPC: Gas-hexadecane PC_320','HPC: Gas-hexadecane PC000','HPC: Gas-hexadecane PC001','HPC: Gas-hexadecane PC010','HPC: Gas-hexadecane PC011','HPC: Gas-hexadecane PC100','HPC: Gas-hexadecane PC101','HPC: Gas-hexadecane PC110','HPC: Gas-hexadecane PC_3001','HPC: Gas-hexadecane PC_3002','HPC: Gas-hexadecane PC_3010','HPC: Gas-hexadecane PC_3011','HPC: Gas-hexadecane PC_3012','HPC: Gas-hexadecane PC_3020','HPC: Gas-hexadecane PC_3021','HPC: Gas-hexadecane PC_3022','HPC: Gas-hexadecane PC_3100','HPC: Gas-hexadecane PC_3101','HPC: Gas-hexadecane PC_3102','HPC: Gas-hexadecane PC_3110','HPC: Gas-hexadecane PC_3120','HPC: Gas-hexadecane PC_3200','HPC: Gas-hexadecane PC_3201','HPC: Gas-hexadecane PC_3202','HPC: Gas-hexadecane PC_3210','HPC: Gas-hexadecane PC_3220','LFICA: Lipoaffinity index composition of A','LFIAB: Lipoaffinity index transition between A and B','LFIA0: Lipoaffinity index distribution of 0.00A','LFIA1: Lipoaffinity index distribution of 0.25A','LFIA2: Lipoaffinity index distribution of 0.50A','LFIA3: Lipoaffinity index distribution of 0.75A','LFIA4: Lipoaffinity index distribution of 1.00A','LFIAB: Lipoaffinity index_3 transition between A and B','LFIAC: Lipoaffinity index_3 transition between A and C','LFI: Lipoaffinity index0','LFI: Lipoaffinity index00','LFI: Lipoaffinity index01','LFI: Lipoaffinity index10','LFI: Lipoaffinity index_301','LFI: Lipoaffinity index_302','LFI: Lipoaffinity index_310','LFI: Lipoaffinity index_320','LFI: Lipoaffinity index000','LFI: Lipoaffinity index001','LFI: Lipoaffinity index010','LFI: Lipoaffinity index011','LFI: Lipoaffinity index100','LFI: Lipoaffinity index101','LFI: Lipoaffinity index110','LFI: Lipoaffinity index_3001','LFI: Lipoaffinity index_3002','LFI: Lipoaffinity index_3010','LFI: Lipoaffinity index_3011','LFI: Lipoaffinity index_3012','LFI: Lipoaffinity index_3020','LFI: Lipoaffinity index_3021','LFI: Lipoaffinity index_3022','LFI: Lipoaffinity index_3100','LFI: Lipoaffinity index_3101','LFI: Lipoaffinity index_3102','LFI: Lipoaffinity index_3110','LFI: Lipoaffinity index_3120','LFI: Lipoaffinity index_3200','LFI: Lipoaffinity index_3201','LFI: Lipoaffinity index_3202','LFI: Lipoaffinity index_3210','LFI: Lipoaffinity index_3220']

#############################################################################
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 11:16:50 2016

@author: jessime

A couple of functions for autodetecting if code is being run in the notebook.
The appropriate tqdm progressbar will be returned.
"""


def _is_kernel():
    if 'IPython' not in sys.modules:
        # IPython hasn't been imported, definitely not
        return False
    from IPython import get_ipython
    # check for `kernel` attribute on the IPython instance
    return getattr(get_ipython(), 'kernel', None) is not None


def my_tqdm():
    return tqdm_notebook if _is_kernel() else tqdm


def my_trange():
    return tnrange if _is_kernel() else trange


######################################################################
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 13:10:49 2016

@author: jessime
"""


class Reader():
    """Fixes any compatibility issues a fasta file might have with this code.

    Parameters
    ----------
    infasta : str (default=None)
        Name of input fasta file to be manipulated
    outfasta : str (default=None)
        location to store extracted data from infasta
    names : iter (default=None)
        Common style names to use in header lines

    Attributes
    ----------
    data : list
        Raw lines of the infasta file
        Note: This is different than the data attribute in other classes

    Examples
    --------
    Putting the sequence on one line instead of breaking it every 80 chars.
    Making sure the whole sequence is capitalized.
    Restructuring the name line to work with GENCODE's naming.
    """

    def __init__(self, infasta=None, outfasta=None, names=None):
        self.infasta = infasta
        self.outfasta = outfasta
        self.names = names

        self.data = None

    def _read_data(self):
        """Sets data to stripped lines from the fasta file
        """
        with open(self.infasta) as infasta:
            self.data = [l.strip() for l in infasta.readlines()]

    def _upper_seq_per_line(self):
        """Sets data to upper case, single line sequences for each header
        """
        new_data = []
        seq = ''
        for i, line in enumerate(self.data):
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
        self.data = new_data

    def get_lines(self):
        self._read_data()
        self._upper_seq_per_line()
        return self.data

    def get_seqs(self):
        clean_data = self.get_lines()
        seqs = clean_data[1::2]
        return seqs

    def get_headers(self):
        clean_data = self.get_lines()
        headers = clean_data[::2]
        return headers

    def get_data(self, tuples_only=False):
        clean_data = self.get_lines()
        headers = clean_data[::2]
        seqs = clean_data[1::2]
        tuples = zip(headers, seqs)
        if tuples_only:
            return tuples
        else:
            return tuples, headers, seqs

    def supply_basic_header(self):
        """Convert headerlines to GENCODE format with only common name and length"""
        new_fasta = []

        if self.names is None:
            self.names = iter(self.get_headers())
        for i, line in enumerate(self.data):
            if line[0] == '>':
                name = next(self.names).strip('>')
                length = len(self.data[i + 1])
                new_fasta.append('>||||{}||{}|'.format(name, length))
            else:
                new_fasta.append(line)
        return new_fasta

    def save(self):
        """Write self.data to a new fasta file"""
        with open(self.outfasta, 'w') as outfasta:
            for line in self.data:
                outfasta.write(line + '\n')


##############################################################################
# !/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description
-----------
Generates a kmer count matrix of m rows by n columns,
where m is the number of transcripts in a fasta file and n is 4^kmer.

Examples
--------
The default settings produce a binary, normalized numpy file:
    $ python kmer_counts.py /path/to/rnas.fa -o /path/to/out.npy

To get a human readable csv file, set the nonbinary flag:
    $ python kmer_counts.py /path/to/rnas.fa -o /path/to/out.csv -nb

If you want to add default labels, also set the label flag:
    $ python kmer_counts.py /path/to/rnas.fa -o /path/to/out.csv -nb -lb

You can change also change the size of the kmer you're using, and prevent normalization:
    $ python kmer_counts.py /path/to/rnas.fa -o /path/to/out.npy -k 4 -nc -ns

Notes
-----
For more sophisticated options, you cannot use the command-line, but need python instead.
To label the axes of the matrix, for example, you can call BasicCounter('/path/rnas.fa').to_csv(names)

Issues
------
Any issues can be reported to https://github.com/CalabreseLab #TODO

---
"""
import Bio.SeqIO as Seq


class BasicCounter:
    """Generates overlapping kmer counts for a fasta file

    Parameters
    ----------
    infasta : str (default=None)
        Full path to fasta file to be counted
    outfile : str (default=None)
        Full path to the counts file to be saved
    k : int (default=6)
        Size of kmer to be counted
    binary : bool (default=True)
        Saves as numpy array if True, else saves as csv
    mean : bool, np.array, str (default=True)
        Set the mean to 0 for each kmer/column of the count matrix.
        If str, provide path to a previously calculated mean array.
    std : bool or str (default=True)
        Set the std. dev. to 1 for each kmer/column of the count matrix
        If str, provide path to a previously calculated std array.
    leave : bool (default=True)
        Set to False if get_counts is used within another tqdm loop
    silent : bool (default=False)
        Set to True to turn off tqdm progress bar

    Attributes
    ----------
    counts : None
        Stores the ndarray of kmer counts
    kmers : list
        str elements of all kmers of size k
    map : dict
        Mapping of kmers to column values
    """

    def __init__(self, infasta=None, k=6,
                 binary=True, mean=False, std=False,
                 leave=True, silent=False, label=False):
        self.infasta = infasta
        self.seqs = None
        if infasta is not None:
            self.seqs = Reader(infasta).get_seqs()
        # self.outfile = outfile
        self.k = k
        self.binary = binary
        self.mean = mean
        if isinstance(mean, str):
            self.mean = np.load(mean)
        self.std = std
        if isinstance(std, str):
            self.std = np.load(std)
        self.leave = leave
        self.silent = silent
        self.label = label

        self.counts = None
        self.kmers = [''.join(i) for i in product('AGTC', repeat=k)]
        self.map = {k: i for k, i in zip(self.kmers, range(4 ** k))}

        if len(self.seqs) == 1 and self.std is True:
            err = ('You cannot standardize a single sequence. '
                   'Please pass the path to an std. dev. array, '
                   'or use raw counts by setting std=False.')
            raise ValueError(err)

    # def occurrences(self, row, seq):
    #     """Counts kmers on a per kilobase scale"""
    #     counts = defaultdict(int)
    #     length = len(seq)
    #     # increment = 1000 / length
    #     increment = 1

    #     for c in range(length - self.k + 1):
    #         kmer = seq[c:c + self.k]
    #         counts[kmer] += increment
    #     sk = length-self.k + 1
    #     ak = 1/(4 ** (4-self.k))
    #     for kmer, n in counts.items():
    #         if kmer in self.map:
    #             # print('count is {}'.format(n))
    #             # print('sk is {}'.format(sk))
    #             # print('ak is {}'.format(ak))
    #             NUM = n * ak / sk # pmlipred paper used

    #             # print('NUM is {}'.format(NUM))
    #             row[self.map[kmer]] = NUM
    #     return row
    # old version
    def occurrences(self, row, seq):
        """Counts kmers on a per kilobase scale"""
        counts = defaultdict(int)
        length = len(seq)
        increment = 1000 / length
        for c in range(length - self.k + 1):
            kmer = seq[c:c + self.k]
            counts[kmer] += increment
        for kmer, n in counts.items():
            if kmer in self.map:
                row[self.map[kmer]] = n
        return row

    def _progress(self):
        """Determine which iterator to loop over for counting."""
        if self.silent:
            return self.seqs

        if not self.leave:
            tqdm_seqs = my_tqdm()(self.seqs, desc='Kmers', leave=False)
        else:
            tqdm_seqs = my_tqdm()(self.seqs)

        return tqdm_seqs

    def center(self):
        """mean center counts by column"""
        if self.mean is True:
            self.mean = np.mean(self.counts, axis=0)
        self.counts -= self.mean

    def standardize(self):
        """divide out the standard deviations from columns of the count matrix"""
        if self.std is True:
            self.std = np.std(self.counts, axis=0)
        self.counts /= self.std

    # def get_counts(self):
    #     """Generates kmer counts for a fasta file"""
    #     self.counts = np.zeros([len(self.seqs), 4 ** self.k], dtype=np.float32)
    #     seqs = self._progress()
    #     for i, seq in enumerate(seqs):
    #         self.counts[i] = self.occurrences(self.counts[i], seq)
    #     if self.mean is not False:
    #         self.center()
    #     if self.std is not False:
    #         self.standardize()

    #     seqname = []
    #     for seq in Seq.parse(self.infasta, 'fasta'):
    #         seqid = seq.id
    #         seqname.append(seqid)
    #     if self.k == 3:
    #         columnanmes = ['KMAAA: Transcript k-mer AAA content','KMAAG: Transcript k-mer AAG content','KMAAT: Transcript k-mer AAT content','KMAAC: Transcript k-mer AAC content','KMAGA: Transcript k-mer AGA content','KMAGG: Transcript k-mer AGG content','KMAGT: Transcript k-mer AGT content','KMAGC: Transcript k-mer AGC content','KMATA: Transcript k-mer ATA content','KMATG: Transcript k-mer ATG content','KMATT: Transcript k-mer ATT content','KMATC: Transcript k-mer ATC content','KMACA: Transcript k-mer ACA content','KMACG: Transcript k-mer ACG content','KMACT: Transcript k-mer ACT content','KMACC: Transcript k-mer ACC content','KMGAA: Transcript k-mer GAA content','KMGAG: Transcript k-mer GAG content','KMGAT: Transcript k-mer GAT content','KMGAC: Transcript k-mer GAC content','KMGGA: Transcript k-mer GGA content','KMGGG: Transcript k-mer GGG content','KMGGT: Transcript k-mer GGT content','KMGGC: Transcript k-mer GGC content','KMGTA: Transcript k-mer GTA content','KMGTG: Transcript k-mer GTG content','KMGTT: Transcript k-mer GTT content','KMGTC: Transcript k-mer GTC content','KMGCA: Transcript k-mer GCA content','KMGCG: Transcript k-mer GCG content','KMGCT: Transcript k-mer GCT content','KMGCC: Transcript k-mer GCC content','KMTAA: Transcript k-mer TAA content','KMTAG: Transcript k-mer TAG content','KMTAT: Transcript k-mer TAT content','KMTAC: Transcript k-mer TAC content','KMTGA: Transcript k-mer TGA content','KMTGG: Transcript k-mer TGG content','KMTGT: Transcript k-mer TGT content','KMTGC: Transcript k-mer TGC content','KMTTA: Transcript k-mer TTA content','KMTTG: Transcript k-mer TTG content','KMTTT: Transcript k-mer TTT content','KMTTC: Transcript k-mer TTC content','KMTCA: Transcript k-mer TCA content','KMTCG: Transcript k-mer TCG content','KMTCT: Transcript k-mer TCT content','KMTCC: Transcript k-mer TCC content','KMCAA: Transcript k-mer CAA content','KMCAG: Transcript k-mer CAG content','KMCAT: Transcript k-mer CAT content','KMCAC: Transcript k-mer CAC content','KMCGA: Transcript k-mer CGA content','KMCGG: Transcript k-mer CGG content','KMCGT: Transcript k-mer CGT content','KMCGC: Transcript k-mer CGC content','KMCTA: Transcript k-mer CTA content','KMCTG: Transcript k-mer CTG content','KMCTT: Transcript k-mer CTT content','KMCTC: Transcript k-mer CTC content','KMCCA: Transcript k-mer CCA content','KMCCG: Transcript k-mer CCG content','KMCCT: Transcript k-mer CCT content','KMCCC: Transcript k-mer CCC content']
    #     else:
    #         columnanmes = self.kmers
    #     df = DataFrame(data=self.counts, index=seqname, columns=columnanmes,dtype='double')
    #     # df.to_csv(self.outfile)
    #     return df

    def get_counts(self):
        """Generates kmer counts for a fasta file"""
        self.counts = np.zeros([len(self.seqs), 4 ** self.k], dtype=np.float32)
        seqs = self._progress()
        for i, seq in enumerate(seqs):
            self.counts[i] = self.occurrences(self.counts[i], seq)
        if self.mean is not False:
            self.center()
        if self.std is not False:
            self.standardize()

        seqname = []
        for seq in Seq.parse(self.infasta, 'fasta'):
            seqid = seq.id
            seqname.append(seqid)
        if self.k == 3:
            columnanmes = ['KMAAA: Transcript k-mer AAA content','KMAAG: Transcript k-mer AAG content','KMAAT: Transcript k-mer AAT content','KMAAC: Transcript k-mer AAC content','KMAGA: Transcript k-mer AGA content','KMAGG: Transcript k-mer AGG content','KMAGT: Transcript k-mer AGT content','KMAGC: Transcript k-mer AGC content','KMATA: Transcript k-mer ATA content','KMATG: Transcript k-mer ATG content','KMATT: Transcript k-mer ATT content','KMATC: Transcript k-mer ATC content','KMACA: Transcript k-mer ACA content','KMACG: Transcript k-mer ACG content','KMACT: Transcript k-mer ACT content','KMACC: Transcript k-mer ACC content','KMGAA: Transcript k-mer GAA content','KMGAG: Transcript k-mer GAG content','KMGAT: Transcript k-mer GAT content','KMGAC: Transcript k-mer GAC content','KMGGA: Transcript k-mer GGA content','KMGGG: Transcript k-mer GGG content','KMGGT: Transcript k-mer GGT content','KMGGC: Transcript k-mer GGC content','KMGTA: Transcript k-mer GTA content','KMGTG: Transcript k-mer GTG content','KMGTT: Transcript k-mer GTT content','KMGTC: Transcript k-mer GTC content','KMGCA: Transcript k-mer GCA content','KMGCG: Transcript k-mer GCG content','KMGCT: Transcript k-mer GCT content','KMGCC: Transcript k-mer GCC content','KMTAA: Transcript k-mer TAA content','KMTAG: Transcript k-mer TAG content','KMTAT: Transcript k-mer TAT content','KMTAC: Transcript k-mer TAC content','KMTGA: Transcript k-mer TGA content','KMTGG: Transcript k-mer TGG content','KMTGT: Transcript k-mer TGT content','KMTGC: Transcript k-mer TGC content','KMTTA: Transcript k-mer TTA content','KMTTG: Transcript k-mer TTG content','KMTTT: Transcript k-mer TTT content','KMTTC: Transcript k-mer TTC content','KMTCA: Transcript k-mer TCA content','KMTCG: Transcript k-mer TCG content','KMTCT: Transcript k-mer TCT content','KMTCC: Transcript k-mer TCC content','KMCAA: Transcript k-mer CAA content','KMCAG: Transcript k-mer CAG content','KMCAT: Transcript k-mer CAT content','KMCAC: Transcript k-mer CAC content','KMCGA: Transcript k-mer CGA content','KMCGG: Transcript k-mer CGG content','KMCGT: Transcript k-mer CGT content','KMCGC: Transcript k-mer CGC content','KMCTA: Transcript k-mer CTA content','KMCTG: Transcript k-mer CTG content','KMCTT: Transcript k-mer CTT content','KMCTC: Transcript k-mer CTC content','KMCCA: Transcript k-mer CCA content','KMCCG: Transcript k-mer CCG content','KMCCT: Transcript k-mer CCT content','KMCCC: Transcript k-mer CCC content']
        else:
            columnanmes = self.kmers
        df = DataFrame(data=self.counts, index=seqname, columns=columnanmes,dtype='double')
        # df.to_csv(self.outfile)
        return df

class BasicCounter_2:
    """Generates overlapping kmer counts for a fasta file

    Parameters
    ----------
    infasta : str (default=None)
        Full path to fasta file to be counted
    outfile : str (default=None)
        Full path to the counts file to be saved
    k : int (default=6)
        Size of kmer to be counted
    binary : bool (default=True)
        Saves as numpy array if True, else saves as csv
    mean : bool, np.array, str (default=True)
        Set the mean to 0 for each kmer/column of the count matrix.
        If str, provide path to a previously calculated mean array.
    std : bool or str (default=True)
        Set the std. dev. to 1 for each kmer/column of the count matrix
        If str, provide path to a previously calculated std array.
    leave : bool (default=True)
        Set to False if get_counts is used within another tqdm loop
    silent : bool (default=False)
        Set to True to turn off tqdm progress bar

    Attributes
    ----------
    counts : None
        Stores the ndarray of kmer counts
    kmers : list
        str elements of all kmers of size k
    map : dict
        Mapping of kmers to column values
    """

    def __init__(self, infasta=None, k=6,
                 binary=True, mean=False, std=False,
                 leave=True, silent=False, label=False, encode = None):
        self.infasta = infasta
        self.seqs = None
        self.encode = encode
        if infasta is not None:
            self.seqs = Reader(infasta).get_seqs()

        
        # self.outfile = outfile
        self.k = k
        self.binary = binary
        self.mean = mean
        if isinstance(mean, str):
            self.mean = np.load(mean)
        self.std = std
        if isinstance(std, str):
            self.std = np.load(std)
        self.leave = leave
        self.silent = silent
        self.label = label

        self.counts = None
        self.kmers = [''.join(i) for i in product('01', repeat=k)]
        self.map = {k: i for k, i in zip(self.kmers, range(2 ** k))}

        if len(self.seqs) == 1 and self.std is True:
            err = ('You cannot standardize a single sequence. '
                   'Please pass the path to an std. dev. array, '
                   'or use raw counts by setting std=False.')
            raise ValueError(err)
    # old version
    def occurrences(self, row, seq):
        """Counts kmers on a per kilobase scale"""
        counts = defaultdict(int)
        length = len(seq)
        increment = 1000 / length
        for c in range(length - self.k + 1):
            kmer = seq[c:c + self.k]
            counts[kmer] += increment
        for kmer, n in counts.items():
            if kmer in self.map:
                row[self.map[kmer]] = n
        frequency = list(row)        
        return frequency

    def _progress(self):
        """Determine which iterator to loop over for counting."""
        if self.silent:
            return self.seqs

        if not self.leave:
            tqdm_seqs = my_tqdm()(self.seqs, desc='Kmers', leave=False)
        else:
            tqdm_seqs = my_tqdm()(self.seqs)

        return tqdm_seqs

    def center(self, data):
        """mean center counts by column"""
        if self.mean is True:
            self.mean = np.mean(data, axis=0)
        data -= self.mean
        return data

    def standardize(self):
        """divide out the standard deviations from columns of the count matrix"""
        if self.std is True:
            self.std = np.std(self.counts, axis=0)
        self.counts /= self.std

    def get_counts(self,):
        """Generates kmer counts for a fasta file"""
        self.counts = np.zeros([len(self.seqs), 2 ** self.k], dtype=np.float32)
        seqs = self._progress()
 
        # if self.encode[0].count("1", 0, 4) == 1:
        feature = []
        for i, seq in enumerate(seqs):
            for encode in self.encode:
                seq01 = seq.replace('A', encode[0]).replace('C', encode[1]).replace('G', encode[2]).replace('T', encode[3])
                feature += self.occurrences(self.counts[i], seq01)
        # print('len(self.seqs)')
        # print(len(self.seqs))
        np.set_printoptions(precision=6)
        data_count = np.array(feature).reshape((len(self.seqs), (2 ** self.k) * len(self.encode)))

        if self.mean is not False:
            data_count = self.center(data_count)
        if self.std is not False:
            self.standardize()

        seqname = []
        for seq in Seq.parse(self.infasta, 'fasta'):
            seqid = seq.id
            seqname.append(seqid)
        # if self.k == 3:
        #     columnanmes = ['KMAAA: Transcript k-mer AAA content','KMAAG: Transcript k-mer AAG content','KMAAT: Transcript k-mer AAT content','KMAAC: Transcript k-mer AAC content','KMAGA: Transcript k-mer AGA content','KMAGG: Transcript k-mer AGG content','KMAGT: Transcript k-mer AGT content','KMAGC: Transcript k-mer AGC content','KMATA: Transcript k-mer ATA content','KMATG: Transcript k-mer ATG content','KMATT: Transcript k-mer ATT content','KMATC: Transcript k-mer ATC content','KMACA: Transcript k-mer ACA content','KMACG: Transcript k-mer ACG content','KMACT: Transcript k-mer ACT content','KMACC: Transcript k-mer ACC content','KMGAA: Transcript k-mer GAA content','KMGAG: Transcript k-mer GAG content','KMGAT: Transcript k-mer GAT content','KMGAC: Transcript k-mer GAC content','KMGGA: Transcript k-mer GGA content','KMGGG: Transcript k-mer GGG content','KMGGT: Transcript k-mer GGT content','KMGGC: Transcript k-mer GGC content','KMGTA: Transcript k-mer GTA content','KMGTG: Transcript k-mer GTG content','KMGTT: Transcript k-mer GTT content','KMGTC: Transcript k-mer GTC content','KMGCA: Transcript k-mer GCA content','KMGCG: Transcript k-mer GCG content','KMGCT: Transcript k-mer GCT content','KMGCC: Transcript k-mer GCC content','KMTAA: Transcript k-mer TAA content','KMTAG: Transcript k-mer TAG content','KMTAT: Transcript k-mer TAT content','KMTAC: Transcript k-mer TAC content','KMTGA: Transcript k-mer TGA content','KMTGG: Transcript k-mer TGG content','KMTGT: Transcript k-mer TGT content','KMTGC: Transcript k-mer TGC content','KMTTA: Transcript k-mer TTA content','KMTTG: Transcript k-mer TTG content','KMTTT: Transcript k-mer TTT content','KMTTC: Transcript k-mer TTC content','KMTCA: Transcript k-mer TCA content','KMTCG: Transcript k-mer TCG content','KMTCT: Transcript k-mer TCT content','KMTCC: Transcript k-mer TCC content','KMCAA: Transcript k-mer CAA content','KMCAG: Transcript k-mer CAG content','KMCAT: Transcript k-mer CAT content','KMCAC: Transcript k-mer CAC content','KMCGA: Transcript k-mer CGA content','KMCGG: Transcript k-mer CGG content','KMCGT: Transcript k-mer CGT content','KMCGC: Transcript k-mer CGC content','KMCTA: Transcript k-mer CTA content','KMCTG: Transcript k-mer CTG content','KMCTT: Transcript k-mer CTT content','KMCTC: Transcript k-mer CTC content','KMCCA: Transcript k-mer CCA content','KMCCG: Transcript k-mer CCG content','KMCCT: Transcript k-mer CCT content','KMCCC: Transcript k-mer CCC content']
        # else:
        dictname = {
            '0010': 'Strong H-Bond donors',
            '0110': 'Linear free energy',
            '0101': 'Molar refractivity',
            '1000': 'Lipoaffinity index',
            '0100': 'Gas-hexadecane PC', 
            '0011': 'NH- count',
            '0001': 'Primary or secondary nitrogens'}
        dictname_sh = {'0010': 'SHD', '0110': 'MLF',
                    '0101': 'MRe', 
                    '1000': 'LFI', '0100': 'HPC',
                    '0011': 'CNH',
                '0001': 'PSN'}
        # columnanmes = [prelabel + conlum  for conlum in self.kmers]
        columnanmes = [dictname_sh[encode] + ": " + dictname[encode] + conlum for encode in self.encode for conlum in self.kmers]

        df = DataFrame(data=data_count, index=seqname, columns=columnanmes,dtype='double')
        # df.to_csv(self.outfile)
        # df01 = df.iloc[:, df.columns.str.contains('0')]
        new_feaname = [i for i in columnanmes if i in  list_new_featrures]
        df_new = df[new_feaname]
        # df = round(df.iloc[:, :], 6)
        return df_new
        # return df
        # elif  self.encode[0].count("1", 0, 4) == 2:
        #     feature = []
        #     for i, seq in enumerate(seqs):
        #         for encode in self.encode:
        #             seq01 = seq.replace('A', encode[0]).replace('C', encode[1]).replace('G', encode[2]).replace('T', encode[3])
        #             feature += self.occurrences(self.counts[i], seq01)
        #     # print('len(self.seqs)')
        #     # print(len(self.seqs))
        #     np.set_printoptions(precision=6)
        #     data_count = np.array(feature).reshape((len(self.seqs), (2 ** self.k) * len(self.encode)))

        #     if self.mean is not False:
        #         data_count = self.center(data_count)
        #     if self.std is not False:
        #         self.standardize()

        #     seqname = []
        #     for seq in Seq.parse(self.infasta, 'fasta'):
        #         seqid = seq.id
        #         seqname.append(seqid)
        #     # if self.k == 3:
        #     #     columnanmes = ['KMAAA: Transcript k-mer AAA content','KMAAG: Transcript k-mer AAG content','KMAAT: Transcript k-mer AAT content','KMAAC: Transcript k-mer AAC content','KMAGA: Transcript k-mer AGA content','KMAGG: Transcript k-mer AGG content','KMAGT: Transcript k-mer AGT content','KMAGC: Transcript k-mer AGC content','KMATA: Transcript k-mer ATA content','KMATG: Transcript k-mer ATG content','KMATT: Transcript k-mer ATT content','KMATC: Transcript k-mer ATC content','KMACA: Transcript k-mer ACA content','KMACG: Transcript k-mer ACG content','KMACT: Transcript k-mer ACT content','KMACC: Transcript k-mer ACC content','KMGAA: Transcript k-mer GAA content','KMGAG: Transcript k-mer GAG content','KMGAT: Transcript k-mer GAT content','KMGAC: Transcript k-mer GAC content','KMGGA: Transcript k-mer GGA content','KMGGG: Transcript k-mer GGG content','KMGGT: Transcript k-mer GGT content','KMGGC: Transcript k-mer GGC content','KMGTA: Transcript k-mer GTA content','KMGTG: Transcript k-mer GTG content','KMGTT: Transcript k-mer GTT content','KMGTC: Transcript k-mer GTC content','KMGCA: Transcript k-mer GCA content','KMGCG: Transcript k-mer GCG content','KMGCT: Transcript k-mer GCT content','KMGCC: Transcript k-mer GCC content','KMTAA: Transcript k-mer TAA content','KMTAG: Transcript k-mer TAG content','KMTAT: Transcript k-mer TAT content','KMTAC: Transcript k-mer TAC content','KMTGA: Transcript k-mer TGA content','KMTGG: Transcript k-mer TGG content','KMTGT: Transcript k-mer TGT content','KMTGC: Transcript k-mer TGC content','KMTTA: Transcript k-mer TTA content','KMTTG: Transcript k-mer TTG content','KMTTT: Transcript k-mer TTT content','KMTTC: Transcript k-mer TTC content','KMTCA: Transcript k-mer TCA content','KMTCG: Transcript k-mer TCG content','KMTCT: Transcript k-mer TCT content','KMTCC: Transcript k-mer TCC content','KMCAA: Transcript k-mer CAA content','KMCAG: Transcript k-mer CAG content','KMCAT: Transcript k-mer CAT content','KMCAC: Transcript k-mer CAC content','KMCGA: Transcript k-mer CGA content','KMCGG: Transcript k-mer CGG content','KMCGT: Transcript k-mer CGT content','KMCGC: Transcript k-mer CGC content','KMCTA: Transcript k-mer CTA content','KMCTG: Transcript k-mer CTG content','KMCTT: Transcript k-mer CTT content','KMCTC: Transcript k-mer CTC content','KMCCA: Transcript k-mer CCA content','KMCCG: Transcript k-mer CCG content','KMCCT: Transcript k-mer CCT content','KMCCC: Transcript k-mer CCC content']
        #     # else:
        #     dictname = {
        #         '0010': 'Strong H-Bond donors',
        #         '0110': 'Linear free energy',
        #         '0101': 'Molar refractivity',
        #         '1000': 'Lipoaffinity index',
        #         '0100': 'Gas-hexadecane PC', 
        #         '0011': 'NH- count',
        #         '0001': 'Primary or secondary nitrogens'}
        #     dictname_sh = {'0010': 'SHD', '0110': 'MLF',
        #                 '0101': 'MRe', 
        #                 '1000': 'LFI', '0100': 'HPC',
        #                 '0011': 'CNH',
        #             '0001': 'PSN'}
        #     # columnanmes = [prelabel + conlum  for conlum in self.kmers]
        #     columnanmes = [dictname_sh[encode] + ": " + dictname[encode] + conlum for encode in self.encode for conlum in self.kmers]

        #     df = DataFrame(data=data_count, index=seqname, columns=columnanmes,dtype='double')

        #     return df
class BasicCounter_3:
    """Generates overlapping kmer counts for a fasta file

    Parameters
    ----------
    infasta : str (default=None)
        Full path to fasta file to be counted
    outfile : str (default=None)
        Full path to the counts file to be saved
    k : int (default=6)
        Size of kmer to be counted
    binary : bool (default=True)
        Saves as numpy array if True, else saves as csv
    mean : bool, np.array, str (default=True)
        Set the mean to 0 for each kmer/column of the count matrix.
        If str, provide path to a previously calculated mean array.
    std : bool or str (default=True)
        Set the std. dev. to 1 for each kmer/column of the count matrix
        If str, provide path to a previously calculated std array.
    leave : bool (default=True)
        Set to False if get_counts is used within another tqdm loop
    silent : bool (default=False)
        Set to True to turn off tqdm progress bar

    Attributes
    ----------
    counts : None
        Stores the ndarray of kmer counts
    kmers : list
        str elements of all kmers of size k
    map : dict
        Mapping of kmers to column values
    """

    def __init__(self, infasta=None, k=6,
                 binary=True, mean=False, std=False,
                 leave=True, silent=False, label=False, encode = None):
        self.infasta = infasta
        self.seqs = None
        self.encodes = encode
        if infasta is not None:
            self.seqs = Reader(infasta).get_seqs()

        
        # self.outfile = outfile
        self.k = k
        self.binary = binary
        self.mean = mean
        if isinstance(mean, str):
            self.mean = np.load(mean)
        self.std = std
        if isinstance(std, str):
            self.std = np.load(std)
        self.leave = leave
        self.silent = silent
        self.label = label

        self.counts = None
        self.kmers = [''.join(i) for i in product('012', repeat=k)]
        self.map = {k: i for k, i in zip(self.kmers, range(3 ** k))}

        if len(self.seqs) == 1 and self.std is True:
            err = ('You cannot standardize a single sequence. '
                   'Please pass the path to an std. dev. array, '
                   'or use raw counts by setting std=False.')
            raise ValueError(err)
    # old version
    def occurrences(self, row, seq):
        """Counts kmers on a per kilobase scale"""
        counts = defaultdict(int)
        length = len(seq)
        increment = 1000 / length
        for c in range(length - self.k + 1):
            kmer = seq[c:c + self.k]
            counts[kmer] += increment
        for kmer, n in counts.items():
            # if kmer in self.map and '0' in kmer:
            row[self.map[kmer]] = n
        frequency = list(row)        
        return frequency

    def _progress(self):
        """Determine which iterator to loop over for counting."""
        if self.silent:
            return self.seqs

        if not self.leave:
            tqdm_seqs = my_tqdm()(self.seqs, desc='Kmers', leave=False)
        else:
            tqdm_seqs = my_tqdm()(self.seqs)

        return tqdm_seqs

    def center(self, data):
        """mean center counts by column"""
        if self.mean is True:
            self.mean = np.mean(data, axis=0)
        data -= self.mean
        return data

    def standardize(self):
        """divide out the standard deviations from columns of the count matrix"""
        if self.std is True:
            self.std = np.std(self.counts, axis=0)
        self.counts /= self.std

    def get_counts(self):
        """Generates kmer counts for a fasta file"""
        self.counts = np.zeros([len(self.seqs), 3 ** self.k], dtype=np.float32)
        seqs = self._progress()

        feature = []
        for i, seq in enumerate(seqs):
            for encode in self.encodes:
                seq01 = seq.replace('A', encode[0]).replace('C', encode[1]).replace('G', encode[2]).replace('T', encode[3])
                # print(seq01)
                feature += self.occurrences(self.counts[i], seq01)
        
        np.set_printoptions(precision=6)
        data_count = np.array(feature).reshape((len(self.seqs), (3 ** self.k) * len(self.encodes)))

        if self.mean is not False:
            data_count = self.center(data_count)
        if self.std is not False:
            self.standardize()

        seqname = []
        for seq in Seq.parse(self.infasta, 'fasta'):
            seqid = seq.id
            seqname.append(seqid)

        dictname = {'1020': 'Lipoaffinity index_3','0102': 'Gas-hexadecane PC_3', '1200':'Strong H-Bond acceptors_3','0120':'Potential Hydrogen Bonds_3','1002':'Sum of path lengths starting from oxygens_3','0021':'Topological polar surface area_3'}

        dictname_sh = {'1020': 'LFI','0102': 'HPC','1200':'SHA','0120':'PHB','1002':'SLF','0021':'TPS'}

        # if self.k == 3:
        #     columnanmes = ['KMAAA: Transcript k-mer AAA content','KMAAG: Transcript k-mer AAG content','KMAAT: Transcript k-mer AAT content','KMAAC: Transcript k-mer AAC content','KMAGA: Transcript k-mer AGA content','KMAGG: Transcript k-mer AGG content','KMAGT: Transcript k-mer AGT content','KMAGC: Transcript k-mer AGC content','KMATA: Transcript k-mer ATA content','KMATG: Transcript k-mer ATG content','KMATT: Transcript k-mer ATT content','KMATC: Transcript k-mer ATC content','KMACA: Transcript k-mer ACA content','KMACG: Transcript k-mer ACG content','KMACT: Transcript k-mer ACT content','KMACC: Transcript k-mer ACC content','KMGAA: Transcript k-mer GAA content','KMGAG: Transcript k-mer GAG content','KMGAT: Transcript k-mer GAT content','KMGAC: Transcript k-mer GAC content','KMGGA: Transcript k-mer GGA content','KMGGG: Transcript k-mer GGG content','KMGGT: Transcript k-mer GGT content','KMGGC: Transcript k-mer GGC content','KMGTA: Transcript k-mer GTA content','KMGTG: Transcript k-mer GTG content','KMGTT: Transcript k-mer GTT content','KMGTC: Transcript k-mer GTC content','KMGCA: Transcript k-mer GCA content','KMGCG: Transcript k-mer GCG content','KMGCT: Transcript k-mer GCT content','KMGCC: Transcript k-mer GCC content','KMTAA: Transcript k-mer TAA content','KMTAG: Transcript k-mer TAG content','KMTAT: Transcript k-mer TAT content','KMTAC: Transcript k-mer TAC content','KMTGA: Transcript k-mer TGA content','KMTGG: Transcript k-mer TGG content','KMTGT: Transcript k-mer TGT content','KMTGC: Transcript k-mer TGC content','KMTTA: Transcript k-mer TTA content','KMTTG: Transcript k-mer TTG content','KMTTT: Transcript k-mer TTT content','KMTTC: Transcript k-mer TTC content','KMTCA: Transcript k-mer TCA content','KMTCG: Transcript k-mer TCG content','KMTCT: Transcript k-mer TCT content','KMTCC: Transcript k-mer TCC content','KMCAA: Transcript k-mer CAA content','KMCAG: Transcript k-mer CAG content','KMCAT: Transcript k-mer CAT content','KMCAC: Transcript k-mer CAC content','KMCGA: Transcript k-mer CGA content','KMCGG: Transcript k-mer CGG content','KMCGT: Transcript k-mer CGT content','KMCGC: Transcript k-mer CGC content','KMCTA: Transcript k-mer CTA content','KMCTG: Transcript k-mer CTG content','KMCTT: Transcript k-mer CTT content','KMCTC: Transcript k-mer CTC content','KMCCA: Transcript k-mer CCA content','KMCCG: Transcript k-mer CCG content','KMCCT: Transcript k-mer CCT content','KMCCC: Transcript k-mer CCC content']
        # else:
        columnanmes = [dictname_sh[encode]+ ": " + dictname[encode] + conlum for encode in self.encodes for conlum in self.kmers]
        df = DataFrame(data=data_count, index=seqname, columns=columnanmes,dtype='double')
        # df01 = df.iloc[:, df.columns.str.contains('0')]
        # if self.encodes == ['0120']:
        #     df.to_csv('/public/home/wangyx/LncRNA/resubmission_code_16methods/Case3_rpi_data/02_EncodingFeatures/RPI1807/test.csv')
        # df.to_csv(self.outfile)
        new_feaname = [i for i in columnanmes if i in  list_new_featrures]
        df_new = df[new_feaname]
        # df = round(df.iloc[:, :], 6)
        return df_new
        # return df