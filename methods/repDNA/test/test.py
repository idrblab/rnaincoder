__author__ = 'Fule Liu'
# -*- coding: UTF-8 -*-
import time
import sys
sys.path.append("../")

if __name__ == '__main__':
    """This a test script for repDNA_manual.
    """
    start_time = time.time()
    error = False

    # ###############################################################################################
    # Basic function.

    # Read sequence data from FASTA files
    from util import get_data

    if get_data(open('example.fasta')) != ['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC']:
        print("Error, the basic function get_data1")
        error = True
    if get_data(['GACTGAACTGCACTTTGGTTTCATATTATTTGctc']) != ['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC']:
        print("Error, the basic function get_data2")
        error = True

    seqs = get_data(open('example.fasta'), desc=True)
    seq = seqs[0]
    if seq.name != 'misc_ppid_3700':
        print("Error, the basic function get_data3")
        error = True
    if seq.seq != 'GACTGAACTGCACTTTGGTTTCATATTATTTGCTC':
        print("Error, the basic function get_data4")
        error = True
    if seq.no != 1:
        print("Error, the basic function get_data5")
        error = True

    # Normalization of physicochemical index.
    from util import normalize_index

    phyche_index = [
        [0.026, 0.036, 0.031, 0.033, 0.016, 0.026, 0.014, 0.031, 0.025, 0.025, 0.026, 0.036, 0.017, 0.025, 0.016,
         0.026],
        [0.038, 0.038, 0.037, 0.036, 0.025, 0.042, 0.026, 0.037, 0.038, 0.036, 0.042, 0.038, 0.018, 0.038, 0.025,
         0.038]]
    if normalize_index(phyche_index) \
            != [[0.06, 1.5, 0.78, 1.07, -1.38, 0.06, -1.66, 0.78, -0.08, -0.08, 0.06, 1.5, -1.23, -0.08, -1.38, 0.06],
                [0.5, 0.5, 0.36, 0.22, -1.36, 1.08, -1.22, 0.36, 0.5, 0.22, 1.08, 0.5, -2.37, 0.5, -1.36, 0.5]]:
        print("Error, the basic function normalize_index")
        error = True

    print("Basic function test end!")

    # ######################################################################################
    # Nucleic acid Composition

    # Basic kmer
    from nac import Kmer

    kmer = Kmer(k=2)
    if kmer.make_kmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC']) \
            != [[1, 3, 0, 3, 2, 0, 0, 4, 2, 2, 1, 1, 2, 2, 4, 7]]:
        print("Error, Nucleic acid composition basic kmer1")
        error = True
    kmer = Kmer(k=2, upto=True)
    if kmer.make_kmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC']) \
            != [[7, 7, 6, 15, 1, 3, 0, 3, 2, 0, 0, 4, 2, 2, 1, 1, 2, 2, 4, 7]]:
        print("Error, Nucleic acid composition basic kmer2")
        error = True
    kmer = Kmer(k=2, normalize=True)
    if kmer.make_kmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC']) \
            != [[0.029, 0.088, 0.0, 0.088, 0.059, 0.0, 0.0, 0.118, 0.059, 0.059, 0.029, 0.029, 0.059, 0.059, 0.118,
                 0.206]]:
        print("Error, Nucleic acid composition basic kmer3")
        error = True
    kmer = Kmer(k=2, normalize=True, upto=True)
    if kmer.make_kmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC']) \
            != [[0.2, 0.2, 0.171, 0.429, 0.029, 0.088, 0.0, 0.088, 0.059, 0.0, 0.0, 0.118, 0.059, 0.059, 0.029, 0.029,
                 0.059, 0.059, 0.118, 0.206]]:
        print("Error, Nucleic acid composition basic kmer4")
        error = True

    print("Basic kmer test end!")

    # RevcKmer
    from nac import RevcKmer

    revckmer = RevcKmer(k=2)
    if revckmer.make_revckmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC']) != [[8, 4, 4, 3, 6, 1, 0, 4, 2, 2]]:
        print("Error, Nucleic acid composition revckmer1")
        error = True
    revckmer = RevcKmer(k=2, upto=True)
    if revckmer.make_revckmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC']) != [[22, 13, 8, 4, 4, 3, 6, 1, 0, 4, 2, 2]]:
        print("Error, Nucleic acid composition revckmer2")
        error = True
    revckmer = RevcKmer(k=2, normalize=True)
    if revckmer.make_revckmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC']) \
            != [[0.235, 0.118, 0.118, 0.088, 0.176, 0.029, 0.0, 0.118, 0.059, 0.059]]:
        print("Error, Nucleic acid composition revckmer3")
        error = True
    revckmer = RevcKmer(k=2, normalize=True, upto=True)
    if revckmer.make_revckmer_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC']) \
            != [[0.629, 0.371, 0.235, 0.118, 0.118, 0.088, 0.176, 0.029, 0.0, 0.118, 0.059, 0.059]]:
        print("Error, Nucleic acid composition revckmer4")
        error = True

    print("RevcKmer test end!")

    # Increment of Diversity
    from nac import IDkmer

    idkmer = IDkmer()
    if idkmer.make_idkmer_vec(open('example.fasta'), open('pos.fasta'), open('neg.fasta')) \
            != [[4.54, 29.19, -103.245, -92.99, -144.589, -145.351, -154.0, -154.0, -153.58, -153.58, -147.207,
                 -147.207]]:
        print("Error, Nucleic acid composition IDkmer1")
        error = True
    idkmer = IDkmer(k=2)
    if idkmer.make_idkmer_vec(open('example.fasta'), open('pos.fasta'), open('neg.fasta')) \
            != [[4.54, 29.19, -103.245, -92.99]]:
        print("Error, Nucleic acid composition IDkmer2")
        error = True
    idkmer = IDkmer(k=2, upto=False)
    if idkmer.make_idkmer_vec(open('example.fasta'), open('pos.fasta'), open('neg.fasta')) != [[16.288, 62.27]]:
        print("Error, Nucleic acid composition IDkmer3")
        error = True

    print("IDkmer test end!")

    # ###############################################################################
    # Autocorrelation

    # Dinucleotide-based Auto covariance
    from ac import DAC

    dac = DAC(2)
    if dac.make_dac_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt']) \
            != [[-0.175, -0.185, -0.173, -0.004]]:
        print("Error, Autocorrelation DAC1")
        error = True
    vec = dac.make_dac_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    if len(vec[0]) != 76:
        print("Error, Autocorrelation DAC2")
        error = True

    from util import normalize_index

    phyche_index = [
        [2.26, 3.03, 2.03, 3.83, 1.78, 1.65, 2.00, 2.03, 1.93, 2.61, 1.65, 3.03, 1.20, 1.93, 1.78, 2.26],
        [7.65, 8.93, 7.08, 9.07, 6.38, 8.04, 6.23, 7.08, 8.56, 9.53, 8.04, 8.93, 6.23, 8.56, 6.38, 7.65]]
    if dac.make_dac_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt'],
                        extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True)) \
            != [[-0.175, -0.185, -0.5, -0.504, -0.173, -0.004, 0.009, 0.106]]:
        print("Error, Autocorrelation DAC3")
        error = True

    print("DAC test end!")

    # Dinucleotide-based Cross covariance
    from ac import DCC

    dcc = DCC(2)
    if dcc.make_dcc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt']) \
            != [[-0.141, -0.238, -0.064, -0.047]]:
        print("Error, Autocorrelation DCC1")
        error = True
    vec = dcc.make_dcc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    if len(vec[0]) != 2812:
        print("Error, Autocorrelation DCC2")
        error = True

    from util import normalize_index

    phyche_index = [
        [2.26, 3.03, 2.03, 3.83, 1.78, 1.65, 2.00, 2.03, 1.93, 2.61, 1.65, 3.03, 1.20, 1.93, 1.78, 2.26],
        [7.65, 8.93, 7.08, 9.07, 6.38, 8.04, 6.23, 7.08, 8.56, 9.53, 8.04, 8.93, 6.23, 8.56, 6.38, 7.65]]
    if dcc.make_dcc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt'],
                        extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True)) \
            != [[-0.141, -0.464, -0.535, -0.238, -0.521, -0.512, -0.18, -0.128, -0.477, -0.097, -0.132, -0.403, -0.064,
                 0.027, 0.042, -0.047, 0.103, 0.142, -0.153, -0.216, -0.061, 0.041, -0.021, 0.162]]:
        print("Error, Autocorrelation DCC3")
        error = True

    print("DCC test end!")

    # Dinucleotide-based Auto-cross covariance
    from ac import DACC

    dacc = DACC(2)
    if dacc.make_dacc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt']) \
            != [[-0.175, -0.185, -0.173, -0.004, -0.141, -0.238, -0.064, -0.047]]:
        print("Error, Autocorrelation DACC1")
        error = True
    vec = dacc.make_dacc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    if len(vec[0]) != 2888:
        print("Error, Autocorrelation DACC2")
        error = True

    from util import normalize_index

    phyche_index = [
        [2.26, 3.03, 2.03, 3.83, 1.78, 1.65, 2.00, 2.03, 1.93, 2.61, 1.65, 3.03, 1.20, 1.93, 1.78, 2.26],
        [7.65, 8.93, 7.08, 9.07, 6.38, 8.04, 6.23, 7.08, 8.56, 9.53, 8.04, 8.93, 6.23, 8.56, 6.38, 7.65]]
    if dacc.make_dacc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt'],
                          extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True)) \
            != [[-0.175, -0.185, -0.5, -0.504, -0.173, -0.004, 0.009, 0.106, -0.141, -0.464, -0.535, -0.238, -0.521,
                 -0.512, -0.18, -0.128, -0.477, -0.097, -0.132, -0.403, -0.064, 0.027, 0.042, -0.047, 0.103, 0.142,
                 -0.153, -0.216, -0.061, 0.041, -0.021, 0.162]]:
        print("Error, Autocorrelation DACC3")
        error = True

    print("DACC test end!")

    # Trinucleotide-based Auto covariance
    from ac import TAC

    tac = TAC(2)
    if tac.make_tac_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome']) \
            != [[-0.332, 0.493, 0.319, 0.012]]:
        print("Error, Autocorrelation TAC1")
        error = True
    vec = tac.make_tac_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    if len(vec[0]) != 24:
        print("Error, Autocorrelation TAC2")
        error = True

    from util import normalize_index

    phyche_index = [
        [7.176, 6.272, 4.736, 7.237, 3.810, 4.156, 4.156, 6.033, 3.410, 3.524, 4.445, 6.033, 1.613, 5.087, 2.169, 7.237,
         3.581, 3.239, 1.668, 2.169, 6.813, 3.868, 5.440, 4.445, 3.810, 4.678, 5.440, 4.156, 2.673, 3.353, 1.668, 4.736,
         4.214, 3.925, 3.353, 5.087, 2.842, 2.448, 4.678, 3.524, 3.581, 2.448, 3.868, 4.156, 3.467, 3.925, 3.239, 6.272,
         2.955, 3.467, 2.673, 1.613, 1.447, 3.581, 3.810, 3.410, 1.447, 2.842, 6.813, 3.810, 2.955, 4.214, 3.581,
         7.176]]
    if tac.make_tac_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome'],
                        extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True)) \
            != [[-0.332, 0.493, 0.296, 0.319, 0.012, -0.308]]:
        print("Error, Autocorrelation TAC3")
        error = True

    print("TAC test end!")

    # Trinucleotide-based Cross covariance
    from ac import TCC

    tcc = TCC(2)
    if tcc.make_tcc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome']) \
            != [[-0.299, -0.11, -0.24, 0.001]]:
        print("Error, Autocorrelation TCC1")
        error = True
    vec = tcc.make_tcc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    if len(vec[0]) != 264:
        print("Error, Autocorrelation TCC2")
        error = True

    from util import normalize_index

    phyche_index = [
        [7.176, 6.272, 4.736, 7.237, 3.810, 4.156, 4.156, 6.033, 3.410, 3.524, 4.445, 6.033, 1.613, 5.087, 2.169, 7.237,
         3.581, 3.239, 1.668, 2.169, 6.813, 3.868, 5.440, 4.445, 3.810, 4.678, 5.440, 4.156, 2.673, 3.353, 1.668, 4.736,
         4.214, 3.925, 3.353, 5.087, 2.842, 2.448, 4.678, 3.524, 3.581, 2.448, 3.868, 4.156, 3.467, 3.925, 3.239, 6.272,
         2.955, 3.467, 2.673, 1.613, 1.447, 3.581, 3.810, 3.410, 1.447, 2.842, 6.813, 3.810, 2.955, 4.214, 3.581,
         7.176]]
    if tcc.make_tcc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome'],
                        extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True)) \
            != [[-0.299, 0.356, -0.11, -0.155, -0.561, -0.307, -0.24, 0.199, 0.001, 0.003, 0.002, 0.257]]:
        print("Error, Autocorrelation TCC3")
        error = True

    print("TCC test end!")

    # Trinucleotide-based Auto-cross covariance
    from ac import TACC

    tacc = TACC(2)
    if tacc.make_tacc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome']) \
            != [[-0.332, 0.493, 0.319, 0.012, -0.299, -0.11, -0.24, 0.001]]:
        print("Error, Autocorrelation TACC1")
        error = True
    vec = tacc.make_tacc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    if len(vec[0]) != 288:
        print("Error, Autocorrelation TACC2")
        error = True

    from util import normalize_index

    phyche_index = [
        [7.176, 6.272, 4.736, 7.237, 3.810, 4.156, 4.156, 6.033, 3.410, 3.524, 4.445, 6.033, 1.613, 5.087, 2.169, 7.237,
         3.581, 3.239, 1.668, 2.169, 6.813, 3.868, 5.440, 4.445, 3.810, 4.678, 5.440, 4.156, 2.673, 3.353, 1.668, 4.736,
         4.214, 3.925, 3.353, 5.087, 2.842, 2.448, 4.678, 3.524, 3.581, 2.448, 3.868, 4.156, 3.467, 3.925, 3.239, 6.272,
         2.955, 3.467, 2.673, 1.613, 1.447, 3.581, 3.810, 3.410, 1.447, 2.842, 6.813, 3.810, 2.955, 4.214, 3.581,
         7.176]]
    if tacc.make_tacc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome'],
                          extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True)) \
            != [[-0.332, 0.493, 0.296, 0.319, 0.012, -0.308, -0.299, 0.356, -0.11, -0.155, -0.561, -0.307, -0.24, 0.199,
                 0.001, 0.003, 0.002, 0.257]]:
        print("Error, Autocorrelation TACC3")
        error = True

    print("TACC test end!")

    # ########################################################################################
    # Pseudo Nucleic acid Composition

    # Pseudo dinucleotide composition
    from psenac import PseDNC

    psednc = PseDNC()
    vec = psednc.make_psednc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    if vec != [
        [0.023, 0.069, 0.0, 0.069, 0.046, 0.0, 0.0, 0.092, 0.046, 0.046, 0.023, 0.023, 0.046, 0.046, 0.092, 0.16,
         0.0823, 0.0705, 0.0694]]:
        print("Error, Pseudo Nucleic acid Composition PseDNC1")
        error = True
    if len(vec[0]) != 19:
        print("Error, Pseudo Nucleic acid Composition PseDNC2")
        error = True
    psednc = PseDNC(lamada=2, w=0.1)
    vec = psednc.make_psednc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    if vec != [
        [0.021, 0.063, 0.0, 0.063, 0.042, 0.0, 0.0, 0.084, 0.042, 0.042, 0.021, 0.021, 0.042, 0.042, 0.084, 0.148,
         0.152, 0.1301]]:
        print("Error, Pseudo Nucleic acid Composition PseDNC3")
        error = True
    if len(vec[0]) != 18:
        print("Error, Pseudo Nucleic acid Composition PseDNC4")
        error = True

    from util import normalize_index

    phyche_index = [
        [1.019, -0.918, 0.488, 0.567, 0.567, -0.070, -0.579, 0.488, -0.654, -2.455, -0.070, -0.918, 1.603, -0.654,
         0.567, 1.019]]
    vec = psednc.make_psednc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'],
                                 extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True))
    if vec != [
        [0.021, 0.063, 0.0, 0.063, 0.042, 0.0, 0.0, 0.085, 0.042, 0.042, 0.021, 0.021, 0.042, 0.042, 0.085, 0.148,
         0.1515, 0.1292]]:
        print("Error, Pseudo Nucleic acid Composition PseDNC5")
        error = True
    if len(vec[0]) != 18:
        print("Error, Pseudo Nucleic acid Composition PseDNC6")
        error = True

    print("PseDNC test end!")

    # Pseudo k-tupler composition
    from psenac import PseKNC

    pseknc = PseKNC()
    vec = pseknc.make_pseknc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    if vec != [
        [0.0, 0.014, 0.0, 0.0, 0.0, 0.0, 0.0, 0.041, 0.0, 0.0, 0.0, 0.0, 0.014, 0.0, 0.0, 0.027, 0.0, 0.014, 0.0, 0.014,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.014, 0.027, 0.014, 0.014, 0.014, 0.0, 0.0, 0.014, 0.0, 0.0,
         0.014, 0.0, 0.0, 0.0, 0.014, 0.0, 0.0, 0.0, 0.014, 0.0, 0.0, 0.0, 0.027, 0.014, 0.0, 0.0, 0.0, 0.014, 0.027,
         0.014, 0.0, 0.014, 0.014, 0.027, 0.041, 0.55]]:
        print("Error, Pseudo Nucleic acid Composition PsekNC1")
        error = True
    if len(vec[0]) != 65:
        print("Error, Pseudo Nucleic acid Composition PsekNC2")
        error = True
    pseknc = PseKNC(k=2, lamada=1, w=0.05)
    vec = pseknc.make_pseknc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'])
    if vec != [
        [0.026, 0.079, 0.0, 0.079, 0.052, 0.0, 0.0, 0.105, 0.052, 0.052, 0.026, 0.026, 0.052, 0.052, 0.105, 0.183,
         0.1089]]:
        print("Error, Pseudo Nucleic acid Composition PsekNC3")
        error = True
    if len(vec[0]) != 17:
        print("Error, Pseudo Nucleic acid Composition PsekNC4")
        error = True

    phyche_index = [[1.019, -0.918, 0.488, 0.567, 0.567, -0.070, -0.579, 0.488, -0.654, -2.455, -0.070, -0.918, 1.603,
                     -0.654, 0.567, 1.019]]
    from util import normalize_index

    vec = pseknc.make_pseknc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'],
                                 extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True))
    if vec != [
        [0.026, 0.079, 0.0, 0.079, 0.053, 0.0, 0.0, 0.105, 0.053, 0.053, 0.026, 0.026, 0.053, 0.053, 0.105, 0.184,
         0.1066]]:
        print("Error, Pseudo Nucleic acid Composition PsekNC5")
        error = True
    if len(vec[0]) != 17:
        print("Error, Pseudo Nucleic acid Composition PsekNC6")
        error = True

    print("PseKNC test end!")

    # Parallel correlation pseudo dinucleotide composition
    from psenac import PCPseDNC

    pc_psednc = PCPseDNC()
    vec = pc_psednc.make_pcpsednc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt'])
    if vec != [[0.027, 0.08, 0.0, 0.08, 0.053, 0.0, 0.0, 0.106, 0.053, 0.053, 0.027, 0.027, 0.053, 0.053, 0.106, 0.186,
                0.0948]]:
        print("Error, Pseudo Nucleic acid Composition PCPseDNC1")
        error = True
    if len(vec[0]) != 17:
        print("Error, Pseudo Nucleic acid Composition PCPseDNC2")
        error = True
    pc_psednc = PCPseDNC(lamada=2, w=0.05)
    vec = pc_psednc.make_pcpsednc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    if vec != [[0.025, 0.075, 0.0, 0.075, 0.05, 0.0, 0.0, 0.1, 0.05, 0.05, 0.025, 0.025, 0.05, 0.05, 0.1, 0.175, 0.072,
                0.0757]]:
        print("Error, Pseudo Nucleic acid Composition PCPseDNC3")
        error = True
    if len(vec[0]) != 18:
        print("Error, Pseudo Nucleic acid Composition PCPseDNC4")
        error = True

    from util import normalize_index

    phyche_index = [
        [1.019, -0.918, 0.488, 0.567, 0.567, -0.070, -0.579, 0.488, -0.654, -2.455, -0.070, -0.918, 1.603, -0.654,
         0.567, 1.019]]
    vec = pc_psednc.make_pcpsednc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt'],
                                      extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True))
    if vec != [
        [0.025, 0.074, 0.0, 0.074, 0.049, 0.0, 0.0, 0.098, 0.049, 0.049, 0.025, 0.025, 0.049, 0.049, 0.098, 0.172,
         0.0869, 0.0771]]:
        print("Error, Pseudo Nucleic acid Composition PCPseDNC5")
        error = True
    if len(vec[0]) != 18:
        print("Error, Pseudo Nucleic acid Composition PCPseDNC6")
        error = True

    print("PC-PseDNC test end!")

    # Parallel correlation pseudo trinucleotide composition
    from psenac import PCPseTNC

    pc_psetnc = PCPseTNC()
    vec = pc_psetnc.make_pcpsetnc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'],
                                      phyche_index=['Dnase I', 'Nucleosome'])
    if vec != [
        [0.0, 0.027, 0.0, 0.0, 0.0, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0, 0.0, 0.027, 0.0, 0.0, 0.053, 0.0, 0.027, 0.0, 0.027,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.027, 0.053, 0.027, 0.027, 0.027, 0.0, 0.0, 0.027, 0.0, 0.0,
         0.027, 0.0, 0.0, 0.0, 0.027, 0.0, 0.0, 0.0, 0.027, 0.0, 0.0, 0.0, 0.053, 0.027, 0.0, 0.0, 0.0, 0.027, 0.053,
         0.027, 0.0, 0.027, 0.027, 0.053, 0.08, 0.1229]]:
        print("Error, Pseudo Nucleic acid Composition PCPseTNC1")
        error = True
    if len(vec[0]) != 65:
        print("Error, Pseudo Nucleic acid Composition PCPseTNC2")
        error = True
    pc_psetnc = PCPseTNC(lamada=2, w=0.05)
    vec = pc_psetnc.make_pcpsetnc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    if vec != [
        [0.0, 0.024, 0.0, 0.0, 0.0, 0.0, 0.0, 0.073, 0.0, 0.0, 0.0, 0.0, 0.024, 0.0, 0.0, 0.048, 0.0, 0.024, 0.0, 0.024,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.024, 0.048, 0.024, 0.024, 0.024, 0.0, 0.0, 0.024, 0.0, 0.0,
         0.024, 0.0, 0.0, 0.0, 0.024, 0.0, 0.0, 0.0, 0.024, 0.0, 0.0, 0.0, 0.048, 0.024, 0.0, 0.0, 0.0, 0.024, 0.048,
         0.024, 0.0, 0.024, 0.024, 0.048, 0.073, 0.0851, 0.1147]]:
        print("Error, Pseudo Nucleic acid Composition PCPseTNC3")
        error = True
    if len(vec[0]) != 66:
        print("Error, Pseudo Nucleic acid Composition PCPseTNC4")
        error = True

    from util import normalize_index

    phyche_index = [
        [7.176, 6.272, 4.736, 7.237, 3.810, 4.156, 4.156, 6.033, 3.410, 3.524, 4.445, 6.033, 1.613, 5.087, 2.169, 7.237,
         3.581, 3.239, 1.668, 2.169, 6.813, 3.868, 5.440, 4.445, 3.810, 4.678, 5.440, 4.156, 2.673, 3.353, 1.668, 4.736,
         4.214, 3.925, 3.353, 5.087, 2.842, 2.448, 4.678, 3.524, 3.581, 2.448, 3.868, 4.156, 3.467, 3.925, 3.239, 6.272,
         2.955, 3.467, 2.673, 1.613, 1.447, 3.581, 3.810, 3.410, 1.447, 2.842, 6.813, 3.810, 2.955, 4.214, 3.581, 7.176]
    ]
    vec = pc_psetnc.make_pcpsetnc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome'],
                                      extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True))
    if vec != [
        [0.0, 0.023, 0.0, 0.0, 0.0, 0.0, 0.0, 0.07, 0.0, 0.0, 0.0, 0.0, 0.023, 0.0, 0.0, 0.046, 0.0, 0.023, 0.0, 0.023,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.023, 0.046, 0.023, 0.023, 0.023, 0.0, 0.0, 0.023, 0.0, 0.0,
         0.023, 0.0, 0.0, 0.0, 0.023, 0.0, 0.0, 0.0, 0.023, 0.0, 0.0, 0.0, 0.046, 0.023, 0.0, 0.0, 0.0, 0.023, 0.046,
         0.023, 0.0, 0.023, 0.023, 0.046, 0.07, 0.1102, 0.1229]]:
        print("Error, Pseudo Nucleic acid Composition PCPseTNC5")
        error = True
    if len(vec[0]) != 66:
        print("Error, Pseudo Nucleic acid Composition PCPseTNC6")
        error = True

    print("PC-PseTNC test end!")

    # Series correlation pseudo dinucleotide composition
    from psenac import SCPseDNC

    sc_psednc = SCPseDNC()
    vec = sc_psednc.make_scpsednc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt'])
    if vec != [[0.03, 0.09, 0.0, 0.09, 0.06, 0.0, 0.0, 0.12, 0.06, 0.06, 0.03, 0.03, 0.06, 0.06, 0.12, 0.21, -0.0088,
                -0.0093]]:
        print("Error, Pseudo Nucleic acid Composition SCPseDNC1")
        error = True
    if len(vec[0]) != 18:
        print("Error, Pseudo Nucleic acid Composition SCPseDNC2")
        error = True
    sc_psednc = SCPseDNC(lamada=2, w=0.05)
    vec = sc_psednc.make_scpsednc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    if len(vec[0]) != 92:
        print("Error, Pseudo Nucleic acid Composition SCPseDNC3")
        error = True

    from util import normalize_index

    phyche_index = [
        [1.019, -0.918, 0.488, 0.567, 0.567, -0.070, -0.579, 0.488, -0.654, -2.455, -0.070, -0.918, 1.603, -0.654,
         0.567, 1.019]]
    vec = sc_psednc.make_scpsednc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Twist', 'Tilt'],
                                      extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True))
    if vec != [
        [0.03, 0.09, 0.0, 0.09, 0.06, 0.0, 0.0, 0.12, 0.06, 0.06, 0.03, 0.03, 0.06, 0.06, 0.12, 0.209, -0.0088, -0.0093,
         0.0004, -0.0088, 0.0, 0.01]]:
        print("Error, Pseudo Nucleic acid Composition SCPseDNC4")
        error = True
    if len(vec[0]) != 22:
        print("Error, Pseudo Nucleic acid Composition SCPseDNC5")
        error = True

    print("SC-PseDNC test end!")

    # Series correlation pseudo trinucleotide composition
    from psenac import SCPseTNC

    sc_psetnc = SCPseTNC()
    vec = sc_psetnc.make_scpsetnc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome'])
    if len(vec[0]) != 66:
        print("Error, Pseudo Nucleic acid Composition SCPseTNC1")
        error = True
    sc_psetnc = SCPseTNC(lamada=2, w=0.05)
    vec = sc_psetnc.make_scpsetnc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], all_property=True)
    if len(vec[0]) != 88:
        print("Error, Pseudo Nucleic acid Composition SCPseTNC2")
        error = True

    from util import normalize_index

    phyche_index = [
        [7.176, 6.272, 4.736, 7.237, 3.810, 4.156, 4.156, 6.033, 3.410, 3.524, 4.445, 6.033, 1.613, 5.087, 2.169, 7.237,
         3.581, 3.239, 1.668, 2.169, 6.813, 3.868, 5.440, 4.445, 3.810, 4.678, 5.440, 4.156, 2.673, 3.353, 1.668, 4.736,
         4.214, 3.925, 3.353, 5.087, 2.842, 2.448, 4.678, 3.524, 3.581, 2.448, 3.868, 4.156, 3.467, 3.925, 3.239, 6.272,
         2.955, 3.467, 2.673, 1.613, 1.447, 3.581, 3.810, 3.410, 1.447, 2.842, 6.813, 3.810, 2.955, 4.214, 3.581,
         7.176]]
    vec = sc_psetnc.make_scpsetnc_vec(['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC'], phyche_index=['Dnase I', 'Nucleosome'],
                                      extra_phyche_index=normalize_index(phyche_index, is_convert_dict=True))
    if len(vec[0]) != 70:
        print("Error, Pseudo Nucleic acid Composition SCPseTNC3")
        error = True

    print("SC-PseTNC test end!")

    if error is False:
        print("Congratulation! It works. :D")
    else:
        print("Although the results are not expected due to machine precision, it also works. :D")

    print("Test end!")
    print("Used time: %ss" % (time.time() - start_time))