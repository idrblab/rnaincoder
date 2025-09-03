#!/usr/bin/env python
# coding: utf-8

# In[1]:


import string, math, copy
import pickle as pkl
import numpy as np
from Bio import SeqIO
from pandas import DataFrame

# In[2]:


AALetter = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

_Hydrophobicity = {'1': 'RKEDQN', '2': 'GASTPHY', '3': 'CLVIMFW'}
# '1'stand for Polar; '2'stand for Neutral, '3' stand for Hydrophobicity

_NormalizedVDWV = {'1': 'GASTPD', '2': 'NVEQIL', '3': 'MHKFRYW'}
# '1'stand for (0-2.78); '2'stand for (2.95-4.0), '3' stand for (4.03-8.08)

_Polarity = {'1': 'LIFWCMVY', '2': 'CPNVEQIL', '3': 'KMHFRYW'}
# '1'stand for (4.9-6.2); '2'stand for (8.0-9.2), '3' stand for (10.4-13.0)

_Charge = {'1': 'KR', '2': 'ANCQGHILMFPSTWYV', '3': 'DE'}
# '1'stand for Positive; '2'stand for Neutral, '3' stand for Negative

_SecondaryStr = {'1': 'EALMQKRH', '2': 'VIYCWFT', '3': 'GNPSD'}
# '1'stand for Helix; '2'stand for Strand, '3' stand for coil

_SolventAccessibility = {'1': 'ALFCGIVW', '2': 'RKQEND', '3': 'MPSTHY'}
# '1'stand for Buried; '2'stand for Exposed, '3' stand for Intermediate

_Polarizability = {'1': 'GASDT', '2': 'CPNVEQIL', '3': 'KMHFRYW'}
# '1'stand for (0-0.108); '2'stand for (0.128-0.186), '3' stand for (0.219-0.409)

_SurfaceTension = {'1': 'GQDNAHR', '2': 'KTSEC', '3': 'ILMFPWYV'}
# In[3]:


_AATProperty = (
_Hydrophobicity, _NormalizedVDWV, _Polarity, _Charge, _SecondaryStr, _SolventAccessibility, _Polarizability,
_SurfaceTension)

_AATPropertyName = (
'_Hydrophobicity', '_NormalizedVDWV', '_Polarity', '_Charge', '_SecondaryStr', '_SolventAccessibility',
'_Polarizability', '_SurfaceTension')


featurename_dict = {
"_Composition of alanine": "_AACoA: Composition of alanine",
"_Composition of cysteine": "_AACoC: Composition of cysteine",
"_Composition of aspartic acid": "_AACoD: Composition of aspartic acid",
"_Composition of glutamic acid": "_AACoE: Composition of glutamic acid",
"_Composition of phenylalanine": "_AACoF: Composition of phenylalanine",
"_Composition of glycine": "_AACoG: Composition of glycine",
"_Composition of histidine": "_AACoH: Composition of histidine",
"_Composition of isoleucine": "_AACoI: Composition of isoleucine",
"_Composition of lysine": "_AACoK: Composition of lysine",
"_Composition of leucine": "_AACoL: Composition of leucine",
"_Composition of methionine": "_AACoM: Composition of methionine",
"_Composition of asparagine": "_AACoN: Composition of asparagine",
"_Composition of proline": "_AACoP: Composition of proline",
"_Composition of glutamine": "_AACoQ: Composition of glutamine",
"_Composition of arginine": "_AACoR: Composition of arginine",
"_Composition of serine": "_AACoS: Composition of serine",
"_Composition of threonine": "_AACoT: Composition of threonine",
"_Composition of valine": "_AACoV: Composition of valine",
"_Composition of tryptophan": "_AACoW: Composition of tryptophan",
"_Composition of tyrosine": "_AACoY: Composition of tyrosine",
"_Electric charge composition of A": "_CagCA: Electric charge composition of A",
"_Electric charge composition of B": "_CagCB: Electric charge composition of B",
"_Electric charge composition of C": "_CagCC: Electric charge composition of C",
"_Electric charge transition between A and B": "_CgTAB: Electric charge transition between A and B",
"_Electric charge transition between A and C": "_CgTAC: Electric charge transition between A and C",
"_Electric charge transition between B and C": "_CgTBC: Electric charge transition between B and C",
"_Electric charge distribution of 0.00A": "_CgDA0: Electric charge distribution of 0.00A",
"_Electric charge distribution of 0.25A": "_CgDA1: Electric charge distribution of 0.25A",
"_Electric charge distribution of 0.50A": "_CgDA2: Electric charge distribution of 0.50A",
"_Electric charge distribution of 0.75A": "_CgDA3: Electric charge distribution of 0.75A",
"_Electric charge distribution of 1.00A": "_CgDA4: Electric charge distribution of 1.00A",
"_Electric charge distribution of 0.00B": "_CgDB0: Electric charge distribution of 0.00B",
"_Electric charge distribution of 0.25B": "_CgDB1: Electric charge distribution of 0.25B",
"_Electric charge distribution of 0.50B": "_CgDB2: Electric charge distribution of 0.50B",
"_Electric charge distribution of 0.75B": "_CgDB3: Electric charge distribution of 0.75B",
"_Electric charge distribution of 1.00B": "_CgDB4: Electric charge distribution of 1.00B",
"_Electric charge distribution of 0.00C": "_CgDC0: Electric charge distribution of 0.00C",
"_Electric charge distribution of 0.25C": "_CgDC1: Electric charge distribution of 0.25C",
"_Electric charge distribution of 0.50C": "_CgDC2: Electric charge distribution of 0.50C",
"_Electric charge distribution of 0.75C": "_CgDC3: Electric charge distribution of 0.75C",
"_Electric charge distribution of 1.00C": "_CgDC4: Electric charge distribution of 1.00C",
"_Hydrophobicity composition of A": "_HydCA: Hydrophobicity composition of A",
"_Hydrophobicity composition of B": "_HydCB: Hydrophobicity composition of B",
"_Hydrophobicity composition of C": "_HydCC: Hydrophobicity composition of C",
"_Hydrophobicity transition between A and B": "_HyTAB: Hydrophobicity transition between A and B",
"_Hydrophobicity transition between A and C": "_HyTAC: Hydrophobicity transition between A and C",
"_Hydrophobicity transition between B and C": "_HyTBC: Hydrophobicity transition between B and C",
"_Hydrophobicity distribution of 0.00A": "_HyDA0: Hydrophobicity distribution of 0.00A",
"_Hydrophobicity distribution of 0.25A": "_HyDA1: Hydrophobicity distribution of 0.25A",
"_Hydrophobicity distribution of 0.50A": "_HyDA2: Hydrophobicity distribution of 0.50A",
"_Hydrophobicity distribution of 0.75A": "_HyDA3: Hydrophobicity distribution of 0.75A",
"_Hydrophobicity distribution of 1.00A": "_HyDA4: Hydrophobicity distribution of 1.00A",
"_Hydrophobicity distribution of 0.00B": "_HyDB0: Hydrophobicity distribution of 0.00B",
"_Hydrophobicity distribution of 0.25B": "_HyDB1: Hydrophobicity distribution of 0.25B",
"_Hydrophobicity distribution of 0.50B": "_HyDB2: Hydrophobicity distribution of 0.50B",
"_Hydrophobicity distribution of 0.75B": "_HyDB3: Hydrophobicity distribution of 0.75B",
"_Hydrophobicity distribution of 1.00B": "_HyDB4: Hydrophobicity distribution of 1.00B",
"_Hydrophobicity distribution of 0.00C": "_HyDC0: Hydrophobicity distribution of 0.00C",
"_Hydrophobicity distribution of 0.25C": "_HyDC1: Hydrophobicity distribution of 0.25C",
"_Hydrophobicity distribution of 0.50C": "_HyDC2: Hydrophobicity distribution of 0.50C",
"_Hydrophobicity distribution of 0.75C": "_HyDC3: Hydrophobicity distribution of 0.75C",
"_Hydrophobicity distribution of 1.00C": "_HyDC4: Hydrophobicity distribution of 1.00C",
"_Polarity composition of A": "_PorCA: Polarity composition of A",
"_Polarity composition of B": "_PorCB: Polarity composition of B",
"_Polarity composition of C": "_PorCC: Polarity composition of C",
"_Polarity transition between A and B": "_PrTAB: Polarity transition between A and B",
"_Polarity transition between A and C": "_PrTAC: Polarity transition between A and C",
"_Polarity transition between B and C": "_PrTBC: Polarity transition between B and C",
"_Polarity distribution of 0.00A": "_PrDA0: Polarity distribution of 0.00A",
"_Polarity distribution of 0.25A": "_PrDA1: Polarity distribution of 0.25A",
"_Polarity distribution of 0.50A": "_PrDA2: Polarity distribution of 0.50A",
"_Polarity distribution of 0.75A": "_PrDA3: Polarity distribution of 0.75A",
"_Polarity distribution of 1.00A": "_PrDA4: Polarity distribution of 1.00A",
"_Polarity distribution of 0.00B": "_PrDB0: Polarity distribution of 0.00B",
"_Polarity distribution of 0.25B": "_PrDB1: Polarity distribution of 0.25B",
"_Polarity distribution of 0.50B": "_PrDB2: Polarity distribution of 0.50B",
"_Polarity distribution of 0.75B": "_PrDB3: Polarity distribution of 0.75B",
"_Polarity distribution of 1.00B": "_PrDB4: Polarity distribution of 1.00B",
"_Polarity distribution of 0.00C": "_PrDC0: Polarity distribution of 0.00C",
"_Polarity distribution of 0.25C": "_PrDC1: Polarity distribution of 0.25C",
"_Polarity distribution of 0.50C": "_PrDC2: Polarity distribution of 0.50C",
"_Polarity distribution of 0.75C": "_PrDC3: Polarity distribution of 0.75C",
"_Polarity distribution of 1.00C": "_PrDC4: Polarity distribution of 1.00C",
"_Polarizability composition of A": "_PozCA: Polarizability composition of A",
"_Polarizability composition of B": "_PozCB: Polarizability composition of B",
"_Polarizability composition of C": "_PozCC: Polarizability composition of C",
"_Polarizability transition between A and B": "_PzTAB: Polarizability transition between A and B",
"_Polarizability transition between A and C": "_PzTAC: Polarizability transition between A and C",
"_Polarizability transition between B and C": "_PzTBC: Polarizability transition between B and C",
"_Polarizability distribution of 0.00A": "_PzDA0: Polarizability distribution of 0.00A",
"_Polarizability distribution of 0.25A": "_PzDA1: Polarizability distribution of 0.25A",
"_Polarizability distribution of 0.50A": "_PzDA2: Polarizability distribution of 0.50A",
"_Polarizability distribution of 0.75A": "_PzDA3: Polarizability distribution of 0.75A",
"_Polarizability distribution of 1.00A": "_PzDA4: Polarizability distribution of 1.00A",
"_Polarizability distribution of 0.00B": "_PzDB0: Polarizability distribution of 0.00B",
"_Polarizability distribution of 0.25B": "_PzDB1: Polarizability distribution of 0.25B",
"_Polarizability distribution of 0.50B": "_PzDB2: Polarizability distribution of 0.50B",
"_Polarizability distribution of 0.75B": "_PzDB3: Polarizability distribution of 0.75B",
"_Polarizability distribution of 1.00B": "_PzDB4: Polarizability distribution of 1.00B",
"_Polarizability distribution of 0.00C": "_PzDC0: Polarizability distribution of 0.00C",
"_Polarizability distribution of 0.25C": "_PzDC1: Polarizability distribution of 0.25C",
"_Polarizability distribution of 0.50C": "_PzDC2: Polarizability distribution of 0.50C",
"_Polarizability distribution of 0.75C": "_PzDC3: Polarizability distribution of 0.75C",
"_Polarizability distribution of 1.00C": "_PzDC4: Polarizability distribution of 1.00C",
"_Solvent accessibility composition of A": "_SoACA: Solvent accessibility composition of A",
"_Solvent accessibility composition of B": "_SoACB: Solvent accessibility composition of B",
"_Solvent accessibility composition of C": "_SoACC: Solvent accessibility composition of C",
"_Solvent accessibility transition between A and B": "_SATAB: Solvent accessibility transition between A and B",
"_Solvent accessibility transition between A and C": "_SATAC: Solvent accessibility transition between A and C",
"_Solvent accessibility transition between B and C": "_SATBC: Solvent accessibility transition between B and C",
"_Solvent accessibility distribution of 0.00A": "_SADA0: Solvent accessibility distribution of 0.00A",
"_Solvent accessibility distribution of 0.25A": "_SADA1: Solvent accessibility distribution of 0.25A",
"_Solvent accessibility distribution of 0.50A": "_SADA2: Solvent accessibility distribution of 0.50A",
"_Solvent accessibility distribution of 0.75A": "_SADA3: Solvent accessibility distribution of 0.75A",
"_Solvent accessibility distribution of 1.00A": "_SADA4: Solvent accessibility distribution of 1.00A",
"_Solvent accessibility distribution of 0.00B": "_SADB0: Solvent accessibility distribution of 0.00B",
"_Solvent accessibility distribution of 0.25B": "_SADB1: Solvent accessibility distribution of 0.25B",
"_Solvent accessibility distribution of 0.50B": "_SADB2: Solvent accessibility distribution of 0.50B",
"_Solvent accessibility distribution of 0.75B": "_SADB3: Solvent accessibility distribution of 0.75B",
"_Solvent accessibility distribution of 1.00B": "_SADB4: Solvent accessibility distribution of 1.00B",
"_Solvent accessibility distribution of 0.00C": "_SADC0: Solvent accessibility distribution of 0.00C",
"_Solvent accessibility distribution of 0.25C": "_SADC1: Solvent accessibility distribution of 0.25C",
"_Solvent accessibility distribution of 0.50C": "_SADC2: Solvent accessibility distribution of 0.50C",
"_Solvent accessibility distribution of 0.75C": "_SADC3: Solvent accessibility distribution of 0.75C",
"_Solvent accessibility distribution of 1.00C": "_SADC4: Solvent accessibility distribution of 1.00C",
"_Surface tension composition of A": "_SfTCA: Surface tension composition of A",
"_Surface tension composition of B": "_SfTCB: Surface tension composition of B",
"_Surface tension composition of C": "_SfTCC: Surface tension composition of C",
"_Surface tension transition between A and B": "_STTAB: Surface tension transition between A and B",
"_Surface tension transition between A and C": "_STTAC: Surface tension transition between A and C",
"_Surface tension transition between B and C": "_STTBC: Surface tension transition between B and C",
"_Surface tension distribution of 0.00A": "_STDA0: Surface tension distribution of 0.00A",
"_Surface tension distribution of 0.25A": "_STDA1: Surface tension distribution of 0.25A",
"_Surface tension distribution of 0.50A": "_STDA2: Surface tension distribution of 0.50A",
"_Surface tension distribution of 0.75A": "_STDA3: Surface tension distribution of 0.75A",
"_Surface tension distribution of 1.00A": "_STDA4: Surface tension distribution of 1.00A",
"_Surface tension distribution of 0.00B": "_STDB0: Surface tension distribution of 0.00B",
"_Surface tension distribution of 0.25B": "_STDB1: Surface tension distribution of 0.25B",
"_Surface tension distribution of 0.50B": "_STDB2: Surface tension distribution of 0.50B",
"_Surface tension distribution of 0.75B": "_STDB3: Surface tension distribution of 0.75B",
"_Surface tension distribution of 1.00B": "_STDB4: Surface tension distribution of 1.00B",
"_Surface tension distribution of 0.00C": "_STDC0: Surface tension distribution of 0.00C",
"_Surface tension distribution of 0.25C": "_STDC1: Surface tension distribution of 0.25C",
"_Surface tension distribution of 0.50C": "_STDC2: Surface tension distribution of 0.50C",
"_Surface tension distribution of 0.75C": "_STDC3: Surface tension distribution of 0.75C",
"_Surface tension distribution of 1.00C": "_STDC4: Surface tension distribution of 1.00C",
"_vdW volume composition of A": "_vdWCA: vdW volume composition of A",
"_vdW volume composition of B": "_vdWCB: vdW volume composition of B",
"_vdW volume composition of C": "_vdWCC: vdW volume composition of C",
"_vdW volume transition between A and B": "_vWTAB: vdW volume transition between A and B",
"_vdW volume transition between A and C": "_vWTAC: vdW volume transition between A and C",
"_vdW volume transition between B and C": "_vWTBC: vdW volume transition between B and C",
"_vdW volume distribution of 0.00A": "_vWDA0: vdW volume distribution of 0.00A",
"_vdW volume distribution of 0.25A": "_vWDA1: vdW volume distribution of 0.25A",
"_vdW volume distribution of 0.50A": "_vWDA2: vdW volume distribution of 0.50A",
"_vdW volume distribution of 0.75A": "_vWDA3: vdW volume distribution of 0.75A",
"_vdW volume distribution of 1.00A": "_vWDA4: vdW volume distribution of 1.00A",
"_vdW volume distribution of 0.00B": "_vWDB0: vdW volume distribution of 0.00B",
"_vdW volume distribution of 0.25B": "_vWDB1: vdW volume distribution of 0.25B",
"_vdW volume distribution of 0.50B": "_vWDB2: vdW volume distribution of 0.50B",
"_vdW volume distribution of 0.75B": "_vWDB3: vdW volume distribution of 0.75B",
"_vdW volume distribution of 1.00B": "_vWDB4: vdW volume distribution of 1.00B",
"_vdW volume distribution of 0.00C": "_vWDC0: vdW volume distribution of 0.00C",
"_vdW volume distribution of 0.25C": "_vWDC1: vdW volume distribution of 0.25C",
"_vdW volume distribution of 0.50C": "_vWDC2: vdW volume distribution of 0.50C",
"_vdW volume distribution of 0.75C": "_vWDC3: vdW volume distribution of 0.75C",
"_vdW volume distribution of 1.00C": "_vWDC4: vdW volume distribution of 1.00C",
"_Secondary structure composition of A": "_SSCCA: Secondary structure composition of A",
"_Secondary structure composition of B": "_SSCCB: Secondary structure composition of B",
"_Secondary structure composition of C": "_SSCCC: Secondary structure composition of C",
"_Secondary structure transition between A and B": "_SSTAB: Secondary structure transition between A and B",
"_Secondary structure transition between A and C": "_SSTAC: Secondary structure transition between A and C",
"_Secondary structure transition between B and C": "_SSTBC: Secondary structure transition between B and C",
"_Secondary structure distribution of 0.00A": "_SSDA0: Secondary structure distribution of 0.00A",
"_Secondary structure distribution of 0.25A": "_SSDA1: Secondary structure distribution of 0.25A",
"_Secondary structure distribution of 0.50A": "_SSDA2: Secondary structure distribution of 0.50A",
"_Secondary structure distribution of 0.75A": "_SSDA3: Secondary structure distribution of 0.75A",
"_Secondary structure distribution of 1.00A": "_SSDA4: Secondary structure distribution of 1.00A",
"_Secondary structure distribution of 0.00B": "_SSDB0: Secondary structure distribution of 0.00B",
"_Secondary structure distribution of 0.25B": "_SSDB1: Secondary structure distribution of 0.25B",
"_Secondary structure distribution of 0.50B": "_SSDB2: Secondary structure distribution of 0.50B",
"_Secondary structure distribution of 0.75B": "_SSDB3: Secondary structure distribution of 0.75B",
"_Secondary structure distribution of 1.00B": "_SSDB4: Secondary structure distribution of 1.00B",
"_Secondary structure distribution of 0.00C": "_SSDC0: Secondary structure distribution of 0.00C",
"_Secondary structure distribution of 0.25C": "_SSDC1: Secondary structure distribution of 0.25C",
"_Secondary structure distribution of 0.50C": "_SSDC2: Secondary structure distribution of 0.50C",
"_Secondary structure distribution of 0.75C": "_SSDC3: Secondary structure distribution of 0.75C",
"_Secondary structure distribution of 1.00C": "_SSDC4: Secondary structure distribution of 1.00C"
}





# In[4]:
def calculateAA(ProteinSequence, AALetter):
    n = len(AALetter)
    result = {}
    for i in AALetter:
        num = ProteinSequence.count(i)
        rate = round(float(num) / n * 100, 3)
        result[i] = rate

    return result


def StringtoNum(ProteinSequence, AAProperty):
    hardProteinSequence = copy.deepcopy(ProteinSequence)
    for k, m in AAProperty.items():
        for index in m:
            hardProteinSequence = str(hardProteinSequence)
            hardProteinSequence = str.replace(hardProteinSequence, index, k)  ###将序列所有氨基酸残基转变成1或2或3的数字
    TProteinSequence = hardProteinSequence

    return TProteinSequence


# In[5]:


def CalculateComposition(ProteinSequence, AAProperty, AAPName):
    TProteinSequence = StringtoNum(ProteinSequence, AAProperty)

    Result = {}
    Num = len(TProteinSequence)  #####seq的长度
    featurename01 = AAPName +' composition of' +  ' A'
    featurename02 = AAPName + ' composition of' + ' B'
    featurename03 = AAPName + ' composition of' + ' C'

    Result[featurename_dict[featurename01]] = round(float(TProteinSequence.count('1')) / Num, 3)  ####round(x, y), y是x保留小数点几位
    Result[featurename_dict[featurename02]] = round(float(TProteinSequence.count('2')) / Num, 3)
    Result[featurename_dict[featurename03]] = round(float(TProteinSequence.count('3')) / Num, 3)
    return Result


# In[6]:


def CalculateTransition(ProteinSequence, AAProperty, AAPName):
    TProteinSequence = StringtoNum(ProteinSequence, AAProperty)
    Result = {}
    Num = len(TProteinSequence)
    CTD = TProteinSequence

    featurename01 = AAPName + ' transition' + ' between A and B'
    featurename02 = AAPName + ' transition' + ' between A and C'
    featurename03 = AAPName + ' transition' + ' between B and C'

    Result[featurename_dict[featurename01]] = round(float(CTD.count('12') + CTD.count('21')) / (Num - 1), 3)
    Result[featurename_dict[featurename02]] = round(float(CTD.count('13') + CTD.count('31')) / (Num - 1), 3)
    Result[featurename_dict[featurename03]] = round(float(CTD.count('23') + CTD.count('32')) / (Num - 1), 3)
    return Result


# In[7]:


def CalculateDistribution(ProteinSequence, AAProperty, AAPName):
    TProteinSequence = StringtoNum(ProteinSequence, AAProperty)
    Result = {}
    Num = len(TProteinSequence)
    temp = ('1', '2', '3')
    for i in temp:
        num = TProteinSequence.count(i)
        ink = 1
        indexk = 0
        cds = []
        while ink <= num:
            indexk = str.find(TProteinSequence, i, indexk) + 1
            cds.append(indexk)
            ink = ink + 1
        dict_name = {
            '1':'A',
            '2':'B',
            '3':'C'
        }
        featurename01 = AAPName + ' distribution of ' + '0.00' + dict_name[i]
        featurename02 = AAPName + ' distribution of ' + '0.25' + dict_name[i]
        featurename03 = AAPName + ' distribution of ' + '0.50' + dict_name[i]
        featurename04 = AAPName + ' distribution of ' + '0.75' + dict_name[i]
        featurename05 = AAPName + ' distribution of ' + '1.00' + dict_name[i]

        if cds == []:
            Result[featurename_dict[featurename01] ] = 0
            Result[featurename_dict[featurename02] ] = 0
            Result[featurename_dict[featurename03]] = 0
            Result[featurename_dict[featurename04]] = 0
            Result[featurename_dict[featurename05] ] = 0
        else:

            Result[featurename_dict[featurename01]] = round(float(cds[0]) / Num * 100, 3)
            Result[featurename_dict[featurename02]] = round(float(cds[int(math.floor(num * 0.25)) - 1]) / Num * 100, 3)
            Result[featurename_dict[featurename03]] = round(float(cds[int(math.floor(num * 0.5)) - 1]) / Num * 100, 3)
            Result[featurename_dict[featurename04]] = round(float(cds[int(math.floor(num * 0.75)) - 1]) / Num * 100, 3)
            Result[featurename_dict[featurename05]] = round(float(cds[-1]) / Num * 100, 3)

    return Result


# In[8]:

def CalculateCompositionSurfaceTension(ProteinSequence):
    result = CalculateComposition(ProteinSequence, _SurfaceTension, '_Surface tension')
    return result


def CalculateCompositionHydrophobicity(ProteinSequence):
    result = CalculateComposition(ProteinSequence, _Hydrophobicity, '_Hydrophobicity')
    return result


# In[9]:


def CalculateCompositionNormalizedVDWV(ProteinSequence):
    result = CalculateComposition(ProteinSequence, _NormalizedVDWV, '_vdW volume')
    return result


# In[10]:


def CalculateCompositionPolarity(ProteinSequence):
    result = CalculateComposition(ProteinSequence, _Polarity, '_Polarity')
    return result


# In[11]:


def CalculateCompositionCharge(ProteinSequence):
    result = CalculateComposition(ProteinSequence, _Charge, '_Electric charge')
    return result


# In[12]:


def CalculateCompositionSecondaryStr(ProteinSequence):
    result = CalculateComposition(ProteinSequence, _SecondaryStr, '_Secondary structure')
    return result


# In[13]:


def CalculateCompositionSolventAccessibility(ProteinSequence):
    result = CalculateComposition(ProteinSequence, _SolventAccessibility, '_Solvent accessibility')
    return result


# In[14]:


def CalculateCompositionPolarizability(ProteinSequence):
    result = CalculateComposition(ProteinSequence, _Polarizability, '_Polarizability')
    return result


# In[15]:


def CalculateTransitionSurfaceTension(ProteinSequence):
    result = CalculateTransition(ProteinSequence, _SurfaceTension, '_Surface tension')
    result = CalculateTransition(ProteinSequence, _SurfaceTension, '_Surface tension')
    return result


def CalculateTransitionHydrophobicity(ProteinSequence):
    result = CalculateTransition(ProteinSequence, _Hydrophobicity, '_Hydrophobicity')
    return result


# In[16]:


def CalculateTransitionNormalizedVDWV(ProteinSequence):
    result = CalculateTransition(ProteinSequence, _NormalizedVDWV, '_vdW volume')
    return result


# In[17]:


def CalculateTransitionPolarity(ProteinSequence):
    result = CalculateTransition(ProteinSequence, _Polarity, '_Polarity')
    return result


# In[18]:


def CalculateTransitionCharge(ProteinSequence):
    result = CalculateTransition(ProteinSequence, _Charge, '_Electric charge')
    return result


# In[19]:


def CalculateTransitionSecondaryStr(ProteinSequence):
    result = CalculateTransition(ProteinSequence, _SecondaryStr, '_Secondary structure')
    return result


# In[20]:


def CalculateTransitionSolventAccessibility(ProteinSequence):
    result = CalculateTransition(ProteinSequence, _SolventAccessibility, '_Solvent accessibility')
    return result


# In[21]:


def CalculateTransitionPolarizability(ProteinSequence):
    result = CalculateTransition(ProteinSequence, _Polarizability, '_Polarizability')
    return result


# In[22]:


def CalculateDistributionHydrophobicity(ProteinSequence):
    result = CalculateDistribution(ProteinSequence, _Hydrophobicity, '_Hydrophobicity')
    return result


# In[23]:


def CalculateDistributionNormalizedVDWV(ProteinSequence):
    result = CalculateDistribution(ProteinSequence, _NormalizedVDWV, '_vdW volume')
    return result


# In[24]:
def CalculateDistributionSurfaceTension(ProteinSequence):
    result = CalculateDistribution(ProteinSequence, _SurfaceTension, '_Surface tension')
    return result


def CalculateDistributionPolarity(ProteinSequence):
    result = CalculateDistribution(ProteinSequence, _Polarity, '_Polarity')
    return result


# In[25]:


def CalculateDistributionCharge(ProteinSequence):
    result = CalculateDistribution(ProteinSequence, _Charge, '_Electric charge')
    return result


# In[26]:


def CalculateDistributionSecondaryStr(ProteinSequence):
    result = CalculateDistribution(ProteinSequence, _SecondaryStr, '_Secondary structure')
    return result


# In[27]:


def CalculateDistributionSolventAccessibility(ProteinSequence):
    result = CalculateDistribution(ProteinSequence, _SolventAccessibility, '_Solvent accessibility')
    return result


# In[28]:


def CalculateDistributionPolarizability(ProteinSequence):
    result = CalculateDistribution(ProteinSequence, _Polarizability, '_Polarizability')
    return result


# In[29]:


def CalculateC(ProteinSequence):
    result = {}
    result.update(CalculateCompositionPolarizability(ProteinSequence))
    result.update(CalculateCompositionSolventAccessibility(ProteinSequence))
    # result.update(CalculateCompositionSecondaryStr(ProteinSequence))
    result.update(CalculateCompositionCharge(ProteinSequence))
    result.update(CalculateCompositionPolarity(ProteinSequence))
    result.update(CalculateCompositionNormalizedVDWV(ProteinSequence))
    result.update(CalculateCompositionHydrophobicity(ProteinSequence))
    return result


# In[30]:


def CalculateT(ProteinSequence):
    result = {}
    result.update(CalculateTransitionPolarizability(ProteinSequence))
    result.update(CalculateTransitionSolventAccessibility(ProteinSequence))
    # result.update(CalculateTransitionSecondaryStr(ProteinSequence))
    result.update(CalculateTransitionCharge(ProteinSequence))
    result.update(CalculateTransitionPolarity(ProteinSequence))
    result.update(CalculateTransitionNormalizedVDWV(ProteinSequence))
    result.update(CalculateTransitionHydrophobicity(ProteinSequence))
    return result


# In[31]:


def CalculateD(ProteinSequence):
    result = {}
    result.update(CalculateDistributionPolarizability(ProteinSequence))
    result.update(CalculateDistributionSolventAccessibility(ProteinSequence))
    # result.update(CalculateDistributionSecondaryStr(ProteinSequence))
    result.update(CalculateDistributionCharge(ProteinSequence))
    result.update(CalculateDistributionPolarity(ProteinSequence))
    result.update(CalculateDistributionNormalizedVDWV(ProteinSequence))
    result.update(CalculateDistributionHydrophobicity(ProteinSequence))
    return result


# In[32]:


def CalculateCompositionTransitionDistribution(ProteinSequence,methodsID):

    result = {}
    if methodsID == 1:
        _Charge = {'1': 'KR', '2': 'ANCQGHILMFPSTWYV', '3': 'DE'}
        # '1'stand for Positive; '2'stand for Neutral, '3' stand for Negative

        result.update(CalculateCompositionCharge(ProteinSequence))
        result.update(CalculateTransitionCharge(ProteinSequence))
        result.update(CalculateDistributionCharge(ProteinSequence))
    if methodsID == 2:
        _Hydrophobicity = {'1': 'RKEDQN', '2': 'GASTPHY', '3': 'CLVIMFW'}
        # '1'stand for Polar; '2'stand for Neutral, '3' stand for Hydrophobicity

        result.update(CalculateCompositionHydrophobicity(ProteinSequence))
        result.update(CalculateTransitionHydrophobicity(ProteinSequence))
        result.update(CalculateDistributionHydrophobicity(ProteinSequence))
    if methodsID == 3:
        _Polarity = {'1': 'LIFWCMVY', '2': 'CPNVEQIL', '3': 'KMHFRYW'}
        # '1'stand for (4.9-6.2); '2'stand for (8.0-9.2), '3' stand for (10.4-13.0)

        result.update(CalculateCompositionPolarity(ProteinSequence))
        result.update(CalculateTransitionPolarity(ProteinSequence))
        result.update(CalculateDistributionPolarity(ProteinSequence))
    if methodsID == 4:
        _Polarizability = {'1': 'GASDT', '2': 'CPNVEQIL', '3': 'KMHFRYW'}
        result.update(CalculateCompositionPolarizability(ProteinSequence))
        result.update(CalculateTransitionPolarizability(ProteinSequence))
        result.update(CalculateDistributionPolarizability(ProteinSequence))
    if methodsID == 5:
        _SolventAccessibility = {'1': 'ALFCGIVW', '2': 'RKQEND', '3': 'MPSTHY'}
        # '1'stand for Buried; '2'stand for Exposed, '3' stand for Intermediate
        result.update(CalculateCompositionSolventAccessibility(ProteinSequence))
        result.update(CalculateTransitionSolventAccessibility(ProteinSequence))
        result.update(CalculateDistributionSolventAccessibility(ProteinSequence))
    if methodsID == 6:
        _SurfaceTension = {'1': 'GQDNAHR', '2': 'KTSEC', '3': 'ILMFPWYV'}
        # '1'stand for (0-0.108); '2'stand for (0.128-0.186), '3' stand for (0.219-0.409)
        result.update(CalculateCompositionSurfaceTension(ProteinSequence))
        result.update(CalculateTransitionSurfaceTension(ProteinSequence))
        result.update(CalculateDistributionSurfaceTension(ProteinSequence))
    if methodsID == 7:

        _NormalizedVDWV = {'1': 'GASTPD', '2': 'NVEQIL', '3': 'MHKFRYW'}
        # '1'stand for (0-2.78); '2'stand for (2.95-4.0), '3' stand for (4.03-8.08)
        result.update(CalculateCompositionNormalizedVDWV(ProteinSequence))
        result.update(CalculateTransitionNormalizedVDWV(ProteinSequence))
        result.update(CalculateDistributionNormalizedVDWV(ProteinSequence))
    if methodsID == 8:
        _SecondaryStr = {'1': 'EALMQKRH', '2': 'VIYCWFT', '3': 'GNPSD'}
        # '1'stand for Helix; '2'stand for Strand, '3' stand for coil
        result.update(CalculateCompositionSecondaryStr(ProteinSequence))
        result.update(CalculateTransitionSecondaryStr(ProteinSequence))
        result.update(CalculateDistributionSecondaryStr(ProteinSequence))
    if methodsID == 9:
        AALetter = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
        result.update(calculateAA(ProteinSequence, AALetter))

    return result


def convert_to_array(protein_dict, N):
    protein_array = np.zeros(N, dtype=np.float)
    k_all = []
    for i, (k, v) in enumerate(protein_dict.items()):
        protein_array[i] = v
        k_all.append(k)
    return protein_array, k_all

def makecolname(strdata):
    if len(strdata)==1:
        datares = strdata
    else:
        datares = strdata[1:]
    return datares

def process_fasta(filepath,methodsID):
    all_protein_array = []
    seqname = []
    column = []
    seq_seq = []
    for seq in SeqIO.parse(filepath, 'fasta'):
        seqid = seq.id
        seqname.append(seqid)
        seq_seq.append(seq.seq)
    for protein in seq_seq:
        result = CalculateCompositionTransitionDistribution(protein,methodsID)
        N = len(result)
        protein_array, k_all = convert_to_array(result, N)
        column = list(map(makecolname,k_all))

        all_protein_array.append(protein_array)
    if methodsID == 9:
        column = ["AACoA: Composition of alanine","AACoR: Composition of arginine","AACoN: Composition of asparagine","AACoD: Composition of aspartic acid","AACoC: Composition of cysteine","AACoE: Composition of glutamic acid","AACoQ: Composition of glutamine","AACoG: Composition of glycine","AACoH: Composition of histidine","AACoI: Composition of isoleucine","AACoL: Composition of leucine","AACoK: Composition of lysine","AACoM: Composition of methionine","AACoF: Composition of phenylalanine","AACoP: Composition of proline","AACoS: Composition of serine","AACoT: Composition of threonine","AACoW: Composition of tryptophan","AACoY: Composition of tyrosine","AACoV: Composition of valine"]

    df_cod = DataFrame(data=all_protein_array, index=seqname, columns=column)
    return df_cod
