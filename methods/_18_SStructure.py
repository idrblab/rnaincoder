import os
import pandas as pd
import Bio.SeqIO as Seq
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
# pandas2ri.activate()
import RNA
# os.environ['R_HOME'] = "/usr/bin/R" #or whereever your R is installed

def extract_SSfeatures(fastapath):
    # import rpy2.robjects as robjects
    # from rpy2.robjects import pandas2ri
    # pandas2ri.activate()

    robjects.r('''
                extract_SSfeatures <- function(fastapath){
                  library(LncFinder)
                  demo_DNA.seq <- seqinr::read.fasta(fastapath)
                  Seqs <- LncFinder::run_RNAfold(demo_DNA.seq, RNAfold.path = "RNAfold", parallel.cores = 2)
                  result_2 <- LncFinder::extract_features(Seqs, label = NULL, SS.features = TRUE,format = "SS", frequencies.file = "human", parallel.cores = -1)
                  res2 <- result_2[,c(12:19)]
                  return(res2)
                }''')
    sstruc = robjects.r['extract_SSfeatures'](fastapath)
    sstruc.columns = ["SLDLD: Structural logarithm distance to lncRNA of acguD", "SLDPD: Structural logarithm distance to pcRNA of acguD", "SLDRD: Structural logarithm distance acguD ratio", "SLDLN: Structural logarithm distance to lncRNA of acguACGU", "SLDPN: Structural logarithm distance to pcRNA of acguACGU", "SLDRN: Structural logarithm distance acguACGU ratio","SDMFE: Secondary structural minimum free energy", "SFPUS: Secondary structural UP frequency paired-unpaired"]
    return sstruc


def makeEIIP(fastapath):
    # import rpy2.robjects as robjects
    # from rpy2.robjects import pandas2ri
    # pandas2ri.activate()
                                #     spectrum.percent = 0.1,
                                # quantile.probs = seq(0, 1, 0.25)

    robjects.r('''
                makeEIIP <- function(fastapath){
                  library(LncFinder)
                  demo_DNA.seq <- seqinr::read.fasta(fastapath)
                  result_1 <- compute_EIIP(
                                demo_DNA.seq,
                                label = NULL,
                                spectrum.percent = 0,
                                quantile.probs = seq(0, 1, 0.25)
                                )
                  return(result_1)
                }''')
    sstruc = robjects.r['makeEIIP'](fastapath)
    sstruc.columns = ['EipSP: Electron-ion interaction pseudopotential signal peak','EipAP: Electron-ion interaction pseudopotential average power','EiSNR: Electron-ion interaction pseudopotential signal/noise ratio','EiPS0: Electron-ion interaction pseudopotential spectrum 0','EiPS1: Electron-ion interaction pseudopotential spectrum 0.25','EiPS2: Electron-ion interaction pseudopotential spectrum 0.5','EiPS3: Electron-ion interaction pseudopotential spectrum 0.75','EiPS4: Electron-ion interaction pseudopotential spectrum 1']
    return sstruc

def makeORFEucDist(fastapath):
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    # pandas2ri.activate()

                  # cds.seq = seqinr::read.fasta('/public/home/wangyx/LncRNA/smallRNA/methods/Data/gencode.v34.pc_transcripts_test.fa')
                  # lncRNA.seq = seqinr::read.fasta('/public/home/wangyx/LncRNA/smallRNA/methods/Data/gencode.v34.lncRNA_transcripts_test.fa')
                  # cds.seq = seqinr::read.fasta('/public/home/wangyx/LncRNA/smallRNA/methods/Data/gencode.v34.pc_transcripts.fa')
                  # lncRNA.seq = seqinr::read.fasta('/public/home/wangyx/LncRNA/smallRNA/methods/Data/gencode.v34.lncRNA_transcripts.fa')                  
    robjects.r('''
                makeORFEucDist <- function(fastapath){
                  library(LncFinder)
                  cds.seq = seqinr::read.fasta('/public/home/wangyx/LncRNA/smallRNA/methods/Data/gencode.v34.pc_transcripts_test.fa')
                  lncRNA.seq = seqinr::read.fasta('/public/home/wangyx/LncRNA/smallRNA/methods/Data/gencode.v34.lncRNA_transcripts_test.fa')   

                  referFreq <- make_referFreq(
                            cds.seq,
                            lncRNA.seq,
                            k = 6,
                            step = 1,
                            alphabet = c("a", "c", "g", "t"),
                            on.orf = TRUE,
                            ignore.illegal = TRUE
                            )

                  demo_DNA.seq <- seqinr::read.fasta(fastapath)
                  EucDis <- compute_EucDistance(
                                demo_DNA.seq,
                                label = NULL,
                                referFreq,
                                k = 6,
                                step = 1,
                                alphabet = c("a", "c", "g", "t"),
                                on.ORF = TRUE,
                                auto.full = FALSE,
                                parallel.cores = -1
                                )


                LogDistance <- compute_LogDistance(
                                    demo_DNA.seq,
                                    label = NULL,
                                    referFreq,
                                    k = 6,
                                    step = 1,
                                    alphabet = c("a", "c", "g", "t"),
                                    on.ORF = TRUE,
                                    auto.full = FALSE,
                                    parallel.cores = -1
                                    )
                hexamerScore <- compute_hexamerScore(
                                    demo_DNA.seq,
                                    label = NULL,
                                    referFreq,
                                    k = 6,
                                    step = 1,
                                    alphabet = c("a", "c", "g", "t"),
                                    on.ORF = TRUE,
                                    auto.full = FALSE,
                                    parallel.cores = -1
                                    )
                hdata2<-cbind(EucDis,LogDistance)
                result <- cbind(hdata2,hexamerScore)
                return(result)

                }'''
              )
    sstruc = robjects.r['makeORFEucDist'](fastapath)
    return sstruc


    
    
fina_df = makeORFEucDist('/public/home/wangyx/LncRNA/smallRNA/Data/CPPredData_test_human_sorf/Homo38_small_ncRNA_test.fasta')
print(fina_df)
# fina_df = pd.DataFrame(fina_df)
# # fina_df = cal_base_pairs('/public/home/wangyx/LncRNA/smallRNA/data/sequence/RPI369_rna_seq_test.fa')
# fina_df.to_csv('/public/home/wangyx/LncRNA/smallRNA/Data/CPPredData_test_human_sorf/Homo38_small_ncRNA_test.csv')





# fastapath = '/public/home/wangyx/LncRNA/resubmission_code_16methods/Case4_rri_data/01_FastaData/fasta_mtr/mirna.fasta'
# fina_df = extract_SSfeatures(fastapath)
# print(fina_df)
# fina_df = pd.DataFrame(fina_df)
# # fina_df = cal_base_pairs('/public/home/wangyx/LncRNA/smallRNA/data/sequence/RPI369_rna_seq_test.fa')
# fina_df.T.to_csv('/public/home/wangyx/LncRNA/resubmission_code_16methods/Case4_rri_data/01_FastaData/fasta_mtr/mirna.csv')