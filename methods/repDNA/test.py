__author__ = 'Fule Liu'

from nacutil import make_kmer_list


if __name__ == "__main__":
    for e in make_kmer_list("AAA", "ACGT"):
        print(e)

    pass
