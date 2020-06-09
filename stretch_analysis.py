"""
The code uses a pre-prepared consensus file and protein fasta.
The matrix parameter must be predefined in 'immutability.py'.
 Using these finds a mutability score array. this array .
 Uses a range of cutoffs that filter the sequence based on its mutability score.
"""

from immutability02 import *
from numpy import arange
# ANALYSIS
STRETCH_SIZE_RANGE = range(18, 25)
CUTOFF_RANGE= arange(-1.0, 2.5, .05)

def find_continuous_stretches(stretch_length,sequence):
    continuous_stretches = 0
    for k in range(len(sequence) - stretch_length):
        if ("".join(sequence[k:k + stretch_length])).__contains__("N") == False: continuous_stretches += 1
    return continuous_stretches

if __name__=="__main__":
    results=[]

    for i in CUTOFF_RANGE:
        results.append(calc_immutability("covid19_4.17.20_consensus.txt", "NC_045512_proteins_fasta.txt", i))

    print ( "cutoff","sum"," ".join(["s"+str(x) for x in STRETCH_SIZE_RANGE]))
        #COLUMNS:   1- cutoffs used
        #           2- sum of nonN nucleotide count in immutable sequence
        #           3 and on- number of non-N-containing continuous strecthes of varying sizes.
    for idx, i in enumerate(CUTOFF_RANGE):

        continuous_stretches=[find_continuous_stretches(x,results[idx][0]) for x in STRETCH_SIZE_RANGE]

        print(i, len([x for x in results[idx][0] if x != "N"]),
              "\t".join(["{:.2f}".format(x) for x in continuous_stretches]))

