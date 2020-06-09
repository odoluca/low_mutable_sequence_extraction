from Bio.SubsMat import MatrixInfo as MI
from Bio.pairwise2 import align
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast.NCBIXML import parse,read
from Bio.Blast.Applications import NcbitblastnCommandline,NcbiblastxCommandline
from Bio.Data.CodonTable import CodonTable
from numpy import std,mean
from sys import argv
from os.path import isfile
from random import uniform


#aa substitution scoring methods

pam=lambda x,y: MI.pam30[x,y]
blo=lambda x,y: MI.blosum30[x,y]
ran=lambda x,y: uniform(-10,10)

# def imutant()

HKYI={
    'A':{'T': 6.0,
    'C': 3.45,
    'G': 3.45},
    'T':{'A': 5.6,
    'C': 11.52,
    'G': 3.68},
    'C':{'A': 5.6,
    'T': 20.05,
    'G': 3.68},
    'G':{'A': 18.71,
    'T': 6.0,
    'C': 3.45}
} #obtained substitution matrix based on HKY+G model using MEGA-X and 146 FASTAs of Covid-19 colelcted from NCBI until 3/14/20.
GTRGI={
    'A':{'T': 0.062,
    'C': 0.037,
    'G': 0.076},
    'T':{'A': 0.058,
    'C': 0.142,
    'G': 0.054},
    'C':{'A': 0.06,
    'T': 0.248,
    'G': 0.031},
    'G':{'A': 0.116,
    'T': 0.088,
    'C': 0.029}
} #best model obtained (BIC=104559.587) using covid19 data up to 4th of April aligned by NCBI and analysed by MEGA X. For estimating ML values, a tree topology was automatically computed. This analysis involved 320 nucleotide FASTAs. Codon positions included were 1st+2nd+3rd+Noncoding.
codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
#SETTINGS
CUTOFF=1
UNCONSERVED_SCORE=None
#PARAMETERS
aa_score_method=blo
consensus_filename="covid19_4.17.20_consensus.txt"
protein_seqs_filename="NC_045512_proteins_fasta.txt"
SUBSTITUTION_MODEL=GTRGI



def replace(str,char,pos):
    return str[:pos]+char+str[pos+1:]



#TODO: iterate through protein_alignments
#TODO:      iterate through aminoacids
#TODO:          iterate through associated nucleotides
#TODO:              calculate nucleotide conservation probability and add to consensus score


def calc_immutability(consensus_filename,protein_seqs_filename,CUTOFF):
    out = NcbitblastnCommandline(subject=consensus_filename, query=protein_seqs_filename, out="results.xml",
                                 outfmt=5)  # evalue=?
    out()  # NCBI tblastn is run locally to align proteins to the consensus sequence

    consensus_seq = SeqIO.parse(consensus_filename, "fasta").__next__().seq._data
    protein_alignments = [x.alignments[0].hsps[0] for x in parse(open(out.out, "r"))]

    mutability_score = [UNCONSERVED_SCORE for _ in range(len(consensus_seq))]
    for alignment in protein_alignments:
        position=alignment.sbjct_start-1

        # print("\ncalculating for",alignment.sbjct)
        for aa in alignment.sbjct:
            assert aa!=" ", "empty aa in alignment subject"

            if aa!="X" :
                codon0= consensus_seq[position:position + 3]
                # print(aa, end="")
                aa0="none"
                try: aa0 = codon_table[codon0]
                except:
                    position+=3
                    continue    #We do not calculate if there is an ambiguous nucleotide in the consensus seq. We could
                                    #have calculated for the other nucleotides in the codon however it requires more
                                    #detailed calculation.

                # print(position,codon0,aa0,end="|")

                for pos_in_codon in range(3):
                    # codonX=codon0
                    if mutability_score[position + pos_in_codon]!=UNCONSERVED_SCORE:continue # if this position was previously processed, there is no reason to process again. SKIP IT. This might be the case if there are several proteins aligned to here.

                    mutability_score[position + pos_in_codon]=0

                    for putative_nuc in SUBSTITUTION_MODEL[codon0[pos_in_codon]].keys():

                        codonX=replace(codon0, putative_nuc, pos_in_codon)
                        aaX = codon_table[codonX]

                        if  aaX!="_" :#and aa0!=aaX: #To include self conversion sounds rigth way as both matrices include these in their calculations. Asuming any deleterious mutation would not be conserved thus not needed for consideration.

                            try:  aa_score=aa_score_method(aa0,aaX)
                            except: aa_score=aa_score_method(aaX,aa0) #pam30 is a triangular matrix. so we check reversed
                            subs_score=SUBSTITUTION_MODEL[codon0[pos_in_codon]][putative_nuc]
                            # print(aa0,"-->",aaX,":",aa_score*subs_score,end="|")

                            mutability_score[position + pos_in_codon]+= aa_score * subs_score
            position+=3

    print("\nconsevation scores complete:")

    standart_dev=std([x for x in mutability_score if x!=None])
    average=mean([x for x in mutability_score if x!=None])

    printlines=["","","","","",""]
    for i in range(len(consensus_seq)):
        if i%100==0:
            print(i,end=" "*(100-len(str(i))))
            printlines[0]+=str(i)+" "*(100-len(str(i)))
    print()
    for i in range(len(consensus_seq)):
        if i%100==0:
            print("|",end="")
            printlines[1]+="|"
        else:
            print(" ",end="")
            printlines[1] += " "
    print()

    prot="" #TODO: automatically find where protein starts, add all one by one. This code below assumes all has ORF +1.It would not work for all cases
    for i in range(2,len(consensus_seq),3):
        d="   "
        if mutability_score[i]!=None:
            try:
                d=" "+codon_table[consensus_seq[i:i+3]]+" "
            except:
                pass
        prot+=d
    print("  "+prot)
    # printlines[2]="  "+prot

    for i in range(len(consensus_seq)):
        print(consensus_seq[i],end="")
        printlines[3]+=consensus_seq[i]
    print()
    for i in range(len(consensus_seq)):
        d="_"
        if mutability_score[i] !=None:
            d=" "
            if mutability_score[i] < average+standart_dev*2:
                d="░"
            if mutability_score[i] < average+standart_dev:
                d="▒"
            if mutability_score[i] < average -standart_dev:
                d = "▓"
            if mutability_score[i] < average-standart_dev*2:
                d = "█"
        print(d,end="")
        printlines[4]+=d
    print()
    immutable_seq= ""
    for i in range(len(consensus_seq)):
        d="N"
        if mutability_score[i]!=None and mutability_score[i]<CUTOFF:
            d=consensus_seq[i]
        print(d,end="")
        printlines[5]+=d

        immutable_seq +=d

    # print()
    # for i in range(0,30000,100):
    #     print("\n".join(printlines[i:i+100]))

    print("\n<<<REPORT>>>")

    print("settings - Score for unconserved nucleotides:",UNCONSERVED_SCORE,"Cutoff score:",CUTOFF)
    print("max score:",max([x for x in mutability_score if x!=None]),"min score:",min([x for x in mutability_score if x!=None]))
    print("mean score",average,"standart dev:",standart_dev)
    print("number of non-Ns in consensus seq:",len([x for x in consensus_seq if x!="N"]))
    print("number of non-Ns in immutable seq:",len([x for x in immutable_seq if x!="N"]))
    print("scores of coding regions summed:", sum([x for x in mutability_score if x!=None]))
    # print("scores of all regions summed:", sum([x for x in mutability_score]))
    return immutable_seq, mutability_score

if __name__=="__main__":

    for idx, arg in enumerate(argv):
        if arg == "-cseq":
            consensus_filename = argv[idx + 1]
            assert isfile(consensus_filename), "No consensus file find."
        elif arg == "-pseq":
            protein_seqs_filename = argv[idx + 1]
            assert isfile(protein_seqs_filename), "No protein sequence file found."
        elif arg=="-c":
            CUTOFF=float(argv[idx+1])
        elif arg.startswith("-"):
            exit("Cannot recognise parameter:"+arg)


    immutable_seq, mutability_score=calc_immutability(consensus_filename, protein_seqs_filename, CUTOFF)

    import matplotlib.pyplot as plt

    # bins=plt.hist([ x for x in mutability_score if x!=None],20)
    # plt.show()
    n, bins, patches = plt.hist(x=[x for x in mutability_score if x!=None], bins='auto', color='#0504aa',
                                alpha=0.7, rwidth=0.85)
    plt.title("mutability score distribution")
    plt.xlabel("mutability score")
    plt.ylabel("number of residues")
    plt.show()

    from pickle import dump
    dump((immutable_seq, mutability_score), open("../primer_evaluation/immutability_score.data", "wb"))


    from numpy import arange
    t=open("immutable_sequences.txt","w")
    for co in arange(0.5,3.1,0.1):
        immutable_seq, mutability_score = calc_immutability(consensus_filename, protein_seqs_filename, co)
        t.write(co.__str__()+"\n")
        t.write(immutable_seq+"\n")
        t.write("\n")
