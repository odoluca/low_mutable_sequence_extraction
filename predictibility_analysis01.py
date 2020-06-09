from immutability02 import calc_immutability
from consensus import extract_consensus
import datetime as dt
from random import sample,shuffle
from matplotlib import pyplot as plt
from numpy import mean,std
from scipy.stats import norm


#TODO:  get list of aligned FASTAs and extract dates. save to a list
#TODO:  for each date:
#TODO:      save all before a date to a temp file
#TODO:      find a consensus sequence using "consensus.py"
#TODO:      from the rest find new mutations by comparing with consensus sequence, keep as new_mutations_list
#TODO:      create 500 random mutations_lists from new_mutations_list
#TODO:      for each mutation list:
#TODO:          sum the corresponding immobility scores for each mutation in randomlists and record to a list
#TODO:      sum the corresponding immobility scores for each mutation in new_mutation_list
#TODO:      calculate the z-score for the new_mutation_list.

#parameters
dated_alignment_filename="covid19_seq_alignments_wDates.aln"
Refseq_protein_seq_filename="NC_045512_proteins_fasta.txt"
CUTOFF=1
VERBOSE=True


FASTAs=[]

#region extracts Fastas and their dates
with open(dated_alignment_filename,"r") as f:
    new_seq = ""
    for line in f.readlines():
        if line.startswith(">"):

            if new_seq!="": FASTAs.append((new_seq,new_seq_date))
            new_seq_date = dt.date.fromisoformat(line.strip("\n").split("|")[2])
            new_seq=""
        new_seq+=line
FASTAs.sort(key=lambda x:x[1]) #NOW FASTAs CONTAIN SEQUENCES AND DATES
#endregion

#region creates a range of dates
# date_cutoff=dt.date.fromisoformat("2020-03-15")
date_range=[dt.date.fromisoformat("2020-03-15") + dt.timedelta(days=x) for x in range(40,50,1)]
#endregion


report_file=open("predictibiliy_report.txt","w") #creates a report file
report_file.write("cutoff_date,old_data_count,new_data_count,new_mutations,unpredicted_mutations,unpredicted_mut_mutability_scores_sum,z_scores,p-values,mean_of_randoms,std_of_randoms\n")

def shuffle_nonFalse(new_consensus_prot_regions,mutability_score):

    a=[x for x in new_consensus_prot_regions if x!=False]
    shuffle(a)

    b=[]
    p=0
    for i in range(new_consensus_prot_regions.__len__()):
        if new_consensus_prot_regions[i]!=False:
            b.append(a[p])
            p+=1
        else:
            b.append(False)
    return b

def mutability_score_sum( mutated_seq,mutability_score):
    return sum([mutability_score[x] for x in range(len(mutability_score)) if (
            mutability_score[x] != False and mutated_seq[x] == "N")])



mutability_profiles={}
#region process each dates
for date_cutoff in date_range:
    print ("Date cufoff:",date_cutoff)
    with open("temp_old.aln","w") as temp_file:
        temp_file.write("".join([x[0] for x in FASTAs if x[1]<=date_cutoff]))
        old_data_count=len([x[0] for x in FASTAs if x[1]<=date_cutoff])
    with open("temp_new.aln","w") as temp_file:
        temp_file.write("".join([x[0] for x in FASTAs if x[1]>date_cutoff]))
        new_data_count=len([x[0] for x in FASTAs if x[1]>date_cutoff])
    consensus_seq=extract_consensus("temp_old.aln",95)
    immutable_seq,mutability_score=calc_immutability("temp_old_consensus.txt","NC_045512_proteins_fasta.txt",CUTOFF)
    mutability_profiles[date_cutoff]=mutability_score
    new_consensus= extract_consensus("temp_new.aln",100) #includes possibly previously happened mutations too.

    new_consensus_in_mutability_regions = [new_consensus[x] if mutability_score[x] != None else False for x in
                                           range(len(new_consensus))]
    #randomized sequences (randomized only at protein regions, and the rest is False)
    random_consensuses=[shuffle_nonFalse(new_consensus_in_mutability_regions, mutability_score) for x in range(300)] #shuffles to generate a new sequence but only using sequences at protein coding regions
    assert new_consensus_in_mutability_regions.count("N") == random_consensuses[0].count("N"), "issue with random sequence generation"
    assert new_consensus_in_mutability_regions.count(False) == random_consensuses[0].count(False), "issue with random sequence generation"

    #calculated mutability score sum for missed mutations
    unpredicted_new_mut_mutability_score_sum =mutability_score_sum(new_consensus_in_mutability_regions, mutability_score)

    #list of mutability score sums for randomly prepared sequences
    random_mutability_score_sums= [mutability_score_sum(random_list, mutability_score) for random_list in random_consensuses]
    mean_of_random=mean(random_mutability_score_sums)
    std_of_random=std(random_mutability_score_sums)
    z_score=(unpredicted_new_mut_mutability_score_sum-mean_of_random )/std_of_random
    p_value=norm.sf(z_score)*2

    if VERBOSE:
        n, bins, patches = plt.hist(x=[x for x in random_mutability_score_sums if x != None], bins='auto', color='#0504aa',
                                    alpha=0.7, rwidth=0.85)
        plt.title("mutability score sum distribution for random sequences")
        plt.xlabel("mutability score sum")
        plt.ylabel("number of randomized sequences")
        plt.show()

    print("sum of unpredicted mutability scores:",unpredicted_new_mut_mutability_score_sum,"z-score",z_score,"p-value",p_value)


    new_mutation_count=[ (consensus_seq[x]!="N" and new_consensus[x]=="N") for x in range(len(mutability_score)) ].count(True)

    unpredicted_new_mut_count=[ (immutable_seq[x]!="N" and new_consensus[x]=="N") for x in range(len(mutability_score)) ].count(True)



    print(new_mutation_count,unpredicted_new_mut_count)
    report_file.write(date_cutoff.__str__()+","+old_data_count.__str__()+","+new_data_count.__str__()+","+new_mutation_count.__str__()+","+unpredicted_new_mut_count.__str__()+","+unpredicted_new_mut_mutability_score_sum.__str__()+","+z_score.__str__()+","+p_value.__str__()+","+mean_of_random.__str__()+","+std_of_random.__str__()+"\n")
    report_file.flush()
report_file.close()

