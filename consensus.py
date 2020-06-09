
#This code takes a FASTA multiple sequence alignment and discovers only conserved alignment.
alignment_filename = "covid19_4.17.20.aln"

def extract_consensus(alignment_filename,conservation_percent=100):
    sequences=[]

    with open(alignment_filename, "r") as f:
        sequence=""
        for line in f.readlines():
            if line.startswith(">"):
                if sequence!="":
                    sequences.append(sequence)
                sequence=""
            else:
                sequence+=line.strip()
    # print (FASTAs)
    print(sequences.__len__(), "FASTAs found.")
    assert all([x.__len__()==sequences[0].__len__() for x in sequences]),"size check fails!"
    print("size check: OK")

    consensus_seq= ""
    consensus_score=[]
    consensus_threshold=conservation_percent
    for i in range(sequences[0].__len__()):
        letters={"A":0,"T":0,"C":0,"G":0,"N":0,"-":0}
        for s in sequences:
            if s[i] in letters.keys():
                letters[s[i]]+=1
        total_determined=float(letters["A"]+letters["T"]+letters["C"]+letters["G"])
        if total_determined!=0:
            assert consensus_threshold>50,"consensus threshold can not be below 50%"
            if letters["A"]/total_determined>=consensus_threshold/100.0:
                consensus_seq= consensus_seq + "A"
                consensus_score.append(letters["A"]/total_determined)
            elif letters["T"]/total_determined>=consensus_threshold/100.0:
                consensus_seq= consensus_seq + "T"
                consensus_score.append(letters["T"]/total_determined)
            elif letters["C"]/total_determined>=consensus_threshold/100.0:
                consensus_seq= consensus_seq + "C"
                consensus_score.append(letters["C"]/total_determined)
            elif letters["G"]/total_determined>=consensus_threshold/100.0:
                consensus_seq= consensus_seq + "G"
                consensus_score.append(letters["G"]/total_determined)
            else:
                consensus_seq+= "N"
        else:
            consensus_seq+= "N"
    print(consensus_seq)

    assert consensus_seq.__len__() == sequences[0].__len__(), "size check fails!"
    print("size check: OK")
    print("using {} as minimum conservation percent limit".format(conservation_percent))
    new_filename= alignment_filename[:-4] + "_consensus.txt"
    with open(new_filename,"w") as f:
        f.write("> Consensus sequence generated from aligned fasta" + alignment_filename+"\n")
        f.write(consensus_seq)
    print("consensus size: ",consensus_seq.__len__())
    print("consensus saved:",new_filename)
    return consensus_seq

con=extract_consensus(alignment_filename,95)











