dates_file="covid19_seq_info_5.25.20.csv"
original_file="covid19_4.17.20.aln"
newfile="covid19_seq_alignments_wDates.aln"

accession_date_dict=dict()

with open(dates_file,"r") as f:
    for line in f.readlines():
        line=line.split(",")
        accession_date_dict.update({line[0]:line[1][0:10]})

with open(original_file,"r") as o:
    with open(newfile,"w") as n:
        for line in o.readlines():
            if line.startswith(">"):
                accession=line[line.find("|")+1:-3]
                print(accession)
                line=line.strip("\n")+"|"+accession_date_dict[accession]+"\n"
            n.write(line)


