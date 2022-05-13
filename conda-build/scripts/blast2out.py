import sys
from Bio import SeqIO

#Open ORF file
orf_file = sys.argv[1]
orfs = open(orf_file)
orf_iterator = SeqIO.parse(orfs,"fasta")

#Get type of sort
sort_type = sys.argv[2]
min_length = int(sys.argv[3])
uniprot = int(sys.argv[4])

#Variables
results = list()
coords = list()
hits = list()
seq_name = " "
orf_name = ""

#Start reading blast output
for line in sys.stdin:
    s = line.split()

    if(orf_name != s[0] and orf_name != ""):
        if(len(coords)>0):
            #Calculate coverage
            coords.sort()
            orf_cov = 0
            start = -1
            end = -1
            for coord in coords:
                if(start == -1):
                    start = coord[0]
                    end = coord[1]
                    continue
                if(coord[0] < end and end < coord[1]):
                    end = coord[1]
                    continue
                if(coord[0] >= end):
                    orf_cov += end - start + 1
                    start = coord[0]
                    end = coord[1]
            if(start != -1):
                orf_cov += end - start + 1
            coords.clear()
            orf_cov = orf_cov/len(orf_seq)
            results.append((orf_name,orf_cov,sorted(hits),orf_seq))
            hits.clear()
            
    while(orf_name != s[0]):
        record = next(orf_iterator)
        orf_name = record.id
    orf_seq = record.seq

    if(len(orf_seq)>=min_length):
        coords.append((int(s[2]),int(s[3])))
        hits.append(s[1][3:-1] + "_" + s[2] + s[3] + "_" + s[4])

if(len(coords)>0):
    #Calculate coverage
    coords.sort()
    orf_cov = 0
    start = -1
    end = -1
    for coord in coords:
        if(start == -1):
            start = coord[0]
            end = coord[1]
            continue
        if(coord[0] < end and end < coord[1]):
            end = coord[1]
            continue
        if(coord[0] >= end):
            orf_cov += end - start + 1
            start = coord[0]
            end = coord[1]
    if(start != -1):
        orf_cov += end - start + 1
    coords.clear()
    orf_cov = orf_cov/len(orf_seq)
    results.append((orf_name,orf_cov,sorted(hits),orf_seq))
    hits.clear()

#Close ORF file
orfs.close()

#Sort fasta output
if(sort_type == "sl"):
    results.sort(reverse = True, key = lambda x: len(x[3]))
if(sort_type == "cv"):
    results.sort(reverse = True, key = lambda x: x[1])
if(sort_type == "ht"):
    results.sort(reverse= True, key = lambda x: len(x[2]))

#Write output
if(len(results)==0):
    exit(0)

for record in results:
    out = ">" + record[0][:-2] + " " + str(record[1]) + " "
    if(uniprot == 1):
        for hit in record[2]:
            out += hit + ","
    print(out[:-1])
    print(record[3])

results.clear()


