import sys

def calculate_seqcount(filename):
    count = 0

    with open(filename) as f:
        for line in f:
            if line.startswith(">"):
                count += 1

    return count

fasta_file = sys.argv[1]
seq_count = calculate_seqcount(fasta_file)
print(seq_count)

#run on command line, python3 script1.py query_genome.fna
#answer is, 5292
#grep -c "^>" query_genome.fna also confirms 

 

