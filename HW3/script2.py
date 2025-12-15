#!/usr/bin/env python3

import argparse
def calculate_num_matched_sequences(filename, min_identity):
    matched_queries = set()
    with open(filename) as f:
        for line in f:
            fields = line.strip().split("\t")
            query_id = fields[0]
            identity = float(fields[2])

            if identity >= min_identity:
                matched_queries.add(query_id)

    return len(matched_queries), matched_queries

def count_total_queries(query_fasta):
    count = 0
    with open(query_fasta) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count

def main():
    parser = argparse.ArgumentParser(description="Parse BLAST output and count matched query sequences")
    parser.add_argument("blast_file", help="BLAST tabular output file")
    parser.add_argument("min_identity", type=float, help="Minimum percent identity threshold")
    parser.add_argument("--ratio", action="store_true", help="Print fraction of matched query sequences")
    parser.add_argument("--query", default="query_genome.fna", help="Query FASTA file")

    args = parser.parse_args()

    num_matched, matched_set = calculate_num_matched_sequences(
        args.blast_file, args.min_identity
    )

    print(f"Number of matched query sequences: {num_matched}")

    if args.ratio:
        total = count_total_queries(args.query)
        ratio = num_matched / total if total > 0 else 0
        print(f"Matched / Total ratio: {ratio:.4f}")


if __name__ == "__main__":
    main()


# python3 script2.py serratia_search.blastn.tsv 95, for 95% identity 

# answer is: Number of matched query sequences: 4357

#$ python3 script2.py serratia_search.blastn.tsv 95 --ratio
# answer: Number of matched query sequences: 4357
#Matched / Total ratio: 0.8233


