#!/usr/bin/env python3

import sys
def find_lowID_best_hits(blastfile, threshold=80.0):
    best_hits = {}
    with open(blastfile) as f:
        for line in f:
            fields = line.strip().split("\t")
            query_id = fields[0]
            pident = float(fields[2])

            # keep best (highest identity) hit per query
            if query_id not in best_hits or pident > best_hits[query_id]:
                best_hits[query_id] = pident

    # queries whose best hit is below threshold
    low_id_hits = []
    for q, pid in best_hits.items():
        if pid < threshold:
            low_id_hits.append((q, pid))

    return low_id_hits

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 script3.py <blastfile> [threshold]")
        sys.exit(1)

    blastfile = sys.argv[1]
    threshold = float(sys.argv[2]) if len(sys.argv) > 2 else 80.0

    results = find_lowID_best_hits(blastfile, threshold)

    print("QueryID\tBestPercentIdentity")
    for q, pid in results:
        print(f"{q}\t{pid:.2f}")


if __name__ == "__main__":
    main()

# run 'python3 script3.py serratia_search.blastn.tsv 80 > low_identity_hits.tsv'
#saved output in file low_identity_hits.tsv
# total lines  in low_identity_hits.tsv: 15 
