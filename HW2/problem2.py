#!/usr/bin/env python3

genes = []

with open("sequences.txt") as f:
    f.readline() 
    for line in f:
        fields = line.strip().split()
        gene_id = fields[0]
        gdna = fields[1]
        mrna = fields[2]

        genes.append((gene_id, len(gdna), len(mrna)))

# sort by gDNA length
genes.sort(key=lambda x: x[1], reverse=True)

for g in genes:
    print(g[0], g[1], g[2])

#ran python3 problem2.py


