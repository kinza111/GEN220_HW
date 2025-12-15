#!/usr/bin/env python3
import sys

for filename in sys.argv[1:]:
    gc = 0
    total = 0

    with open(filename) as f:
        for line in f:
            if line.startswith(">"):
                continue
            line = line.strip().upper()
            gc += line.count("G") + line.count("C")
            total += len(line)

    gc_percent = (gc / total) * 100
    print(f"{filename}\t{gc_percent:.2f}%")

#answer 
#in the command line run 'python3 problem1.py GCF_000009045.1_ASM904v1_genomic.fna GCF_001399775.1_ASM139977v1_genomic.fna'
#Result: 

#GCF_000009045.1_ASM904v1_genomic.fna	43.51%
#GCF_001399775.1_ASM139977v1_genomic.fna	68.04%

