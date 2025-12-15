#!/usr/bin/env python3

BIN_SIZE = 500
INPUT_FILE = "sequences.txt"
OUTPUT_FILE = "lengths_binned.csv"

lengths = []
with open(INPUT_FILE) as f:
    f.readline()  
    for line in f:
        line = line.strip()
        if not line:
            continue
        fields = line.split()
        gdna = fields[1]              
       lengths.append(len(gdna))     

# summary statistics counts 
if len(lengths) == 0:
    print("No genes found in file.")
    raise SystemExit

avg_len = sum(lengths) / len(lengths)
min_len = min(lengths)
max_len = max(lengths)

print("Gene length summary (gDNA):")
print(f"Count:   {len(lengths)}")
print(f"Average: {avg_len:.2f}")
print(f"Shortest:{min_len}")
print(f"Longest: {max_len}")

#bin lengths by 500 bp and write csv file 
max_bin = (max_len // BIN_SIZE) * BIN_SIZE
bin_counts = {b: 0 for b in range(0, max_bin + BIN_SIZE, BIN_SIZE)}

for L in lengths:
    b = (L // BIN_SIZE) * BIN_SIZE
    bin_counts[b] += 1

with open(OUTPUT_FILE, "w") as out:
    out.write("bin_size,count\n")
    for b in sorted(bin_counts):
        out.write(f"{b},{bin_counts[b]}\n")

# run python3 problem3.py


