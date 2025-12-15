#!/bin/bash -l

#download genomes 

#Thermus aquaticus genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/399/775/GCF_001399775.1_ASM139977v1/GCF_001399775.1_ASM139977v1_genomic.fna.gz
#Bacillus subtilis genome 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz

#unzip file 
gunzip *.gz 


