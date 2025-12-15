# HW4 — AI-assisted GFF gene architecture summary

## Prompt used (AI coding prompt)
I want a Python tool that reads GFF3/GTF annotation files from multiple smallgenomes and reports genome annotation summary statistics:
- exon size distribution (min/mean/max)
- intron size distribution (min/mean/max) computed as gaps between exons within the same transcript/gene when possible
- average number of exons per gene
- intergenic distance distribution between neighboring genes on the same contig/chromosome

The twist added is it  also compute gene density (genes per Mb) and an approximate coding fraction (total CDS bp / genome span bp per contig, using annotation coordinates). It outputs summary table in tsv format. 
## project taask/end goal 
Compare gene structure between a bacterium (few/no introns) and a eukaryote (introns/exons) using GFF annotations, and quantify how gene density and coding region fiffer.  

toy data: 1. Thermus aquaticus, 2. Bacillus subtilis, Saccharomyces cerevisiae

#download dataset in toy_data directory 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/399/775/GCF_001399775.1_ASM139977v1/GCF_001399775.1_ASM139977v1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz
gunzip *.gz 

mv GCF_001399775.1_ASM139977v1_genomic.gff therm_aquaticus.gf
mv GCF_000009045.1_ASM904v1_genomic.gff bacillus_sub.gff
mv GCF_000146045.2_R64_genomic.gff scc_cervisiae.gf

## run python script to perform task 'gff_stats.py'

python3 gff_stats.py --gff toy_data/therm_aquaticus.gff --out thermus_summary.tsv
python3 gff_stats.py --gff toy_data/bacillus_sub.gff --out bacillus_summary.tsv
python3 gff_stats.py --gff toy_data/scc_cervisise.gff --out yeast_summary.tsv

# make one final file to compare and analyze gene densities and coding regions 
head -n 1 thermus_summary.tsv > combined_summary.tsv
tail -n 1 thermus_summary.tsv >> combined_summary.tsv
tail -n 1 bacillus_summary.tsv >> combined_summary.tsv
tail -n 1 yeast_summary.tsv >> combined_summary.tsv






