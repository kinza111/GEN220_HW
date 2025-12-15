#!/bin/bash 

#reference genome 
wget https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_030291735.1/

#cds file 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/291/735/GCF_030291735.1_ASM3029173v1/GCF_030291735.1_ASM3029173v1_cds_from_genomic.fna.gz
mv GCF_030291735.1_ASM3029173v1_cds_from_genomic.fna.gz Ref_genome.fna.gz

#query genome 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/050/750/665/GCF_050750665.1_Serratia_CH31/GCF_050750665.1_Serratia_CH31_cds_from_genomic.fna.gz
mv GCF_050750665.1_Serratia_CH31_cds_from_genomic.fna.gz query_genome.fna.gz

gunzip *.gz 

#blast 

Ref="Ref_genome.fna"
Query="query_genome.fna"

database="Ref_genome_db"
Out_tsv="serratia_search.blastn.tsv"

#make blast database 
makeblastdb -in "$Ref" -dbtype nucl -out "$database"

blastn -query "$Query" -db "$database" \
  -out "$Out_tsv" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

echo "Done. Wrote: $Out_tsv"











