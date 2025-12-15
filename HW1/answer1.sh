#!/bin/bash -l

#dataset 

cd ~/bigdata/Kinza/GEN220/GEN220_data/tabular
ll # list the files 

#copy files to HW1 directory 
cp threatened-species.csv.gz ~/bigdata/Kinza/GEN220/HW1/

# size of threatened-species.csv.gz

ls -l 

# 3.6M, size of threatened-species.csv.gz

#number of lines in the file
wc -l threatened-species.csv.gz
#13577 threatened-species.csv.gz

zcat threatened-species.csv.gz | head -n 1
#taxonid,kingdom_name,phylum_name,class_name,order_name,family_name,genus_name,scientific_name,taxonomic_authority,infra_rank,infra_name,population,category,main_common_name

# count unique phyla and order name 

zcat threatened-species.csv.gz | cut -d',' -f2 | tail -n +2 | sort | uniq | wc -l
#answer is: 4 
zcat threatened-species.csv.gz | cut -d',' -f3 | tail -n +2 | sort | uniq | wc -l
#answer is: 16 

# how many kingdom FUNGI are present in the file?
zcat threatened-species.csv.gz | cut -d',' -f2 | tail -n +2 | grep -w FUNGI | wc -l
#answer is: 627 


