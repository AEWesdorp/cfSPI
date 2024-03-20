#!/bin/bash 
### SELECT FUNGAL NCBI GENOMES FOR SIMULATION ###

mkdir -p ../sim_data/

wget ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt -O ../sim_data/assembly_summary_refseq.txt

grep -e "Penicillium" -e "Aspergillus" -e "Komagataella phaffii" -e "Mucor" -e "Rhizopus" -e "Pneumocystis" -e "Fusarium" -e "Botrytis cinerea" ../sim_data/assembly_summary_refseq.txt | \
grep -v "virus" | cut -f 8 | sort | uniq > ../sim_data/assembly_summary_refseq_ofInterest.txt
