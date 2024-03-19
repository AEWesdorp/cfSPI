#!/bin/bash 
### SELECT FUNGAL NCBI GENOMES FOR SIMULATION ###

mkdir -p ../sim_data/
mkdir -p ../sim_data/genomes_NCBI/

wget ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt -O ../sim_data/assembly_summary_genbank.txt
wget ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt -O ../sim_data/assembly_summary_refseq.txt
wget ftp.ncbi.nlm.nih.gov/genomes/genbank/sim_data/assembly_summary.txt -O ../sim_data/assembly_summary.txt

grep -e "Penicillium" -e "Aspergillus" -e "Komagataella phaffii" -e "Mucor" -e "Rhizopus" -e "Pneumocystis" -e "Fusarium" -e "Botrytis cinerea" ../sim_data/assembly_summary_refseq.txt | \
grep -v "virus" | cut -f 8 | sort | \
uniq > ../sim_data/assembly_summary_refseq_ofInterest.txt

grep -e "Penicillium" -e "Aspergillus" -e "Komagataella phaffii" -e "Mucor" -e "Rhizopus" -e "Pneumocystis" -e "Fusarium" -e "Botrytis cinerea" ../sim_data/assembly_summary_refseq.txt | \
grep -v "virus" | cut -f 7-8 | sort | \
uniq | sed 's| ||g' | perl -ane 'print "$F[1]\t$F[0]\n"' > ../NCBI_genome2taxid.map
