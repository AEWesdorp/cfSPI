#!/bin/bash 
### SIMULATE FUNGAL READS ###
echo 'Seed number: ' $1
seed_nr=$1

mkdir -p ./sim_data/ReSeq_simulation/
while read specie; do
	echo $specie
	sp_dir="$(grep "$specie" ./sim_data/assembly_summary_refseq.txt | cut -d$'\t' -f 20  )"
	sp_nm="$(echo $sp_dir | awk -F/ '{print $NF}')"
	wget "${sp_dir}/${sp_nm}_genomic.fna.gz" -O "./sim_data/genomes_NCBI/$(echo $specie | sed 's/ /_/g' )_genomic.fna.gz" 
	gunzip --force "./sim_data/genomes_NCBI/$(echo $specie | sed 's/ /_/g' )_genomic.fna.gz" 

	for name in $(ls ../ReSeq_statsOnly/ | head -n 1); do #only simulate based on one sample
		mkdir -p ./sim_data/ReSeq_simulation/${name}/
		#for seed_nr in seed; do 
		mkdir -p "./sim_data/ReSeq_simulation/${name}/${name}NCBI$(echo $specie | sed 's/ //g' )Syn${seed_nr}/"

		reseq  illuminaPE -j 32 -R "./sim_data/genomes_NCBI/$(echo $specie | sed 's/ /_/g' )_genomic.fna"  \
		-s ../ReSeq_statsOnly/${name}/${name}.bam.reseq \
		--refBias no --numReads 100000 --seed $seed_nr \
		-1 "./sim_data/ReSeq_simulation/${name}/${name}NCBI$(echo $specie | sed 's/ //g' )Syn${seed_nr}/${name}_NCBI_$(echo $specie | sed 's/ //g' )_seed${seed_nr}.simulatedR1.fastq" \
		-2 "./sim_data/ReSeq_simulation/${name}/${name}NCBI$(echo $specie | sed 's/ //g' )Syn${seed_nr}/${name}_NCBI_$(echo $specie | sed 's/ //g' )_seed${seed_nr}.simulatedR2.fastq"

		gzip -f "./sim_data/ReSeq_simulation/${name}/${name}NCBI$(echo $specie | sed 's/ //g' )Syn${seed_nr}/${name}_NCBI_$(echo $specie | sed 's/ //g' )_seed${seed_nr}.simulatedR1.fastq"
		gzip -f "./sim_data/ReSeq_simulation/${name}/${name}NCBI$(echo $specie | sed 's/ //g' )Syn${seed_nr}/${name}_NCBI_$(echo $specie | sed 's/ //g' )_seed${seed_nr}.simulatedR2.fastq"
		#done
	done
	rm "./sim_data/genomes_NCBI/$(echo $specie | sed 's/ /_/g' )_genomic.fna"
done < ./sim_data/assembly_summary_refseq_ofInterest.txt
