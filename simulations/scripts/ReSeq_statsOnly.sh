#!/bin/bash 
echo 'Human genome: ' $1 

HumanGenome=$1

###read generation using ReSeq 
while read -r line
do 
	echo "$line"
	file_dir=$(echo $line | awk '{print $1}')
	echo $file_dir
	name=$(echo $line | awk '{print $2}')
	echo $name
	
	### GET RESEQ STATS ###
 	mkdir -p ../sim_data/
	mkdir -p ../sim_data/ReSeq_statsOnly/
	mkdir -p ../sim_data/ReSeq_statsOnly/${name}/

	#combine sequencing data (lanes)
	pigz -cd -p 32 ${file_dir}*R1*.fastq.gz > ../sim_data/ReSeq_statsOnly/${name}/${name}_R1.fastq
	pigz -cd -p 32 ${file_dir}*R2*.fastq.gz > ../sim_data/ReSeq_statsOnly/${name}/${name}_R2.fastq

	#mapping human reads 
	bowtie2 -p 32 -X 2000 -x $HumanGenome \
	-1 ../sim_data/ReSeq_statsOnly/${name}/${name}_R1.fastq \
	-2 ../sim_data/ReSeq_statsOnly/${name}/${name}_R2.fastq \
	| samtools sort -m 10G -@ 4 -T ../sim_data/ReSeq_statsOnly/${name}/${name}.tmp \
	-o ../sim_data/ReSeq_statsOnly/${name}/${name}.bam -

	#cleanup 
	#rm ../sim_data/ReSeq_statsOnly/${name}/${name}_R1.fastq
	#rm ../sim_data/ReSeq_statsOnly/${name}/${name}_R2.fastq

	#get stats ReSeq
	reseq illuminaPE -j 32 -r "${HumanGenome}.fna" \
	-b ../sim_data/ReSeq_statsOnly/${name}/${name}.bam --statsOnly --noTiles
	
done < ./sim_data/input_sequencing_files.txt
