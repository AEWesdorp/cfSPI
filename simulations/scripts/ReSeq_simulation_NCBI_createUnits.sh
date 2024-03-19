#!/bin/bash
echo 'Human genome: ' $1
echo 'Seed number: ' $2

HumanGenome=$1
seed_nr=$2

mkdir -p ../config/
CONFIG=../config/units_ReSeq_simulation_NCBI.txt
head -n 1 ./head_units.txt > $CONFIG

while read specie; do
	for name in $(ls ../ReSeq_statsOnly/ | head -n 1); #does only generate config for one sample
	do
		UDI=$(awk -v nm="$name" '{if ($2 == nm ) print $3 }' ../input_sequencing_files.txt )
		FASTQ_DIR2=./synthetic_datasets/ReSeq_simulation/${name}/${name}NCBI$(echo $specie | sed 's/ //g' )Syn${seed_nr}/
		printf "${name}$NCBI$(echo $specie | sed 's| ||g' |  sed 's|[.]||g' | sed 's|[-]||g' )Syn${seed_nr}\t${name}NCBI$(echo $specie | sed 's/ //g' )Syn${seed_nr}\tnone\tpatient\tplasma\tsynthetic\tSRSLY\tSRSLY_dual_index\t${UDI}\t\t\tRunX\t${HumanGenome}\t${FASTQ_DIR2}\n" >> $CONFIG
	done
done < ../sim_data/assembly_summary_refseq_ofInterest.txt
