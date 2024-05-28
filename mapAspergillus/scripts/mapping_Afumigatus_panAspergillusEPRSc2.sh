## TaxID of Aspergillus fumigatus
taxid="746128" #A. fumigatus

for sample in $(cut -f 1  ../../cfspi/config/samples.txt | grep -v sample_name); do 
	echo ${sample}
	./extract_kraken2_taxids_edit.py -r ../../output/cfspi/samples/results/kraken2_report/after_host_mapping/${sample}_EPRSc2_conf0.4.report -t ${taxid}  \
		--include-children --output_file ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_taxIDs.txt

	awk 'NR == FNR { keywords[$1]=1; next; } { if ($3 in keywords) print $2; }'         ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_taxIDs.txt ../../output/cfspi/samples/results/kraken2_output/after_host_mapping/${sample}_EPRSc2_conf0.4.output > ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_nms.txt

	nr_lines=$( wc -l  ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_nms.txt | cut -d" " -f 1 )
	echo $nr_lines

	##mapping for reads classified as A. fumigatus
	if [ $nr_lines != "0" ]; then
		echo "mapping reads classified to this taxa"
		filterbyname.sh in1=../../output/cfspi/samples/results/host_mapping/${sample}_unmapped_host_r1.fq in2=../../output/cfspi/samples/results/host_mapping/${sample}_unmapped_host_r2.fq out1=../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_R1.fastq out2=../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_R2.fastq names=../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_nms.txt include=TRUE overwrite=TRUE

		for AspGen in $(ls ../resources/EPRSc2_Aspergillus/*fna ); do
			echo $AspGen
			sh_AspGen=$(echo $AspGen | sed 's|.*FungiDB-46_||g' | sed 's|_Genome_cleaned_final.fna||g')
			mp_AspGen=$(echo $AspGen | sed 's|.fna$||g')

			bowtie2 -p 8 -x ${mp_AspGen} -1 ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_R1.fastq -2 ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_R2.fastq > ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_${sh_AspGen}.sam;
			samtools view -b ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_${sh_AspGen}.sam -o ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_${sh_AspGen}.bam;
			samtools sort ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_${sh_AspGen}.bam > ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_${sh_AspGen}_srt.bam
			samtools index ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_${sh_AspGen}_srt.bam

			rm ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_${sh_AspGen}.sam
			rm ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_${sh_AspGen}.bam
		done

		rm ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_R1.fastq 
		rm ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_R2.fastq 
	fi

	rm ../../output/mapAspergillus/${sample}_EPRSc2_conf0.4_${taxid}_nms.txt
done
