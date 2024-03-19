### cfDNA read simulations

To assess the accuracy of Aspergillus classification and determine the Limit of Significant detection (LoSD), we simulated Illumina reads resembling cell-free DNA (cfDNA) from fungal genomes sourced from the NCBI RefSeq. This process generated a total of 87 simulated datasets, comprising 54 Aspergillus genomes, 7 Penicillium genomes, and 25 genomes from other pathogenic fungi. For specifics regarding the number of simulated reads per dataset, please refer to Supplementary Table 2. 

Instructions and code for executing these simulations are outlined below. Simulated files can be found here: **XXX**

Generate tab separated *input_sequencing_files.txt* file. 
First column contains path to folder containing raw paired-end sequencing data. Raw files should be named "*R1*.fastq.gz" and "*R2*.fastq.gz"
Second column contains sample_name (e.g. sample1)
Third column containd UDI (e.g. UDI01)

```bash
##prepare for simulations
PATH_HumanGenome=/home/user/human_genome #path to human genome (without .fna or .fa); make sure the human genome is indexed using bowtie2-build
SEED_NR=1 #set seed; for example 1 

##start ReSeq simulation
cd ./scripts

#obtain stats ReSeq
sh ./ReSeq_statsOnly.sh $PATH_HumanGenome

#select NCBI refseq genomes for simulation; modify if needed - "-e taxid name"
sh ./genomeSelection_NCBI.sh

#simulate microbial reads from NCBI refseq genomes
sh ./ReSeq_simulation_NCBI.sh $SEED_NR

#generate units files for Snakemake pipeline
sh ./ReSeq_simulation_NCBI_createUnits.sh $PATH_HumanGenome $SEED_NR
```
