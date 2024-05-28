# cfSPI-pipeline
This repository contains the **cfDNA Single-strand Pathogen Identification (cfSPI) pipeline**, an open-source data analysis Snakemake workflow designed for detecting pathogenic microbial species. The cfSPI pipeline processes paired-end Illumina sequencing data and has been optimized for detecting *Aspergillus* species.
In short, duplicates were removed (using nubeam), after which high-quality sequencing data was generated (using fastp) by default removal of low quality reads and usage of a low complexity filter as well as by adapter removal and removal of short (<35bp) reads (using AdapterRemoval). After subtraction of host sequences by mapping to the human reference genome using Bowtie2, the remaining paired-end reads were taxonomically classified using kraken2.

## Install snakemake in a new/clean conda environment
Install snakemake in your environment (installation suggestion see snakemake website: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

## Create a samplesheet and configfile 
Create a samplesheet and a configfile in folder `configs/` accordingly. 
An example samplesheet can be found at `config/samples.txt` and an example configfil at `config/config_samples.yaml`. Within `config/config_samples.yaml`, ensure to specify the `run_name`, `outdir` and `db_dir`. 

## Running the pipeline by submitting jobs via [slurm](https://slurm.schedmd.com/documentation.html) scheduler:
1. Start a screen session. 
2. Request an interactive node for submitting jobs (long enough for all jobs to finish, e.g. 48 hours), with 16G mem, 2 cores.
3. While running snakemake, simply run command:
`snakemake --configfile ./config/config_samples.yaml  --snakefile workflow/Snakefile_IFD  --profile ./profile/slurm --conda-frontend conda --use-conda `. (Resources defined in each rule will be used. If not defined, default resources defined in the `profile/slurm/config.yaml` will be used. )

More information See [snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) and page [snakemake slurm](https://snakemake.readthedocs.io/en/stable/executing/cluster.html#executing-on-slurm-clusters). 

## Trouble shooting
- Input fasta folder always should contains following files to start with `*R1*.fastq.gz, *R2*.fastq.gz`
- The current pipeline version includes adapter sequence information utilized by the SRSLY Claret Kit, the KAPA Kit, and 384 IDT UMI's. If another library preparation method is employed, kindly update the `config/config_samples.yaml` file and append the adapter index sequences to the `resources/adapter_indexes/` directory."
