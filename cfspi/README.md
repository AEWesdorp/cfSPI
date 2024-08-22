# cfSPI-pipeline
This repository contains the **cfDNA Single-strand Pathogen Identification (cfSPI) pipeline**, an open-source data analysis Snakemake workflow designed for detecting pathogenic microbial species. The cfSPI pipeline processes paired-end Illumina sequencing data and has been optimized for detecting *Aspergillus* species.

In short, duplicates were removed (using nubeam), after which high-quality sequencing data was generated (using fastp) by default removal of low quality reads and usage of a low complexity filter as well as by adapter removal and removal of short (<35bp) reads (using AdapterRemoval). After subtraction of host sequences by mapping to the human reference genome using Bowtie2, the remaining paired-end reads were taxonomically classified using kraken2.

## Getting started
To use the cfSPI-pipeline, follow these steps: 

#### Prerequisites
1. Ensure you have `conda` installed on your system.

#### Installation
1. Clone the GitHub Repository:
    ```bash
    git clone https://github.com/AEWesdorp/cfSPI.git
    ```
2. Create and Activate a Conda Environment:
    ```bash
    # Create a new empty environment called "cfspi_env"
    conda create -c conda-forge -c bioconda -n cfspi_env snakemake bowtie2
    # Activate the environment "cfspi_env"
    conda activate cfspi_env
    ```
    
## Create a combined indexed version of the human reference genome 
To mitigate potential false positives arising from incomplete host read subtraction, we implemented a dual-mapping strategy using `bowtie2-2.5.1`. This involved aligning reads to a combined host genome version comprising GRCh38.p14 and CHM13v2.

We obtained the GRCh38.p14 'Genome sequence (FASTA)' from the [NCBI RefSeq](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) and the CHM13v2 genome via this [GitHub repository](https://github.com/marbl/CHM13). To index these genomes, we executed the following commands:
```bash
cat GCF_000001405.40_GRCh38.p14_genomic.fna chm13v2.0.fa > chm13v2.0_PLUS_GCF_000001405.40_GRCh38.p14_genomic.fna
bowtie2-build chm13v2.0_PLUS_GCF_000001405.40_GRCh38.p14_genomic.fna chm13v2.0_PLUS_GCF_000001405.40_GRCh38.p14_genomic
```
This process yielded an indexed dual-genome, **`chm13v2.0_PLUS_GCF_000001405.40_GRCh38.p14_genomic`**.

Note: If you plan to use the `Snakemake_IFD_plus` or are specifically interested in mapping to either `chm13v2.0.fa` or `GCF_000001405.40_GRCh38.p14_genomic.fna` using the default `Snakemake_IFD` pipeline, you can index these reference genomes with the following commands:
```bash
bowtie2-build chm13v2.0.fa chm13v2.0
bowtie2-build GCF_000001405.40_GRCh38.p14_genomic.fna GCF_000001405.40_GRCh38.p14_genomic
```

## Create a samplesheet and configfile 
Create a samplesheet and a configfile in folder `configs/` accordingly. 
An example samplesheet can be found at **`./config/samples.tsv`** and an example configfile at **`./config/config_samples.yaml`**. 

In **`./config/samples.tsv`**, ensure to specify the following: 
- *sample_name*: Specify name of your sample.
- *library_prep*: Specify the library preparation you have used (options: `SRSLY` or `KAPA`)
- *adapter_type*: Specify which adapter were used during Illumina library contstruction (options: `SRSLY_dual_index`, `KAPA_single_index` or `IDT384UMI_dual`)
- *UDI*: Set the UDI of the sample.
- *path_to_R1_R2*: Please provide the absolute directory path where the raw sequencing files are stored (path to `*_R1_*.fastq.gz`, `*_R2_*.fastq.gz`). 
  
In **`./config/config_samples.yaml`**, ensure to specify the following:

General Settings:
- *units*: Specify name of the samplesheet (for example, `./config/samples.tsv`). Relative path recommended.  
- *output_folder*: Please provide a unique name. An output folder with this name will be created within the specified *outdir*.
- *outdir*: Please provide the absolute path to the output directory where the *output_folder* will be generated and the results will be stored.

Reference Genome Settings:
- *reference_genome*: Indicate the reference genome to be used. We recommend using the dual-genome `chm13v2.0_PLUS_GCF_000001405.40_GRCh38.p14_genomic`
- *reference_genome_dir*: Please provide the absolute path to the directory where the reference genome(s) are stored.
  
Kraken2 Classification Settings:
- *database*: Specify the name of the database(s) to be used for kraken2 classification. Each *database* name should correspond to a directory within the *database_dir*. the *database* folder should be containing all necessary database files (e.g., `hash.k2d`, `opts.k2d`, `taxo.k2d`). Please note that our kraken2 database(s) are **not** included in this GitHub repository, but that information on the databases can be found in `./../databases/'; instructions for database construction are provided in the manuscript. Alternatively, you can access standard prebuilt Kraken2 index databases [here](https://benlangmead.github.io/aws-indexes/k2). 
- *database_dir*: Provide the absolute path to the directory where the kraken2 database(s) are located.
- *k2_threshold*: Set the threshold value for Kraken2 classification.

## Running the cfspi-pipeline on an interactive node
1. Start a screen session. 
2. Request an interactive node for for running the jobs (long enough to finish all jobs of one liquid biopsy sample, e.g. 24 hours), with 450G mem, 16 cores. 
3. Move to the `cfspi/` sub-directory within the cloned Git directory where your workflow resides.
4. Activate your conda environment.
      ```bash
    conda activate cfspi_env
    ```
5. Run the snakemake pipeline.
   ```bash
   snakemake --configfile ./config/config_samples.yaml  --snakefile workflow/Snakefile_IFD  --cores all --conda-frontend conda --use-conda
   ```

## Running the cfspi-pipeline by submitting jobs via [slurm](https://slurm.schedmd.com/documentation.html) scheduler:
1. Start a screen session. 
2. Request an interactive node for submitting jobs (long enough for all jobs to finish, e.g. 48 hours), with 16G mem, 2 cores.
3. Move to the `cfspi/` sub-directory within the cloned Git directory where your workflow resides.
4. Activate your conda environment.
    ```bash
    conda activate cfspi_env
    ```
5. Run the snakemake pipeline.
   ```bash
   snakemake --configfile ./config/config_samples.yaml  --snakefile workflow/Snakefile_IFD  --profile ./profile/slurm --conda-frontend conda --use-conda
   ```
   (Resources defined in each rule will be used. If not defined, default resources defined in the `profile/slurm/config.yaml` will be used.)

More information See [snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) and page [snakemake slurm](https://snakemake.readthedocs.io/en/stable/executing/cluster.html#executing-on-slurm-clusters). 

## Trouble shooting
- Input fasta folder always should contains following files to start with `*_R1_*.fastq.gz, *_R2_*.fastq.gz`
- The current pipeline version includes adapter sequence information utilized by the SRSLY Claret Kit, the KAPA Kit, and 384 IDT UMI's. If another library preparation method is employed, kindly update the `workflow/rules/trim.smk` file and append the adapter index sequences to the `resources/adapter_indexes/` directory.
