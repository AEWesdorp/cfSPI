# cfSPI

cfSPI is a cfDNA Single-strand Pathogen Identification workflow, designed for paired-end Illumina sequencing data. This repository is associated with the publication: "NGS-based *Aspergillus* detection in plasma and lung lavage of children with invasive pulmonary aspergillosis."

For details on the pipeline, see: [cfSPI](https://github.com/AEWesdorp/cfSPI/tree/main/cfspi).

## Table of Contents
1. [Introduction](#introduction)
2. [Repository Content](#repository-content)
3. [Installation](#installation)
4. [Usage](#usage)
5. [License](#license)
6. [Contact](#contact)

## Introduction

cfSPI is a cfDNA Single-strand Pathogen Identification workflow, designed for paired-end Illumina sequencing data. This repository was created by Emmy Wesdorp from [De Ridder lab](https://www.deridderlab.nl/) at the Center of Molecular Medicine, University Medical Center Utrecht, the Netherlands.

You can find analysis pipelines and scripts for generating figures related to the paper: "NGS-based *Aspergillus* detection in plasma and lung lavage of children with invasive pulmonary aspergillosis."

## Repository Content
The various parts of the analyses are organized into different folders within the main directory. Each folder contains a ` README.md`  file with specific details relevant to the analyses conducted/information provided within that folder.

### cfSPI-pipeline
Details on the pipeline, which processes paired-end Illumina sequencing data to identify pathogenic species, optimized for detecting *Aspergillus* species: [cfSPI](https://github.com/AEWesdorp/cfSPI/tree/main/cfspi).

### databases
General information on the nine hash-table databases created for this project: [databases](https://github.com/AEWesdorp/cfSPI/tree/main/databases).

### figures
Details on data processing and figure generation: [figures](https://github.com/AEWesdorp/cfSPI/tree/main/figures).

### mapAspergillus
Details on the alignment of *A. fumigatus* classified reads to diverse *Aspergillus* genomes for quality assurance: [mapAspergillus](https://github.com/AEWesdorp/cfSPI/tree/main/mapAspergillus).

### simulations
Details on simulating Illumina reads resembling cell-free DNA (cfDNA) from fungal genomes sourced from the NCBI RefSeq: [simulations](https://github.com/AEWesdorp/cfSPI/tree/main/simulations).

## Installation

To install cfSPI, follow these steps:

1. Clone the repository:
    ```bash
    git clone https://github.com/AEWesdorp/cfSPI.git
    cd cfSPI
    ```

2. Install `snakemake` and `R` in your environment.
3. Further dependencies are described in each folder in this repository separately.

## Usage
A `README.md` file could be found in each folder concerning relevant analyses.

## License
This project is licensed under the GNU GENERAL PUBLIC LICENSE. See the LICENSE file for more details.

## Contact
Please contact the authors and create an issue on github to get help.
