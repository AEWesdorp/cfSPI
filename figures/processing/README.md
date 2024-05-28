# Preprocessing for figure making 

In this directory, we conduct preprocessing, such as handling the output of the cfSPI, to prepare the data for figure creation.

- `DB_stats.ipynb`: This notebook provides statistics on the kraken2 hash-table databases. The resulting Rdata files can be found here `../../output/preprocessing/DBstats/`
- `samples.DBs_thresholds.ipynb`: This notebook processes the sample output of the cfSPI pipeline. The resulting Rdata files can be found here `../../output/preprocessing/samples/`
- `simulations.DBs_thresholds.ipynb`: This notebook processes the simulation output of the cfSPI pipeline. The resulting Rdata files can be found here `../../output/preprocessing/simulations/`
- `calculations_LOSD_species.ipynb`: In this notebook, we compute the LoSD (Limit of Detection) in MPM (Mapped Per Million) at the species level. The resulting Rdata files can be found here `../../output/preprocessing/LoSD/`
- `calculations_LOSD_genus.ipynb`: This notebook computes the LoSD in MPM at the genus level. The resulting Rdata files can be found here `../../output/preprocessing/LoSD/`
- `calculations.Proof_of_principle.ipynb`: Here, we calculate the mean Fisher's tests. The resulting Rdata files can be found here `../../output/preprocessing/samples/`

#### Required
- R 4.2.0
- dplyr 1.1.2
- khroma 1.11.0 
- ggplot2 3.4.4
- ggpattern 1.0.1
- reshape2 1.4.4
- tidyverse 2.0.0
- stringr 1.5.0
- patchwork 1.1.2
- Hmisc 5.1
- rstatix 0.7.2
- readr 2.1.4
- Rsamtools 2.14.0
- seqinr 4.2-27
- ggpubr 0.6.0
- ggbreak 0.1.1
- gghighlight 0.4.0
- ggnewscale 0.4.9
- ggh4x 0.2.6 
