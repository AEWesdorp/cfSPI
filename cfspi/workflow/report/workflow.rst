This workflow performs pathogen detection on (single stranded enriched) cell-free DNA data.
Steps:
1. First, quality check was applied on to illumina sequencing fastq data with fastp low complexity low_complexity_filter
2. Adapter removal is then applied with AdapterRemoval.
3. Then, create fastqc report to make sure bases are trimmed properly.
4. Afterward, reads are mapped to host genome.
5. Take the unmapped reads, use kraken2 to analyse pathogen composition. 
--> Directory of the viral db locally: 
6. Calculate ratio of top 2 species identified to the total amount of unmapped reads (for both singleton and paired). Calculate the domain level ratio. Note the level of "unclassified reads"
7. download taxon_id.fasta of the top 2 species. build as reference.
--> https://github.com/npbhavya/Kraken2-output-manipulation

8. extract reads mapped to each species, map them to pathogen genome reference. & Dedup PCR products to count unique mappings.
9. Calculate ratio of top 2 species identified to the total amount of unmapped (unique!) reads (for both singleton and paired).
10. 




7. Generate report (containing results from 6.)
