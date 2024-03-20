### Databases 

Kraken2 for classification depends on a hash-table database. In **XXX** we sought to identify the most suitable database for the purpose of screening for potential pathogenic *Aspergillus*. The nine hash-table databases we created are available through ***XXX***. Hence, it is crucial to mention that:
- The **RS** also known as **uR.7**; 
- **EPRSc2** also known as **cRE.21**; 
- **EPRSFv46** also known as **uRE.21**; 
- **EPRSFv46DM** also known as **dRE.21**; 
- **EPRSFv64DM** also kown as **dRE.31**; 
- **EPRSFv64** also known as **uRE.31**; 
- **EPRSFv46MCAspDM** also known as **dREM.258**; 
- **EPRSFv64MCAspDM** also known as **dREM.260**; 

For general information on the database compositon we would like to refer to the manuscript itself. Supplementary Table 5 thereby offers details on the construction of each database. Here, we provide detailed information about the content of a Kraken 2 database in the **_inspect.txt* files, obtained via ```kraken2-inspect --db```. Thereby, for example, the ***RS_inspect.txt*** file contains the details of the **RS** database, the ***EPRSc2_inspect.txt*** file contains the details of the **EPRSc2**, and so forth. 

Also we provide the *seqid2taxid.map* files utilized in creating these nine databases. Thereby, for example, we used the ***RS_seqid2taxid.map*** for the construction of the **RS** database, the ***EPRSc2_seqid2taxid.map*** for the construction of the **EPRSc2**, and so forth. 
