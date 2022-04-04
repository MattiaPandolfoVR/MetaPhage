# MetaPhage dependencies

To ensure the reproducibility of results, a conda environment YAML
file is provided to:

1. Generate the conda environment
2. Generate a Docker/Singularity container based on the conda environment


The environment contains all the tools used in the pipeline and Nextflow itself:

 
```mermaid
graph TD;
 style input fill:#ff9,stroke:#333,stroke-width:2px
 classDef miner fill:#f99,stroke:#333,stroke-width:2px
 input(INPUT_READS) --> FASTP;
 FASTP --> SEQ_SCREEN;
 SEQ_SCREEN --> KRAKEN2;
 KRAKEN2 --> KRONA_plot;
 SEQ_SCREEN --> ASSEMBLY;
 ASSEMBLY --> QUAST;
 ASSEMBLY --> VirFinder:::miner;
 ASSEMBLY --> VirSorter:::miner;
 ASSEMBLY --> VIBRANT:::miner;
 ASSEMBLY --> Phigaro:::miner;
 VirFinder --> CD-HIT;
 VirSorter --> CD-HIT ;
 VIBRANT --> CD-HIT ;
 Phigaro --> CD-HIT ;
 CD-HIT --> BACKMAPPING;
 BACKMAPPING --> BAMTOCOUNTS;
 CD-HIT --> PRODIGAL;
 PRODIGAL --> vConTACT2;
 vConTACT2 --> GraphAnalyzer;
 GraphAnalyzer --> R_Diversity_Scripts;
 BAMTOCOUNTS --> R_Diversity_Scripts;
 FASTP --> REPORT;
 QUAST --> REPORT;
 R_Diversity_Scripts --> REPORT; 
```
