# conda install r-virfinder==1.1 -c bioconda
# conda install r-seqinr==3.1_3 -c bioconda
# Rscript VS_extraction.R sint_A_reads_metaspades_scaffolds.fasta

library(VirFinder)
 
args <- commandArgs(trailingOnly = TRUE) 
filein <- args[1] 
            
predResult <- VF.pred(filein) 
subsetted <- subset(predResult, pvalue < 0.005)

options(width=10000) 
predResult <- predResult[order(predResult$pvalue),]
subsetted <- subsetted[order(subsetted$pvalue),]

sink("results.txt")
print(predResult)
sink()

write.csv(subsetted,"virus.csv")