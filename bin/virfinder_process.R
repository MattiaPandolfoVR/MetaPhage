# conda install r-virfinder==1.1 -c bioconda
# conda install r-seqinr==3.1_3 -c bioconda
# Rscript VS_analysis.R sint_A_reads_metaspades_scaffolds.fasta virus.CSV

library(seqinr)
 
args <- commandArgs(trailingOnly = TRUE) 
fileinFASTA <- args[1] 
fileinCSV <- args[2]
            
dfFASTA <- read.csv(file = fileinCSV, header=TRUE)

seqFASTA <- read.fasta(file = fileinFASTA, forceDNAtolower = FALSE, 
            seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)

sub <- seqFASTA[names(seqFASTA) %in% dfFASTA$name]
write.fasta(sequences = sub, names = names(sub), as.string=TRUE,
            nbchar = 60, file.out = "viral_sequences.fasta")
