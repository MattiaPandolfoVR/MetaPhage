library(VirFinder)
require(Biostrings)
require(parallel)
library(seqinr)

# Take the input passed by Nextflow, which are
# ${scaffold}
# ${task.cpus}
# $workflow.projectDir
args <- commandArgs(trailingOnly = TRUE) 
fileinFASTA <- args[1] 
taskcpus <- args[2]
workdir <- args[3]

# Check if the parVF_pred script is in the bin folder
if(!exists("parVF_pred", mode="function")) source(file.path(workdir, "bin/Rscript/parVF_pred.R"))

# In predResult, store parVF_pred results, runned with ${scaffold} and ${task.cpus} values from NF
predResult <- parVF.pred(fileinFASTA, taskcpus)

# Subset the result for pvalues < 0.005
subsetted <- subset(predResult, pvalue < 0.005)
options(width=10000) 

# Order the result by increasing pvalues (smaller first)
predResult <- predResult[order(predResult$pvalue),]

# Order the subset by increasing pvalues (smaller first)
subsetted <- subsetted[order(subsetted$pvalue),]

# Fix rownames of predResult to not repeat
rownames(predResult) <- 1:nrow(predResult)

# Store the predResult output in results.txt file 
sink("results.txt")
print(predResult)
sink()

# Create a df of the subset output
dfFASTA <- data.frame(subsetted)
# Fix name problem creating empy fasta files
dfFASTA$name <- sub(" .*", "", dfFASTA$name)

seqFASTA <- read.fasta(file = fileinFASTA, forceDNAtolower = FALSE, 
                       seqtype = "DNA", as.string = TRUE, set.attributes = FALSE)

sub <- seqFASTA[names(seqFASTA) %in% dfFASTA$name]
names(sub) = sub("^", "virfinder_", names(sub))
write.fasta(sequences = sub, names = names(sub), as.string=TRUE,
            nbchar = 60, file.out = "viral_sequences.fasta")