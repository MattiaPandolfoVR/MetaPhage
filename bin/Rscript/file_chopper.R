shhh <- suppressPackageStartupMessages
shhh(library(seqinr))
shhh(library(tidyverse))

# Reading input
args <- commandArgs(trailingOnly = TRUE)
############################# vOTU CONSENSUS FILE ##############################
nucleotides <- args[1]
fna <- read.fasta(nucleotides, 
                  set.attributes = FALSE, 
                  as.string=TRUE, forceDNAtolower = FALSE)
############################# vOTU PROTEINS FILE ###############################
aminoacids <- args[2]
faa <- read.fasta(aminoacids, 
                  set.attributes = FALSE, 
                  as.string=TRUE, seqtype = "AA")
############################# vOTU COORDS FILE #################################
coords <- args[3]
gff <- read.delim(coords, 
                  header=F, comment.char="#")

# for each vOTU in the nucleotide multifasta, create a single fasta file.
for(i in 1:length(fna)){
  write.fasta(sequences = fna[[i]], 
              names = names(fna)[i], 
              file.out = paste0(names(fna)[i], ".fasta"), 
              open = "w")
}
# for each vOTU, create a .faa file containing the relative protein sequences
faa_groups <- split(x = faa, 
                    f = str_extract(names(faa), "vOTU_[0-9]+"), 
                    drop = FALSE)
for(i in 1:length(faa_groups)){
  mapply(names(faa_groups[[i]]), 
         unlist(faa_groups[[i]]), 
         FUN = function(name, sequence){
           write.fasta(sequences = sequence, 
                       names = name, 
                       file.out = paste0(substr(names(faa_groups[i]), 1,
                                                nchar(names(faa_groups[i]))),
                                         ".faa"), open = "a")
           })
}

# for each vOTU, create a .gff file containing the relative coordinates
gff_groups <- unique(gff$V1)
for(i in 1:length(gff_groups)){
  df_gff <- subset(gff, V1 == gff_groups[i])
  writeLines("##gff-version  3", paste0(gff_groups[i], ".gff"))
  write.table(df_gff, paste0(gff_groups[i], ".gff"), 
              quote = FALSE, sep = "\t", col.names=FALSE, 
              row.names=FALSE, append=TRUE)
}