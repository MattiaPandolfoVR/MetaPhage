# Coded by Mattia Pandolfo (mattia.pandolfo@univr.it)
shhh <- suppressPackageStartupMessages
shhh(library(DT))
shhh(library(readr))
shhh(library(dplyr))
shhh(library(plyr))
shhh(library(seqinr))
shhh(library(tidyverse))
shhh(library(gtools))

# Reading input
args <- commandArgs(trailingOnly = TRUE)
# check arguments
if (length(args) < 1) {
 stop("Usage: checkv_table.R <quality_summary> <>\n")
}
############################### QUALITY SUMMARY ################################
file_qs <- args[1]
if (! file.exists(file_qs)){
  stop("ERROR: quality summary file not found in: ", file_qs, "\n")
}
qual_sum <- read.delim(file_qs, sep = "\t", check.names = F)
# change characters to factors
qual_sum <- mutate_if(qual_sum, is.character, as.factor)
df <- qual_sum
# rename contigs_id column to vOTU
colnames(df)[1] <- "ViralOTU"
# create the table with DT package
qs <- DT::datatable(df, filter="top", rownames = FALSE,
                      extensions = c('Buttons','Scroller'), escape = F, 
                      options = list(scrollY = "650px",
                                     scrollX = TRUE,
                                     scroller = TRUE,
                                     dom = 'BlPfrtip',
                                     buttons = list('copy',
                                                    'print',
                                                    list(extend = 'collection',
                                                         buttons = c('csv',
                                                                     'excel',
                                                                     'pdf'),
                                                         text = 'Download as:')
                                     )
                      )
)
qs$width = 1500
htmlwidgets::saveWidget(qs, "./checkv_table.html",
                        selfcontained = TRUE, libdir = NULL)
