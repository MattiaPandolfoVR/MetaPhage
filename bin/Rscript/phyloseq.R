shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(plotly))
shhh(library(phyloseq))
shhh(library(metagenomeSeq))
shhh(source("bin/Rscript/filter&CSSnormalize.R"))

# Reading input
args <- commandArgs(trailingOnly = TRUE)
################################### COUNT TABLE ################################
file_count <- args[1]
count <- read.csv(file_count, sep = '\t', row.names = 1)
################################# TAXONOMY TABLE ###############################
file_taxo <- args[2]
taxo <- read.csv(file_taxo, sep = '\t', row.names = 1)
################################### METADATA ###################################
metadata <- args[3]                                                             
metadata <- read.csv(metadata, sep = ',', row.names = 1)
metadata$Sample <- rownames(metadata)

# Phyloseq object creation
ps0 <- phyloseq(otu_table(count, taxa_are_rows = TRUE),
                tax_table(as.matrix(taxo)),
                sample_data(metadata))
# save phyloseq object
saveRDS(ps0, "phyloseq.rds")

####################### FILTER & CSS NORMALIZE #################################
ps <- filter_CSSnormalize(ps0)

# save filtered and CSS normalized phyloseq object
saveRDS(ps, "phyloseq_filt_css_norm.rds")