shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(plotly))
shhh(library(phyloseq))
shhh(library(metagenomeSeq))

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

################################## FILTERING ###################################
# Filter uncharacterized taxas
ps <- subset_taxa(ps0, !is.na(Phylum) & !Phylum %in% c("", "Unclassified"))
# relative abundance filter
FSr  = transform_sample_counts(ps, function(x) x / sum(x))
FSfr = filter_taxa(FSr, function(x) mean(x) < 0.00005, TRUE)
rmtaxa = taxa_names(FSfr)
alltaxa = taxa_names(ps)
myTaxa = alltaxa[!alltaxa %in% rmtaxa]
physeqaFS <- prune_taxa(myTaxa, ps)
ps = filter_taxa(physeqaFS, function(x) sum(x >= 1) >= (2), TRUE)

################################ NORMALIZATION #################################
# CSS
ps_m <- phyloseq::phyloseq_to_metagenomeSeq(ps)
# normalized count matrix
ps_norm <- metagenomeSeq::MRcounts(ps_m, norm = TRUE, log = TRUE)
# abundance to css-normalized
phyloseq::otu_table(ps) <- phyloseq::otu_table(ps_norm, taxa_are_rows = T)
# Restore sample_data rownames
row.names(ps@sam_data) <- c(1:nrow(ps@sam_data))
# save normalized phyloseq object
saveRDS(ps, "phyloseq_css_norm.rds")