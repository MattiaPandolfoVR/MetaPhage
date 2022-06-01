shhh <- suppressPackageStartupMessages
shhh(library(plotly))
shhh(library(matrixStats))
shhh(library(dplyr))
shhh(library(readr))
shhh(library(phyloseq))
shhh(library(metagenomeSeq))

# Reading input
args <- commandArgs(trailingOnly = TRUE)
############################### COUNT TABLE ####################################
file_count <- args[1]
count <- read.csv(file_count, sep = '\t', row.names = 1)
############################## TAXONOMY TABLE ##################################
file_taxo <- args[2]
taxo <- read.csv(file_taxo, sep = '\t', row.names = 1)
################################ METADATA ######################################
file_meta <- args[3]
metadata <- read.csv(file_meta, sep = ',', row.names = 1)
metadata$Sample <- rownames(metadata)
violin_var = args[4]

# Phyloseq object creation
ps0 <- phyloseq(otu_table(count, taxa_are_rows = TRUE),
                tax_table(as.matrix(taxo)),
                sample_data(metadata))
ps_r <- ps0
# Phyloseq object creation
ps <- ps0

####################### FILTER & CSS NORMALIZE #################################
################################## FILTERING ###################################
# Filter uncharacterized taxas
ps_filter <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "Unclassified"))
# relative abundance filter
FSr  = transform_sample_counts(ps_filter, function(x) x / sum(x))
FSfr = filter_taxa(FSr, function(x) mean(x) > 0.0001, TRUE)
keeptaxa = taxa_names(FSfr)
alltaxa = taxa_names(ps_filter)
myTaxa = alltaxa[alltaxa %in% keeptaxa]
physeqaFS <- prune_taxa(myTaxa, ps_filter)
ps_filter = filter_taxa(physeqaFS, function(x) sum(x >= 1) >= (2), TRUE)

################################ NORMALIZATION #################################
# CSS
ps_m <- phyloseq::phyloseq_to_metagenomeSeq(ps_filter)
# normalized count matrix
ps_norm <- metagenomeSeq::MRcounts(ps_m, norm = TRUE, log = TRUE)
# abundance to css-normalized
phyloseq::otu_table(ps_filter) <- phyloseq::otu_table(ps_norm, taxa_are_rows = T)
# Restore sample_data rownames
row.names(ps_filter@sam_data) <- c(1:nrow(ps_filter@sam_data))
ps <- ps_filter

################################## PLOTS #######################################
# Metadata processing
df_new_meta = data.frame(Sample = metadata$Sample, metadata[[violin_var]])
colnames(df_new_meta)[2] <- paste(violin_var)
# Dataframes creation
df_tot = data.frame()
for(sample in colnames(ps@otu_table[1,1:ncol(ps@otu_table)])){
  df_new = data.frame(sample, rownames(ps@otu_table), ps@otu_table[,sample])
  colnames(df_new) <- c("Sample","viralOTU","Abundance")
  df_tot <- rbind(df_tot, df_new)
}
rownames(df_tot) <- 1:nrow(df_tot)
df_tot$temp <- df_new_meta[,2][match(df_tot$Sample, df_new_meta$Sample)]
colnames(df_tot)[4] <- paste(violin_var)
df_tot <- df_tot[order(df_tot$Sample),]
ps_rich <- estimate_richness(ps_r)
df_rich <- data.frame(Sample = as.character(rownames(ps_rich)),
                      Richness = ps_rich$Observed)
df_rich <- df_rich[order(df_rich$Sample),]
df_tot$Richness <- df_rich[,2][match(df_tot$Sample, df_rich$Sample)]

v1 <- df_tot %>%
  plot_ly(
    x = ~get(violin_var),
    y = ~Richness,
    split = ~get(violin_var),
    type = "violin", bandwidth = 50,
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    )
  )
htmlwidgets::saveWidget(v1, "violin_plot_richness.html",
                        selfcontained = TRUE, libdir = NULL)
