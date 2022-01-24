shhh <- suppressPackageStartupMessages
shhh(library(heatmaply))
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
########################### HEATMAP VARIABLE ###################################
heatmap_var <- args[4]

# Phyloseq object creation
ps0 <- phyloseq(otu_table(count, taxa_are_rows = TRUE),
                tax_table(as.matrix(taxo)),
                sample_data(metadata))

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

# Data-frames creation
df_merged <- merge(x = ps@otu_table, y = ps@tax_table, by = 0 )
rownames(df_merged) = df_merged$Row.names
df_merged$Row.names <- NULL
df_merged$Accession <- NULL
df_merged$Status <- NULL
df_merged$VC <- NULL
df_merged$Level <- NULL
df_merged$Weight <- NULL
df_merged$Host <- NULL
df_merged$BaltimoreGroup <- NULL
df_merged$Realm <- NULL
df_merged$Kingdom <- NULL
df_meta <- metadata[metadata$Sample %in% colnames(count),]
# Fix labels in rowside and colside columns
colside <- as.data.frame(df_meta[,c("Sample",heatmap_var)])
#colside <- colside[order(colside$Study_group),]
colside <- colside[order(colside[,2]),]
colside$Sample <- NULL
# create data_frame to plot
df_heatmap <- df_merged[,c(match(rownames(colside),colnames(df_merged)),
                           (ncol(df_merged)-6):ncol(df_merged))]
# Heatmap creation
h1 <- heatmaply(df_heatmap,
                main = "Per sample vOTUs abundance and Viral taxonomy",
                xlab = "Samples",
                ylab = "vOTUs",
                row_dend_left = FALSE,
                col_side_colors = colside,
                dendrogram = "row",
                show_dendrogram = FALSE,
                scale="row",
                showticklabels = c(TRUE,FALSE),
                k_col = NA, k_row = NA,
                row_text_angle = 0,
                column_text_angle = 45,
                label_names = c("vOTU", "Sample", "Abundance"),
                plot_method = c("ggplot"),
                margins = c(100,10,50,10),
                hide_colorbar = FALSE) %>%
  plotly::layout(showlegend = FALSE,
                 annotations = list(
                   visible = FALSE))
htmlwidgets::saveWidget(h1, "heatmap.html", 
                        selfcontained = TRUE, libdir = NULL)
