# Coded by Mattia Pandolfo (mattia.pandolfo@univr.it)
shhh <- suppressPackageStartupMessages
shhh(library(heatmaply))
shhh(library(phyloseq))
shhh(library(metagenomeSeq))

# Reading input
args <- commandArgs(trailingOnly = TRUE)
# check arguments
if (length(args) < 4) {
  stop("Usage: heatmap.R <count_table> <taxonomy_table> <metadata> <heatmap_var>\n")
}
################################## COUNT TABLE #################################
file_count <- args[1]
if (! file.exists(file_count)){
  stop("ERROR: count table file not found in: ", file_count, "\n")
}
count <- read.delim(file_count, row.names=1, sep = "\t", check.names = F)
################################ TAXONOMY TABLE ################################
file_taxo <- args[2]
if (! file.exists(file_taxo)){
  stop("ERROR: taxonomy table file not found in: ", file_taxo, "\n")
}
# change "O", "n.a." and "uncharacterize" to "NA" in taxonomy table
taxo <- read.delim(file_taxo, row.names=1, sep = "\t", check.names = F, na.strings = c("n.a.", "O", "Unclassified"))
################################### METADATA ###################################
file_meta <- args[3]
if (! file.exists(file_meta)){
  stop("ERROR: metadata file not found in: ", file_meta, "\n")
}
metadata <- read.delim(file_meta, row.names=1, sep = ",", check.names = F)
metadata$Sample <- rownames(metadata)
metadata <- metadata %>%
  select(Sample, everything())
############################### HEATMAP VARIABLE ###############################
heatmap_var <- args[4]
if (heatmap_var == FALSE){
  stop("ERROR: heatmap_var parameter not found. You can set it in your project config.file (e.g. heatmap_var = \"metadata_variable1\" \n")
}

# Phyloseq object creation
ps <- phyloseq(otu_table(count, taxa_are_rows = TRUE),
               tax_table(as.matrix(taxo)),
               sample_data(metadata))

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
colnames(ps@tax_table)
# remove unwanted columns
ps@tax_table <- ps@tax_table[, -(2:10)]

# Data-frames creation
df_merged <- merge(x = ps@otu_table, y = ps@tax_table, by = 0 )
rownames(df_merged) = df_merged$Row.names
df_merged$Row.names <- NULL
df_meta <- metadata[metadata$Sample %in% colnames(count),]
# Fix labels in rowside and colside columns
colside <- as.data.frame(df_meta[,c("Sample",heatmap_var)])
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
                scale= 'column',
                show_dendrogram = FALSE,
                showticklabels = c(TRUE,FALSE),
                k_col = NA, k_row = NA,
                row_text_angle = 0,
                column_text_angle = 45,
                label_names = c("vOTU", "Sample", "log_CSS_normalized_counts"),
                plot_method = c("ggplot"),
                margins = c(100,10,50,10),
                hide_colorbar = FALSE) %>%
  plotly::layout(showlegend = FALSE,
                 annotations = list(
                   visible = FALSE))
htmlwidgets::saveWidget(h1, "./heatmap.html", 
                        selfcontained = TRUE, libdir = NULL)