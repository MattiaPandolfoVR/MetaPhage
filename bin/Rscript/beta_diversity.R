shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(data.table))
shhh(library(plotly))
shhh(library(phyloseq))
shhh(library(metagenomeSeq))

# Reading input
args <- commandArgs(trailingOnly = TRUE)
# check arguments
if (length(args) < 4) {
  stop("Usage: summary_report.R <count_table> <taxonomy_table> <metadata> <beta_var>\n")
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
################################# BETA VARIABLE ################################
beta_var <- "Study_group"
  #args[4]
if (beta_var == FALSE){
  stop("ERROR: beta_var parameter not found. You can set it in your project config.file (e.g. beta_var = \"metadata_variable1\" \n")
}

# Phyloseq object creation
ps <- phyloseq(otu_table(count, taxa_are_rows = TRUE),
               tax_table(as.matrix(taxo)),
               sample_data(metadata))

####################### FILTER & CSS NORMALIZE #################################
################################## FILTERING ###################################
# Filter uncharacterized taxas
ps_filter <- subset_taxa(ps, !is.na(Phylum))
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

################################# 2D PLOTS #####################################
# Beta-Diversity measurments
# Bray-Curtis
set.seed(12345)
b_ord = ordinate(ps, method = "PCoA", "bray")
b_ord_df <- as.data.frame(b_ord$vectors)
# Jaccard
set.seed(12345)
j_ord = ordinate(ps, method = "PCoA", "jaccard")
j_ord_df <- as.data.frame(j_ord$vectors)
# Plots axis labels
x <- list(
  title = "PCoA 1"
)
y <- list(
  title = "PCoA 2"
)
# check if data are suitable for 2d plots
if (ncol(b_ord_df)<2 || ncol(j_ord_df) <2){
  message ("Sample numerosity in your data doesn't hallow 2D plots creation!")
} else {
  # Bray-Curtis
  # create a column for the user metadata of interest and populate it
  b_ord_df[,beta_var] <- sample_data(ps) %>%
    data.frame() %>%
    select(all_of(beta_var)) %>%
    mutate_if(is.factor,as.character)
  # 2D plot
  p_2d_b = plot_ly(b_ord_df, x= ~Axis.1, y= ~Axis.2,
                   color = ~b_ord_df[,beta_var], colors = "Dark2",
                   mode = "markers",type = "scatter",
                   text  = row.names(b_ord_df),
                   alpha = 0.8, size = 800) %>% layout(xaxis = x, yaxis = y)
  htmlwidgets::saveWidget(p_2d_b, "./beta2D_bray.html",
                          selfcontained = TRUE, libdir = NULL)
  # Jaccard
  # create a column for the user metadata of interest and populate it
  j_ord_df[,beta_var] <- sample_data(ps) %>%
    data.frame() %>%
    select(all_of(beta_var)) %>%
    mutate_if(is.factor,as.character)
  # 2D plot
  p_2d_j = plot_ly(j_ord_df, x= ~Axis.1, y= ~Axis.2,
                   color = ~j_ord_df[,beta_var], colors = "Dark2",
                   mode = "markers",type = "scatter",
                   text  = row.names(j_ord_df),
                   alpha = 0.8, size = 800) %>% layout(xaxis = x, yaxis = y)
  htmlwidgets::saveWidget(p_2d_j, "./beta2D_jaccard.html",
                          selfcontained = TRUE, libdir = NULL)
}
################################# 3D PLOTS #####################################
# check if data are suitable for 3d plots
if (ncol(b_ord_df)<3 || ncol(j_ord_df)<3){
  message ("Sample numerosity in your data doesn't hallow 3D plots creation!")
} else {
  # Bray-Curtis
  p_3d_b <- plot_ly(b_ord_df, x = ~Axis.1, y = ~Axis.2, z = ~Axis.3,
                    color = ~b_ord_df[,beta_var], colors = "Dark2",
                    text = row.names(b_ord_df),
                    mode = "markers", type = "scatter3d")
  p_3d_b <- p_3d_b %>% layout(scene = list(xaxis = list(title = 'PCoA 1'),
                                           yaxis = list(title = 'PCoA 2'),
                                           zaxis = list(title = 'PCoA 3')))
  htmlwidgets::saveWidget(p_3d_b, "./beta3D_bray.html",
                          selfcontained = TRUE, libdir = NULL)
  
  # Jaccard 
  p_3d_j <- plot_ly(j_ord_df, x = ~Axis.1, y = ~Axis.2, z = ~Axis.3,
                    color = ~j_ord_df[,beta_var], colors = "Dark2",
                    text = row.names(j_ord_df),
                    mode = "markers", type = "scatter3d")
  p_3d_j <- p_3d_j %>% layout(scene = list(xaxis = list(title = 'PCoA 1'),
                                           yaxis = list(title = 'PCoA 2'),
                                           zaxis = list(title = 'PCoA 3')))
  htmlwidgets::saveWidget(p_3d_b, "./beta3D_jaccard.html",
                          selfcontained = TRUE, libdir = NULL)
}

