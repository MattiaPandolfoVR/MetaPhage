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
################################### BETA_VAR ###################################
beta_var <- args[4]

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

################################# 2D PLOTS #####################################
# Beta-Diversity measurments
# Bray-Curtis
set.seed(12345)
b_ord = ordinate(ps, method = "PCoA", "bray")
b_comps <- b_ord$vectors
# Jaccard
set.seed(12345)
j_ord = ordinate(ps, method = "PCoA", "jaccard", )
j_comps <- j_ord$vectors
# Plots axis labels
x <- list(
  title = "PCoA 1"
)
y <- list(
  title = "PCoA 2"
)
# Bray-Curtis
b_ord_df = as.data.frame(b_comps)
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
j_ord_df = as.data.frame(j_comps)
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

################################# 3D PLOTS #####################################
# check if data are suitable for 3d plots
# Bray-Curtis
if (ncol(b_comps)<3){
  message ("Sample numerosity in your data doesn't hallow 3D plots creation!")
} else {
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
