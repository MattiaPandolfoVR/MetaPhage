# Coded by Mattia Pandolfo (mattia.pandolfo@univr.it)
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
  stop("Usage: beta_diversity.R <count_table> <taxonomy_table> <metadata> <beta_var>\n")
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
beta_var <- args[4]
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
colnames(b_ord_df) <- gsub(pattern = "Axis.", replacement = "Bray_axis.", x=colnames(b_ord_df))
b_ord_df$Sample <- rownames(b_ord_df)
# Jaccard
set.seed(12345)
j_ord = ordinate(ps, method = "PCoA", "jaccard")
j_ord_df <- as.data.frame(j_ord$vectors)
colnames(j_ord_df) <- gsub(pattern = "Axis.", replacement = "Jac_axis.", x=colnames(j_ord_df))
j_ord_df$Sample <- rownames(j_ord_df)

#check if data are suitable for 2d plots
if (ncol(b_ord_df)<2 || ncol(j_ord_df) <2){
 message ("Sample numerosity in your data does not hallow 2D plots creation!")
} else {
  # how many values in beta var?
  # Bray-Curtis
  # create a column for the user metadata of interest and populate it
  b_ord_df[,beta_var] <- sample_data(ps) %>%
    data.frame() %>%
    select(all_of(beta_var)) %>%
    mutate_if(is.factor,as.character)
  b_ord_df <- b_ord_df[order(b_ord_df[,beta_var]), ]
  # Jaccard
  # create a column for the user metadata of interest and populate it
  j_ord_df[,beta_var] <- sample_data(ps) %>%
    data.frame() %>%
    select(all_of(beta_var)) %>%
    mutate_if(is.factor,as.character)
  j_ord_df <- j_ord_df[order(j_ord_df[,beta_var]), ]
  # Merge the single betas in one
  df_beta <- merge(b_ord_df, j_ord_df, by = c("Sample", beta_var))
  df_beta[,beta_var] <- as.factor(df_beta[,beta_var])
  n_beta <- length(levels(df_beta[,beta_var]))
  # Generate how many true and false to plot the correct graph with buttons
  if (n_beta == 1) {
    items_TOT = TRUE
    items_TOT_inv = FALSE
  } else {
    # T and F for Bray-curtis
    items_TRUE <- paste(rep(TRUE, n_beta))
    items_FALSE <- paste(rep(FALSE, n_beta))
    items_TOT <- as.logical(c(items_TRUE,items_FALSE))
    # F and T for Jaccard
    items_TOT_inv <- as.logical(c(items_FALSE,items_TRUE))
  }
  # create item lists (for 2D plots)
  items2d <- list(
    list(label="Bray_Curtis", method = "update", args=list(list(visible=items_TOT),
                                                           list(title = "2D-Beta diversity<br>Bray_Curtis"))),
    list(label="Jaccard", method = "update", args=list(list(visible=items_TOT_inv),
                                                       list(title = "2D-Beta diversity<br>Jaccard")))
  )
  # Plot 2D-FIGURE
  fig_beta_2d <- plot_ly(data = df_beta, x = ~Bray_axis.1, y = ~Bray_axis.2, color = ~df_beta[,beta_var],
                         colors = "Dark2",  visible=T, type = "scatter", mode = "markers", groupclick = "toggleitem",
                         hovertext = paste("Sample: ", df_beta$Sample),
                         alpha = 0.8, size = 800) %>%
    add_trace(x = ~Jac_axis.1, y = ~Jac_axis.2, color = ~df_beta[,beta_var],
              colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_beta$Sample)) %>%
    layout(
      title = "2D-Beta diversity<br>Select the metric from the dropdown menu",
      xaxis = list(title = "PCoA 1"),
      yaxis = list(title = "PCoA 2"),
      legend= list(title = list(text= beta_var),
                   traceorder = "normal",
                   itemsizing = "constant"),
      autosize = TRUE,
      boxmode = "group",
      margin = list(autosize = TRUE, width = 1400, height = 900),
      updatemenus = list(
        list(buttons = items2d)
      )
    )
  htmlwidgets::saveWidget(fig_beta_2d, "./beta2D.html",
                         selfcontained = TRUE, libdir = NULL)
}
################################ 3D PLOTS #####################################
# check if data are suitable for 3d plots
if (ncol(b_ord_df)<3 || ncol(j_ord_df)<3){
  message ("Sample numerosity in your data doesn't hallow 3D plots creation!")
} else {
  # how many values in beta var?
  # Bray-Curtis
  # create a column for the user metadata of interest and populate it
  b_ord_df[,beta_var] <- sample_data(ps) %>%
    data.frame() %>%
    select(all_of(beta_var)) %>%
    mutate_if(is.factor,as.character)
  b_ord_df <- b_ord_df[order(b_ord_df[,beta_var]), ]
  # Jaccard
  # create a column for the user metadata of interest and populate it
  j_ord_df[,beta_var] <- sample_data(ps) %>%
    data.frame() %>%
    select(all_of(beta_var)) %>%
    mutate_if(is.factor,as.character)
  j_ord_df <- j_ord_df[order(j_ord_df[,beta_var]), ]
  # Merge the single betas in one
  df_beta <- merge(b_ord_df, j_ord_df, by = c("Sample", beta_var))
  df_beta[,beta_var] <- as.factor(df_beta[,beta_var])
  n_beta <- length(levels(df_beta[,beta_var]))
  # Generate how many true and false to plot the correct graph with buttons
  if (n_beta == 1) {
    items_TOT = TRUE
    items_TOT_inv = FALSE
  } else {
    # T and F for Bray-curtis
    items_TRUE <- paste(rep(TRUE, n_beta))
    items_FALSE <- paste(rep(FALSE, n_beta))
    items_TOT <- as.logical(c(items_TRUE,items_FALSE))
    # F and T for Jaccard
    items_TOT_inv <- as.logical(c(items_FALSE,items_TRUE))
  }
  # create item lists (for 3D plots)
  items3d <- list(
    list(label="Bray_Curtis", method = "update", args=list(list(visible=items_TOT),
                                                           list(title = "3D-Beta diversity<br>Bray_Curtis"))),
    list(label="Jaccard", method = "update", args=list(list(visible=items_TOT_inv),
                                                       list(title = "3D-Beta diversity<br>Jaccard")))
  )
  # Plot 3D-FIGURE
  fig_beta_3d <- plot_ly(data = df_beta, x = ~Bray_axis.1, y = ~Bray_axis.2, z = ~Bray_axis.3,
                         color = ~df_beta[,beta_var], colors = "Dark2",
                         text = paste("Sample: ", df_beta$Sample),
                         type = "scatter3d", mode = "markers",
                         visible=T, hovertext = paste("Sample: ", df_beta$Sample),
                         alpha = 0.8) %>%
    add_trace(x = ~Jac_axis.1, y = ~Jac_axis.2, z = ~Jac_axis.3, color = ~df_beta[,beta_var],
              colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_beta$Sample)) %>%
    layout(
      title = "3D-Beta diversity<br>Select the metric from the dropdown menu",
      xaxis = list(title = "PCoA 1"),
      yaxis = list(title = "PCoA 2"),
      yaxis = list(title = "PCoA 3"),
      legend= list(title = list(text= beta_var),
                   traceorder = "normal",
                   itemsizing = "constant"),
      autosize = TRUE,
      boxmode = "group",
      margin = list(autosize = TRUE, width = 1400, height = 900),
      updatemenus = list(
        list(buttons = items3d)
      )
    )
  htmlwidgets::saveWidget(fig_beta_3d, "./beta3D.html",
                          selfcontained = TRUE, libdir = NULL)
}
