shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(plyr))
shhh(library(plotly))
shhh(library(phyloseq))
shhh(library(data.table))

# Reading input
args <- commandArgs(trailingOnly = TRUE)
############################### COUNT TABLE ####################################
file_count <- args[1]
count <- read.csv(file_count, sep = '\t', row.names = 1)
############################## TAXONOMY TABLE ##################################
file_taxo <- args[2]
taxo <- read.csv(file_taxo, sep = '\t', row.names = 1)
# Extract the scaffold id and the last 2 column from taxonomy table
taxo <- taxo[,c("Host","Family","Genus")]
################################ METADATA ######################################
metadata <- args[3]                                                             
metadata <- read.csv(metadata, sep = ',', row.names = 1)
metadata$Sample <- rownames(metadata)
############################### META_VARS ######################################
alpha_var1 <- args[4]
alpha_var2 <- args[5]
# Check if only one variable is passed
if(alpha_var2 == FALSE){
  alpha_var2 <- alpha_var1
}
# Phyloseq object creation
ps <- phyloseq(otu_table(count, taxa_are_rows = TRUE),
               tax_table(as.matrix(taxo)),
               sample_data(metadata))
# Remove NAs from taxonomy table
ps <- subset_taxa(ps, !is.na(Family) & !is.na(Genus))
# Restore sample_data rownames
row.names(ps@sam_data) <- c(1:nrow(ps@sam_data))
# Calculate alpha_diversity measurments
df_div <- estimate_richness(ps)
df_div$Sample <- rownames(df_div)
# Transform ps@sample_data in a df
sample_df <- as.matrix(ps@sam_data)
sample_df <- as.data.frame(sample_df)
# Merge the two data-frames
df_alpha <- merge(sample_df, df_div, by="Sample")
# And melt it
df_melt <- melt(df_alpha, na.rm = FALSE)

################################# PLOTS ########################################
# Observed subplot
# create annotations for y axis
y_Observed <- list(title = "Observed")
Observed <- df_melt %>% filter(variable == "Observed") %>% 
  plot_ly(
    x = ~get(alpha_var1), 
    y = ~value, 
    color = ~get(alpha_var2),
    colors = "Dark2",
    type = "box",
    legendgroup = ~get(alpha_var2),
    showlegend = TRUE) %>% 
  subplot(nrows = 1, shareX = TRUE, shareY =TRUE, titleX = FALSE) %>% 
  layout(title = "Alpha-diversity by different indexes",
         boxmode = "group",
         groupclick = "toggleitem",
         yaxis = y_Observed)
# Chao1 subplot
y_Chao1 <- list(title = "Chao1")
Chao1 <- df_melt %>% filter(variable == "Chao1") %>%
  plot_ly(
    x = ~get(alpha_var1), 
    y = ~value, 
    color = ~get(alpha_var2),
    colors = "Dark2",
    type = "box",
    legendgroup = ~get(alpha_var2),
    showlegend = FALSE) %>% 
  subplot(nrows = 1, shareX = TRUE, shareY =TRUE, titleX = FALSE) %>% 
  layout(boxmode = "group", yaxis = y_Chao1)
# se.chao1 subplot
y_se.Chao1 <- list(title = "se.Chao1")
se.Chao1 <- df_melt %>% filter(variable == "se.chao1") %>%
  plot_ly(
    x = ~get(alpha_var1), 
    y = ~value, 
    color = ~get(alpha_var2),
    colors = "Dark2",
    type = "box",
    legendgroup = ~get(alpha_var2),
    showlegend = FALSE) %>% 
  subplot(nrows = 1, shareX = TRUE, shareY =TRUE, titleX = FALSE) %>% 
  layout(boxmode = "group", yaxis = y_se.Chao1)
# ACE subplot
y_ACE <- list(title = "ACE")
ACE <- df_melt %>% filter(variable == "ACE") %>%
  plot_ly(
    x = ~get(alpha_var1), 
    y = ~value, 
    color = ~get(alpha_var2),
    colors = "Dark2",
    type = "box",
    legendgroup = ~get(alpha_var2),
    showlegend = TRUE) %>% 
  subplot(nrows = 1, shareX = TRUE, shareY =TRUE, titleX = FALSE) %>% 
  layout(boxmode = "group", groupclick = "toggleitem", yaxis = y_ACE)
# se.ACE subplot
y_se.ACE <- list(title = "se.ACE")
se.ACE <- df_melt %>% filter(variable == "se.ACE") %>%
  plot_ly(
    x = ~get(alpha_var1), 
    y = ~value, 
    color = ~get(alpha_var2),
    colors = "Dark2",
    type = "box",
    legendgroup = ~get(alpha_var2),
    showlegend = FALSE) %>% 
  subplot(nrows = 1, shareX = TRUE, shareY =TRUE, titleX = FALSE) %>% 
  layout(boxmode = "group", yaxis = y_se.ACE)
# Shannon subplot
y_Shannon <- list(title = "Shannon")
Shannon <- df_melt %>% filter(variable == "Shannon") %>%
  plot_ly(
    x = ~get(alpha_var1), 
    y = ~value, 
    color = ~get(alpha_var2),
    colors = "Dark2",
    type = "box",
    legendgroup = ~get(alpha_var2),
    showlegend = FALSE) %>% 
  subplot(nrows = 1, shareX = TRUE, shareY =TRUE, titleX = FALSE) %>% 
  layout(boxmode = "group", yaxis = y_Shannon)
# Simpson subplot
y_Simpson <- list(title = "Simpson")
Simpson <- df_melt %>% filter(variable == "Simpson") %>%
  plot_ly(
    x = ~get(alpha_var1), 
    y = ~value, 
    color = ~get(alpha_var2),
    colors = "Dark2",
    type = "box",
    legendgroup = ~get(alpha_var2),
    showlegend = TRUE) %>% 
  subplot(nrows = 1, shareX = TRUE, shareY =TRUE, titleX = FALSE) %>% 
  layout(boxmode = "group", groupclick = "toggleitem", yaxis = y_Simpson)
# InvSimpson subplot
y_InvSimpson <- list(title = "InvSimpson")
InvSimpson <- df_melt %>% filter(variable == "InvSimpson") %>%
  plot_ly(
    x = ~get(alpha_var1), 
    y = ~value, 
    color = ~get(alpha_var2),
    colors = "Dark2",
    type = "box",
    legendgroup = ~get(alpha_var2),
    showlegend = FALSE) %>% 
  subplot(nrows = 1, shareX = TRUE, shareY =TRUE, titleX = FALSE) %>% 
  layout(boxmode = "group", yaxis = y_InvSimpson)
# Fisher subplot
y_Fisher <- list(title = "Fisher")
Fisher <- df_melt %>% filter(variable == "Fisher") %>%
  plot_ly(
    x = ~get(alpha_var1), 
    y = ~value, 
    color = ~get(alpha_var2),
    colors = "Dark2",
    type = "box",
    legendgroup = ~get(alpha_var2),
    showlegend = FALSE) %>% 
  subplot(nrows = 1, shareX = TRUE, shareY =TRUE, titleX = FALSE) %>% 
  layout(boxmode = "group", yaxis = y_Fisher)
# Combine subplots
p_OCS <- subplot(nrows = 1,
                 Observed,Chao1,se.Chao1,
                 titleX = TRUE, titleY = TRUE)
p_ASS <- subplot(nrows = 1,
                 ACE,se.ACE,Shannon,
                 titleX = TRUE, titleY = TRUE)
p_SIF <- subplot(nrows = 1,
                 Simpson,InvSimpson,Fisher,
                 titleX = TRUE, titleY = TRUE)
htmlwidgets::saveWidget(p_OCS, "alpha_div_Observed_Chao1_se.Chao1.html",
                        selfcontained = TRUE, libdir = NULL)
htmlwidgets::saveWidget(p_ASS, "alpha_div_ACE_se.ACE_Shannon.html",
                        selfcontained = TRUE, libdir = NULL)
htmlwidgets::saveWidget(p_SIF, "alpha_div_Simpson_InvSimpson_Fisher.html",
                        selfcontained = TRUE, libdir = NULL)