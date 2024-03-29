shhh <- suppressPackageStartupMessages
shhh(library(reshape2))
shhh(library(Biostrings))
shhh(library(gtools))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(plotly))
shhh(library(matrixStats))
shhh(library(readr))
shhh(library(tidyverse))

############################## EXTRACT SAMPLE FUNC. ############################
extractSampleID <- function(rawSample, stripPattern) {
  sampleID <- basename(rawSample) %>% str_remove(stripPattern)
  return(sampleID)
}

# Reading input
args <- commandArgs(trailingOnly = TRUE)
# check arguments
if (length(args) < 5) {
  stop("Usage: summary_report.R <count_table> <vOTUs_consensus> <metadata> <sum_viol_var> <path>\n")
}
################################## COUNT TABLE #################################
file_count <- args[1]
if (! file.exists(file_count)){
  stop("ERROR: count table file not found in: ", file_count, "\n")
}
count <- read.delim(file_count, sep = "\t", check.names = F)

################################ MULTIFASTA FILE ###############################
file_fasta <- args[2]
if (! file.exists(file_fasta)){
  stop("ERROR: vOTUs consensus fasta file not found in: ", file_fasta, "\n")
}
fasta <- readDNAStringSet(file_fasta)
################################### METADATA ###################################
file_meta <- args[3]
if (! file.exists(file_meta)){
  stop("ERROR: metadata file not found in: ", file_meta, "\n")
}
metadata <- read.delim(file_meta, row.names=1, sep = ",", check.names = F)
metadata$Sample <- rownames(metadata)
metadata <- metadata %>%
  select(Sample, everything())
################################ VIOLIN VARIABLE ###############################
violin_color <- args[4]
if (violin_color == FALSE){
  stop("ERROR: sum_viol_var parameter not found. You can set it in your project config.file (e.g. sum_viol_var = \"metadata_variable1\" \n")
}
################################# FILES PATH ###################################
file_path <- args[5]
if (! file.exists(file_path)){
  stop("ERROR: files path not found in: ", file_path, "\n")
}
# create the violin single plot folder
single_path = file.path(file_path, "report/plots")
dir.create(single_path, showWarnings = FALSE)
single_path = file.path(single_path, "single_violins")
dir.create(single_path, showWarnings = FALSE)

############################# SUMMARY DATAFRAME ################################
# create a dataframe with fasta sequences length
fasta_range = data.frame(ViralOTU = fasta@ranges@NAMES, 
                         length = fasta@ranges@width)
#match count_table vOTU with the fasta sequence
df_count <- merge(fasta_range, count, by= "ViralOTU")
# sort the table numerically
df_count <- df_count[mixedorder(as.character(df_count$ViralOTU)),]
rm(fasta, fasta_range)
# melt df_count
melt_count <- melt(df_count, id.vars = c("ViralOTU", "length"), na.rm = FALSE)
# use ddply to summarize on each sample 
df <- ddply(melt_count, .variables = ~variable, function(x){
  summary(x[x$value>0,"length"])[c("Min.", "Mean", "Max.")]
})
# rename df columns
names(df)[names(df) == "Min."] <- "min_vOTU_length"
names(df)[names(df) == "Max."] <- "max_vOTU_length"
names(df)[names(df) == "Mean"] <- "avg_vOTU_length"
# change rownames of df with sample IDs
rownames(df) <- df$variable
df$variable <- NULL
# create the summary data_frame with vOTUs total count and number of vOTU
# per sample
df_sample = as.data.frame(apply(df_count[,-c(1,2)], 2, sum))
colnames(df_sample) = "vOTUs_tot_abundance"
# number of unique vOTU per sample
dfloop <- df_count[,-c(2)]
dfloop$ViralOTU <- NULL
df_votus <- as.data.frame(colSums(dfloop != 0))
colnames(df_votus) = "vOTUs_number"
rm(df_count, melt_count, dfloop)
# merge the three df
df$Sample <- rownames(df)
df_sample$Sample <- rownames(df_sample)
df_votus$Sample <- rownames(df_votus)
df_final <- join_all(list(df_votus,df_sample,df), by = 'Sample', type = 'full')
head(df_final)
rm(df, df_sample, df_votus)

############################## SINGLE VIOLINS ##################################
# Metadata processing
df_new_meta = data.frame(Sample = metadata$Sample, metadata[[violin_color]])
colnames(df_new_meta)[2] <- paste(violin_color)
# Dataframes creation
#df_tot = data.frame()
# Create a list to hold the plot objects
pltList <- list()
for(sample in colnames(count[1,2:ncol(count)])){
  df_tot = data.frame()
  df_new = data.frame(count$ViralOTU, count[sample], sample)
  colnames(df_new) <- c("viralOTU","Abundance","Sample")
  df_tot <- rbind(df_tot, df_new)
  df_tot$temp <- df_new_meta[,2][match(df_tot$Sample, df_new_meta$Sample)]
  colnames(df_tot)[4] <- paste(violin_color)
  # single violin generation
  # Create plot name
  pltName <- paste(sample, '_violin', sep = '' )
  # Store a plot in the list using the name as an index
  pltList[[pltName]] <- df_tot %>%
    plot_ly(
      x = ~Sample,
      y = ~log10(Abundance),
      type = 'violin',
      box = list(
        visible = T
      ), 
      meanline = list(
        visible = T
      ),
      jitter = 1,
      text = ~viralOTU
    )
  violin = paste0(paste(single_path, pltName, sep = "/"), "_plot.html")
  violin
  htmlwidgets::saveWidget(pltList[[pltName]], violin,
                          selfcontained = TRUE, libdir = NULL)
}

################################ VIOLIN FILES ##################################
# read the coordinates files and create a data.frame of paths
violins_path <- list.files(path = single_path, 
                           pattern = ".+_violin_plot.html", 
                           full.names=T)
df_violin_fp <- as.data.frame(violins_path)
# match the coordinates to the vOTU in taxo table
violins_names = unique( extractSampleID(violins_path, "_violin_plot.html") )
matching_violins = violins_names[violins_names %in% df_final$Sample]
# and to the path of the relative .gff
violins_names_fp = as.data.frame(
  grep(paste(matching_violins,collapse="|"), df_violin_fp$violins_path, value=TRUE)
)
colnames(violins_names_fp) <- c("violins_path")
# modify the matching .gff file path with a relative path
violins_fp = paste0("../", str_extract(violins_names_fp$violins_path,
                                       "report/plots/single_violins/.+_violin_plot.html"))
# create the link table
df_violins = tibble("Sample" = matching_violins, file = violins_names_fp) %>%
  mutate(file = str_replace_all(file,
                                '([^;]*)plots/single_violins/', ''),
         path_violins = file.path(violins_fp),
         Violin_plots = paste0('<a target=_blank href=',
                               path_violins, '>', file,'</a>'))
df_violins$file <- NULL
df_violins$path_violins <- NULL
# Create data.frame with links
df_dt <- join_all(list(df_final,df_violins), by = 'Sample', type = 'full')
df_dt <- df_dt[,c(2,1,3,4,5,6,7)]
# Create DT version of the data.frame
dt_table <- DT::datatable(df_dt, filter="top", rownames = FALSE,
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
) %>% DT::formatRound(columns=c(2:6), digits=2)
dt_table$width = 1500
htmlwidgets::saveWidget(dt_table, "vOTUs_summary.html",
                        selfcontained = TRUE, libdir = NULL)