# Coded by Mattia Pandolfo (mattia.pandolfo@univr.it)
shhh <- suppressPackageStartupMessages
shhh(library(DT))
shhh(library(readr))
shhh(library(dplyr))
shhh(library(plyr))
shhh(library(seqinr))
shhh(library(tidyverse))
shhh(library(gtools))
shhh(library(magrittr))

# Reading input
args <- commandArgs(trailingOnly = TRUE)
# check arguments
if (length(args) != 2) {
  stop("Usage: kraken_files.R <outdir_path> <metadata>")
}
############################# FILES AND PATHS ##################################
file_paths <- args[1]
# Check if kraken_path and krona_path exist
kraken_path <- file.path(file_paths, "report/taxonomy/kraken2")
if (! file.exists(kraken_path)){
  stop("ERROR: report/taxonomy/kraken2 folder not found in: ", kraken_path, "\n")
}
krona_path <-  file.path(file_paths, "report/taxonomy/krona")
if (! file.exists(krona_path)){
  warning("WARNING: report/taxonomy/krona folder not found in: ", krona_path, "\n")
}
############################### METADATA #######################################
file_meta <- args[2]
if (! file.exists(file_meta)){
  stop("ERROR: metadata file not found in: ", file_meta, "\n")
}
metadata <- read.delim(file_meta, row.names=1, sep = ",", check.names = F)
metadata$Sample <- rownames(metadata)
metadata <- metadata %>%
  select(Sample, everything())
################################ KRAKEN FILES ##################################
# read the kraken files and create a data.frame of paths
krakens_path <- list.files(path = kraken_path, 
                           pattern = "_report.txt", 
                           recursive = TRUE,
                           full.names=T)
# create a df of kraken files path
df_kraken_fp <- as.data.frame(krakens_path)
# match the kraken files to the metadata sample names
krakens_names = unique(basename(krakens_path)) 
krakens_names <- sapply(strsplit(krakens_names,"_report.txt"), `[`, 1)
matching_krakens = krakens_names[krakens_names %in% metadata$Sample]
# and create a dataframe of theese
krakens_names_fp = as.data.frame(
  grep(paste(matching_krakens,collapse="|"), df_kraken_fp$krakens_path, value=TRUE)
)
colnames(krakens_names_fp) <- c("krakens_path")
# modify the matching kraken file path with a relative path: like
krakens_fp = paste0("./", str_extract(krakens_names_fp$krakens_path,
                                       "taxonomy/kraken2/.+/.+_report.txt"))
# create the link table
df_krakens = tibble("Sample" = matching_krakens, file = krakens_fp) %>%
  mutate(file = str_replace_all(file,
                                '([^;]*)taxonomy/kraken2/.+/', ''),
         path_krakens = file.path(krakens_fp),
         Kraken_Report = paste0('<a target=_blank href=',
                                path_krakens, '>', file,'</a>'))
df_krakens$file <- NULL
df_krakens$path_krakens <- NULL

################################ KRONA FILES ###################################
# check if the krona folder exist
if(file.exists(krona_path)){
  # read the krona files and create a data.frame of paths
  kronas_path <- list.files(path = krona_path, 
                            ".+_krak_krona_abundancies.html", 
                            recursive = TRUE,
                            full.names=T)
  
  df_krona_fp <- as.data.frame(kronas_path)
  # match the krona files to the metadata sample names
  kronas_names = unique(basename(kronas_path)) 
  kronas_names <- sapply(strsplit(kronas_names,"_krak_krona_abundancies.html"), `[`, 1)
  matching_kronas = kronas_names[kronas_names %in% metadata$Sample]
  # and create a dataframe of theese
  kronas_names_fp <- filter(df_krona_fp,
                            matching_kronas %in% matching_kronas)
  # modify the matching krona file path with a relative path
  kronas_fp = paste0("./", str_extract(kronas_names_fp$kronas_path,
                                        "taxonomy/krona/.+/.+_krak_krona_abundancies.html"))
  # create the link table
  df_kronas = tibble("Sample" = matching_kronas, file = kronas_fp) %>%
    mutate(file = str_replace_all(file,
                                  '([^;]*)taxonomy/krona/.*\\d+/', ''),
           path_kronas = file.path(kronas_fp),
           Krona_Plot = paste0('<a target=_blank href=',
                               path_kronas, '>', file,'</a>'))
  df_kronas$file <- NULL
  df_kronas$path_kronas <- NULL
  # Merge the two df
  df_kk <- join_all(list(df_krakens,df_kronas), by = 'Sample', type = 'full')
  # Create DT version of the data.frame
  dt_KK <- DT::datatable(df_kk, filter="top", rownames = FALSE,
                         extensions = c('Buttons','Scroller'), escape = F, 
                         options = list(scrollY = "650px",
                                        scrollX = TRUE,
                                        scroller = TRUE,
                                        dom = 'BlPfrtip',
                                        buttons = list('copy', 'print')
                         )
  )
  dt_KK$width = 1500
  htmlwidgets::saveWidget(dt_KK, "./kraken_files.html",
                          selfcontained = TRUE, libdir = NULL)
} else {
  # Create DT version of the data.frame with kraken files paths only
  dt_K <- DT::datatable(df_krakens, filter="top", rownames = FALSE,
                        extensions = c('Buttons','Scroller'), escape = F, 
                        options = list(scrollY = "650px",
                                       scrollX = TRUE,
                                       scroller = TRUE,
                                       dom = 'BlPfrtip',
                                       buttons = list('copy', 'print')
                        )
  )
  dt_K$width = 1200
  htmlwidgets::saveWidget(dt_K, "./kraken_files.html",
                          selfcontained = TRUE, libdir = NULL)
}
