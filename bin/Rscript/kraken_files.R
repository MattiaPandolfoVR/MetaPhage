#!/usr/bin/env Rscript
# A script to generate the HTML table with the links
# to Kraken/Krona files

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
if (length(args) != 2) {
  stop("Usage: kraken_files.R <input_path> <metadata>")
}

############################# FILES AND PATHS ##################################
file_paths <- args[1]

cat(" Input: ", file_paths, "\n")
kraken_path <- file.path(file_paths, "taxonomy/kraken2")
krona_path <-  file.path(file_paths, "taxonomy/krona")
if (! file.exists(kraken_path)){
  stop("ERROR: taxonomy/kraken2 not found in: ", kraken_path, "\n")
}
if (! file.exists(krona_path)){
  cat("WARNING: taxonomy/krona not found in: ", krona_path, "\n")
}
# Check if kraken_path and krona_path exist

############################### METADATA #######################################
metadata <- args[2]    
cat(" Metadata: ", metadata, "\n") 
if (! file.exists(metadata)){
  stop("ERROR: metadata not found in: ", metadata, "\n")
}                                                        
meta <- read.csv(metadata, sep = ',')
head(meta$Sample)
#deb#head(meta)
################################ KRAKEN FILES ##################################
# read the kraken files and create a data.frame of paths
krakens_path <- list.files(path = kraken_path, 
                           pattern = "_report.txt", 
                           recursive = TRUE,
                           full.names=T)

# [1] "/Users/telatina/Downloads/MetaPhage//taxonomy/kraken2/SRR8652861/SRR8652861_report.txt"
# [2] "/Users/telatina/Downloads/MetaPhage//taxonomy/kraken2/SRR8652914/SRR8652914_report.txt"
# [3] "/Users/telatina/Downloads/MetaPhage//taxonomy/kraken2/SRR8652969/SRR8652969_report.txt"
cat("Raw Files: ", length(krakens_path),"\n")
df_kraken_fp <- as.data.frame(krakens_path)
cat("Files: ", length(df_kraken_fp$krakens_path),"\n")

#df_kraken_fp
#                                                                             krakens_path
#1  /Users/telatina/Downloads/MetaPhage//taxonomy/kraken2/SRR8652861/SRR8652861_report.txt
#2  /Users/telatina/Downloads/MetaPhage//taxonomy/kraken2/SRR8652914/SRR8652914_report.txt
#3  /Users/telatina/Downloads/MetaPhage//taxonomy/kraken2/SRR8652969/SRR8652969_report.txt
#4  /Users/telatina/Downloads/MetaPhage//taxonomy/kraken2/SRR8653084/SRR8653084_report.txt


# match the kraken files to the metadata sample names
krakens_names = unique(basename(krakens_path)) 

# [1] "SRR8652861" "SRR8652914" "SRR8652969" "SRR8653084" "SRR8653090" "SRR8653218" "SRR8653221" "SRR8653245"
# [9] "SRR8653247" "SRR8653248"

krakens_names <- sapply(strsplit(krakens_names,"_report.txt"), `[`, 1)

cat("Kraken names: ", length(krakens_names), " ", head(krakens_names), "\n")
matching_krakens = krakens_names[krakens_names %in% meta$Sample]
cat("Kraken matched: ", length(matching_krakens), " ", head(matching_krakens),"\n")

#matching_krakens
#[1] "SRR8652861" "SRR8652914" "SRR8652969" "SRR8653084" "SRR8653090" ...


# and create a dataframe of theese

#krakens_names_fp =  df_kraken_fp[match(matching_krakens,
#                                       str_extract(krakens_path, "SRR\\d+")),]




krakens_names_fp = as.data.frame(
  grep(paste(matching_krakens,collapse="|"), df_kraken_fp$krakens_path, value=TRUE)
)
colnames(krakens_names_fp) <- c("krakens_path")
cat("Kraken absolute: ", length(krakens_names_fp$krakens_path), "\n")
#> krakens_names_fp
#[1] "/Users/telatina/Downloads/MetaPhage//taxonomy/kraken2/SRR8652861/SRR8652861_report.txt"
#[2] "/Users/telatina/Downloads/MetaPhage//taxonomy/kraken2/SRR8652914/SRR8652914_report.txt"


# modify the matching kraken file path with a relative path: like
# [1] "../taxonomy/kraken2/SRR8652861/SRR8652861_report.txt" ...

krakens_fp = paste0("../", str_extract(krakens_names_fp$krakens_path,
                                       "taxonomy/kraken2/.+/.+_report.txt"))

cat("Kraken relative: ", length(krakens_fp), "\n")
# create the link table
df_krakens = tibble("Sample" = matching_krakens, file = krakens_fp) %>%
  mutate(file = str_replace_all(file,
                                '([^;]*)taxonomy/kraken2/.+/', ''),
         path_krakens = file.path(krakens_fp),
         Kraken_Report = paste0('<a target=_blank href=',
                                path_krakens, '>', file,'</a>'))
df_krakens$file <- NULL
df_krakens$path_krakens <- NULL

################################ KRONA FILES ##################################
# check if the krona folder exist
if(file.exists(krona_path)){
  cat("Create DT version of the data.frame with kraken AND krona files\n")
  # read the krona files and create a data.frame of paths
  kronas_path <- list.files(path = krona_path, 
                            ".+_krak_krona_abundancies.html", 
                            recursive = TRUE,
                            full.names=T)
  
  df_krona_fp <- as.data.frame(kronas_path)
  kronas_path
  # match the krona files to the metadata sample names
  kronas_names = unique(basename(kronas_path)) 
  kronas_names
  kronas_names <- sapply(strsplit(kronas_names,"_krak_krona_abundancies.html"), `[`, 1)
  
  
  matching_kronas = kronas_names[kronas_names %in% meta$Sample]
  
  # and create a dataframe of theese
  kronas_names_fp <- filter(df_krona_fp,
                            matching_kronas %in% matching_kronas)
  
  # modify the matching krona file path with a relative path
  
  kronas_fp = paste0("../", str_extract(kronas_names_fp$kronas_path,
                                        "taxonomy/krona/.+/.+_krak_krona_abundancies.html"))
  # create the link table
  df_kronas = tibble("Sample" = matching_kronas, file = kronas_fp) %>%
    mutate(file = str_replace_all(file,
                                  '([^;]*)taxonomy/krona/SRR\\d+/', ''),
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
  htmlwidgets::saveWidget(dt_KK, "kraken_files.html",
                          selfcontained = TRUE, libdir = NULL)
} else {
  # Create DT version of the data.frame with kraken files paths only
  cat("Create DT version of the data.frame with kraken files paths only\n")
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
  htmlwidgets::saveWidget(dt_K, "kraken_files.html",
                          selfcontained = TRUE, libdir = NULL)
}
