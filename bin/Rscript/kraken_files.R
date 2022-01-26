shhh <- suppressPackageStartupMessages
shhh(library(DT))
shhh(library(readr))
shhh(library(dplyr))
shhh(library(plyr))
shhh(library(seqinr))
shhh(library(tidyverse))
shhh(library(gtools))

# Reading input
args <- commandArgs(trailingOnly = TRUE)
############################# FILES AND PATHS ##################################
file_paths <- args[1]
kraken_path <- file.path(file_paths, "taxonomy/kraken2")
krona_path <-  file.path(file_paths, "taxonomy/krona")

############################### METADATA #######################################
metadata <- args[2]                                                             
meta <- read.csv(metadata, sep = ',')

################################ KRAKEN FILES ##################################
# read the kraken files and create a data.frame of paths
krakens_path <- list.files(path = kraken_path, 
                           pattern = "SRR\\d+_report.txt", 
                           recursive = TRUE,
                           full.names=T)
df_kraken_fp <- as.data.frame(krakens_path)
# match the kraken files to the metadata sample names
krakens_names = unique(str_extract(krakens_path, "SRR\\d+"))
matching_krakens = krakens_names[krakens_names %in% meta$Sample]
# and create a dataframe of theese
krakens_names_fp =  df_kraken_fp[match(matching_krakens,
                                       str_extract(krakens_path, "SRR\\d+")),]
# modify the matching kraken file path with a relative path
krakens_fp = paste0("../", str_extract(krakens_names_fp,
                                       "taxonomy/kraken2/SRR\\d+/SRR\\d+_report.txt"))
# create the link table
df_krakens = tibble("Sample" = matching_krakens, file = krakens_fp) %>%
  mutate(file = str_replace_all(file,
                                '([^;]*)taxonomy/kraken2/SRR\\d+/', ''),
         path_krakens = file.path(krakens_fp),
         Kraken_Report = paste0('<a target=_blank href=',
                                path_krakens, '>', file,'</a>'))
df_krakens$file <- NULL
df_krakens$path_krakens <- NULL

################################ KRONA FILES ##################################
# check if the krona folder exist
if(file.exists(krona_path)){
  # read the krona files and create a data.frame of paths
  kronas_path <- list.files(path = krona_path, 
                            "SRR\\d+_krak_krona_abundancies.html", 
                            recursive = TRUE,
                            full.names=T)
  df_krona_fp <- as.data.frame(kronas_path)
  # match the krona files to the metadata sample names
  kronas_names = unique(str_extract(kronas_path, "SRR\\d+"))
  matching_kronas = kronas_names[kronas_names %in% meta$Sample]
  # and create a dataframe of theese
  kronas_names_fp =  df_krona_fp[match(matching_kronas,
                                      str_extract(kronas_path, "SRR\\d+")),]
  # modify the matching krona file path with a relative path
  kronas_fp = paste0("../", str_extract(kronas_names_fp,
                                        "taxonomy/krona/SRR\\d+/SRR\\d+_krak_krona_abundancies.html"))
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
  htmlwidgets::saveWidget(dt_KK, "kraken_files.html",
                          selfcontained = TRUE, libdir = NULL)
}
