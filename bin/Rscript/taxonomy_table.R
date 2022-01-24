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
############################## TAXONOMY TABLE ##################################
file_df <- args[1]
df <- read.csv(file_df, sep = '\t')
# mutate the taxonomy data.frame
df <- mutate_if(df, is.character, as.factor)
################################## FILES #######################################
file_fna <- args[2]
fna <- read.fasta(file_fna, set.attributes = FALSE, as.string=TRUE,
                  forceDNAtolower = FALSE)
file_faa <- args[3]
faa <- read.fasta(file_faa, set.attributes = FALSE, as.string=TRUE,
                  forceDNAtolower = FALSE)
file_gff <- args[4]
gff <- read.delim(file_gff, header=F, comment.char="#")
################################## PATHS #######################################
file_paths <- args[5]
fna_path <- file.path(file_paths, "cd-hit/vOTUs_consensus")
faa_path <- file.path(file_paths, "cd-hit/vOTUs_proteins")
gff_path <- file.path(file_paths, "cd-hit/vOTUs_coords")
graph_path <- file.path(file_paths, "taxonomy/vcontact2/single-views_megahit")
kraken_path <- file.path(file_paths, "taxonomy/kraken2")
krona_path <-  file.path(file_paths, "taxonomy/krona")

############################ MINER COMPARISON ##################################
file_min <- args[6]
min_com <- read.csv(file_min, sep = '')
############################### METADATA #######################################
metadata <- args[7]                                                             
meta <- read.csv(metadata, sep = ',')

############################## CONSENSUS FILES #################################
# read the consensus fasta files and create a data.frame of paths
consensus_path <- list.files(path = fna_path, 
                               pattern = "vOTU*", 
                               full.names=T)
df_fasta_fp <- as.data.frame(consensus_path)
# match the consensus to the vOTU in taxo table
matching_cons = names(fna)[names(fna) %in% df$ViralOTU]
# and to the path of the relative .fasta
match_cons_fp =  df_fasta_fp[match(matching_cons, 
                           str_extract(consensus_path, 
                                       "vOTU_\\d+")),]
# modify the matching fasta file path with a relative path
consensus_fp = paste0("../", str_extract(match_cons_fp,
                              "cd-hit/vOTUs_consensus/vOTU_\\d+.fasta"))
# create the link table
df_consensus = tibble("ViralOTU" = matching_cons, file= consensus_fp) %>%
  mutate(file = str_replace_all(file,
                                '([^;]*)cd-hit/vOTUs_consensus/', ''),
         path_fasta = file.path(consensus_fp),
         Consensus = paste0('<a target=_blank href=',
                            path_fasta, '>', file,'</a>'))
df_consensus$file <- NULL
df_consensus$path_fasta <- NULL

################################ PROTEIN FILES #################################
# read the protein files and create a data.frame of paths
proteins_path <- list.files(path = faa_path, 
                              pattern = "vOTU*", 
                              full.names=T)
df_prot_fp <- as.data.frame(proteins_path)
# match the proteins to the vOTU in taxo table
prot_names = unique(str_extract(names(faa), "vOTU_[0-9]+"))
matching_prots = prot_names[prot_names %in% df$ViralOTU]
# and to the path of the relative .faa
matching_prots_fp =  df_prot_fp[match(matching_prots,
                                     str_extract(proteins_path, "vOTU_\\d+")),]
# modify the matching .faa file path with a relative path
protein_fp = paste0("../", str_extract(matching_prots_fp,
                                "cd-hit/vOTUs_proteins/vOTU_\\d+.faa"))
# create the link table
df_proteins = tibble("ViralOTU" = matching_prots, file = protein_fp) %>%
  mutate(file = str_replace_all(file, 
                                '([^;]*)cd-hit/vOTUs_proteins/', ''),
         path_prot = file.path(protein_fp),
         Proteins = paste0('<a target=_blank href=',
                           path_prot, '>', file,'</a>'))
df_proteins$file <- NULL
df_proteins$path_prot <- NULL

################################# GFF FILES ####################################
# read the coordinates files and create a data.frame of paths
coords_path <- list.files(path = gff_path, 
                              pattern = "vOTU*", 
                              full.names=T)
df_coord_fp <- as.data.frame(coords_path)
# match the coordinates to the vOTU in taxo table
coord_names = unique(gff$V1)
matching_coords = coord_names[coord_names %in% df$ViralOTU]
# and to the path of the relative .gff
matching_coords_fp =  df_coord_fp[match(matching_coords,
                                     str_extract(coords_path, "vOTU_\\d+")),]
# modify the matching .gff file path with a relative path
coords_fp = paste0("../", str_extract(matching_coords_fp,
                                "cd-hit/vOTUs_coords/vOTU_\\d+.gff"))
# create the link table
df_coords = tibble("ViralOTU" = matching_coords, file = coords_fp) %>%
  mutate(file = str_replace_all(file,
                                '([^;]*)cd-hit/vOTUs_coords/', ''),
         path_coords = file.path(coords_fp),
         Coordinates = paste0('<a target=_blank href=',
                              path_coords, '>', file,'</a>'))
df_coords$file <- NULL
df_coords$path_coords <- NULL

################################ GRAPHS FILES ##################################
# read the coordinates files and create a data.frame of paths
graphs_path <- list.files(path = graph_path, 
                          pattern = "vOTU*", 
                          full.names=T)
df_graphs_fp <- as.data.frame(graphs_path)
# match the coordinates to the vOTU in taxo table
graphs_names = unique(str_extract(graphs_path, "vOTU_[0-9]+"))
matching_graphs = graphs_names[graphs_names %in% df$ViralOTU]
# and to the path of the relative .gff
graphs_names_fp =  df_graphs_fp[match(matching_graphs,
                                        str_extract(graphs_path, "vOTU_\\d+")),]
# modify the matching .gff file path with a relative path
graphs_fp = paste0("../", str_extract(graphs_names_fp,
                                     "taxonomy/vcontact2/single-views_megahit/vOTU_\\d+.html"))
# create the link table
df_graphs = tibble("ViralOTU" = matching_graphs, file.z = graphs_fp) %>%
  mutate(file.z = str_replace_all(file.z,
                                '([^;]*)taxonomy/vcontact2/single-views_megahit/', ''),
         path_graphs = file.path(graphs_fp),
         Neighbors = paste0('<a target=_blank href=',
                              path_graphs, '>', file.z,'</a>'))
df_graphs$file.z <- NULL
df_graphs$path_graphs <- NULL

############################## PER vOTU MINERS #################################
df_mincomp <- merge(df[,c("ViralOTU","Status")],
                    min_com[, c("ViralOTU",
                                "Phigaro",
                                "Vibrant",
                                "Virfinder",
                                "Virsorter")],
                    by= "ViralOTU")
df_mincomp$Status <- NULL
df_mincomp$Phigaro <- gsub(1,"Yes", df_mincomp$Phigaro)
df_mincomp$Phigaro <- gsub(0,"No", df_mincomp$Phigaro)
df_mincomp$Vibrant <- gsub(1,"Yes", df_mincomp$Vibrant)
df_mincomp$Vibrant <- gsub(0,"No", df_mincomp$Vibrant)
df_mincomp$Virfinder <- gsub(1,"Yes", df_mincomp$Virfinder)
df_mincomp$Virfinder <- gsub(0,"No", df_mincomp$Virfinder)
df_mincomp$Virsorter <- gsub(1,"Yes", df_mincomp$Virsorter)
df_mincomp$Virsorter <- gsub(0,"No", df_mincomp$Virsorter)
df_mincomp <- mutate_if(df_mincomp, is.numeric, as.factor)

# merge the data.frames in a single one, containing for each vOTU the links
# to the relative files
df_tot <-Reduce(function(x,y) merge(x = x, y = y, by = "ViralOTU"), 
       list(df, df_consensus, df_proteins, df_coords, df_mincomp))
df_tot <- merge(df_tot, df_graphs, by = "ViralOTU", all=TRUE)
df_tot[is.na(df_tot)] <- "n.a."
# change characters to factors
df_tot <- mutate_if(df_tot, is.character, as.factor)
# reorder and cleanup the dataframe
df_tot <- df_tot[mixedorder(as.character(df_tot$ViralOTU)),]
df_tot$BaltimoreGroup <- NULL
df_tot$Realm <- NULL
df_tot$Kingdom <- NULL
df_tot$Class <- NULL
df_tot$Order <- NULL
df_tot$Subfamily <- NULL
colnames(df_tot)[2] <- "Inherit from" 
df_tot <- df_tot[,c(1,2,8,9,10,11,19,12,13,14,15,16,17,18,3,4,5,6,7)]
rownames(df_tot) <- 1:length(rownames(df_tot))
# create the table with DT package
taxo <- DT::datatable(df_tot, filter="top", rownames = FALSE,
                      extensions = c('Buttons','Scroller'), escape = F, 
                      options = list(scrollY = "650px",
                                     scrollX = TRUE,
                                     scroller = TRUE,
                                     columnDefs = list(list(visible=FALSE,
                                                            targets=c(7,8,9,
                                                                      14,15,16,
                                                                      17,18))),
                                     dom = 'BlPfrtip',
                                     buttons = list('copy',
                                                    'print',
                                                    list(extend = 'collection',
                                                         buttons = c('csv',
                                                                     'excel',
                                                                     'pdf'),
                                                         text = 'Download as:'),
                                                    list(extend = 'colvis',
                                                         columns = c(1,2,3,4,5,
                                                                     6,10,11,12,
                                                                     13,14,15,16,
                                                                     17,18)),
                                                    list(extend = 'colvis',
                                                         text = "Fasta files",
                                                         columns = c(7,8,9))
                                     )
                      )
)
taxo$width = 1500
htmlwidgets::saveWidget(taxo, "taxonomy_table.html",
                        selfcontained = TRUE, libdir = NULL)

################################ KRAKEN FILES ##################################
# read the coordinates files and create a data.frame of paths
krakens_path <- list.files(path = kraken_path, 
                           pattern = "SRR\\d+_report.txt", 
                           recursive = TRUE,
                           full.names=T)
df_kraken_fp <- as.data.frame(krakens_path)
# match the coordinates to the vOTU in taxo table
krakens_names = unique(str_extract(krakens_path, "SRR\\d+"))
matching_krakens = krakens_names[krakens_names %in% meta$Sample]
# and to the path of the relative .gff
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
# read the coordinates files and create a data.frame of paths
kronas_path <- list.files(path = krona_path, 
                          "SRR\\d+_krak_krona_abundancies.html", 
                          recursive = TRUE,
                          full.names=T)
df_krona_fp <- as.data.frame(kronas_path)
# match the coordinates to the vOTU in taxo table
kronas_names = unique(str_extract(kronas_path, "SRR\\d+"))
matching_kronas = kronas_names[kronas_names %in% meta$Sample]
# and to the path of the relative .gff
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
