# Coded by Mattia Pandolfo (mattia.pandolfo@univr.it)
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
# check arguments
if (length(args) < 7) {
  stop("Usage: taxonomy_table.R <taxonomy_table> <vOTUs_consensus> <vOTUs_proteins> <vOTUs_coords> <outdir_path> <min_comp> <metadata>\n")
}
################################ TAXONOMY TABLE ################################
file_taxo <- args[1]
if (! file.exists(file_taxo)){
  stop("ERROR: taxonomy table file not found in: ", file_taxo, "\n")
}
taxo <- read.delim(file_taxo, sep = "\t", check.names = F)
# mutate the taxonomy data.frame
taxo <- mutate_if(taxo, is.character, as.factor)
# ########################## .fna, .faa, .gff FILES ##############################
file_fna <- args[2]
if (! file.exists(file_fna)){
  stop("ERROR: .fna file not found in: ", file_fna, "\n")
}
fna <- read.fasta(file_fna, set.attributes = FALSE, as.string=TRUE,
                  forceDNAtolower = FALSE)
file_faa <- args[3]
if (! file.exists(file_faa)){
  stop("ERROR: .faa file not found in: ", file_faa, "\n")
}
faa <- read.fasta(file_faa, set.attributes = FALSE, as.string=TRUE,
                  forceDNAtolower = FALSE)
file_gff <- args[4]
if (! file.exists(file_gff)){
  stop("ERROR: .gff file not found in: ", file_gff, "\n")
}
gff <- read.delim(file_gff, header=F, comment.char="#")
#################################### PATHS #####################################
file_paths <- args[5]
fna_path <- file.path(file_paths, "report/files/vOTUs_consensus")
if (! file.exists(fna_path)){
  cat("WARNING: report/files/vOTUs_consensus folder not found in: ", fna_path, "\n")
}
faa_path <- file.path(file_paths, "report/files/vOTUs_proteins")
if (! file.exists(faa_path)){
  cat("WARNING: report/files/vOTUs_proteins folder not found in: ", faa_path, "\n")
}
gff_path <- file.path(file_paths, "report/files/vOTUs_coords")
if (! file.exists(gff_path)){
  cat("WARNING: report/files/vOTUs_coords folder not found in: ", gff_path, "\n")
}
graph_path <- file.path(file_paths, "report/taxonomy/vcontact2/single-views_megahit")
if (! file.exists(graph_path)){
  cat("WARNING: report/taxonomy/vcontact2/single-views_megahit folder not found in: ", graph_path, "\n")
}
############################### MINER COMPARISON ###############################
file_min <- args[6]
if (! file.exists(file_min)){
  stop("ERROR: miner comparison file not found in: ", file_min, "\n")
}
min_com <- read.delim(file_min, sep = "", check.names = F)
################################### METADATA ###################################
file_meta <- args[7]
if (! file.exists(file_meta)){
  stop("ERROR: metadata file not found in: ", file_meta, "\n")
}
metadata <- read.delim(file_meta, row.names=1, sep = ",", check.names = F)
metadata$Sample <- rownames(metadata)
metadata <- metadata %>%
  select(Sample, everything())

############################## CONSENSUS FILES ################################
# read the consensus fasta files and create a data.frame of paths
consensus_path <- list.files(path = fna_path,
                               pattern = "vOTU*",
                               full.names=T)
df_fasta_fp <- as.data.frame(consensus_path)
# match the consensus to the vOTU in taxo table
matching_cons = names(fna)[names(fna) %in% taxo$ViralOTU]
# and to the path of the relative .fasta
match_cons_fp =  df_fasta_fp[match(matching_cons,
                           str_extract(consensus_path,
                                       "vOTU_\\d+")),]
# modify the matching fasta file path with a relative path
consensus_fp = paste0("./", str_extract(match_cons_fp,
                              "files/vOTUs_consensus/vOTU_\\d+.fasta"))
# create the link table
df_consensus = tibble("ViralOTU" = matching_cons, file= consensus_fp) %>%
  mutate(file = str_replace_all(file,
                                '([^;]*)files/vOTUs_consensus/', ''),
         path_fasta = file.path(consensus_fp),
         Consensus = paste0('<a target=_blank href=',
                            path_fasta, '>', file,'</a>'))
df_consensus$file <- NULL
df_consensus$path_fasta <- NULL

############################### PROTEIN FILES #################################
# read the protein files and create a data.frame of paths
proteins_path <- list.files(path = faa_path,
                              pattern = "vOTU*",
                              full.names=T)
df_prot_fp <- as.data.frame(proteins_path)
# match the proteins to the vOTU in taxo table
prot_names = unique(str_extract(names(faa), "vOTU_[0-9]+"))
matching_prots = prot_names[prot_names %in% taxo$ViralOTU]
# and to the path of the relative .faa
matching_prots_fp =  df_prot_fp[match(matching_prots,
                                     str_extract(proteins_path, "vOTU_\\d+")),]
# modify the matching .faa file path with a relative path
protein_fp = paste0("./", str_extract(matching_prots_fp,
                                "files/vOTUs_proteins/vOTU_\\d+.faa"))
# create the link table
df_proteins = tibble("ViralOTU" = matching_prots, file = protein_fp) %>%
  mutate(file = str_replace_all(file,
                                '([^;]*)files/vOTUs_proteins/', ''),
         path_prot = file.path(protein_fp),
         Proteins = paste0('<a target=_blank href=',
                           path_prot, '>', file,'</a>'))
df_proteins$file <- NULL
df_proteins$path_prot <- NULL

################################ GFF FILES ####################################
# read the coordinates files and create a data.frame of paths
coords_path <- list.files(path = gff_path,
                              pattern = "vOTU*",
                              full.names=T)
df_coord_fp <- as.data.frame(coords_path)
# match the coordinates to the vOTU in taxo table
coord_names = unique(gff$V1)
matching_coords = coord_names[coord_names %in% taxo$ViralOTU]
# and to the path of the relative .gff
matching_coords_fp =  df_coord_fp[match(matching_coords,
                                     str_extract(coords_path, "vOTU_\\d+")),]
# modify the matching .gff file path with a relative path
coords_fp = paste0("./", str_extract(matching_coords_fp,
                                "files/vOTUs_coords/vOTU_\\d+.gff"))
# create the link table
df_coords = tibble("ViralOTU" = matching_coords, file = coords_fp) %>%
  mutate(file = str_replace_all(file,
                                '([^;]*)files/vOTUs_coords/', ''),
         path_coords = file.path(coords_fp),
         Coordinates = paste0('<a target=_blank href=',
                              path_coords, '>', file,'</a>'))
df_coords$file <- NULL
df_coords$path_coords <- NULL

############################### GRAPHS FILES ##################################
# read the coordinates files and create a data.frame of paths
graphs_path <- list.files(path = graph_path,
                          pattern = "vOTU*",
                          full.names=T)
df_graphs_fp <- as.data.frame(graphs_path)
# match the coordinates to the vOTU in taxo table
graphs_names = unique(str_extract(graphs_path, "vOTU_[0-9]+"))
matching_graphs = graphs_names[graphs_names %in% taxo$ViralOTU]
# and to the path of the relative .gff
graphs_names_fp =  df_graphs_fp[match(matching_graphs,
                                        str_extract(graphs_path, "vOTU_\\d+")),]
# modify the matching .gff file path with a relative path
graphs_fp = paste0("./", str_extract(graphs_names_fp,
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
df_mincomp <- merge(taxo[,c("ViralOTU","Status")],
                    min_com[, c("ViralOTU",
                                "Deepvirfinder",
                                "Phigaro",
                                "Vibrant",
                                "Virfinder",
                                "Virsorter2")],
                    by= "ViralOTU")
df_mincomp$Status <- NULL
df_mincomp$Deepvirfinder <- gsub(1,"Yes", df_mincomp$Deepvirfinder)
df_mincomp$Deepvirfinder <- gsub(0,"No", df_mincomp$Deepvirfinder)
df_mincomp$Phigaro <- gsub(1,"Yes", df_mincomp$Phigaro)
df_mincomp$Phigaro <- gsub(0,"No", df_mincomp$Phigaro)
df_mincomp$Vibrant <- gsub(1,"Yes", df_mincomp$Vibrant)
df_mincomp$Vibrant <- gsub(0,"No", df_mincomp$Vibrant)
df_mincomp$Virfinder <- gsub(1,"Yes", df_mincomp$Virfinder)
df_mincomp$Virfinder <- gsub(0,"No", df_mincomp$Virfinder)
df_mincomp$Virsorter2 <- gsub(1,"Yes", df_mincomp$Virsorter2)
df_mincomp$Virsorter2 <- gsub(0,"No", df_mincomp$Virsorter2)
df_mincomp <- mutate_if(df_mincomp, is.numeric, as.factor)

# merge the data.frames in a single one, containing for each vOTU the links
# to the relative files
df_tot <-Reduce(function(x,y) merge(x = x, y = y, by = "ViralOTU"), 
       list(taxo, df_consensus, df_proteins, df_coords, df_mincomp))
colnames(df_tot)[2] <- "Inherited_from"
df_tot <- merge(df_tot, df_graphs, by = "ViralOTU", all=TRUE)
df_tot[is.na(df_tot)] <- "n.a."
# change characters to factors
df_tot <- mutate_if(df_tot, is.character, as.factor)
# reorder and cleanup the dataframe
df_tot <- df_tot[mixedorder(as.character(df_tot$ViralOTU)),]
df_tot_ord <-  subset(x = df_tot, select = c(ViralOTU, Inherited_from, Host,
                                               Phylum, Family, Genus, Neighbors,
                                               Consensus, Proteins, Coordinates,
                                               Accession, Status, VC, Level, Weight,
                                               Deepvirfinder, Phigaro, Vibrant,
                                               Virfinder, Virsorter2
                                               ))
rownames(df_tot_ord) <- 1:length(rownames(df_tot_ord))
# create the table with DT package
taxo <- DT::datatable(df_tot_ord, filter="top", rownames = FALSE,
                      extensions = c('Buttons','Scroller'), escape = F, 
                      options = list(scrollY = "650px",
                                     scrollX = TRUE,
                                     scroller = TRUE,
                                     columnDefs = list(list(visible=FALSE,
                                                            targets=c(7,8,9,
                                                                      10,11,12))),
                                     dom = 'BlPfrtip',
                                     buttons = list('copy','print',
                                                    list(extend = 'collection',
                                                         buttons = c('csv',
                                                                     'excel',
                                                                     'pdf'),
                                                         text = 'Download as:'),
                                                    list(extend = 'colvis',
                                                         columns = c(1,2,3,4,5,
                                                                     6,7,8,9,10,11,12,
                                                                     13,14,15,16,17,18,19)),
                                                    list(extend = 'colvis',
                                                         text = "Fasta files",
                                                         columns = c(7,8,9))
                                     )
                      )
)
taxo$width = 1500
htmlwidgets::saveWidget(taxo, "./taxonomy_table.html",
                        selfcontained = TRUE, libdir = NULL)
