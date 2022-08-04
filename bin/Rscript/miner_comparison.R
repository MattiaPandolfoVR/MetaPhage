# Coded by Mattia Pandolfo (mattia.pandolfo@univr.it)
shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(plyr))
shhh(library(UpSetR))
shhh(library(plotly))
shhh(library(reshape2))
#Reading input
args <- commandArgs(trailingOnly = TRUE)
# check arguments
if (length(args) < 3) {
  stop("Usage: miner_comparison.R <vc_min_comp> <vOTUs_min_comp> <params.minlen>\n")
}
###################### VIRAL CONTIGS MINER COMPARISON MATRIX ###################
vc_matrix <- args[1]
if (! file.exists(vc_matrix)){
  stop("ERROR: viral contigs miner comparison matrix file not found in: ",
       vc_matrix, "\n")
}
dfv <- read.delim(vc_matrix, header = TRUE, sep = "", check.names = F)
######################### vOTUs MINER COMPARISON MATRIX ########################
mc_matrix <- args[2]
if (! file.exists(mc_matrix)){
  stop("ERROR: miner comparison matrix file not found in: ", mc_matrix, "\n")
}
dfm <- read.delim(mc_matrix, header = TRUE, sep = "", check.names = F)
minlen <- args[3]
if (minlen == FALSE){
  minlen = 0
}
############################ upSetR MINER COMPARISON ###########################
# Remove all vOTUs which length is under the minlen threshold
dfm <- dfm[dfm$Length >= minlen,]
dfno0 <- dfm
dfno0$Cluster <- NULL
dfno0$Header <- NULL
dfno0$Length <- NULL
# Remove columns of phage-miners full of 0s
dfno0 <- dfno0[, !apply(dfno0 == 0, 2, all)]
# set the columns to be printed in upSetR
set = colnames(dfno0)[-1]
# and colour their bars in lightblue
colours = rep("#0277bd", length(set))

# base upSet plot
png(filename="./votus_miner_comparison_mqc.png", width = 1400, height= 800)
upset(dfno0, sets = set,
      sets.bar.color = colours,
      mainbar.y.label = "vOTUs intersection",
      sets.x.label = "vOTUs per miner",
      order.by = "freq",
      mb.ratio = c(0.6, 0.4),
      number.angles = 0,
      point.size = 3.5,
      line.size = 2,
      text.scale = 2)
dev.off()
######################## LENGTH DISTRIBUTION VIOLIN-PLOT #######################
# remove number headers
dfv_k <- dfv %>%
  filter(grepl(pattern = "^k", dfv$ContigsID))
# melt df
df_melt <- melt(dfv_k, id.vars = c("ContigsID", "Length"), na.rm = TRUE,
                variable.name = "Miner", value.name = "Presence")
df_melt_no0 <- df_melt[df_melt$Presence != 0,]
df_melt_no0$Presence[df_melt_no0$Presence > 1] <- 1
# threshold line
tresh_line <- function(y = 0, color = "red") {
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color, dash="dot")
  )
}

# Plot
vir_dist_plot <- df_melt_no0 %>%
  plot_ly(
    y = ~Length,
    hovertext = paste("viral contig: ", df_melt_no0$ContigsID),
    hoveron = "points+kde",
    type = 'violin',
    split = ~Miner,
    box = list(
      visible = T
    ),
    meanline = list(
      visible = T
    )
  )
vir_dist_plot <- vir_dist_plot %>%
  layout(
    title = "Viral contigs distribution plot<br>Red dashed line represent the contig minimum length threshold.<br>Viral contigs under this threshold are discarded after dereplication.",
    shapes = tresh_line(minlen),
    yaxis = list(
      rangemode = "nonnegative"
    )
  )
################################### PLOT GRID ##################################
htmlwidgets::saveWidget(vir_dist_plot, "./vc_distribution_plot_mqc.html",
                       selfcontained = TRUE, libdir = NULL)

############################ VC MINER COMPARISON TABLE #########################
# change characters to factors
df_dfv <- dfv_k
#Some data manipulation for better visualization
df_dfv$Deepvirfinder[df_dfv$Deepvirfinder > 1] <- 1
df_dfv$Phigaro[df_dfv$Phigaro > 1] <- 1
df_dfv$Vibrant[df_dfv$Vibrant > 1] <- 1
df_dfv$Virfinder[df_dfv$Virfinder > 1] <- 1
df_dfv$Virsorter2[df_dfv$Virsorter2 > 1] <- 1
# change miner value from numeric range to character
df_dfv$Deepvirfinder <- gsub(1,"Yes", df_dfv$Deepvirfinder)
df_dfv$Deepvirfinder <- gsub(0,"No", df_dfv$Deepvirfinder)
df_dfv$Phigaro <- gsub(1,"Yes", df_dfv$Phigaro)
df_dfv$Phigaro <- gsub(0,"No", df_dfv$Phigaro)
df_dfv$Vibrant <- gsub(1,"Yes", df_dfv$Vibrant)
df_dfv$Vibrant <- gsub(0,"No", df_dfv$Vibrant)
df_dfv$Virfinder <- gsub(1,"Yes", df_dfv$Virfinder)
df_dfv$Virfinder <- gsub(0,"No", df_dfv$Virfinder)
df_dfv$Virsorter2 <- gsub(1,"Yes", df_dfv$Virsorter2)
df_dfv$Virsorter2 <- gsub(0,"No", df_dfv$Virsorter2)
# change characters to factors
df_dfv <- mutate_if(df_dfv, is.character, as.factor)
# rename contigs_id column to vOTU
colnames(df_dfv)[1] <- "ViralContig"
# create the table with DT package
vc_t <- DT::datatable(df_dfv, filter="top", rownames = FALSE,
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
)
vc_t$width = 1500
htmlwidgets::saveWidget(vc_t, "./vc_miner_comparison_table.html",
                        selfcontained = TRUE, libdir = NULL)

