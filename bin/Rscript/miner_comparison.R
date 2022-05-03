shhh <- suppressPackageStartupMessages
shhh(library(UpSetR))

# Reading input
args <- commandArgs(trailingOnly = TRUE)
########################## MINER COMPARISON MATRIX #############################
mc_matrix <- args[1]
df <- read.csv(mc_matrix, sep = '')

df$Cluster <- NULL
df$Header <- NULL
df$Length <- NULL
# base upSet plot
png(filename="./miner_comparison.png", width = 1400, height= 800)
upset(df, sets = c("Phigaro","Vibrant","Virfinder","Virsorter"),
      sets.bar.color = c("#0277bd","#0277bd","#0277bd","#0277bd"),
      mainbar.y.label = "vOTUs intersection", sets.x.label = "vOTUs per miner",
      order.by = "freq",
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      point.size = 3.5,
      line.size = 2, 
      text.scale = 2)
dev.off()
