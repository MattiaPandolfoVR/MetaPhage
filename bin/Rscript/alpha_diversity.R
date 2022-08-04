# Coded by Mattia Pandolfo (mattia.pandolfo@univr.it)
shhh <- suppressPackageStartupMessages
shhh(library(dplyr))
shhh(library(plyr))
shhh(library(plotly))
shhh(library(phyloseq))
shhh(library(data.table))
shhh(library(reshape2))

# Reading input
args <- commandArgs(trailingOnly = TRUE)
# check arguments
if (length(args) < 5) {
  stop("Usage: alpha_diversity.R <count_table> <taxonomy_table> <metadata> <alpha_var1> <alpha_var2>\n")
}
################################## COUNT TABLE #################################
file_count <- args[1]
if (! file.exists(file_count)){
  stop("ERROR: count table file not found in: ", file_count, "\n")
}
count <- read.delim(file_count, row.names = 1, sep = "\t", check.names = F)
################################ TAXONOMY TABLE ################################
file_taxo <- args[2]
if (! file.exists(file_taxo)){
  stop("ERROR: taxonomy table file not found in: ", file_taxo, "\n")
}
# change "O", "n.a." and "uncharacterize" to "NA" in taxonomy table
taxo <- read.delim(file_taxo, row.names = 1, sep = "\t", check.names = F, na.strings = c("n.a.", "O", "Unclassified"))
# Extract the scaffold id and the last 2 column from taxonomy table
taxo <- taxo[,c("Host","Family","Genus")]
################################### METADATA ###################################
file_meta <- args[3]
if (! file.exists(file_meta)){
  stop("ERROR: metadata file not found in: ", file_meta, "\n")
}
metadata <- read.delim(file_meta, row.names = 1, sep = ",", check.names = F)
metadata$Sample <- rownames(metadata)
metadata <- metadata %>%
  select(Sample, everything())
################################ ALPHA VARIABLES ###############################
alpha_var1 <- args[4]
alpha_var1_c = 1
alpha_var2 <- args[5]
alpha_var2_c = 1
# Check if only one variable is passed
if(alpha_var1 == FALSE && alpha_var2 == FALSE){
  stop("ERROR: alpha_var1 and alpha_var2 parameters not found. Please set at lease one of the two in your project config.file (e.g. alpha_var1 = \"metadata_variable1\" \n")
} else if(alpha_var1 == FALSE) {
  alpha_var1_c = FALSE
  alpha_var1 <- alpha_var2
} else if(alpha_var2 == FALSE) {
  alpha_var2_c = FALSE
  alpha_var2 <- alpha_var1
}

################################## PHYLOSEQ ####################################
ps <- phyloseq(otu_table(count, taxa_are_rows = TRUE),
               tax_table(as.matrix(taxo)),
               sample_data(metadata))
# Remove NAs at phylum level
# Restore sample_data rownames
row.names(ps@sam_data) <- c(1:nrow(ps@sam_data))
# Calculate alpha_diversity measurments
df_div <- estimate_richness(ps)
row.names(df_div) <- ps@sam_data$Sample
df_div$Sample <- row.names(df_div)
df_div <- df_div %>%
  select(Sample, everything())
row.names(df_div) <- c(1:nrow(df_div))
# Transform ps@sample_data in a df
sample_df <- as.matrix(ps@sam_data)
sample_df <- as.data.frame(sample_df)
# Merge the two data-frames
df_alpha <- merge(sample_df, df_div, by="Sample")
df_alpha_min <- subset(x = df_alpha, select = c(Sample, get(alpha_var1), get(alpha_var2), Observed,
                                                Chao1, se.chao1, ACE, se.ACE, Shannon, Simpson, InvSimpson, Fisher))
df_alpha_min_dt <- setDT(df_alpha_min)

# create item lists (for dichotomous and not dichotomous variables)
items_dicho <- list(
  list(label="Observed", method = "update", args=list(list(visible=c(T,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)),
                                                      list(title = "Alpha diversity<br>Observed"))),
  list(label="Chao1", method = "update", args=list(list(visible=c(F,F,T,T,F,F,F,F,F,F,F,F,F,F,F,F,F,F)),
                                                   list(title = "Alpha diversity<br>Chao1"))),
  list(label="se.Chao1", method = "update", args=list(list(visible=c(F,F,F,F,T,T,F,F,F,F,F,F,F,F,F,F,F,F)),
                                                      list(title = "Alpha diversity<br>se.Chao1"))),
  list(label="ACE", method = "update", args=list(list(visible=c(F,F,F,F,F,F,T,T,F,F,F,F,F,F,F,F,F,F)),
                                                 list(title = "Alpha diversity<br>ACE"))),
  list(label="se.ACE", method = "update", args=list(list(visible=c(F,F,F,F,F,F,F,F,T,T,F,F,F,F,F,F,F,F)),
                                                    list(title = "Alpha diversity<br>se.ACE"))),
  list(label="Shannon", method = "update", args=list(list(visible=c(F,F,F,F,F,F,F,F,F,F,T,T,F,F,F,F,F,F)),
                                                     list(title = "Alpha diversity<br>Shannon"))),
  list(label="Simpson", method = "update", args=list(list(visible=c(F,F,F,F,F,F,F,F,F,F,F,F,T,T,F,F,F,F)),
                                                     list(title = "Alpha diversity<br>Simpson"))),
  list(label="InvSimpson", method = "update", args=list(list(visible=c(F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,T,F,F)),
                                                        list(title = "Alpha diversity<br>InvSimpson"))),
  list(label="Fisher", method = "update", args=list(list(visible=c(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,T,T)),
                                                    list(title = "Alpha diversity<br>Fisher")))
)

items_nodicho <- list(
  list(label="Observed", method = "update", args=list(list(visible=c(T,F,F,F,F,F,F,F,F)),
                                                      list(title = "Alpha diversity<br>Observed"))),
  list(label="Chao1", method = "update", args=list(list(visible=c(F,T,F,F,F,F,F,F,F)),
                                                   list(title = "Alpha diversity<br>Chao1"))),
  list(label="se.Chao1", method = "update", args=list(list(visible=c(F,F,T,F,F,F,F,F,F)),
                                                      list(title = "Alpha diversity<br>se.Chao1"))),
  list(label="ACE", method = "update", args=list(list(visible=c(F,F,F,T,F,F,F,F,F)),
                                                 list(title = "Alpha diversity<br>ACE"))),
  list(label="se.ACE", method = "update", args=list(list(visible=c(F,F,F,F,T,F,F,F,F)),
                                                    list(title = "Alpha diversity<br>se.ACE"))),
  list(label="Shannon", method = "update", args=list(list(visible=c(F,F,F,F,F,T,F,F,F)),
                                                     list(title = "Alpha diversity<br>Shannon"))),
  list(label="Simpson", method = "update", args=list(list(visible=c(F,F,F,F,F,F,T,F,F)),
                                                     list(title = "Alpha diversity<br>Simpson"))),
  list(label="InvSimpson", method = "update", args=list(list(visible=c(F,F,F,F,F,F,F,T,F)),
                                                        list(title = "Alpha diversity<br>InvSimpson"))),
  list(label="Fisher", method = "update", args=list(list(visible=c(F,F,F,F,F,F,F,F,T)),
                                                    list(title = "Alpha diversity<br>Fisher")))
)

# if one of the alpha_var is a dichotomous variable, use it to color the plot
if(alpha_var1_c == 1 && alpha_var2_c == 1){
  if(nrow(unique(df_alpha_min_dt[,2])) == 2){ #alpha_var1 is dichotomous
    # plot
    fig_alpha <- plot_ly(data=df_alpha_min_dt, x = ~get(alpha_var2), y = ~Observed, color = ~get(alpha_var1),
                         colors = "Dark2",  visible=T, type = "box", boxpoints = "all", groupclick = "toggleitem",
                         hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~Chao1, color = ~get(alpha_var1),
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~se.chao1, color = ~get(alpha_var1),
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~ACE, color = ~get(alpha_var1),
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$SampleE)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~se.ACE, color = ~get(alpha_var1),
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~Shannon, color = ~get(alpha_var1),
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~Simpson, color = ~get(alpha_var1),
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~InvSimpson, color = ~get(alpha_var1),
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~Fisher, color = ~get(alpha_var1),
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      layout(
        title = "Alpha diversity<br>Select the index from the dropdown menu",
        xaxis = list(title= paste(alpha_var1)),
        yaxis = list(title = "Species Richness"),
        legend= list(title=list(text= alpha_var1),
                     traceorder= "normal",
                     itemsizing= "constant"),
        autosize = TRUE,
        boxmode = "group",
        margin = list(autosize = TRUE, width = 1400, height = 900),
        updatemenus = list(
          list(
            buttons = items_dicho)
        )
      )
    htmlwidgets::saveWidget(fig_alpha, "./alpha_diversities.html",
                            selfcontained = TRUE, libdir = NULL)
    } else if(nrow(unique(df_alpha_min_dt[,3])) == 2){ #alpha_var2 is dichotomous
      # plot
      fig_alpha <- plot_ly(data=df_alpha_min_dt, x = ~get(alpha_var1), y = ~Observed, color = ~get(alpha_var2),
                           colors = "Dark2",  visible=T, type = "box", boxpoints = "all", groupclick = "toggleitem",
                           hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
        add_boxplot(x = ~get(alpha_var1), y = ~Chao1, color = ~get(alpha_var2),
                    colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
        add_boxplot(x = ~get(alpha_var1), y = ~se.chao1, color = ~get(alpha_var2),
                    colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
        add_boxplot(x = ~get(alpha_var1), y = ~ACE, color = ~get(alpha_var2),
                    colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$SampleE)) %>%
        add_boxplot(x = ~get(alpha_var1), y = ~se.ACE, color = ~get(alpha_var2),
                    colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
        add_boxplot(x = ~get(alpha_var1), y = ~Shannon, color = ~get(alpha_var2),
                    colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
        add_boxplot(x = ~get(alpha_var1), y = ~Simpson, color = ~get(alpha_var2),
                    colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
        add_boxplot(x = ~get(alpha_var1), y = ~InvSimpson, color = ~get(alpha_var2),
                    colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
        add_boxplot(x = ~get(alpha_var1), y = ~Fisher, color = ~get(alpha_var2),
                    colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
        layout(
          title = "Alpha diversity<br>Select the index from the dropdown menu",
          xaxis = list(title= paste(alpha_var1)),
          yaxis = list(title = "Species Richness"),
          legend= list(title=list(text= alpha_var2),
                       traceorder= "normal",
                       itemsizing= "constant"),
          autosize = TRUE,
          boxmode = "group",
          margin = list(autosize = TRUE, width = 1400, height = 900),
          updatemenus = list(
            list(
              buttons = items_dicho)
          )
        )
      htmlwidgets::saveWidget(fig_alpha, "./alpha_diversities.html",
                              selfcontained = TRUE, libdir = NULL)
      } else { #alpha_var1 nor alpha_var2 are dichotomous
        # plot
        fig_alpha <- plot_ly(data=df_alpha_min_dt, x = ~get(alpha_var1), y = ~Observed, color = "#66C2A5",
                             colors = "Dark2",  visible=T, type = "box", boxpoints = "all", groupclick = "toggleitem",
                             hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
          add_boxplot(x = ~get(alpha_var1), y = ~Chao1, color = "#FC8D62",
                      colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
          add_boxplot(x = ~get(alpha_var1), y = ~se.chao1, color = "#8DA0CB",
                      colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
          add_boxplot(x = ~get(alpha_var1), y = ~ACE, color = "#E78AC3",
                      colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$SampleE)) %>%
          add_boxplot(x = ~get(alpha_var1), y = ~se.ACE, color = "#A6D854",
                      colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
          add_boxplot(x = ~get(alpha_var1), y = ~Shannon, color = "#FFD92F",
                      colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
          add_boxplot(x = ~get(alpha_var1), y = ~Simpson, color = "#E5C494",
                      colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
          add_boxplot(x = ~get(alpha_var1), y = ~InvSimpson, color = "#B3B3B3",
                      colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
          add_boxplot(x = ~get(alpha_var1), y = ~Fisher, color = "#FFD92F",
                      colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
          layout(
            title = "Alpha diversity<br>Select the index from the dropdown menu",
            xaxis = list(title= paste(alpha_var1)),
            yaxis = list(title = "Species Richness"),
            showlegend = FALSE,
            autosize = TRUE,
            boxmode = "group",
            margin = list(autosize = TRUE, width = 1400, height = 900),
            updatemenus = list(
              list(
                buttons = items_nodicho)
            )
          )
        htmlwidgets::saveWidget(fig_alpha, "./alpha_diversities.html",
                                selfcontained = TRUE, libdir = NULL)
      }
} else if(alpha_var1_c == 1 && alpha_var2_c != 1){
  # plot
  fig_alpha <- plot_ly(data=df_alpha_min_dt, x = ~get(alpha_var1), y = ~Observed, color = "#66C2A5",
                       colors = "Dark2",  visible=T, type = "box", boxpoints = "all", groupclick = "toggleitem",
                       hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
    add_boxplot(x = ~get(alpha_var1), y = ~Chao1, color = "#FC8D62",
                colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
    add_boxplot(x = ~get(alpha_var1), y = ~se.chao1, color = "#8DA0CB",
                colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
    add_boxplot(x = ~get(alpha_var1), y = ~ACE, color = "#E78AC3",
                colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$SampleE)) %>%
    add_boxplot(x = ~get(alpha_var1), y = ~se.ACE, color = "#A6D854",
                colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
    add_boxplot(x = ~get(alpha_var1), y = ~Shannon, color = "#FFD92F",
                colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
    add_boxplot(x = ~get(alpha_var1), y = ~Simpson, color = "#E5C494",
                colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
    add_boxplot(x = ~get(alpha_var1), y = ~InvSimpson, color = "#B3B3B3",
                colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
    add_boxplot(x = ~get(alpha_var1), y = ~Fisher, color = "#FFD92F",
                colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
    layout(
      title = "Alpha diversity<br>Select the index from the dropdown menu",
      xaxis = list(title= paste(alpha_var1)),
      yaxis = list(title = "Species Richness"),
      showlegend = FALSE,
      autosize = TRUE,
      boxmode = "group",
      margin = list(autosize = TRUE, width = 1400, height = 900),
      updatemenus = list(
        list(
          buttons = items_nodicho)
      )
    )
  htmlwidgets::saveWidget(fig_alpha, "./alpha_diversities.html",
                          selfcontained = TRUE, libdir = NULL)
  } else if(alpha_var1_c != 1 && alpha_var2_c == 1){
    # plot
    fig_alpha <- plot_ly(data=df_alpha_min_dt, x = ~get(alpha_var2), y = ~Observed, color = "#66C2A5",
                         colors = "Dark2",  visible=T, type = "box", boxpoints = "all", groupclick = "toggleitem",
                         hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~Chao1, color = "#FC8D62",
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~se.chao1, color = "#8DA0CB",
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~ACE, color = "#E78AC3",
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$SampleE)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~se.ACE, color = "#A6D854",
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~Shannon, color = "#FFD92F",
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~Simpson, color = "#E5C494",
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~InvSimpson, color = "#B3B3B3",
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      add_boxplot(x = ~get(alpha_var2), y = ~Fisher, color = "#FFD92F",
                  colors = "Dark2",  visible=F, hovertext = paste("Sample: ", df_alpha_min_dt$Sample)) %>%
      layout(
        title = "Alpha diversity<br>Select the index from the dropdown menu",
        xaxis = list(title= paste(alpha_var2)),
        yaxis = list(title = "Species Richness"),
        showlegend = FALSE,
        autosize = TRUE,
        boxmode = "group",
        margin = list(autosize = TRUE, width = 1400, height = 900),
        updatemenus = list(
          list(
            buttons = items_nodicho)
        ))
    htmlwidgets::saveWidget(fig_alpha, "./alpha_diversities.html",
                            selfcontained = TRUE, libdir = NULL)
    }