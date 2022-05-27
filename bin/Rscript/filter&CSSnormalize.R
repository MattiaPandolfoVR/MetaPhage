filter_CSSnormalize <- function(ps){
  ################################## FILTERING ###################################
  # Filter uncharacterized taxas
  ps_filter <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "Unclassified"))
  # relative abundance filter
  FSr  = transform_sample_counts(ps_filter, function(x) x / sum(x))
  FSfr = filter_taxa(FSr, function(x) mean(x) > 0.0001, TRUE)
  keeptaxa = taxa_names(FSfr)
  alltaxa = taxa_names(ps_filter)
  myTaxa = alltaxa[alltaxa %in% keeptaxa]
  physeqaFS <- prune_taxa(myTaxa, ps_filter)
  ps_filter = filter_taxa(physeqaFS, function(x) sum(x >= 1) >= (2), TRUE)
  
  ################################ NORMALIZATION #################################
  # CSS
  ps_m <- phyloseq::phyloseq_to_metagenomeSeq(ps_filter)
  # normalized count matrix
  ps_norm <- metagenomeSeq::MRcounts(ps_m, norm = TRUE, log = TRUE)
  # abundance to css-normalized
  phyloseq::otu_table(ps_filter) <- phyloseq::otu_table(ps_norm, taxa_are_rows = T)
  # Restore sample_data rownames
  row.names(ps_filter@sam_data) <- c(1:nrow(ps_filter@sam_data))
  return(ps_filter)
}