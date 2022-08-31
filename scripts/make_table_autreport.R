#Gut-Microbiome Report
setwd("C:/Users/ctata/Documents/Lab/M3/Aut-Reports")

set.seed(100)
library(ggplot2)

#functions

percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

#Go through by sibling Pairs

#Select all classes that have at least >0.1% abundance in either sample
ps <- readRDS("ps_DeSeq_pass_min_postDD_min0.03.Rda") 

#reformatting Pair IDs for matching
map <- sample_data(ps)
pair_ids <- gsub("_N", "", map$Natural.sibling.included.in.present.study)
pair_ids <- gsub("_A", "", pair_ids)
sample_data(ps)$Pair <- pair_ids
pair_ids <- unique(pair_ids)

#assigning taxonomy (phylo object is organzied only by sequence)
taxlevel_tab <- as.data.frame(tax_table(ps))
#taxa_with_levels <- assignTaxonomy(as.character(taxlevel_tab$Sequence), "~/Data/silva_nr_v132_train_set.fa.gz", multithread = TRUE, verbose=TRUE)
taxa_with_levels <- readRDS("taxawithlevels.RDS")
taxa_names(ps) <- taxlevel_tab$Sequence
#rownames(otu_table(ps)) < taxlevel_tab$Sequence
tax_table(ps)<- taxa_with_levels         

ps <- tax_glom(ps, "Class")

ps  = transform_sample_counts(ps, function(x) x / sum(x) )

#perc_otu <- otu_table(ps) / rep(colSums(otu_table(ps)), each = nrow(otu_table(ps)))
#otu_table(ps) <- otu_table(perc_otu, taxa_are_rows = T)

#for timepoint 1
ps_time1 <- subset_samples(ps, Within.study.sampling.date == "Timepoint 1")
ps_time2 <- subset_samples(ps, Within.study.sampling.date == "Timepoint 2")
ps_time3 <- subset_samples(ps, Within.study.sampling.date == "Timepoint 3")


makelist <- function(pair_id, pstime){
  print(pair_id)
  #pair_ps <- subset_samples(ps, Pair %in% pair_id) #Subset_samples breaks in function
  pair_toKeep <- sample_data(ps_time1)$Pair == pair_id
  pair_ps <- prune_samples(pair_toKeep, pstime)
  
  toKeep <- apply(data.frame(otu_table(pair_ps)), 1, function(taxa) return(any(taxa > 0.001)))
  pair_ps_filt <- prune_taxa(toKeep, pair_ps)
  
  tax_label <- paste("Kingdom:", gsub("k__", "",  tax_table(pair_ps_filt)[,1]),
                     "Phylum:", gsub("p__", "",  tax_table(pair_ps_filt)[,2]),
                     "Class:", gsub("c__", "",  tax_table(pair_ps_filt)[,3]))
  
  table_print <- otu_table(pair_ps_filt)
  table_print <- as.data.frame(table_print)
  rownames(table_print) <- tax_label
  table_print[,1] <- percent(x= table_print[,1])
  table_print[,2] <- percent(x= table_print[,2])
  table_print$Class <- rownames(table_print)
  table_print <- data.frame(table_print[,3], table_print[,1], table_print[,2], row.names = NULL)
  colnames(table_print) <- c("Bacterial Class > 0.1%","Your Child with Autism", "His/Her Sibling")
  return(table_print)
  
}
output_list_1 <- lapply(pair_ids, makelist, ps_time1)
output_list_2 <- lapply(pair_ids, makelist, ps_time2)
output_list_3 <- lapply(pair_ids, makelist, ps_time3)


makeplot <- function(pair_id, output_list){
  pair1_time1 <- as.data.frame(output_list[[pair_id]])
  pair1_time1 <- data.frame(pair1_time1[,2], pair1_time1[,3], row.names = pair1_time1[,1])
  colnames(pair1_time1) <- c("Your Child with Autism", "His/Her Sibling")
  
  test <- pair1_time1
  test <- t(test)
  test <- melt(test, id.vars = "Class" )
  test$value <- as.numeric(gsub("%", "", test$value))
  colnames(test) <- c("Person", "Class", "Percentage")
  p <- ggplot(data = test, aes(x = Person,  y = Percentage, fill = Class)) +
    geom_bar(stat = "identity")
  return(p)
}

tmp <- seq(length(output_list_1))
sib_plots_1 <- lapply(tmp, makeplot, output_list_1)

