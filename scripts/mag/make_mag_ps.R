library(phyloseq)

setwd("C:/Users/kitikomp/Documents/Lab/M3")
tab <- read.delim("data/mtg/contig_abund_table.csv", sep=",", row.names = 1)
colnames(tab) <- gsub("X", "", colnames(tab))
mapping <- read.delim("data/mapping.csv", sep= ",", row.names = 1)

ps <- phyloseq(otu_table(tab, taxa_are_rows = T), sample_data(mapping))
saveRDS(ps, "data/mtg/ps_mag_age_filtered.rds")
