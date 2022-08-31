setwd("C:/Users/kitikomp/Documents/Lab/M3")
source('scripts/clean_mapping_ml.R')
ps_16s <- readRDS("data/16s/ps_dds_nooutliers.rds")
ps_mtg <- readRDS("data/mtg/ps_rle_nooutliers.rds")
ps_mtt <- readRDS("data/mtt/ps_rle_nooutliers.rds")
ps_metabol <- readRDS("data/metabol/ps_rle_nooutliers.rds")

dim(ps_mtg@otu_table)


ps1 <- ps_mtg
ps2 <- ps_mtt
sample_set <- intersect(sample_names(ps1), sample_names(ps2))
print(paste("Number of samples overlapping: ", length(sample_set)))
ps1 <- prune_samples(sample_set, ps1)
ps2 <- prune_samples(sample_set, ps2)

write.table(ps1@otu_table, "results/halla/mtg_tab.txt", sep = "\t", quote = F)
write.table(ps2@otu_table, "results/halla/mtt_tab.txt", sep = "\t", quote = F)
