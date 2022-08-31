library(vegan)
setwd("C:/Users/kitikomp/Documents/Lab/M3/")

#Goal: remove outlying samples based on PCoA plots using the chosen distance metrics and normalizations
#Goal: remove genes that are not present in greater than 10% of samples


lookForOutliers <- function(ps, dist_metric){
  dist <- vegdist(t(ps@otu_table), method = dist_metric, na.rm = F)
  ord <- ordinate(ps, distance = dist, method = "PCoA")
  plot_ordination(ps, ord, type = "samples", color = "phenotype", label = "host_name")
}

ps_16s <- readRDS("data/16s/ps_dds.rds")
lookForOutliers(ps_16s, dist_metric = "bray")
#outliers <- c("075_N")
#ps_16s <- prune_samples(!(ps_16s@sam_data$host_name %in% c(outliers)), ps_16s)
#lookForOutliers(ps_16s, dist_metric = "bray")
#saveRDS(ps_16s, "data/16s/ps_dds_nooutliers.rds")
write.table(ps_16s@otu_table, "data/16s/asv_table_nooutliers.txt", sep = "\t", quote = F)
write.table(ps_16s@sam_data, "data/16s/mapping_nooutliers.txt", sep = "\t", quote = F)


ps_archae <- readRDS("data/archae/ps_dds.rds")
lookForOutliers(ps_archae, dist_metric = "horn")


ps_mtg <- readRDS("data/mtg/ps_rle.rds")
lookForOutliers(ps_mtg, dist_metric = "bray")
outliers <- c("038_A")
ps_mtg <- prune_samples(!(ps_mtg@sam_data$host_name %in% c(outliers)), ps_mtg)
lookForOutliers(ps_mtg, dist_metric = "manhattan")
#ps_mtg <- prune_taxa(apply(ps_mtg@otu_table, 1, function(x) return(sum(x > 0) > nsamples(ps_mtg) * .1)), ps_mtg)
saveRDS(ps_mtg, "data/mtg/ps_mtg_rle_nooutliers.rds")

# adju counts 
counts <- asinh(ps_mtg@otu_table) * 10
ps_mtg <- phyloseq(otu_table(counts, taxa_are_rows = T), sample_data(ps_mtg@sam_data), tax_table(ps_mtg@tax_table))
saveRDS(ps_mtg, "data/mtg/ps_mtg_rle_nooutliers_adjcounts.rds")

# write table
write.table(ps_mtg@otu_table, "data/mtg/asv_table_nooutliers.txt", sep = "\t", quote = F)
write.table(ps_mtg@sam_data, "data/mtg/mapping_nooutliers.txt", sep = "\t", quote = F)

ps_mtt <- readRDS("data/mtt/ps_rle.rds")
lookForOutliers(ps_mtt, dist_metric = "manhattan")
#outliers <- c("038_A")
#ps_mtt <-  prune_samples(!(ps_mtt@sam_data$host_name %in% c(outliers)), ps_mtt)
#lookForOutliers(ps_mtt, dist_metric = "bray")

saveRDS(ps_mtt, "data/mtt/ps_mtt_rle_nooutliers.rds")

# adj counts 
counts <- asinh(ps_mtt@otu_table) * 10
ps_mtt <- phyloseq(otu_table(counts, taxa_are_rows = T), sample_data(ps_mtt@sam_data), tax_table(ps_mtt@tax_table))
saveRDS(ps_mtt, "data/mtt/ps_mtt_rle_nooutliers_adjcounts.rds")

write.table(ps_mtt@otu_table, "data/mtt/asv_table_nooutliers.txt", sep = "\t", quote = F)
write.table(ps_mtt@sam_data, "data/mtt/mapping_nooutliers.txt", sep = "\t", quote = F)

ps_metabol <- readRDS("data/mbx/ps_rle.rds")
lookForOutliers(ps_metabol, dist_metric = "manhattan")
outliers <- c("177_A", "038_A", "134_N") #144_A and 134_A are questionable, I'm keeping them for now as a judgement call
ps_metabol <- prune_samples(!(ps_metabol@sam_data$host_name %in% c(outliers)), ps_metabol)
lookForOutliers(ps_metabol, dist_metric = "manhattan")
taxa_names(ps_metabol) <- data.frame(ps_metabol@tax_table)$BIOCHEMICAL
hist(apply(ps_metabol@otu_table, 2, max), breaks = 20)
saveRDS(ps_metabol, "data/mbx/ps_mbx_rle_nooutliers.rds")


# We need lower counts to be able to do topic models
counts <- asinh(ps_metabol@otu_table) * 10
ps_mbx <- phyloseq(otu_table(counts, taxa_are_rows = T), sample_data(ps_metabol@sam_data), tax_table(ps_metabol@tax_table))
saveRDS(ps_mbx, "data/mbx/ps_mbx_rle_nooutliers_adjcounts.rds")

#ps_metabol@otu_table <- otu_table(apply(counts, 1, as.integer), taxa_are_rows = T)
#taxa_names(ps_metabol) <- data.frame(ps_metabol@tax_table)$BIOCHEMICAL
#saveRDS(ps_metabol, "data/metabol/ps_mbx_rle_nooutliers_adj_counts.rds")
write.table(ps_metabol@otu_table, "data/metabol/asv_table_nooutliers.txt", sep = "\t", quote = F)
write.table(ps_metabol@sam_data, "data/metabol/mapping_nooutliers.txt", sep = "\t", quote = F)
write.table(ps_metabol@tax_table, "data/metabol/metabolite_nooutliers.txt", sep = "\t", quote = F, col.names = NA)


# Fix mapping file metabolomics
ps_16s <- readRDS("data/16s/ps_dds.rds")
ps_mtg <- readRDS("data/mtg/ps_rle.rds")
ps_mtt <- readRDS("data/mtt/ps_rle.rds")
ps_metabol <- readRDS("data/mbx/ps_mbx_rle_nooutliers_adjcounts.rds")

library(plyr)
map_16s <- ps_16s@sam_data
map_16s$omic <- "16s"
map_mtg <- ps_mtg@sam_data
map_mtg$omic <- "mtg"
map_mtt <- ps_mtt@sam_data
map_mtt$omic <- "mtt"
map_metabol <- ps_metabol@sam_data
map_metabol$omic <- "metabol"

map_master <- data.frame(rbind.fill(map_16s, map_mtg, map_mtt, map_metabol))
map_master$biospecimen_name <- gsub("F1", "P1", map_master$biospecimen_name)
map_master_one_host <- map_master[!duplicated(map_master$host_name), ]

dim(map_master_one_host)




map <- ps_metabol@sam_data
rownames(map) <- map$biospecimen_name
rownames(map) <- gsub("F1", "P1", rownames(map))

map_master <- 
map_tmp <- map_master[map_master$biospecimen_name %in% rownames(map), ]
map_tmp <- map_tmp[!duplicated(map_tmp$biospecimen_name), ]
rownames(map_tmp) <- map_tmp$biospecimen_name
sample_data(ps_metabol) <- sample_data(map_tmp[rownames(map), ])


# add kegg id
ps_tmp <- readRDS("data/mbx/ps_clean.rds")
df <- data.frame(tax_table(ps_tmp))
kegg_ids <- df[taxa_names(ps_metabol), "KEGG"]
tmp <- data.frame(ps_metabol@tax_table)
tmp$KEGG <- kegg_ids
ps_metabol@tax_table <- tax_table(as.matrix(tmp))

saveRDS(ps_metabol, "data/mbx/ps_mbx_rle_nooutliers_adjcounts_fixedmapping.rds")

