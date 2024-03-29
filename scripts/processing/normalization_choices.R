library(pbapply)
library(vegan)
library(DESeq2)
library(metagenomeSeq)
library(phyloseq)
library(reshape2)
library(ggplot2)
library(edgeR)
library(compositions)
library(RColorBrewer)

distance_metrics <- c('manhattan', 'euclidean', 'canberra', 'clark', 'bray', 'kulczynski', 'jaccard', 
                      'altGower')

deSeqNorm <- function(ps, scale = F){
  if(scale){
    otu_table(ps) <- ps@otu_table / 100
  }
  ps_dds <- phyloseq_to_deseq2(ps, ~ phenotype)
  ps_dds <- estimateSizeFactors(ps_dds, type = "poscounts")
  ps_dds <- estimateDispersions(ps_dds)
  abund <- getVarianceStabilizedData(ps_dds)
  abund <- abund + abs(min(abund)) #don't allow deseq to return negative counts
  ps_deSeq <- phyloseq(otu_table(abund, taxa_are_rows = T), sample_data(ps), tax_table = tax_table(ps))
  return(ps_deSeq)
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

clrNorm <- function(ps){
  tab <- ps@otu_table
  tab <- tab + 1
  tab <- apply(tab, 2, function(x) {
    geo_mean <- gm_mean(x)
    return(log(x / geo_mean))
  })
  tab <- tab + abs(min(tab))
  ps <- phyloseq(otu_table(tab, taxa_are_rows = T), sample_data(ps), tax_table = tax_table(ps))
  return(ps)
}

tmmNorm <- function(ps){
  if(taxa_are_rows(ps)){
    norm_factors <- edgeR::calcNormFactors(data.frame(ps@otu_table), method = "TMM")
    seqtab_tmm <- sweep(ps@otu_table, 2, norm_factors, '/')
  }else{
    norm_factors <- edgeR::calcNormFactors(t(ps@otu_table), method = "TMM")
    seqtab_tmm <- sweep(t(ps@otu_table), 2, norm_factors, '/')
  }
  ps <- phyloseq(otu_table(seqtab_tmm, taxa_are_rows = T), sample_data = sample_data(ps), tax_table = tax_table(ps))
  return(ps)
}

rleNorm <- function(ps, scale = F){
  if(scale){
    s = 10000
  }else{
    s = 1
  }
  if(taxa_are_rows(ps)){
    norm_factors <- edgeR::calcNormFactors(ps@otu_table, method = "RLE")
    seqtab_mtg_tmm <- sweep(ps@otu_table, 2, norm_factors^2, '/')
  }else{
    norm_factors <- edgeR::calcNormFactors(t(ps@otu_table), method = "RLE")
    seqtab_mtg_tmm <- sweep(t(ps@otu_table), 2, norm_factors^2, '/')
  }
  ps <- phyloseq(otu_table(seqtab_mtg_tmm, taxa_are_rows = T), sample_data = sample_data(ps), tax_table = tax_table(ps))
  return(ps)
}

cssNorm <- function(ps){
  p_metag<-phyloseq_to_metagenomeSeq(ps)
  norm_df <- MRcounts(p_metag, norm=T)
  otu_table(ps) <- otu_table(norm_df, taxa_are_rows = T)
  return(ps)
}

#Calculate the average distance between replicates. Compare that to the average distance between random individuals.
#The lower this ratio, the better the normalization method is. This will also depend on the distance metric selected.
getInOutRatio <- function(ps, dists, dist_metric){
  mean_replicate_dists <- pbsapply(unique(ps@sam_data$host_name), function(host_name){
    samples <- sample_names(ps)[ps@sam_data$host_name == host_name]
    mean_dist <- mean(unlist(dists[samples, samples]))
    return(mean_dist)
  }) # one element per subject
  mean_nonreplicate_dists <- pbsapply(unique(ps@sam_data$host_name), function(host_name){
    return(mean(unlist(dists[ps@sam_data$host_name == host_name, ps@sam_data$host_name != host_name])))
  }) # one element per subject
  
  return(mean_replicate_dists / mean_nonreplicate_dists) #we want this ratio to be low
}

getInOutRatio_sibs <- function(ps, dists, dist_metric){
  mean_replicate_dists <- pbsapply(unique(ps@sam_data$familyID), function(famID){
    samples <- sample_names(ps)[ps@sam_data$familyID == famID]
    mean_dist <- mean(unlist(dists[samples, samples]))
    return(mean_dist)
  }) # one element per subject
  mean_nonreplicate_dists <- pbsapply(unique(ps@sam_data$familyID), function(famID){
    return(mean(unlist(dists[ps@sam_data$familyID == famID, ps@sam_data$familyID != famID])))
  }) # one element per subject
  
  return(mean_replicate_dists / mean_nonreplicate_dists) #we want this ratio to be low
}

getRatioList <- function(ps, comparison){
  meanInOutRatios <- list()
  otu = t(otu_table(ps))
  for(dist_metric in distance_metrics){
    print(dist_metric)
    if(dist_metric == "morisita"){
      otu = round(otu)
    }
    dists <- as.data.frame(as.matrix(vegdist(otu, method = dist_metric)))
    if(comparison == "self"){
      meanInOutRatios[[dist_metric]] <- mean(getInOutRatio(ps, dists, dist_metric))
    }
    if(comparison == "sib"){
      meanInOutRatios[[dist_metric]] <- mean(getInOutRatio_sibs(ps, dists, dist_metric))
    }
  }
  return(meanInOutRatios)
}


makeBoxplot <- function(ratios_dds, ratios_css, ratios_rle, ratios_tmm){
  tmp <- matrix(c( unlist(ratios_dds), unlist(ratios_css), unlist(ratios_rle), unlist(ratios_tmm)), nrow = 4, byrow = T)
  rownames(tmp) <- c("dds", "css", "rle", "tmm")
  colnames(tmp) <- distance_metrics
  tmp <- melt(tmp)
  colnames(tmp) <- c("normalization", "distance_metric", "ratio")
  
  library(reshape2)
  ggplot(tmp, aes(x = normalization, y = ratio, fill = distance_metric, color = distance_metric)) +
    geom_boxplot() + ylim(0, .6) + geom_boxplot(size = 1.7) + 
    scale_color_brewer(palette = "Paired") + scale_fill_brewer(palette = "Paired")+
    theme_bw()+
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14))+
    ylab("Intra/Inter Sibling Ratio")+
    xlab("Normalization method")
  #from the zoomed out view, gower is definitely out as a distance metric. It is too dependent on normalization method
  
  #ggplot(tmp, aes(x = normalization, y = ratio, fill = distance_metric, color = distance_metric)) +
  #  geom_boxplot(size = 1.5) + ylim(0, 10)
  # we should not use gower, canberra, OR clark with rle
  
  #ggplot(tmp, aes(x = normalization, y = ratio, fill = distance_metric, color = distance_metric)) +
  #  geom_boxplot(size = 1.5) + ylim(0, 1)
  
}

runNormalizations <- function(ps_raw, omic, scale = F){
  # Normalize in 3 ways
  
  ps_rle <- rleNorm(ps_raw, scale)
  saveRDS(ps_rle, paste0("../data/", omic, "/ps_rle.rds"))
  ps_css <- cssNorm(ps_raw)
  saveRDS(ps_css, paste0("../data/", omic, "/ps_css.rds"))
  ps_dds <- deSeqNorm(ps_raw, scale)
  saveRDS(ps_dds, paste0("../data/", omic, "/ps_dds.rds"))
  ps_tmm <- tmmNorm(ps_raw)
  saveRDS(ps_tmm, paste0("../data/", omic, "/ps_tmm.rds"))
  #ps_clr <- clrNorm(ps_raw)
  #saveRDS(ps_clr, paste0("../data/", omic, "/ps_clr.rds"))
}

runRatioTest <- function(omic, comparison){
  #Get the average distance between replicates / average distance between non-replicates
  #ps_raw <- readRDS(paste0("../data/", omic, "/ps_clean.rds"))
  #famIDs_keep <- names(table(ps_raw@sam_data$familyID)>1)[table(ps_raw@sam_data$familyID)>1]
  #ps_raw <- prune_samples(ps_raw@sam_data$familyID %in% famIDs_keep, ps_raw)
  #meanInOutRatios_raw <- getRatioList(ps_raw, comparison = comparison)
  #saveRDS(meanInOutRatios_raw, paste0("output/normalizations/", omic, "/mean_replicate_nonreplicate_distances_raw.rds"))
  
  ps_rle <- readRDS(paste0("../data/", omic, "/ps_rle.rds"))
  famIDs_keep <- names(table(ps_rle@sam_data$familyID)>1)[table(ps_rle@sam_data$familyID)>1]
  ps_rle <- prune_samples(ps_rle@sam_data$familyID %in% famIDs_keep, ps_rle)
  meanInOutRatios_rle <- getRatioList(ps_rle, comparison = comparison)
  saveRDS(meanInOutRatios_rle, paste0("../results/summary_statistics/normalizations/", omic, "/mean_replicate_nonreplicate_distances_rle.rds"))
  
  ps_css <- readRDS(paste0("../data/", omic, "/ps_css.rds"))
  famIDs_keep <- names(table(ps_css@sam_data$familyID)>1)[table(ps_css@sam_data$familyID)>1]
  ps_css <- prune_samples(ps_css@sam_data$familyID %in% famIDs_keep, ps_css)
  meanInOutRatios_css <- getRatioList(ps_css, comparison = comparison)
  saveRDS(meanInOutRatios_css, paste0("../results/summary_statistics/normalizations/", omic, "/mean_replicate_nonreplicate_distances_css.rds"))
  
  ps_dds <- readRDS(paste0("../data/", omic, "/ps_dds.rds"))
  famIDs_keep <- names(table(ps_dds@sam_data$familyID)>1)[table(ps_dds@sam_data$familyID)>1]
  ps_dds <- prune_samples(ps_dds@sam_data$familyID %in% famIDs_keep, ps_dds)
  meanInOutRatios_dds <- getRatioList(ps_dds, comparison = comparison)
  saveRDS(meanInOutRatios_dds, paste0("../results/summary_statistics/normalizations/", omic, "/mean_replicate_nonreplicate_distances_dds.rds"))
  
  ps_tmm <- readRDS(paste0("../data/", omic, "/ps_tmm.rds"))
  famIDs_keep <- names(table(ps_tmm@sam_data$familyID)>1)[table(ps_tmm@sam_data$familyID)>1]
  ps_tmm <- prune_samples(ps_tmm@sam_data$familyID %in% famIDs_keep, ps_tmm)
  meanInOutRatios_tmm <- getRatioList(ps_tmm, comparison = comparison)
  saveRDS(meanInOutRatios_tmm, paste0("../results/summary_statistics/normalizations/", omic, "/mean_replicate_nonreplicate_distances_tmm.rds"))
  
  #ps_clr <- readRDS(paste0("../data/", omic, "/ps_clr.rds"))
  #famIDs_keep <- names(table(ps_clr@sam_data$familyID)>1)[table(ps_clr@sam_data$familyID)>1]
  #ps_clr <- prune_samples(ps_clr@sam_data$familyID %in% famIDs_keep, ps_clr)
  #meanInOutRatios_clr <- getRatioList(ps_clr, comparison = comparison)
  #saveRDS(meanInOutRatios_clr, paste0("../results/summary_statistics/normalizations/", omic, "/mean_replicate_nonreplicate_distances_clr.rds"))
}

displayPlot <- function(omic){
  # Evaluate normalizations for criteria
  #1. replicates should be close to each other (16s)
  #ratios_raw <- readRDS(paste0("../results/summary_statistics/normalizations/", omic, "/mean_replicate_nonreplicate_distances_raw.rds"))
  distance_metrics <- c('manhattan', 'euclidean', 'canberra', 'clark', 'bray', 'kulczynski', 'jaccard', 
                        'altGower')
  
  
  ratios_dds <- readRDS(paste0("../results/summary_statistics/normalizations/", omic, "/mean_replicate_nonreplicate_distances_dds.rds"))
  ratios_css <- readRDS(paste0("../results/summary_statistics/normalizations/", omic, "/mean_replicate_nonreplicate_distances_css.rds"))
  ratios_rle <- readRDS(paste0("../results/summary_statistics/normalizations/", omic, "/mean_replicate_nonreplicate_distances_rle.rds"))
  ratios_tmm <- readRDS(paste0("../results/summary_statistics/normalizations/", omic, "/mean_replicate_nonreplicate_distances_tmm.rds"))
  #ratios_clr <- readRDS(paste0("../results/summary_statistics/normalizations/", omic, "/mean_replicate_nonreplicate_distances_clr.rds"))
  
  
  p <- makeBoxplot(ratios_dds[distance_metrics], ratios_css[distance_metrics], ratios_rle[distance_metrics], ratios_tmm[distance_metrics])
  if(omic == "metabol"){
    omic = "mbx"
  }
  p <- p + ggtitle(omic)
  return(p)
}
setwd("~/Lab/M3/scripts")

omic = "16s"
ps_16s <- readRDS("../data/16s/ps_not_norm_mapping_cleaned.rds") #resaved 16s/ps_not_norm_mapping_cleaned.rds as ps_clean.rds
runNormalizations(ps_16s, "16s")
runRatioTest("16s", comparison = "self")
p <- displayPlot("16s")
p
ggsave(paste0("../results/summary_statistics/normalizations/", omic, "normtest.pdf"))
# we should not use gower, canberra, OR clark with rle
#overall dds minimizes the replicate/nonreplicate ratios, which was our criteria for selection. 
#It emerges the winner for 16s. We already think that morisita or horn are going to be the best distance metrics as well

omic = "mtg"
ps_mtg <- readRDS(paste0("../data/", omic, "/ps_clean.rds"))
runNormalizations(ps_mtg, omic)
runRatioTest(omic, comparison = "sib")
p <- displayPlot(omic)
p
ggsave(paste0("../results/summary_statistics/normalizations/", omic, "normtest.pdf"))

omic = "mtt"
ps_mtt <- readRDS(paste0("../data/", omic, "/ps_clean.rds"))
runNormalizations(ps_mtt, omic)
runRatioTest(omic, comparison = "sib")
p <- displayPlot(omic)
p
ggsave(paste0("../results/summary_statistics/normalizations/", omic, "normtest.pdf"))

omic = "metabol"
ps_metabol <- readRDS(paste0("../data/", omic, "/ps_clean.rds"))

# Going back to original scale data
mbx <- read.csv("~/Lab/M3/data/mbx/metabolomics_data_orgi_scale.csv", row.names = 1)
mbx[is.na(mbx)] <- 0
# need to make the numbers smaller or they overflow the integer value
# Scale each metabolite to be median 1
#mbx <- t(apply(mbx, 1, function(x) return(x / (sd(x)))))
mbx <- mbx[, ps_metabol@sam_data$biospecimen_id]
#mbx <- apply(mbx, 2, function(x) return(sapply(x, as.integer)))
colnames(mbx) <- ps_metabol@sam_data$biospecimen_name
colnames(mbx) <- gsub("F", "P", colnames(mbx))
ps_mbx <- phyloseq(otu_table(mbx, taxa_are_rows = T), sample_data(ps_metabol@sam_data), tax_table(ps_metabol@tax_table))
hist(apply(ps_mbx@otu_table, 1, max), main = "max count present for that metabolite")
hist(apply(ps_mbx@otu_table, 2, max), main = "max count present in each sample", breaks = 20)
hist(sample_sums(ps_mbx))

# mean scale each metabolite on top of rle
counts <- ps_mbx@otu_table

ps_mbx_real <- readRDS("../data/mbx/ps_mbx_rle_nooutliers_adjcounts_fixedmapping.rds")
ps_mbx_real@otu_table[1:10, 1:10]
counts[1:10, 1:10]

counts_adj <- apply(counts, 2, function(x) return(x / median(x)))
counts_adj[1:10, 1:10]

runNormalizations(ps_mbx, omic, scale = F)
runRatioTest(omic, comparison = "sib")
p <- displayPlot(omic)
p
ggsave(paste0("../results/summary_statistics/normalizations/", omic, "normtest.pdf"))

p <- ggarrange(displayPlot("16s"), 
               displayPlot("mtg"), 
               displayPlot("mtt"), 
               displayPlot("metabol"), 
               common.legend = T)
p
ggsave("../results/summary_statistics/normalizations/norm_justification.pdf")


omic = "mag"
ps_mag <- readRDS("../data/mtg/ps_mag_fixed_mapping.rds")
runNormalizations(ps_mag, omic)
runRatioTest(omic, comparison = "sib")
p <- displayPlot(omic)
p
ggsave(paste0("../results/summary_statistics/normalizations/", omic, "normtest.pdf"))

p <- ggarrange(displayPlot("16s"), 
               displayPlot("mtg"), 
               displayPlot("mtt"), 
               displayPlot("metabol"), 
               common.legend = T)
p
ggsave("../results/summary_statistics/normalizations/normalization_justification.pdf", width = 10, height = 6)

