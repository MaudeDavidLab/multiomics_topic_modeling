library(cowplot)
library(ggplotify)
library(ape)
library(pheatmap)
library(vegan)
library(ggplot2)
########################################
#setwd("C:/Users/kitikomp/Documents/Lab/multiomics/")
## After normalization choice results ##



#Filter metabolites to include only those in database
#metabolites_not_in_db = read.csv("data/metabol/metabolites_not_in_database.csv")
#metabolites_in_db_directly = taxa_names(ps_metabol)[!(taxa_names(ps_metabol) %in% metabolites_not_in_db$BIOCHEMICAL)]
#keep1 = !(metabolites_not_in_db$IN_MICROBIAL_DATABASE %in% c("FALSE", 'FALSE; not in pubchem', 'FALSE, not in pubchem'))
#keep2 = !is.na(metabolites_not_in_db$IN_MICROBIAL_DATABASE)
#keep = keep1 & keep2
#metabolites_in_db_indirectly = as.character(metabolites_not_in_db$BIOCHEMICAL[keep])
#metabolites_keep = c(metabolites_in_db_directly, metabolites_in_db_indirectly)
#metabolites_keep = metabolites_keep[!grepl("X - ", metabolites_keep)]
#ps_metabol_filt <- prune_taxa(metabolites_keep, ps_metabol)



runMantelTest <- function(dist1, dist2, samples, null_test = F){
  if(null_test){
    dist1 <- dist1[rownames(dist1) %in% samples, colnames(dist1) %in% samples]
    dist2 <- dist2[rownames(dist2) %in% samples, colnames(dist2) %in% samples]
    print(paste("Row names match: ", sum(rownames(dist1) == rownames(dist2))))
    return(mantel(dist1, dist2, method = 'pearson'))
  }else{
    dist1 = dist1[samples, samples]
    dist2 = dist2[samples, samples]
    print(paste("Row names match: ", sum(rownames(dist1) == rownames(dist2))))
    return(mantel(dist1, dist2 ))
  }
  
}
getHeatmaps <- function(ps_16s, ps_mtg, ps_mtt, ps_metabol, method = "bray"){
  
  dist_16s <- as.data.frame(as.matrix(vegdist(t(otu_table(ps_16s)), method = method )))
  dist_mtg <- as.data.frame(as.matrix(vegdist(t(otu_table(ps_mtg)), method = method )))
  dist_mtt <- as.data.frame(as.matrix(vegdist(t(otu_table(ps_mtt)), method = method )))
  dist_metabol <- as.data.frame(as.matrix(vegdist(t(otu_table(ps_metabol)), method = method )))
  #dist_mag <- as.data.frame(as.matrix(vegdist(otu_table(ps_mag), method = method )))
  #dist_metabol_filt <- as.data.frame(as.matrix(vegdist(t(otu_table(ps_metabol_filt)), method = 'bray')))
  
  samples_16s_mtg <- intersect(sample_names(ps_16s), sample_names(ps_mtg))
  samples_16s_mtt <- intersect(sample_names(ps_16s), sample_names(ps_mtt))
  samples_16s_metabol <- intersect(sample_names(ps_16s), sample_names(ps_metabol))
  #samples_16s_mag <- intersect(sample_names(ps_16s), sample_names(ps_mag))
  
  samples_mtg_mtt <- intersect(sample_names(ps_mtg), sample_names(ps_mtt))
  samples_mtg_metabol <- intersect(sample_names(ps_mtg), sample_names(ps_metabol))
 # samples_mtg_mag <- intersect(sample_names(ps_mtg), sample_names(ps_mag))
  
  samples_mtt_metabol <- intersect(sample_names(ps_mtt), sample_names(ps_metabol))
  #samples_mtt_mag <- intersect(sample_names(ps_mtt), sample_names(ps_mag))
  
 # samples_metabol_mag <- intersect(sample_names(ps_metabol), sample_names(ps_mag))
  
  df <- data.frame(omic1 = rep("", 5), omic2 = rep("", 5), statistic = rep(0, 5), signif = rep(0, 5))
  
  mantel_16s_mtg <- runMantelTest(dist_16s, dist_mtg, samples_16s_mtg)
  mantel_16s_mtt <- runMantelTest(dist_16s, dist_mtt, samples_16s_mtt)
  mantel_16s_metabol <- runMantelTest(dist_16s, dist_metabol, samples_16s_metabol)
  #mantel_16s_mag <- runMantelTest(dist_16s, dist_mag, samples_16s_mag)
  
  #mantel_16s_metabol_filt <- runMantelTest(dist_16s, dist_metabol_filt, samples_16s_metabol)
  mantel_mtg_mtt <- runMantelTest(dist_mtg, dist_mtt, samples_mtg_mtt)
  mantel_mtg_metabol <- runMantelTest(dist_mtg, dist_metabol, samples_mtg_metabol)
 # mantel_mtg_mag <- runMantelTest(dist_mtg, dist_mag, samples_mtg_mag)
  #mantel_mtg_metabol_filt <- runMantelTest(dist_mtg, dist_metabol_filt, samples_mtg_metabol)
  mantel_mtt_metabol <- runMantelTest(dist_mtt, dist_metabol, samples_mtt_metabol)
 # mantel_mtt_mag <- runMantelTest(dist_mtt, dist_mag, samples_mtt_mag)
  #mantel_mtt_metabol_filt <- runMantelTest(dist_mtt, dist_metabol_filt, samples_mtt_metabol)
  #mantel_metabol_metabol_filt <- runMantelTest(dist_metabol, dist_metabol_filt, sample_names(ps_metabol))
  #mantel_metabol_mag <- runMantelTest(dist_metabol, dist_mag, samples_metabol_mag)
  

  c("16s", "mtg", "mtt", "metabol", "mag")
  pval_heatmap <- matrix(c(0, mantel_16s_mtg$signif, mantel_16s_mtt$signif, mantel_16s_metabol$signif,
                           mantel_16s_mtg$signif, 0,  mantel_mtg_mtt$signif, mantel_mtg_metabol$signif, 
                           mantel_16s_mtt$signif, mantel_mtg_mtt$signif, 0, mantel_mtt_metabol$signif,  
                           mantel_16s_metabol$signif, mantel_mtg_metabol$signif, mantel_mtt_metabol$signif, 0),
                         byrow = T, ncol =4)
  #mantel_16s_metabol_filt$signif, mantel_mtg_metabol_filt$signif, mantel_mtt_metabol_filt$signif, mantel_metabol_metabol_filt$signif, 0)
  
  z_heatmap <- matrix(c(0, mantel_16s_mtg$statistic,  mantel_16s_mtt$statistic, mantel_16s_metabol$statistic,
                        mantel_16s_mtg$statistic, 0,  mantel_mtg_mtt$statistic, mantel_mtg_metabol$statistic, 
                        mantel_16s_mtt$statistic, mantel_mtg_mtt$statistic, 0, mantel_mtt_metabol$statistic, 
                        mantel_16s_metabol$statistic, mantel_mtg_metabol$statistic, mantel_mtt_metabol$statistic, 0),
                        byrow = T, ncol =4)
  #mantel_16s_metabol_filt$statistic, mantel_mtg_metabol_filt$statistic, mantel_mtt_metabol_filt$statistic, mantel_metabol_metabol_filt$statistic, 0)
  diag(z_heatmap) <- rep(1, 4)
  rownames(z_heatmap) <- c("16s", "mtg",  "mtt", "metabol")
  colnames(z_heatmap) <- c("16s", "mtg",  "mtt", "metabol")
  rownames(pval_heatmap) <- c("16s", "mtg",  "mtt", "metabol")
  colnames(pval_heatmap) <- c("16s", "mtg", "mtt", "metabol")
  
  z_heatmap <- round(z_heatmap, 3)
  z_heatmap_labels <- z_heatmap
  z_heatmap_labels[pval_heatmap < .05 & z_heatmap < 1] <- paste(z_heatmap[pval_heatmap < .05 & z_heatmap < 1], "*")
  
  stat = pheatmap(z_heatmap,
                  display_numbers= z_heatmap_labels,
                  color = diverging_hcl(100, palette = "PurpleGreen", rev = T, alpha = 0.8),
                  fontsize_number=15,
                  cellheight=20,
                  number_color = "black",
                  cluster_cols = F,
                  cluster_rows= F,
                  fontsize_col = 14,
                  fontsize_row = 14)
  #p = pheatmap(pval_heatmap, display_numbers=T, cluster_rows = F, cluster_col = F,
  #             fontsize = 18, fontsize_row = 14, fontsize_col = 14,  color = sequential_hcl(100, palette = "Heat2")) #1-p val. Red is more correlated
  stat = as.ggplot(stat) + ggtitle(paste("Omic similarity"))
  #p = as.ggplot(p) + ggtitle("Mantel test p-value")
  return(stat)
}


#t = ggdraw() + draw_label("best")
#plot_grid(as.ggplot(stat), as.ggplot(p), nrow = 2)
#as expected, mtg tracks well with both 16s and mtt, while 16s tracks less well with mtt (one step further)
#metabolomics doesn't track with anything


#CCA
#tmp <- phyloseq(otu_table(ps_16s@otu_table, taxa_are_rows = T), as.data.frame(t(ps_mtg@otu_table)))
#samples <- samples_16s_mtg
#df1 <- data.frame(t(otu_table(ps_16s@otu_table)))
#df1 <- df1[samples, ]
#df2 <- data.frame(t(otu_table(ps_mtg@otu_table)))
#df2 <- df2[samples, ]
#ccamodel <- cca(df1~., df2)
#plot(ccamodel)



