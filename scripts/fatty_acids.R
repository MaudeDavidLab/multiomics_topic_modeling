# This script is meant to answer the question, "do we see more x because the gut bacteria are producing more x?"
library(KEGGREST)
library(ggplot2)

ps_mtt <- readRDS("C:/Users/ctata/Documents/Lab/M3/PS_kegg_fun_MTT.rds")
print(nsamples(ps_mtt))
ps_mtt <- addNumericMetadata(ps_mtt)
ps_mtt <- prune_samples(sample_sums(ps_mtt) > 100000, ps_mtt)
print(nsamples(ps_mtt))

set.seed(0)
ps_mtt <- rarefy_even_depth(ps_mtt, sample.size = min(sample_sums(ps_mtt)))
print(nsamples(ps_mtt))
saveRDS(ps_mtt, "ps_metatranscriptome.rds")

ord <- ordinate(ps_mtt, "PCoA", "bray")
ps_mtt@sam_data$depth <- sample_sums(ps_mtt)
plot_ordination(ps_mtt, ord, type = "sample", color = "")


ps_mtg <- readRDS("ps_metagenome.rds")
###########################################################################################################################

gene_list <- c("K00665", "K00667", "K11533", "K00647", "K09458", "K11539", "K13370", "K00059", "K01716", "K02372", "K22540",
               "K00208", "K02371", "K10780", "K00209", "K07512")


getEnzymeKos <- function(enzyme){
  kos <- NA
  tryCatch({
    kos <- keggGet(enzyme)[[1]]$NAME
  }, error = function(e){
    print(paste("Enzyme error"))
  })
  return(kos)
}

ko_names <- lapply(gene_list, getEnzymeKos)
ko_names <- unlist(lapply(ko_names, function(x) return(strsplit(x, split = ", ")[[1]][1])))
names(ko_names) <- gene_list

#########################################################################################################

plotBoxplots <- function(ps_selected){
  data = data.frame(t(ps_selected@otu_table))
  data$phenotype <- ps_selected@sam_data$phenotype
  data$sampleID <- rownames(data)
  data = melt(data)
  ggplot(data = data, aes(x = phenotype, y = value)) + geom_boxplot(aes(fill = phenotype)) + facet_wrap(~variable)

}

doWilcoxTest <- function(ps_aut, ps_nt){
  for(i in seq(1,12)){
    print(wilcox.test(as.numeric(ps_aut@otu_table[i,]), as.numeric(ps_nt@otu_table[i,])))
  }
}

pairSamples <- function(ps_aut, ps_nt){
  familyIds <- intersect(ps_aut@sam_data$familyID, ps_nt@sam_data$familyID)
  ps_aut <- prune_samples(ps_aut@sam_data$familyID %in% familyIds, ps_aut)
  ps_aut <- prune_samples(!duplicated(ps_aut@sam_data$familyID), ps_aut)
  ps_nt <- prune_samples(ps_nt@sam_data$familyID %in% familyIds, ps_nt)
  ps_nt <- prune_samples(!duplicated(ps_nt@sam_data$familyID), ps_nt)
  ps_aut@sam_data$familyID == ps_nt@sam_data$familyID
  return(list(ps_aut, ps_nt))
}

plotDiffs <- function(ps_aut, ps_nt){
  plots <- list()
  diffs_list <- list()
  par(mfrow = c(3, 4))
  for(i in seq(1,12)){
    diffs <- ps_aut@otu_table[i, ] - ps_nt@otu_table[i, ]
    #p <- hist(diffs, breaks = 40, main = ko_names[i])
    #plots[[i]] <- p
    diffs_list[[i]] <- diffs
   # p
  }
  return(diffs_list)
}


binTest <- function(diffs, aut){
  if(aut){
    observed_successes = sum(diffs > 0) #if the difference is high, enriched in aut
  }else{
    observed_successes = sum(diffs < 0) #if the difference in low, enriched in nt
  }
  num_sample = length(diffs)
  pval <- 1 - sum(dbinom(seq(0, observed_successes-1), p = 0.5, size = num_sample))
  print(pval)
  return(pval)
}



runTests <- function(ps, aut){
  ps_selected <- prune_taxa(gene_list, ps)
  ps_aut <- prune_samples(ps_selected@sam_data$phenotype == "A", ps_selected)
  ps_nt <- prune_samples(ps_selected@sam_data$phenotype == "N", ps_selected)
  tmp <- pairSamples(ps_aut, ps_nt)
  ps_aut <- tmp[[1]]
  ps_nt <- tmp[[2]]
  diffs_list <- plotDiffs(ps_aut, ps_nt)
  pvals <- unlist(lapply(diffs_list, binTest, aut))
  return(pvals)
}

pvals_mtt <- runTests(ps_mtt, aut=T)
pvals_mtg <- runTests(ps_mtg, aut=T)
df_aut <- data.frame(taxa_names(ps_selected), as.character(ko_names[taxa_names(ps_selected)]), pvals_mtg, pvals_mtt)

pvals_mtt <- runTests(ps_mtt, aut = F)
pvals_mtg <- runTests(ps_mtg, aut= F)
df_nt <- data.frame(taxa_names(ps_selected), as.character(ko_names[taxa_names(ps_selected)]), pvals_mtg, pvals_mtt)



#are there more fatty acid genes in general in one or the other
ps_selected <- prune_taxa(gene_list, ps_mtg)
ps_aut <- prune_samples(ps_selected@sam_data$phenotype == "A", ps_selected)
ps_nt <- prune_samples(ps_selected@sam_data$phenotype == "N", ps_selected)
tmp <- pairSamples(ps_aut, ps_nt)
ps_aut <- tmp[[1]]
ps_nt <- tmp[[2]]

aut_fattyacid <- colSums(ps_aut@otu_table)
nt_fattyacid <- colSums(ps_nt@otu_table)
par(mfrow = c(1, 1))
hist(aut_fattyacid - nt_fattyacid, breaks = 40)
binTest(aut_fattyacid - nt_fattyacid, aut = F)
wilcox.test(aut_fattyacid, nt_fattyacid)
#There's not more fatty acid biosynthesis genes at large in one grou por another


############################
############################
### Fatty Acid Metabolites
ps_metabolite <- readRDS("ps_metabolite.rds")
sat_fatty_acids <- taxa_names(ps_metabolite)[grep("ate.*\\(.*0)", taxa_names(ps_metabolite))]
unsat_fatty_acids <- taxa_names(ps_metabolite)[grep("ate.*\\(.*:.*[^0])", taxa_names(ps_metabolite))]
carnitines <- taxa_names(ps_metabolite)[grep("carnitine", taxa_names(ps_metabolite) )]
carnitines <- c("carnitine")

ps_selected <- prune_taxa(carnitines, ps_metabolite)
ps_aut <- prune_samples(ps_selected@sam_data$phenotype == "A", ps_selected)
ps_nt <- prune_samples(ps_selected@sam_data$phenotype == "N", ps_selected)
tmp <- pairSamples(ps_aut, ps_nt)
ps_aut <- tmp[[1]]
ps_nt <- tmp[[2]]
binTest(sample_sums(ps_aut) - sample_sums(ps_nt), aut = F)
hist(sample_sums(ps_aut) - sample_sums(ps_nt), breaks = 40)

aut <- ps_selected@sam_data$phenotype == "A"
nt <- ps_selected@sam_data$phenotype == "N"
aut_fattyacid <- sample_sums(ps_selected)[aut]
nt_fattyacid <- sample_sums(ps_selected)[nt]
wilcox.test(aut_fattyacid, nt_fattyacid)
mean(aut_fattyacid)
mean(nt_fattyacid)

par(mfrow = c(1, 2))
hist(unsat[aut] / sat[aut], breaks = 60, xlim = c(0,12))
hist(unsat[nt] / sat[nt], breaks = 30)
wilcox.test(unsat[aut] / sat[aut], unsat[nt] / sat[nt])
#There are more fatty acids IN GENERAL in neurotypicals
# more saturated fatty acids in NT
# NO difference in unsaturated fatty acids
# Options
#1. NT are absorbing less saturated fatty acids
#2. NT are producing more saturated fatty acids


#I have no idea how to tell if the metagenome and metatranscriptome support that microbes are actually producing more
