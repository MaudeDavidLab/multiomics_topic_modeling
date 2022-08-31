ps <- readRDS("~/../Downloads/ps_rle.rds")
ps_5d <- prune_taxa(taxa_names(ps) == "5-dodecenoate (12:1n7)" , ps)
aut_samples <- ps@otu_table[, ps@sam_data$phenotype == "A"]
nt_samples <- ps@otu_table[, ps@sam_data$phenotype == "N"]
wilcox.test(as.numeric(aut_samples), as.numeric(nt_samples))


ps <- readRDS("~/Lab/M3/data/16s/ps_DeSeq_pass_min_postDD_min0.03.rds")
sample_names(ps) <- ps@sam_data$Biospecimen.Name
samples <- intersect(sample_names(ps_5d), sample_names(ps))

ps_16s <- prune_samples(samples, ps)
ps_5d <- prune_samples(samples, ps_5d)

tab_16s <- ps_16s@otu_table
tab_5d <- ps_5d@otu_table

tab_16s <- tab_16s[ , samples]
tab_5d <- tab_5d[, samples]
colnames(tab_16s) == colnames(tab_5d)

pvals <- apply(tab_16s, 1, function(x){
  return(cor.test(x, as.numeric(tab_5d))$p.value)
})

values <- apply(tab_16s, 1, function(x){
  return(cor.test(x, as.numeric(tab_5d))$statistic)
})

asvs <- rownames(tab_16s)[p.adjust(pvals, "fdr") < .05]
values[p.adjust(pvals, "fdr") < .05]
ps_16s@tax_table[asvs, ]
