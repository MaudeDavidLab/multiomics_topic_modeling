
library(KEGGREST)
library("EnrichmentBrowser")
library(clusterProfiler)
dir = "C:/Users/kitikomp/Documents/Lab/M3/"
setwd(dir)




############################
#log2 ratio of means########
############################
### MTG ####################
ps <- readRDS("data/mtg/ps_rle_nooutliers.rds")
tmp <- apply(ps@otu_table, 1, function(x) {
  aut <- ps@sam_data$phenotype == "A"
  nt <- ps@sam_data$phenotype == "N"
  return(log2(mean(x[aut]) / mean(x[nt])))
  
})
df <- data.frame(tmp)
write.table(df, file = "results/gsea/mtg_log2mean_ranklist.csv", row.names = T, sep= ",")

### MTT ####################
ps <- readRDS("data/mtt/ps_rle_nooutliers.rds")
tmp <- apply(ps@otu_table, 1, function(x) {
  aut <- ps@sam_data$phenotype == "A"
  nt <- ps@sam_data$phenotype == "N"
  return(log2(mean(x[aut]) / mean(x[nt])))
  
})
df <- data.frame(tmp)
write.table(df, file = "results/gsea/mtt_log2mean_ranklist.csv", row.names = T, sep= ",", quote = F)


##################################
#### Get gene list  ##############
##################################
org <- keggList("organism")
prokaryotes <- grepl("Prokaryotes", org[,4])
org <- org[prokaryotes, "organism"]

pathways <- lapply(org, function(code){
  tryCatch(
    {
      tmp <- keggLink("pathway", code)
      pathway_numbers <- gsub("path:", "", as.character(tmp))
      pathway_numbers <- gsub(code, "", as.character(pathway_numbers))
      pathway_numbers <- unique(pathway_numbers)
      return(pathway_numbers)
    }, error=function(cond){
      return(NULL)
    })
  
})

pathways_final <- unique(unlist(pathways))
pathways_final <- paste("ko", pathways_final, sep = "")
saveRDS(pathways_final, "results/gsea/prokaryote_pathway_ids")
pathway_ids <- readRDS("results/gsea/prokaryote_pathway_ids")

pathways_final <- readRDS("results/gsea/prokaryote_pathway_ids")

library(KEGGREST)

for(pathway in pathways_final){
  tmp <- keggGet(pathway)
  genes <- names(tmp[[1]]$ORTHOLOGY)
  file = 'results/gsea/prokaryote_gene_set.gmt'
  line <- paste(pathway, paste(genes, collapse = ","), sep = ",")
  write(line, 'results/gsea/prokaryote_gene_set.gmt', append = T)
}









getActualDiffs <- function(ps){
  seqtab <- t(as.data.frame(ps@otu_table)) #sample by ASV
  #For each sibling pair, find actual difference
  counts <- table(ps@sam_data$familyID)
  remove <- names(counts)[counts < 2]
  remove <- c(remove, "45") #both autism
  fam_ids <- unique(ps@sam_data$familyID)
  fam_ids <- fam_ids[!(fam_ids %in% remove)]
  
  actual_diffs <- lapply(fam_ids, function(fam_id){
    sib_aut <- sample_names(ps)[ps@sam_data$familyID == fam_id & ps@sam_data$phenotype == 'A']
    sib_nt <- sample_names(ps)[ps@sam_data$familyID == fam_id & ps@sam_data$phenotype == 'N']
    print(fam_id)
    return(seqtab[sib_aut, ] - seqtab[sib_nt, ])
  })
  actual_diffs <- do.call(rbind, actual_diffs)
}

findSig_Binomial_removezeros <- function(actual_diffs, pos, type = "ko", thresh = 0.05){
  pvals <- c()
  num_samples <- c()
  ko_ids <- c()
  for(ko in colnames(actual_diffs)){
    diffs <- actual_diffs[, ko]
    num_sample <- sum(diffs != 0)
    if(num_sample > 5){
      print(paste("Num_samples: ", num_sample))
      if(pos){
        observed_successes <- sum(diffs > 0)
      }else{
        observed_successes <- sum(diffs < 0)
      }
      pval <- 1 - sum(dbinom(seq(0, observed_successes-1), p = 0.5, size = num_sample))  
      
      pvals <- c(pvals, pval)
      num_samples <- c(num_samples, num_sample)
      ko_ids <- c(ko_ids, ko)
    }
  }
  pvals_adj <- p.adjust(pvals, method = "fdr")
  #print(pvals)
  #hist(pvals)
  #ko_ids <- colnames(actual_diffs)
  df <- data.frame(ko_ids, pvals_adj, pvals, num_samples)
  #ko_ids <- colnames(actual_diffs)[which(pvals < thresh)]
  #df <- data.frame(ko_ids, pvals[pvals< thresh])
  return(df)
}


############### Binomial model  ############
actual_diffs <- getActualDiffs(ps)
sig_kos_pos <- findSig_Binomial_removezeros(actual_diffs, pos = T)
sig_kos_neg <- findSig_Binomial_removezeros(actual_diffs, pos = F)

sig_kos_pos$importance <- 1 - sig_kos_pos$pvals
ord <- sig_kos_pos[order(sig_kos_pos$importance, decreasing = T), ]
ord_skinny <- ord %>% select(-c("pvals_adj", "pvals", "num_samples"))



library(KEGGREST)
gene_symbols <- sapply(ord_skinny$ko_ids, function(x) {
  try({
    tmp <- keggGet(x)
    return(tmp[[1]]$NAME)
  }) 
  return(NA)
})

ord_skinny$gene_symbol = gene_symbols[ord_skinny$ko_ids]
write.table(ord_skinny, file = "results/gsea/mtg_binomial_ranklist.csv", row.names = F, sep= ",")


gene_set <- list()
pathways <- keggList("pathway")
print(length(pathways)) #548
i=0
for(path in names(pathways)){
  print(i, end = "\t")
  try({
    tmp <- keggGet(path)
    ko_pathway <- tmp[[1]]$KO_PATHWAY
    tmp <- keggGet(ko_pathway)
    genes <- names(tmp[[1]]$ORTHOLOGY)
    gene_set[[ko_pathway]] = genes
  })
  i = i + 1

}

df <- rbind.fill(lapply(gene_set,function(y){as.data.frame(t(y),stringsAsFactors=FALSE)}))

ddd <- data.frame(I(unlist(lapply(gene_set,paste,collapse=","))))

rownames(df) <- names(gene_set)
write.table(ddd, "kegg_gene_set.gmt", col.names = F, quote = F, na = "", sep = ",")

