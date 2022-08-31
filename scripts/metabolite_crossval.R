library(KEGGREST)

metabolite_nominations <- read.delim("C:/Users/ctata/Documents/Lab/M3/metabolite_res_combined.tsv", row.names=1)
kegg_compound_ids <- metabolite_nominations$KEGG
kegg_compound_ids <- as.character(kegg_compound_ids)
kegg_compound_ids <- kegg_compound_ids[kegg_compound_ids != ""]
kegg_compound_ids <- kegg_compound_ids[!is.na(kegg_compound_ids)]

getEnzymeKos <- function(enzyme){
  kos <- NA
  tryCatch({
    kos <- names(keggGet(enzyme)[[1]]$ORTHOLOGY)
  }, error = function(e){
    print(paste("Enzyme error"))
  })
  return(kos)
}




#Lots of the compounds don't have modules associated with them
kos_associated_metabolite <- list()
for(id in kegg_compound_ids){
  print(id)
  tryCatch({
    paths <- gsub("map", "ko", names(keggGet(id)[[1]]$PATHWAY))
    kos_associated_metabolite[[id]] <- paths
    #enzymes <- keggGet(id)[[1]]$ENZYME
    #kos_list <- lapply(enzymes, getEnzymeKos)
    #kos <- unlist(kos_list)
    #kos <- kos[!is.na(kos)]
    
    #modules <- keggGet(id)[[1]]$MODULE
    #orth <- lapply(names(modules), function(mod) return(keggGet(mod)[[1]]$ORTHOLOGY))
    #print(length(orth))
    #kos_associated_metabolite[[id]] <- unlist(strsplit(names(unlist(orth)), ","))
    #kos_associated_metabolite[[id]] <- kos
  }, error = function(e){
    print("Compound Error")
  })
}



ko_res <- read.delim("C:/Users/ctata/Documents/Lab/M3/ko_res.tsv", row.names = 1)
sig_KOs <- rownames(ko_res)
kos_associated_KO <- list()
for(id in sig_KOs){
  print(id)
  tryCatch({
    paths <- names(keggGet(id)[[1]]$PATHWAY)
    kos_associated_KO[[id]] <- paths
  }, error = function(e){
    print("Compound Error")
  })
}



enzymes_associated_KO <- list()
for(id in sig_KOs){
  print(id)
  tryCatch({
    tmp <- keggGet(id)[[1]]
    enzymes <- strsplit(gsub("]", "", strsplit(tmp$DEFINITION, "EC:")[[1]][2]), " ")[[1]]
    enzymes_associated_KO[[id]] <- enzymes
  }, error = function(e){
    print("Error")
  })
}
enzymes_ko <- unlist(enzymes_associated_KO)
enzymes_ko <- enzymes_ko[!is.na(enzymes_ko)]

mtt_res <- read.delim("C:/Users/ctata/Documents/Lab/M3/mtt_res.tsv", row.names = 1)
sig_mtt <- rownames(mtt_res)
kos_associated_mtt <- list()
for(id in sig_mtt){
  print(id)
  tryCatch({
    paths <- names(keggGet(id)[[1]]$PATHWAY)
    kos_associated_mtt[[id]] <- paths
  }, error = function(e){
    print("Compound Error")
  })
}


enzymes_associated_mtt <- list()
for(id in sig_mtt){
  print(id)
  tryCatch({
    tmp <- keggGet(id)[[1]]
    enzymes <- strsplit(gsub("]", "", strsplit(tmp$DEFINITION, "EC:")[[1]][2]), " ")[[1]]
    enzymes_associated_mtt[[id]] <- enzymes
  }, error = function(e){
    print("Error")
  })
}
enzymes_mtt <- unlist(enzymes_associated_mtt)
enzymes_mtt <- enzymes_mtt[!is.na(enzymes_mtt)]



res <- read.delim("res_combined_mtg_mtt.tsv")
enzymes_list <- list()
for(id in res$ko_names){
  print(id)
  enzymes_list[[id]] <- NA
  tryCatch({
    tmp <- keggGet(id)[[1]]
    enzymes <- strsplit(gsub("]", "", strsplit(tmp$DEFINITION, "EC:")[[1]][2]), " ")[[1]]
    enzymes_list[[id]] <- enzymes
  }, error = function(e){
    print("Error")
  })
}
enzyme_list_collapse <- lapply(enzymes_list, function(x) return(paste(x, collapse = ", ")))
enzyme_list_collapse <- unlist(enzyme_list_collapse)
res$enzymes <- enzyme_list_collapse[res$ko_names]


write.table(res, "res_combined_mtg_mtt_enzymes.tsv", sep = "\t", quote = F, row.names = F)



pathway_res <- read.delim("C:/Users/ctata/Documents/Lab/M3/pathway_res.tsv", row.names = 1)
pathway_names <- rownames(pathway_res)
kos_from_pathways <- lapply(pathway_names, function(path_id){
  kos <- keggGet(path_id)[[1]]$ORTHOLOGY
  return(names(kos))
})

sig_kos <- c(sig_kos, unlist(kos_from_pathways))

################################################
#######    Pathways   ##########################
################################################
mtt_mtg_res <- read.delim("C:/Users/ctata/Documents/Lab/M3/res_combined_mtg_mtt.tsv")
paths_for_kos  <- pblapply(mtt_mtg_res$ko_names, function(id){
  tryCatch({
    paths <- names(keggGet(id)[[1]]$PATHWAY)
    return(paths)
  }, error = function(e){
    print("Error")
  })
})
saveRDS(paths_for_kos, "paths_for_kos_combined_mtg_mtt.rds")
paths_for_mtt_mtg <- unlist(readRDS("paths_for_kos_combined_mtg_mtt.rds"))
paths_for_mtt_mtg <- paths_for_mtt_mtg[paths_for_mtt_mtg != "Error"]
paths_for_mtt_mtg <- table(paths_for_mtt_mtg)
paths_for_mtt_mtg <- paths_for_mtt_mtg[order( paths_for_mtt_mtg)]

piph_res <- read.delim("C:/Users/ctata/Documents/Lab/M3/piph_res_combined.tsv")
paths_for_kos  <- pblapply(as.character(piph_res$X), function(id){
  tryCatch({
    paths <- names(keggGet(id)[[1]]$PATHWAY)
    return(paths)
  }, error = function(e){
    print("Error")
  })
})
saveRDS(paths_for_kos, "paths_for_kos_combined_piph.rds")

paths_for_piph <- unlist(readRDS("paths_for_kos_combined_piph.rds"))
paths_for_piph <- paths_for_piph[paths_for_piph != "Error"]
paths_for_piph <- table(paths_for_piph)
paths_for_piph <- paths_for_piph[order(paths_for_piph)]


paths_for_metabolites  <- pblapply(as.character(kegg_compound_ids), function(id){
  tryCatch({
    paths <- names(keggGet(id)[[1]]$PATHWAY)
    return(paths)
  }, error = function(e){
    print("Error")
  })
})
saveRDS(paths_for_metabolites, "paths_for_metabolites_combined.rds")

paths_for_metabolites <- unlist(readRDS("paths_for_metabolites_combined.rds"))
paths_for_metabolites <- paths_for_metabolites[paths_for_metabolites != "Error"]
paths_for_metabolites <- gsub("map", "ko", paths_for_metabolites)
paths_for_metabolites <- table(paths_for_metabolites)
paths_for_metabolites <- paths_for_metabolites[order(paths_for_metabolites)]


###############
## Intersect
paths_keep <- intersect(names(paths_for_mtt_mtg), names(paths_for_piph))
paths_keep <- intersect(paths_keep, names(paths_for_metabolites))
paths_keep <- unique(paths_keep)



############################################
# Each compound is a part of some pathways
# Each KO is a part of some pathways
# Find significantly different compounds and KOs independently, then look for pathway overlap
# Score = # of times a metabolite's pathway overlaps with any KO's pathway / # of pathway for that metabolite
matches_metagenome <- lapply(kos_associated_metabolite, function(kos){
  num_sig_pathways <- sapply(kos, function(ko){
    return(as.numeric(table(unlist(kos_associated_KO))[ko])) #number of times pathway was hit by significant kos
  })
  return(num_sig_pathways)
})

matches_mtt<- lapply(kos_associated_metabolite, function(kos){
  num_sig_pathways <- sapply(kos, function(ko){
    return(as.numeric(table(unlist(kos_associated_mtt))[ko])) #number of times pathway was hit by significant kos
  })
  return(num_sig_pathways)
})


keep <- sapply(matches_mtt, function(x) return(length(x) > 0))
matches_mtt <- matches_mtt[keep]

keep <- sapply(matches_metagenome, function(x) return(length(x) > 0))
matches_metagenome <- matches_metagenome[keep]

paths_per_compound <- unlist(lapply(kos_associated_metabolite, length))
overlap_ratio_metagenome <- unlist(matches_metagenome) / paths_per_compound
overlap_ratio_metagenome[order(overlap_ratio_metagenome, decreasing = T)]

overlap_ratio_mtt <- unlist(matches_mtt) / paths_per_compound
overlap_ratio_mtt[order(overlap_ratio_mtt, decreasing = T)]

compounds <- c("C06573", "C06573", "C00624", "C15987", "C08580") #both metagenome and metatranscriptome

getPathwayName <- function(path_id){
  path_id <- gsub("ko", "", path_id)
  path_name <- NA
  tryCatch({
    entry <- keggGet(paste("map", path_id, sep = ""))
    path_name <- entry[[1]]$NAME
    class_name <- entry[[1]]$CLASS
  }, error = function(e){
    print(e)
  })
  return(path_name)
}
path_names <- lapply(matches_mtt, function(paths){
  sapply(names(paths), getPathwayName)
})
paste(names(path_names[[6]]), path_names[[6]], sep = ":")

sapply(compounds, function(comp){
  sapply(kos_associated_metabolite[[comp]], getPathwayName)
}) 




### Demo code for permutation test examples
df <- data.frame(a = c(0.5,2,3,4,6), b = c(2, 2.5, 3.5, 4.5, 7))
ggplot(melt(df), aes(fill = variable, y = value, x = variable)) + geom_boxplot() + geom_jitter(position = position_jitter(width = .00))

wilcox.test(df$a, df$b)
t.test(df$a, df$b)


ggplot(melt(df), aes(fill = variable, y = value, x = variable)) + geom_boxplot() + geom_jitter(position = position_jitter(width = .00))