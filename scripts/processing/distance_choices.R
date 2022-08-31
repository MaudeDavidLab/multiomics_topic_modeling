
runMantelTest <- function(dist1, dist2, samples, null_test = F){
  if(null_test){
    dist1 <- dist1[rownames(dist1) %in% samples, colnames(dist1) %in% samples]
    dist2 <- dist2[rownames(dist2) %in% samples, colnames(dist2) %in% samples]
    print(paste("Row names match: ", sum(rownames(dist1) == rownames(dist2))))
    return(mantel(dist1, dist2))
  }else{
    dist1 = dist1[samples, samples]
    dist2 = dist2[samples, samples]
    print(paste("Row names match: ", sum(rownames(dist1) == rownames(dist2))))
    return(mantel(dist1, dist2))
  }
}



singleDistTest <- function(ps1, ps2, dist1, dist2, output_folder){
  dist_mat1 <- as.data.frame(as.matrix(vegdist(t(otu_table(ps1)), method = dist1)))
  dist_mat2 <- as.data.frame(as.matrix(vegdist(t(otu_table(ps2)), method = dist2)))
  samples <- intersect(sample_names(ps1), sample_names(ps2))
  mantel_res <- runMantelTest(dist_mat1, dist_mat2, samples)
  return(mantel_res)
}

runDistanceTests <- function(output_folder){
  distance_metrics <- c('manhattan', 'euclidean', 'canberra', 'clark', 'bray', 'kulczynski', 'jaccard', 'gower', 'altGower',
                        'mountford', 'raup', 'binomial', 'chao', 'cao', 'mahalanobis')
  for(dist_metric in distance_metrics){
    dist_16s <- as.data.frame(as.matrix(vegdist(t(otu_table(ps_16s)), method = dist_metric)))
    dist_mtg <- as.data.frame(as.matrix(vegdist(t(otu_table(ps_mtg)), method = dist_metric)))
    dist_mtt <- as.data.frame(as.matrix(vegdist(t(otu_table(ps_mtt)), method = dist_metric)))
    dist_metabol <- as.data.frame(as.matrix(vegdist(t(otu_table(ps_metabol)), method = dist_metric)))
    
    mantel_16s_mtg <- runMantelTest(dist_16s, dist_mtg, samples_16s_mtg)
    mantel_16s_mtt <- runMantelTest(dist_16s, dist_mtt, samples_16s_mtt)
    mantel_16s_metabol <- runMantelTest(dist_16s, dist_metabol, samples_16s_metabol)
    mantel_mtg_mtt <- runMantelTest(dist_mtg, dist_mtt, samples_mtg_mtt)
    mantel_mtg_metabol <- runMantelTest(dist_mtg, dist_metabol, samples_mtg_metabol)
    mantel_mtt_metabol <- runMantelTest(dist_mtt, dist_metabol, samples_mtt_metabol)
    
    saveRDS(mantel_16s_mtg, paste(output_folder, dist_metric, '_16s_mtg.rds', sep = ""))
    saveRDS(mantel_16s_mtt, paste(output_folder, dist_metric, '_16s_mtt.rds', sep = ""))
    saveRDS(mantel_16s_metabol, paste(output_folder, dist_metric, '_16s_metabol.rds', sep = ""))
    saveRDS(mantel_mtg_mtt, paste(output_folder, dist_metric, '_mtg_mtt.rds', sep = ""))
    saveRDS(mantel_mtg_metabol, paste(output_folder, dist_metric, '_mtg_metabol.rds', sep = ""))
    saveRDS(mantel_mtt_metabol, paste(output_folder, dist_metric, '_mtt_metabol.rds', sep = ""))
  }
  
}

printTests <- function(dist_metric, folder){
  print("16S vs mtg")
  print(readRDS(paste(folder, "/", dist_metric, '_16s_mtg.rds', sep = ""))$signif)
  print("16S vs mtt")
  print(readRDS(paste(folder, "/",dist_metric, '_16s_mtt.rds', sep = ""))$signif)
  print("16S vs metabol")
  print(readRDS(paste(folder, "/",dist_metric, '_16s_metabol.rds', sep = ""))$signif)
  print("mtg vs mtt")
  print(readRDS(paste(folder,"/", dist_metric, '_mtg_mtt.rds', sep = ""))$signif)
  print("mtg vs metabol")
  print(readRDS(paste(folder,"/", dist_metric, '_mtg_metabol.rds', sep = ""))$signif)
  print("mtt vs metabol")
  print(readRDS(paste(folder, "/",dist_metric, '_mtt_metabol.rds', sep = ""))$signif)
}

runNullVerify <- function(ps1, ps2, samples, dist_metric){
  set.seed(123)
  permuted_1 <- t(otu_table(ps1))
  permuted_1 <- permuted_1[sample(seq(1, nrow(permuted_1)), nrow(permuted_1)), ]
  permuted_2 <- t(otu_table(ps2))
  permuted_2 <- permuted_2[sample(seq(1, nrow(permuted_2)), nrow(permuted_2)), ]
  
  dist_1_null <- as.data.frame(as.matrix(vegdist(permuted_1, method = dist_metric)))
  dist_2_null <- as.data.frame(as.matrix(vegdist(permuted_2, method = dist_metric)))
  mantel_null <- runMantelTest(dist_1_null, dist_2_null, samples, null_test = T)
  return(mantel_null)
}

printNullTests <- function(dist_metric){
  print("16S vs mtg")
  print(runNullVerify(ps_16s, ps_mtg, samples_16s_mtg, dist_metric))
  print("16S vs mtt")
  print(runNullVerify(ps_16s, ps_mtt, samples_16s_mtt, dist_metric))
  print("16S vs metabol")
  print(runNullVerify(ps_16s, ps_metabol, samples_16s_metabol, dist_metric))
  print("mtg vs mtt")
  print(runNullVerify(ps_mtg, ps_mtt, samples_mtg_mtt, dist_metric))
  print("mtg vs metabol")
  print(runNullVerify(ps_mtg, ps_metabol, samples_mtg_metabol, dist_metric))
  print("mtt vs metabol")
  print(runNullVerify(ps_mtt, ps_metabol, samples_mtt_metabol, dist_metric))
  
}

#Average distance between samples in any dataset comparison
samples_16s_mtg <- intersect(sample_names(ps_16s), sample_names(ps_mtg))
samples_16s_mtt <- intersect(sample_names(ps_16s), sample_names(ps_mtt))
samples_16s_metabol <- intersect(sample_names(ps_16s), sample_names(ps_metabol))
samples_mtg_mtt <- intersect(sample_names(ps_mtg), sample_names(ps_mtt))
samples_mtg_metabol <- intersect(sample_names(ps_mtg), sample_names(ps_metabol))
samples_mtt_metabol <- intersect(sample_names(ps_mtt), sample_names(ps_metabol))

folder = "../results/summary_statistics/distance_choices/mantel_tests/rle/"
runDistanceTests(output_folder = folder)
# We need some criteria for picking a distance metric. I'm declaring that 
#1. 16s has to match mtg
#2. mtg has to match mtt
# Starred metrics match these criteria

#Print pvalues : low value means the two omics are significantly related
printTests("euclidean", folder = folder)
printTests("canberra", folder = folder) 
printTests("clark", folder = folder) 
printTests('bray', folder = folder) 
printTests('kulczynski', folder = folder) 
printTests('jaccard', folder = folder) 
printTests('gower', folder = folder) 
printTests('altGower', folder = folder)
printTests('morisita', folder = folder) 
printTests('horn', folder = folder)
#printTests('mountford', folder = 'output/mantel_tests/rle/')
#printTests('raup', folder = 'output/mantel_tests/rle/')
#There was an error in binomial so the rest didn't get calculated.

# Just to be super sure, we can permute the samples in the table to make sure the result is NOT significant when it shouldn't be
printNullTests('horn') # none are significant, which is good.
printNullTests('morisita')# none are significant, which is good.


## Confirm distance metric and normalization choices
dist_16s <- as.data.frame(as.matrix(vegdist(t(otu_table(ps_16s)), method = "horn")))
dist_mtg <- as.data.frame(as.matrix(vegdist(t(otu_table(ps_mtg)), method = "horn")))
dist_mtt <- as.data.frame(as.matrix(vegdist(t(otu_table(ps_mtt)), method = "horn")))
dist_metabol <- as.data.frame(as.matrix(vegdist(t(otu_table(ps_metabol)), method = "euclidean")))

getSampleOverlap <- function(df1, df2){
  return(intersect(rownames(df1), rownames(df2)))
}
mantel_16s_mtg <- runMantelTest(dist_16s, dist_mtg, getSampleOverlap(dist_16s, dist_mtg))
mantel_16s_mtt <- runMantelTest(dist_16s, dist_mtt, getSampleOverlap(dist_16s, dist_mtt))
mantel_16s_metabol <- runMantelTest(dist_16s, dist_metabol, getSampleOverlap(dist_16s, dist_metabol))
mantel_mtg_mtt <- runMantelTest(dist_mtg, dist_mtt, getSampleOverlap(dist_mtg, dist_mtt))
mantel_mtg_metabol <- runMantelTest(dist_mtg, dist_metabol, getSampleOverlap(dist_mtg, dist_metabol))
mantel_mtt_metabol <- runMantelTest(dist_mtt, dist_metabol, getSampleOverlap(dist_mtt, dist_metabol))

output_folder = "../results/summary_statistics/distance_choices/mantel_tests/"
saveRDS(mantel_16s_mtg, paste(output_folder, 'best_16s_mtg.rds', sep = ""))
saveRDS(mantel_16s_mtt, paste(output_folder,  'best_16s_mtt.rds', sep = ""))
saveRDS(mantel_16s_metabol, paste(output_folder,  'best_16s_metabol.rds', sep = ""))
saveRDS(mantel_mtg_mtt, paste(output_folder,  'best_mtg_mtt.rds', sep = ""))
saveRDS(mantel_mtg_metabol, paste(output_folder,  'best_mtg_metabol.rds', sep = ""))
saveRDS(mantel_mtt_metabol, paste(output_folder,  'best_mtt_metabol.rds', sep = ""))

dist_metabol <- as.data.frame(as.matrix(vegdist(t(otu_table(ps_metabol)), method = "euclidean"))) 
mantel_res = singleDistTest(ps_16s, ps_mtg, "horn", "horn")
mantel_res$signif
mantel_res = singleDistTest(ps_16s, ps_mtg, "horn", "horn")
mantel_res$signif



tmp <- readRDS(paste0(output_folder, 'best_mtg_mtt.rds'))

