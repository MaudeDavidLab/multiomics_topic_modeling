

addNumericMetadata <- function(ps){
  mapping_ml_numeric <- read.csv("data/mapping_ml_numeric_full.csv")
  rownames(mapping_ml_numeric) <- mapping_ml_numeric$biospecimen_id
  sample_data(ps) <- mapping_ml_numeric
  familyID <- ps@sam_data$familyID
  host_name <- ps@sam_data$host_name
  phenotype <- ps@sam_data$phenotype
  sample_data(ps) <- ps@sam_data[ , !(sample_variables(ps) %in% c("mother_age", "father_age", "age", "racial_group",
                                                                  "timepoint", "host_name", "biospecimen_id", "familyID","min_time_antibiotics"))]
  
  #Normalize
  df <- ps@sam_data
  na_count <- apply(df, 2, function(x) return(sum(is.na(x))))
  df <- df[ , na_count < nsamples(ps) * .5]
  for(i in 1:ncol(df)){
    x <- as.numeric(df[,i][[1]])
    df[ , i] <- as.numeric(x) / max(x, na.rm = T) #make everything between 0 and 1
  }
  
  
  
  sample_data(ps) <- df
  sample_data(ps)$familyID <- familyID
  sample_data(ps)$host_name <- host_name
  sample_data(ps)$phenotype <- phenotype
  return(ps)
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
normalize_rle <- function(df, addOne = F){
  #assume variable by sample
  df <- df + 1
  geoMeans <- apply(df, 1, gm_mean) #geometric mean of each taxa
  norm_factor <- apply(df, 2, function(sample){
    median(sample / geoMeans)
  })
  return(t(t(df) / norm_factor))
}

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

getKOName <- function(ko_id){
  tryCatch({
    entry <- keggGet(ko_id)[[1]]
    def <- entry$DEFINITION
    return(def)
  }, error = function(e){
    print(e)
  })
}

findSig_Binomial <- function(actual_diffs, pos, type = "ko", thresh = 0.05){
  pvals <- c()
  pvals <- apply(actual_diffs, 2, function(diffs){
    if(pos){
      observed_successes <- sum(diffs > 0)
    }else{
      observed_successes <- sum(diffs < 0)
    }
    pval <- 1 - sum(dbinom(seq(0, observed_successes), p = 0.5, size = nrow(actual_diffs))) 
    #pval <- sum(dbinom(seq(observed_successes, nrow(actual_diffs)), p = 0.5, size = nrow(actual_diffs))) # what we saw or crazier
    return(pval)
  })
  pvals <- p.adjust(pvals, method = "fdr")
  print(pvals)
  hist(pvals)
  ko_ids <- colnames(actual_diffs)
  df <- data.frame(ko_ids, pvals)
  #ko_ids <- colnames(actual_diffs)[which(pvals < thresh)]
  #df <- data.frame(ko_ids, pvals[pvals< thresh])
  return(df)
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


findSig_Binomial_old <- function(null_dist, actual_diffs, pos = T, type = "ko", thresh = 0.05){
  pvals <- c()
  for(ko in colnames(null_dist)){
    if(pos){
      null_prob <- sum(null_dist[ ,ko] > 0)/nrow(null_dist) #probability that you observe something above 0 (i.e prob. that one person has more given randomness)
      observed_successes <- sum(actual_diffs[,ko] > 0) #what we actually observed
    }else if(!pos){
      null_prob <- sum(null_dist[ ,ko] < 0)/nrow(null_dist) #probability that you observe something above 0 (i.e prob. that one person has more given randomness)
      observed_successes <- sum(actual_diffs[,ko] < 0) #what we actually observed
    }
    
    pval <- 1 - sum(dbinom(seq(0, observed_successes-1), p = null_prob, size = length(actual_diffs[,ko]))) # 1 - less crazy = what we saw or crazier
    print(paste("Null prob: ", null_prob))
    pvals <- c(pvals, pval)
  }
  p_adj <- p.adjust(pvals, method = "fdr")
  ko_ids <- colnames(null_dist)[which(p_adj < thresh)]
  df <- data.frame(ko_ids, p_adj[p_adj< thresh])
  return(df)
}


#if(type == "pathway"){
#  path_names <- sapply(ko_ids, getPathwayName)
#sapply(ko_ids, function(ko_id){
#  plotNullHist(null_dist[, ko_id], actual_diffs[,ko_id], ko = ko_id, path = getPathwayName(ko_id))
#})
# }
# if(type == "metabolite"){
#   path_names <- metabolite2pathway[ko_ids, ]
#sapply(ko_ids, function(ko_id){
#  plotNullHist(null_dist[, ko_id], actual_diffs[,ko_id], ko = ko_id, path = metabolite2pathway[ko_id, ])
#})
#}
# if(type == "ko"){
#   path_names <- sapply(ko_ids, getKOName)
#sapply(ko_ids, function(ko_id){
#  plotNullHist(null_dist[, ko_id], actual_diffs[,ko_id], ko = ko_id, path = getKOName(ko_id))
#})
# }


getNullDist <- function(ps, iter = 10000){
  seqtab <- ps@otu_table
  seqtab <- t(seqtab)
  null_dist <- list()
  for(i in seq(1, iter)){
    if(i %% 100 == 0){
      print(i)
    }
    sample_pair <- sample(unique(ps@sam_data$host_name), 2)
    while(substr(sample_pair[1], 1, 3) == substr(sample_pair[2], 1, 3)){
      sample_pair <- sample(unique(ps@sam_data$host_name), 2)
    }
    sample_names <- sample_names(ps)[ps@sam_data$host_name %in% sample_pair]
    df_tmp <- seqtab[as.character(sample_names), ]
    null_dist[[i]] <-  (df_tmp[1,] - df_tmp[2,])
  }
  
  null_dist <- do.call(rbind, null_dist)
  colnames(null_dist) <- colnames(seqtab)
  return(null_dist)
}


#We need to plot them on the same frequency scale, ideally still in a histogram
plotNullHist <- function(g, actual_diffs_vec = c(0), ko = NA, path = NA){
  g <- as.numeric(g)
  actual_diffs_vec <- as.numeric(actual_diffs_vec)
  minimum <-  min(min(actual_diffs_vec, min(g)))
  maximum <- max(max(actual_diffs_vec, max(g)))
  par(mar = c(5,5,2,5))
  h <- hist(g, breaks = 40,  main = paste(ko, ":", path))
  
  xfit <- seq(min(g), max(g), length = 40) 
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g)) 
  yfit <- yfit * diff(h$mids[1:2]) * length(g) 
  lines(xfit, yfit, col = "black", lwd = 2,  xlim =c(minimum, maximum))
  for(x in actual_diffs_vec){
    abline(v = x, col = "red")
  }
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

linearModeling_weightMetadata <- function(ps, type, reg = 2.3, vars = NA){
  if(is.na(vars)){
    vars <- sample_variables(ps)
  }
  if(!("familyID" %in% vars)){
    vars <- c(vars, "familyID")
  }
  if(!("phenotype" %in% vars)){
    vars <- c(vars, "phenotype")
  }
  if(!("host_name" %in% vars)){
    vars <- c(vars, "host_name")
  }

  sample_data(ps) <- ps@sam_data[ , vars]

  ps <- filterNAs(ps)
  df <- ps@sam_data
  df <- df %>% select(-c("familyID", "phenotype", "host_name")) #don't count in metadata dist
  meta_dist <- as.matrix(dist(df, method = "euclidean"))
  ps <- prune_samples(rownames(meta_dist), ps)
  seqtab <- t(ps@otu_table)
  
  preds <- c()
  pos_kos <- c()
  neg_kos <- c()
  for(sam_id in rownames(df)){
    print(sam_id)
    fam_id <- ps@sam_data[sam_id, "familyID"][[1]]
    train_data <- seqtab[ps@sam_data$familyID != fam_id, ]
    weights <- meta_dist[sam_id, rownames(train_data)]
    y <- ps@sam_data[rownames(train_data), "phenotype"][[1]]
    y <- ifelse(y == "A", 1, 0)
    
    data <- as.data.frame(cbind(train_data, y))
    model <- cv.glmnet(x = train_data, y = as.factor(y),
                       family = "binomial", type.measure = "auc", weights = 1/(weights + .01))
    l <- reg * sd(model$lambda) + mean(model$lambda)
    coefs <- coefficients(model, s = l)
    pos_coefs <- rownames(coefs)[as.numeric(coefs) > 0]
    neg_coefs <- rownames(coefs)[as.numeric(coefs) < 0]
    print(pos_coefs)
    pos_kos <- c(pos_kos, pos_coefs)
    neg_kos <- c(neg_kos, neg_coefs)
    pred <-predict(model, newx = t(ps@otu_table[ , sam_id]),
                   s = l, type = "response")
    preds <- c(preds, pred)
    #print(paste("prediction: ", pred))
    #print(paste("Actual: ", ps@sam_data[sam_id, "phenotype"][[1]]))
  }
  preds_num <- ifelse(preds <0.5, 0, 1)
  phenotype_num <- ifelse(ps@sam_data$phenotype == "A", 1, 0)
  sum(preds_num == phenotype_num) / nsamples(ps)
  
  roc_obj <- roc(response = phenotype_num, predictor = preds)
  print(roc_obj)
  table(pos_kos)
  table(neg_kos)
  sig_kos_pos_linear <- data.frame(table(pos_kos))
  sig_kos_neg_linear <- data.frame(table(neg_kos))
  
  #if(type == "pathway"){
  #  sig_kos_pos_linear <- data.frame(table(pos_kos), sapply(names(table(pos_kos)), getPathwayName))
  #  sig_kos_neg_linear <- data.frame(table(neg_kos), sapply(names(table(neg_kos)), getPathwayName))
  #}
  #if(type == "metabolite"){
  #  sig_kos_pos_linear <- data.frame(table(pos_kos), pathway = sapply(names(table(pos_kos)), function(x) return(as.character(metabolite2pathway[x, ]))))
  #  sig_kos_neg_linear <- data.frame(table(neg_kos), pathway = sapply(names(table(neg_kos)), function(x) return(as.character(metabolite2pathway[x, ]))))
  #}
  #if(type == "ko"){
  #  sig_kos_pos_linear <- data.frame(table(pos_kos), sapply(names(table(pos_kos)), getKOName))
  #  sig_kos_neg_linear <- data.frame(table(neg_kos),  sapply(names(table(neg_kos)), getKOName))
  #}
  return(list(aut = sig_kos_pos_linear, nt = sig_kos_neg_linear))
}

filterNAs <- function(ps){
  df <- sample_data(ps)
  df <- df[ , colSums(is.na(df)) < 100] #delete variables that everyone has as NA
  df <- na.omit(df) #delete samples that have too many NAs
  sample_data(ps) <- df
  return(ps)
}

runPermanova <- function(ps_tmp){
  
  ps_permanova <- filterNAs(ps_tmp)
  familyID <- ps_permanova@sam_data$familyID
  ps_permanova@sam_data <- ps_permanova@sam_data %>% select(-c("host_name", "familyID"))
  ps_permanova@sam_data$phenotype <- ifelse( ps_permanova@sam_data$phenotype == "A", 1, 0)
  sample_data(ps_permanova) <- data.frame(data.matrix(ps_permanova@sam_data)) #convert to numeric in case characters
  dists <- vegdist(t(ps_permanova@otu_table), method = "bray")

  sam_dat <- ps_permanova@sam_data
  sam_dat$familyID <- familyID
  
  form <- formula(paste("dists" , paste(sample_variables(ps_permanova), collapse = "+"), sep = " ~ "))
  perm <- how(nperm = 10000)
  setBlocks(perm) <- with(sam_dat, familyID)
  perm <- adonis(form, data = data.frame(ps_permanova@sam_data), permutations = perm)
  perm_tab <- as.data.frame(perm$aov.tab)
  perm_tab <- perm_tab[!(rownames(perm_tab) %in% c("Total", "Residuals")), ]
  
  items <- rownames(perm_tab)
  disp_p_vals <- list()
  for(item in items){
    p <- permutest(betadisper(dists, group = ps_permanova@sam_data[ , item][[1]]))
    disp_p_vals[[item]] <-  p$tab$`Pr(>F)`
  }
  disp_p_vals <- unlist(lapply(disp_p_vals, function(x) return(x[1])))
  
  perm_sig <- cbind(perm_tab, disp_p_vals)
  colnames(perm_sig)[7] <- "disp_pval"
  return(perm_sig)
}


trainRf <- function(input, map){
  folds <- groupKFold(map$familyID, k = 5) 
  train_control <- trainControl(method="cv", index = folds, repeats = 3, search = "grid")
  set.seed(10)
  tunegrid <- expand.grid(.mtry=c(1:15))
  model <- train(x = input, y = as.factor(map$phenotype), trControl=train_control, method="rf", tuneGrid = tunegrid, metric = "Accuracy")
  # summarize results
  return(model)
}