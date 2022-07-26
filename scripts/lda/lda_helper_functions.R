getSimulatedWord <- function(sample_topic_dist, topic_taxa_dist){
  topic <- sample(seq(length(sample_topic_dist)), prob = sample_topic_dist, size = 1)
  taxa <- sample(seq(ncol(topic_taxa_dist)), prob = topic_taxa_dist[topic, ], size = 1)
  return(taxa)
}

getSimulatedSample <- function(sample_topic_dist, topic_taxa_dist, n){
  sim_count <- sapply(seq(n), function(x) return(getSimulatedWord(sample_topic_dist, topic_taxa_dist)))
  sim_count <- table(sim_count)
  missing <- seq(1, ncol(topic_taxa_dist))[!seq(1, ncol(topic_taxa_dist)) %in% names(sim_count)]
  tmp <- rep(0, length(missing))
  names(tmp) <- missing
  sim_count <- c(sim_count, tmp)
  return(sim_count[order(as.numeric(names(sim_count)))])
}

getSimulatedDataset <- function(ps, model, n = 100, divisor = 1){
  thetas <- model@gamma
  betas <- exp(model@beta)
  samples <- sample_names(ps)
  
  rownames(thetas) <- samples
  sample_counts <- colSums(ps@otu_table)
  
  #sample_counts[sample]
  if(n == "sample_counts"){
    sim_list <- lapply(samples, function(sample) return(getSimulatedSample(thetas[sample, ], betas, n = sample_counts[sample]/divisor)))
  }else{
    sim_list <- lapply(samples, function(sample) return(getSimulatedSample(thetas[sample, ], betas, n = n)))
  }
  
  
  sim_seq_tab <- do.call("rbind", sim_list)
  rownames(sim_seq_tab) <- sample_names(ps)
  colnames(sim_seq_tab) <- taxa_names(ps)
  return(sim_seq_tab)
}

filterByPrevalence <- function(ps, filter_thresh){
  data <- round(ps@otu_table)
  if(taxa_are_rows(ps)){
    keep <- apply(data, 1, function(x) return(sum(x>0) > filter_thresh*nsamples(ps)))
    data <- data[keep, ]
    ps_new <- phyloseq(otu_table(data, taxa_are_rows = T), sample_data(ps@sam_data), tax_table(ps@tax_table))
  }else{
    keep <- apply(data, 2, function(x) return(sum(x > 0)> filter_thresh*nsamples(ps)))
    data <- data[ , keep]
    ps_new <- phyloseq(otu_table(data, taxa_are_rows = F), sample_data(ps@sam_data), tax_table(ps@tax_table))
  }
  return(ps_new)
}


procrustesTest <- function(data_true, data_sim, k = 50, scale = FALSE, title = ""){
  scores_true <- svd(scale(data_true, scale = scale))$u[ , seq(k)]
  loadings_true <-  svd(scale(data_true, scale = scale))$v[ , seq(k)]
  evals_true <-  svd(scale(data_true, scale = scale))$d[k]
  
  scores_sim <- svd(scale(data_sim, scale = scale))$u[ , seq(k)]
  loadings_sim <-  svd(scale(data_sim, scale = scale))$v[ , seq(k)]
  evals_sim <-  svd(scale(data_sim, scale = scale))$d[k]
  
  aligned_scores <- procrustes(scores_true, scores_sim)$Yrot
  aligned_loadings <- procrustes(loadings_true, loadings_sim)$Yrot
  
  tmp <- procrustes(as.matrix(scores_true), as.matrix(scores_sim))
  
  plot(tmp, main = title, type= "text")
  plot(tmp, kind = 2, main = title)
  protest(X = as.matrix(scores_true), Y = as.matrix(scores_sim), scores = "sites", permutations = 999)
  
  
  keep <- tmp[[1]][,2] > -0.1 & tmp[[1]][,1] < 0.2
  data_true <- data_true[keep, ]
  data_sim <- data_sim[keep, ]
  tmp[[1]] <- tmp[[1]][tmp[[1]][,2] > -0.1 & tmp[[1]][,1] < 0.2, ]
  plot(tmp, main = title, type= "text")
  return(tmp)
}

plotNumberTopics <- function(directory, prefix, title){
  likelihoods <- list(rep(NA, 9))
  for(file in list.files(directory)){
    i <- sub(".*?_", "", file)
    i <- sub(".*?_", "", i)
    i <- sub("topics.*", "", i)
    #i <- sub("topics_3000iter_0.05filter.rds", "", i)
    i <- as.numeric(i)-1
    model <- readRDS(paste0(directory, file))
    likelihoods[[i]] <- model@loglikelihood
    rm(model)
  }
  
  log_likelihoods <- unlist(likelihoods)
  df <- data.frame(numTopics = seq(2, length(log_likelihoods)+1), likelihoods = log_likelihoods)
  p <- ggplot(df, aes(x = numTopics, y = likelihoods))+
    geom_point()+
    geom_line()+
    xlab("Number of Topics") + ylab("Log Likelihood")+
    theme_bw()
  return(p)
}

quantileTest <- function(data_true, data_sim, type ="samples"){
  if(type == "samples"){
    quantiles_true <- apply(data_true, 1, quantile, probs = seq(0, .98, 0.05))
    quantiles_sim <- apply(data_sim, 1, quantile, probs = seq(0, .98, 0.05))
  }else if(type == "features"){
    quantiles_true <- apply(data_true, 2, quantile, probs = seq(0, .98, 0.05))
    quantiles_sim <- apply(data_sim, 2, quantile, probs = seq(0, .98, 0.05))
  }
  
  melted_true <- melt(quantiles_true)
  melted_true$type <- rep("true", nrow(melted_true))
  colnames(melted_true) <- c("quantile", "sample", "quantile_value", "type")
  
  melted_sim <- melt(quantiles_sim)
  melted_sim$type <- rep("sim", nrow(melted_sim))
  colnames(melted_sim) <- c("quantile", "sample", "quantile_value", "type")
  
  df <- data.frame(cbind(melted_true, melted_sim))
  line <- lm(quantile_value ~ quantile_value.1, df)
  p <- ggplot(df, aes(x = quantile_value, y = quantile_value.1)) +
    geom_point() +
    geom_smooth(method='lm', formula= y~x)+
    ggplot2::annotate("text", size = 4, x = min(df$quantile_value) + sd(df$quantile_value) + 2, y = max(df$quantile_value.1) + sd(df$quantile_value.1), label = paste0("R2 = ", round(summary(line)[[9]], 3)))+
    theme(axis.text = element_text(size = 6.5))+
    theme_bw()+
    xlab("Quantile Values True") + ylab("Quantile Values Simulated")+
    ylim(c(0,  max(df$quantile_value.1) + sd(df$quantile_value.1)))
  print(p)
  return(list(p=p, line=line))
  
}

marginalDistributionCorrelation <- function(data_true, data_sim){
  # feature by feature correlation matrix
  f_corr_true <- cor(data_true, method = "spearman")
  f_corr_sim <- cor(data_sim, method = "spearman")
  c <- cor(c(f_corr_true), c(f_corr_sim), method = "spearman") #correlation between the pairwise marginal distribution
  plot(c(f_corr_true), c(f_corr_sim))
  return(c)
}

getModelFit <- function(ps, sim_counts, title = ""){
  data_true = t(ps@otu_table)
  data_sim = sim_counts
  pro_res <- procrustesTest(data_true, data_sim, scale = F, k = 30, title = title)
  p3 <- quantileTest(data_true, data_sim)
  return(list(pro_res, p3))
}


getCorrelationsMetadata <- function(ps, thetas, ending,  n ){
  
  rownames(thetas) == sample_names(ps@sam_data)
  
  df <- data.frame(cbind(thetas, ps@sam_data))
  df$stool_freq <- as.numeric(ps@sam_data$stool_freq)
  df$cat <- as.integer(ps@sam_data$cat)
  df$phenotype <- as.integer(ifelse(ps@sam_data$phenotype == "A", 1, 0))
  df$csection <- as.integer(as.logical(ps@sam_data$csection))
  df$age <- as.numeric(df$age)
  colnames(df)[1:n] <- colnames(thetas)
  topics <- colnames(df)[1:n]
  
  control_variables <- c("whole_grain", "fermented_vegetables", "dairy", "fruit", "meal_home_prep", "meal_ready_eat", "meat", "olive_oil", "seafood", "sweetened_drink", "vegetable", "restaurant", "sugary_food", "starchy_food", "dairy_freq", "fat_oil_freq", "vegetable_freq", "fruit_freq", "phenotype", "probiotic", "csection", "vitamin_B", "vitamin_D", "age")
  
  correlations <- lapply(seq(1, n), function(i){
    topic <- colnames(df)[i]
    print(topic)
    topic_vec <- df[, i]
    if(sum(is.na(topic_vec)) == 0){
      cor_res <- lapply(control_variables, function(x){
        print(x)
        tmp <- df[!is.na(df[ , x]), ] 
        var_vec <- as.numeric(tmp[,x])#variable vector
        topic_vec <- tmp[ , i]
        return(cor.test(var_vec, topic_vec, method = "spearman"))
      })
      pvals <- unlist(lapply(cor_res, function(x) return(x$p.value)))
      names(pvals) <- control_variables
      rvals <- unlist(lapply(cor_res, function(x) return(x$estimate)))
      names(rvals) <- control_variables
      p_adj <- p.adjust(pvals, method = "fdr")
      print(sum(p_adj < .05) / length(p_adj))
      print(control_variables[p_adj < 0.05 & abs(rvals) > 0.1])
      return(list(rvals = rvals, p_adj = p_adj))
    }else{
      return(list(rvals = NA, p_adj = NA))
    }
    
  })
 # names(correlations) <- topics
  
  rvals_mat <- do.call(rbind, lapply(correlations, function(x){ return(x$rvals)}))
  rownames(rvals_mat) <- topics
  pvals_mat <- do.call(rbind, lapply(correlations, function(x){ return(x$p_adj)}))
  rownames(pvals_mat) <- topics
  return(list(rvals = rvals_mat, p_adj = pvals_mat))
}


getImportantFactors <- function(ps, data, iter = 5){
  imp <- c()
  for(i in seq(1, iter)){
    print(i)
    set.seed(i)
    group <- data.frame(t(ps@otu_table)[rownames(data), ])
    group$phenotype <- as.factor(ifelse(data$phenotype == 1, "A", "N"))
    out <- Boruta(group %>% select(-c("phenotype")), group$phenotype)
    out
    df_plot<- group[ , names(out$finalDecision[out$finalDecision == "Confirmed" | out$finalDecision == "Tentative"]), drop = F]
    imp <- c(imp, colnames(df_plot))
  }
  df_plot <- data.frame(t(ps@otu_table)[rownames(data), ])
  #keep <- names(table(imp))[table(imp) > 25]
  #df_plot<- group[, keep, drop = F]
  df_plot$phenotype <- as.factor(ifelse(data$phenotype == 1, "A", "N"))
  return(list(df_plot= df_plot, imp = imp))
}

keepMostImp <- function(df_plot,imp,num){
  phenotype <- df_plot$phenotype
  keep <- names(table(imp))[table(imp) >= num]
  df_plot <- df_plot[, keep, drop = F]
  df_plot$phenotype <- phenotype
  return(df_plot)
}
