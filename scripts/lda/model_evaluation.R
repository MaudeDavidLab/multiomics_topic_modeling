library(phyloseq)

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
  return(c)
}

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
  
  sim_list <- lapply(samples, function(sample){
    total_counts <- sample_counts[sample]
    
    # For every sample, for every word that should be in that sample, pick a topic. Count the number of times each topic is picked
    topics <- rmultinom(n = 1, size = total_counts, prob = thetas[sample, ])
    
    # Pick a word for each topic directly from the number of times that topic was selected above
    words <- sapply(seq(1, length(topics)), function(i){
      counts_for_topic <- topics[i]
      feature_counts <- rmultinom(n = 1, size = counts_for_topic, prob = betas[i, ])
    })
    words <- rowSums(words)
    return(words)
  })
  
  sim_seq_tab <- do.call("rbind", sim_list)
  rownames(sim_seq_tab) <- sample_names(ps)
  colnames(sim_seq_tab) <- taxa_names(ps)
  return(sim_seq_tab)
}


plotQuantileCorrelations <- function(ps, directory, title, sim = F, divisor = 1){
  files <- list.files(directory)
  files <- files[grepl("model", files)]
  r2s_samples <- c()
  r2s_features <- c()
  pmdc <- c()
  topics <- c()
  for(file in files){
    i <- sub(".*?_", "", file)
    i <- sub(".*?_", "", i)
    i <- sub("topics.*", "", i)
    #i <- sub("topics_3000iter_0.05filter.rds", "", i)
    i <- as.numeric(i)
    print(i)
    topics <- c(topics, i)
    model <- readRDS(paste0(directory, file))
    # likelihoods[[i]] <- model@loglikelihood
    if(sim){
      sim_counts <- getSimulatedDataset(ps, model, n = "sample_counts", divisor = divisor)
      saveRDS(sim_counts, paste0(directory, "/sim_counts_", i, "topics.rds"))
    }else{
      sim_counts <- readRDS(paste0(directory, "/sim_counts_", i, "topics.rds"))
    }
    # sample quantiles
    tmp <- quantileTest(t(ps@otu_table), sim_counts, type = "samples")
    line <- tmp$line
    r2s_samples <- c(r2s_samples, summary(line)[[9]])
    
    # feature quantiles
    tmp <- quantileTest(t(ps@otu_table), sim_counts, type = "features")
    line <- tmp$line
    r2s_features <- c(r2s_features, summary(line)[[9]])
    
    pairwise_marg_dist_corr <- marginalDistributionCorrelation(t(ps@otu_table), sim_counts)
    pmdc <- c(pmdc, pairwise_marg_dist_corr)
    
    rm(model)
    rm(sim_counts)
  }
  df <- data.frame(numTopics = topics, r2s_samples = r2s_samples, r2s_features = r2s_features, pmdc = pmdc)
  p_samples <- ggplot(df, aes(x = numTopics, y = r2s_samples)) +
    geom_point()+ 
    geom_line()+
    theme_bw()+
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14))+
    xlab("Number of Topics") + ylab("R^2")+
    ggtitle(paste(title, "Quantile Correlation Samples"))
  
  p_features <- ggplot(df, aes(x = numTopics, y = r2s_features)) +
    geom_point()+ 
    geom_line()+
    theme_bw()+
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14))+
    xlab("Number of Topics") + ylab("R^2")+
    ggtitle(paste(title, "Quantile Correlation Features"))
  
  p_pmdc <- ggplot(df, aes(x = numTopics, y = pmdc)) +
    geom_point()+ 
    geom_line()+
    theme_bw()+
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14))+
    xlab("Number of Topics") + ylab("R^2")+
    ggtitle(paste(title, "Pairwise Marginal Dist. Correlation"))
  return(list(p_samples = p_samples, p_features = p_features, p_pmdc = p_pmdc))
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

ps_mtg <- readRDS("ps_mtg_rle_nooutliers_adjcounts.rds")
ps_mtg_filt <- filterByPrevalence(ps_mtg, 0.1)
p_mtg <- plotQuantileCorrelations(ps_mtg_filt, "models_mtg_50000/", title = "MTG", sim = T, divisor = 1)
p <- ggarrange(p_mtg$p_samples, p_mtg$p_features, p_mtg$p_pmdc)
ggsave("models_mtg_50000/quantile_correlation_num_topics.pdf", p, width = 6, height = 3.5)


ps_mtt <- readRDS("ps_mtt_rle_nooutliers_adjcounts.rds")
ps_mtt_filt <- filterByPrevalence(ps_mtt, 0.1)
p_mtt <- plotQuantileCorrelations(ps_mtt_filt, "models_mtt_50000/", title = "MTT", sim = T, divisor = 1)
p <- ggarrange(p_mtt$p_samples, p_mtt$p_features, p_mtt$p_pmdc)
ggsave("models_mtt_50000/quantile_correlation_num_topics.pdf", p, width = 6, height = 3.5)


ps_mbx <- readRDS("ps_mbx_rle_nooutliers_adjcounts_fixedmapping.rds")
ps_mbx_filt <- filterByPrevalence(ps_mbx, 0.1)
p_mbx <- plotQuantileCorrelations(ps_mbx_filt, "models_mbx_50000/", title = "MBX", sim = T, divisor = 1)
p <- ggarrange(p_mbx$p_samples, p_mbx$p_features, p_mbx$p_pmdc)
ggsave("models_mbx_50000/quantile_correlation_num_topics.pdf", p, width = 6, height = 3.5)