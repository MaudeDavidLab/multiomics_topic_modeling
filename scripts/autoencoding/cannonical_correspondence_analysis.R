library(vegan)
library(phyloseq)
library(dplyr)
library(lme4)
library(caret)
library(DESeq2)
library(ape)
library(reshape2)
library(randomForest)
library(caret)
library(dplyr)
library(KEGGREST)
library(pROC)
library(glmnet)
library(patchwork)


setwd("C:/Users/kitikomp/Documents/Lab/M3")
source('scripts/clean_mapping_ml.R')
ps_16s <- readRDS("data/16s/ps_dds_nooutliers.rds")
ps_mtg <- readRDS("data/mtg/ps_rle_nooutliers.rds")
ps_mtt <- readRDS("data/mtt/ps_rle_nooutliers.rds")
ps_metabol <- readRDS("data/metabol/ps_rle_nooutliers.rds")


getCoinertiaPs <- function(ps1, ps2){
  # Same samples
  sample_set <- intersect(sample_names(ps1), sample_names(ps2))
  print(paste("Number of samples overlapping: ", length(sample_set)))
  ps1 <- prune_samples(sample_set, ps1)
  ps2 <- prune_samples(sample_set, ps2)
  
  #Ordinate each independently
  ord1 <- ordinate(ps1, method = "PCoA", distance = "horn")
  plot_ordination(ps1, ord1, type = "samples", color = "phenotype")
  
  ord2 <- ordinate(ps2, method = "PCoA", distance = "horn")
  plot_ordination(ps2, ord2, type = "samples", color = "phenotype")
  
  
  #Co-inertia analysis
  tmp <- data.frame(t(ps2@otu_table))
  ps_cca <- phyloseq(otu_table(ps1@otu_table, taxa_are_rows = T), sample_data = sample_data(tmp))
  ord_cca <- ordinate(ps_cca, method = "CCA", distance = "euclidean")
  plot_ordination(ps_cca, ord_cca, type = "samples")
  
  
  sample_by_pc <- ord_cca$CA$u
  ps_red <- phyloseq(otu_table(sample_by_pc, taxa_are_rows = F), sample_data = sample_data(ps1@sam_data))
  return(ps_red)
}



testCorrelations <- function(data, mapping_numerical, component){
  corrs <- apply(mapping_numerical, 2, function(x){
    return(cor.test(x, as.numeric(data[, component]), method = "spearman"))
    
  })
  corrs
  return(corrs)
}

testCategorial <- function(data, mapping_char_use, component){
  tests <- apply(mapping_char_use, 2, function(x){
    option1 <- unique(x)[1]
    option2 <- unique(x)[2]
    keep <- !is.na(x)
    x <- x[keep]
    projection <- data[keep , component]

    return(wilcox.test(projection[x==option1], projection[x==option2]))
  } )
  return(tests)
}

plotSignificantVariables <- function(data, corrs, component1, component2, data_type = "num", title = "", sig_alpha = 0.07){
  pvals <- as.numeric(lapply(corrs, function(x) return(x$p.value)))
  pvals_adj <- p.adjust(pvals, 'fdr')
  if(data_type == "num"){
    rs <- lapply(corrs, function(x) return(x$estimate[[1]]))
  }
  if(data_type == "cat"){
    rs <- lapply(corrs, function(x) return(x$statistic[[1]]))
  }
  plots <- list()
  keep <-  pvals_adj < sig_alpha
  if(sum(keep) > 0){
    pvals_keep <- pvals_adj[keep]
    rs_keep <- rs[keep]
    if(data_type == "num"){
      statistics <- paste(paste("R2 = ", round(as.numeric(rs_keep), 3)), paste("pval_adj = ", round(as.numeric(pvals_keep), 3)))
    }
    if(data_type == "Cat"){
      statistics <- paste(paste("W = ", round(as.numeric(rs_keep), 3)), paste("pval_adj = ", round(as.numeric(pvals_keep), 3)))
    }
    
    label_df <- data.frame(statistics = statistics, variable = names(rs_keep))
    
    for( variable in names(rs_keep)){
      print(variable)
      data_use <- data[!is.na(data[, variable]), ]
      p1 <- ggplot(data = data_use, mapping = aes_string(x = component1, y = component2, label = "host_name", color = variable)) +
        geom_point( size = 3) +
        geom_text(aes_string(label = "host_name"), size = 3) + stat_ellipse(type = "norm") +
        ggtitle(paste0("Correlated with ", component1, " : ", variable))
      
      
      label <- label_df$statistics[label_df$variable == variable]
      print(label)
      p2 <- ggplot(data = data_use, mapping = aes_string(x = component1, y = variable)) + geom_point()+
        geom_smooth(method = "lm") + 
        ggtitle(variable)
      if(data_type == "num"){
        y_max <- max(data_use[, variable])
        y_min <- min(data_use[, variable])
        p2 <- p2 + geom_text(x = min(data_use[ , component1]) - 0.1, y = y_max + 0.4, label = label,  hjust = 0, size = 3) + ylim(y_min - 0.25, y_max + 0.55)
      }
      
      print(p1 + p2)
      plots[[variable]] = p1 + p2
    }
  }
  
  return(plots)
}

cleanMapping <- function(mapping){
  mapping_use <- mapping %>% select(-c("biospecimen_id", "biospecimen_name", "host_name", "timepoint", "familyID",
                                       "racial_group", "specific_food_allergy"))
  mapping_use$age <- as.numeric(mapping_use$age)
  mapping_character <- mapping_use[, sapply(mapping_use, class) == 'character' | sapply(mapping_use, class) == 'logical']
  mapping_character <- mapping_character %>% select(-c("stool_freq"))
  if("date" %in% colnames(mapping_character)){
    mapping_character <- mapping_character %>% select(-c("date"))
  }
  mapping_character$phenotype <- ifelse(mapping_character$phenotype == "A", T, F)
  mapping_character <- data.frame(apply(mapping_character, 2, as.logical))
  mapping_numerical <- mapping_use[, sapply(mapping_use, class) == 'numeric']
  mapping_character <- mapping_character %>% select(-c("diarrhea", "constipation", "env_tobacco", "other_GI_symptoms")) #not enough examples
  return(list(mapping_numerical, mapping_character))
}

testCoinertia <- function(ps_red){
  mapping <- data.frame(ps_red@sam_data)
  tmp <- cleanMapping(mapping)
  mapping_numerical = tmp[[1]]
  mapping_char = tmp[[2]]
  #split variables into numerical and categorical
  
  
  #make plotting data
  data <- data.frame(ps_red@otu_table)
  data <- cbind(data, mapping_numerical, mapping_char)
  data$host_name <- ps_red@sam_data$host_name
  
  #visualize
  ggplot(data = data, mapping = aes(x = CA1, y = CA2, label = host_name)) +
    geom_point(aes(color = phenotype), size = 3) + geom_text(aes(label = host_name), size = 3)
  
  #Calculate the numerical test scores
  corrs_num_ca1 <- testCorrelations(data, mapping_numerical, "CA1")
  corrs_num_ca2 <- testCorrelations(data, mapping_numerical, "CA2")
  corrs_num_ca3 <- testCorrelations(data, mapping_numerical, "CA3")
  
  #Visualize some of those interesting variables
  ca1 <- plotSignificantVariables(data, corrs = corrs_num_ca1, component1 = "CA1", component2 = "CA2", data_type = "num", sig_alpha = 0.07)
  ca2 <- plotSignificantVariables(data, corrs = corrs_num_ca2, component1 = "CA2", component2 = "CA3", data_type = "num", sig_alpha = 0.07)
  ca3 <- plotSignificantVariables(data, corrs = corrs_num_ca3, component1 = "CA3", component2 = "CA4", data_type = "num", sig_alpha = 0.07)
  
  
  #Calculate the categorical test scores
  tests_char_ca1 <- testCategorial(data, mapping_char, "CA1")
  tests_char_ca2 <- testCategorial(data, mapping_char, "CA2")
  tests_char_ca3 <- testCategorial(data, mapping_char, "CA3")
  
  #Visualize some of those interesting variables
  ca4 <- plotSignificantVariables(data, corrs = tests_char_ca1, component1 = "CA1", component2 = "CA2", data_type = "cat", sig_alpha = 0.07)
  ca5 <- plotSignificantVariables(data, corrs = tests_char_ca2, component1 = "CA2", component2 = "CA3", data_type = "cat", sig_alpha = 0.07)
  ca6 <- plotSignificantVariables(data, corrs = tests_char_ca3, component1 = "CA3", component2 = "CA4", data_type = "cat", sig_alpha = 0.07)
  
  return(list(numerical_ca1 = ca1, numerical_ca2 = ca2, numerical_ca3 = ca3, categorical_ca1 = ca4, categorical_ca2 = ca5, categorical_ca3 = ca6))
  
}

testOneOmic <- function(ps){
  mapping <- data.frame(ps@sam_data)
  tmp <- cleanMapping(mapping)
  mapping_numerical = tmp[[1]]
  mapping_char = tmp[[2]]
  
  ord <- ordinate(ps, method = "PCoA", distance = "bray")
  
  data <- data.frame(ord$vectors)
  data <- cbind(data, mapping_numerical, mapping_char)
  data$host_name <- ps@sam_data$host_name
  
  
  #Calculate the numerical test scores
  corrs_num_ca1 <- testCorrelations(data, mapping_numerical, "Axis.1")
  corrs_num_ca2 <- testCorrelations(data, mapping_numerical, "Axis.2")
  corrs_num_ca3 <- testCorrelations(data, mapping_numerical, "Axis.3")
  
  #Visualize some of those interesting variables
  ca1 <- plotSignificantVariables(data, corrs_num_ca1, "Axis.1", "Axis.2", data_type = "num", sig_alpha = 0.07)
  ca2 <- plotSignificantVariables(data, corrs_num_ca2, "Axis.2", "Axis.3", data_type = "num", sig_alpha = 0.07)
  ca3 <- plotSignificantVariables(data, corrs_num_ca3, "Axis.3", "Axis.4", data_type = "num", sig_alpha = 0.07)
  
  
  #Calculate the categorical test scores
  tests_char_ca1 <- testCategorial(data, mapping_char, "Axis.1")
  tests_char_ca2 <- testCategorial(data, mapping_char, "Axis.2")
  tests_char_ca3 <- testCategorial(data, mapping_char, "Axis.3")
  
  #Visualize some of those interesting variables
  ca4 <- plotSignificantVariables(data, tests_char_ca1, "Axis.1", "Axis.2", data_type = "cat", sig_alpha = 0.07)
  ca5 <- plotSignificantVariables(data, tests_char_ca2, "Axis.2", "Axis.3", data_type = "cat", sig_alpha = 0.07)
  ca6 <- plotSignificantVariables(data, tests_char_ca3, "Axis.3", "Axis.4", data_type = "cat", sig_alpha = 0.07)
  return(list(numerical_ca1 = ca1, numerical_ca2 = ca2, numerical_ca3 = ca3, categorical_ca1 = ca4, categorical_ca2 = ca5, categorical_ca3 = ca6))
}

printOutPlots <- function(p, variable_pair){
  for(i in seq(1, length(p))){
    if(length(p[[i]]) > 0){
      if(i %in% c(1,2,3)){
        filename <- paste0("results/coinertia/", variable_pair, "_coinertia_numerical_ca", i, ".png")
        ggsave(filename, wrap_plots(wrap_plots(p[[i]], ncol = 1)) + plot_annotation(title = variable_pair), height = 15, width = 10)
      }
      if(i %in% c(4,5,6)){
        filename <- paste0("results/coinertia/", variable_pair, "_coinertia_categorical_ca", i-3, ".png")
        ggsave(filename, wrap_plots(wrap_plots(p[[i]], ncol = 1)) + plot_annotation(title = variable_pair), height = 15, width = 10)
      }
      
    }
  }
}

variable_pair <- "mtg_mtt"
plot(x=1, y=1)
text(x=1, y=1, "MTG and MTT")
ps_red <- getCoinertiaPs(ps_mtg, ps_mtt)
p <- testCoinertia(ps_red)
printOutPlots(p, variable_pair)

#min_time_antibiotics
#meal_home_prep
#vegetable
#fat_oil_freq
#vegetable_freq
#fruit_freq

variable_pair <- "mtg_metabol"
plot(x=1, y=1)
text(x=1, y=1, "MTG and Metabol")
ps_red <- getCoinertiaPs(ps_mtg, ps_metabol)
p <- testCoinertia(ps_red)
printOutPlots(p, variable_pair)


#fat_oil_freq
#vegetable_freq
#fruit_freq
#dietary_restriction

variable_pair <- "mtt_metabol"
plot(x=1, y=1)
text(x=1, y=1, "MTT and Metabol")
ps_red <- getCoinertiaPs(ps_mtt, ps_metabol)
p <- testCoinertia(ps_red)
printOutPlots(p, variable_pair)

variable_pair <- "16s_metabol"
plot(x=1, y=1)
text(x=1, y=1, "16s and Metabol")
ps_red <- getCoinertiaPs(ps_16s, ps_metabol)
p <- testCoinertia(ps_red)
printOutPlots(p, variable_pair)

#sex
#dairy
#vegetable
#infant_diet
#age
#meats_and_seafood
#dairy_freq
#fat_oil_freq
#vegetable_freq
#fruit_freq
#fat_oil_freq

variable_pair <- "16s_mtg"
plot(x=1, y=1)
text(x=1, y=1, "16s and MTG")
ps_red <- getCoinertiaPs(ps_16s, ps_mtg)
p <- testCoinertia(ps_red)
printOutPlots(p, variable_pair)
#vegetable
#infant_diet
#age
#fat_oil_freq
#vegetable_freq
#fruit_freq
#whole_grain
#fruit
#restaurant
#sugary_food
#fruit_freq

variable_pair <- "16s_mtt"
plot(x=1, y=1)
text(x=1, y=1, "16s and MTT")
ps_red <- getCoinertiaPs(ps_16s, ps_mtt)
p <- testCoinertia(ps_red)
printOutPlots(p, variable_pair)
#seafood
#vegetable
#infant_diet
#fat_oil_freq
#vegetable_freq
#fruit_freq
#dog
#recently_ill



#MTG reflects fat/oil, vegetables, and fruit when combined with other omics
#By itself, MTG reflects fat/oil and age.
variable_pair <- "mtg"
plot(x=1, y=1)
text(x=1, y=1, "MTG")
p <- testOneOmic(ps_mtg)
printOutPlots(p, variable_pair)

#16s reflects infant_diet, on top of fat/oil, vegetables, and fruit
#By itself, 16s reflects...everything
variable_pair <- "16s"
plot(x=1, y=1)
text(x=1, y=1, "16s")
p <- testOneOmic(ps_16s)
printOutPlots(p, variable_pair)


#Can't attribute anything in particular to MTT
variable_pair <- "mtt"
plot(x=1, y=1)
text(x=1, y=1, "MTT")
p <- testOneOmic(ps_mtt)
printOutPlots(p, variable_pair)

#Metabol reflects
#issues with mapping file
map_1 <- data.frame(ps_mtg@sam_data)
map_2 <- data.frame(ps_mtt@sam_data)
map_1 <- map_1 %>% select(-c("date"))
map_2 <- map_2 %>% select(-c("date"))
map_tmp <- rbind(map_1, map_2)
keep <- sample_names(ps_metabol)[sample_names(ps_metabol) %in% rownames(map_tmp)]
ps_metabol <- prune_samples(keep, ps_metabol)
map_tmp <- map_tmp[sample_names(ps_metabol), ]
sample_data(ps_metabol) <- map_tmp
write.csv(map_tmp, "data/metabol/mapping_nooutliers2.csv")

variable_pair <- "metabol"
plot(x=1, y=1)
text(x=1, y=1, "Metabol")
p <- testOneOmic(ps_metabol)
printOutPlots(p, variable_pair)


#Conclusions: By looking at the omics together using a co-inertia analysis, we are able to triangulate what lifestyle
#variables are reflected by each omic. This is important for researchers to know, because these lifestyle variables
#manifest as confounding variables in gut microbiome studies, leading to potentially incorrect or biased conclusions if
#they are not controlled for. 

###############################################
#################  MARA  ######################
###############################################
#Use only autism samples to look at MARA
plotMara <- function(ps_red){
  ps_red <- prune_samples(ps_red@sam_data$phenotype == "A", ps_red)
  mapping <- data.frame(ps_red@sam_data)
  tmp <- cleanMapping(mapping)
  mapping_numerical = tmp[[1]]
  mapping_char = tmp[[2]]
  #split variables into numerical and categorical
  #make plotting data
  data <- data.frame(ps_red@otu_table)
  data <- cbind(data, mapping_numerical, mapping_char)
  data$host_name <- ps_red@sam_data$host_name
  
  #visualize
  ggplot(data = data, mapping = aes(x = CA1, y = CA2, label = host_name)) +
    geom_point(aes(color = phenotype), size = 3) + geom_text(aes(label = host_name), size = 3)
  
  #Calculate the numerical test scores
  corrs_num_ca1 <- testCorrelations(data, mapping_numerical, "CA1")
  corrs_num_ca2 <- testCorrelations(data, mapping_numerical, "CA2")
  corrs_num_ca3 <- testCorrelations(data, mapping_numerical, "CA3")
  
  ca1 <- plotSignificantVariables(data, corrs = corrs_num_ca1["MARA"], component1 = "CA1", component2 = "CA2", data_type = "num", sig_alpha = 1)
  ca2 <- plotSignificantVariables(data, corrs = corrs_num_ca2["MARA"], component1 = "CA2", component2 = "CA3", data_type = "num", sig_alpha = 1)
  ca3 <- plotSignificantVariables(data, corrs = corrs_num_ca3["MARA"], component1 = "CA3", component2 = "CA4", data_type = "num", sig_alpha = 1)
  
  return(ca1$MARA / ca2$MARA / ca3$MARA)
}

ps_16s@sam_data$MARA <- as.numeric(ps_16s@sam_data$MARA )
ps_mtg@sam_data$MARA <- as.numeric(ps_mtg@sam_data$MARA )
ps_mtt@sam_data$MARA <- as.numeric(ps_mtt@sam_data$MARA )
ps_metabol@sam_data$MARA <- as.numeric(ps_metabol@sam_data$MARA )

variable_pair <- "mtg_mtt"
plot(x=1, y=1)
text(x=1, y=1, "MTG and MTT")
ps_red <- getCoinertiaPs(ps_mtg, ps_mtt)
p <- plotMara(ps_red)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_mara", ".png")
ggsave(filename, p, height = 15, width = 10)


variable_pair <- "mtg_metabol"
plot(x=1, y=1)
ps_red <- getCoinertiaPs(ps_mtg, ps_metabol)
p <- plotMara(ps_red)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_mara", ".png")
ggsave(filename, p, height = 15, width = 10)

variable_pair <- "mtt_metabol"
plot(x=1, y=1)
text(x=1, y=1, "MTT and Metabol")
ps_red <- getCoinertiaPs(ps_mtt, ps_metabol)
p <- plotMara(ps_red)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_mara", ".png")
ggsave(filename, p, height = 15, width = 10)

variable_pair <- "16s_metabol"
plot(x=1, y=1)
text(x=1, y=1, "16s and Metabol")
ps_red <- getCoinertiaPs(ps_16s, ps_metabol)
p <- plotMara(ps_red)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_mara", ".png")
ggsave(filename, p, height = 15, width = 10)


variable_pair <- "16s_mtg"
plot(x=1, y=1)
text(x=1, y=1, "16s and MTG")
ps_red <- getCoinertiaPs(ps_16s, ps_mtg)
p <- plotMara(ps_red)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_mara", ".png")
ggsave(filename, p, height = 15, width = 10)


variable_pair <- "16s_mtt"
plot(x=1, y=1)
text(x=1, y=1, "16s and MTT")
ps_red <- getCoinertiaPs(ps_16s, ps_mtt)
p <- plotMara(ps_red)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_mara", ".png")
ggsave(filename, p, height = 15, width = 10)


variable_pair <- "mtg"
plot(x=1, y=1)
text(x=1, y=1, "MTG")
p <- plotMara(ps_mtg)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_mara", ".png")
ggsave(filename, p, height = 15, width = 10)

variable_pair <- "16s"
plot(x=1, y=1)
text(x=1, y=1, "16s")
p <- plotMara(ps_16s)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_mara", ".png")
ggsave(filename, p, height = 15, width = 10)

variable_pair <- "mtt"
plot(x=1, y=1)
text(x=1, y=1, "MTT")
p <- plotMara(ps_mtt)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_mara", ".png")
ggsave(filename, p, height = 15, width = 10)




map_1 <- data.frame(ps_mtg@sam_data)
map_2 <- data.frame(ps_mtt@sam_data)
map_1 <- map_1 %>% select(-c("date"))
map_2 <- map_2 %>% select(-c("date"))
map_tmp <- rbind(map_1, map_2)
keep <- sample_names(ps_metabol)[sample_names(ps_metabol) %in% rownames(map_tmp)]
ps_metabol <- prune_samples(keep, ps_metabol)
map_tmp <- map_tmp[sample_names(ps_metabol), ]
sample_data(ps_metabol) <- map_tmp

variable_pair <- "metabol"
plot(x=1, y=1)
text(x=1, y=1, "Metabol")
p <- plotMara(ps_metabol)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_mara", ".png")
ggsave(filename, p, height = 15, width = 10)




####################################
#############  Anxiety  ############
####################################
plotAnxiety <- function(ps_red){
  ps_anxiety <- readRDS("data/16s/ps_DeSeq_pass_min_postDD_min0.03.rds")
  anxiety <- ps_anxiety@sam_data$Recent.anxiety..caretaker.reported.
  anxiety_values <- rep(0, length(anxiety))
  anxiety_values[anxiety == "No elevated anxiety"] = 0
  anxiety_values[anxiety == "Somewhat elevated anxiety"] = 1
  anxiety_values[anxiety == "Elevated anxiety"] = 2
  names(anxiety_values) <- ps_anxiety@sam_data$Biospecimen.Name
  
  samples <- sample_names(ps_red)[sample_names(ps_red) %in% names(anxiety_values)]
  ps_red <- prune_samples(samples, ps_red)
  ps_red@sam_data$anxiety <- anxiety_values[samples]
  
  mapping <- data.frame(ps_red@sam_data)
  tmp <- cleanMapping(mapping)
  mapping_numerical = tmp[[1]]
  mapping_char = tmp[[2]]
  #split variables into numerical and categorical
  #make plotting data
  data <- data.frame(ps_red@otu_table)
  data <- cbind(data, mapping_numerical, mapping_char)
  data$host_name <- ps_red@sam_data$host_name
  
  #visualize
  ggplot(data = data, mapping = aes(x = CA1, y = CA2, label = host_name)) +
    geom_point(aes(color = phenotype), size = 3) + geom_text(aes(label = host_name), size = 3)
  
  #Calculate the numerical test scores
  corrs_num_ca1 <- testCorrelations(data, mapping_numerical, "CA1")
  corrs_num_ca2 <- testCorrelations(data, mapping_numerical, "CA2")
  corrs_num_ca3 <- testCorrelations(data, mapping_numerical, "CA3")
  
  ca1 <- plotSignificantVariables(data, corrs = corrs_num_ca1["anxiety"], component1 = "CA1", component2 = "CA2", data_type = "num", sig_alpha = 1)
  ca2 <- plotSignificantVariables(data, corrs = corrs_num_ca2["anxiety"], component1 = "CA2", component2 = "CA3", data_type = "num", sig_alpha = 1)
  ca3 <- plotSignificantVariables(data, corrs = corrs_num_ca3["anxiety"], component1 = "CA3", component2 = "CA4", data_type = "num", sig_alpha = 1)
  
  return(ca1$anxiety / ca2$anxiety / ca3$anxiety)
}


variable_pair <- "mtg_mtt"
plot(x=1, y=1)
text(x=1, y=1, "MTG and MTT")
ps_red <- getCoinertiaPs(ps_mtg, ps_mtt)
p <- plotAnxiety(ps_red)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_anxiety", ".png")
ggsave(filename, p, height = 15, width = 10)


variable_pair <- "mtg_metabol"
plot(x=1, y=1)
ps_red <- getCoinertiaPs(ps_mtg, ps_metabol)
p <- plotAnxiety(ps_red)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_anxiety", ".png")
ggsave(filename, p, height = 15, width = 10)

variable_pair <- "mtt_metabol"
plot(x=1, y=1)
text(x=1, y=1, "MTT and Metabol")
ps_red <- getCoinertiaPs(ps_mtt, ps_metabol)
p <- plotAnxiety(ps_red)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_anxiety", ".png")
ggsave(filename, p, height = 15, width = 10)

variable_pair <- "16s_metabol"
plot(x=1, y=1)
text(x=1, y=1, "16s and Metabol")
ps_red <- getCoinertiaPs(ps_16s, ps_metabol)
p <- plotAnxiety(ps_red)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_anxiety", ".png")
ggsave(filename, p, height = 15, width = 10)


variable_pair <- "16s_mtg"
plot(x=1, y=1)
text(x=1, y=1, "16s and MTG")
ps_red <- getCoinertiaPs(ps_16s, ps_mtg)
p <- plotAnxiety(ps_red)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_anxiety", ".png")
ggsave(filename, p, height = 15, width = 10)


variable_pair <- "16s_mtt"
plot(x=1, y=1)
text(x=1, y=1, "16s and MTT")
ps_red <- getCoinertiaPs(ps_16s, ps_mtt)
p <- plotAnxiety(ps_red)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_anxiety", ".png")
ggsave(filename, p, height = 15, width = 10)


variable_pair <- "mtg"
plot(x=1, y=1)
text(x=1, y=1, "MTG")
p <- plotAnxiety(ps_mtg)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_anxiety", ".png")
ggsave(filename, p, height = 15, width = 10)

variable_pair <- "16s"
plot(x=1, y=1)
text(x=1, y=1, "16s")
p <- plotAnxiety(ps_16s)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_anxiety", ".png")
ggsave(filename, p, height = 15, width = 10)

variable_pair <- "mtt"
plot(x=1, y=1)
text(x=1, y=1, "MTT")
p <- plotAnxiety(ps_mtt)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_anxiety", ".png")
ggsave(filename, p, height = 15, width = 10)




map_1 <- data.frame(ps_mtg@sam_data)
map_2 <- data.frame(ps_mtt@sam_data)
map_1 <- map_1 %>% select(-c("date"))
map_2 <- map_2 %>% select(-c("date"))
map_tmp <- rbind(map_1, map_2)
keep <- sample_names(ps_metabol)[sample_names(ps_metabol) %in% rownames(map_tmp)]
ps_metabol <- prune_samples(keep, ps_metabol)
map_tmp <- map_tmp[sample_names(ps_metabol), ]
sample_data(ps_metabol) <- map_tmp

variable_pair <- "metabol"
plot(x=1, y=1)
text(x=1, y=1, "Metabol")
p <- plotAnxiety(ps_metabol)
p <- p + plot_annotation(variable_pair)
filename <- paste0("results/coinertia/", variable_pair, "_anxiety", ".png")
ggsave(filename, p, height = 15, width = 10)





#phenotype nor MARA score was directly related to any combination of these omics, however, dietary preferences were.

trainClassifier <- function(ps){
  set.seed(2)
  folds <- groupKFold(group = ps@sam_data$familyID, k = 7)
  control <- trainControl(method='repeatedcv', 
                          number=4, 
                          repeats=3,
                          allowParallel = F)
  
  tunegrid <- expand.grid(.mtry=c(1:20))
  
  #ps@sam_data$phenotype <- ifelse(ps@sam_data$phenotype == "A", 1, 0) 
  index_train = folds[[1]]
  x_train <- ps@otu_table[index_train, ]
  x_test <- ps@otu_table[-index_train, ]
  metadata <- data.frame(ps@sam_data) %>% select(-c( biospecimen_id, familyID,biospecimen_name,
                                                     host_name, timepoint, age, sex, racial_group, MARA,
                                                     env_tobacco, gluten_allergy, diarrhea, constipation))
  if("date" %in% colnames(metadata)){
    metadata <- metadata %>% select(-c(date))
  }
  
  data_train <- metadata[index_train, ]
  data_test <- metadata[-index_train, ]
  model_metadata <- train(phenotype ~ ., 
                 data = data_train, 
                 method='rf', 
                 metric='Accuracy', 
                 tuneGrid=tunegrid, 
                 trControl=control,
                 na.action = na.exclude)
  best_accuracy_meta <- max(model_metadata$results$Accuracy)
  
  
  data_train <- data.frame(cbind(x_train, metadata[index_train, ]))
  model_all <- train(phenotype ~ ., 
                          data = data_train, 
                          method='rf', 
                          metric='Accuracy', 
                          tuneGrid=tunegrid, 
                          trControl=control,
                          na.action = na.exclude)
  best_accuracy_all <- max(model_all$results$Accuracy)
  
  print(paste("Best accuracy meta: ", best_accuracy_meta))
  print(paste("Best accuracy all: ", best_accuracy_all))
}

#Classifiers
trainClassifier(getCoinertiaPs(ps_16s, ps_mtg))
trainClassifier(getCoinertiaPs(ps_16s, ps_mtt))
trainClassifier(getCoinertiaPs(ps_16s, ps_metabol))

trainClassifier(getCoinertiaPs(ps_mtg, ps_mtt))
trainClassifier(getCoinertiaPs(ps_mtg, ps_metabol))
trainClassifier(getCoinertiaPs(ps_mtt, ps_metabol))



#split the phyloseq objects into training and testing to make our lives easier later on

ps_train <- phyloseq(otu_table(ps@otu_table[index_train, ], taxa_are_rows = F), ps@sam_data[index_train, ])
ps_test <- phyloseq(otu_table(ps@otu_table[-index_train, ], taxa_are_rows = F), ps@sam_data[-index_train, ])

set.seed(1)
data_train = data.frame(x_train)
data_train$hive_role = ps_train@sam_data$bee_type
control <- trainControl(method='repeatedcv', 
                        number=3, 
                        repeats=3,
                        allowParallel = F)

tunegrid <- expand.grid(.mtry=c(1:20)) #mtry is the depth of each decision tree. We'll be trying out models where each tree is 3 to 20 splits deep
rf <- train(hive_role ~., 
            data= data_train, 
            method='rf', 
            metric='Accuracy', 
            tuneGrid=tunegrid, 
            trControl=control)
print(rf)












