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

samples <- intersect(sample_names(ps_mtg), sample_names(ps_mtt))
mtg <- ps_mtg@otu_table[ , samples]
mtt <- ps_mtt@otu_table[, samples]

genes <- intersect(rownames(mtg), rownames(mtt))
mtg <- mtg[genes, ]
mtt <- mtt[genes, ]

pvals <- c()
for(i in seq(1, nrow(mtg))){
  pval <- cor.test(as.numeric(mtg[i, ]), as.numeric(mtt[i, ]), method = "spearman")$p.value
  pvals <- c(pvals, pval)
  print(i)
}

keep <- p.adjust(pvals, "fdr") < .05


#Differential abundance
mtg_filt <- mtg[keep, ]
mtt_filt <- mtt[keep, ]
metadata <- ps_mtg@sam_data[samples, ]

pvals <- c()
for(i in seq(1, ncol(mtg))){
  pval <- wilcox.test(as.numeric(mtg[ ,metadata$phenotype == "A"]), as.numeric(mtg[ , metadata$phenotype == "N"]))$p.value
  pvals <- c(pvals, pval)
  print(i)
}