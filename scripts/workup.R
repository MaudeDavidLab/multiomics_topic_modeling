library(vegan)
library(phyloseq)
library(dplyr)
library(lme4)
library(caret)
library(fossil)
library(DESeq2)
library(ape)
library(reshape2)
library(randomForest)
library(caret)
library(dplyr)
library(KEGGREST)
library(pROC)
library(glmnet)



dir = "C:/Users/kitikomp/Documents/Lab/M3/"
setwd(dir)

############################################################################################################################
############################################################################################################################
#######################################################   Setup functions  #################################################
############################################################################################################################
############################################################################################################################
#load functions.R

############################################################################################################################
############################################################################################################################
#######################################################   16S  #############################################################
############################################################################################################################
############################################################################################################################

ps <- readRDS("ps_not_norm_comp_pass_min_postDD_age_filtered.rds")
hist(sample_sums(ps))

#Rarefy only the incredibly abundant samples. There are just way too many counts
abundant_samples <- sample_names(ps)[sample_sums(ps) > 600000]
rarefied_samples <- rrarefy(t(ps@otu_table[ , abundant_samples]), 400000)
ps@otu_table[ , abundant_samples] <- rarefied_samples
hist(sample_sums(ps))

months <- unlist(lapply(ps@sam_data$Biospecimen.Date.Collected, function(x) return(substring(x, 6, 7))))
ps@sam_data$month <- as.factor(months)

#age buckets
ps@sam_data$age_buckets <- cut(ps@sam_data$Age..months., breaks = c(0, 24, 24*2, 24*4, 24*6))

### write piphillan
seqtab <- ps@otu_table
seqs <- as.character(ps@tax_table[rownames(seqtab), 2])
headers <- paste(">otu", seq(1, nrow(seqtab)), sep = "")
sample_names <- paste("Sample1", seq(1, ncol(seqtab)), sep = "")
colnames(seqtab) <- sample_names
rownames(seqtab) <- gsub(">", "", headers)
write.csv(seqtab, "seqtab_piphillan.csv", quote = F)

fasta <- paste(headers, seqs, collapse = "\n", sep = "\n")
write(fasta, "repseqs.fasta")


#Make piphillan object
piph_table <- read.delim("C:/Users/ctata/Documents/Lab/M3/data/piphillan/results/ko_abund_table_unnorm.txt", row.names=1)
colnames(piph_table) <- sample_names(ps)
ps_piph <- phyloseq(otu_table(piph_table, taxa_are_rows = T), sample_data(ps))
ps_piph <- addNumericMetadata(ps_piph)
saveRDS(ps_piph, "ps_piph.rds")


####################################
#### Taxonomic Class ###############
####################################
seqs <- as.character(ps@tax_table[,2])
taxa <- assignSpecies(seqs, "../public_data/rdp_species_assignment_16.fa.gz")

#####################################
##### Metadata analysis #############
#####################################

ps_tmp <- addNumericMetadata(ps_rle)


#########################################
### Select normalization based        ###
### on individual timepoint proximity ###
#########################################
#Check this function:
# m <- matrix(c(x, y), nrow = 2, byrow = T)
# vegdist(m, "bray) == bray_curtis_dist(x, y)
bray_curtis_dist <- function(x, y){
  return(sum(abs(x-y)) / sum(x + y))
}

getMeanOutDists <- function(host_name, ps, otu_table){
  out_samples <- otu_table[ , ps@sam_data$Host.Name != host_name]
  in_samples <- otu_table[ , ps@sam_data$Host.Name == host_name]

  dists2 <- apply(in_samples, 2, function(in_sample){
    dists <- apply(out_samples, 2, function(out_sample){
      bray.curtis(in_sample, out_sample)
    })
    mean(dists)
  })
  mean_out_dists <- mean(dists2)
  return(mean_out_dists)
}

getMeanInDists <- function(ps, host_name, otu_table){
  samples <- otu_table[ , ps@sam_data$Host.Name == host_name]
  mean_dist <- mean(vegdist(t(samples), method = 'bray'))
  return(mean_dist)
}

#Calculate the average distance between replicates
getInOutRatio <- function(otu_table){
  mean_dists <- pbsapply(unique(ps@sam_data$Host.Name), function(n){
    return(getMeanInDists(ps, n, otu_table))
  }) # one element per subject

  mean_out_dists <- pbsapply(unique(ps@sam_data$Host.Name), function(host_name) return(getMeanOutDists(host_name, ps, otu_table))) # one element per subject

  return(mean_dists / mean_out_dists) #we want this ratio to be low
}



deSeqNorm <- function(ps){
  ps_dds <- phyloseq_to_deseq2(ps, ~ phenotype)
  ps_dds <- estimateSizeFactors(ps_dds, type = "poscounts")
  ps_dds <- estimateDispersions(ps_dds)
  abund <- getVarianceStabilizedData(ps_dds)
  abund <- abund + abs(min(abund)) #don't allow deseq to return negative counts
  ps_deSeq <- phyloseq(otu_table(t(abund), taxa_are_rows = F), sample_data(ps))
  return(ps_deSeq)
}

saveRDS(ps_deSeq, "ps_deseq.rds")
cssNorm <- function(ps){
  p_metag<-phyloseq_to_metagenomeSeq(ps)
  norm_df <- MRcounts(p_metag, norm=T)
  otu_table(ps) <- otu_table(norm_df, taxa_are_rows = T)
  return(ps)
}
ps_css <- cssNorm(ps_tmp)
ps_deseq <- deSeqNorm(ps)
ps_asinh <- phyloseq(asinh(ps@otu_table), sample_data(ps@sam_data))
ps_rle <- phyloseq(otu_table(normalize_rle(ps@otu_table, addOne = T)), sample_data(ps@sam_data))

rle <- getInOutRatio(ps_rle@otu_table)
asinh <- getInOutRatio(ps_asinh@otu_table)
#Both these options seem acceptable
#css <- getInOutRatio(ps_css@otu_table)
#deseq <- getInOutRatio(ps_deseq@otu_table)


##########################################
######### Variance across replicates  ####
##########################################

getMeanReplicateDistance <- function(names, ps){
  mean_dists <- sapply(names, function(sample_name){
    sample_ids <- rownames(ps@sam_data)[ps@sam_data$Host.Name == sample_name]
    dists <- vegdist(ps@otu_table[, sample_ids], method = "bray")
    return(mean(dists))
  })
  return(mean_dists)
}


aut_names <- unique(ps_rle@sam_data$Host.Name[ps_rle@sam_data$phenotype == "A"])
mean_aut_dists <- getMeanReplicateDistance(aut_names, ps_rle)

nt_names <- unique(ps_rle@sam_data$Host.Name[ps_rle@sam_data$phenotype == "N"])
mean_nt_dists <- getMeanReplicateDistance(nt_names, ps_rle)

######################
#ordinate samples##
######################
ord <- ordinate(ps_rle, "NMDS", "bray")
plot_ordination(ps_rle, ord, type = "sample", color = "Family.group.ID..Biospecimen.", label = "Biospecimen.Name" )

######################
#### PRUNE OUTLIERS ###
######################
#We see that we have some massive outliers, both with an without normalization. The read counts are not disparate between
#these 3 and the rest of the population. I'm not going to look more into why these are outliers, instead I'll just remove them
#Possible reasons: illness, antibiotics, age
outliers <- c("M3_043_P3_N_V1", "M3_186_P1_N_V1", "M3_097_P2_N_V1")
colSums(ps@otu_table[, ps@sam_data$Biospecimen.Name %in% outliers])
ps <- prune_samples(!(ps@sam_data$Biospecimen.Name %in% outliers), ps)

########################
### WRITE FASTA ########
########################
ps@tax_table
headers <- paste(">", rownames(ps@tax_table), sep = "")
seqs <- as.character(ps@tax_table[, 2])
write(paste(headers, seqs, sep = "\n", collapse = "\n"), file = "repseqs.fasta")

#############################
#ordinate of samples again ##
#############################

ord <- ordinate(ps, "NMDS", "bray")
ps@sam_data$Family.group.ID..Biospecimen. <- as.factor(ps@sam_data$Family.group.ID..Biospecimen.)
plot_ordination(ps, ord, type = "sample", color = "Family.group.ID..Biospecimen." ) + facet_wrap(~age_buckets)
#We're trying to see that roughly timepoints and individuals are in the same area
plot_ordination(ps, ord, type = "sample", color = "phenotype" ) + facet_wrap(~age_buckets)


#Theres a correlation between age and phenotype. Does the microbiome look affected by age?
form <- as.formula(paste("dists" , paste(colnames(ps@sam_data), collapse = "+"), sep = " ~ "))
dists <- phyloseq::distance(ps, method = "bray", type = "samples")
adonis2(dists ~ age_buckets, data = data.frame(ps@sam_data))
#The microbiome is affected by age, as is the phenotype. This is something we should control for




#############################
### PCoA on Metadata ########
#############################
df <- select(ps_tmp@sam_data, -c("familyID", "host_name"))
meta_dist <- dist(ps_tmp@sam_data, method = "euclidean")
ord <- ordinate(ps_tmp, "PCoA", distance = meta_dist)
plot_ordination(ps_tmp, ord, type = "samples", color = "phenotype", label = "host_name")

#Seems like there is a greater than random split here - kids with autism have some lifestyle characteristics
#different than their neurotypical counterparts

pcoa_obj <- pcoa( meta_dist)
rownames(pcoa_obj$vectors) <- rep("", nrow(df))
biplot(pcoa_obj, df, expand = 10)


##############################
##### PERMANOVA  #############
##############################


perm_sig <- runPermanova(ps_tmp)
sig_vars <- rownames(perm_sig)[perm_sig$`Pr(>F)` < .05 & perm_sig$disp_pval > .05]
saveRDS(sig_vars, "sig_vars.rds")
#Phenotype seems like it's supposed to be significantly associated with microbiome structure. There is hope.
#cat, sex, gluten allergy, whole_grain, dairy, fruit, lactose intolerance, multivitamin, dietary restriction,
#starchy food, meats and seafood, dairy freq, fat and oil freq




#######################
####   CCA   ##########
#######################
sig_vars <- readRDS("sig_vars.rds")
ps_cca <- ps_tmp
for(var in sig_vars){
  ps_cca <- prune_samples(!is.na(ps_cca@sam_data[ , var][[1]]), ps_cca)
}

form <- formula(paste("ps_cca" , paste(sig_vars, collapse = "+"), sep = " ~ "))
ord_cca <- ordinate(ps_cca, formula = form, method = "CCA")
p0 <- plot_ordination(ps_cca, ord_cca, type = "split", label = "host_name", color = "phenotype")

# Now add the environmental variables as arrows
arrowmat = vegan::scores(ord_cca, disp = "bp")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = CCA1, yend = CCA2, x = 0, y = 0, shape = NULL, color = NULL,
                label = labels)
extend <- 7
label_map = aes(x = extend * CCA1, y = extend * CCA2, shape = NULL, color = NULL,
                label = labels)
# Make a new graphic
arrowhead = arrow(length = unit(0.55, "npc"))
p1 = p0 + geom_segment(arrow_map, size = 1, data = arrowdf, color = "gray",
                       arrow = arrowhead) + geom_text(label_map, size = 3, data = arrowdf)
p1




##########################################
#########  Metadata between phenotypes ###
##########################################

df <- select(ps_tmp@sam_data, -c("familyID", "host_name"))
df$phenotype <- as.factor(df$phenotype) #0 is autism

plots <- c()
for(var in sig_vars){
  print(var)
  p <- ggplot(data.frame(df), aes_string(x = "phenotype" , y = var)) + geom_boxplot() + geom_jitter()
  plot(p)
}
#Conclusion: we will need to add these important factors to any model because they are correlate with phenotype
# and also microbial structure. A microbe or KO that is different between NT and AUT maybe be due to differences
# in dairy consumption, and so will likely not be a good candidate for mouse feeding. Likewise, a baseline model
# will include these metadata variables. We'd like our taxa or KOs to give us MORE information on top of this.


#####################
### Embed  ##########
#####################

embed <- function(seqtab, best_hits, repseqs_file, out_dir, id_thresh = 99, eval_thresh = 1*10^(-29), saveData = T){
  #seqtab has samples for rows and ASVs for columns

  #Filter best_hits table to only include hits that pass the e-value threshold
  best_hits <- best_hits[best_hits$evalue < 1*10^(-29), ]
  best_hits <- best_hits[best_hits$piden >= id_thresh, ]

  #Match up the query id in best_hits to the full ASV in colnames(seqtab) using the repseqs.fasta file
  fasta <- read.table(repseqs_file)
  headers <- gsub(">", "", fasta[seq(1, nrow(fasta), by = 2), 1])
  query_seqs <- as.character(fasta[seq(2, nrow(fasta), by = 2), 1])
  fasta_df <- data.frame(query_seqs, row.names = headers)
  best_hits$full_query_seq <- as.character(fasta_df[rownames(best_hits), 1])
  best_hits <- best_hits[as.character(best_hits$full_query_seq) %in% colnames(seqtab), ]
  print(dim(best_hits))

  #Drop any ASVs from the table that don't have near enough hits in the transformation matrix
  seqtab_hits <- seqtab[ , as.character(best_hits$full_query_seq)]
  print(dim(seqtab_hits))
  #Assign the id of each ASV's nearest hit in the embedding transformation table.
  seqtab_out <- seqtab_hits
  colnames(seqtab_out) <- best_hits$query_seq
  colnames(seqtab_hits) <- best_hits$hit_id


  qual_vec_file = "../quality_vectors_final/data/embed/embed_.07_100dim.txt" #need to figure out how to include data object in package
  qual_vecs <- read.table(qual_vec_file)
  qual_vecs_ord <- qual_vecs[colnames(seqtab_hits), ]

  embedded <- as.matrix(seqtab_hits) %*% as.matrix(qual_vecs_ord)
  if(saveData){
    embedded_out_file = paste(out_dir, "embedded_.07_100dim_",id_thresh, sep = "")
    saveRDS(embedded, paste(embedded_out_file, ".rds", sep = ""))
    write.table(embedded, paste(embedded_out_file, ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)

    seqtab_file = paste(out_dir, "seqtab_.07_100dim_", id_thresh, sep = "")
    write.table(seqtab_out, paste(seqtab_file, ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
  }

  return(embedded)
}

#change rownames to full ASVs
seqtab <- ps_tmp@otu_table
rownames(seqtab) <- as.character(ps@tax_table[, 2])
seqtab <- t(seqtab)

best_hits <- read.delim("embed/best_hits.tsv", header=FALSE, row.names = 1)
colnames(best_hits) <- c("hit_id", "query_seq", "hit_seq", "evalue", "bitscore", "length", "piden")

repseqs_file = "repseqs.fasta"
out_dir = "embed/"

embedded <- embed(seqtab, best_hits, repseqs_file, out_dir)


################################################
#### Linear Mixed Models  ######################
#### Controlling for Lifestyle Variables #######
################################################

getResiduals <- function(y, df){
  df$y <- y
  form <- formula(paste(paste("y" , paste(sig_vars[sig_vars!="phenotype"], collapse = "+"), sep = " ~ "), " + (1|familyID)"))
  #form <- formula(y ~  dairy + fruit + gluten_allergy + whole_grain + (1|familyID) + phenotype )
  lme_fit <- lmer(form, data = df)
  res <- summary(lme_fit)
  resids <- residuals(lme_fit)
  return(resids)
}

getSig <- function(y, df){
  df$y <- y
  form <- formula("y ~ phenotype")
  fit <- lm(form, data = df)
  pval <- summary(fit)$coefficients["phenotype" , 4]
  return(pval)
}

sig_vars <- readRDS("sig_vars.rds")
familyID <- ps_tmp@sam_data$familyID
df <- data.frame(ps_tmp@sam_data[ , sig_vars])
df$familyID <- as.factor(familyID)
#df <- data.frame(apply(df, 2, function(x) as.factor(x)))
df <- na.omit(df)
ps_prune <- ps_tmp
ps_prune <- prune_samples(sample_names(ps_prune) %in% rownames(df), ps_prune)

resids <- list()
keep <- apply(ps_prune@otu_table, 1, function(x) return(var(x) > 1))
ps_prune <- prune_taxa(keep, ps_prune)
for(i in 1:ntaxa(ps_prune)){
  print(i)
  y <- as.numeric(ps_prune@otu_table[i, rownames(df)])
  tryCatch({
    res <- getResiduals(y, df)
    #pvals[i] <- pval
    resids[[i]] <- res
  }, error=function(e, i){
    print(e)
    })
}

pvals <- rep(NA, length(resids))
for(i in 1:length(resids)){
  y <- resids[[i]]
  pval <- getSig(y, df)
  pvals[i] <- pval
}


p_adj <- p.adjust(pvals)
sum(p_adj < .05, na.rm = T)
ps_tmp@tax_table[rownames(ps_tmp@otu_table)[which(p_adj < .05)]]
#  "ASV_1008" "ASV_1016" "ASV_1017" "ASV_1020" "ASV_1297" "ASV_606"  "ASV_982"  "ASV_987"



#############################################
#### Mixed models with embeddings  ##########
#############################################
familyID <- ps_tmp@sam_data$familyID
df$familyID <- familyID
df <- na.omit(df)
ps_embed <- phyloseq(otu_table(t(embedded), taxa_are_rows = T), sample_data(df))

resids <- list()
for(i in 1:ntaxa(ps_embed)){
  print(i)
  y <- as.numeric(ps_embed@otu_table[i, rownames(df)])
  tryCatch({
    res <- getResiduals(y, data.frame(df))
    #pvals[i] <- pval
    resids[[i]] <- res
  }, error=function(e, i){
    print(e)
  })
}

pvals <- rep(NA, length(resids))
for(i in 1:length(resids)){
  y <- resids[[i]]
  pval <- getSig(y, data.frame(df))
  pvals[i] <- pval
}
hist(pvals)


########### Sequential control for metadata
y <- as.numeric(ps_tmp@otu_table["ASV_1016", ])
form <- formula(paste(paste("y" , paste(sig_vars[sig_vars != "phenotype"], collapse = "+"), sep = " ~ "), " + (1|familyID)"))
df$y <- y
fit <- lmer(form, data = df)

y <- residuals(fit)
df <- df[names(y), ]
df$y <- y
fit <- lmer(formula("y ~ phenotype + (1|familyID)"), data = df)
print(summary(fit))

###################################
#########   Use Deseq  ############
###################################







#####################
### Random Forest ###
#####################




input <- cbind(embedded[rownames(df), ], data.frame(df) %>% select( -c(familyID, phenotype)))
metadata_model = trainRf(data.frame(df) %>% select( -c(familyID, phenotype)), map = data.frame(df))
embed_model = trainRf(input, map = data.frame(df))

input <- cbind(t(ps_tmp@otu_table[ , rownames(df)]), data.frame(df) %>% select( -c(familyID, phenotype)))
asv_model = trainRf(input, map = data.frame(df))

######################
### Kmer models ######
######################
kmer_7 <- readRDS("taxa_kmer_7.rds")
kmer_7 <- kmer_7[colnames(seqtab), ]
asv_kmer_7 <- as.matrix(seqtab) %*% as.matrix(kmer_7)
kmer_7_model <- trainRf(asv_kmer_7)

############################################################################################################################
############################################################################################################################
#######################################################   Metagenomes  #####################################################
############################################################################################################################
############################################################################################################################

ps_kegg <- readRDS("data/mtg/ps_rle_nooutliers.rds")
hist(as.numeric(sample_sums(ps_kegg)))

#Our options are drop the samples with high counts (7) or use percents
# a full 66% of the variance in the dataset is along the axis of the outliers after normalization. I think this means we need to
#take out the outliers, because they're just way too different and will obscure everything


#####################################
##### Get pathway names #############
#####################################


path_names <- sapply(taxa_names(ps_kegg), function(x) return(getPathwayName(x)))

#Keep only prokaryote pathways
pro1 <- path_names[grepl("ko00", names(path_names))]
pro2 <- path_names[grepl("ko01", names(path_names))]
pro3 <- path_names[grepl("ko02", names(path_names))]
prokar_paths <- c(pro1, pro2, pro3)
keep <- taxa_names(ps_kegg) %in% names(prokar_paths)
ps_kegg <- prune_taxa(keep, ps_kegg)


#We should also look at the human associated pathways, but this might be from contamination

ord <- ordinate(ps_kegg, "PCoA", "bray")
plot_ordination(ps_kegg, ord, type = "sample", color = "phenotype", label = "Biospecimen.Name" )
outliers <- ord$vectors[,1] > 0.1
ps_kegg <- prune_samples(!outliers, ps_kegg)
ord <- ordinate(ps_kegg, "PCoA", "bray")
plot_ordination(ps_kegg, ord, type = "sample", color = "phenotype", label = "Biospecimen.Name" )

###########################
## Assume pre-normalized ##
###########################
#otu_table(ps_kegg) <- normalize_rle(ps_kegg@otu_table)
#ord <- ordinate(ps_kegg, "PCoA", "bray")
#plot_ordination(ps_kegg, ord, type = "sample", color = "phenotype", label = "Biospecimen.Name" )


###########################
##### Numeric metadata ####
###########################
months <- as.numeric(unlist(lapply(ps_kegg@sam_data$Biospecimen.Date.Collected, function(x) return(substring(x, 6, 7)))))
months <- abs(months - 12) # we want december to be considered 0
sample_data(ps_kegg)$months <- months
ps_meta <- addNumericMetadata(ps_kegg)
sample_data(ps_meta)$months <- ps_kegg@sam_data[sample_names(ps_meta), "months" ]
ps_kegg <- ps_meta

saveRDS(ps_kegg, "ps_kegg_propath_metadata.rds")

ps_kegg <- readRDS("ps_kegg_propath_metadata.rds")
###############################
######  Permutation test  #####
###############################
#To create a null distribution, pick random pairs so long as they are not siblings


############### Binomial model  ############
null_dist <- getNullDist(ps_kegg, iter = 100000)
actual_diffs <- getActualDiffs(ps_kegg)
sig_kos_pos <- findSig_Binomial(null_dist, actual_diffs, pos = T)
sig_kos_neg <- findSig_Binomial(null_dist, actual_diffs, pos = F)

#ko_ids                                  path_names p_adj.p_adj...0.05.
#ko00190 ko00190                   Oxidative phosphorylation          0.01139605
#ko00310 ko00310                          Lysine degradation          0.01139605
#ko00520 ko00520 Amino sugar and nucleotide sugar metabolism          0.03583506
#ko00950 ko00950          Isoquinoline alkaloid biosynthesis         0.06160240

#ko00626 ko00626  Naphthalene degradation         0.07939011
#ko01524 ko01524 Platinum drug resistance         0.09333663

#we figure the neurotypical sibling is the best example of what your microbiome should look like given your similarity in metadata

######################################
###### Bayesian approach  ############
######################################
#somehow we want to weight the vote that the microbiome composition of each sample gets based on its proximity in the
#metadata to the query sample


#Prior will be the probability of autism given metadata
df <- ps_kegg@sam_data
df <- df  %>% select(-c("host_name", "phenotype", "familyID"))
df_aut <- ps_kegg@sam_data[ps_kegg@sam_data$phenotype == 1]
df_aut <- df_aut %>% select(-c("host_name"))
avg_aut <- rowMeans(df_aut, na.rm = T)

df_nt <- ps_kegg@sam_data[ps_kegg@sam_data$phenotype == 0]
df_nt <- df_nt %>% select(-c("host_name"))
avg_nt <- rowMeans(df_nt, na.rm = T)


nb <- naive_bayes(df, as.factor(ps_kegg@sam_data$phenotype))
prior_aut <- predict(nb, df, "prob")[,2]
sum(predict(nb, df, "class") == ps_kegg@sam_data$phenotype) / nsamples(ps_kegg)

nb_asv <- naive_bayes(t(ps_kegg@otu_table), as.factor(ps_kegg@sam_data$phenotype) )
sum(predict(nb_asv, t(ps_kegg@otu_table)) == ps_kegg@sam_data$phenotype) / nsamples(ps_kegg)

input <- cbind(df, t(ps_kegg@otu_table))
nb_combo <- naive_bayes(input, as.factor(ps_kegg@sam_data$phenotype))
sum(predict(nb_combo, input) == ps_kegg@sam_data$phenotype) / nsamples(ps_kegg)

#########################################################
######## Linear modeling - weight by metadata  ##########
#########################################################


#ko00254 ko00254   53                        Aflatoxin biosynthesis
#ko00310 ko00310   48                            Lysine degradation
#ko00514 ko00514   12          Other types of O-glycan biosynthesis
#ko00860 ko00860  151          Porphyrin and chlorophyll metabolism


#ko00073         ko00073    1          Cutin, suberine and wax biosynthesis
#ko00944         ko00944    1             Flavone and flavonol biosynthesis
#ko00965         ko00965    7                         Betalain biosynthesis
#ko01524         ko01524   47                      Platinum drug resistance


#Area under the curve: 0.7755 2.5 * sd(model$lambda) + mean(model$lambda)
#Area under the curve: 0.7503 2.3 * sd(model$lambda) + mean(model$lambda)
#Area under the curve: 0.7046  2 * sd(model$lambda) + mean(model$lambda)



############################
##  Permanova ##############
############################
perm_sig <- runPermanova(ps_tmp)
sig_vars <- rownames(perm_sig)[perm_sig$`Pr(>F)` < .05 & perm_sig$disp_pval > .05]

#############################
### Metagenome classifier ###
#############################
input <- cbind(t(ps_tmp@otu_table), ps_tmp@sam_data %>% select( -c(familyID, phenotype)))

input <- ps_tmp@sam_data %>% select( -c(familyID, phenotype, host_name, multivitamin, probiotic, dietary_restriction,
                                        dietary_supplement, vitamin_B, vitamin_D, sex))
input <- na.omit(input)
map <- ps_tmp@sam_data[rownames(input)]
metadata_model = trainRf(input, map)
imp <- randomForest(input, as.factor(map$phenotype), mtry = 6, importance = T)
varImpPlot(imp)


input_combo <- cbind(input, t(ps_tmp@otu_table)[rownames(input), ])
combo_model <- trainRf(input_combo, map)

asv_model <- trainRf(t(ps_tmp@otu_table)[rownames(input), ], map)

# What key metabolites are missing if you aren't eating fruit, or dairy, or vegetables.
# Can an intervention simply be: high fiber diet?




##########################################
###### Differential pathway analysis  ####
##########################################

hist(as.numeric(ps_tmp@otu_table[5,]))
aut_samples <- ps_tmp@sam_data$phenotype == 1
nt_samples <- ps_tmp@sam_data$phenotype == 0



apply(ps_tmp@otu_table, 1, function(x){
  return(t.test(as.numeric( x[aut_samples]), as.numeric(x[nt_samples]))$p.val)
})

apply(ps_tmp@otu_table, 1, function(x){
  return(wilcox.test(as.numeric( x[aut_samples]), as.numeric(x[nt_samples]))$p.val)
})


hist(apply(ps_tmp@otu_table, 1, var), breaks = 20)
filt <- ps_tmp@otu_table[apply(ps_tmp@otu_table, 1, var) > .00002, ]
apply(filt, 1, function(x){
  return(t.test(as.numeric( x[aut_samples]), as.numeric(x[nt_samples]))$p.val)
})

linearModeling <- function(ps){
  input <- na.omit(data.frame(ps@sam_data))
  pathway_meta_ass <- data.frame()
  for(path in taxa_names(ps)){
    y <- as.numeric(ps@otu_table[path, rownames(input)])
    input$y <- y
    path_name <- getPathwayName(path)
    m <- lm(y ~ fruit + vegetable_freq + gluten_allergy + stool_freq + dairy_freq + recently_ill + lactose_intolerance + phenotype, data = data.frame(input))
    tmp <- summary(m)$coefficients
    keep <- data.frame(tmp[tmp[,4]<.1 & rownames(tmp) != "(Intercept)", ,drop = FALSE])
    print(path_name)
    if(nrow(keep) > 0){
      path_name_rep <- rep(path_name, nrow(keep))
      path_id <- rep(path, nrow(keep))
      keep$var <- rownames(keep)
      keep <- cbind(keep, path_name_rep, path_id)
      pathway_meta_ass <- rbind(pathway_meta_ass, keep)
    }
  }
  return(pathway_meta_ass)
}
pathway_meta_ass <- linearModeling(ps_kegg)
path_ids <- as.character(pathway_meta_ass[pathway_meta_ass$var == "phenotype", ]$path_id)
lapply(path_ids, function(path_id) return(boxplot(as.numeric(data[path_id, aut_samples]), as.numeric(data[path_id, nt_samples]), main = path_id)))


#################################################
############# Biplot ############################
#################################################

ord <- ordinate(ps_tmp, method = "CCA", distance = "euclidean")
p <- plot_ordination(ps_tmp, ord, type = "taxa", color = "Pathway")

plot_data <- data.frame(ord$CA$v)
plot_data$path_name <- as.character(ps_tmp@tax_table[rownames(plot_data), ])
to_label <- (plot_data$CA1 < -5 & plot_data$CA2 < -10) | plot_data$path_name == "Sesquiterpenoid and triterpenoid biosynthesis"
plot_data$path_name[!to_label] = ""
label_map = aes(x = 1.2 * CA1, y = 1.2 * CA2, shape = NULL, color = NULL,
                label = path_name)
p <- ggplot(plot_data) +
  geom_point(aes(x = CA1, y = CA2)) +
  geom_text_repel(label_map, size = 4, box.padding = 2.5)
p

# Now add the environmental variables as arrows
arrowmat = envfit(ord, env = df, na.rm = T)
pvals <- arrowmat$vectors[[4]]

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat$vectors[[1]]), arrowmat$vectors[[1]])
arrowdf <- arrowdf[pvals < .2 | rownames(arrowdf) == "phenotype", ]
# Define the arrow aesthetic mapping
inflate = 10
arrow_map = aes(xend = inflate * CA1, yend = inflate*CA2, x = 0, y = 0, shape = NULL, color = NULL,
                label = labels)
label_map = aes(x = inflate * CA1, y = inflate * CA2, shape = NULL, color = NULL,
                label = labels)
arrowhead = arrow(length = unit(0.05, "npc"))
p1 = p + geom_segment(arrow_map, size = 0.5, data = arrowdf, color = "gray",
                      arrow = arrowhead) + geom_text_repel(label_map, size = 3, data = arrowdf, color ="purple")
p1




#####################################
### Does adding microbiome data   ###
### improve performance - linear  ###
#####################################




###############################################
##############  Metatranscriptomics  ##########
###############################################
ps_trans <- readRDS("C:/Users/ctata/Documents/Lab/M3/PS_kegg_pathway_MTT.rds")
"ko00909"

otu_table(ps_trans) <- normalize_rle(ps_trans@otu_table)
ord <- ordinate(ps_trans, "PCoA", "bray")
plot_ordination(ps_trans, ord, type = "sample", color = "phenotype", label = "Biospecimen.Name" )

ps_trans <- addNumericMetadata(ps_trans)

last_bac_gene <- "ko03430"
remove <- taxa_names(ps_trans)[(which(taxa_names(ps_trans) == "ko03430") + 1) : ntaxa(ps_trans)] #removing all the human associated pathways
ps_trans <- prune_taxa(!(taxa_names(ps_trans) %in% remove), ps_trans)


pathway_meta_ass <- linearModeling(ps_trans)
aut_samples <- sample_names(ps_trans)[ps_trans@sam_data$phenotype == 1]
nt_samples <- sample_names(ps_trans)[ps_trans@sam_data$phenotype == 0]
data <- data.frame(ps_trans@otu_table)

path_ids <- as.character(pathway_meta_ass[pathway_meta_ass$var == "phenotype", ]$path_id)
lapply(path_ids, function(path_id) return(boxplot(as.numeric(data[path_id, aut_samples]), as.numeric(data[path_id, nt_samples]), main = path_id)))

boxplot(as.numeric(data["ko01502", aut_samples]), as.numeric(data["ko01502", nt_samples]))
boxplot(as.numeric(data["ko01502", aut_samples]), as.numeric(data["ko01502", nt_samples]))


pvals_wilcox <- apply(ps_trans@otu_table, 1, function(x){
  return(wilcox.test(as.numeric( x[aut_samples]), as.numeric(x[nt_samples]))$p.val)
})

#######################################
#### linear modeling  #################
#######################################

#Does adding in microbiome pathway information change what we know based on metadata?



######################################
##### Dirichlet multinomial  #########
######################################






#ASVs
#Kos
#pathways normalized

#linear model for each
#Linear models with only permanova variables included
#binomial model


#Metabolite Model
  #linear model
  #binomial model
  #* linear model with permanova variables

#Metagenome Model
  #

############################################################################################################################
############################################################################################################################
#######################################################   Metabolomics  #####################################################
############################################################################################################################
############################################################################################################################

metabolomics <- read.csv("data/metabol/metabolomics_3.csv")
rownames(metabolomics) <- metabolomics$biochemical
metabolite2pathway <- data.frame(metabolomics$sub_pathway, row.names = metabolomics$biochemical)

metabolomics <- metabolomics[ , 9:ncol(metabolomics)]
sample_manifest <- read.csv("data/metabol/metabolomics_sample_manifest.csv", row.names = 1)
colnames(metabolomics) <- as.character(sample_manifest[colnames(metabolomics), "CLIENT_SAMPLE_ID"])

#read in as a tmp mapping file
mapping <- read.delim("data/metabol/mapping_metabolomics.csv", sep = ",")
rownames(mapping) <- mapping$CLIENT_SAMPLE_ID

ps_metabolite <- phyloseq(otu_table(metabolomics, taxa_are_rows = T), sample_data = sample_data(mapping))
ps_metabolite <- addNumericMetadata(ps_metabolite)
saveRDS(ps_metabolite, "data/metabol/ps_metabolite_tmp.rds")

ord <- ordinate(ps_metabolite, "PCoA", "bray")
plot_ordination(ps_metabolite, ord, type = "sample", color = "phenotype")

#### Binomial Model Metabolomics ######
null_dist <- getNullDist(ps_metabolite, iter = 100000)
actual_diffs <- getActualDiffs(ps_metabolite)
bin_pos_metabolite <- findSig_Binomial_old(null_dist, actual_diffs, pos = T, type = "metabolite")
bin_neg_metabolite <- findSig_Binomial_old(null_dist, actual_diffs, pos = F, type = "metabolite")



sapply(bin_neg_metabolite$ko_ids, function(ko_id){
  plotNullHist(null_dist[, ko_id], actual_diffs[,ko_id], ko = ko_id, path = metabolite2pathway[ko_id, ])
})



#pos
#alpha-lipoate


#### LOO linear model Metabolomics
res_linearModel <- linearModeling_weightMetadata(ps_metabolite, type = "metabolite", reg = 2)
line_pos_metabolite <- res_linearModel[[1]]
line_neg_metabolite <- res_linearModel[[2]]

#write metabolite results
res <- list(bin_neg_metabolite$ko_ids, line_neg_metabolite$neg_kos)
df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df) <- c( "Bin_metabolite", "Line_metabolite")
for(i in seq(1, length(res))){
  for(var in res[[i]]){
    df[var, i] <- 1
  }
}
df <- cbind(getMetaboliteInfo(rownames(df)), df)
df <- df[rownames(df) != "NA", ]
write.table(df, "metabolite_res.tsv", sep = "\t", col.names = NA, row.names = T)


### Permanova
perm_metabolite <- runPermanova(ps_metabolite)
perm_metabolite[perm_metabolite$`Pr(>F)` < .05, ]
sig_vars_metabolite <- rownames(perm)[perm$`Pr(>F)` < .07 & perm$disp_pval > .05]

#### LOO linear model
res_linearModel <- linearModeling_weightMetadata(ps_metabolite, type = "metabolite", reg = 2.3, vars = sig_vars_metabolite)
line_pos_perm_metabolite <- res_linearModel[[1]]
line_neg_perm_metabolite <- res_linearModel[[2]]


#### Classifier

meta_impute <- apply(ps_metabolite@sam_data %>% select(-c(familyID, phenotype, host_name)), 2, function(col){
  col[is.na(col)] = mean(col, na.rm = T)
  return(col)
})

library(ROCR)
plotModel <- function(model, y){
  preds <- predict(model, meta_impute, y, type = "prob")
  perf <- performance(predictions(preds, y), "tpr", "fpr")
  roc_obj <- roc(perf)
  plot(roc_obj)
  print(auc(roc_obj))
}

metadata_model <- trainRf(meta_impute, map = ps_metabolite@sam_data)
model <- randomForest(meta_impute, ps_metabolite@sam_data$phenotype, mtry = metadata_model[[1]]$best)
plotModel(model, y = ps_metabolite@sam_data$phenotype)

metabolite_input <- t(ps_metabolite@otu_table)
metabolite_model <- trainRf(metabolite_input, map = ps_metabolite@sam_data)

metabolite_meta_input <- cbind(t(ps_metabolite@otu_table), meta_impute)
metabolite_meta_model <- trainRf(metabolite_meta_input, map = ps_metabolite@sam_data)
#  18    0.6146804  0.2294162


############################################################################################################################
############################################################################################################################
#######################################################   KOs  #####################################################
############################################################################################################################
############################################################################################################################

ps_kegg <- readRDS("C:/Users/ctata/Documents/Lab/M3/PS_kegg_fun_MTG_age_filtered.rds")
ps_kegg <- addNumericMetadata(ps_kegg)


ord <- ordinate(ps_kegg, "PCoA", "bray")
plot_ordination(ps_kegg, ord, type = "sample", color = "phenotype")

ps_kegg_norm <- phyloseq(otu_table(normalize_rle(ps_kegg@otu_table)), sample_data(ps_kegg@sam_data))
saveRDS(ps_kegg_norm, "ps_metagenome.rds")


ord <- ordinate(ps_kegg, "PCoA", "bray")
plot_ordination(ps_kegg, ord, type = "sample", color = "phenotype")
#not sure about the noramlization

#### Binomial Model ######
null_dist <- getNullDist(ps_kegg, iter = 10000)
actual_diffs <- getActualDiffs(ps_kegg)
bin_pos_ko <- findSig_Binomial(null_dist, actual_diffs, pos = T, type = "ko")
bin_neg_ko <- findSig_Binomial(null_dist, actual_diffs, pos = F, type = "ko")

#### LOO linear model
res_linearModel <- linearModeling_weightMetadata(ps_kegg, type = "ko", reg = 2.3)
line_pos_ko <- res_linearModel[[1]]
line_neg_ko <- res_linearModel[[2]]
line_pos_ko_filt <- line_pos_ko[line_pos_ko$Freq > 30, ]
line_neg_ko_filt <- line_neg_ko[line_neg_ko$Freq > 30, ]

line_neg_ko$def <- sapply(line_neg_ko$neg_kos, getKOName)



#write ko results
res <- list(bin_neg_ko$ko_ids, line_neg_ko_filt$neg_kos)
df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df) <- c( "Bin_ko", "Line_ko")
for(i in seq(1, length(res))){
  for(var in res[[i]]){
    df[var, i] <- 1
  }
}
write.table(df, "ko_res.tsv", sep = '\t', col.names = NA, row.names = T)






### Permanova
perm_ko <- runPermanova(ps_kegg)
write.csv(perm_ko[perm_ko$`Pr(>F)` < .05, ], "perm_ko.csv", quote = F)

sig_vars_kegg <- rownames(perm)[perm$`Pr(>F)` < .07 & perm$disp_pval > .05]

#### LOO linear model
res_linearModel <- linearModeling_weightMetadata(ps_kegg, type = "ko", reg = 2, vars = sig_vars_kegg)
line_pos_perm_ko <- res_linearModel[[1]]
line_neg_perm_ko <- res_linearModel[[2]]




#### Classifier
library(ROCR)
plotModel <- function(model, data, y){
  preds <- predict(model, data, y, type = "prob")
  preds_obj <- prediction(preds[,2], y)
  perf <- performance(preds_obj, "tpr", "fpr")
  plot(perf)
  print(performance(preds_obj, measure = "auc")@y.values[[1]])
}

meta_impute <- apply(ps_kegg@sam_data %>% select(-c(familyID, phenotype, host_name)), 2, function(col){
  col[is.na(col)] = mean(col, na.rm = T)
  return(col)
})

#ko_input <- t(ps_kegg@otu_table)
#ko_model <- trainRf(ko_input, map = ps_kegg@sam_data)


ko_meta_input <- cbind(t(ps_kegg@otu_table), meta_impute)
ko_meta_model <- trainRf(ko_meta_input, map = ps_kegg@sam_data)
ko_meta_model <- randomForest(ko_meta_input, ps_kegg@sam_data$phenotype, mtry = ko_meta_model$bestTune[[1]])
plotModel(model = ko_meta_model, data = ko_meta_input, y = ps_kegg@sam_data$phenotype)


############################################################################################################################
############################################################################################################################
#######################################################   MTT  #####################################################
############################################################################################################################
############################################################################################################################

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

#### Binomial Model ######
null_dist <- getNullDist(ps_mtt, iter = 10000)
actual_diffs <- getActualDiffs(ps_mtt)
bin_pos_mtt <- findSig_Binomial(null_dist, actual_diffs, pos = T, type = "ko")
bin_neg_mtt <- findSig_Binomial(null_dist, actual_diffs, pos = F, type = "ko")

#### LOO linear model
res_linearModel <- linearModeling_weightMetadata(ps_mtt, type = "ko", reg = 2.3)
line_pos_mtt <- res_linearModel[[1]]
line_neg_mtt <- res_linearModel[[2]]
line_pos_mtt_filt <- line_pos_mtt[line_pos_mtt$Freq > 30, ]
line_neg_mtt_filt <- line_neg_mtt[line_neg_mtt$Freq > 30, ]


#write transcriptome results
res <- list(bin_neg_mtt$ko_ids, line_neg_mtt_filt$neg_kos)
df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df) <- c( "Bin_mtt", "Line_mtt")
for(i in seq(1, length(res))){
  for(var in res[[i]]){
    df[var, i] <- 1
  }
}
write.table(df, "mtt_res.tsv", sep = '\t', col.names = NA, row.names = T)


### Permanova
perm_mtt <- runPermanova(ps_mtt)
write.csv(perm_mtt[perm_mtt$`Pr(>F)` < .05, ], "perm_mtt.csv", quote = F)




############################################################################################################################
############################################################################################################################
#######################################################   Permutation Metagenome  ##########################################
############################################################################################################################
############################################################################################################################
ps_metabolite <- readRDS("ps_metabolite.rds")
otu_table(ps_metabolite) <- log(ps_metabolite@otu_table)
variances <- apply(ps_metabolite@otu_table, 1, function(x) return(var(x)))
#ps_metabolite <- prune_taxa(variances > 0.5, ps_metabolite)
actual_diffs <- getActualDiffs(ps_metabolite)
bin_pos_metabolite <- findSig_Binomial_removezeros(actual_diffs, pos = T, type = "ko")
bin_neg_metabolite <- findSig_Binomial_removezeros(actual_diffs, pos = F, type = "ko")
View(bin_neg_metabolite[order(bin_neg_metabolite$pvals), ])

bin_neg_metabolite <- bin_neg_metabolite[bin_neg_metabolite$pvals < .05, ]
bin_pos_metabolite <- bin_pos_metabolite[bin_pos_metabolite$pvals < .05, ]

#write metabolite results
res <- list(bin_pos_metabolite$ko_ids, line_pos_metabolite$neg_kos)
df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df) <- c( "Bin_metabolite", "Line_metabolite")
for(i in seq(1, length(res))){
  for(var in res[[i]]){
    df[var, i] <- 1
  }
}
df <- cbind(getMetaboliteInfo(rownames(df)), df)
df <- df[rownames(df) != "NA", ]
write.table(df, "metabolite_res_pos.tsv", sep = "\t", col.names = NA, row.names = T)


## Put together pos/neg metabolites
neg <- read.delim("metabolite_res_neg.tsv", header = T, row.names = 1)
pos <- read.delim("metabolite_res_pos.tsv", header = T, row.names = 1)
df <- rbind(neg, pos)
df$direction <- c(rep(-1, nrow(neg)), rep(1, nrow(pos)))
df <- df %>% select(-c("Bin_metabolite", "Line_metabolite"))
write.table(df, "metabolite_res_combined.tsv", sep = "\t", col.names = NA, row.names = T)


ps_mtt <- readRDS("ps_metatranscriptome.rds")
#null_dist <- getNullDist(ps_mtt, iter = 10000)
actual_diffs <- getActualDiffs(ps_mtt)
bin_pos_mtt<- findSig_Binomial_removezeros(actual_diffs = actual_diffs, pos = F, type = "ko")
bin_neg_mtt<- findSig_Binomial_removezeros(actual_diffs = actual_diffs, pos = F, type = "ko")
bin_pos_mtt <- bin_pos_mtt[bin_pos_mtt$pvals < .05, ]
bin_neg_mtt <- bin_neg_mtt[bin_neg_mtt$pvals < .05, ]

#write transcriptome results
res <- list(bin_neg_mtt$ko_ids, line_neg_mtt_filt$neg_kos)
df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df) <- c( "Bin_mtt", "Line_mtt")
for(i in seq(1, length(res))){
  for(var in res[[i]]){
    df[var, i] <- 1
  }
}
write.table(df, "mtt_res_neg.tsv", sep = '\t', col.names = NA, row.names = T)


## Put together pos/neg metabolites
neg <- read.delim("mtt_res_neg.tsv", header = T, row.names = 1)
pos <- read.delim("mtt_res_pos.tsv", header = T, row.names = 1)
df <- rbind(neg, pos)
df$direction <- c(rep(-1, nrow(neg)), rep(1, nrow(pos)))
df <- df %>% select(-c("Bin_mtt", "Line_mtt"))
write.table(df, "mtt_res_combined.tsv", sep = '\t', col.names = NA, row.names = T)
df_mtt <- df


ps_mgg <- readRDS("ps_metagenome.rds")
actual_diffs <- getActualDiffs(ps_mgg)
bin_pos_ko <- findSig_Binomial_removezeros(actual_diffs, pos = T, type = "ko")
bin_neg_ko <- findSig_Binomial_removezeros(actual_diffs, pos = F, type = "ko")
bin_pos_ko <- bin_pos_ko[bin_pos_ko$pvals < .05, ]
bin_neg_ko <- bin_neg_ko[bin_neg_ko$pvals < .05, ]

#write transcriptome results
res <- list(bin_neg_ko$ko_ids, line_neg_ko_filt$neg_kos)
df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df) <- c( "Bin_ko", "Line_ko")
for(i in seq(1, length(res))){
  for(var in res[[i]]){
    df[var, i] <- 1
  }
}
write.table(df, "mtg_res_neg.tsv", sep = '\t', col.names = NA, row.names = T)

## Put together pos/neg metagneome
neg <- read.delim("mtg_res_neg.tsv", header = T, row.names = 1)
pos <- read.delim("mtg_res_pos.tsv", header = T, row.names = 1)
df <- rbind(neg, pos)
df$direction <- c(rep(-1, nrow(neg)), rep(1, nrow(pos)))
df <- df %>% select(-c("Bin_ko", "Line_ko"))
write.table(df, "mtg_res_combined.tsv", sep = '\t', col.names = NA, row.names = T)
df_mtg <- df


mtg_rownames <- rownames(df_mtg)
mtt_rownames <- rownames(df_mtt)
mtg_direction <- df_mtg$direction
mtt_direction <- df_mtt$direction
df_combo <- data.frame(ko_names = c(mtg_rownames, mtt_rownames), direction =  c(mtg_direction, mtt_direction))

df_combo$mtg <- rep(0, nrow(df_combo))
df_combo$mtg[df_combo$ko_names %in% rownames(df_mtg)] <- 1
df_combo$mtt <- rep(0, nrow(df_combo))
df_combo$mtt[df_combo$ko_names%in% rownames(df_mtt)] <- 1

messed_ids <- intersect(rownames(df_mtg)[df_mtg$direction == 1], rownames(df_mtt)[df_mtt$direction == -1])
messed_ids <- c(messed_ids, intersect(rownames(df_mtg)[df_mtg$direction == -1], rownames(df_mtt)[df_mtt$direction == 1]))
df_combo <- df_combo[!(rownames(df_combo) %in% messed_ids), , drop = FALSE]
write.table(df_combo, "res_combined_mtg_mtt.tsv",  sep = '\t',  row.names = F)


############################################################################################################################
########################################################################################################################
#######################################################   Piphillan  #####################################################
############################################################################################################################
############################################################################################################################
ps_piph <- readRDS("ps_piph.rds")
#keep only taxa that are present in all timepoints for each sample
ps_piph <- prune_samples(!duplicated(ps_piph@sam_data$host_name), ps_piph)

actual_diffs <- getActualDiffs(ps_piph)
bin_pos_piph <- findSig_Binomial_removezeros(actual_diffs, pos = T, type = "ko")
bin_neg_piph <- findSig_Binomial_removezeros(actual_diffs, pos = F, type = "ko")
bin_pos_piph <- bin_pos_piph[bin_pos_piph$pvals < .05, ]
bin_neg_piph <- bin_neg_piph[bin_neg_piph$pvals < .05, ]

#### LOO linear model
res_linearModel <- linearModeling_weightMetadata(ps_piph, type = "ko", reg = 2.3)
line_pos_piph <- res_linearModel[[1]]
line_neg_piph <- res_linearModel[[2]]
line_pos_piph_filt <- line_pos_piph[line_pos_piph$Freq > 30, ]
line_neg_piph_filt <- line_neg_piph[line_neg_piph$Freq > 30, ]

writeResults <- function(bin_res, line_res, outname){
  res <- list(bin_res$ko_ids, line_res$neg_kos)
  df <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(df) <- c( "Bin_ko", "Line_ko")
  for(i in seq(1, length(res))){
    for(var in res[[i]]){
      df[var, i] <- 1
    }
  }
  write.table(df, outname, sep = '\t', col.names = NA, row.names = T)
}
writeResults(bin_neg_piph, line_neg_piph_filt, "piph_res_neg.tsv")
writeResults(bin_pos_piph, line_pos_piph_filt, "piph_res_pos.tsv")

## Put together pos/neg metagneome
neg <- read.delim("piph_res_neg.tsv", header = T, row.names = 1)
pos <- read.delim("piph_res_pos.tsv", header = T, row.names = 1)
df <- rbind(neg, pos)
df$direction <- c(rep(-1, nrow(neg)), rep(1, nrow(pos)))
df <- df %>% select(-c("Bin_ko", "Line_ko"))
write.table(df, "piph_res_combined.tsv", sep = '\t', col.names = NA, row.names = T)





############################################################################################################################
############################################################################################################################
#######################################################   Pathways  #####################################################
############################################################################################################################
############################################################################################################################
ps_pathway <- readRDS("ps_kegg_pathway_rlenorm_agefilt.rds") #already normalized and numerized

######kegg

#### Binomial Model ######
null_dist <- getNullDist(ps_pathway, iter = 10000)
actual_diffs <- getActualDiffs(ps_pathway)
bin_pos_path <- findSig_Binomial(null_dist, actual_diffs, pos = T, type = "pathway")
bin_neg_path <- findSig_Binomial(null_dist, actual_diffs, pos = F, type = "pathway")

#### LOO linear model
res_linearModel <- linearModeling_weightMetadata(ps_pathway, type = "pathway", reg = 2)
line_pos_path <- res_linearModel[[1]]
line_neg_path <- res_linearModel[[2]]


################## Combine results ##########



#write pathway results
res <- list(bin_neg_path$ko_ids, line_neg_path$neg_kos)
df <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(df) <- c( "Bin_pathway", "Line_pathway")
for(i in seq(1, length(res))){
  for(var in res[[i]]){
    df[var, i] <- 1
  }
}
write.table(df, "pathway_res.tsv", sep = '\t', col.names = NA, row.names = T)


res <- list(bin_neg_ko$ko_ids, line_neg_ko$neg_kos, bin_neg_metabolite$ko_ids, line_neg_metabolite$neg_kos)
df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(df) <- c("Bin_ko", "Line_ko", "Bin_metabolite", "Line_metabolite")
for(i in seq(1, length(res))){
  for(var in res[[i]]){
    df[var, i] <- 1
  }
}

best <- df[rowSums(df, na.rm = T) >=2, ]
colnames(best) <- c("Metagenome_Binomial", "Metagenome_Linear", "Metabolome_Binomial", "Metabolome_Linear")
best <- best[order(rowSums(best, na.rm = T), decreasing = T), ]

######## get ids ###########
getMetaboliteInfo <- function(ids){
  metabolomics_SG <- read.csv("C:/Users/ctata/Documents/Lab/M3/metabolomics_SG.csv", header=FALSE)
  metabolite_ids <- metabolomics_SG[14:nrow(metabolomics_SG), 1:13]
  rownames(metabolite_ids) <- metabolite_ids[,2]
  metabolite_ids <- metabolite_ids[ , 3:13]
  colnames(metabolite_ids) <- c("SUPER_PATHWAY", "SUB_PATHWAY", "COMP_ID", "PLATFORM", "CHEMICAL_ID", "RI", "MASS", "PUBCHEM", "CAS", "KEGG", "HMDB")
  return(metabolite_ids[ids, ])
}



###### Keep
keep <- metabolite_ids[rownames(best),]
colnames(keep) <- c("SUPER_PATHWAY", "SUB_PATHWAY", "COMP_ID", "PLATFORM", "CHEMICAL_ID", "RI", "MASS", "PUBCHEM", "CAS", "KEGG", "HMDB")

final_metabolites <- cbind(keep, best)
final_metabolites <- final_metabolites[rownames(final_metabolites)!="NA", ]
View(final_metabolites)
write.table(final_metabolites, "metabolite_res.tsv", sep = '\t', col.names = NA, row.names = T)







#Correlate with phenotype
mapping <- read.csv("G:/Shared drives/M3-shared/Data/Metabolomics/M3_metabolome_meta_data.csv")
ps_metabolite <- readRDS("ps_metabolite.rds")
sample_names(ps_metabolite) %in% as.character(mapping$Biospecimen.Barcode)
rownames(mapping) <- mapping$Biospecimen.Barcode
mapping <- mapping[sample_names(ps_metabolite), ]
phenotypes <- mapping[, 84 : (84 +14)]
phenotypes <- phenotypes[!is.na(phenotypes$Language.ability.and.use), ]

language <- list("Able to speak fluently" = 0, "Phrase speech" = 1, "Single word speech" = 2, "Little to no speech" = 3)
conversation <- list("Able to have conversation" = 0, "Difficulty with conversation" = 2,"Limited conversation ability" = 3, "Cannot have a conversation" = 4)
speech <- list("Understands nearly all words"  = 0, "Understands most words" = 1, "Understands many words" = 2, "Understands about half of words" = 3, "Understands few or no words"  = 4)
timing <- list("Regularly" = 0, "Sometimes" = 1, "Rarely" = 2, "Never" = 3, "No opportunity to play with other children" = NA)
eye_contact <- list("Consistent eye contact"=0, "Some eye contact"=1, "Little of no eye contact" = 2)
rep_motion <- list("Never" = 0, "Sometimes" = 1, "Regularly" = 2)
#development scale is not clearly ordered
sleep <- list("Healthy sleep pattern" = 0, "Some sleep difficulties" = 1, "Constant sleep difficulties" = 2)
sounds <- list("Not bothered by typical sounds" = 0, "Sensitive to typical sounds" = 1, "Highly sensitive to typical sounds" = 2 )
self_harm <- list("No self-injurious behavior" = 0 , "Mild self-harming behavior" = 1, "Dangerous or frequent self-harming behavior" =2)
gi <- list("No issues" = 0, "Sometimes" = 1, "Continuous" = 2)
imitate <- list("Imitates actions or gestures of others" = 0, "Imitates others when prompted" = 1, "Does not imitate others" =2)

phenotypes_numeric <- phenotypes
phenotypes_numeric$language <- unlist(language[phenotypes$Language.ability.and.use])
phenotypes_numeric$conversation <- unlist(conversation[phenotypes$Conversation.ability])
phenotypes_numeric$speech <- unlist(speech[phenotypes$Understands.speech])
phenotypes_numeric$play_imagine_alone <- unlist(timing[phenotypes$Plays.imaginatively.when.alone])
phenotypes_numeric$play_imagine_others <- unlist(timing[phenotypes$Plays.imaginatively.with.others])
phenotypes_numeric$play_group_others <- unlist(timing[phenotypes$Plays.in.a.group.with.others])
phenotypes_numeric$eye_contact <- unlist(eye_contact[phenotypes$Eye.contact.finding])
phenotypes_numeric$repetitive_motion <- unlist(rep_motion[phenotypes$Repetitive.motion])
phenotypes_numeric$shows_objects_others <- unlist(timing[phenotypes$Picks.up.objects.to.show.to.others])
phenotypes_numeric$sleep <- unlist(sleep[phenotypes$Sleep.pattern.finding])
phenotypes_numeric$sounds <- unlist(sounds[phenotypes$Response.to.typical.sounds])
phenotypes_numeric$self_harm <- unlist(self_harm[phenotypes$Self.injurious.behavior.finding])
phenotypes_numeric$gi <- unlist(gi[phenotypes$Gastrointestinal.problems..M3.])
phenotypes_numeric$imitate_behavior <- unlist(imitate[phenotypes$Imitation.behavior])
phenotypes_numeric <- phenotypes_numeric %>% select(c("language", "conversation", "speech", "play_imagine_alone", "play_imagine_others",
                                "play_group_others", "eye_contact", "repetitive_motion", "shows_objects_others", "sleep",
                                "sounds", "self_harm", "gi", "imitate_behavior")) #14


sample_data(ps_metabolite) <- phenotypes_numeric

## Do any of the metabolites of interest correlate with a particular phenotype
## Spearman cor(x, y, method = "spearman")
nominated_metabolites <- read.csv("C:/Users/ctata/Documents/Lab/M3/Nominated_metabolites.csv")
metabolites <- as.character(nominated_metabolites$metabolon_name)
metabolites <- metabolites[metabolites != ""]
metabolites <- metabolites[metabolites %in% taxa_names(ps_metabolite)]

ps_metabolite_pruned <- prune_taxa(metabolites, ps_metabolite)
ps_metabolite_pruned@otu_table <- log(ps_metabolite_pruned@otu_table)
correlations <- apply(ps_metabolite_pruned@otu_table, 1, function(metabolite_counts){
  cor(as.numeric(unlist(metabolite_counts)), ps_metabolite_pruned@sam_data, use = "na.or.complete")
})
rownames(correlations) <- colnames(ps_metabolite_pruned@sam_data)
library(pheatmap)
tiff("metabolite_corr_phenotype.tif", width = 12, height = 6, units = 'in', res = 300)
pheatmap(correlations, angle_col = 90)
dev.off()

trace(pheatmap::draw_colnames, edit = T) #Edit pheatmap functions for myself :O



##### Boxplots
aut <- prune_samples(ps_metabolite@sam_data$phenotype == "A", ps_metabolite )
nt <- prune_samples(ps_metabolite@sam_data$phenotype == "N", ps_metabolite)
metabolite_nominations <- read.delim("C:/Users/ctata/Documents/Lab/M3/metabolite_res_combined.tsv", row.names=1)
octanoate <- taxa_names(aut)[grepl("octanoate", taxa_names(aut))]
hexanoate <- taxa_names(aut)[grepl("hexanoate", taxa_names(aut))]
decanoate <- taxa_names(aut)[grepl("decanoate", taxa_names(aut))]
oates <- c(octanoate, hexanoate, decanoate)

boxplot(as.numeric(aut@otu_table[octanoate[1],]), as.numeric(nt@otu_table[octanoate[1], ]))
boxplot(as.numeric(aut@otu_table[octanoate[2],]), as.numeric(nt@otu_table[octanoate[2], ]))
boxplot(as.numeric(aut@otu_table[octanoate[3],]), as.numeric(nt@otu_table[octanoate[3], ]))
boxplot(as.numeric(colSums(aut@otu_table[oates,])), as.numeric(colSums(nt@otu_table[oates, ])))


