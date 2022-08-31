setwd("C:/Users/ctata/Documents/Lab/M3")
library(dplyr)
library(phyloseq)
library(glmnet)
library(cvAUC)
library(lsa)
library(DESeq2)


###Read in phyloseq obj
ps <- readRDS("ps_not_norm_comp_pass_min_postDD.Rda")
print(nsamples(ps))
seqRead_thresh <- mean(sample_sums(ps)) + 3*sd(sample_sums(ps))
ps <- prune_samples(sample_sums(ps) < seqRead_thresh, ps)
print(nsamples(ps))
seqRead_thresh <- mean(sample_sums(ps)) - 3*sd(sample_sums(ps))
ps <- prune_samples(sample_sums(ps) > seqRead_thresh, ps)
print(nsamples(ps))


#Deal with mapping file
map_ps <- ps@sam_data
map_ps$sampleID <- substr(map_ps$Host.Name, 1, 5)
map_ps <- map_ps[!duplicated(map_ps$sampleID), ] #use this to create clean ml mapping file
metadata_keep <- c("Within.study.sampling.date..Biospecimen." , "Family.group.ID..Biospecimen.",
                   "Pet.dog..Biospecimen.", "Pet.cat..Biospecimen." ,
                   "Minimum.time.since.antibiotics..Biospecimen.", "Stool.frequency..Biospecimen.",
                   "Born.by.caesarian.section..Biospecimen.", "Other.GI.symptoms..M3...Biospecimen.",
                   "Biological.sex..Biospecimen.", "Gluten.allergy..Biospecimen.",
                   "Non.celiac.gluten.sensitivity..Biospecimen." , 
                   "Whole.grain..consumption.frequency...Biospecimen.", "Fermented.vegetable..consumption.frequency...Biospecimen.",               
                   "Dairy..consumption.frequency...Biospecimen." ,"Fruit..consumption.frequency...Biospecimen.",
                   "Meals.prepared.at.home..consumption.frequency...Biospecimen.","Ready.to.eat.meals..consumption.frequency...Biospecimen.",
                   "Meat..consumption.frequency...Biospecimen.", "Olive.oil.used.in.cooking..M3...Biospecimen." ,
                   "Seafood..consumption.frequency...Biospecimen.","Sweetened.drink..consumption.frequency...Biospecimen." ,
                   "Vegetable..consumption.frequency...Biospecimen.", "Restaurant.prepared.meals..consumption.frequency...Biospecimen.",
                   "Sugary.food..consumption.frequency...Biospecimen." ,"Diet.during.infancy..Biospecimen." ,
                   "Lactose.intolerance..Biospecimen." , "Born.prematurely..Biospecimen.",
                   "Environmental.tobacco.smoke.exposure..Biospecimen.","Vitamin.B.complex.supplement..consumption.frequency...Biospecimen.",
                   "Vitamin.D..consumption.frequency...Biospecimen.", "Recently.ill..Biospecimen.",
                   "Diarrhea..Biospecimen.","Constipation..Biospecimen.",
                   "Starchy.food..consumption.frequency...longitudinal...Biospecimen.","Bread..consumption.frequency...longitudinal...Biospecimen.",
                   "Natural.father.current.age..months...Biospecimen." ,"Natural.mother.current.age..months...Biospecimen." ,
                   "Flu.like.symptoms..Biospecimen.", "Fever..Biospecimen.", "phenotype", "Age..months.")

map_ps <- map_ps[ , metadata_keep]
rename_cols <- c("timepoint", "familyID",
                  "dog", "cat",
                  "time_since_antibiotics", "stool_freq", 
                  "csection", "other_GI_symptoms", 
                  "sex", "gluten_allergy",
                  "nonceliac_sensitivity",
                   "whole_grain", "fermented_vegetables",
                   "dairy", "fruit",
                   "meal_home_prep", "meal_ready_eat", 
                   "meat", "olive_oil", 
                   "seafood", "sweetened_drink",
                   "vegetable", "restaurant", 
                   "sugary_food", "infant_diet",
                   "lactose_intolerance", "prematurely_born", 
                   "env_tobacco",  "vitamin_B",
                   "vitamin_D", "recently_ill", 
                   "diarrhea", "constipation", 
                   "starchy_food",  "bread", 
                  "father_age", "mother_age", 
                  "flu_symptoms", "fever",
                  "phenotype", "age")

colnames(map_ps) <- rename_cols

embed <- read.table("embed/embedded_100.txt", sep = "\t")
timepoints <- substr(rownames(embed), 7, 8)
embed_1 <- embed[timepoints == 1, ]
embed_2 <- embed[timepoints == 2, ]
embed_3 <- embed[timepoints == 3, ]
rownames(embed_1) <- substr(rownames(embed_1), 1, 5)
rownames(embed_2) <- substr(rownames(embed_2), 1, 5)
rownames(embed_3) <- substr(rownames(embed_3), 1, 5)
ps_embed_1 <- phyloseq(otu_table(embed_1, taxa_are_rows = F), sample_data = sample_data(map))
ps_embed_2 <- phyloseq(otu_table(embed_2, taxa_are_rows = F), sample_data = sample_data(map))
ps_embed_3 <- phyloseq(otu_table(embed_3, taxa_are_rows = F), sample_data = sample_data(map))

seqtab <- read.table("seqtab.txt", sep = "\t")
timepoints <- substr(rownames(seqtab), 7, 8)
seqtab_1 <- seqtab[timepoints == 1, ]
seqtab_2 <- seqtab[timepoints == 2, ]
seqtab_3 <- seqtab[timepoints == 3, ]
rownames(seqtab_1) <- substr(rownames(seqtab_1), 1, 5)
rownames(seqtab_2) <- substr(rownames(seqtab_2), 1, 5)
rownames(seqtab_3) <- substr(rownames(seqtab_3), 1, 5)
ps_seqtab_1 <- phyloseq(otu_table(asinh(seqtab_1), taxa_are_rows = F), sample_data = sample_data(map))
ps_seqtab_2 <- phyloseq(otu_table(asinh(seqtab_2), taxa_are_rows = F), sample_data = sample_data(map))
ps_seqtab_3 <- phyloseq(otu_table(asinh(seqtab_3), taxa_are_rows = F), sample_data = sample_data(map))


#deseq normalize
deSeqNorm <- function(ps){
  library(DESeq2)
  ps_dds <- phyloseq_to_deseq2(ps, ~ phenotype )
  ps_dds <- estimateSizeFactors(ps_dds, type = "poscounts")
  ps_dds <- estimateDispersions(ps_dds)
  abund <- getVarianceStabilizedData(ps_dds)
  abund <- abund + abs(min(abund)) #don't allow deseq to return negative counts
  ps_deSeq <- phyloseq(otu_table(t(abund), taxa_are_rows = F), sample_data(ps))
  return(ps_deSeq)
}

ps_deseq_1 <- deSeqNorm(phyloseq(otu_table(seqtab_1, taxa_are_rows = F), sample_data = sample_data(map)))
ps_deseq_2 <- deSeqNorm(phyloseq(otu_table(seqtab_2, taxa_are_rows = F), sample_data = sample_data(map)))
ps_deseq_3 <- deSeqNorm(phyloseq(otu_table(seqtab_3, taxa_are_rows = F), sample_data = sample_data(map)))


###################################
#########MACHINE LEARNING##########
##################################



################GLMNET################
#####################

library(ROCR)

getMLInput <- function(otu, mapping, test_pairs, include_meta = T, meta_only = F){
  set.seed(0)


  test <- mapping$familyID %in% test_pairs
  train <- !(mapping$familyID %in% test_pairs)
  
  y_train <- mapping$phenotype[train]
  y_test <- mapping$phenotype[test]
  y_train <- ifelse(y_train == "A", 1, -1)
  y_test <- ifelse(y_test == "A", 1, -1)
  
  familyID <- mapping$familyID
  
  mapping <- mapping %>% select(-c("familyID", "phenotype"))
  father_age <- as.numeric(mapping$father_age)
  mapping <- apply(mapping, 2, function(x){
    x <- unlist(x)
    tmp <- (x - mean(x)) / sd(x)
    return( tmp)
  })
  data <- data.frame(otu)
  #data <- data[rownames(data) %in% rownames(mapping), ]
  

  if(include_meta){
    data <- cbind(data, mapping)
  }
  if(meta_only){
    data <- mapping
    tmp_df <- apply(mapping, 2, function(x){
      return(as.numeric(x) * father_age)
    })
    data <- cbind(data, tmp_df)
  }
  
  data_train <- data[train, ]
  data_test <- data[test, ]
  return(list(data_train, data_test, y_train, y_test, familyID[train], familyID[test]))
}


getMLInput_full <- function(otu, mapping, include_meta = T, meta_only = F){
  set.seed(0)
  
  y <- mapping$phenotype
  y <- ifelse(y == "A", 1, -1)
  
  mapping <- mapping %>% select(-c("familyID", "phenotype"))
  
  mapping <- apply(mapping, 2, function(x){
    x <- unlist(x)
    tmp <- (x - mean(x)) / sd(x)
    return( tmp)
  })
  data <- data.frame(otu)
  
  if(include_meta){
    data <- cbind(data, mapping)
  }
  if(meta_only){
    data <- mapping
    
  }
  
  return(list(data, y))
}




#Logistic regression
predictAut <- function(otu, mapping, test_pairs, col = "black", include_meta= T, meta_only = F){
  tmp <- getMLInput(otu, mapping, test_pairs, include_meta = include_meta, meta_only = meta_only)
  data_train <- tmp[[1]]
  data_test <- tmp[[2]]
  y_train <- tmp[[3]]
  y_test <- tmp[[4]]
  
  cv <- cv.glmnet(x = as.matrix(data_train), y = as.factor(y_train), type.measure = "deviance", family = "binomial")
  lambda <- cv$lambda.min
  if(max(cv$lambda) == lambda){
    lambda <- cv$lambda[30]
  }
  
  acc <- sum((predict(cv, newx = as.matrix(data_test), s = lambda) > 0) == y_test) / length(y_test)
  plot(cv)

  preds_prob <- predict(cv, newx = as.matrix(data_test), s = lambda, type = "response") #5, #12
  pred <- prediction(preds_prob, y_test)
  
  # calculate probabilities for TPR/FPR for predictions
  perf <- performance(pred,"tpr","fpr")
  print(paste("AUC: ", unlist(performance(pred,"auc")@y.values))) # shows calculated AUC for model
  
  return(list(preds_prob = preds_prob, y_test = y_test, cv = cv))
}


include_meta = T
set.seed(0)
ps_seqtab_use <- ps_seqtab_1
ps_embed_use <- ps_embed_1
ps_deseq_use <- ps_deseq_1
if(include_meta){
  sample_data(ps_seqtab_use) <- ps_seqtab_use@sam_data[complete.cases(ps_seqtab_use@sam_data), ]
  sample_data(ps_embed_use) <- ps_embed_use@sam_data[complete.cases(ps_embed_use@sam_data), ]
  sample_data(ps_deseq_use) <- ps_deseq_use@sam_data[complete.cases(ps_deseq_use@sam_data), ]
}
pairs <- unique(ps_seqtab_use@sam_data$familyID)
increment <- length(pairs) / 5
folds <- list(a = pairs[1: increment], b = pairs[(increment + 1): (2*increment)], 
              c = pairs[(2*increment + 1):(3*increment)], d = pairs[(3*increment + 1): (4*increment)])


otu_table(ps_embed_use) <- otu_table(apply(ps_embed_use@otu_table, 2, function(x) return( (x- mean(x)) / sd(x))), taxa_are_rows = F)
otu_table(ps_seqtab_use) <- ps_seqtab_use@otu_table[, colSums(ps_seqtab_use@otu_table) > 0]
otu_table(ps_seqtab_use) <- otu_table(apply(ps_seqtab_use@otu_table, 2, function(x) return( (x- mean(x)) / sd(x))), taxa_are_rows = F)
otu_table(ps_deseq_use) <- ps_deseq_use@otu_table[, colSums(ps_deseq_use@otu_table) > 0]
otu_table(ps_deseq_use) <- otu_table(apply(ps_deseq_use@otu_table, 2, function(x) return( (x- mean(x)) / sd(x))), taxa_are_rows = F)

preds_asinh <- list()
preds_embed <- list()
preds_pca <- list()
preds_meta <- list()

test_labels <- list()
test_labels_meta <- list()


for(fold in names(folds)){
  test_pairs <- folds[[fold]]
  
  #Raw
  print("ASIN")
  asin <- predictAut(otu = ps_seqtab_use@otu_table, mapping = ps_seqtab_use@sam_data, test_pairs,
                     include_meta = include_meta, meta_only = F)
  preds_asinh[[fold]] <- asin$preds_prob
  test_labels[[fold]] <- asin$y_test
  
  #Embed
  print("EMBED")
  embed <- predictAut(otu = ps_embed_use@otu_table,  mapping = ps_embed_use@sam_data, test_pairs,
                      include_meta = include_meta, meta_only = F)
  preds_embed[[fold]] <- embed$preds_prob
  
  ##pca
  print("PCA")
  pca <- prcomp(ps_seqtab_use@otu_table)
  data <- as.matrix(ps_seqtab_use@otu_table) %*% as.matrix(pca$rotation)
  pca <- predictAut(otu = data,  mapping = ps_seqtab_use@sam_data, test_pairs, include_meta = include_meta)
  preds_pca[[fold]] <- pca$preds_prob
  
  #meta only
  print("META")
  meta <- predictAut(otu = ps_embed_use@otu_table,  mapping = ps_embed_use@sam_data, test_pairs, include_meta = F, meta_only= T)
  preds_meta[[fold]] <- meta$preds_prob
}

plotROC <- function(preds_asinh, preds_embed, preds_pca, preds_meta, test_labels){
  
  if(!is.na(preds_embed)){
    auc_embed <- cvAUC(preds_embed, test_labels)
    red_light <-  rgb(1,0,0,alpha=0.2) 
    plot(auc_embed$perf, col = red_light)
    plot(auc_embed$perf, col = "red", avg = "vertical", add = T, lwd = 5)
    
    embed_auc_print <- mean(unlist(mapply(function(x, y) return( performance(prediction(x, y), "auc")@y.values),
                                          preds_embed, test_labels)))
  }
  if(!is.na(preds_pca)){
    auc_pca <- cvAUC(preds_pca, test_labels)
    green_light <-  rgb(0,1,0,alpha=0.2) 
    plot(auc_pca$perf, col = green_light, add = T)
    plot(auc_pca$perf, col = "green", avg = "vertical", add = T, lwd = 5)
    pca_auc_print <- mean(unlist(mapply(function(x, y) return( performance(prediction(x, y), "auc")@y.values),
                                        preds_pca, test_labels)))
  } 
  
  
  if(!is.na(preds_asinh)){
    auc_asinh <- cvAUC(preds_asinh, test_labels)
    red_light <-  rgb(1,0,0,alpha=0.2) 
    plot(auc_asinh$perf, col = red_light)
    plot(auc_asinh$perf, col = "red", avg = "vertical", add = T, lwd = 5)
    
    asinh_auc_print <- mean(unlist(mapply(function(x, y) return( performance(prediction(x, y), "auc")@y.values),
                                          preds_asinh, test_labels)))
  }

  if(!is.na(preds_meta)){
    auc_meta <- cvAUC(preds_meta, test_labels_meta)
    col_light <- rgb(0.5, 0, 0.5, alpha = 0.2)
    col_dark <- rgb(0.5, 0, 0.5, alpha = 1)
    plot(auc_meta$perf, col = col_light, add = T)
    plot(auc_meta$perf, col = col_dark, avg = "vertical", add = T, lwd = 5)
    lines(c(0,1),c(0,1),col = "black", lty = 4, lwd = 5 )
    
    meta_auc_print <- mean(unlist(mapply(function(x, y) return( performance(prediction(x, y), "auc")@y.values),
                                         preds_meta, test_labels_meta)))
  }
  
  
  legend("bottomright", legend = c(paste("Embed AUC: ", round(embed_auc_print, 3)),
                                   paste("PCA AUC: ", round(pca_auc_print, 3)),
                                   paste("Asinh AUC: ", round(asinh_auc_print, 3)),
                                   paste("Meta only AUC: ", round(meta_auc_print, 3))),
         fill = c("red", "green", "blue", "purple"))
  
}



plotROC(preds_asinh, preds_embed, preds_pca, preds_meta)


##########################
##########################
## Look at coefficients ##

printCoefficients <- function(tmp){
  data = tmp[[1]]
  y = tmp[[2]]
  cv <- cv.glmnet(x = as.matrix(data), y = as.factor(y), type.measure = "deviance", family = "binomial")
  print(coef(cv, s = cv$lambda.min))
}

tmp <- getMLInput_full(otu = ps_embed_use@otu_table, mapping = ps_embed_use@sam_data, include_meta = T)
printCoefficients(tmp)

tmp <- getMLInput_full(otu = ps_seqtab_use@otu_table, mapping = ps_seqtab_use@sam_data, include_meta = T)
printCoefficients(tmp)

tmp <- getMLInput_full(otu = ps_embed_use@otu_table, mapping = ps_embed_use@sam_data, meta_only = T)
printCoefficients(tmp)


###RANDOM FOREST##################
library(randomForest)
#Logistic regression
predictAut_rf <- function(otu, mapping, test_pairs, col = "black", include_meta= T, meta_only = F){
  tmp <- getMLInput(otu, mapping, test_pairs, include_meta = include_meta, meta_only = meta_only)
  data_train <- tmp[[1]]
  data_test <- tmp[[2]]
  y_train <- tmp[[3]]
  y_test <- tmp[[4]]
  
  ######### Use cross-validation to pick parameters from the train set ######
  control <- trainControl(method = "oob", number = 5, search = "random")
  seed <- 7
  metric = "Accuracy"
  set.seed(seed)
  
  rf_default <- train(x = data_train, y = as.factor(y_train),
                      method="rf", metric=metric, tuneLength = 50, trControl=control)
  print(rf_default$bestTune)
  rf <- randomForest(x = data_train, y = as.factor(y_train), mtry = as.numeric(rf_default$bestTune), ntree = 500)
  preds_prob <- predict(rf, newdata = as.matrix(data_test), type = "prob")[,2]
  pred <- prediction(preds_prob, y_test)
  
  # calculate probabilities for TPR/FPR for predictions
  perf <- performance(pred,"tpr","fpr")
  print(paste("AUC: ", unlist(performance(pred,"auc")@y.values))) # shows calculated AUC for model
  
  return(list(preds_prob = preds_prob, y_test = y_test, rf = rf))
}

preds_asinh <- list()
preds_embed <- list()
preds_pca <- list()
preds_meta <- list()

test_labels <- list()
test_labels_meta <- list()


for(fold in names(folds)){
  test_pairs <- folds[[fold]]
  
  #Raw
  #print("ASIN")
  #asin <- predictAut_rf(otu = ps_seqtab_use@otu_table, mapping = ps_seqtab_use@sam_data, test_pairs,
  #                   include_meta = include_meta, meta_only = F)
  #preds_asinh[[fold]] <- asin$preds_prob
  #test_labels[[fold]] <- asin$y_test
  
  #Embed
  print("EMBED")
  embed <- predictAut_rf(otu = ps_embed_use@otu_table,  mapping = ps_embed_use@sam_data, test_pairs,
                      include_meta = include_meta, meta_only = F)
  preds_embed[[fold]] <- embed$preds_prob
  test_labels[[fold]] <- embed$y_test
  
  ##pca
  print("PCA")
  pca <- prcomp(ps_seqtab_use@otu_table)
  data <- as.matrix(ps_seqtab_use@otu_table) %*% as.matrix(pca$rotation)
  pca <- predictAut_rf(otu = data,  mapping = ps_seqtab_use@sam_data, test_pairs, include_meta = include_meta)
  preds_pca[[fold]] <- pca$preds_prob
  
  #meta only
  print("META")
  meta <- predictAut_rf(otu = ps_embed_use@otu_table,  mapping = ps_embed_use@sam_data, test_pairs, include_meta = F, meta_only= T)
  preds_meta[[fold]] <- meta$preds_prob
  test_labels_meta[[fold]] <- meta$y_test
}

plotROC(preds_asinh = NA, preds_embed = preds_embed, preds_pca = preds_pca, preds_meta = preds_meta, test_labels = test_labels)


################################################
#### That's cool, but really we'd  #############
#### like to have a single model  ##############
#### for later use #############################
################################################


#roc_embed <- roc(embed$y_test, embed$preds_prob, smoothed = TRUE, plot = TRUE)
###########plot
plotROC_noCrossVal <- function(asin, embed, pca, meta){
  pred <- prediction(embed$preds_prob, embed$y_test)
  perf <- performance(pred,"tpr","fpr")
  plot(perf, col = "red")
  auc_embed <- performance(pred, measure = "auc")@y.values[[1]]
  
  
  pred <- prediction(pca$preds_prob, pca$y_test)
  perf <- performance(pred,"tpr","fpr")
  lines(perf@x.values[[1]], perf@y.values[[1]], col = "green")
  auc_pca <- performance(pred, measure = "auc")@y.values[[1]]
  
  
  pred <- prediction(asin$preds_prob, asin$y_test)
  perf <- performance(pred,"tpr","fpr")
  lines(perf@x.values[[1]], perf@y.values[[1]], col = "blue")
  auc_asin <- performance(pred, measure = "auc")@y.values[[1]]
  
  pred <- prediction(meta$preds_prob, asin$y_test)
  perf <- performance(pred,"tpr","fpr")
  lines(perf@x.values[[1]], perf@y.values[[1]], col = "purple")
  auc_meta <- performance(pred, measure = "auc")@y.values[[1]]
  
  
  legend("bottomright", legend = c(paste("Embed AUC: ", round(auc_embed, 3)),
                                   paste("PCA AUC: ", round(auc_pca, 3)),
                                   paste("DESeq AUC: ", round(auc_asin, 3)),
                                   paste("Meta only AUC: ", round(auc_meta, 3))),
         fill = c("red", "green", "blue", "purple"))
}

library(caret)
set.seed(0)
test_pairs <- sample( unique(ps_seqtab_use@sam_data$familyID), 20)


##### With metadata
include_meta = T
asin_rf_meta <- predictAut_rf(otu = ps_seqtab_use@otu_table, mapping = ps_seqtab_use@sam_data, test_pairs,
                      include_meta = include_meta, meta_only = F)

embed_rf_meta <- predictAut_rf(otu = ps_embed_use@otu_table,  mapping = ps_embed_use@sam_data, test_pairs,
                       include_meta = include_meta, meta_only = F)
test_labels <- ifelse(ps_embed_use@sam_data$phenotype == "A", 1, -1)


pca <- prcomp(ps_seqtab_use@otu_table)
data <- as.matrix(ps_seqtab_use@otu_table) %*% as.matrix(pca$rotation)
pca_rf_meta <- predictAut_rf(otu = data,  mapping = ps_seqtab_use@sam_data, test_pairs, include_meta = include_meta)

meta_rf_meta <- predictAut_rf(otu = ps_embed_use@otu_table,  mapping = ps_embed_use@sam_data, test_pairs, meta_only = T)

deseq_rf_meta <- predictAut_rf(otu = ps_deseq_use@otu_table,  mapping = ps_deseq_use@sam_data, test_pairs, include_meta = include_meta)
plotROC_noCrossVal(deseq_rf_meta, embed_rf_meta, pca_rf_meta, meta_rf_meta)



############Without metadata
include_meta = F
asin_rf <- predictAut_rf(otu = ps_seqtab_use@otu_table, mapping = ps_seqtab_use@sam_data, test_pairs,
                      include_meta = include_meta, meta_only = F)

deseq_rf <- predictAut_rf(otu = ps_deseq_use@otu_table, mapping = ps_deseq_use@sam_data, test_pairs,
                         include_meta = include_meta, meta_only = F)


embed_rf <- predictAut_rf(otu = ps_embed_use@otu_table,  mapping = ps_embed_use@sam_data, test_pairs,
                       include_meta = include_meta, meta_only = F)
test_labels <- ifelse(ps_embed_use@sam_data$phenotype == "A", 1, -1)


pca <- prcomp(ps_seqtab_use@otu_table)
data <- as.matrix(ps_seqtab_use@otu_table) %*% as.matrix(pca$rotation)
pca_rf <- predictAut_rf(otu = data,  mapping = ps_seqtab_use@sam_data, test_pairs, include_meta = include_meta)

meta_rf <- predictAut_rf(otu = ps_embed_use@otu_table,  mapping = ps_embed_use@sam_data, test_pairs, meta_only = T)
plotROC_noCrossVal(deseq_rf, embed_rf, pca_rf, meta_rf)





####Check out variable importance w/ and w/o metadata
par(mar=c(2, 4, 2, 4))
varImpPlot(meta_rf$rf)
varImpPlot(embed_rf_meta$rf)
varImpPlot(embed_rf$rf)
varImpPlot(deseq_rf_meta$rf)
varImpPlot(deseq_rf$rf)



###### Barplot for Age
library(ggplot2)
map_tmp <- ps_seqtab_use@sam_data
ggplot(map_tmp, aes(x = phenotype, y = age, fill = phenotype)) + geom_boxplot()


### Drop age and redo models
map_tmp <- map_tmp[ , !(colnames(map_tmp) %in% c("age"))]
sample_data(ps_seqtab_use) <- map_tmp
sample_data(ps_embed_use) <- map_tmp
include_meta = T
asin <- predictAut_rf(otu = ps_seqtab_use@otu_table, mapping = ps_seqtab_use@sam_data, test_pairs,
                      include_meta = include_meta, meta_only = F)

embed <- predictAut_rf(otu = ps_embed_use@otu_table,  mapping = ps_embed_use@sam_data, test_pairs,
                       include_meta = include_meta, meta_only = F)
test_labels <- ifelse(ps_embed_use@sam_data$phenotype == "A", 1, -1)


pca <- prcomp(ps_seqtab_use@otu_table)
data <- as.matrix(ps_seqtab_use@otu_table) %*% as.matrix(pca$rotation)
pca <- predictAut_rf(otu = data,  mapping = ps_seqtab_use@sam_data, test_pairs, include_meta = include_meta)

meta <- predictAut_rf(otu = ps_embed_use@otu_table,  mapping = ps_embed_use@sam_data, test_pairs, meta_only = T)
plotROC_noCrossVal(asin, embed, pca, meta)




########################## 
### glmnet ##############
#########################



####Include meta
include_meta = T
asin_meta <- predictAut(otu = ps_seqtab_use@otu_table, mapping = ps_seqtab_use@sam_data, test_pairs,
                      include_meta = include_meta, meta_only = F)

deseq_meta <- predictAut(otu = ps_deseq_use@otu_table, mapping = ps_deseq_use@sam_data, test_pairs,
                   include_meta = include_meta, meta_only = F)

embed_meta <- predictAut(otu = ps_embed_use@otu_table,  mapping = ps_embed_use@sam_data, test_pairs,
                       include_meta = include_meta, meta_only = F)
test_labels <- ifelse(ps_embed_use@sam_data$phenotype == "A", 1, -1)

pca_trans <- prcomp(ps_seqtab_use@otu_table)
data <- as.matrix(ps_seqtab_use@otu_table) %*% as.matrix(pca_trans$rotation)
pca_meta <- predictAut(otu = data,  mapping = ps_seqtab_use@sam_data, test_pairs, include_meta = include_meta)

meta_meta <- predictAut(otu = ps_embed_use@otu_table,  mapping = ps_embed_use@sam_data, test_pairs, meta_only = T)
plotROC_noCrossVal(deseq_meta, embed_meta, pca_meta, meta)

plotVariables(deseq_meta, lambda = deseq_meta$cv$lambda[30])
plotVariables(embed_meta, lambda = embed_meta$cv$lambda[30])
plotVariables(meta_meta, lambda = meta_meta$cv$lambda.min)

### Don't include meta
include_meta = F
asin <- predictAut(otu = ps_seqtab_use@otu_table, mapping = ps_seqtab_use@sam_data, test_pairs,
                   include_meta = include_meta, meta_only = F)

deseq <- predictAut(otu = ps_deseq_use@otu_table, mapping = ps_deseq_use@sam_data, test_pairs,
                    include_meta = include_meta, meta_only = F)

embed <- predictAut(otu = ps_embed_use@otu_table,  mapping = ps_embed_use@sam_data, test_pairs,
                    include_meta = include_meta, meta_only = F)
test_labels <- ifelse(ps_embed_use@sam_data$phenotype == "A", 1, -1)

pca <- prcomp(ps_seqtab_use@otu_table)
data <- as.matrix(ps_seqtab_use@otu_table) %*% as.matrix(pca$rotation)
pca <- predictAut(otu = data,  mapping = ps_seqtab_use@sam_data, test_pairs, include_meta = include_meta)

meta <- predictAut(otu = ps_embed_use@otu_table,  mapping = ps_embed_use@sam_data, test_pairs, meta_only = T)
plotROC_noCrossVal(deseq, embed, pca, meta)


#### Variables

plotVariables <- function(obj, lambda){
  par(mar=c(10, 4, 5, 2))
  obj_coef <- as.matrix(coef(obj$cv, s = lambda))
  barplot(obj_coef[abs(obj_coef) > .01 ,], col = ifelse(obj_coef[abs(obj_coef) > .01 ,] > 0, "red", "blue"), las = 2)
  
}
plotVariables(deseq, lambda = deseq$cv$lambda[30])
plotVariables(embed, lambda = embed$cv$lambda[30])
plotVariables(meta, lambda = meta$cv$lambda.min)

######### Which bugs are those ASVs?
deseq_coef <- as.matrix(coef(deseq$cv, s = deseq$cv$lambda[30]))
asvs <- names(deseq_coef[abs(deseq_coef) > .01 ,])
asvs <- asvs[grepl( "ASV", asvs)]
fasta <- read.table("repseqs.fasta")
headers <- gsub(">", "", as.character(fasta[seq(1,nrow(fasta), by = 2), 1]))
seqs <- as.character(fasta[seq(2, nrow(fasta), by = 2), 1])
library(dada2)
seqs[headers %in% asvs]
tax <- assignTaxonomy(seqs[headers %in% asvs], "silva_nr_v132_train_set.fa.gz")
species <- assignSpecies(seqs[headers %in% asvs], "rdp_species_assignment_16.fa.gz")

tax_mat <- data.frame(tax)
coef_vals <- deseq_coef[abs(deseq_coef) > .01 & grepl("ASV", rownames(deseq_coef)),]
tax_mat[ , 'coef_vals'] <- coef_vals
rownames(tax_mat) <- asvs
rownames(species) <-  asvs
write.table(tax_mat, "imp_taxa_asinh_vals.csv", sep = ",", row.names = T, col.names = T)
write.table(species, "imp_species_asinh.csv", sep = ",", row.names = T, col.names = T)


# Hypothesis: The asv model is only good as far as the asv serves as a proxy for some other metadata variable
#ASV 74 and vegetable frequency
cor.test(as.numeric(ps_seqtab_use@otu_table[ , "ASV_74"]), as.numeric(ps_seqtab_use@sam_data$vegetable))
cor.test(as.numeric(ps_seqtab_use@otu_table[ , "ASV_2778"]), as.numeric(ps_seqtab_use@sam_data$nonceliac_sensitivity))

#####################
#### Permanova ######
sample_data(ps_seqtab_1) <- ps_seqtab_1@sam_data[complete.cases(ps_seqtab_1@sam_data), ] 
dists <- distance(ps_seqtab_1, method = "bray", type = "samples")
dists <- distance(ps_embed_use, method = "euclidean", type = "samples")
#Read in qual vec distance matrix
library(vegan)
library(lsa)
#ps_use <- ps_asinh
#ps_tmp <- prune_samples(sample_names(ps_use)[!is.na(ps_use@sam_data$Milk..Cheese)], ps_use)
#dists <- rowWise_cosine_dist(ps_tmp@otu_table)
#dists <- dists + abs(min(dists))
#dists <- as.dist(dists)
#data = as.data.frame(ps_tmp@sam_data)
#data = data.frame(sample_data(ps_tmp))
map<- data.frame(ps_embed_use@sam_data)
form <- as.formula(paste("dists" , paste(colnames(ps_embed_use@sam_data), collapse = "+"), sep = " ~ "))

perm <- adonis(form, data = map)
perm

