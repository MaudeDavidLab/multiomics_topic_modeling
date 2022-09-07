library(tm)
library(slam)
library(topicmodels)
library(reshape2)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
input_file = as.character(args[1])
output_file = as.character(args[2])
filter_thresh = as.numeric(args[3])
num_topics = as.numeric(args[4])
iter=as.numeric(args[5])

print(paste0("input file: " , input_file))
print(paste0("output file: " , output_file))
print(paste0("filter threshold: " , filter_thresh))
print(paste0("number of topics:"  , num_topics))
print(paste0("number of iterations: " , iter))

seed_num = 2014

getLDAInput <- function(ps, filter_thresh = 0.05){
  data.final <- t(round(ps@otu_table))
  print(dim(data.final))
  print(filter_thresh)
	#Filter out less common ASVs
  keep <- apply(data.final, 2, function(x) return(sum(x>0) > filter_thresh*nrow(data.final)))
  print(sum(keep))
  data.final <- data.final[ , keep]

  dtm=as.simple_triplet_matrix(data.final)
  return(dtm)
}

#Read in data 
ps <- readRDS(input_file)
print("making input")
dtm <- getLDAInput(ps, filter_thresh = filter_thresh)
print("finished making input")

# label using molecule names
matched <- match(colnames(dtm), taxa_names(ps))
colnames(dtm) == taxa_names(ps)[matched]
colnames(dtm) <- as.character(ps@tax_table[matched, 1])
print("Added labels")

#Train model
print("Beginning training...")
print(num_topics)
model <- LDA(dtm, k = num_topics, method = "Gibbs", control = list(seed = seed_num, burnin = 1000, thin = 100, iter = iter, verbose=10))
print("Finished training...")
saveRDS(model, output_file)



#####################################
############# 16s ###################
#####################################

# Read in data 
#ps_16s <- readRDS("ps_16s_dds.rds")
#taxa_names(ps_16s) <- ps_16s@tax_table[,2]

# get full taxa table
#ps_not_norm_comp <- readRDS("ps_DeSeq_pass_min_postDD_min0.03.rds")
#tax_table(ps_16s) <- ps_not_norm_comp@tax_table


# Prepare input 
#dtm <- getLDAInput(ps_16s, filter_thresh = 0.05)
#colnames(dtm) <- paste(data.frame(ps_16s@tax_table)[colnames(dtm), "Phylum"], 
#                       data.frame(ps_16s@tax_table)[colnames(dtm), "Family"],
#                       data.frame(ps_16s@tax_table)[colnames(dtm), "Species"])
# Train model
#model_16s = LDA(dtm, k = 3, method = "Gibbs", control = list(seed = seed_num, burnin = 1000, thin = 100, iter = 1000))
#saveRDS(model_16s, "model_gibbs_16s_3topics.rds")


####################################
##########  MTG  ###################
####################################

# Read in data
#ps_mtg <- readRDS("ps_mtg_rle_nooutliers.rds")
#print("making input")
#dtm <- getLDAInput(ps_mtg, filter_thresh = 0.1)
#print("finished making input")

# label using molecule names
#matched <- match(colnames(dtm), taxa_names(ps_mtg))
#colnames(dtm) == taxa_names(ps_mtg)[matched]
#colnames(dtm) <- as.character(ps_mtg@tax_table[matched, 1])
#print("Added labels")

#Train model 
#print("Beginning training...")
#model_mtg <- LDA(dtm, k = 4, method = "Gibbs", control = list(seed = seed_num, burnin = 1000, thin = 100, iter = 1000))
#print("Finished training...")
#saveRDS(model_mtg, "model_gibbs_mtg_4topics_2.rds")




