library(KEGGREST)
library(Matrix)
library(aod)

genes <- keggList("ko")
compounds <- keggList("cpd")
gene_names <- names(genes)
compound_names <- names(compounds)
compound_names <- gsub("cpd:", "", compound_names)



getReactionsFromGene <- function(gene_name){
  try({
    gene <- keggGet(gene_name)
    reaction_string <- gene[[1]]$DBLINKS[1]
    reaction_string <- gsub("RN: ", "", reaction_string)
    reactions <- strsplit(reaction_string, " ")[[1]]
    return(reactions)
  })
}

getReactionsFromMetabolite <- function(metabolite_name){
  try({ 
     metabolite <- keggGet(metabolite_name)
     reactions <- metabolite[[1]]$REACTION[[1]]
     return(reactions)
  })
}

getMetabolitesFromReaction <- function(reaction){
  try({
    tmp <- keggGet(reaction)
    reaction_equation <- tmp[[1]]$EQUATION
    reaction_equation <- gsub("<=>", "+", reaction_equation)
    reaction_equation <- gsub(" ", "", reaction_equation)
    metabolites <- strsplit(reaction_equation, "\\+")[[1]]
    return(metabolites)
  })
}

getGenesFromReaction <- function(reaction){
  try({
    tmp <- keggGet(reaction)
    genes <- names(tmp[[1]]$ORTHOLOGY)
    return(genes)
  })
}

possible_transitions <- Matrix(matrix(rep(0, length(compound_names) * length(gene_names)),
                                      ncol = length(compound_names)), sparse = T)
rownames(possible_transitions) <- gsub("ko:", "", gene_names)
colnames(possible_transitions) <- compound_names

print(length(gene_names))
i=1
for(gene_name in gene_names[721:length(gene_names)]){
  print(i, end = "\t")
  reactions <- getReactionsFromGene(gene_name)
  metabolite_vec <- sapply(reactions, function(reaction){
    metabolites <- getMetabolitesFromReaction(reaction)
    return(metabolites)
  })
  metabolite_vec <- as.vector(unlist(metabolite_vec))
  if(!(TRUE %in% grepl("Error", metabolite_vec))){
    #check for the presence of an error message
    metabolite_vec <- paste("C", sub(".*?C", "", metabolite_vec), sep = "")
    metabolite_vec <-sub("\\(.*\\)", "", metabolite_vec)
    metabolite_vec <- unique(metabolite_vec)
    metabolite_vec <- metabolite_vec[metabolite_vec %in% compound_names]
    if(length(metabolite_vec) > 1){
      indicies <- match(metabolite_vec, compound_names)
      possible_transitions[i,indicies] <- 1
    }
    
  }else{
    print("error, moving on")
  }
  i = i + 1
}
saveRDS(possible_transitions, "possible_transitions_genes_compouns.rds")


# Do wald test to get ranked list of metabolites
setwd("C:/Users/kitikomp/Documents/Lab/M3")
ps_metabol <- readRDS("data/metabol/ps_rle_nooutliers_fixedmapping.rds")
df <- as.data.frame(t(ps_metabol@otu_table))
df$phenotype <- ifelse(ps_metabol@sam_data$phenotype == "A", 1, 0)
colnames(df) <- gsub("-", ".", colnames(df))
colnames(df) <- gsub(",", ".", colnames(df))
colnames(df) <- gsub("\\(", ".", colnames(df))
colnames(df) <- gsub("\\)", ".", colnames(df))
colnames(df) <- gsub(" ", ".", colnames(df))
colnames(df) <- gsub("\\*", ".", colnames(df))
colnames(df) <- gsub("\\^", ".", colnames(df))
colnames(df) <- gsub("\\+", ".", colnames(df))
colnames(df) <- gsub("\\/", ".", colnames(df))
colnames(df) <- gsub("\\[", ".", colnames(df))
colnames(df) <- gsub("\\:", ".", colnames(df))
colnames(df) <- gsub("\\;", ".", colnames(df))
colnames(df) <- gsub("\\]", ".", colnames(df))

colnames(df) <- gsub('[[:digit:]]+', '', colnames(df))
stats <- sapply(colnames(df), function(x){
  f <- formula(paste("phenotype ~ ", x))
  line <- glm(f, data = df, family = "binomial")
  res <- wald.test(b = coef(line), Sigma = vcov(line), Terms = 1)
  stat <- res$result$chi2[1]
  return(stat)
})


ord <- order(stats, decreasing = T)
df <- data.frame(taxa_names(ps_metabol)[ord],
                 stat = stats[ord],
                 kegg_id = data.frame(tax_table(ps_metabol))[ord, "KEGG"])


df <- df[df$kegg_id!= "", ]
df$kegg_id <- sub(".+?C", "", df$kegg_id)
df$kegg_id <- sub(",.*", "", df$kegg_id)
df <- df[grepl("C", df$kegg_id), ]
df <- df[!duplicated(df$kegg_id), ]
rownames(df) <- df$kegg_id

possible_transitions <- readRDS("results/page_rank/possible_transitions_genes_compouns.rds")
zeros <- Matrix(matrix(rep(0, ncol(possible_transitions)^2), ncol = ncol(possible_transitions)), sparse = T)
zeros_2 <- Matrix(matrix(rep(0, nrow(possible_transitions)^2), ncol = nrow(possible_transitions)), sparse = T)

possible_transitions_t <- t(possible_transitions)

tmp1 <- rbind(possible_transitions, zeros)
tmp2 <- rbind(zeros_2, possible_transitions_t)
transitions <- cbind(tmp1, tmp2)





