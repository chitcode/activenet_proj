library(EBSeq)
library(igraph)
library(stringr)


genes_network <- read.csv("~/UW_classes/BMI826/project/data/protein_interactions.txt",
                          sep = "\t",stringsAsFactors = FALSE, header = FALSE)

genes_network <- genes_network[,c('V2','V5','V6','V15')]
names(genes_network) <- c("species",'protein1','protein2','score')

str_split <- function(s,sep='-'){
  splits = strsplit(s,sep)[[1]]
  if(length(splits) >= 2)
    return(splits[2])
  else
    return(s)
}

genes_network$species = sapply(genes_network$species, function(x) str_split(x,'_'))
genes_network$protein1 = sapply(genes_network$protein1, function(x) str_split(x,':'))
genes_network$protein2 = sapply(genes_network$protein2, function(x) str_split(x,':'))
genes_network$score = sapply(genes_network$score, function(x) str_split(x,':'))
genes_network$score <- as.numeric(genes_network$score)


genes_network <- genes_network[genes_network$species == 'HUMAN',c("protein1","protein2","score")]

protein.graph <- graph.data.frame(genes_network,directed = FALSE)
vcount(protein.graph)


lungbrain_data <- read.csv("~/UW_classes/BMI826/project/data/lungbrain_data_normed.csv", stringsAsFactors = FALSE)

row.names(lungbrain_data) <- lungbrain_data[,1]
lungbrain_data <- lungbrain_data[,2:dim(lungbrain_data)[2]]

conditions = sapply(names(lungbrain_data), function(x) str_sub(x,-2,-1))
lt_bt <- conditions == 'LT' |  conditions == 'BT'

data_lt_bt <- lungbrain_data[,lt_bt]

EBOut <- EBTest(Data=as.matrix(data_lt_bt), Conditions=conditions[lt_bt], 
                sizeFactors=rep(1,dim(data_lt_bt)[2]), maxround=5)


EBDERes = GetDEResults(EBOut, FDR=0.05)
dim(EBDERes$PPMat)

p_values <- EBOut$PPDE


# Cleaning the network (removing the nodes not present in the RNA sequencing data)
genes_remove <- setdiff(genes_network$protein1,names(p_values))
G1 <- delete_vertices(protein.graph,genes_remove)

p_mean <- mean(p_values)
p_sd <- sd(p_values)
z_scores <- (p_values-p_mean)/p_sd

# Building random Z score for corresponding k sized network
# let's build k upto 100 using Monte-Carlo simulation

max_k <- length(z_scores)
mu_k <- numeric(max_k)
sigma_k <- numeric(max_k)
for (k in 1:max_k){
  sim_count <- 1000
  mc_z <- numeric(sim_count)
  for (i in 1:sim_count){
    mc_z[i] <- sum(sample(z_scores,size = k,replace = TRUE))/sqrt(k)
  }
  mu_k[k] <- mean(mc_z)
  sigma_k[k] <- sd(mc_z)
  
  if(k %% 200 == 0){
    cat(".")
  }
}


score_subgrpup <- function(group_genes){
  #getting the z-score for the genes
  k <- length(group_genes)
  z_A = sum(group_genes)/sqrt(k)
  
  s_A <- (z_A - mu_k[k])/sigma_k[k]
  return(s_A)
}

simulated_annealing <- function(G,N,score_func,Temp_high,cool){
  nodes_count <- length(V(G))
  nodes_active <- rep(1,nodes_count)#1* (runif(nodes_count) >= 0.5)
  subnet_score <- score_func(z_scores[nodes_active])
  temperature = Temp_high
  
  scores = numeric(N)
  for(count in 1:N){
    #what's the new temperature
    temperature <- temperature * cool
    if(temperature > 0.001){
      node_selected <- sample(1:nodes_count,1)
      changed_state <- nodes_active
      changed_state[node_selected] <- 1* (!changed_state[node_selected])#fliping the node state
      
      changed_subnet_score <- score_func(z_scores[changed_state])
      print(sum(nodes_active))
      print(subnet_score)
      print(sum(changed_state))
      print(changed_subnet_score)
      p = exp((changed_subnet_score-subnet_score)/temperature)
      if((changed_subnet_score > subnet_score) | (runif(1) < p)){
        nodes_active <- changed_state
        subnet_score <- changed_subnet_score
      }
      
      scores[count] = subnet_score
    }else{
      break()
    }
    if(count %% 100 ==0){
      cat(".")
    }
  }
  return(list(nodes_active,scores))
}

result <- simulated_annealing(G1,100,score_subgrpup,1,0.99)

