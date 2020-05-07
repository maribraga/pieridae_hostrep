
library(picante)
library(ggtree)
library(tidyverse)
library(patchwork)
library(viridis)
library(bipartite)
library(ggraph)
library(tidygraph)
#library(igraph)

# library(data.tree)
# ?Traverse

load("./inference/char_hist.RData")


ggtree(tree) + geom_tiplab(size=2) + geom_nodelab(size=2) + theme_tree() + xlim(c(0,110))

# To match node/tip names and numbers
labels <- tibble(label = c(rev(tree$tip.label), paste0("Index_",67:131)), number = 1:131)

# Get a list of realized (high pp) host repertoires per age
pt <- 90 
# gather networks in a list
list_m_at_ages_50h <- list_m_at_ages
list_m_at_ages_50h[[9]] <- ext_net_50h

rep50_list <- list()
for(m in 1:length(list_m_at_ages_50h)){
  matrix <- list_m_at_ages_50h[[m]]
  for(i in 1:nrow(matrix)){
    for(j in 1:ncol(matrix)){
      if(matrix[i,j] < pt/100){
        matrix[i,j] = 0
      } else{
        matrix[i,j] = 1
      }
    }
  }
  df <- as.data.frame(matrix) # so that it doesn't become a vector
  rep50_list[[m]] <- df
}

# Get first and last age at which each node and tip is present 
duration <- list() 
for(i in 1:nrow(labels)){
  node <- labels %>% filter(number == i) %>% pull(label)
  times <- c()
  
  for(a in 1:length(ages)){
    if(node %in% colnames(rep50_list[[a]])){
      times <- c(times,ages[a])
    }
  }
  duration[[i]] <- times
}

duration[[130]] <- duration[[131]] <- NA


# Get the parent of each node
parents <- history_dat_bl1 %>%
  select(node_index,parent_index) %>% 
  distinct() %>% 
  arrange(node_index)

# List of vectors with the sequence of nodes from tip to root for each tip
traverse <- list()
for(t in 1:Ntip(tree)){
  n = t
  lineage <- c(n)
  while(!is.na(n)){
    parent <- parents %>% filter(node_index == n) %>% pull(parent_index)
    lineage <- c(lineage, parent)
    n = parent
  }
  traverse[[t]] <- lineage
}


# Make a list with repertoire by age for each lineage (tip to root)
rep_by_lineage <- list()

for(t in 1:Ntip(tree)){
  
  rep_by_age <- list()
  lineage <- traverse[[t]]
  times <- duration[[t]]
  
  lin_pos <- 1
  n <- lineage[lin_pos]
  label <- labels %>% filter(number == n) %>% pull(label) 
  
  for(a in length(ages):1){
    
    while(!(label %in% colnames(rep50_list[[a]]))){
      lin_pos <- lin_pos + 1
      n <- lineage[lin_pos]
      label <- labels %>% filter(number == n) %>% pull(label)
    }
    
    rep_by_age[[a]] <- select(rep50_list[[a]], any_of(label))
  }
  rba <- bind_cols(rep_by_age)
  colnames(rba) <- ages
  rep_by_lineage[[t]] <- rba
}

View(rep_by_lineage[[10]])

# It works, but most of the columns (repertoires) are repeated, because of the tree structure 






