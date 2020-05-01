#'---
#'title: "Pieridae host repertoire - character history"
#'author: "Mariana Braga"
#'date: "`r format(Sys.time(), '%d %B, %Y')`"
#'output: github_document
#'---

#'-------------
#'
#' Script 2 for empirical study performed in Braga et al. 2020
#' *Evolution of butterfly-plant networks revealed by Bayesian inference of host repertoire*.

#+ include = FALSE

library(ape)
library(phytools) # needed?
library(dispRity)
library(ggtree)
library(picante)

library(tidyverse)
library(patchwork)
library(wesanderson)

library(MCMCpack)
library(coda)
library(kdensity)
#library(phylotate) # needed?
#library(data.table) # needed in script sourced

library(bipartite)
library(ggraph)
library(tidygraph)
library(igraph)
#library(ndtv)        # needed?
#library(intergraph)  # needed?

#' ### Data
#' First we read in the phylogenetic trees for butterflies and plants.
#' Then, we read in the interaction matrix and remove plants (rows)
#' that are not hosts to any butterfly.
#'
/*#### INPUT DATA ####

# _Read trees ----
*/
#' **Trees**
#' 
tree <- read.tree("./data/bphy_pie_ladder.phy")
host_tree <- read.tree("./data/angio_pie_50tips_ladder.phy")

#+ include = FALSE
## for figure 2
# host_tree_bl1 <- read.tree("./data/angio_pie_50tips_bl1.phy")
# 
# ggtree(host_tree_bl1) + 
#   theme_tree2() + 
#   geom_tiplab(align=TRUE, linesize=.5) + 
#   xlim(0, 25)


#' **Extant network**
/*# _Extant network ----
*/
#+ results='hide'
ext_net_50h <- as.matrix(read.csv("./data/incidence_pieridae.csv", header = T, row.names = 1))
identical(colnames(ext_net_50h), tree$tip.label)
identical(rownames(ext_net_50h), host_tree$tip.label)

#+
ext_net <- ext_net_50h[which(rowSums(ext_net_50h) != 0),]
dim(ext_net)
  
#' ### Character history
#' 

/*##### CHARACTER HISTORY #####
  
# Read in .history.txt file ----
*/

#' **Read in .history.txt files**
#'
#' These files can get quite big, so I can't upload them in Github.
#' Also, you might want to thin out these files to speed up their parsing. 
#' This is how you can do it.
#' 
#' (Note that in the analyses in the paper we did not thin the histories, 
#' so we'll show results with all sampled histories).

colclasses <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))
#+ eval = FALSE
history_dat_time = read.table("./inference/out.2.real.pieridae.2s.history.txt", sep="\t", header=T, colClasses = colclasses)
 
# define burnin and sampling interval  
it_seq <- seq(20000,200000, 100)  

# Time-calibrated tree
history_dat_time <- filter(history_dat_time, iteration %in% it_seq) %>% 
  mutate(node_index = node_index + 1)

write.table(history_dat_time,"./inference/history.time.txt", sep="\t", quote = F, row.names = F)

# Tree with branch lengths = 1
history_dat_bl1 = read.table("./inference/out.3.bl1.pieridae.2s.history.txt", sep="\t", header=T, colClasses = colclasses)

history_dat_bl1 <- filter(history_dat_bl1, iteration %in% it_seq) %>% 
  mutate(node_index = node_index + 1)

write.table(history_dat_bl1,"./inference/history.bl1.txt", sep="\t", quote = F, row.names = F)

#' The next time you run this script you just need to read in the thinned histories.
#' One file with sampled histories when using the time-cklibrated host tree,
#' and one using the transformed host tree (all branch lengths = 1).
#' 

#+
history_dat_time = read.table("./inference/history.time.txt", sep="\t", header=T, colClasses = colclasses)
history_dat_bl1 = read.table("./inference/history.bl1.txt", sep="\t", header=T, colClasses = colclasses)


#' **Calculate effective rate of evolution**
#' 
/*# Calculate effective rate of evolution ----
*/
  
tree_length <- sum(tree$edge.length)+86
n_events_time <- group_by(history_dat_time,iteration) %>% 
  summarise(n = n()) %>% 
  summarise(mean = mean(n)) %>% 
  pull(mean)
n_events_bl1 <- group_by(history_dat_bl1,iteration) %>% 
  summarise(n = n()) %>% 
  summarise(mean = mean(n)) %>% 
  pull(mean)

(rate_time <- n_events_time/tree_length)
(rate_bl1 <- n_events_bl1/tree_length)

#' In both analyses, the rate of host repertoire evolution was near 1 event every 10 Myr (along each branch).
#'

/*
## (States at nodes) ----

#+ eval = FALSE
source("functions_ancestral_states.R")

hosts <- host_tree$tip.label
nodes <- 67:131

pp <- make_matrix_nodes(history_dat, nodes, c(2))
row.names(pp) <- paste0("Index_",nodes)
colnames(pp) <- host_tree$tip.label

assign(paste0("pp_",name), pp)

graph <- graph_from_incidence_matrix(pp, weighted = TRUE)
el <- get.data.frame(graph, what = "edges") %>% 
  mutate(to = factor(to, levels = hosts),
         from = factor(from, levels = paste0("Index_",nodes))) %>% 
  rename(p = weight)

#' Read probability matrix
#' 

#el <- readRDS("./inference/states_at_nodes_bl1.rds")
el <- readRDS("./inference/states_at_nodes_time.rds")

# all interactions
ggplot(el, aes(x = to, y = from)) +
  geom_tile(aes(fill = p)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_gradient(low = "white", high = "black") +
  labs(fill = "Posterior\nprobability") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())

# only high probability
filter(el, p > 0.9) %>% 
  ggplot(aes(x = to, y = from)) + 
  geom_tile(aes(fill = p)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_gradient(low = "grey50", high = "black") +
  labs(fill = "Posterior\nprobability") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())

*/
 

/*## States at ages ----

# _Get posteriors at ages ----
*/

#' **Ancestral networks**
#'
#' First load the script `functions_ancestral_states.R` which contains the functions 
#' to calculate the posterior probabilities for each plant-butterfly 
#' interaction at given times in the past.
#' 
#+ eval = FALSE
source("functions_ancestral_states.R")

#' Then choose the times in the past at which you want to reconstruct the network.
#' 
ages = c(80,70,60,50,40,30,20,10,0)

#' This part is slow, so I won't run it here.
#' 
#+ eval = FALSE
list_m_at_ages = list()
for (i in 1:(length(ages)-1)) {
  age = ages[i]
  list_m_at_ages[[i]] = t(make_matrix_at_age( history_dat_time, age, s_hit=c(2) ))
}

#' I have done this before and saved the list as a .rds file.
#' 
list_m_at_ages <- readRDS("./inference/list_m_at_ages_time.rds")
  
length(list_m_at_ages)
head(list_m_at_ages[[1]])

#' Then we add the extant network in the list with all ancestral networks.
#'
list_m_at_ages[[9]] <- ext_net

#' We can build two kinds of ancestral networks, binary and quantitative.
#' For binary networks, we need to choose a probability threshold above which 
#' we consider the interaction has enough support. Interactions above the threshold 
#' are coded as 1, and the remaining interactions are coded as 0.
#'   

#' **Binary networks with a probability threshold**
#' 
/*# _(Binary networks with a probability threshold) ----
*/

#' After defining the probability threshold, we transform the probabilities into 1s and 0s.
#'    
# probability threshold
pt <- 90 
# gather networks in a list
net_list <- list()
for(m in 1:length(list_m_at_ages)){
  matrix <- list_m_at_ages[[m]]
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
  df = df[ rowSums(df)!=0, ]  # if only one row/column is left
  df = df[ ,colSums(df)!=0 ]
  net_list[[m]] <- df
}

net_list[[1]]


/*# __Calculate modularity with bipartite (STOCHASTIC STEP!)----
*/
#' *Calculate modularity (stochastic step!)*
#' 

#' This doesn't work for the network at 80Ma because it only has one host, 
#' so we add it manually before calculatig for the other networks.
#' 

all_mod <- tibble(name = c(rownames(net_list[[1]]),colnames(net_list[[1]])), 
                  age = 80,
                  original_module = 1)

for(i in 2:length(net_list)){
  set.seed(5)
  mod <- computeModules(net_list[[i]])
  assign(paste0("mod_",ages[i]),mod)
  mod_list <- listModuleInformation(mod)[[2]]
  nmod <- length(mod_list)
  
  for(m in 1:nmod){
    members <- unlist(mod_list[[m]])
    mtbl <- tibble(name = members, 
                   age = rep(ages[i], length(members)),
                   original_module = rep(m, length(members)))
    
    all_mod <- bind_rows(all_mod, mtbl)
  }
}

# check modules for some network
plotModuleWeb(mod_60, labsize = 0.4)


  
# __(Make tidygraphs with modules) ----

minp <- 10

list_tgraphs <- list()
for(n in 1:length(net_list)){
  # get names from nets
  net <- net_list[[n]]
  rnames <- rownames(net)
  cnames <- colnames(net)
  
  # get weighted graph from list_m_at_ages - after setting minimum weight
  wnet <- list_m_at_ages[[n]]
  for(i in 1:nrow(wnet)){
    for(j in 1:ncol(wnet)){
      if(wnet[i,j] < minp/100){
        wnet[i,j] = 0
      }
    }
  }
  wnet = wnet[ rowSums(wnet)!=0, ]
  wnet = wnet[ ,colSums(wnet)!=0 ]
  
  wgraph <- as_tbl_graph(wnet, directed = F) %>% 
    left_join(filter(all_mod, age == ages[n]), by = "name") %>% 
    select(type, name, original_module)
  
  #assign(paste0("wgraph_",ages[n]), wgraph)
  
  list_tgraphs[[n]] <- wgraph
}

/*
  
# __Match modules across ages ----

# if fixed on excel

all_mod_edited <- read.csv(paste0(path_bip,"quant_modules_bipartite_seed5.csv"), header = T, stringsAsFactors = F)
for(n in 1:length(ages)){
  list_tgraphs[[n]] <- list_tgraphs[[n]] %>% 
    activate(what = "nodes") %>% 
    left_join(filter(all_wmod_edited, age == ages[n]) %>% select(name, module), by = "name")
}





  
/*# _Weighted networks with a low probability threshold ----
*/

#' **Weighted networks**
#' 
#' To build weighted networks we use the posterior probabilities as weights for each interaction.
#' But since many interactions have really small probabilities, we can set a minimum probability,
#' below which the weight is set to 0.
#' 

# probability threshold
lpt <- 10 
# gather networks in a list
list_wnets <- list()
for(m in 1:length(list_m_at_ages)){
  matrix <- list_m_at_ages[[m]]
  for(i in 1:nrow(matrix)){
    for(j in 1:ncol(matrix)){
      if(matrix[i,j] < lpt/100){
        matrix[i,j] = 0
      }
    }
  }
  matrix = matrix[ rowSums(matrix)!=0, ]
  matrix = matrix[ ,colSums(matrix)!=0 ]
  list_wnets[[m]] <- matrix
  #assign(paste0("net_",pt,"p_",ages[m]), matrix)
}

list_wnets[[1]]


# __Calculate weighted modularity with bipartite (STOCHASTIC STEP!)----

all_wmod <- tibble()
#set.seed(5) # set.seed(5) produces 10 modules

for(i in 1:length(list_wnets)){
  set.seed(5)
  wmod <- computeModules(list_wnets[[i]])
  assign(paste0("wmod_",ages[i]),wmod)
  wmod_list <- listModuleInformation(wmod)[[2]]
  nwmod <- length(wmod_list)
  
  for(m in 1:nwmod){
    members <- unlist(wmod_list[[m]])
    mtbl <- tibble(name = members, 
                   age = rep(ages[i], length(members)),
                   original_module = rep(m, length(members)))
                   
    all_wmod <- bind_rows(all_wmod, mtbl)
  }
}

# check modules for last wnet
plotModuleWeb(wmod, labsize = 0.4)

#path_bip <- paste0("./net_structure/bipartite/")
#write.csv(all_wmod, paste0(path_bip,"quant_modules_bipartite_seed5.csv"), row.names = F)


# __Make tidygraphs with bipartite modules ----

list_tgraphs <- list()
for(n in 1:length(list_wnets)){
  
  # get weighted graph
  wnet <- list_wnets[[n]]
  # get names
  rnames <- rownames(wnet)
  cnames <- colnames(wnet)
  
  wgraph <- as_tbl_graph(wnet, directed = F) %>% 
    left_join(filter(all_wmod, age == ages[n]), by = "name") %>% 
    select(type, name, original_module)
  
  #assign(paste0("wgraph_",ages[n]), wgraph)
  
  list_tgraphs[[n]] <- wgraph
}

#save.image("~/Box/Mari+Michael/network_evol/pieridae/graphs and trees.RData")

# __Match modules across ages ----

# if fixed on excel

all_wmod_edited <- read.csv(paste0(path_bip,"quant_modules_bipartite_seed5.csv"), header = T, stringsAsFactors = F)
for(n in 1:length(ages)){
  list_tgraphs[[n]] <- list_tgraphs[[n]] %>% 
    activate(what = "nodes") %>% 
    left_join(filter(all_wmod_edited, age == ages[n]) %>% select(name, module), by = "name")
}

# # try to fix it here
# 
# (degree_mod <- list_tgraphs[[9]] %>% 
#    activate(what = "nodes") %>% 
#    mutate(degree = centrality_degree()) %>% 
#    as_tibble() %>%
#    group_by(type, original_module) %>% 
#    top_n(1, degree) %>% 
#    group_by(original_module) %>% 
#    top_n(1, degree) %>% 
#    arrange(original_module))
# 
# list_modules_order <- list(
#   c('M1','M2'),
#   c('M1','M2'),
#   c('M9','M2','M3','M1'), # M9 = Pseudopontia
#   c('M2','M9','M3','M1'),
#   c('M4','M9','M2','M1','M3'),
#   c('M2','M7','M3','M4','M1','M9'), # M7 = Nepheronia
#   c('M5','M7','M4','M3','M2','M1','M6'), # Nepheronia is in M2 and Pseudopontia is in M4; M5 = Asteraceae, M6 = Rhamnaceae, M7 = Kricogonia
#   c('M2','M4','M6','M1','M5','M8','M7','M3') # new: M8 = Aporia
# )
# 
# 
# all_wmod$module <- NA
# 
# for(i in 1:nrow(all_wmod)){
#   a <- which(ages == all_wmod$age[i])
#   m <- all_wmod$original_module[i]
#   all_wmod$module[i] <- list_modules_order[[a]][m]
# }
# 
# View(all_wmod)
    

# __Make ggtree ----
# Read modified .tre file to read index as node labels
tree <- treeio::read.newick(paste0(path_data,"tree_nodelab.tre"), node.label = "label")
# Add its root time
tree$root.time <- max(tree.age(tree)$ages)


## Slice the tree at ages and create data frame with module info

list_subtrees <- list()
list_tip_data <- list()

for(i in 1:(length(ages)-1)){
  subtree <- slice.tree(tree, age = ages[[i]], "acctran")
  list_subtrees[[i]] <- subtree
  #assign(paste0("tree_",i), subtree)
  
  graph <- list_tgraphs[[i]]
  mod_from_graph <- tibble(module = activate(graph,nodes) %>% filter(type == TRUE) %>% pull(module),
                           label = activate(graph,nodes) %>% filter(type == TRUE) %>% pull(name))
  # extra step just to check that tip labels and graph node names match
  tip_data <- tibble(label = subtree$tip.label) %>% 
    inner_join(mod_from_graph) 
  list_tip_data[[i]] <- tip_data
  #assign(paste0("tip_data_",i), tip_data)
}
list_subtrees[[9]] <- tree
list_tip_data[[9]] <- tibble(label = tree$tip.label) %>% 
  inner_join(filter(all_wmod_edited, age == 0), by = c("label" = "name"))

# __Plot ggtree and ggraph ----


# brewer.pal(12, "Paired")

mod_levels <- c(paste0('M',1:12),paste0('T',1:5))
custom_pal <- c("#b4356c","#1b1581","#e34c5b","#fca33a","#fbeba9","#fdc486",
                "#802b72","#f8c4cc","#c8d9ee","#82a0be","#00a2bf","#006e82", 
                paste0("grey",c(10,30,50,70,90)))
tip_size = c(5.5,5,5,4,3,2,2,2,2)
node_size = c(5,5,4,3,3,3,3,3,3)

for(i in 1:length(ages)){
  
  subtree <- list_subtrees[[i]]
  ggt <- ggtree(subtree, ladderize = F) %<+% list_tip_data[[i]] +
    geom_tippoint(aes(color = factor(module, levels = mod_levels)), size = tip_size[i]) + 
    #geom_tiplab(size = 1) +
    geom_rootedge(rootedge = 1) +
    scale_color_manual(values = custom_pal,na.value = "grey70", drop = F) +
    #scale_color_viridis(discrete = T, na.value = "grey70", option = "plasma", end = 0.9) +
    xlim(c(0,tree$root.time)) +
    #labs(color = "")
    theme(legend.position = "none")
  
  assign(paste0("ggt_",ages[[i]]), ggt)
  
  graph <- list_tgraphs[[i]]
  laybip = layout_as_bipartite(graph)
  laybip = laybip[,c(2,1)]
  
  #ggn <- 
    
    ggraph(graph, layout = laybip) +
    geom_edge_link(aes(width = weight), color = "grey50") + 
    geom_node_point(aes(shape = type, color = factor(module, levels = mod_levels)), size = 2) +
    scale_shape_manual(values = c("square","circle")) +
    #scale_edge_color_grey(start = 0.8, end = 0.5) +
    scale_color_manual(values = custom_pal, na.value = "grey70", drop = F) +
    scale_edge_width("Probability", range = c(0.1,1)) +
    labs(title = paste0(ages[[i]]," Ma"), shape = "", color = "Module") +    # CHANGE
    theme_void() #+
    #theme(legend.position = "none")

  
  # ggn <- ggraph(graph, layout = "stress") +
  #   geom_edge_link(aes(width = weight), color = "grey50") +  # , color = weight
  #   geom_node_point(aes(shape = type, color = factor(module, levels = mod_levels)), size = node_size[i]) + 
  #   scale_shape_manual(values = c("square","circle")) +
  #   #scale_edge_color_grey(start = 0.8, end = 0.5) +
  #   scale_color_manual(values = custom_pal, na.value = "grey70", drop = F) +
  #   #scale_color_viridis(discrete = T, na.value = "grey70", option = "plasma", end = 0.9) +
  #   scale_edge_width("Probability", range = c(0.3,1)) +
  #   labs(title = paste0(ages[[i]]," Ma"), shape = "", color = "Module") +    # CHANGE
  #   theme_void() + 
  #   theme(legend.position = "none") #(legend.margin = margin(0,1,0,1,"cm"))
  
  assign(paste0("ggn_",ages[[i]]), ggn)
}

# plot subtree and network per age
ggt_20 + ggn_20

# plot all subtrees
ggt_80 + ggt_70 + ggt_60 + 
  ggt_50 + ggt_40 + ggt_30 + 
  ggt_20 + ggt_10 + ggt_0 +
  plot_layout(nrow = 3,ncol = 3, heights = c(1,2,3), guides = 'keep')

# plot all networks
(ggn_80 | ggn_70 | ggn_60) /
(ggn_50 | ggn_40 ) 

(ggn_30 | ggn_20) 

(ggn_10 | ggn_0) +
  plot_layout(widths = c(2,3))


# ___Figure 3 ----

(ga <- ggt_80 + ggn_80 + plot_layout(widths = c(2,3)))

(gb <- ggt_70 + ggn_70 + plot_layout(widths = c(2,3)))

(gc <- ggt_60 + ggn_60 + plot_layout(widths = c(2,3)))

(gd <- ggt_50 + ggn_50 + plot_layout(widths = c(2,3)))

(ge <- ggt_40 + ggn_40 + plot_layout(widths = c(2,3)))

(gf <- ggt_30 + ggn_30 + plot_layout(widths = c(2,3)))

(gg <- ggt_20 + ggn_20 + plot_layout(widths = c(2,3)))

(gh <- ggt_10 + ggn_10 + plot_layout(widths = c(2,3)))

(gi <- ggt_0 + ggn_0 + plot_layout(widths = c(2,3)))

(ga / gb / gc / gd) |
  (ge / gf / gg) |
  (gh / gi) + 
  plot_layout(width = c(1,2,3))


# ___Figure 2 ----

edge_list <- get.data.frame(list_tgraphs[[9]], what = "edges") %>%
  inner_join(all_wmod_edited %>% filter(age == 0) %>% select(name, module), by = c("from" = "name")) %>%
  inner_join(all_wmod_edited %>% filter(age == 0) %>% select(name, module), by = c("to" = "name")) %>%
  mutate(Module = ifelse(module.x == module.y, module.x, NA)) #%>% factor())

phylob <- tree$tip.label
phylop <- host_tree$tip.label

plot_net <- edge_list %>% mutate(
  to = factor(to, levels = phylob),
  from = factor(from, levels = phylop))

ggplot(plot_net, aes(x = from, y = to, fill = factor(Module, levels = mod_levels))) +
  geom_tile() +
  theme_bw() +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_manual(values = custom_pal, na.value = "grey70", drop = T) +
  labs(fill = "Module") +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())


# host tree with modules

host_tip_mod <- tibble(label = host_tree$tip.label) %>% 
  inner_join(filter(all_wmod_edited, age == 0), by = c("label" = "name"))

ggtree(host_tree, ladderize = F) %<+% host_tip_mod +
  geom_tippoint(aes(color = factor(module, levels = mod_levels)), size = 2, shape = "square") + 
  scale_color_manual(values = custom_pal,na.value = "grey70", drop = F)


*/