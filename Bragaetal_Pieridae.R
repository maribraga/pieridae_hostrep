#'---
#'title: "Pieridae host repertoire"
#'author: "Mariana Braga"
#'date: "4/22/2020"
#'output: github_document
#'---

#' Script for empirical study performed in Braga et al. 2020
#' Evolution of butterfly-plant networks revealed by Bayesian inference of host repertoire

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


#### INPUT DATA ####

# _Read trees ----
path_data <- "./data/"
tree <- read.tree(paste0(path_data,"bphy_pie_ladder.phy"))
host_tree <- read.tree(paste0(path_data,"angio_pie_50tips_ladder.phy"))

## for figure 2
host_tree_bl1 <- read.tree(paste0(path_data,"angio_pie_50tips_bl1.phy"))

ggtree(host_tree_bl1) + 
  theme_tree2() + 
  geom_tiplab(align=TRUE, linesize=.5) + 
  xlim(0, 25)


# _Extant network ----

ext_net_50h <- as.matrix(read.csv(paste0(path_data,"incidence_pieridae.csv"), header = T, row.names = 1))
identical(colnames(ext_net_50h), tree$tip.label)
identical(rownames(ext_net_50h), host_tree$tip.label)

ext_net <- ext_net_50h[which(rowSums(ext_net_50h) != 0),]

# JUMP to Character history


##### PARAMETER ESTIMATES #####

# _Read log files - output from MCMC ----
#path_logs <- "./Rev_analyses/server/output-NEW/"
path_logs <- "./Rev_analyses/for_paper/"

chain1 <- read.table(paste0(path_logs,"out.2.real.pieridae.2s.log"), header = TRUE)[,c(1,5,7:9)]
chain3 <- read.table(paste0(path_logs,"out.3.bl1.pieridae.2s.log"), header = TRUE)[,c(1,5,7:9)]
chain5 <- read.table(paste0(path_logs,"out.2.b0.pieridae.2s.log"), header = TRUE)[,c(1,5,8:9)]
colnames(chain1) <- colnames(chain3) <- c("generation","clock","beta", "lambda[02]", "lambda[20]")
colnames(chain5) <- c("generation","clock", "lambda[02]", "lambda[20]")

postr <- filter(chain1, generation >= 20000 & generation <= 200000)
postb <- filter(chain3, generation >= 20000 & generation <= 200000)
post0 <- filter(chain5, generation >= 20000 & generation <= 200000)

posterior <- bind_rows(postr,postb,post0) %>% 
  mutate(tree = c(rep("time", 3601),rep("bl1", 3601),rep("zero", 3601)))


# _Mean estimates ----
means <- group_by(posterior, tree) %>% 
  summarise_all(mean)
#  tree  generation  clock  beta  `lambda[02]` `lambda[20]`
#   bl1      110000 0.0186  1.48       0.0268        0.973
#  time      110000 0.0200  2.10       0.0347        0.965
#  zero      110000 0.0197    NA       0.0334        0.967


# _Density ----

nrow = 5001
priors <- tibble(
  generation = 1:nrow,
  clock = rexp(nrow, rate = 10),
  beta = rexp(nrow, rate = 1),
  lambdas = dbeta(seq(0,1,0.0002), 1, 3), # `lambda[02]`
  tree = rep("prior", nrow)
)

kd_beta <- kdensity(x = priors$beta, kernel='gamma', support=c(0,Inf), bw = 0.05)
kd_clock <- kdensity(x = priors$clock, kernel='gamma', support=c(0,Inf), bw = 0.01)
kd_lambda <- kdensity(x = priors$lambdas, kernel='gamma', support=c(0,Inf), bw = 0.01)

x10 = seq(0,10,0.002)
x1 = seq(0,1,0.0002)

# plot(x10, kd_beta(x10))
# plot(x1, kd_clock(x1))
# plot(x1, kd_lambda(x1))

plot_priors <- tibble(
  param = c(rep("beta", nrow), rep("clock", nrow), rep("lambda", nrow)),
  x = c(x10, x1, x1),
  y = c(kd_beta(x10), kd_clock(x1), kd_lambda(x1))
) 

all_post <- mutate(posterior, beta = case_when(is.na(beta) ~ 0, TRUE ~ beta),
                   tree = factor(tree, levels = c("zero","time","bl1")))


palwes <- wes_palette("Zissou1", 4, type = "continuous")
pal <- c(palwes[4],palwes[3],palwes[1])

ggbeta <- ggplot(all_post, aes(beta)) +
  geom_line(aes(x,y), data = filter(plot_priors, param == "beta")) +
  geom_density(aes(group = tree, fill = tree, col = tree), alpha = 0.7) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  xlim(c(0,7)) +
  theme_bw()

ggclock <- ggplot(all_post, aes(clock)) +
  geom_line(aes(x,y), data = filter(plot_priors, param == "clock")) +
  geom_density(aes(group = tree, fill = tree, col = tree), alpha = 0.7) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  xlim(c(0,0.1)) +
  theme_bw()

gg02 <- ggplot(all_post, aes(`lambda[02]`)) +
  geom_line(aes(x,y), data = filter(plot_priors, param == "lambda")) +
  geom_density(aes(group = tree, fill = tree, col = tree), alpha = 0.7) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  xlim(c(0,0.1)) +
  theme_bw()

gg20 <- ggplot(all_post, aes(`lambda[20]`)) +
  geom_line(aes(x,y), data = filter(plot_priors, param == "lambda")) +
  geom_density(aes(group = tree, fill = tree, col = tree), alpha = 0.7) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  xlim(c(0.9,1)) +
  theme_bw()

ggclock + ggbeta + gg02 + gg20 + plot_layout(ncol = 4, guides = 'collect')

# __Wrong prior density at 0 ----
# nrow = 5001
# priors <- tibble(
#   generation = 1:nrow,
#   clock = rexp(nrow, rate = 10),
#   beta = rexp(nrow, rate = 1),
#   `lambda[02]` = dbeta(seq(0,1,0.0002), 1, 3),
#   `lambda[20]` = dbeta(seq(0,1,0.0002), 1, 3),
#   tree = rep("prior", nrow)
# )
# 
# all_post <- bind_rows(mutate(posterior, beta = case_when(is.na(beta) ~ 0, TRUE ~ beta)), 
#                       priors)  %>%
#   mutate(tree = factor(tree, levels = c("prior","zero","time","bl1")))
# 
# palwes <- wes_palette("Zissou1", 4, type = "continuous")
# pal <- c("grey",palwes[4],palwes[3],palwes[1])
# 
# ggbeta <- ggplot(all_post, aes(beta)) +
#   geom_line(aes(), data = filter(plot_priors, param == "beta")) +
#   geom_density(aes(group = tree, fill = tree, col = tree), alpha = 0.7) +
#   scale_fill_manual(values = pal) +
#   scale_color_manual(values = pal) +
#   xlim(c(0,7)) +
#   theme_bw()
# 
# ggclock <- ggplot(all_post, aes(clock)) +
#   geom_density(aes(group = tree, fill = tree, col = tree), alpha = 0.7) +
#   scale_fill_manual(values = pal) +
#   scale_color_manual(values = pal) +
#   xlim(c(0,0.1)) +
#   theme_bw()
# 
# gg02 <- ggplot(all_post, aes(`lambda[02]`)) +
#   geom_density(aes(group = tree, fill = tree, col = tree), alpha = 0.7) +
#   scale_fill_manual(values = pal) +
#   scale_color_manual(values = pal) +
#   xlim(c(0,0.1)) +
#   theme_bw()
# 
# gg20 <- ggplot(all_post, aes(`lambda[20]`)) +
#   geom_density(aes(group = tree, fill = tree, col = tree), alpha = 0.7) +
#   scale_fill_manual(values = pal) +
#   scale_color_manual(values = pal) +
#   xlim(c(0.9,1)) +
#   theme_bw()
# 
# ggclock + ggbeta + gg02 + gg20 + plot_layout(ncol = 4, guides = 'collect')

# __Old way of doing density plots with kdensity ----
# prior_beta <- data.frame(sim = NA, beta = rexp(5001, rate = 1))
# prior_clock <- data.frame(sim = NA, clock = rexp(5001, rate = 10))
# 
# # get density function for all parameters
# for(m in c("time", "bl1", "zero")){
#   for(p in colnames(posterior)[c(2:5)]){
#     post <- filter(posterior, tree == m) %>% dplyr::select(p)
#     if(p %in% c("beta", "clock")) {
#       kdens <- kdensity(x = post[,1], kernel='gamma', support=c(0,Inf), bw = 0.01)
#     } else{
#       kdens <- kdensity(x = post[,1], kernel='gamma', support=c(0,Inf), bw = 0.001)
#     }
#     assign(paste0("kd_",p,"_",m,"_"),kdens)
# 
#   }
# }
# 
# 
# for(m in c("time", "bl1", "zero")){
#   for(p in colnames(posterior)[c(2:5)]){
#     post <- filter(posterior, tree == m) %>% dplyr::select(p)
#     if(!(p == "beta" & m == "zero")) {
#       if(p == "beta") {
#         kdens <- kdensity(x = post[,1], kernel='gamma', support=c(0,Inf), bw = 0.01)
#       } else {
#         kdens <- kdensity(x = post[,1], kernel='gamma', support=c(0,Inf), bw = 0.001)
#         }
#       assign(paste0("kd_",p,"_",m),kdens)
#     }
#   }
# }
# 
# 
# # create tibbles and plot
# z = seq(0,5,0.001)
# y = seq(0,1,0.0002)
# 
# # beta
# dens_beta <- tibble(x = z, prior = dexp(x = z, rate=1), time = kd_beta_time(x), bl1 = kd_beta_bl1(x))
# 
# ggplot(dens_beta) +
#   geom_line(aes(x, prior), col = "grey50", alpha = 0.8) +
#   geom_line(aes(x, time), col = "blue", alpha = 0.8) +
#   geom_line(aes(x, bl1), col = "orange", alpha = 0.8) +
#   labs(x = expression("Estimated phylogenetic-distance power, " ~ beta),
#        y = "Density") +
#   theme_light()
# 
# # clock
# dens_clock <- tibble(x = y, prior = dexp(x = z, rate=10), time = kd_clock_time(x), bl1 = kd_clock_bl1(x), 
#                      b0 = kd_clock_zero(x))
# ggplot(dens_clock) +
#   geom_line(aes(x, prior), col = "grey50", alpha = 0.8) +
#   geom_line(aes(x, time), col = "blue", alpha = 0.8) +
#   geom_density(aes(clock), data = postr) +
#   #geom_line(aes(x, bl1), col = "orange", alpha = 0.8) +
#   #geom_line(aes(x, b0), col = "red", alpha = 0.8) +
#   scale_x_continuous(limits = c(0,0.1)) +
#   labs(x = expression("Estimated rate of host-repertoire evolution, " ~ mu)) +
#   theme_light()
# 
# # rates
# prior_rates <- dbeta(y, 1, 3)
# dens_02 <- tibble(x = y, prior = prior_rates, time = `kd_lambda[02]_time`(x), bl1 = `kd_lambda[02]_bl1`(x), 
#                   b0 = `kd_lambda[02]_zero`(x))
# ggplot(dens_02) +
#   geom_line(aes(x, prior), col = "grey50", alpha = 0.8) +
#   geom_line(aes(x, time), col = "blue", alpha = 0.8) +
#   geom_line(aes(x, bl1), col = "orange", alpha = 0.8) +
#   geom_line(aes(x, b0), col = "red", alpha = 0.8) +
#   scale_x_continuous(limits = c(0,0.08)) +
#   theme_light()
# 
# dens_20 <- tibble(x = y, prior = prior_rates, time = `kd_lambda[20]_time`(x), bl1 = `kd_lambda[20]_bl1`(x), 
#                   b0 = `kd_lambda[20]_zero`(x))
# ggplot(dens_20) +
#   geom_line(aes(x, prior), col = "grey50", alpha = 0.8) +
#   geom_line(aes(x, time), col = "blue", alpha = 0.8) +
#   geom_density(aes(`lambda[20]`), data = postr) +
#   #geom_line(aes(x, bl1), col = "orange", alpha = 0.8) +
#   #geom_line(aes(x, b0), col = "red", alpha = 0.8) +
#   scale_x_continuous(limits = c(0.8,1)) +
#   theme_light()



# _Bayes factor ----

d_prior <- dexp(x=0, rate=1)

kd_beta_time <- kdensity(x = filter(posterior, tree == "time") %>% pull(beta), 
                         kernel='gamma', 
                         support=c(0,Inf), 
                         bw = 0.02)
kd_beta_bl1 <- kdensity(x = filter(posterior, tree == "bl1") %>% pull(beta), 
                         kernel='gamma', 
                         support=c(0,Inf), 
                         bw = 0.02)
max_time = kd_beta_time(0)
max_bl1 = kd_beta_bl1(0)

BF_time <- d_prior/max_time
# 7.029484
BF_bl1 <- d_prior/max_bl1
# 65545525



##### CHARACTER HISTORY #####

# Read in .history.txt file ----

path_out <- "./output/"

nhosts <- Ntip(host_tree)
hosts <- host_tree$tip.label

name1 <- "out.3.bl1.pieridae.2s.history.txt"
name2 <- "out.2.real.pieridae.2s.history.txt"
name3 <- "out.2.b0.pieridae.2s.history.txt"

name <- name1

colclasses <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))
history_dat = read.table(paste0(path_out,name), sep="\t", header=T, colClasses = colclasses)

history_dat <- filter(history_dat, iteration >= 20000 & iteration <= 200000) %>% 
  mutate(node_index = node_index + 1)


# Calculate effective rate of evolution ----
tree_length <- sum(tree$edge.length)+86
n_events <- group_by(history_dat,iteration) %>% 
  summarise(n = n()) %>% 
  summarise(mean = mean(n)) %>% 
  pull(mean)

rate <- n_events/tree_length
rate
# bl1 = 0.0950814 (out 3), out 2 = 0.0958
#---


## (States at nodes) ----

source("./code/functions_ancestral_states.R")

nodes <- c(84,128,129,130,131) #67:131

# 2 state model
pp <- make_matrix_nodes(history_dat, nodes, c(2))
row.names(pp) <- paste0("Index_",nodes)
colnames(pp) <- host_tree$tip.label

assign(paste0("pp_",name), pp)

graph <- graph_from_incidence_matrix(pp, weighted = TRUE)
el <- get.data.frame(graph, what = "edges") %>% 
  mutate(to = factor(to, levels = hosts),
         from = factor(from, levels = paste0("Index_",nodes))) %>% 
  rename(p = weight)

ppbl1 <- ggplot(el, aes(x = to, y = from)) +
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




## States at ages ----

# _Get posteriors at ages ----

source("./code/functions_ancestral_states.R")

ages = c(80,70,60,50,40,30,20,10,0)
#ages = c(70)

list_m_at_ages = list()
for (i in 1:(length(ages)-1)) {
  age = ages[i]
  list_m_at_ages[[i]] = t(make_matrix_at_age( history_dat, age, s_hit=c(2) ))
  #assign(paste0("m_",age), list_m_at_ages[[i]])
}

list_m_at_ages[[9]] <- ext_net

View(list_m_at_ages[[1]])
View(list_m_at_ages[[9]])

# _(Binary networks with a probability threshold) ----

# probability threshold
pt <- 80 
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
  matrix = matrix[ rowSums(matrix)!=0, ]
  matrix = matrix[ ,colSums(matrix)!=0 ]
  net_list[[m]] <- matrix
  #assign(paste0("net_",pt,"p_",ages[m]), matrix)
}

View(net_list[[1]])


# __(Output matrices for MODULAR) ----

# folder that will be read by MODULAR - change when repeating with different host tree
path_net <- "/Users/mari/Google\ Drive/P4\ -\ networks/host_tree_bl1/modular/"    
#path_net <- "/Users/mari/Google\ Drive/P4\ -\ networks/host_tree_time/modular/"    

# export networks as .txt files without row and column names
for(n in 1:length(net_list)){
  net <- net_list[[n]]
  name <- paste0(path_net,"net_",pt,"p_",ages[n],"Ma.txt")
  write.table(net, name, row.names = FALSE, col.names = FALSE)
}


# __(Make tidygraphs with MODULAR modules) ----

# called in the beginning as well
#path_mod <- paste0("./net_structure/modular/")
minp <- 10

for(n in 1:length(net_list)){
  # get names from nets
  net <- net_list[[n]]
  rnames <- rownames(net)
  cnames <- colnames(net)
  
  # get weighted graph from list_m_at_ages - after setting minimum weight
  matrix <- list_m_at_ages[[n]]
  for(i in 1:nrow(matrix)){
    for(j in 1:ncol(matrix)){
      if(matrix[i,j] < minp/100){
        matrix[i,j] = 0
      }
    }
  }
  matrix = matrix[ rowSums(matrix)!=0, ]
  matrix = matrix[ ,colSums(matrix)!=0 ]
  
  #graph <- graph_from_incidence_matrix(matrix, weighted = TRUE)
  
  nCs <- length(cnames)
  nRs <- length(rnames)
  levels <- c(paste0("C",1:nCs), paste0("R",1:nRs))
  
  # mod <- read.table(paste0(path_mod,"bl1/MEMBERS_net_",pt,"p_",ages[n],"Ma.txt"), header = TRUE) %>% 
  #   mutate(Node = factor(Node, levels = levels), module = Module) %>% 
  #   arrange(Node) %>% 
  #   mutate(name = c(cnames,rnames))
  
  
  mod <- read.csv(paste0(path_mod, "bl1/modules_pieridae.csv"), header = TRUE) %>% 
    filter(age == ages[n])

  tgraph <- as_tbl_graph(matrix) %>% 
    left_join(mod, by = "name") %>% 
    select(type, name, module) %>% 
    mutate(#type = case_when(type == FALSE ~ "butterfly",
           #                 type == TRUE ~ "plant"),
           degree = centrality_degree()) %>% 
    activate(edges) %>% 
    mutate(highp = case_when(weight >= 0.8 ~ TRUE,
                             weight < 0.8 ~ FALSE))
  
  # nodes <- get.data.frame(graph, what = "vertices") %>% 
  #   left_join(mod, by = "name")
  # V(graph)$module <- nodes$Module
  
  assign(paste0("tgraph_",ages[n]), tgraph)
}




# OR
# _Weighted networks with a low probability threshold ----

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

dim(list_wnets[[1]])
dim(list_wnets[[9]])

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


# NETWORK STRUCTURE ----

# Sizes over time ----

net_size <- tibble()

for(i in 1:length(ages)){
  ns <- list_tgraphs[[i]] %>% activate(what = "nodes") %>% as_tibble()
  es <- list_tgraphs[[i]] %>% activate(what = "edges") %>% as_tibble()
  n.high <- length(which(ns$type == TRUE))
  n.low <- length(which(ns$type == FALSE))
  n.mod <- length(unique(ns$module))
  n.int <- nrow(es)
  net_size <- bind_rows(net_size, tibble(age=ages[i], 
                                         n.high=n.high, 
                                         n.low=n.low, 
                                         n.mod=n.mod, 
                                         n.int=n.int,
                                         i.high=n.int/n.high,
                                         i.low=n.int/n.low,
                                         p.low=n.low/n.high,
                                         connectance=n.int/(n.high*n.low)))
}

net_size <- pivot_longer(net_size, 2:9, names_to = "index", values_to = "value")
net_size_n <- filter(net_size, index %in% c("n.high","n.low","n.mod","n.int"))
net_size_r <- filter(net_size, index %in% c("i.high","i.low","p.low","connectance"))

ggplot(net_size_n, aes(age,value, col = index)) +
  geom_line() +
  geom_point() +
  scale_x_reverse() +
  theme_bw()

ggplot(net_size_r, aes(age,value, col = index)) +
  geom_line() +
  geom_point() +
  scale_x_reverse() +
  theme_bw()


# Null models ----
# r00_samp (random), c0_samp (keeps column sums), r2d (keeps marginal sums)

Nulls_r00 <- list()
Nulls_c0 <- list()
Nulls_r2d <- list()

nit <- 1000

for(i in 1:length(ages)){
  null_r00 <- vegan::nullmodel(list_wnets[[i]], "r00_samp") 
  sim_r00 <- simulate(null_r00, nsim=nit, seed = 1)
  Nulls_r00[[i]] <- sim_r00
  
  null_c0 <- vegan::nullmodel(list_wnets[[i]], "c0_samp") 
  sim_c0 <- simulate(null_c0, nsim=nit, seed = 1)
  Nulls_c0[[i]] <- sim_c0
  
  if(ages[i] == 0){
    count <- list_wnets[[i]]
  } else {
    count <- round(list_wnets[[i]]*100)
  }
  
  null_r2d <- bipartite::nullmodel(count, N=nit, method=1)
  Nulls_r2d[[i]] <- null_r2d
  
}  

# second round of null models

Nulls_random <- list()
Nulls_colsums <- list()
Nulls_allsums <- list()

nit <- 1000

for(i in 1:length(ages)){
  
  if(ages[i] == 0){
    count <- list_wnets[[i]]
  } else {
    count <- round(list_wnets[[i]]*100)
  }
  
  null_r00 <- vegan::nullmodel(count, "r00_both") 
  sim_r00 <- simulate(null_r00, nsim=nit, seed = 1)
  Nulls_random[[i]] <- sim_r00
  
  null_c0 <- vegan::nullmodel(count, "c0_both") 
  sim_c0 <- simulate(null_c0, nsim=nit, seed = 1)
  Nulls_colsums[[i]] <- sim_c0
  
  null_swap <- vegan::nullmodel(count, "quasiswap_count") 
  sim_swap <- simulate(null_swap, nsim=nit, seed = 1)
  Nulls_allsums[[i]] <- sim_swap
  
} 

  
# _Modularity ----

# slow!

# Qnull <- tibble()
# 
# for(i in 1:length(ages)){
# 
#   sim_random <- Nulls_random[[i]]
#   sim_colsums <- Nulls_colsums[[i]]
#   sim_allsums <- Nulls_allsums[[i]]
# 
#   for(j in 1:nit){
#     Qrandom <- computeModules(sim_random[,,j])@likelihood
#     Qcolsums <- computeModules(sim_colsums[,,j])@likelihood
#     Qallsums <- computeModules(sim_allsums[,,j])@likelihood
# 
#     Qnull <- bind_rows(Qnull, tibble(age = ages[i], sim=j, allsums=Qallsums, random=Qrandom, colsums=Qcolsums))
#   }
# }

Qnull <- Qnull %>% pivot_longer(3:5, names_to = "model", values_to = "Q")


# observed 
Qobs <- tibble()

for(i in 1:length(ages)){
  q <- get(paste0("wmod_",ages[i]))
  Qobs<- bind_rows(Qobs, tibble(age = ages[i], Q = q@likelihood))
}


# _Nestedness ----

Nnull <- tibble()

for(i in 1:length(ages)){
  
  sim_random <- Nulls_random[[i]]
  sim_colsums <- Nulls_colsums[[i]]
  sim_allsums <- Nulls_allsums[[i]]  
  
  if(ages[i] == 0){
    for(j in 1:nit){
      Nrandom <- networklevel(sim_random[,,j],index="NODF")
      Ncolsums <- networklevel(sim_colsums[,,j],index="NODF")
      Nallsums <- networklevel(sim_allsums[,,j],index="NODF")
      Nnull <- bind_rows(Nnull, tibble(age = ages[i], sim=j, 
                                       colsums=Ncolsums, random=Nrandom, allsums=Nallsums))
    }
  } else {
    for(j in 1:nit){
      Nrandom <- networklevel(sim_random[,,j],index="weighted NODF")
      Ncolsums <- networklevel(sim_colsums[,,j],index="weighted NODF")
      Nallsums <- networklevel(sim_allsums[,,j],index="weighted NODF")
      Nnull <- bind_rows(Nnull, tibble(age = ages[i], sim=j,
                                       colsums=Ncolsums, random=Nrandom, allsums=Nallsums))
    }
  }
}

Nnull <- Nnull %>% pivot_longer(3:5, names_to = "model", values_to = "N")


# observed
Nobs <- tibble()

for(i in 1:length(ages)){
  if(ages[i] == 0){
    nodf <- networklevel(list_wnets[[i]],index="NODF")
    Nobs<- bind_rows(Nobs, tibble(age = ages[i], nodf = nodf))
  } else {
    nodf <- networklevel(list_wnets[[i]],index="weighted NODF")
    Nobs<- bind_rows(Nobs, tibble(age = ages[i], nodf = nodf))
  }
}


# _Density plots ----

# modularity
for(g in ages){
  
  gg <- ggplot(filter(Qnull, age == g)) +
    geom_density(aes(Q, group = model, color = model)) +
    scale_color_grey() +
    geom_vline(xintercept = filter(Qobs, age == g)$Q, col = "blue") +
    labs(title = paste0(g," Ma"), y = NULL) +
    theme_bw() + 
    theme(axis.ticks.y = element_line(linetype = "blank"), 
          axis.text.y = element_blank())
  
  assign(paste0("mod_plot_",g), gg)
}

(mod_plot_80 | mod_plot_70 | mod_plot_60) /
  (mod_plot_50 | mod_plot_40 | mod_plot_30) /
  (mod_plot_20 | mod_plot_10 | mod_plot_0) +
  plot_layout(guides = 'collect')


# nestedness
for(g in ages){
  gg <- ggplot(filter(Nnull, age == g)) +
    geom_density(aes(N, group = model, color = model)) +
    scale_color_grey() +
    geom_vline(xintercept = filter(Nobs, age == g)$nodf, col = "blue") +
    labs(title = paste0(g," Ma"), y = NULL) +
    theme_bw() + 
    theme(axis.ticks.y = element_line(linetype = "blank"), 
          axis.text.y = element_blank())
  
  assign(paste0("nodf_plot_",g), gg)
}

(nodf_plot_80 | nodf_plot_70 | nodf_plot_60) /
  (nodf_plot_50 | nodf_plot_40 | nodf_plot_30) /
  (nodf_plot_20 | nodf_plot_10 | nodf_plot_0) +
  plot_layout(guides = 'collect')


# _Z-scores ----

Qzscore <- Qnull %>% 
  group_by(age, model) %>% 
  summarize(mean = mean(Q),
            sd = sd(Q)) %>% 
  left_join(Qobs) %>% 
  mutate(z = (Q - mean)/sd) %>% 
  left_join(Qnull %>% 
              left_join(rename(Qobs, Qobs = Q)) %>% 
              group_by(age, model) %>% 
              summarise(p = sum(Q > Qobs)/nit))


Nzscore <- Nnull %>% 
  group_by(age, model) %>% 
  summarize(mean = mean(N),
            sd = sd(N)) %>% 
  left_join(Nobs) %>% 
  mutate(z = (nodf - mean)/sd) %>% 
  left_join(Nnull %>% 
              left_join(rename(Nobs, Nobs = nodf)) %>% 
              group_by(age, model) %>% 
              summarise(p = sum(N > Nobs)/nit))


# colors from parameter estimates plot
palwes <- wes_palette("Zissou1", 4, type = "continuous")
pal <- c(palwes[4],palwes[3],palwes[1])

plot_qz <- ggplot(Qzscore) +
  geom_line(aes(age, z, group = model, col = model)) +
  geom_point(aes(age, z, group = model, col = model),
             data = filter(Qzscore, p <= 0.05),
             size = 2, alpha = 0.7) +
  scale_color_manual(values = pal) +
  scale_x_reverse() +
  labs(title = "Modularity, Q", y = "Z-score", x = "Millions of years ago, Ma", col = "Null model") +
  theme_bw()

plot_nz <- ggplot(Nzscore) +
  geom_line(aes(age, z, group = model, col = model)) +
  geom_point(aes(age, z, group = model, col = model),
             data = filter(Nzscore, p <= 0.05),
             size = 2, alpha = 0.7) +
  scale_color_manual(values = pal) +
  scale_x_reverse() +
  labs(title = "Nestedness, N", y = "Z-score", x = "Millions of years ago, Ma", col = "Null model") +
  theme_bw()

plot_qz / plot_nz + plot_layout(guides = 'collect')



# Checking null networks

net0 <- Nulls_r00[[7]][,,1]
net1 <- Nulls_c0[[7]][,,1]
net2 <- Nulls_r2d[[7]][[1]]
net <- list_wnets[[7]]

net0 <- Nulls_random[[8]][,,10]
net1 <- Nulls_colsums[[8]][,,10]
net2 <- Nulls_allsums[[8]][,,10]
net <- list_wnets[[8]]

par(mfrow = c(2,2))

visweb(net, type = "none", labsize = 0.1) #prednames = F, preynames = F)
visweb(net2, type = "none", labsize = 0.1)
visweb(net1, type = "none", labsize = 0.1)
visweb(net0, type = "none", labsize = 0.1)

visweb(net, type = "diagonal", labsize = 0.1) #prednames = F, preynames = F)
visweb(net2, type = "diagonal", labsize = 0.1)
visweb(net1, type = "diagonal", labsize = 0.1)
visweb(net0, type = "diagonal", labsize = 0.1)

visweb(net, type = "nested", labsize = 0.1) 
visweb(net2, type = "nested", labsize = 0.1)
visweb(net1, type = "nested", labsize = 0.1)
visweb(net0, type = "nested", labsize = 0.1)


plotModuleWeb(wmod_80)



# Phylogenetic diversity ----

# _Butterflies ----

butterflies <- setdiff(all_wmod_edited$name, host_tree$tip.label)
mod_matrix <- filter(all_wmod_edited, name %in% butterflies) %>% 
  frame2webs(c("module","name","age"))

mod_matrix <- mod_matrix[rev(names(mod_matrix))]

list_but_pd <- list()
for(i in 1:length(ages)){
  bpd <- ses.pd(mod_matrix[[i]], list_subtrees[[i]], null.model="taxa.labels") 
  list_but_pd[[i]] <- bpd
}

str(list_but_pd)


# _Plants ----

mod_matrix_hosts <- filter(all_wmod_edited, name %in% host_tree$tip.label) %>% 
  frame2webs(c("module","name","age"))

mod_matrix_hosts <- mod_matrix_hosts[rev(names(mod_matrix_hosts))]

list_host_pd <- list()
for(i in 1:length(ages)){
  
  out <- setdiff(host_tree$tip.label,colnames(mod_matrix_hosts[[i]]))
  host_subtree <- drop.tip(host_tree, out)
  ppd <- ses.pd(mod_matrix_hosts[[i]], host_subtree, null.model="taxa.labels") 
  list_host_pd[[i]] <- ppd
}

str(list_host_pd)


mod_pd <- tibble()
for(i in 1:length(ages)){
  age <- ages[i]
  ppd <- list_host_pd[[i]]
  bpd <- list_but_pd[[i]]
  
  btbl <- tibble(age=age, type="butterfly", module=rownames(bpd), 
                 PDz=bpd$pd.obs.z, PDp=bpd$pd.obs.p, ntaxa=bpd$ntaxa)
  ptbl <- tibble(age=age, type="plant", module=rownames(ppd), 
                 PDz=ppd$pd.obs.z, PDp=ppd$pd.obs.p, ntaxa=ppd$ntaxa)
  
  mod_pd <- bind_rows(mod_pd, btbl, ptbl)
}


gpdb <- ggplot(filter(mod_pd, type == "butterfly")) +
  geom_line(aes(age,PDz, group=factor(module, levels = mod_levels), 
                col=factor(module, levels = mod_levels))) +
  geom_point(aes(age,PDz, group=factor(module, levels = mod_levels), 
                 col=factor(module, levels = mod_levels)),
             data = filter(mod_pd, PDp <= 0.05, type == "butterfly"),
             size = 2) +
  scale_color_manual(values = custom_pal,na.value = "grey70", drop = F) +
  scale_x_reverse() +
  labs(title = "Butterflies", color = "Module", x = "Millions of years ago, Ma", y = "PD z-score") +
  theme_bw()

gpdp <- ggplot(filter(mod_pd, type == "plant")) +
  geom_line(aes(age,PDz, group=factor(module, levels = mod_levels), 
                col=factor(module, levels = mod_levels))) +
  geom_point(aes(age,PDz, group=factor(module, levels = mod_levels), 
                 col=factor(module, levels = mod_levels)),
             data = filter(mod_pd, PDp <= 0.05, type == "plant"),
             size = 2) +
  scale_color_manual(values = custom_pal,na.value = "grey70", drop = F) +
  scale_x_reverse() +
  labs(title = "Plants", color = "Module", x = "Millions of years ago, Ma", y = "PD z-score") +
  theme_bw()


gpdb / gpdp + plot_layout(guides = 'collect')


# Number of nodes and modules through time ----

gnb <- ggplot() + 
  geom_point(data=filter(mod_pd, type == "butterfly"), 
             aes(age,ntaxa, col = factor(module, levels = mod_levels)),
             alpha = 0.7) + 
  geom_line(data=filter(mod_pd, type == "butterfly"),
            aes(age,ntaxa, col = factor(module, levels = mod_levels))) + 
  #geom_point(data=filter(mod_pd, type == "plant"), aes(age,ntaxa, col = factor(module, levels = mod_levels)), shape = "square") + 
  scale_color_manual(values = custom_pal,na.value = "grey70", drop = F) +
  #scale_shape_manual(values = c("square","circle")) +
  scale_x_reverse() +
  labs(color = "Module", x = "Millions of years ago, Ma", y = "Number of butterfly taxa") +
  theme_bw()
            
gnp <- ggplot() + 
  geom_point(data=filter(mod_pd, type == "plant"), 
              aes(age,ntaxa, col = factor(module, levels = mod_levels)),
              alpha = 0.7) + 
  geom_line(data=filter(mod_pd, type == "plant"),
            aes(age,ntaxa, col = factor(module, levels = mod_levels))) + 
  scale_color_manual(values = custom_pal,na.value = "grey70", drop = F) +
  scale_x_reverse() +
  labs(color = "Module", x = "Millions of years ago, Ma", y = "Number of plant taxa") +
  theme_bw()

gnb / gnp + plot_layout(guides = 'collect')





# -------- END of main script -------- #


#set.seed(2)

cz <- czvalues(wmod_0, weighted = T, level = "lower")
plot(cz[[1]], cz[[2]], pch=16, xlab="c", ylab="z", cex=0.8, xlim=c(0,1), las=1)
abline(v=0.62) # threshold of Olesen et al. 2007
abline(h=2.5)   # dito
text(cz[[1]], cz[[2]], names(cz[[1]]), pos=4, cex=0.7)

# example for computing a c- or z-threshold:
czobs <- czvalues(wmod_20, weighted = T, level = "lower")
nulls <- nullmodel(list_wnets[[7]], N=10) # this should be larger, of course
null.mod.list <- sapply(nulls, computeModules)
null.cz <- lapply(null.mod.list, czvalues)
# compute 95
null.cs <- sapply(null.cz, function(x) x$c) # c-values across all species in nulls
quantile(null.cs, 0.95) 
# this could now serve as thresholds for identifying particularly uncommonly high c-values
# and analogously for z, of course




# ----------- Plot all nets with ggplot ---------- #

tree <- list_subtrees[[2]]
graph <- list_tgraphs[[2]] %>% 
     activate(what = "nodes") %>%
     mutate(degree = centrality_degree(weights = weight))

mod <- all_wmod_edited %>% filter(age == ages[1])

edge_list <- get.data.frame(graph, what = "edges") %>%
  inner_join(mod %>% select(name, module), by = c("from" = "name")) %>%
  inner_join(mod %>% select(name, module), by = c("to" = "name")) %>%
  mutate(Module = ifelse(module.x == module.y, module.x, NA)) #%>% factor())

# by phylo
phylob <- tree$tip.label
phylop <- host_tree$tip.label

plot_net <- edge_list %>% mutate(
  to = factor(to, levels = phylob),
  from = factor(from, levels = phylop))

# by degree
korder1 <- graph %>% as_tibble() %>% filter(type == TRUE) %>% arrange(degree) %>% pull(name)
korder2 <- graph %>% as_tibble() %>% filter(type == FALSE) %>% arrange(desc(degree)) %>% pull(name)

plot_net <- edge_list %>% mutate(
  to = factor(to, levels = korder1),
  from = factor(from, levels = korder2))


ggplot(plot_net, aes(x = from, y = to, fill = factor(Module, levels = mod_levels), alpha = weight)) +
  geom_tile() +
  theme_bw() +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_manual(values = custom_pal, na.value = "grey70", drop = F) +
  labs(fill = "Module") +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none")

# ---------------------------------------- #






ggn_80 + geom_node_text(aes(label = name), repel = T)


#Plot tree for each age

plot(tree_10, cex = 0.5)
axisPhylo()


# tests

ages
tree <- read.tree(paste0(path_data,"tree_index.tre"), node.label = paste0("Index",67:131))
tree$Nnode

# Trying to set node labels to match RevBayes output
tree <- read.newick(paste0(path_data,"tree_nodelab.tre"), node.label = "label")
ggtree(tree) + geom_nodelab() + geom_tiplab() + theme_tree2()
#same as
#plot(tree, show.node.label = T)
#
# node labels are read if using treeio (read.newick instead of read.tree)
# but treeSlice ignores it!
root_age <- node.depth.edgelength(tree)[1]
for(i in ages){
  subtree <- treeSlice(tree, root_age - i, orientation = "rootwards")
  assign(paste0("tree_",i), subtree)
}

plot(tree_80)

# Trying to set node labels to match RevBayes output, but nodes are labeled from the root in R and from tips in rb
tree2 <- makeNodeLabel(tree)
tree2$node.label
tree2$node.label <- paste0("Index",131:67)
plot(tree2, show.node.label = T)
# ---







# # get weighted edge lists from list_m_at_ages (to plot probabilities) 
# el_ages <- tibble()
# for (i in 1:length(list_m_at_ages)){
#   graph <- graph_from_incidence_matrix(list_m_at_ages[[i]], weighted = TRUE)
#   nodes <- get.data.frame(graph, what = "vertices")
#   node_list <- full_join(nodes, nodes_mod) %>% 
#     dplyr::select(name, mod2) %>% 
#     dplyr::rename(mod = mod2)
#   
#   el <- get.data.frame(graph, what = "edges") %>% 
#     inner_join(node_list %>% dplyr::select(name, mod), by = c("to" = "name")) %>%
#     mutate(age = ages[i],
#            to = factor(to, levels = hosts),
#            from = factor(from, levels = c(rev(tree$tip.label), paste0("Index_",67:131)))) %>% 
#     rename(p = weight)
#   
#   el_ages <- bind_rows(el_ages,el)
# }
