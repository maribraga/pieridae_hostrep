#'---
#'title: "Pieridae host repertoire - network evolution"
#'author: "Mariana Braga"
#'date: "`r format(Sys.time(), '%d %B, %Y')`"
#'output: github_document
#'---

#'-------------
#'
#' Script 3 for empirical study performed in Braga et al. 2020
#' *Evolution of butterfly-plant networks revealed by Bayesian inference of host repertoire*.
#'
#' (You must have finished script 2 first)

#+ include = FALSE

library(ape)
library(ggtree)
library(picante)

library(tidyverse)
library(patchwork)
library(wesanderson)

library(bipartite)
library(ggraph)
library(tidygraph)
library(igraph)

load("./inference/char_hist.RData")

/*# NETWORK STRUCTURE ----
*/

/*# Sizes over time ----
*/
  
#' ## Network size
#' 
#' First we will look at how network size changed over time in Pieridae.  
#' Let's look at the number of butterflies and plants, the number of interactions and modules.
#' We can also calculate some ratios. Here we have: *i.high* = mean number of interactions per butterfly,
#' *i.low* = mean number of interactions per host plant, *p.low* = number of hosts / number of butterflies,
#' and *connectance* = the proportion of possible interactions that are realized.

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

gg_size_n <- ggplot(net_size_n, aes(age,value, col = index)) +
  geom_line() +
  geom_point() +
  scale_x_reverse() +
  theme_bw()

gg_size_r <- ggplot(net_size_r, aes(age,value, col = index)) +
  geom_line() +
  geom_point() +
  scale_x_reverse() +
  theme_bw()

#+ net_size, fig.width = 8, fig.height = 3
gg_size_n + gg_size_r 


/*
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
*/