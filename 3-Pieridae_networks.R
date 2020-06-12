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
#' This is a continuation of script 2 - Character history, so make sure you complete that one first.

#+ include = FALSE

library(picante)
library(tidyverse)
library(patchwork)
library(viridis)
library(bipartite)
library(ggraph)
library(tidygraph)
#library(igraph)

load("./inference/char_hist.RData")


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
net_size90 <- tibble()

for(i in 1:length(ages)){
  ns <- list_tgraphs[[i]] %>% activate(what = "nodes") %>% as_tibble()
  es <- list_tgraphs[[i]] %>% activate(what = "edges") %>% as_tibble()
  ns90 <- list_tgraphs[[i]] %>% activate(what = "nodes") %>% as_tibble() %>% filter(!is.na(module))
  es90 <- list_tgraphs[[i]] %>% activate(what = "edges") %>% as_tibble() %>% filter(weight >= 0.9)
  
  n.high <- length(which(ns$type == TRUE))
  n.low <- length(which(ns$type == FALSE))
  #n.mod <- length(unique(ns$module))
  n.int <- nrow(es)
  n90.high <- length(which(ns90$type == TRUE))
  n90.low <- length(which(ns90$type == FALSE))
  #n90.mod <- 
  n90.int <- nrow(es90)  
  
  net_size <- bind_rows(net_size, tibble(age=ages[i], 
                                         n.high=n.high, 
                                         n.low=n.low,
                                         n.int=n.int,
                                         i.high=n.int/n.high))
  
  net_size90 <- bind_rows(net_size90, tibble(age=ages[i],
                                           n.high=n90.high, 
                                           n.low=n90.low, 
                                           n.int=n90.int,
                                           i.high=n90.int/n90.high))
}

net_size <- pivot_longer(net_size, 2:5, names_to = "index", values_to = "value") %>% 
  mutate(prob = "all")
net_size90 <- pivot_longer(net_size90, 2:5, names_to = "index", values_to = "value") %>% 
  mutate(prob = "high")
net_size <- bind_rows(net_size,net_size90)

# 2-color palette
pal_2c <- viridis_pal()(3)[2:1]

ggplot(net_size, aes(age,value,col = prob)) +
  geom_line(aes(group=prob)) +
  geom_point() +
  scale_color_manual(values = pal_2c) +
  facet_grid(factor(index, levels = c("n.int","n.high","n.low","i.high"))~., scales = "free") +
  scale_x_reverse() +
  theme_bw()


/* # Effective partners over time ----
*/
  
#' ### Effective partners
#' 

net_ep <- tibble()

for(i in 1:length(ages)){
  
  net <- net_list[[i]]
  wnet <- list_wnets[[i]]
  
  m_ep <- mean(specieslevel(net, index = "effective partners", level = "higher")$effective.partners)
  mw_ep <- mean(specieslevel(wnet, index = "effective partners", level = "higher")$effective.partners)
  
  net_ep <- bind_rows(net_ep, tibble(age=ages[i], network = c('high-binary','all-weighted'), EP = c(m_ep, mw_ep)))
  
}

#+ fig.width = 5.5, fig.height = 2
ggplot(net_ep, aes(age,EP,col = network)) +
  geom_line(aes(group=network)) +
  geom_point() +
  scale_color_manual(values = pal_2c) +
  scale_x_reverse() +
  theme_bw()



/*# BINARY networks ---- 
*/
  
#' ## Network structure - Binary networks
#' 
/* # fix this later
#' We are gonna generate networks using two null models.
#' `r00` only keeps the size and number of interactions of the original network,
#' while `r2dtable` also keeps the marginal sums, i.e. the number of interactions 
#' that each butterfly/plant has.
*/
/*# Null models ----
*/

#' ### Null models
#' 

# number of null networks to be generated
nit <- 1000

#+ eval = FALSE  
Nulls_r00 <- list()
#Nulls_r2d <- list()

# gotta skip network at 80 Ma because it only has one host
for(i in 2:length(ages)){
  null_r00 <- vegan::nullmodel(net_list[[i]], "r00") 
  sim_r00 <- simulate(null_r00, nsim=nit, seed = 1)
  Nulls_r00[[i]] <- sim_r00
  
  #null_r2d <- vegan::nullmodel(net_list[[i]], "r2dtable") 
  #sim_r2d <- simulate(null_r2d, nsim=nit, seed = 1)
  #Nulls_r2d[[i]] <- sim_r2d
  
}  

  
/*# _Modularity ----
*/

#' ### Modularity
#' 
#' Let's calculate the modularity of each null network,
#' except for 80 and 70 Ma. They are too small.
  
source("my_compute_module.R")  

#+ eval = FALSE  
Qnull <- tibble()

for(i in 3:length(ages)){

  sim_random <- Nulls_r00[[i]]
  #sim_allsums <- Nulls_r2d[[i]]

  for(j in 1:nit){
    Qrandom <- mycomputeModules(sim_random[,,j])@likelihood
    #Qallsums <- mycomputeModules(sim_allsums[,,j])@likelihood

    Qnull <- bind_rows(Qnull, tibble(age = ages[i], sim=j, random=Qrandom)) #allsums=Qallsums,
  }
}

Qnull <- Qnull %>% pivot_longer(3, names_to = "model", values_to = "Q")

/*
saveRDS(Qnull, "./networks/Qnull.rds")
*/

#+
# Modularity of binary null networks  
Qnull <- readRDS("./networks/Qnull.rds")

    
# observed modularity for binary networks
Qobs <- tibble()

for(i in 3:length(ages)){
  q <- get(paste0("mod_",ages[i]))
  Qobs<- bind_rows(Qobs, tibble(age = ages[i], Q = q@likelihood))
}


# _Nestedness ----

#+ eval = FALSE
# do not calculate for networks at 80 and 70 Ma. They are too small. 
Nnull <- tibble()

for(i in 3:length(ages)){
  
  sim_random <- Nulls_r00[[i]]
  #sim_allsums <- Nulls_r2d[[i]]  
  
  for(j in 1:nit){
    Nrandom <- networklevel(sim_random[,,j],index="NODF")
    #Nallsums <- networklevel(sim_allsums[,,j],index="NODF")
    Nnull <- bind_rows(Nnull, tibble(age = ages[i], sim=j, random=Nrandom)) #, allsums=Nallsums))
  }
}

Nnull <- Nnull %>% pivot_longer(3, names_to = "model", values_to = "N")

/*
saveRDS(Nnull, "./networks/Nnull.rds")
*/
  
#+ 
# Nestedness of binary null networks
Nnull <- readRDS("./networks/Nnull.rds")


# observed nestedness for binary networks
Nobs <- tibble()

for(i in 3:length(ages)){
    nodf <- networklevel(net_list[[i]],index="NODF")
    Nobs<- bind_rows(Nobs, tibble(age = ages[i], nodf = nodf))
}


/*# _Z-scores ----
*/
#' ### Z-scores
  
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


/*# WEIGHTED networks ---- 
*/
#' ## Network structure - Weighted networks
#' 

/*# Null models ----
*/
#' ### Null models
#'

Nulls_rb <- list()
#Nulls_colsums <- list()
#Nulls_allsums <- list()

nit <- 1000

for(i in 1:(length(ages)-1)){
  
  count <- round(list_wnets[[i]]*100)
  
  null_rb <- vegan::nullmodel(count, "r00_both") 
  sim_rb <- simulate(null_rb, nsim=nit, seed = 1)
  Nulls_rb[[i]] <- sim_rb
  
  # null_c0 <- vegan::nullmodel(count, "c0_both") 
  # sim_c0 <- simulate(null_c0, nsim=nit, seed = 1)
  # Nulls_colsums[[i]] <- sim_c0
  
  # null_swap <- vegan::nullmodel(count, "quasiswap_count") 
  # sim_swap <- simulate(null_swap, nsim=nit, seed = 1)
  # Nulls_allsums[[i]] <- sim_swap
  
} 

  
/*# _Modularity ----
*/
#' ### Modularity
#' 

#+ eval = FALSE
Qwnull <- tibble()

for(i in 1:(length(ages)-1)){
  
  simw_random <- Nulls_rb[[i]]

  for(j in 1:nit){
    Qwrandom <- computeModules(simw_random[,,j])@likelihood
    Qwnull <- bind_rows(Qwnull, tibble(age = ages[i], sim=j, random=Qwrandom))
  }
}

Qwnull <- Qwnull %>% pivot_longer(3, names_to = "model", values_to = "Q")

/*
saveRDS(Qwnull, "./networks/Qwnull.rds")
*/
  
#+
# Modularity of weighted null networks 
Qwnull <- readRDS("./networks/Qwnull.rds")

# observed
Qwobs <- tibble()

for(i in 1:(length(ages)-1)){
  q <- get(paste0("wmod_",ages[i]))
  Qwobs<- bind_rows(Qwobs, tibble(age = ages[i], Q = q@likelihood))
}


/*# _Nestedness ----
*/
#' ### Nestedness
#'   

#+ eval = FALSE
Nwnull <- tibble()

for(i in 1:(length(ages)-1)){
  
  sim_rb <- Nulls_rb[[i]]

  for(j in 1:nit){
    Nrb <- networklevel(sim_rb[,,j],index="weighted NODF")

    Nwnull <- bind_rows(Nwnull, tibble(age = ages[i], sim=j,random=Nrb))
  }
}

Nwnull <- Nwnull %>% pivot_longer(3, names_to = "model", values_to = "N")

/*
saveRDS(Nwnull, "./networks/Nwnull.rds")
*/
  
#+ 
# Nestedness of weighted null networks
Nwnull <- readRDS("./networks/Nwnull.rds")

# observed
Nwobs <- tibble()

for(i in 1:(length(ages)-1)){
  wnodf <- networklevel(list_wnets[[i]],index="weighted NODF")
  Nwobs<- bind_rows(Nwobs, tibble(age = ages[i], nodf = wnodf))
}

    
/*# _Z-scores ----
*/
#' ### Z-scores
#' 
    
Qwzscore <- Qwnull %>%
  group_by(age, model) %>%
  summarize(mean = mean(Q),
            sd = sd(Q)) %>%
  left_join(Qwobs) %>%
  mutate(z = (Q - mean)/sd) %>%
  left_join(Qwnull %>%
              left_join(rename(Qwobs, Qobs = Q)) %>%
              group_by(age, model) %>%
              summarise(p = sum(Q > Qobs)/nit)) %>% 
  bind_rows(filter(Qzscore, age == 0))

Nwzscore <- Nwnull %>% 
  group_by(age, model) %>% 
  summarize(mean = mean(N),
            sd = sd(N)) %>% 
  left_join(Nwobs) %>% 
  mutate(z = (nodf - mean)/sd) %>% 
  left_join(Nwnull %>% 
              left_join(rename(Nwobs, Nobs = nodf)) %>% 
              group_by(age, model) %>% 
              summarise(p = sum(N > Nobs)/nit)) %>% 
  bind_rows(filter(Nzscore, age == 0))

Qz <- bind_rows(filter(Qzscore, model == 'random') %>% mutate(network = 'binary'), 
                filter(Qwzscore, model == 'random')%>% mutate(network = 'weighted'))

Nz <- bind_rows(filter(Nzscore, model == 'random') %>% mutate(network = 'binary'), 
                filter(Nwzscore, model == 'random')%>% mutate(network = 'weighted'))

plot_qz <- ggplot(Qz) +
  geom_line(aes(age, z, group = network, col = network)) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Qz, p <= 0.05),
             size = 2, alpha = 0.7) +
  scale_color_manual(values = pal_2c) +
  scale_x_reverse() +
  labs(title = "Modularity, Q", y = "Z-score", x = "Millions of years ago, Ma") +
  theme_bw()

plot_nz <- ggplot(Nz) +
  geom_line(aes(age, z, group = network, col = network)) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Nz, p <= 0.05),
             size = 2, alpha = 0.7) +
  scale_color_manual(values = pal_2c) +
  scale_x_reverse() +
  labs(title = "Nestedness, N", y = "Z-score", x = "Millions of years ago, Ma") +
  theme_bw()

#+ zscore, fig.width = 6, fig.height = 4.5, warning = FALSE
plot_qz / plot_nz + plot_layout(guides = 'collect')


/*
# plot_qz <- ggplot(Qzscore) +
#   geom_line(aes(age, z, group = model, col = model)) +
#   geom_point(aes(age, z, group = model, col = model),
#              data = filter(Qzscore, p <= 0.05),
#              size = 2, alpha = 0.7) +
#   scale_color_manual(values = pal) +
#   scale_x_reverse() +
#   labs(title = "Modularity, Q", y = "Z-score", x = "Millions of years ago, Ma", col = "Null model") +
#   theme_bw()
# 
# plot_nz <- ggplot(Nzscore) +
#   geom_line(aes(age, z, group = model, col = model)) +
#   geom_point(aes(age, z, group = model, col = model),
#              data = filter(Nzscore, p <= 0.05),
#              size = 2, alpha = 0.7) +
#   scale_color_manual(values = pal) +
#   scale_x_reverse() +
#   labs(title = "Nestedness, N", y = "Z-score", x = "Millions of years ago, Ma", col = "Null model") +
#   theme_bw()
# 
# #+ zscore, fig.width = 6.6, fig.height = 5.3, warning = FALSE
# plot_qz / plot_nz + plot_layout(guides = 'collect')
*/

# Phylogenetic diversity ----

# _Butterflies ----

# __Binary networks ----
butterflies <- setdiff(all_mod_edited$name, host_tree$tip.label)
mod_matrix <- filter(all_mod_edited, name %in% butterflies) %>% 
  frame2webs(c("module","name","age"))

mod_matrix <- mod_matrix[rev(names(mod_matrix))]

list_but_pd <- list()
for(i in 2:length(ages)){
  bpd <- ses.pd(mod_matrix[[i]], list_subtrees[[i]], null.model="taxa.labels") 
  list_but_pd[[i]] <- bpd
}

# __Weighted networks ----
butterflies <- setdiff(all_wmod_edited$name, host_tree$tip.label)
wmod_matrix <- filter(all_wmod_edited, name %in% butterflies) %>% 
  frame2webs(c("module","name","age"))

wmod_matrix <- wmod_matrix[rev(names(wmod_matrix))]

list_but_pdw <- list()
for(i in 2:length(ages)){
  bpdw <- ses.pd(wmod_matrix[[i]], list_subtreesw[[i]], null.model="taxa.labels") 
  list_but_pdw[[i]] <- bpdw
}

# _Plants ----

# __Binary networks ----
mod_matrix_hosts <- filter(all_mod_edited, name %in% host_tree$tip.label) %>% 
  frame2webs(c("module","name","age"))

mod_matrix_hosts <- mod_matrix_hosts[rev(names(mod_matrix_hosts))]

list_host_pd <- list()
for(i in 2:length(ages)){
  
  out <- setdiff(host_tree$tip.label,colnames(mod_matrix_hosts[[i]]))
  host_subtree <- drop.tip(host_tree, out)
  ppd <- ses.pd(mod_matrix_hosts[[i]], host_subtree, null.model="taxa.labels") 
  list_host_pd[[i]] <- ppd
}

# __Weighted networks ----
wmod_matrix_hosts <- filter(all_wmod_edited, name %in% host_tree$tip.label) %>% 
  frame2webs(c("module","name","age"))

wmod_matrix_hosts <- wmod_matrix_hosts[rev(names(wmod_matrix_hosts))]

list_host_pdw <- list()
for(i in 2:length(ages)){
  
  out <- setdiff(host_tree$tip.label,colnames(wmod_matrix_hosts[[i]]))
  host_subtree <- drop.tip(host_tree, out)
  ppd <- ses.pd(wmod_matrix_hosts[[i]], host_subtree, null.model="taxa.labels") 
  list_host_pdw[[i]] <- ppd
}

#-- 

mod_pd <- tibble()
for(i in 2:length(ages)){
  age <- ages[i]
  ppd <- list_host_pd[[i]]
  bpd <- list_but_pd[[i]]
  
  btbl <- tibble(age=age, type="butterfly", module=rownames(bpd), 
                 PDz=bpd$pd.obs.z, PDp=bpd$pd.obs.p, ntaxa=bpd$ntaxa)
  ptbl <- tibble(age=age, type="plant", module=rownames(ppd), 
                 PDz=ppd$pd.obs.z, PDp=ppd$pd.obs.p, ntaxa=ppd$ntaxa)
  
  mod_pd <- bind_rows(mod_pd, btbl, ptbl)
}

wmod_pd <- tibble()
for(i in 2:length(ages)){
  age <- ages[i]
  ppd <- list_host_pdw[[i]]
  bpd <- list_but_pdw[[i]]
  
  btbl <- tibble(age=age, type="butterfly", module=rownames(bpd), 
                 PDz=bpd$pd.obs.z, PDp=bpd$pd.obs.p, ntaxa=bpd$ntaxa)
  ptbl <- tibble(age=age, type="plant", module=rownames(ppd), 
                 PDz=ppd$pd.obs.z, PDp=ppd$pd.obs.p, ntaxa=ppd$ntaxa)
  
  wmod_pd <- bind_rows(wmod_pd, btbl, ptbl)
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
  labs(title = "Butterflies - binary networks", color = "Module", x = "Millions of years ago, Ma", y = "PD z-score") +
  theme_bw() 

gpdbw <- ggplot(filter(wmod_pd, type == "butterfly")) +
  geom_line(aes(age,PDz, group=factor(module, levels = wmod_levels), 
                col=factor(module, levels = mod_levels))) +
  geom_point(aes(age,PDz, group=factor(module, levels = wmod_levels), 
                 col=factor(module, levels = mod_levels)),
             data = filter(wmod_pd, PDp <= 0.05, type == "butterfly"),
             size = 2) +
  scale_color_manual(values = custom_pal,na.value = "grey70", drop = F) +
  scale_x_reverse() +
  labs(title = "Butterflies - weighted networks", color = "Module", x = "Millions of years ago, Ma", y = "PD z-score") +
  theme_bw() +
  theme(legend.position = "none")

# gpdp <- ggplot(filter(mod_pd, type == "plant")) +
#   geom_line(aes(age,PDz, group=factor(module, levels = mod_levels), 
#                 col=factor(module, levels = mod_levels))) +
#   geom_point(aes(age,PDz, group=factor(module, levels = mod_levels), 
#                  col=factor(module, levels = mod_levels)),
#              data = filter(mod_pd, PDp <= 0.05, type == "plant"),
#              size = 2) +
#   scale_color_manual(values = custom_pal,na.value = "grey70", drop = F) +
#   scale_x_reverse() +
#   labs(title = "Plants", color = "Module", x = "Millions of years ago, Ma", y = "PD z-score") +
#   theme_bw()

#+ pd, fig.width = 6.6, fig.height = 5.3, warning = FALSE
#gpdb / gpdp + plot_layout(guides = 'collect')
gpdb / gpdbw + plot_layout(guides = 'collect')


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
  labs(title = "Butterflies - Binary networks", color = "Module", x = "Millions of years ago, Ma", y = "Number of butterfly taxa") +
  theme_bw() +
  theme(legend.position = "none")

gnbw <- ggplot() + 
  geom_point(data=filter(wmod_pd, type == "butterfly"), 
             aes(age,ntaxa, col = factor(module, levels = wmod_levels)),
             alpha = 0.7) + 
  geom_line(data=filter(wmod_pd, type == "butterfly"),
            aes(age,ntaxa, col = factor(module, levels = wmod_levels))) + 
  #geom_point(data=filter(mod_pd, type == "plant"), aes(age,ntaxa, col = factor(module, levels = mod_levels)), shape = "square") + 
  scale_color_manual(values = custom_pal,na.value = "grey70", drop = F) +
  #scale_shape_manual(values = c("square","circle")) +
  scale_x_reverse() +
  labs(title = "Butterflies - Weighted networks", color = "Module", x = "Millions of years ago, Ma", y = "Number of butterfly taxa") +
  theme_bw() +
  theme(legend.position = "none")

gnp <- ggplot() + 
  geom_point(data=filter(mod_pd, type == "plant"), 
             aes(age,ntaxa, col = factor(module, levels = mod_levels)),
             alpha = 0.7) + 
  geom_line(data=filter(mod_pd, type == "plant"),
            aes(age,ntaxa, col = factor(module, levels = mod_levels))) + 
  scale_color_manual(values = custom_pal,na.value = "grey70", drop = F) +
  scale_x_reverse() +
  labs(title = "Plants - Binary networks", color = "Module", x = "Millions of years ago, Ma", y = "Number of plant taxa") +
  theme_bw() +
  theme(legend.position = "none")

gnpw <- ggplot() + 
  geom_point(data=filter(wmod_pd, type == "plant"), 
              aes(age,ntaxa, col = factor(module, levels = wmod_levels)),
              alpha = 0.7) + 
  geom_line(data=filter(wmod_pd, type == "plant"),
            aes(age,ntaxa, col = factor(module, levels = wmod_levels))) + 
  scale_color_manual(values = custom_pal,na.value = "grey70", drop = F) +
  scale_x_reverse() +
  labs(title = "Plants - Weighted networks", color = "Module", x = "Millions of years ago, Ma", y = "Number of plant taxa") +
  theme_bw()


#+ nmod, fig.width = 6.6, fig.height = 5.3, warning = FALSE
gnb + gnp + gnbw + gnpw + plot_layout(guides = 'collect')




/*
  
# _Density plots ----

# Binary networks
# modularity
for(g in ages[3:9]){
  
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

# nestedness
for(g in ages[3:9]){
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

#+ z_dens, fig.width = 11, fig.height = 3
(nodf_plot_60 / mod_plot_60 |
    nodf_plot_50 / mod_plot_50 |
    nodf_plot_40 / mod_plot_40 |
    nodf_plot_30 / mod_plot_30 |
    nodf_plot_20 / mod_plot_20 |
    nodf_plot_10 / mod_plot_10 |
    nodf_plot_0 / mod_plot_0) +
  plot_layout(guides = 'collect')


# Weighted networks
# modularity
for(g in ages[1:8]){
  
  gg <- ggplot(filter(Qwnull, age == g)) +
    geom_density(aes(Q, group = model, color = model)) +
    scale_color_grey() +
    geom_vline(xintercept = filter(Qwobs, age == g)$Q, col = "blue") +
    labs(title = paste0(g," Ma"), y = NULL) +
    theme_bw() + 
    theme(axis.ticks.y = element_line(linetype = "blank"), 
          axis.text.y = element_blank())
  
  assign(paste0("wmod_plot_",g), gg)
}

# nestedness
for(g in ages[1:8]){
  gg <- ggplot(filter(Nwnull, age == g)) +
    geom_density(aes(N, group = model, color = model)) +
    scale_color_grey() +
    geom_vline(xintercept = filter(Nwobs, age == g)$nodf, col = "blue") +
    labs(title = paste0(g," Ma"), y = NULL) +
    theme_bw() + 
    theme(axis.ticks.y = element_line(linetype = "blank"), 
          axis.text.y = element_blank())
  
  assign(paste0("wnodf_plot_",g), gg)
}

#+ z_dens, fig.width = 11, fig.height = 3
(wnodf_plot_80 / wmod_plot_80 |
    wnodf_plot_70 / wmod_plot_70 |
    wnodf_plot_60 / wmod_plot_60 |
    wnodf_plot_50 / wmod_plot_50 |
    wnodf_plot_40 / wmod_plot_40 |
    wnodf_plot_30 / wmod_plot_30 |
    wnodf_plot_20 / wmod_plot_20 |
    wnodf_plot_10 / wmod_plot_10) +
  plot_layout(guides = 'collect')



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

*/
  