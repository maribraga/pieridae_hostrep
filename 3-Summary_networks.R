#'---
#'title: "Pieridae host repertoire - Structure of summary networks"
#'author: "Mariana Braga"
#'date: "`r format(Sys.time(), '%d %B, %Y')`"
#'output: github_document
#'---

#'-------------
#'
#' Script 3 for analyses performed in Braga et al. 2021 Ecology Letters
#' *Phylogenetic reconstruction of ancestral ecological networks through time for pierid butterflies and their host plants*.
#'
#' This is a continuation of script 2 - Character history, so make sure you complete that one first.

#+ include = FALSE
load("./inference/char_hist_complete.RData")

library(evolnets)
library(tidyverse)
library(patchwork)
library(bipartite)
library(ggraph)
library(tidygraph)
#library(igraph)


#' In this script we will investigate how modularity and nestedness in the interactions
#' between Pieridae and their host plants changed over time. For that we will calculate
#' z-scores for Q and NODF. 
#' 

/*# Extant network ----
*/
#' ### Extant network
#' 
#' First, we'll calculate the z-scores for the extant network
#' 
  
/*# _Null models ----
*/
#' **Null models**
#'

# number of null networks  
nit <- 1000

#+ eval = FALSE
null_model <- vegan::nullmodel(ext_net, "r00") 
Nulls_ext <- simulate(null_model, nsim=nit, seed = 1)

/*# __Modularity ----
*/
#' **Modularity**
#' 

#+ eval = FALSE
# this takes a lot of time!
Qnull_ext <- tibble()

for(j in 1:nit){
  Qrandom <- mycomputeModules(Nulls_ext[,,j])@likelihood
  Qnull_ext <- bind_rows(Qnull_ext, tibble(age = 0, sim=j, Q=Qrandom))
}

/*
saveRDS(Qnull_ext, "./networks/Qnull_ext.rds")
*/

#+
# quick way
# Null distribution of modularity for extant network
Qnull_ext <- readRDS("./networks/Qnull_ext.rds")

# Observed modularity in the extant network
Qobs_ext <- mycomputeModules(ext_net)@likelihood


/*# __Nestedness ----
*/
#' **Nestedness**
#' 
  
#+ eval = FALSE
# slow way
Nnull_ext <- tibble()
 
for(j in 1:nit){
  Nrandom <- networklevel(Nulls_ext[,,j],index="NODF")
  Nnull_ext <- bind_rows(Nnull_ext, tibble(age = 0, sim=j, NODF=Nrandom))
}

/*
saveRDS(Nnull_ext, "./networks/Nnull_ext.rds")
*/
  
#+ 
# quick way
# Null distribution of nestedness for extant network
Nnull_ext <- readRDS("./networks/Nnull_ext.rds")

# Observed nestedness in the extant network
Nobs_ext <- networklevel(ext_net,index="NODF")


/*# _Z-scores ----
*/
#' **Z-scores and p**
  
Qzscore_ext <- Qnull_ext %>% 
  summarize(mean = mean(Q),
            sd = sd(Q)) %>% 
  mutate(z = (Qobs_ext - mean)/sd) 

Qp_ext <- Qnull_ext %>% 
  summarise(p = sum(Q > Qobs_ext)/nit) %>% 
  pull()


Nzscore_ext <- Nnull_ext %>%
  summarize(mean = mean(NODF),
            sd = sd(NODF)) %>% 
  mutate(z = (Nobs_ext - mean)/sd) 


Np_ext <- Nnull_ext %>% 
  summarise(p = sum(NODF > Nobs_ext)/nit) %>% 
  pull()

  
  
/*# Ancestral networks ----
*/
#' ### Ancestral networks
#' 
#' Now let's calculate the z-scores for the ancestral networks. 
#' I'll go through all the steps with the `weighted_net_50` network. For binary networks,
#' follow the steps done for the extant network above.
#' 
  
/*# _Null models ----
*/
#' **Null models**
#'
#' We will generate 1000 null networks that will be used to produce null distributions
#' for Q and NODF for each summary network at each age.

nit <- 1000
Nulls <- list()

for(i in 1:(length(ages)-1)){
  
  count <- round(weighted_net_50[[i]]*100) # transform the probabilities into counts 
                                           # to use the null model 'r00_both'
  null_rb <- vegan::nullmodel(count, "r00_both") 
  sim_rb <- simulate(null_rb, nsim=nit, seed = 1)
  Nulls[[i]] <- sim_rb
  
} 

/*
saveRDS(Nulls, "./networks/Nulls_50pp.rds")
*/
  


/*# _Modularity ----
*/
#' **Modularity**
#' 
#' Now we calculate modularity for each null network

#+ eval = FALSE
# slow! takes hours
Qwnull50 <- tibble()

for(i in 1:(length(ages)-1)){
  
  simw_random <- Nulls[[i]]

  for(j in 1:nit){
    Qwrandom <- computeModules(simw_random[,,j])@likelihood
    Qwnull50 <- bind_rows(Qwnull50, tibble(age = ages[i], sim=j, Q=Qwrandom))
  }
}

/*
saveRDS(Qwnull50, "./networks/Qwnull50.rds")
*/
  
#+
# Null distribution of modularity for the summary network with probability threshold of 0.5 
Qwnull50 <- readRDS("./networks/Qwnull50.rds")

# Observed modularity (uses an object produced in script 2)
Qwobs50 <- tibble()

for(i in 1:(length(ages)-1)){
  q <- get(paste0("wmod50_",ages[i]))
  Qwobs50<- bind_rows(Qwobs50, tibble(age = ages[i], Q = q@likelihood))
}


/*# _Nestedness ----
*/
#' **Nestedness**
#'   

#+ eval = FALSE
Nwnull50 <- tibble()

for(i in 1:(length(ages)-1)){
  
  sim_rb <- Nulls[[i]]

  for(j in 1:nit){
    Nrb <- networklevel(sim_rb[,,j],index="weighted NODF")  # important difference between binary and weighted networks
    Nwnull50 <- bind_rows(Nwnull50, tibble(age = ages[i], sim=j, NODF=Nrb))
  }
}

/*
saveRDS(Nwnull50, "./networks/Nwnull50.rds")
*/
  
#+ 
# Nestedness of weighted null networks
Nwnull50 <- readRDS("./networks/Nwnull50.rds")

# observed
Nwobs50 <- tibble()

for(i in 1:(length(ages)-1)){
  wnodf <- networklevel(weighted_net_50[[i]],index="weighted NODF")
  Nwobs50<- bind_rows(Nwobs50, tibble(age = ages[i], NODF = wnodf))
}

    
/*# _Z-scores ----
*/
#' **Z-scores**
#' 

#+ message = FALSE      
Qwzscore50 <- Qwnull50 %>%
  group_by(age) %>%
  summarize(mean = mean(Q),
            sd = sd(Q)) %>%
  left_join(Qwobs50) %>%
  mutate(z = (Q - mean)/sd) %>%
  left_join(Qwnull50 %>%
              left_join(rename(Qwobs50, Qobs = Q)) %>%
              group_by(age) %>%
              summarise(p = sum(Q > Qobs)/nit))

Nwzscore50 <- Nwnull50 %>% 
  group_by(age) %>% 
  summarize(mean = mean(NODF),
            sd = sd(NODF)) %>% 
  left_join(Nwobs50) %>% 
  mutate(z = (NODF - mean)/sd) %>% 
  left_join(Nwnull50 %>% 
              left_join(rename(Nwobs50, Nobs = NODF)) %>% 
              group_by(age) %>% 
              summarise(p = sum(NODF > Nobs)/nit))

