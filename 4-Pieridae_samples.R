#'---
#'title: "Pieridae host repertoire - networks per sample"
#'author: "Mariana Braga"
#'date: "`r format(Sys.time(), '%d %B, %Y')`"
#'output: github_document
#'---

#'-------------
#'
#' Script 4 for empirical study performed in Braga et al. 2020
#' *Evolution of butterfly-plant networks revealed by Bayesian inference of host repertoire*.
#'

#+ include = FALSE

library(ape)
library(tidyverse)
library(patchwork)
library(viridis)
library(bipartite)
#library(ggraph)
#library(tidygraph)
#library(igraph)


/*# Get samples from the posterior ----
*/
  
#' ## Read in necessary files: trees and character history
#' 

tree <- read.tree("./data/bphy_pie_ladder.phy")
host_tree <- read.tree("./data/angio_pie_50tips_ladder.phy")

colclasses <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))
history_dat_bl1 = read.table("./inference/history.bl1.txt", sep="\t", header=T, colClasses = colclasses)

it_seq <- seq(20000,200000, 1800)  
history <- filter(history_dat_bl1, iteration %in% it_seq[-1])
nsamp <- length(unique(history$iteration))

ages = c(80,70,60,50,40,30,20,10,0)

source("functions_ancestral_states.R")

# list of arrays: sample (# iterations in history) x butterfly x plant
samples_ages = list()
for (i in 1:(length(ages)-1)) {
  age = ages[i]
  samples_ages[[i]] = make_matrix_samples_at_age( history, age, s_hit=c(2), tree, host_tree )
}


/*# Structure of sampled networks ----
*/
  
#' ### Structure of sampled networks
#' 

/*# _Modularity ----
*/

#+ eval = FALSE

source("my_compute_module.R")

# Calculate Q for observed networks (samples) AND save modules

Qsamples <- tibble()
Mod_samples <- tibble()

for(a in 1:(length(ages)-1)){
  for(i in 1:nsamp){
    net <- samples_ages[[a]][i,,]
    set.seed(2)                    ### remember to check if we set seed every time
    mod <- mycomputeModules(net)
    
    q <- mod@likelihood
    Qsamples <- bind_rows(Qsamples, tibble(age = ages[a], sample = i, Q=q))
    
    mod_list <- listModuleInformation(mod)[[2]]
    nmod <- length(mod_list)
    for(m in 1:nmod){
      members <- unlist(mod_list[[m]])
      mtbl <- tibble(name = members, 
                     age = rep(ages[a], length(members)),
                     sample = rep(i, length(members)),
                     original_module = rep(m, length(members)))
      
      Mod_samples <- bind_rows(Mod_samples, mtbl)
    }
  }
}

#saveRDS(Qsamples, "./networks/Qsamples.rds")
#saveRDS(Mod_samples, "./networks/Mod_samples.rds")

Qsamples <- readRDS("./networks/Qsamples.rds")        # tibble with Q values
Mod_samples <- readRDS("./networks/Mod_samples.rds")  # tibble with module information



/*# _Nestedness ----
*/

# Calculate NODF for observed networks (samples)

Nsamples <- tibble()

for(a in 1:(length(ages)-1)){
  for(i in 1:nsamp){
    net <- samples_ages[[a]][i,,]
    nodf <- networklevel(net, index="NODF")
    Nsamples <- bind_rows(Nsamples, tibble(age = ages[a], sample = i, N=nodf))
  }
}

#saveRDS(Nsamples, "./networks/Nsamples.rds")

Nsamples <- readRDS("./networks/Nsamples.rds")




/*# Null models ----
*/
  
#' ### Null models
#' 
#' `Nulls_age` are lists of length 100 (samples), where each element is an array of 100 null networks (buts x hosts x null)  

# number of null networks to be generated
nnull <- 100

#+ eval = FALSE  
Nulls_10 <- list()

for(i in 1:nsamp){
  net <- samples_ages[[8]][i,,]
  net = net[ rowSums(net)!=0, ]
  net = net[ ,colSums(net)!=0 ]
  
  null <- vegan::nullmodel(net, "r00") 
  sim <- simulate(null, nsim=nnull, seed = 1)
  Nulls_10[[i]] <- sim
  
}  


#' ### Nestedness and Modularity
#' 
#' Let's calculate the nestedness and modularity of each null network, per age
#' 
  
/*# _Nestedness ----
*/

#+ eval = FALSE

Nnull <- tibble()
for(a in 1:(length(ages)-1)){
  
  Nulls_age = paste0("Nulls_",ages[a])
  
  for(i in 1:nsamp){
    
    sim <- get(Nulls_age)[[i]]
    
    for(j in 1:nnull){
      Nrandom <- networklevel(sim[,,j],index="NODF")
      Nnull <- bind_rows(Nnull, tibble(age = ages[a], sample = i, sim=j, Nrandom=Nrandom))
    }
  }
}

#' Once the calculations are done in the server, read the resulting tibble
#' 

Nnull <- readRDS("./networks/Nnull_samples.rds")



/*# _Modularity ----
*/

#+ eval = FALSE

source("my_compute_module.R")

Qnull <- tibble()

for(a in 1:(length(ages)-1)){
  Nulls_age = paste0("Nulls_",ages[a])

  for(i in 1:nsamp){
    sim <- get(Nulls_age)[[i]]
    
    for(j in 1:nnull){
      Qrandom <- mycomputeModules(sim[,,j])@likelihood
      Qnull <- bind_rows(Qnull, tibble(age = ages[a], sample=i, sim=j, Qrandom=Qrandom))
    }
  }
}

saveRDS(Qnull, file = "Qnull_samples.rds") #(if everything was done in one job)

#' After running on the server, we can continue from here
#'

#(If ages are ran separately, 8 files will be produced)
Qnull_80 <- readRDS("./networks/Qnull_samples_80.rds")
Qnull_70 <- readRDS("./networks/Qnull_samples_70.rds")
Qnull_60 <- readRDS("./networks/Qnull_samples_60.rds")
Qnull_50 <- readRDS("./networks/Qnull_samples_50.rds")
Qnull_40 <- readRDS("./networks/Qnull_samples_40.rds")
Qnull_30 <- readRDS("./networks/Qnull_samples_30.rds")
Qnull_20 <- readRDS("./networks/Qnull_samples_20.rds")
Qnull_10 <- readRDS("./networks/Qnull_samples_10.rds")

Qnull <- bind_rows(Qnull_80,Qnull_70,Qnull_60,Qnull_50,
                   Qnull_40,Qnull_30,Qnull_20,Qnull_10)



/*# _Z-scores ----
*/
#' ### Z-scores
  
Qzscore <- Qnull %>% 
  group_by(age, sample) %>% 
  summarize(mean = mean(Qrandom),
            sd = sd(Qrandom)) %>% 
  left_join(Qsamples) %>% 
  mutate(z = (Q - mean)/sd) %>% 
  left_join(Qnull %>% 
              left_join(Qsamples) %>% 
              group_by(age, sample) %>% 
              summarise(p = sum(Qrandom >= Q)/nnull))


Nzscore <- Nnull %>% 
  group_by(age, sample) %>% 
  summarize(mean = mean(Nrandom),
            sd = sd(Nrandom)) %>% 
  left_join(Nsamples) %>% 
  mutate(z = (N - mean)/sd) %>% 
  left_join(Nnull %>% 
              left_join(Nsamples) %>% 
              group_by(age, sample) %>% 
              summarise(p = sum(Nrandom >= N)/nnull))

# 
# plot_qz <- ggplot(Qzscore) +
#   geom_point(aes(age, z), size = 1, col = "grey") +
#   geom_point(aes(age, z),
#              data = filter(Qzscore, p <= 0.05),
#              size = 2, alpha = 0.7) +
#   scale_x_reverse() +
#   labs(title = "Modularity, Q", y = "Z-score", x = "Millions of years ago, Ma") +
#   theme_bw()
# 
# plot_nz <- ggplot(Nzscore) +
#   geom_point(aes(age, z), size = 1, col = "grey") +
#   geom_point(aes(age, z),
#              data = filter(Nzscore, p <= 0.05),
#              size = 2, alpha = 0.7) +
#   scale_x_reverse() +
#   labs(title = "Nestedness, N", y = "Z-score", x = "Millions of years ago, Ma") +
#   theme_bw()
# 
# #+ zscore, fig.width = 6, fig.height = 4.5, warning = FALSE
# plot_qz / plot_nz + plot_layout(guides = 'collect')


#' Compare with summary networks (binary and weighted)
#' 

Nobs <- readRDS("./networks/Observed_Ns.rds")
Qobs <- readRDS("./networks/Observed_Qs.rds")

Nbnull <- readRDS("./networks/Nbnull.rds") %>%  filter(age != 0)
Qbnull <- readRDS("./networks/Qbnull.rds") %>%  filter(age != 0)

Nwnull <- readRDS("./networks/Nwnull.rds")
Qwnull <- readRDS("./networks/Qwnull.rds")


nit = 1000

# binary nets
Qbz <- Qbnull %>% 
  group_by(age, model) %>% 
  summarize(mean = mean(Q),
            sd = sd(Q)) %>% 
  left_join(Qobs) %>% 
  mutate(z = (Qbin - mean)/sd) %>% 
  left_join(Qbnull %>% 
              left_join(Qobs) %>% 
              group_by(age, model) %>% 
              summarise(p = sum(Q >= Qbin)/nit))


Nbz <- Nbnull %>% 
  group_by(age, model) %>% 
  summarize(mean = mean(N),
            sd = sd(N)) %>% 
  left_join(Nobs) %>% 
  mutate(z = (Nbin - mean)/sd) %>% 
  left_join(Nbnull %>% 
              left_join(Nobs) %>% 
              group_by(age, model) %>% 
              summarise(p = sum(N >= Nbin)/nit))


# weighted nets
Qwz <- Qwnull %>% 
  group_by(age, model) %>% 
  summarize(mean = mean(Q),
            sd = sd(Q)) %>% 
  left_join(Qobs) %>% 
  mutate(z = (Qwei - mean)/sd) %>% 
  left_join(Qwnull %>% 
              left_join(Qobs) %>% 
              group_by(age, model) %>% 
              summarise(p = sum(Q >= Qwei)/nit))


Nwz <- Nwnull %>% 
  group_by(age, model) %>% 
  summarize(mean = mean(N),
            sd = sd(N)) %>% 
  left_join(Nobs) %>% 
  mutate(z = (Nwei - mean)/sd) %>% 
  left_join(Nwnull %>% 
              left_join(Nobs) %>% 
              group_by(age, model) %>% 
              summarise(p = sum(N >= Nwei)/nit))

# merge to plot

# 2-color palette
pal_2c <- viridis_pal()(3)[2:1]

Qz <- bind_rows(filter(Qbz, model == 'random') %>% mutate(network = 'binary'), 
                filter(Qwz, model == 'random')%>% mutate(network = 'weighted'))

Nz <- bind_rows(filter(Nbz, model == 'random') %>% mutate(network = 'binary'), 
                filter(Nwz, model == 'random')%>% mutate(network = 'weighted'))

plot_qz <- ggplot(Qz) +
  geom_point(aes(age, z), col = "grey", alpha = 0.7, data = Qzscore) +
  geom_line(aes(age, z, group = network, col = network)) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Qz, p <= 0.05),
             size = 2, alpha = 0.7) +
  scale_color_manual(values = pal_2c) +
  scale_x_reverse() +
  labs(title = "Modularity, Q", y = "Z-score", x = "Millions of years ago, Ma") +
  theme_bw()

plot_nz <- ggplot(Nz) +
  geom_point(aes(age, z), col = "grey", alpha = 0.7, data = Nzscore) +
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


viol <- ggplot(Nz) +
  geom_violin(aes(age, z, group = age), col = "grey", data = Nzscore) +
  geom_line(aes(age, z, group = network, col = network)) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Nz, p <= 0.05),
             size = 2, alpha = 0.7) +
  scale_color_manual(values = pal_2c) +
  scale_x_reverse() +
  labs(title = "Nestedness, N", y = "Z-score", x = "Millions of years ago, Ma") +
  theme_bw()

point <- ggplot(Nz) +
  geom_point(aes(age, z), col = "grey", alpha = 0.7, data = Nzscore) +
  geom_line(aes(age, z, group = network, col = network)) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Nz, p <= 0.05),
             size = 2, alpha = 0.7) +
  scale_color_manual(values = pal_2c) +
  scale_x_reverse() +
  labs(title = "Nestedness, N", y = "Z-score", x = "Millions of years ago, Ma") +
  theme_bw()

viol / point

ggplot(Qzscore) +
  geom_density(aes(z)) +
  facet_grid(age ~ .) +
  theme_bw()


/*# _Pp of N and Q ----
*/
#' ### Posterior probability of network structure
#' 

ppN <- Nzscore %>% 
  group_by(age) %>% 
  filter(p <= 0.05) %>% 
  summarise(pp = n()/nsamp)
  
ppQ <- Qzscore %>% 
 group_by(age) %>% 
 filter(p <= 0.05) %>% 
 summarise(pp = n()/nsamp) 

plot(ppN)
plot(ppQ)

# > ppN
# # A tibble: 8 x 2
#     age    pp
#   <dbl> <dbl>
# 1    10 1    
# 2    20 1    
# 3    30 0.95 
# 4    40 0.580
# 5    50 0.13 
# 6    60 0.06 
# 7    70 0.03 
# 8    80 0.01 
# > ppQ
# # A tibble: 7 x 2
#     age    pp
#   <dbl> <dbl>
# 1    10  0.99
# 2    20  0.95
# 3    30  0.96
# 4    40  0.32
# 5    50  0.26
# 6    60  0.15
# 7    70  0.03


