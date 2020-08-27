library(ape)
library(dplyr)
require(stringr)
require(data.table)
library(vegan)
library(bipartite)

load("./Null_nets_samples.RData")
source("my_compute_module.R")


Qnull <- tibble()

a = 8
Nulls_age = paste0("Nulls_",ages[a])

for(i in 1:nsamp){

  sim <- get(Nulls_age)[[i]]

  for(j in 1:nnull){
    set.seed(2)         ### forgot to add this before
    Qrandom <- mycomputeModules(sim[,,j])@likelihood
    Qnull <- bind_rows(Qnull, tibble(age = ages[a], sample=i, sim=j, Qrandom=Qrandom))
  }
}

saveRDS(Qnull, file = paste0("Qnull_samples_",ages[a],".rds"))

