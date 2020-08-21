library(ape)
library(dplyr)
require(stringr)
require(data.table)
library(vegan)
library(bipartite)

load("./Null_nets_samples.RData")

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

saveRDS(Nnull, file = "Nnull_samples.rds")
