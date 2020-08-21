
library(ggtree)
library(tidyverse)
library(patchwork)
library(ggraph)
library(tidygraph)
library(bipartite)
library(igraph)


tree <- read.tree("./data/bphy_pie_ladder.phy")
host_tree <- read.tree("./data/angio_pie_50tips_ladder.phy")

colclasses <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))
history_dat_bl1 = read.table("./inference/history.bl1.txt", sep="\t", header=T, colClasses = colclasses)

it_seq <- seq(20000,200000, 1800)  
history <- filter(history_dat_bl1, iteration %in% it_seq[-1])
length(unique(history$iteration))

ages = c(80,70,60,50,40,30,20,10,0)

source("functions_ancestral_states.R")

# list of arrays: sample (# iterations in history) x butterfly x plant
samples_ages = list()
for (i in 1:(length(ages)-1)) {
  age = ages[i]
  samples_ages[[i]] = make_matrix_samples_at_age( history, age, s_hit=c(2), tree, host_tree )
}


# Modularity

nodes_mod <- tibble()
qsamples <- tibble()

nsamples <- 100
for (a in 1:(length(ages)-1)) {
  for(i in 1:nsamples){
    matrix <- samples_ages[[a]][i,,]
    matrix = matrix[ rowSums(matrix)!=0, ]
    matrix = matrix[ ,colSums(matrix)!=0 ]
    set.seed(5)
    lpa <- DIRT_LPA_wb_plus(matrix, mini=4, reps=50)
    mod <- convert2moduleWeb(matrix, lpa)
    assign(paste0("mod_",ages[a],"_",i),mod)
    mod_list <- listModuleInformation(mod)[[2]]
    nmod <- length(mod_list)
    for(m in 1:nmod){
      members <- unlist(mod_list[[m]])
      mtbl <- tibble(name = members, 
                     age = rep(ages[a], length(members)),
                     sample = rep(i, length(members)),
                     original_module = rep(m, length(members)))
      
      nodes_mod <- bind_rows(nodes_mod, mtbl)
      
      qtbl <- tibble(age = ages[a],
                     sample = i,
                     q = mod@likelihood)
      qsamples <- bind_rows(qsamples, qtbl)
    }
  }
}

saveRDS(qsamples, "qsamples.rds")

Qobs <- readRDS("./networks/Observed_Qs.rds")

plot_qsamples <- qsamples %>% 
  left_join(Qobs, by = "age")

#pal_2c <- viridis_pal()(3)[2:1]

ggQ_samples <- ggplot(plot_qsamples) +
  geom_violin(aes(age,q, group = age), col = "grey") +
  geom_point(aes(age,Qbin), col = "#21908CFF") +
  geom_point(aes(age,Qwei), col = "#440154FF") +
  labs(title = "Modularity per sample, binary (green), weighted (purple)",
       x = "Age (time before present in million years)",
       y = "Modularity") +
  scale_x_reverse() +
  theme_bw()


plotModuleWeb(mod_60_1, labsize = 0.4)
plotModuleWeb(mod_60_2, labsize = 0.4)
plotModuleWeb(mod_60_3, labsize = 0.4)

plotModuleWeb(mod_10_1, labsize = 0.4)
plotModuleWeb(mod_10_2, labsize = 0.4)
plotModuleWeb(mod_10_3, labsize = 0.4)


# Nestedness

Nsamples <- tibble()

nsamples <- 100
for (a in 1:(length(ages)-1)) {
  for(i in 1:nsamples){
    matrix <- samples_ages[[a]][i,,]
    matrix = matrix[ rowSums(matrix)!=0, ]
    matrix = matrix[ ,colSums(matrix)!=0 ]
    nodf <- networklevel(matrix,index="NODF")
    
    ntbl <- tibble(age = ages[a],
                   sample = i,
                   n = nodf)
    Nsamples <- bind_rows(Nsamples, ntbl)
    }
  }
}

saveRDS(Nsamples, "Nsamples.rds")

Nobs <- readRDS("./networks/Observed_Ns.rds")

plot_Nsamples <- Nsamples %>% 
  left_join(Nobs, by = "age")

#pal_2c <- viridis_pal()(3)[2:1]

ggN_samples <- ggplot(plot_Nsamples) +
  geom_violin(aes(age,n, group = age), col = "grey") +
  geom_point(aes(age,Nbin), col = "#21908CFF") +
  geom_point(aes(age,Nwei), col = "#440154FF") +
  labs(title = "Nestedness per sample, binary (green), weighted (purple)",
       x = "Age (time before present in million years)",
       y = "Nestedness") +
  scale_x_reverse() +
  theme_bw()


ggQ_samples / ggN_samples


# wgraph <- as_tbl_graph(wnet, directed = F) %>% 
#   left_join(filter(all_mod_edited, age == ages[n]), by = "name") %>% 
#   select(type, name, module)





# Maximizing modularity ----

# increasing the number of reps per number of modules

ext_net_50h <- as.matrix(read.csv("./data/incidence_pieridae.csv", header = T, row.names = 1))
ext_net <- ext_net_50h[which(rowSums(ext_net_50h) != 0),]

list_m_at_ages <- readRDS("./inference/list_m_at_ages_bl1.rds")
length(list_m_at_ages)
list_m_at_ages[[9]] <- ext_net

pt <- 90 
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

net30 <- as.matrix(net_list[[6]])
net10 <- as.matrix(net_list[[8]])
  
set.seed(5)
lpa <- DIRT_LPA_wb_plus(ext_net, mini=8, reps=100)
mdef <- convert2moduleWeb(ext_net, lpa)
plotModuleWeb(mdef, labsize = 0.4)
mdef@likelihood
#0.6432
length(unique(lpa[[1]]))
# 13 modules

set.seed(5)
lpa4 <- DIRT_LPA_wb_plus(ext_net, mini=4, reps=100)
mdef4 <- convert2moduleWeb(ext_net, lpa4)
plotModuleWeb(mdef4, labsize = 0.4)
mdef4@likelihood
#0.647616
length(unique(lpa4[[1]]))
# 10 modules

A = LPA_wb_plus(ext_net) # 26 modules
qlpa <- convert2moduleWeb(ext_net, A)
plotModuleWeb(qlpa, labsize = 0.4)

aa = 13
B = LPA_wb_plus(ext_net, aa)
lpab <- convert2moduleWeb(ext_net, B)
plotModuleWeb(lpab, labsize = 0.4)
lpab@likelihood

#-----------------------------------------------#
function (MATRIX, mini = 4, reps = 10) 
{
  A = LPA_wb_plus(MATRIX)
  mods = length(unique(A[[1]]))
  if ((mods - mini) > 0) {
    for (aa in mini:mods) {
      for (bb in 1:reps) {
        B = LPA_wb_plus(MATRIX, aa)
        if (B[[3]] > A[[3]]) 
          A = B
      }
    }
  }
  return(list(Row_labels = A[[1]], Col_labels = A[[2]], modularity = A[[3]]))
}
#-----------------------------------------------#



set.seed(5)
mod <- computeModules(net30)
plotModuleWeb(mod, labsize = 0.4)
mod@likelihood

set.seed(5)
l30 <- DIRT_LPA_wb_plus(net30, mini=4, reps=100)
m30 <- convert2moduleWeb(net30, l30)
plotModuleWeb(m30, labsize = 0.4)
m30@likelihood
#0.59875
length(unique(l30[[1]]))
# 4 modules


# 500 reps -- different seeds still don't converge

set.seed(5)
lpa500_5 <- DIRT_LPA_wb_plus(ext_net, mini=4, reps=500)
m500_5 <- convert2moduleWeb(ext_net, lpa500_5)
plotModuleWeb(m500_5, labsize = 0.4)
m500_5@likelihood
# 0.64864
length(unique(lpa500_5[[1]]))
# 8 modules

set.seed(2)
lpa500_2 <- DIRT_LPA_wb_plus(ext_net, mini=4, reps=500)
m500_2 <- convert2moduleWeb(ext_net, lpa500_2)
plotModuleWeb(m500_2, labsize = 0.4)
m500_2@likelihood
# 0.64832
length(unique(lpa500_2[[1]]))
# 10 modules

