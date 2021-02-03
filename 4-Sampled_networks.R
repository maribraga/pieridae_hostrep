#'---
#'title: "Pieridae host repertoire - Sampled networks"
#'author: "Mariana Braga"
#'date: "`r format(Sys.time(), '%d %B, %Y')`"
#'output: github_document
#'---

#'-------------
#'
#' Script 4 for analyses performed in Braga et al. 2021
#' *Evolution of butterfly-plant networks over time, as revealed by Bayesian inference of host repertoire*.
#'

#+ include = FALSE

library(evolnets)
library(ape)
library(tidyverse)
library(patchwork)
library(bipartite)

#library(viridis)
#library(RColorBrewer)
#library(wesanderson)
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
  
Qzsamples <- Qnull %>% 
  group_by(age, sample) %>% 
  summarize(mean = mean(Qrandom),
            sd = sd(Qrandom)) %>% 
  left_join(Qsamples) %>% 
  mutate(z = (Q - mean)/sd) %>% 
  left_join(Qnull %>% 
              left_join(Qsamples) %>% 
              group_by(age, sample) %>% 
              summarise(p = sum(Qrandom >= Q)/nnull))


Nzsamples <- Nnull %>% 
  group_by(age, sample) %>% 
  summarize(mean = mean(Nrandom),
            sd = sd(Nrandom)) %>% 
  left_join(Nsamples) %>% 
  mutate(z = (N - mean)/sd) %>% 
  left_join(Nnull %>% 
              left_join(Nsamples) %>% 
              group_by(age, sample) %>% 
              summarise(p = sum(Nrandom >= N)/nnull))


/*# _Pp of N and Q ----
*/
#' ### Posterior probability of network structure
#' 

ppN <- Nzsamples %>% 
  group_by(age) %>% 
  filter(p <= 0.05) %>% 
  summarise(pp = n()/nsamp) %>% 
  mutate(y = 16)
  
ppQ <- Qzsamples %>% 
 group_by(age) %>% 
 filter(p <= 0.05) %>% 
 summarise(pp = n()/nsamp)  %>% 
  mutate(y = 13)


/*# __Plot Z-scores samples and consensus ----
*/
#' ### Plot Z-scores samples and consensus

# Read consensus z-scores

Qz <- readRDS("./networks/Qz.rds")
Nz <- readRDS("./networks/Nz.rds")

pal_3c <- brewer.pal(n = 9, name = 'Blues')[c(9,7,5)]
wes_orange <- wes_palette("Zissou1",7, type = "continuous")[6]

# first version #----------------------------------------------------#
plot_qz <- ggplot(Qz) +
  geom_point(aes(age, z), col = wes_orange, alpha = 0.7, data = Qzsamples) +
  stat_summary(aes(age, z), fun = "mean", geom = "line", col = wes_orange, data = Qzsamples) +
  geom_text(aes(age, y, label = pp), data = ppQ) +
  geom_line(aes(age, z, group = network, col = network)) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Qz, p <= 0.05),
             size = 2, alpha = 0.7) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Qz, p >= 0.95),
             size = 2, shape = 1) +
  scale_color_manual(values = pal_3c) +
  scale_x_reverse() +
  labs(title = "Modularity, Q", y = "Z-score", x = "Millions of years ago, Ma") +
  theme_bw()

plot_nz <- ggplot(Nz) +
  geom_point(aes(age, z), col = wes_orange, alpha = 0.7, data = Nzsamples) +
  stat_summary(aes(age, z), fun = "mean", geom = "line", col = wes_orange, data = Nzsamples) +
  geom_text(aes(age, y, label = pp), data = ppN) +
  geom_line(aes(age, z, group = network, col = network)) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Nz, p <= 0.05),
             size = 2, alpha = 0.7) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Nz, p >= 0.95),
             size = 2, shape = 1) +
  scale_color_manual(values = pal_3c) +
  scale_x_reverse() +
  labs(title = "Nestedness, N", y = "Z-score", x = "Millions of years ago, Ma") +
  theme_bw()

#+ zscore, fig.width = 6, fig.height = 4.5, warning = FALSE
plot_qz / plot_nz + plot_layout(guides = 'collect')
#----------------------------------------------------#

# second version #----------------------------------------------------#
plot_qz <- ggplot(Qz) +
  geom_point(aes(age, z), col = wes_orange, shape = 1, data = Qzsamples) +
  geom_point(aes(age, z), col = wes_orange, alpha = 0.8, data = filter(Qzsamples, p <= 0.05)) +
  stat_summary(aes(age, z), fun = "mean", geom = "line", col = wes_orange, data = Qzsamples) +
  geom_text(aes(age, y, label = pp), data = ppQ) +
  geom_line(aes(age, z, group = network, col = network)) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Qz, p <= 0.05),
             size = 2, alpha = 0.7) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Qz, p >= 0.95),
             size = 2, shape = 1) +
  scale_color_manual(values = pal_3c) +
  scale_x_reverse() +
  labs(title = "Modularity, Q", y = "Z-score", x = "Millions of years ago, Ma") +
  theme_bw()

plot_nz <- ggplot(Nz) +
  geom_point(aes(age, z), col = wes_orange, alpha = 0.7, data = Nzsamples) +
  stat_summary(aes(age, z), fun = "mean", geom = "line", col = wes_orange, data = Nzsamples) +
  geom_text(aes(age, y, label = pp), data = ppN) +
  geom_line(aes(age, z, group = network, col = network)) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Nz, p <= 0.05),
             size = 2, alpha = 0.7) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Nz, p >= 0.95),
             size = 2, shape = 1) +
  scale_color_manual(values = pal_3c) +
  scale_x_reverse() +
  labs(title = "Nestedness, N", y = "Z-score", x = "Millions of years ago, Ma") +
  theme_bw()

#+ zscore, fig.width = 6, fig.height = 4.5, warning = FALSE
plot_qz / plot_nz + plot_layout(guides = 'collect')

#----------------------------------------------------#

# Violin plot

violN <- ggplot(Nz) +
  geom_violin(aes(age, z, group = age), col = "white", fill = wes_orange, alpha = 0.5, data = Nzsamples) +
  stat_summary(aes(age, z), fun = "mean", geom = "line", col = wes_orange, data = Nzsamples) +
  stat_summary(aes(age, z), fun = "mean", geom = "point", col = wes_orange, data = filter(Nzsamples, age < 40)) +
  geom_text(aes(age, y, label = pp), data = ppN) +
  geom_line(aes(age, z, group = network, col = network)) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Nz, p <= 0.05),
             size = 2, alpha = 0.7) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Nz, p >= 0.95),
             size = 2, shape = 1) +
  scale_color_manual(values = pal_3c) +
  scale_x_reverse() +
  labs(title = "Nestedness, N", y = "Z-score", x = "Millions of years ago, Ma") +
  theme_bw()

violQ <- ggplot(Qz) +
  geom_violin(aes(age, z, group = age), col = "white", fill = wes_orange, alpha = 0.5, data = Qzsamples) +
  stat_summary(aes(age, z), fun = "mean", geom = "line", col = wes_orange, data = Qzsamples) +
  stat_summary(aes(age, z), fun = "mean", geom = "point", col = wes_orange, data = filter(Qzsamples, age < 40)) +
  geom_text(aes(age, y, label = pp), data = ppQ) +
  geom_line(aes(age, z, group = network, col = network)) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Qz, p <= 0.05),
             size = 2, alpha = 0.7) +
  geom_point(aes(age, z, group = network, col = network),
             data = filter(Qz, p >= 0.95),
             size = 2, shape = 1) +
  scale_color_manual(values = pal_3c) +
  scale_x_reverse() +
  labs(title = "Modularity, Q", y = "Z-score", x = "Millions of years ago, Ma") +
  theme_bw()

violQ / violN

ggplot(Qzsamples) +
  geom_density(aes(z)) +
  facet_grid(age ~ .) +
  theme_bw()


# Figure S3 - raw Q and NODF ------------

Qs <- bind_rows(Qobs %>% mutate(network = '0.9'),
                Qwobs50 %>% mutate(network = '0.5'),
                Qwobs %>% mutate(network = '0.1'))
Ns <- bind_rows(Nobs %>% mutate(network = '0.9'),
                Nwobs50 %>% mutate(network = '0.5'),
                Nwobs %>% mutate(network = '0.1'))


rawN <- ggplot() +
  geom_violin(aes(age, N, group = age), col = "white", fill = wes_orange, alpha = 0.5, data = Nsamples) +
  stat_summary(aes(age, N), fun = "mean", geom = "point", col = wes_orange, data = Nsamples) +
  geom_point(aes(age, nodf, group = network, col = network), data = filter(Ns, age > 0)) +
  geom_point(aes(age, nodf), data = filter(Ns, age == 0)) +
  scale_color_manual(values = pal_3c) +
  scale_x_reverse() +
  labs(title = "Nestedness, NODF", y = "NODF", x = "Millions of years ago, Ma") +
  theme_bw()

rawQ <- ggplot(Qs) +
  geom_violin(aes(age, Q, group = age), col = "white", fill = wes_orange, alpha = 0.5, data = Qsamples) +
  stat_summary(aes(age, Q), fun = "mean", geom = "point", col = wes_orange, data = Qsamples) +
  geom_point(aes(age, Q, group = network, col = network), data = filter(Qs, age > 0)) +
  geom_point(aes(age, Q), data = filter(Qs, age == 0)) +
  scale_color_manual(values = pal_3c) +
  scale_x_reverse() +
  labs(title = "Modularity, Q", y = "Q", x = "Millions of years ago, Ma") +
  theme_bw()

rawQ / rawN




/*# Pairwise module membership ----
*/
#' ### Pairwise module membership
#' 

Mod_samples <- readRDS("./networks/Mod_samples.rds")  # tibble with module information

ages = c(80,70,60,50,40,30,20,10)
nsamp = 100

pair_mod_matrix <- list()
pair_mod_tbl <- list()

# Find sampled networks at 80Ma that got several modules with the same hosts (connectance = 1)
out <- Mod_samples %>% group_by(age, sample) %>% distinct(name) %>% summarize(u=n()) %>% 
  left_join(Mod_samples %>% group_by(age, sample) %>% summarize(n=n())) %>% 
  mutate(problem = case_when(u != n ~ "YES", u == n ~ "NO")) %>% 
  filter(problem == "YES")

good_samples_at_80 <- setdiff(1:100,out$sample)

for(a in 1:length(ages)){
  taxa = Mod_samples %>% filter(age == ages[a]) %>% distinct(name)
  ntaxa = nrow(taxa)
  
  heat <- matrix(data=0, nrow = ntaxa, ncol = ntaxa)
  rownames(heat) <- colnames(heat) <- taxa$name
  
  tbl <- tibble(row = taxa$name, col = taxa$name) %>% 
    complete(row,col) %>% 
    mutate(freq = 0)
  
  for(i in 1:nsamp){
    
    if(a == 1){
      if(i %in% good_samples_at_80){
        table <- Mod_samples %>% filter(age == ages[a], sample == i)
        mods <- unique(table$original_module)
        
        for(m in 1:length(mods)){
          module <- filter(table, original_module == m)
          
          for(r in 1:nrow(heat)){
            for(c in 1:ncol(heat)){
              if(rownames(heat)[r] %in% module$name & colnames(heat)[c] %in% module$name){
                heat[r,c] = heat[r,c] + 1
                tbl <- mutate(tbl, freq = case_when(row == rownames(heat)[r] & col == colnames(heat)[c] ~ freq + 1,
                                                    TRUE ~ freq))
              }
            }
          }
        }
      }
      
    } else{
      table <- Mod_samples %>% filter(age == ages[a], sample == i)
      mods <- unique(table$original_module)
      
      for(m in 1:length(mods)){
        module <- filter(table, original_module == m)
        
        for(r in 1:nrow(heat)){
          for(c in 1:ncol(heat)){
            if(rownames(heat)[r] %in% module$name & colnames(heat)[c] %in% module$name){
              heat[r,c] = heat[r,c] + 1
              tbl <- mutate(tbl, freq = case_when(row == rownames(heat)[r] & col == colnames(heat)[c] ~ freq + 1,
                                                  TRUE ~ freq))
            }
          }
        }
      }
    }
  }
  pair_mod_matrix[[a]] <- heat
  pair_mod_tbl[[a]] <- tbl
}





#saveRDS(pair_mod_tbl,"./networks/pair_mod_tbl.rds")
#pair_mod_tbl <- readRDS("./networks/pair_mod_tbl.rds")

# plot with base R heatmap()
#heatmap(pair_mod_matrix[[3]], cexRow = 0.5, cexCol = 0.5)
#heatmap(pair_mod_matrix[[5]], cexRow = 0.5, cexCol = 0.5)
#heatmap(pair_mod_matrix[[5]], Rowv = NA, Colv = NA, scale = "column", cexRow = 0.5, cexCol = 0.5)


# Add module information to plot with ggplot2

all_wmod50_edited <- read.csv("./networks/all_wmod50_bl1.csv", header = T, stringsAsFactors = F)

wmod_levels50 <- c(paste0('M',1:12))
custom_palw50 <- c("#8a1c4c","#1b1581","#e34c5b","#fca33a","#fbeba9","#fdc486",
                 "#b370a8","#f8c4cc","#c8d9ee","#82a0be","#00a2bf","#006e82")


# Loop across all ages

for(i in 1:length(ages)){
  
  edge_list <- as.tbl(pair_mod_tbl[[i]]) %>%
  #left_join(all_wmod50_edited %>% filter(age == ages[i]) %>% select(name, module), by = c("row" = "name")) %>%
  #left_join(all_wmod50_edited %>% filter(age == ages[i]) %>% select(name, module), by = c("col" = "name")) %>%
  inner_join(all_wmod50_edited %>% filter(age == ages[i]) %>% select(name, module), by = c("row" = "name")) %>%
  inner_join(all_wmod50_edited %>% filter(age == ages[i]) %>% select(name, module), by = c("col" = "name")) %>%
  mutate(Module = ifelse(module.x == module.y, module.x, NA),
         Frequency = freq/max(freq)) 

  order <- edge_list %>% arrange(module.x) %>% pull(row) %>% unique()

  plot_net <- edge_list %>% mutate(
   row = factor(row, levels = order),
   col = factor(col, levels = order))
  
  gg <- ggplot(plot_net, aes(x = row, y = reorder(col,desc(col)), fill = factor(Module, levels = wmod_levels50), alpha = Frequency)) +
  geom_tile() +
  theme_bw() +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_manual(values = custom_palw50, na.value = "grey20", drop = F) +
  scale_alpha(range = c(min(edge_list$Frequency),max(edge_list$Frequency))) +
  labs(fill = "Module") +
  theme(
    #axis.text.x = element_text(angle = 270, hjust = 0, size = 6),
    #axis.text.y = element_text(size = 6),
    axis.text = element_blank(),     
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())
  
  assign(paste0("heat_mod50_",ages[[i]]), gg)
}


# define layout
design <- c(patchwork::area(1,1,1,1),
            patchwork::area(2,1,2,1),
            patchwork::area(3,1,3,1),
            patchwork::area(4,1,4,1),
            patchwork::area(1,2,2,3),
            patchwork::area(3,2,4,3),
            patchwork::area(1,4,2,5),
            patchwork::area(3,4,4,5))

# design <- c(patchwork::area(1,1,1,1),
#             patchwork::area(1,2,1,2),
#             patchwork::area(1,3,1,3),
#             patchwork::area(1,4,1,4),
#             patchwork::area(2,1,3,2),
#             patchwork::area(2,3,3,4),
#             patchwork::area(4,1,6,2),
#             patchwork::area(4,3,6,4))

#+ fig3, fig.width = 20, fig.height = 15, warning = F
# plot!
heat_mod50_80 +
  heat_mod50_70 + 
  heat_mod50_60 + 
  heat_mod50_50 + 
  heat_mod50_40 + 
  heat_mod50_30 + 
  heat_mod50_20 + 
  heat_mod50_10 + 
  plot_layout(guides = 'collect', design = design)

