Pieridae host repertoire - character history
================
Mariana Braga
01 May, 2020

-----

Script 2 for empirical study performed in Braga et al. 2020 *Evolution
of butterfly-plant networks revealed by Bayesian inference of host
repertoire*.

### Data

First we read in the phylogenetic trees for butterflies and plants.
Then, we read in the interaction matrix and remove plants (rows) that
are not hosts to any butterfly.

**Trees**

``` r
tree <- read.tree("./data/bphy_pie_ladder.phy")
host_tree <- read.tree("./data/angio_pie_50tips_ladder.phy")
```

**Extant network**

``` r
ext_net_50h <- as.matrix(read.csv("./data/incidence_pieridae.csv", header = T, row.names = 1))
identical(colnames(ext_net_50h), tree$tip.label)
identical(rownames(ext_net_50h), host_tree$tip.label)
```

``` r
ext_net <- ext_net_50h[which(rowSums(ext_net_50h) != 0),]
dim(ext_net)
```

    ## [1] 33 66

### Character history

**Read in .history.txt files**

These files can get quite big, so I can’t upload them in Github. Also,
you might want to thin out these files to speed up their parsing. This
is how you can do it. (Note that in the analyses in the paper we did not
thin the histories, so we’ll show results with all sampled histories).

``` r
colclasses <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))
```

``` r
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
```

The next time you run this script you just need to read in the thinned
histories. One file with sampled histories when using the
time-cklibrated host tree, and one using the transformed host tree (all
branch lengths = 1).

``` r
history_dat_time = read.table("./inference/history.time.txt", sep="\t", header=T, colClasses = colclasses)
history_dat_bl1 = read.table("./inference/history.bl1.txt", sep="\t", header=T, colClasses = colclasses)
```

**Calculate effective rate of evolution**

``` r
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
```

    ## [1] 0.09766287

``` r
(rate_bl1 <- n_events_bl1/tree_length)
```

    ## [1] 0.09510045

In both analyses, the rate of host repertoire evolution was near 1 event
every 10 Myr (along each branch).

**Ancestral networks**

First load the script `functions_ancestral_states.R` which contains the
functions to calculate the posterior probabilities for each
plant-butterfly interaction at given times in the past.

``` r
source("functions_ancestral_states.R")
```

Then choose the times in the past at which you want to reconstruct the
network.

``` r
ages = c(80,70,60,50,40,30,20,10,0)
```

This part is slow, so I won’t run it here.

``` r
list_m_at_ages = list()
for (i in 1:(length(ages)-1)) {
  age = ages[i]
  list_m_at_ages[[i]] = t(make_matrix_at_age( history_dat_time, age, s_hit=c(2) ))
}
```

I have done this before and saved the list as a .rds file.

``` r
list_m_at_ages <- readRDS("./inference/list_m_at_ages_time.rds")
  
length(list_m_at_ages)
```

    ## [1] 8

``` r
head(list_m_at_ages[[1]])
```

    ##                Pseudopontia    Index_70   Index_129
    ## Amborellaceae   0.002499306 0.003610108 0.001666204
    ## Cabombaceae     0.006942516 0.011108026 0.004443210
    ## Schisandraceae  0.015551236 0.021938350 0.010552624
    ## Annonaceae      0.021938350 0.030824771 0.020549847
    ## Siparunaceae    0.017772841 0.031380172 0.016106637
    ## Saururaceae     0.011385726 0.026103860 0.009997223

Then we add the extant network in the list with all ancestral networks.

``` r
list_m_at_ages[[9]] <- ext_net
```

We can build two kinds of ancestral networks, binary and quantitative.
For binary networks, we need to choose a probability threshold above
which we consider the interaction has enough support. Interactions above
the threshold are coded as 1, and the remaining interactions are coded
as 0.

**Binary networks with a probability threshold**

After defining the probability threshold, we transform the probabilities
into 1s and 0s.

``` r
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
```

    ##          Pseudopontia Index_70 Index_129
    ## Fabaceae            1        1         1

*Calculate modularity (stochastic step\!)*

This doesn’t work for the network at 80Ma because it only has one host,
so we add it manually before calculatig for the other networks.

``` r
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
```

![](Pieridae_histories_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
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
```
