Pieridae host repertoire - character history
================
Mariana Braga
29 April, 2020

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
you might want to thin out these files to speed up their parsing.

This is how you can do it.

``` r
colclasses <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))
# history_dat_time = read.table("./inference/out.2.real.pieridae.2s.history.txt", sep="\t", header=T, colClasses = colclasses)
# 
## define burnin and sampling interval  
#it_seq <- seq(20000,200000, 100)  
#
## Time-calibrated tree
#history_dat_time <- filter(history_dat_time, iteration %in% it_seq) %>% 
#  mutate(node_index = node_index + 1)
#
#write.table(history_dat_time,"./inference/history.time.txt", sep="\t", quote = F, row.names = F)
#
## Tree with branch lengths = 1
#history_dat_bl1 = read.table("./inference/out.3.bl1.pieridae.2s.history.txt", sep="\t", header=T, colClasses = colclasses)
#
#history_dat_bl1 <- filter(history_dat_bl1, iteration %in% it_seq) %>% 
#  mutate(node_index = node_index + 1)
#
#write.table(history_dat_bl1,"./inference/history.bl1.txt", sep="\t", quote = F, row.names = F)
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

Then choose the times in the past at which you want to reconstruct the
network.

``` r
#ages = c(80,70,60,50,40,30,20,10,0)
ages = c(80,70,0)

list_m_at_ages = list()
for (i in 1:(length(ages)-1)) {
  age = ages[i]
  list_m_at_ages[[i]] = t(make_matrix_at_age( history_dat_time, age, s_hit=c(2) ))
  #assign(paste0("m_",age), list_m_at_ages[[i]])
}

length(list_m_at_ages)
```

    ## [1] 2

``` r
head(list_m_at_ages[[1]])
```

    ##                Pseudopontia   Index_70   Index_129
    ## Amborellaceae   0.001665741 0.00388673 0.001665741
    ## Cabombaceae     0.006662965 0.01332593 0.003331483
    ## Schisandraceae  0.014436424 0.02276513 0.009439200
    ## Annonaceae      0.021654636 0.03109384 0.021654636
    ## Siparunaceae    0.017212660 0.03220433 0.018323154
    ## Saururaceae     0.010549695 0.02665186 0.008328706

Then we add the extant network in the list with all ancestral networks.

``` r
#list_m_at_ages[[9]] <- ext_net
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

head(net_list[[1]])
```

    ##             Pseudopontia Index_70 Index_129
    ## Capparaceae            0        0         1
    ## Fabaceae               1        1         1

**Weighted networks**

To build weighted networks we use the posterior probabilities as weights
for each interaction. But since many interactions have really small
probabilities, we can set a minimum probability, below which the weight
is set to 0.

``` r
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
```

    ## [1] 5 3

``` r
#dim(list_wnets[[9]])
```
