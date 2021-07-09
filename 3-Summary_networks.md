Pieridae host repertoire - Structure of summary networks
================
Mariana Braga
09 July, 2021

------------------------------------------------------------------------

Script 3 for analyses performed in Braga et al. 2021 Ecology Letters
*Phylogenetic reconstruction of ancestral ecological networks through
time for pierid butterflies and their host plants*.

This is a continuation of script 2 - Character history, so make sure you
complete that one first.

In this script we will investigate how modularity and nestedness in the
interactions between Pieridae and their host plants changed over time.
For that we will calculate z-scores for Q and NODF.

### Extant network

First, we’ll calculate the z-scores for the extant network

**Null models**

``` r
# number of null networks  
nit <- 1000
```

``` r
null_model <- vegan::nullmodel(ext_net, "r00") 
Nulls_ext <- simulate(null_model, nsim=nit, seed = 1)
```

**Modularity**

``` r
# this takes a lot of time!
Qnull_ext <- tibble()

for(j in 1:nit){
  Qrandom <- mycomputeModules(Nulls_ext[,,j])@likelihood
  Qnull_ext <- bind_rows(Qnull_ext, tibble(age = 0, sim=j, Q=Qrandom))
}
```

``` r
# quick way
# Null distribution of modularity for extant network
Qnull_ext <- readRDS("./networks/Qnull_ext.rds")

# Observed modularity in the extant network
Qobs_ext <- mycomputeModules(ext_net)@likelihood
```

**Nestedness**

``` r
# slow way
Nnull_ext <- tibble()
 
for(j in 1:nit){
  Nrandom <- networklevel(Nulls_ext[,,j],index="NODF")
  Nnull_ext <- bind_rows(Nnull_ext, tibble(age = 0, sim=j, NODF=Nrandom))
}
```

``` r
# quick way
# Null distribution of nestedness for extant network
Nnull_ext <- readRDS("./networks/Nnull_ext.rds")

# Observed nestedness in the extant network
Nobs_ext <- networklevel(ext_net,index="NODF")
```

**Z-scores and p**

``` r
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
```

### Ancestral networks

Now let’s calculate the z-scores for the ancestral networks. I’ll go
through all the steps with the `weighted_net_50` network. For binary
networks, follow the steps done for the extant network above.

**Null models**

We will generate 1000 null networks that will be used to produce null
distributions for Q and NODF for each summary network at each age.

``` r
nit <- 1000
Nulls <- list()

for(i in 1:(length(ages)-1)){
  
  count <- round(weighted_net_50[[i]]*100) # transform the probabilities into counts 
                                           # to use the null model 'r00_both'
  null_rb <- vegan::nullmodel(count, "r00_both") 
  sim_rb <- simulate(null_rb, nsim=nit, seed = 1)
  Nulls[[i]] <- sim_rb
  
} 
```

**Modularity**

Now we calculate modularity for each null network

``` r
# slow! takes hours
Qwnull50 <- tibble()

for(i in 1:(length(ages)-1)){
  
  simw_random <- Nulls[[i]]

  for(j in 1:nit){
    Qwrandom <- computeModules(simw_random[,,j])@likelihood
    Qwnull50 <- bind_rows(Qwnull50, tibble(age = ages[i], sim=j, Q=Qwrandom))
  }
}
```

``` r
# Null distribution of modularity for the summary network with probability threshold of 0.5 
Qwnull50 <- readRDS("./networks/Qwnull50.rds")

# Observed modularity (uses an object produced in script 2)
Qwobs50 <- tibble()

for(i in 1:(length(ages)-1)){
  q <- get(paste0("wmod50_",ages[i]))
  Qwobs50<- bind_rows(Qwobs50, tibble(age = ages[i], Q = q@likelihood))
}
```

**Nestedness**

``` r
Nwnull50 <- tibble()

for(i in 1:(length(ages)-1)){
  
  sim_rb <- Nulls[[i]]

  for(j in 1:nit){
    Nrb <- networklevel(sim_rb[,,j],index="weighted NODF")  # important difference between binary and weighted networks
    Nwnull50 <- bind_rows(Nwnull50, tibble(age = ages[i], sim=j, NODF=Nrb))
  }
}
```

``` r
# Nestedness of weighted null networks
Nwnull50 <- readRDS("./networks/Nwnull50.rds")

# observed
Nwobs50 <- tibble()

for(i in 1:(length(ages)-1)){
  wnodf <- networklevel(weighted_net_50[[i]],index="weighted NODF")
  Nwobs50<- bind_rows(Nwobs50, tibble(age = ages[i], NODF = wnodf))
}
```

**Z-scores**

``` r
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


Qwzscore50
```

    ## # A tibble: 8 x 6
    ##     age   mean      sd      Q      z     p
    ##   <dbl>  <dbl>   <dbl>  <dbl>  <dbl> <dbl>
    ## 1    10 0.523  0.0142  0.689  11.7   0    
    ## 2    20 0.516  0.0161  0.655   8.58  0    
    ## 3    30 0.492  0.0189  0.661   8.95  0    
    ## 4    40 0.496  0.0240  0.598   4.26  0    
    ## 5    50 0.445  0.0344  0.456   0.319 0.389
    ## 6    60 0.293  0.0410  0.472   4.35  0    
    ## 7    70 0.169  0.0226  0.176   0.340 0.362
    ## 8    80 0.0836 0.00849 0.0227 -7.17  1

``` r
Nwzscore50
```

    ## # A tibble: 8 x 6
    ##     age  mean     sd  NODF      z     p
    ##   <dbl> <dbl>  <dbl> <dbl>  <dbl> <dbl>
    ## 1    10  4.06  0.509  2.63 -2.80  1    
    ## 2    20  4.42  0.597  1.93 -4.17  1    
    ## 3    30  6.39  0.997  7.17  0.783 0.214
    ## 4    40  7.97  1.65   5.50 -1.50  0.934
    ## 5    50 12.5   3.62   0    -3.46  1    
    ## 6    60 23.9   9.50   3.23 -2.18  0.996
    ## 7    70 30.9  14.2    0    -2.17  0.98 
    ## 8    80 21.7   8.69   0    -2.49  1
