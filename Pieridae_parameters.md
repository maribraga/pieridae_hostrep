Pieridae host repertoire
================
Mariana Braga
4/22/2020

-----

Script 1 for empirical study performed in Braga et al. 2020 *Evolution
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

ext_net <- ext_net_50h[which(rowSums(ext_net_50h) != 0),]
```

### Outputs from RevBayes

**Parameter estimates**

Now we’ll read in the `log files` from the analyses in RevBayes. They
contain all the sampled parameter values during MCMC.

``` r
##### PARAMETER ESTIMATES #####

chain1 <- read.table("./inference/out.2.real.pieridae.2s.log", header = TRUE)[,c(1,5,7:9)]
chain2 <- read.table("./inference/out.3.real.pieridae.2s.log", header = TRUE)[,c(1,5,7:9)]
chain3 <- read.table("./inference/out.3.bl1.pieridae.2s.log", header = TRUE)[,c(1,5,7:9)]
chain4 <- read.table("./inference/out.2.bl1.pieridae.2s.log", header = TRUE)[,c(1,5,7:9)]
colnames(chain1) <- colnames(chain2) <- colnames(chain3) <- colnames(chain4) <- 
  c("generation","clock","beta", "lambda[02]", "lambda[20]")

postt <- filter(chain1, generation >= 20000 & generation <= 200000)
postt2 <- filter(chain2, generation >= 20000 & generation <= 200000)
postb <- filter(chain3, generation >= 20000 & generation <= 200000)
postb2 <- filter(chain4, generation >= 20000 & generation <= 200000)
```

  - *Convergence*

We’ll use the Gelman and Rubin’s convergence diagnostic

``` r
gelman.diag(mcmc.list(as.mcmc(postt), as.mcmc(postt2)))
```

    ## Potential scale reduction factors:
    ## 
    ##            Point est. Upper C.I.
    ## generation        NaN        NaN
    ## clock               1       1.01
    ## beta                1       1.02
    ## lambda[02]          1       1.00
    ## lambda[20]          1       1.00
    ## 
    ## Multivariate psrf
    ## 
    ## 1.01

``` r
gelman.diag(mcmc.list(as.mcmc(postb), as.mcmc(postb2)))
```

    ## Potential scale reduction factors:
    ## 
    ##            Point est. Upper C.I.
    ## generation        NaN        NaN
    ## clock            1.01       1.03
    ## beta             1.00       1.00
    ## lambda[02]       1.00       1.00
    ## lambda[20]       1.00       1.00
    ## 
    ## Multivariate psrf
    ## 
    ## 1.01

We are good to go.

  - *Mean estimates*

<!-- end list -->

``` r
posterior <- bind_rows(postt,postb) %>% 
  mutate(tree = c(rep("time", 3601),rep("bl1", 3601)))

means <- group_by(posterior, tree) %>% 
  select(-generation) %>% 
  summarise_all(mean)
means
```

    ## # A tibble: 2 x 5
    ##   tree   clock  beta `lambda[02]` `lambda[20]`
    ##   <chr>  <dbl> <dbl>        <dbl>        <dbl>
    ## 1 bl1   0.0186  1.48       0.0268        0.973
    ## 2 time  0.0200  2.10       0.0347        0.965

  - *Density plots*

<!-- end list -->

``` r
# prior distributions  
plot_priors <- readRDS("./inference/priors.rds")

all_post <- mutate(posterior, beta = case_when(is.na(beta) ~ 0, TRUE ~ beta),
                   tree = factor(tree, levels = c("time","bl1")))

pal <- c("#08519c","#6baed6")

ggbeta <- ggplot(all_post, aes(beta)) +
  geom_line(aes(x,y), data = filter(plot_priors, param == "beta")) +
  geom_density(aes(group = tree, fill = tree, col = tree), alpha = 0.7) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  xlim(c(0,7)) +
  theme_bw()

ggclock <- ggplot(all_post, aes(clock)) +
  geom_line(aes(x,y), data = filter(plot_priors, param == "clock")) +
  geom_density(aes(group = tree, fill = tree, col = tree), alpha = 0.7) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  xlim(c(0,0.1)) +
  theme_bw()

gg02 <- ggplot(all_post, aes(`lambda[02]`)) +
  geom_line(aes(x,y), data = filter(plot_priors, param == "lambda")) +
  geom_density(aes(group = tree, fill = tree, col = tree), alpha = 0.7) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  xlim(c(0,0.1)) +
  theme_bw()

gg20 <- ggplot(all_post, aes(`lambda[20]`)) +
  geom_line(aes(x,y), data = filter(plot_priors, param == "lambda")) +
  geom_density(aes(group = tree, fill = tree, col = tree), alpha = 0.7) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  xlim(c(0.9,1)) +
  theme_bw()
```

``` r
ggclock + ggbeta + gg02 + gg20 + plot_layout(ncol = 4, guides = 'collect')
```

    ## Warning: Removed 4500 row(s) containing missing values (geom_path).

    ## Warning: Removed 1500 row(s) containing missing values (geom_path).

    ## Warning: Removed 4500 row(s) containing missing values (geom_path).
    
    ## Warning: Removed 4500 row(s) containing missing values (geom_path).

![](Pieridae_parameters_files/figure-gfm/densities-1.png)<!-- -->

**Bayes factor**

``` r
d_prior <- dexp(x=0, rate=1)

kd_beta_time <- kdensity(x = filter(posterior, tree == "time") %>% pull(beta), 
                         kernel='gamma', 
                         support=c(0,Inf), 
                         bw = 0.02)
```

    ## Warning: namespace 'extraDistr' is not available and has been replaced
    ## by .GlobalEnv when processing object ''

``` r
kd_beta_bl1 <- kdensity(x = filter(posterior, tree == "bl1") %>% pull(beta), 
                         kernel='gamma', 
                         support=c(0,Inf), 
                         bw = 0.02)
max_time = kd_beta_time(0)
max_bl1 = kd_beta_bl1(0)

(BF_time <- d_prior/max_time)
```

    ## [1] 7.029484

``` r
(BF_bl1 <- d_prior/max_bl1)
```

    ## [1] 65545525
