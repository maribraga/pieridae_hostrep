#'---
#'title: "Pieridae host repertoire - parameters"
#'author: "Mariana Braga"
#'date: "`r format(Sys.time(), '%d %B, %Y')`"
#'output: github_document
#'---

#'-------------
#'
#' Script 1 for analyses performed in Braga et al. 2021
#' *Phylogenetic reconstruction of ancestral ecological networks through time for pierid butterflies and their host plants*,
#' Ecology Letters.

#+ include = FALSE
library(tidyverse)
library(patchwork)
library(MCMCpack)
library(coda)
library(kdensity)

  
#' ### Parameter Estimates
#' 
#' First we'll read in the `log files` from the analyses in RevBayes. 
#' They contain all the sampled parameter values during MCMC.
#' The first two chains used the time-calibrated host tree,
#' while the last two used the transformed host tree (all branch lengths = 1).

/*##### PARAMETER ESTIMATES #####

# _Read log files - output from MCMC ----
*/
chain1 <- read.table("./inference/out.2.real.pieridae.2s.log", header = TRUE)[,c(1,5,7:9)]
chain2 <- read.table("./inference/out.3.real.pieridae.2s.log", header = TRUE)[,c(1,5,7:9)]
chain3 <- read.table("./inference/out.3.bl1.pieridae.2s.log", header = TRUE)[,c(1,5,7:9)]
chain4 <- read.table("./inference/out.2.bl1.pieridae.2s.log", header = TRUE)[,c(1,5,7:9)]
colnames(chain1) <- colnames(chain2) <- colnames(chain3) <- colnames(chain4) <- 
  c("generation","clock","beta", "lambda[02]", "lambda[20]")

postt <- filter(chain1, generation >= 20000)
postt2 <- filter(chain2, generation >= 20000)
postb <- filter(chain3, generation >= 20000)
postb2 <- filter(chain4, generation >= 20000)

#' - *Convergence*
#' 
#' We'll use the Gelman and Rubin's convergence diagnostic

/*# _Convergence ----
*/
# Test convergence when measuring anagenetic distance between hosts (with time-calibrated host tree)  
gelman.diag(mcmc.list(as.mcmc(postt), as.mcmc(postt2)))

# Test convergence when measuring cladogenetic distance between hosts (with transformed host tree)  
gelman.diag(mcmc.list(as.mcmc(postb), as.mcmc(postb2)))

#' We are good to go.
#'
#' - *Mean estimates*
#' 
posterior <- bind_rows(postt,postb) %>% 
  mutate(tree = c(rep("time", 3601),rep("bl1", 3601)))

/*# _Mean estimates ----
*/
means <- group_by(posterior, tree) %>% 
  dplyr::select(-generation) %>% 
  summarise_all(mean)
means

#' - *Density plots*

/*# _Density ----
*/

#+ eval = FALSE    
# slow chunk! Read RDS file below
# calculate prior distributions  
nrow = 1e6

priors <- tibble(
  generation = 1:nrow,
  clock = rexp(nrow, rate = 10),
  beta = rexp(nrow, rate = 1),
  lambda = rdirichlet(nrow,c(1,1))[,1],
  tree = rep("prior", nrow)
)

kd_beta <- kdensity(x = priors$beta[1:5000], kernel='gamma', support=c(0,Inf), bw = 0.05)
kd_clock <- kdensity(x = priors$clock[1:5000], kernel='gamma', support=c(0,Inf), bw = 0.01)
kd_lambda <- kdensity(x = priors$lambda, kernel='gamma', support=c(0,Inf), bw = 0.01)

x10 = seq(0,10,0.002)
x1 = seq(0,1,0.0002)


plot_priors <- tibble(
  param = c(rep("beta", length(x1)), rep("clock", length(x1)), rep("lambda", length(x1))),
  x = c(x10, x1, x1),
  y = c(kd_beta(x10), kd_clock(x1), kd_lambda(x1))
) 

#+ 
# fast solution: read prior distributions from file
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

#+ densities, fig.width=18, fig.height=4, warning = FALSE
ggclock + ggbeta + gg02 + gg20 + plot_layout(ncol = 4, guides = 'collect')


#' **Bayes factor**
#' 
#' With Bayes factors we check which model has more support, 
#' the full model, where beta > 0, or the simplified model where beta = 0.
#' 
/*# _Bayes factor ----
*/

#+ warning = FALSE
d_prior <- dexp(x=0, rate=1)

kd_beta_time <- kdensity(x = filter(posterior, tree == "time") %>% pull(beta), 
                         kernel='gamma', 
                         support=c(0,Inf), 
                         bw = 0.02)
kd_beta_bl1 <- kdensity(x = filter(posterior, tree == "bl1") %>% pull(beta), 
                         kernel='gamma', 
                         support=c(0,Inf), 
                         bw = 0.02)
max_time = kd_beta_time(0)
max_bl1 = kd_beta_bl1(0)

(BF_time <- d_prior/max_time)
(BF_bl1 <- d_prior/max_bl1)

#' According to the calculated Bayes factors,
#' both analyses support the full model, where beta > 0.
#' 
#' Now we are ready to move on to Character history inference.
#'  