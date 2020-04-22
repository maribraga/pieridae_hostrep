#'---
#'title: "Pieridae host repertoire"
#'author: "Mariana Braga"
#'date: "4/22/2020"
#'output: github_document
#'---

#' Script for empirical study performed in Braga et al. 2020
#' Evolution of butterfly-plant networks revealed by Bayesian inference of host repertoire

library(ape)
library(phytools) # needed?
library(dispRity)
library(ggtree)
library(picante)

library(tidyverse)
library(patchwork)
library(wesanderson)

library(MCMCpack)
library(coda)
library(kdensity)
#library(phylotate) # needed?
#library(data.table) # needed in script sourced

library(bipartite)
library(ggraph)
library(tidygraph)
library(igraph)
#library(ndtv)        # needed?
#library(intergraph)  # needed?


#### INPUT DATA ####

# _Read trees ----
path_data <- "./data/"
tree <- read.tree(paste0(path_data,"bphy_pie_ladder.phy"))
host_tree <- read.tree(paste0(path_data,"angio_pie_50tips_ladder.phy"))

## for figure 2
host_tree_bl1 <- read.tree(paste0(path_data,"angio_pie_50tips_bl1.phy"))

ggtree(host_tree_bl1) + 
  theme_tree2() + 
  geom_tiplab(align=TRUE, linesize=.5) + 
  xlim(0, 25)

