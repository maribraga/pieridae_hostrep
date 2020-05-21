source("functions_ancestral_states.R")

library(dplyr)
library(ape)

col_classes <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))
dat_full = read.table( file="./inference/history.bl1.txt", sep="\t", header=T, colClasses=col_classes)
tree <- read.tree("./data/bphy_pie_ladder.phy")
host_tree <- read.tree("./data/angio_pie_50tips_ladder.phy")
# make sure to get history and modules for the same analysis (bl1 or time)
all_mod_edited <- read.csv("./networks/all_mod_bl1.csv", header = T, stringsAsFactors = F)

# thin matrix to 10% of original size (for speed)
f_thin = 0.5
it = unique(dat_full$iteration)
it_thin = it[ seq(1,length(it),length.out=f_thin*length(it)) ]
dat = dat_full[dat_full$iteration %in% it_thin, ]

# record ages
ages = rev(sort(unique(all_mod_edited$age)))
n_ages = length(ages)

# collect posterior samples of graphs for each age
# each element in graphs is an array [1801 iterations (if thin .5), 131 butterflies, 50 hosts]
graphs = list()
for (i in 1:n_ages) {
    graphs[[i]] = make_matrix_samples_at_age(dat=dat, age=ages[i], s_hit=c(2), tree, host_tree)
}

# get all subgraphs (and probs) for module
all_mod_prob = compute_all_module_probs(graphs, all_mod_edited)

# save output
# saveRDS(all_mod_prob, file="all_mod_prob_bl1_thin50.rds")

# load output
# all_mod_prob = readRDS(file="all_mod_prob.rds")


# The output is a 3-dimensional list; the elements are age x module x pattern;
# You could always print the entire object.
# print(all_mod_prob)

# ... or print all modules x patterns for age "50"
# print(all_mod_prob[[ "50" ]])  #!!!! check M3 here

# ... or print all patterns for a module "M2" at age "50"
# print(all_mod_prob[[ "50" ]][[ "M2" ]])

# Printing a specific pattern for a given module and age is trickier, since we
# distinguish sampled graphs for modules using a vector-string representation of
# the subgraph, so it can be used as a key in a list. Alternatively, you can
# access list elements using an integer as an index (per usual). Module
# subgraph patterns are not sorted in order of probability, but we could
# do that if needed.
# 
# ... in any case, to print the 1st pattern for module "M2" at age "50"
# print(all_mod_prob[[ "50" ]][[ "M2" ]][[ 1 ]])



## Checking M3 at 50Ma - it's not the thinning

## Thin = 0.1
# print(all_mod_prob[[ "50" ]][[ "M3" ]])
# $`0`
# $`0`$count
# [1] 312
# 
# $`0`$graph
#           Brassicaceae
# Index_123            0
# 
# $`0`$str
# [1] "0"
# 
# $`0`$prob
# [1] 0.8642659
# 
# 
# $`1`
# $`1`$count
# [1] 49
# 
# $`1`$graph
#           Brassicaceae
# Index_123            1
# 
# $`1`$str
# [1] "1"
# 
# $`1`$prob
# [1] 0.1357341


## Thin = 0.5

# > print(all_mod_prob[[ "50" ]][[ "M3" ]])
# $`0`
# $`0`$count
# [1] 1572
# 
# $`0`$graph
#           Brassicaceae
# Index_123            0
# 
# $`0`$str
# [1] "0"
# 
# $`0`$prob
# [1] 0.8728484
# 
# 
# $`1`
# $`1`$count
# [1] 229
# 
# $`1`$graph
#           Brassicaceae
# Index_123            1
# 
# $`1`$str
# [1] "1"
# 
# $`1`$prob
# [1] 0.1271516


