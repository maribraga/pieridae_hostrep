source("functions_ancestral_states.R")

library(dplyr)
library(ape)

col_classes <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))
dat_full = read.table( file="./inference/history.time.txt", sep="\t", header=T, colClasses=col_classes)
tree <- read.tree("./data/bphy_pie_ladder.phy")
host_tree <- read.tree("./data/angio_pie_50tips_ladder.phy")
all_mod_edited <- read.csv("./networks/all_mod_bl1.csv", header = T, stringsAsFactors = F)

dat = dat_full[dat_full$iteration %in% unique(dat_full$iteration)[1:100], ]
#dat = dat_full
ages = c(80,70,60,50,40,30,20,10,0)
#ages = c(60)
n_ages = length(ages)



graphs = list()
for (i in 1:n_ages) {
    graphs[[i]] = make_matrix_samples_at_age(dat=dat, age=ages[i], s_hit=c(2), tree, host_tree)
    # modules = make_modules_for_age(dat)
}

# get all subgraphs (and probs) for module
all_mod_prob = compute_all_module_probs(graphs, all_mod_edited)

# The output is a 3-dimensional list; the elements are age x module x pattern;
# You could always print the entire object.
# print(all_mod_prob)

# ... or print all modules x patterns for age "50"
# print(all_mod_prob[[ "50" ]])

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