source("functions_ancestral_states.R")

library(dplyr)
library(ape)
library(data.table)

col_classes <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))
dat_full = read.table( file="./inference/history.bl1.txt", sep="\t", header=T, colClasses=col_classes)
tree <- read.tree("./data/bphy_pie_ladder.phy")
host_tree <- read.tree("./data/angio_pie_50tips_bl1.phy")
# make sure to get history and modules for the same analysis (bl1 or time)
all_mod_edited <- read.csv("./networks/all_mod_bl1.csv", header = T, stringsAsFactors = F)


# NOTE: this manually increments various node index values in the history file
#       may or may not be necessary (check w/ Mari) 
# - not necessary if you get history.bl1.txt from the zip file
# node_idx = c("node_index","parent_index","child1_index","child2_index")
# dat_full[,node_idx] = dat_full[,node_idx] + 1

# this file is used to compare values between the new `graphs` object and
# Mari's pre-computed values
list_m_at_ages <- readRDS("./inference/list_m_at_ages_bl1.rds")

# thin matrix to 10% of original size (for speed)
f_thin = 0.1
it = unique(dat_full$iteration)
it = it[ (it >= 2e4 & it <= 2e5) ]
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

# NOTE: uncomment these lines to compare marginal probs from graphs against
#       the marginal probs in list_m_at_ages
# g_check = t( colSums(graphs[[1]],dim=1) / length(it_thin) )
# row_idx = c("Fabaceae")
# col_idx = c("Index_70","Index_129")
# g_check[row_idx,col_idx]
# list_m_at_ages[[1]][row_idx,col_idx]

# Get all subgraphs (and probs) for module
all_mod_prob = compute_all_module_probs(graphs, all_mod_edited)

# Get the probability for all subgraph patterns for a module that share
# a percent of edges (p) that exactly match the most probable module, which is
# assumed to be the pre-computed module in `all_mod_edited'. We call this
# the 'soft-match' posterior probability.
#
# For example, suppose we wanted to compute the soft-match posterior probability
# for the module
#       Y Z
#     A 0 1
#     B 1 1
#     C 1 1
#     D 1 0
#
# We can encode this as G1 = (AZ, BY, BZ, CY, CZ, DY), i.e. a subgraph with 8
# potential interactions. Now suppose G1 has posterior probability 0.6. If we
# set p=100, then the only posterior samples with an exact match for the 6 edges
# and 2 non-edges will be included in our 'soft-match' metric. Now suppose
# G2 = (AZ, BY, BZ, CY, DY) has pp=0.3 and G3 = (AZ, BY, CY, DY) has pp=0.1.
# If we set p=75 for the soft-match then 75% of the interactions exactly match
# G1, so the soft-match-75 to G1 woudl be pp_G1 + pp_G2 = 0.9.

all_mod_prob_90 = get_match_prob(all_mod_prob, p=0.9)


# save output
# saveRDS(all_mod_prob, file="networks/all_mod_prob_bl1_thin10.rds")

# load output
# all_mod_prob = readRDS(file="networks/all_mod_prob_bl1_thin10.rds")


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






