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


# 
# ###
# # NOTE: Here is an example for a `module_list` variable
# # that we will want to generate automatically. 
# 
# module_template = graphs[[1]][1,,]
# module_template[ 1:length(module_template) ] = 0
# 
# module_example = module_template
# module_example_edges = matrix(c(
#     "Index_70", "Fabaceae",
#     "Index_70", "Brassicaceae"), ncol=2, byrow=T)
# for (i in 1:nrow(module_example_edges)) {
#     e = module_example_edges[i,]
#     module_example[ c(e[1]), c(e[2]) ] = 1
# }
# 
# module_list = list( module_example_edges )
# 
# for (i in length(module_list)) {
#     module = module_list[[i]]
#     p = compute_prob_subgraph_edges( graphs[[1]], module )
#     print(p)
# }
# 
# 
# for (k in 1:dim(graphs[[1]])[1]) {
#     g = graphs[[1]][k,,]
#     ni = rownames(g)
#     nj = colnames(g)
#     for (i in 1:dim(g)[1]) {
#         for (j in 1:dim(g)[2]) {
#             if (g[i,j] == 1) {
#                 cat(k,"--",ni[i],"--",nj[j],":",i,j,"\n")
#             }
#         }
#     }
#     cat("\n")
# }
#         
#     