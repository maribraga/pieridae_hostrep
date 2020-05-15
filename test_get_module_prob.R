source("functions_ancestral_states.R")

library(dplyr)
library(ape)

col_classes <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))
dat = read.table( file="./inference/history.time.txt", sep="\t", header=T, colClasses=col_classes)
tree <- read.tree("./data/bphy_pie_ladder.phy")
host_tree <- read.tree("./data/angio_pie_50tips_ladder.phy")

# ages = c(80,70,60,50,40,30,20,10,0)
ages = c(60)
n_ages = length(ages)

graphs = list()
for (i in 1:n_ages) {
    graphs[[i]] = make_matrix_samples_at_age(dat=dat, age=ages[i], s_hit=c(2), tree, host_tree)
    # modules = make_modules_for_age(dat)
}
    
    