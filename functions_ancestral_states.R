require(stringr)
require(data.table)
#require(igraph)

# find lineages (and their repertoires) that exist during the specified age (time slice)
make_dat_timeslice = function(dat, age) {
    
    # reduce dat to only include relevant branches
    dat2 = dat[dat$branch_start_time>=age & dat$branch_end_time<=age, ]
    
    nodes = unique(dat2$node_index)
    
    dat3 = rbind( dat2[ dat2$transition_type=="no_change", ],
                  dat2[ dat2$transition_type=="anagenetic" & dat2$transition_time>=age, ] )
    
    ret = list()
    for (i in 1:length(nodes)) {
      if(!(nodes[i] %in% dat3$node_index)){               # if changes happened only after age
        parent = dat[ dat$node_index==nodes[i], 12][1]    # get state at parent node           ## option would be to get start.state instead of end.state
        dat4 = dat[ dat$node_index==parent, ]
        dat4$node_index = nodes[i]                        # fix node_index - back to child nodes
        if(nrow(dat4)==1){                                # when type is no_change, time is NA and which.min doesn't work
          ret[[i]] = dat4
        } else{
          ret[[i]] = dat4[ which.min(dat4$transition_time), ]
        }
      } else {
      dat4 = dat3[ dat3$node_index==nodes[i], ]
      if(nrow(dat4)==1){
        ret[[i]] = dat4
      } else{
      ret[[i]] = dat4[ which.min(dat4$transition_time), ]  # get state at minimum, not maximum time (that is greater than age)
      }
      }
    }
    return(rbindlist(ret))
}

# make matrix with posterior probabilities of interactions at specified ages (time slices)
make_matrix_at_age = function(dat, age, s_hit=c(2), drop_empty=T) {
    
  iterations = sort(unique(dat$iteration))
  n_iter = length(iterations)
  
  # get dimensions
  n_host_tip = length( str_split( dat$start_state[1], "" )[[1]] )
  n_parasite_lineage = length(unique(dat$node_index))
  n_parasite_tip = (n_parasite_lineage + 1) / 2
  m = matrix(data = 0, nrow = n_parasite_lineage, ncol = n_host_tip)
    
    for (it in iterations) {
        
        # get dataset for iteration
        dat_it = dat[ dat$iteration==it, ]
        
        # extract relevant branches
        dat_it = make_dat_timeslice(dat_it, age)
    
        # add edges ( parasite x host )
        for (i in 1:nrow(dat_it)) {
            n_idx = dat_it$node_index[i]
            s = as.numeric( str_split (dat_it$end_state[i], "")[[1]] )
            s_idx = s %in% s_hit
            #print(n_idx)
            #print(s_idx)
            m[ n_idx, s_idx ] = m[ n_idx, s_idx ] + 1
        }
    }
    
    row.names(m) <- c(rev(tree$tip.label), paste0("Index_",(n_parasite_tip+1):n_parasite_lineage))
    colnames(m) <- host_tree$tip.label
  
    # remove empty rows/columns
    if (drop_empty) {
        m = m[ rowSums(m)!=0, ]
        m = m[ ,colSums(m)!=0 ]
    }
    
    # convert to probability
    m = m * (1/n_iter)
    
    return(m)
}

# make matrix with posterior probabilities of interactions at specified ages (time slices)
make_matrix_samples_at_age = function(dat, age, s_hit=c(2), drop_empty=T) {
    
    iterations = sort(unique(dat$iteration))
    n_iter = length(iterations)
  
    # get dimensions
    n_host_tip = length( str_split( dat$start_state[1], "" )[[1]] )
    n_parasite_lineage = length(unique(dat$node_index))
    n_parasite_tip = (n_parasite_lineage + 1) / 2
    
    m_names = list( 1:n_iter,
                    c(rev(tree$tip.label), paste0("Index_",(n_parasite_tip+1):n_parasite_lineage)),
                    host_tree$tip.label )

    m = array(0, dim=c(n_iter, n_host_tip, n_parasite_tip), dimnames=m_names)
    #m = matrix(data = 0, nrow = n_parasite_lineage, ncol = n_host_tip)
    
    for (it_idx in 1:length(iterations)) {
        it = iterations[it_idx]
        
        # get dataset for iteration
        dat_it = dat[ dat$iteration==it, ]
        
        # extract relevant branches
        dat_it = make_dat_timeslice(dat_it, age)
    
        # add edges ( parasite x host )
        for (i in 1:nrow(dat_it)) {
            n_idx = dat_it$node_index[i]
            s = as.numeric( str_split (dat_it$end_state[i], "")[[1]] )
            s_idx = s %in% s_hit
            #print(n_idx)
            #print(s_idx)
            m[it_idx, n_idx, s_idx ] = 1
        }
    }
  
    # remove empty rows/columns
    #if (drop_empty) {
    #    m = m[ rowSums(m)!=0, ]
    #    m = m[ ,colSums(m)!=0 ]
    #}
    
    # convert to probability
    #m = m * (1/n_iter)
    
    return(m)
}

compute_prob_subgraph_pattern = function(graphs, pattern, tol=0) {
    n_it = dim(graphs)[1]
    min_score = length(graphs[1,,]) - tol
    n_match = 0
    for (i in 1:n_it) {
        score = sum(c(graphs[i,,]==pattern))
        if (score >= min_score) {
            n_match = n_match + 1
        }
    }
    return(n_match/n_it)
}

# make matrix with posterior probabilities of interactions at internal nodes
make_matrix_nodes = function(dat, nodes, state) { 
  
  dat <- filter(dat, node_index %in% nodes)
  iterations = sort(unique(dat$iteration))
  n_iter = length(iterations)
  
  # get dimensions
  n_host_tip = length( str_split( dat$start_state[1], "" )[[1]] )
  n_parasite_lineage = length(unique(dat$node_index))
  #n_parasite_tip = (n_parasite_lineage + 1) / 2
  g = matrix(data = 0, nrow = n_parasite_lineage, ncol = n_host_tip)
  
  for (it in iterations) {
    dat_it = dat[ dat$iteration == it, ]
    ret = list()
    for (i in 1:length(nodes)) {
      dat2 = dat_it[ dat_it$node_index == nodes[i], ]
      if(nrow(dat2)==1){
        ret[[i]] = dat2
      } else{
        ret[[i]] = dat2[ which.min(dat2$transition_time), ]
      }
    }
    
    ret <- rbindlist(ret)
    
    for (r in 1:nrow(ret)) {
      #n_idx = ret$node_index[r]
      s = as.numeric( str_split (ret$end_state[r], "")[[1]] )
      s_idx = s %in% state
      #print(n_idx)
      #print(s_idx)
      g[ r, s_idx ] = g[ r, s_idx ] + 1
    }
  }
  
  # convert to probability
  g = g * (1/n_iter)
  
  return(g)
}
