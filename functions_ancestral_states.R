require(stringr)
require(data.table)
#require(igraph)

# find lineages (and their repertoires) that exist during the specified age (time slice)
make_dat_timeslice = function(dat, age, tree, host_tree) {
    
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

# make matrix with marginal posterior probabilities of interactions at specified ages (time slices)
make_matrix_at_age = function(dat, age, s_hit=c(2), tree, host_tree, drop_empty=T) {
    
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
        dat_it = make_dat_timeslice(dat_it, age, tree, host_tree)
    
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

# make matrix with posterior probabilities of repertoires at specified ages (time slices)
make_matrix_samples_at_age = function(dat, age, s_hit=c(2), tree, host_tree, drop_empty=T) {
    
    iterations = sort(unique(dat$iteration))
    n_iter = length(iterations)
  
    # get dimensions
    n_host_tip = length( str_split( dat$start_state[1], "" )[[1]] )
    n_parasite_lineage = length(unique(dat$node_index))
    n_parasite_tip = (n_parasite_lineage + 1) / 2
    
    m_names = list( 1:n_iter,
                    c(rev(tree$tip.label), paste0("Index_",(n_parasite_tip+1):n_parasite_lineage)),
                    host_tree$tip.label )

    m = array(0, dim=c(n_iter, n_parasite_lineage, n_host_tip), dimnames=m_names)
    #m = matrix(data = 0, nrow = n_parasite_lineage, ncol = n_host_tip)
    
    for (it_idx in 1:length(iterations)) {
        it = iterations[it_idx]
        
        # get dataset for iteration
        dat_it = dat[ dat$iteration==it, ]
        
        # extract relevant branches
        dat_it = make_dat_timeslice(dat_it, age, tree, host_tree)
    
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


compute_all_module_probs = function(graphs, modules) {
    
    # initialize the module-age-prob list
    mod_prob_list = list()
    
    # initialize and re-order age vector
    ages = as.character(rev(sort(unique(modules$age))))
    n_ages = length(ages)
    
    # get the full set of row (insect) and column (plant) names
    g_row_names = rownames(graphs[[1]][1,,])
    g_col_names = colnames(graphs[[1]][1,,])
    
    # compute for each age
    for (i in 1:n_ages) {
        
        # initialize the module probs for this age
        mod_prob_list[[ ages[i] ]] = list()

        # get the list of modules for each age
        modules_for_age = modules[ modules$age == ages[i], ]
        module_names_for_age = sort(unique( modules_for_age$module ))

        # get the number of iterations
        n_it = dim(graphs[[i]])[1]

        # process each each module m at age i
        for (m_idx in 1:length(module_names_for_age))
        {
            # get the module name
            m_name = module_names_for_age[m_idx]
            
            # extract the module of interest from the module dataframe
            m = modules_for_age[ modules_for_age$module == m_name, ]

            # initialize info for a particular module+age
            mod_prob_list[[ ages[i] ]][[ m_name ]] = list()

            # find num. matches across iterations for the module at age i
            for (it in 1:n_it)
            {
                
                # get the reduced set of names for m (the module-age pair)
                m_row_names = intersect( m$name, g_row_names )
                m_col_names = intersect( m$name, g_col_names )
                
                # get posterior sample of graph for iteration it at age i
                g_it = graphs[[i]][it,,]
                
                # now, extract the subgraph for module m for iteration it at age i
                m_it = matrix(g_it[ m_row_names, m_col_names],
                              nrow=length(m_row_names),
                              ncol=length(m_col_names))
                colnames(m_it)=m_col_names
                rownames(m_it)=m_row_names
                
                # because the row and col names are constant, each m_it will have
                # the same dimensions; this is useful, because it means we can
                # then convert m_it into a string representation, and use that
                # string to identify unique subgraph structures among m_it
                
                # so, first, let's create our string representation of m_it
                m_it_flat = as.character(c(m_it))
                s_it_flat = paste( m_it_flat, collapse="" )
                
                # then we create a new record for the subgraph for
                # m_it if no such previous subgraph has been previously recorded
                subgraph_exists = is.null( mod_prob_list[[ ages[i] ]][[ m_name ]][[ s_it_flat ]] )
                
                # ... create new record if needed!
                if (subgraph_exists) {
                    mod_prob_list[[ ages[i] ]][[ m_name ]][[ s_it_flat ]] = list()
                    mod_prob_list[[ ages[i] ]][[ m_name ]][[ s_it_flat ]]$count = 0
                    mod_prob_list[[ ages[i] ]][[ m_name ]][[ s_it_flat ]]$graph = m_it
                    mod_prob_list[[ ages[i] ]][[ m_name ]][[ s_it_flat ]]$str = s_it_flat
                }
                
                # increment the count for the subgraph pattern
                m_it_count = mod_prob_list[[ ages[i] ]][[ m_name ]][[ s_it_flat ]]$count
                mod_prob_list[[ ages[i] ]][[ m_name ]][[ s_it_flat ]]$count = m_it_count + 1
            }
            
            
            # revisit each module
            for (k in 1:length(mod_prob_list[[ ages[i] ]][[ m_name ]])) {
                
                # now compute the probability of the module subgraph pattern
                m_it_prob = mod_prob_list[[ ages[i] ]][[ m_name ]][[ k ]]$count / n_it
                mod_prob_list[[ ages[i] ]][[ m_name ]][[ k ]]$prob = m_it_prob
            }
            
        }
    }

    return(mod_prob_list)
}


compute_all_module_probs_more = function(graphs, modules, tol=0, tol_row=0, tol_col=0) {
    
    # initialize the module-age-prob list
    mod_prob_list = list()
    
    # initialize and re-order age vector
    ages = as.character(rev(sort(unique(modules$age))))
    n_ages = length(ages)
    
    # get the full set of row (insect) and column (plant) names
    g_row_names = rownames(graphs[[1]][1,,])
    g_col_names = colnames(graphs[[1]][1,,])
    
    # compute for each age
    for (i in 1:n_ages) {
        
        # initialize the module probs for this age
        mod_prob_list[[ ages[i] ]] = list()

        # get the list of modules for each age
        modules_for_age = modules[ modules$age == ages[i], ]
        module_names_for_age = sort(unique( modules_for_age$module ))

        # get the number of iterations
        n_it = dim(graphs[[i]])[1]

        # process each each module m at age i
        for (m_idx in 1:length(module_names_for_age))
        {
            # get the module name
            m_name = module_names_for_age[m_idx]
            
            # extract the module of interest from the module dataframe
            m = modules_for_age[ modules_for_age$module == m_name, ]

            # initialize info for a particular module+age
            mod_prob_list[[ ages[i] ]][[ m_name ]] = list()

            # find num. matches across iterations for the module at age i
            for (it in 1:n_it)
            {
                
                # get the reduced set of names for m (the module-age pair)
                m_row_names = intersect( m$name, g_row_names )
                m_col_names = intersect( m$name, g_col_names )
                
                # get posterior sample of graph for iteration it at age i
                g_it = graphs[[i]][it,,]
                
                # now, extract the subgraph for module m for iteration it at age i
                m_it = matrix(g_it[ m_row_names, m_col_names],
                              nrow=length(m_row_names),
                              ncol=length(m_col_names))
                rownames(m_it)=m_row_names
                colnames(m_it)=m_col_names
                
                # make module w/ for module row-set with expanded cols
                m_it_row = matrix(g_it[ m_row_names, ],
                              nrow=length(m_row_names),
                              ncol=length(g_col_names))
                rownames(m_it_row)=m_row_names
                colnames(m_it_row)=g_col_names
                m_it_row[ m_row_names, m_col_names ] = m_it_row[ m_row_names, m_col_names ] - m_it[ m_row_names, m_col_names ]
                n_row_excess = sum(m_it_row)
                
                
                # make module w/ for module col-set with expanded rows
                m_it_col = matrix(g_it[ m_row_names, ],
                              nrow=length(g_row_names),
                              ncol=length(m_col_names))
                rownames(m_it_col)=g_row_names
                colnames(m_it_col)=m_col_names
                m_it_col[ m_row_names, m_col_names ] = m_it_col[ m_row_names, m_col_names ] - m_it[ m_row_names, m_col_names ]
                n_col_excess = sum(m_it_col)
                
                # what is total excess among row/col relationships?
                n_excess = n_row_excess + n_col_excess
                
                # because the row and col names are constant, each m_it will have
                # the same dimensions; this is useful, because it means we can
                # then convert m_it into a string representation, and use that
                # string to identify unique subgraph structures among m_it
                
                # so, first, let's create our string representation of m_it
                m_it_flat = as.character(c(m_it))
                s_it_flat = paste( m_it_flat, collapse="" )
                
                # then we create a new record for the subgraph for
                # m_it if no such previous subgraph has been previously recorded
                subgraph_exists = is.null( mod_prob_list[[ ages[i] ]][[ m_name ]][[ s_it_flat ]] )
                
                # ... create new record if needed!
                if (subgraph_exists) {
                    mod_prob_list[[ ages[i] ]][[ m_name ]][[ s_it_flat ]] = list()
                    mod_prob_list[[ ages[i] ]][[ m_name ]][[ s_it_flat ]]$count = 0
                    mod_prob_list[[ ages[i] ]][[ m_name ]][[ s_it_flat ]]$graph = m_it
                    mod_prob_list[[ ages[i] ]][[ m_name ]][[ s_it_flat ]]$str = s_it_flat
                }
                
                # increment the count for the subgraph pattern
                m_it_count = mod_prob_list[[ ages[i] ]][[ m_name ]][[ s_it_flat ]]$count
                mod_prob_list[[ ages[i] ]][[ m_name ]][[ s_it_flat ]]$count = m_it_count + 1
            }
            
            
            # revisit each module
            for (k in 1:length(mod_prob_list[[ ages[i] ]][[ m_name ]])) {
                
                # now compute the probability of the module subgraph pattern
                m_it_prob = mod_prob_list[[ ages[i] ]][[ m_name ]][[ k ]]$count / n_it
                mod_prob_list[[ ages[i] ]][[ m_name ]][[ k ]]$prob = m_it_prob
            }
            
        }
    }

    return(mod_prob_list)
}



# make matrix with marginal posterior probabilities of interactions at internal nodes
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



### GRAVEYARD OF DEADLY CODE


# compute_module_probs = function(graphs, modules) {
#     
#     
#     # age x mod x (graph, prob)
#     mod_prob_list = list()
#     
#     modules$prob = 0
#     ages = rev(sort(unique(modules$age)))
#     n_ages = length(ages)
#     
#     g_row_names = rownames(graphs[[1]][1,,])
#     g_col_names = colnames(graphs[[1]][1,,])
#     
#     # compute for each age
#     for (i in 1:1) { #n_ages) {
#     
#         mod_prob_list[[ ages[i] ]] = list()
#         
#         modules_for_age = modules[ modules$age == ages[i], ]
#         module_names_for_age = sort(unique( modules_for_age$module ))
# 
#         n_it = dim(graphs[[i]])[1]
#         
#         # compute for each module at age i
#         for (m_idx in 1:length(module_names_for_age))
#         {
#             m_name = module_names_for_age[m_idx]
#             mod_prob_list[[ ages[i] ]][[ m_name ]] = list()
#             
#             # create temp. variable for the module
#             #print( module_names_for_age[m_idx])
#             m = modules_for_age[ modules_for_age$module == m_name, ]
#             
#             
#             # get row/col idx names for module
#             m_row_names = intersect( m$name, g_row_names )
#             m_col_names = intersect( m$name, g_col_names )
#             
#             n_match = 0
#             
#             # find num. matches across iterations for the module at age i
#             for (it in 1:n_it)
#             {
#                 # get graph sample with index it at age i
#                 g_it = graphs[[i]][it,,]
#                 #print(g_it)
#                 print(m_row_names)
#                 print(m_col_names)
#                 print(g_it[ m_row_names, m_col_names] )
#                 #break
#                 # all nodes in m must be fully connected to match
#                 match = TRUE
#                 
#                 for (j in 1:length(m_row_names))
#                 {
#                     s_j = m_row_names[j]
#                     for (k in 1:length(m_col_names))
#                     {
#                         s_k = m_col_names[k]
#                         #cat(it, ages[i], s_j, s_k,  g_it[ c(s_j), c(s_k) ], "\n")
#                         if ( g_it[ c(s_j), c(s_k) ] == 0 ) {
#                             #print(g_it[ c(s_j), c(s_k) ])
#                             match = FALSE
#                         }
#                         if (!match) break
#                     }
#                     if (!match) break
#                 }
#                 if (match)
#                 {
#                     n_match = n_match + 1
#                 }
#             }
#             m$prob = n_match/n_it
#             #print(m)
#             modules_for_age[ modules_for_age$module == module_names_for_age[m_idx], ] = m
#         }
#         modules[ modules$age == ages[i], ] = modules_for_age
#     }
#     
#     return(modules)
# }
# 


# compute_prob_subgraph_pattern_old = function(graphs, pattern, tol=0) {
#     n_it = dim(graphs)[1]
#     min_score = length(graphs[1,,]) - tol
#     print(min_score)
#     print(length(pattern))
#     n_match = 0
#     for (i in 1:n_it) {
#         score = sum(c(graphs[i,,]==pattern))
#         print(score)
#         if (score >= min_score) {
#             n_match = n_match + 1
#         }
#     }
#     return(n_match/n_it)
# }


# compute_prob_subgraph_edges = function(graphs, edges, tol=0) {
#     n_it = dim(graphs)[1]
#     print(n_it)
#     n_match = 0
#     for (i in 1:n_it) {
#         #print(i)
#         match = T
#         for (j in 1:nrow(edges)) {
#             e = edges[j,]
#             #print(e)
#             if ( graphs[i,e[1],e[2]] == 0 ) {
#                 match = F
#                 break
#             }
#         }        
#         if (match) {
#             cat("it = ", i, " n_match = ", n_match, "\n")
#             n_match = n_match + 1
#         }
#     }
#     return(n_match/n_it)
# }


get_match_prob = function(x, p=0.9) {
  frac = p
  ret = x
  for (i in 1:length(x)) {
    for (j in 1:length(x[[i]])) {
      # get subgraph w/ highest prob
      prob = 0
      idx = 0
      for (k in 1:length(x[[i]][[j]])) {
        g = x[[i]][[j]][[k]]
        if (g$prob > prob) {
          prob = g$prob
          idx = k
        }
      }
      gmax = x[[i]][[j]][[idx]]
      #print(gmax)
      n_int = length(gmax$graph)
      ret[[i]][[j]] = gmax
      ret[[i]][[j]]$prob = 0
      ret[[i]][[j]]$graph = list()
      kk = 1
      for (k in 1:length(x[[i]][[j]])) {
        gk = x[[i]][[j]][[k]]
        diff = gmax$graph - gk$graph
        count_small = sum(diff[ diff < 0 ])
        count_large = sum(diff[ diff > 0 ])
        gk_frac  = 1.0 - ((count_small + count_large) / n_int)
        if (gk_frac >= frac) {
          ret[[i]][[j]]$prob = ret[[i]][[j]]$prob + gk$prob
          ret[[i]][[j]]$graph[[kk]] = gk$graph
          kk = kk + 1
          #print(gk$graph)
        }
      }
    }
  }
  return(ret)
}
 