library(gRbase)

# Auxiliary functions

count.duplicates <- function(X){
  x <- do.call('paste', c(X, sep = '\r'))
  ox <- order(x)
  rl <- rle(x[ox])
  cbind(X[ox[cumsum(rl$lengths)],,drop=FALSE],count = rl$lengths)  
}

is_numeric0 = function(x){
  is.numeric(x) && length(x) == 0L
}

pa = function (set, object){
  
  # find the parents of nodes "set" in graph "object"
  
  amat <- as(object,"matrix")
  if (is_ugMAT(amat)) 
    return(NULL)
  pa <- names(which(amat[, set[1]] > 0))
  pa <- setdiff(pa, set)
  if (length(pa)) 
    as.numeric(pa)
  else NULL
}


library(gRbase)
library(ggm)

chain_comp = function(ess){
  
  # Given an EG with adjacency matrix "ess" returns its set of chain components Tau
  
  ess = as(ess,"igraph")
  
  amat = as.matrix(get.adjacency(ess))  # if the argument is ess.obs (the essential graph, a graph object)
  wmat = matrix(as.integer((amat + t(amat)) > 1), nrow = nrow(amat))
  wg = graph.adjacency(wmat, mode = "undirected")
  cc = clusters(wg)
  neworder <- sort.list(cc$membership, na.last = NA)
  a = matrix(0, nrow = length(cc$csize), ncol = length(cc$csize))
  b = cumsum(cc$csize)
  wmat = amat[neworder, neworder]
  
  for(i in 1: length(cc$csize)){
    for(j in 1: length(cc$csize)){
      if(j != i){
        a[i,j] = as.integer(sum(wmat[(max(b[i-1],0)+1):b[i],
                                     (max(b[j-1],0)+1):b[j]]) > 0)
      }
    }
  }
  
  rownames(a) = colnames(a) = as.character(1:length(b))
  
  chainorder = topOrder(a)
  vertorder = c()
  chainsize = c()
  
  for(k in 1:length(b)){
    vertorder = c(vertorder, which(cc$membership == chainorder[k]))
    chainsize = c(chainsize, cc$csize[chainorder[k]])
  }
  
  q = list(vert.order=vertorder,chain.size=chainsize)
  order = q$vert.order
  size = c(1,q$chain.size)
  
  Tau = list(0)
  for(i in 2:(length(size))){
    Tau[[i-1]] = order[sum(size[1:(i-1)]):(sum(size[2:i]))]
  }
  
  return(Tau)
  
}

# Main functions

marg_A = function(A, Y_tab, card_y, tau, pa_tau){
  
  ################################################################################################
  ## Function to compute the marginal likelihood of a complete categorical model for nodes in A ##
  ################################################################################################
  
  ###########
  ## INPUT ##
  ###########
  
  ## A      : set of nodes (variables)
  ## Y_tab  : contingency table of counts
  ## card_y : (q,1) vector with number of levels of each categorical variable
  ## tau    : chain component of A
  ## pa_tau : parents of nodes in chain component tau
  
  ############
  ## OUTPUT ##
  ############
  
  ## (log) marginal likelihood of a complete categorical model for variables in A
  
  card_A       = prod(card_y[A])
  card_A_compl = prod(card_y[setdiff(tau,A)])
  
  card_tau    = prod(card_y[tau])
  card_pa_tau = prod(card_y[pa_tau])
  
  card_AX = card_A*card_pa_tau
  
  a = 1/card_tau
  
  n_A_obs  = aggregate(Y_tab$count, by = as.list(Y_tab[c(A,pa_tau)]), FUN = sum)$x
  n_A_zero = rep(0, card_AX - length(n_A_obs))
  n_A      = c(n_A_obs, n_A_zero)
  
  n_pa = sum(n_A)
  
  if(length(pa_tau) > 0){
    
    n_pa_obs  = aggregate(Y_tab$count, by = as.list(Y_tab[pa_tau]), FUN = sum)$x
    n_pa_zero = rep(0, card_pa_tau - length(n_pa_obs))
    n_pa     = c(n_pa_obs, n_pa_zero)
    
  }
  
  a_A = card_A_compl*a
  a_tau = card_tau*a
  
  out_A = sum(lgamma(a_tau) - lgamma(n_pa + a_tau)) + sum(lgamma(n_A + a_A) - lgamma(a_A))
  
  return(out_A)
  
}
        
        
marg_ess = function(ess, Y_tab, card_y){
  
  #####################################################################
  ## Function to compute the marginal likelihood of a categorical EG ##
  #####################################################################
  
  ###########
  ## INPUT ##
  ###########
  
  ## ess    : adjacency matrix of the input EG
  ## Y_tab  : contingency table of counts
  ## card_y : (q,1) vector with number of levels of each categorical variable
  
  ############
  ## OUTPUT ##
  ############
  
  ## (log) marginal likelihood of ess
  
  Tau = chain_comp(ess)
  
  return(sum(unlist(lapply(Tau, m_tau, ess = ess, Y_tab = Y_tab, card_y))))
  
}


m_tau = function(tau, ess, Y_tab, card_y){
  
  ##################################################################################
  ## Function to compute the marginal likelihood of a decomposable UG model G_tau ##
  ##################################################################################
  
  ###########
  ## INPUT ##
  ###########
  
  ## tau    : chain component
  ## ess    : EG from which tau is taken
  ## card_y : (q,1) vector with number of levels of each categorical variable
  ## Y_tab  : contingency table of counts
  
  ############
  ## OUTPUT ##
  ############
  
  ## (log) marginal likelihood of a decomposable UG model G_tau
  
  tau_c = length(tau)
  
  pa_tau = as.numeric(pa(tau, as(ess, "graphNEL")))
  
  G_tau = subGraph(as.character(tau), as(ess, "graphNEL"))
  
  if(tau_c <= 2){
    
    marg_tau = marg_A(tau, Y_tab, card_y, tau, pa_tau)
    
  } else{
    
    C_set = lapply(mpd(G_tau)$cliques, as.numeric)
    S_set = lapply(mpd(G_tau)$separators, as.numeric)
    
    # exclude numeric0 elements in S_set
    
    S_set = S_set[-unlist(lapply(S_set, is_numeric0))]
    
    marg_sep = 0
    
    if(length(S_set) > 0){
      
      marg_sep = sum(sapply(S_set, marg_A, Y_tab = Y_tab, card_y = card_y, tau = tau, pa_tau = pa_tau))
    
    }
    
      marg_tau = sum(sapply(C_set, marg_A, Y_tab, card_y, tau, pa_tau)) - marg_sep
      
  }
  
  return(marg_tau)
  
} 
