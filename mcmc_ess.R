# MCMC scheme

library(pcalg)
library(igraph)
library(gRbase)
library(graph)
library(ggm)

source("move.r")
source("marg_ess.r")

count.duplicates <- function(X){
  x <- do.call('paste', c(X, sep = '\r'))
  ox <- order(x)
  rl <- rle(x[ox])
  cbind(X[ox[cumsum(rl$lengths)],,drop=FALSE],count = rl$lengths)  
}

n.edge <- function(A){
  
  # A: adjacency matrix of a graph
  length(which(A[lower.tri(A)] == 1 | t(A)[lower.tri(A)] == 1))
  
}

post_ess_graphs = function(Y, m, T){
  
  ##############################################################
  ## Function to move from one EG to another in the EG space) ##
  ##############################################################
  
  ###########
  ## INPUT ##
  ###########
  
  ## Y : (n,q) data matrix
  ## m : maximum number of edges in the EG space (sparsity constraint)
  ## T : number of MCMC iterations
  
  ############
  ## OUTPUT ##
  ############
  
  ## chain.G : a collection of EGs visited by the MCMC chain
  
  ## represented as a (q*q, T) matrix with vectorized adjacency matrices of the EGs (one for each column)
  
  
  q = ncol(Y)
  
  isFALSE = function(x){identical(x, FALSE)}
  
  card_y = sapply(1:q, function(i) length(unique(Y[,i])))
  
  Y_tab = count.duplicates(Y)
  
  # initial graph (empty graph)
  
  A_0 = matrix(0,q,q)
  
  colnames(A_0) = rownames(A_0) = 1:q
  
  # matrix collecting the output of the MCMC (each column is the by-column-vectorized adjacency matrix of the accepted graph)
  
  A_chain = matrix(0, nrow = q*q, ncol = T)
  A_chain[,1] = A_0
  
  A_old = A_0
    
  colnames(A_old) = rownames(A_old) = 1:q
    
  m_old = marg_ess(A_old, Y_tab, card_y)
  
  colnames(A_0) = rownames(A_0) = 1:q
  
  # MCMC iterations
  
  p = 1.5/(2*q-2)
  
  for(t in 2:T){
    
    # proposed move
    
    prop = move(A_old, m, q)
    
    A_new = prop[[1]]
    m_new = marg_ess(ess = A_new, Y_tab, card_y = card_y)
    
    type.operator = prop$O_new
    
    # Distinguish 3 cases:
    
    if(type.operator == 1 | type.operator == 3){
      
      # (1) Insert an edge
      
      logprior = log(p/(1-p))
      
    }else{
      
      if(type.operator == 2 | type.operator == 4){
        
        # (2) Delete an edge
        
        logprior = log((1-p)/p)
        
      }else{
        
        # (3) Reverse a directed edge, Make a v-structure or Remove a v-structure
        
        logprior = log(1)
        
      }
      
    }
    
    
    # acceptance ratio
    ratio = min(0, m_new - m_old + logprior)
    
    if(log(runif(1,0,1)) < ratio){ # accept move
      
      A_old = A_new
      m_old = m_new
      
    }
    
    # store chain value
    
    A_chain[,t] = A_old
    
    if(t%%10 == 0) {print(paste0("Iteration ", t))}
    
  }
  
  return(list(chain.G = A_chain))
  
}
