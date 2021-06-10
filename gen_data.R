library(pcalg)
library(gRbase)
library(mvtnorm)

gen_dag_data = function(i, n, q, p){
  
  #########################################################
  ## Function to simulate DAG-dependent categorical data ##
  #########################################################
  
  ###########
  ## INPUT ##
  ###########
  
  ## i : seed for random generation
  ## n : number of observations
  ## q : number of variables
  ## p : probability of edge inclusion
  
  ############
  ## OUTPUT ##
  ############
  
  ## D : the true DAG 
  ## Y : a (n,q) data matrix
  
  set.seed(i)
  
  true_dag = randomDAG(q, prob = p)
  
  A = as(true_dag, "matrix")
  
  B = A*matrix(runif(q*q, 0.1, 1), q, q)*sample(c(-1, 1), size = q*q, replace = TRUE)
  
  B_1 = t(B)
  
  I = diag(rep(0.1, q))
  
  Sigma_cond = diag(rep(1, q))
  
  S = solve(I - B_1)%*%Sigma_cond%*%t(solve(I - B_1))
  
  mu = c(rep(0, q))
  
  library(mvtnorm)
  
  Y = data.frame(rmvnorm(n, mu, S))
  
  # generate random cut-offs uniformly in [-1, 1]
  
  cut_offs = rep(0.1, q)
  
  for(j in 1:q){
    
    gamma_j = runif(1, quantile(Y[,j], 0.05), quantile(Y[,j], 0.95))
    
    Y[,j][Y[,j] >= gamma_j] = 1; Y[,j][Y[,j] < gamma_j] = 0
    
  }
  
  return(list(D = true_dag, Y = Y))
  
}
