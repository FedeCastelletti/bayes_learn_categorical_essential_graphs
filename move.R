# Auxiliary functions

source("operators.r")

n.edge = function(A){
  length(which(A[lower.tri(A)] == 1 | t(A)[lower.tri(A)] == 1))
}

names   = c("action","test","x","y","z")
types   = paste0("perfect.",c("iu","du","id","dd","mv","rv","rd"))
actions = c("iu","du","id","dd","mv","rv","rd")

move <- function(A, m, q = q){
  
  ##############################################################
  ## Function to move from one EG to another in the EG space) ##
  ##############################################################
  
  ###########
  ## INPUT ##
  ###########
  
  ## A : adjacency matrix of the input EG
  ## m : maximum number of edges in the graph space (sparsity constraint)
  ## q : number of nodes in the graph space
  
  ############
  ## OUTPUT ##
  ############
  
  ## A_new : a new EG (direct successor of input EG)
  ## O_new : the type of operator applied to obtain the new EG
  ## nodes : nodes involved in the move
  
  # Seven types of operators are:
  
  # iu : insert an undirected edge
  # du : delete an undirected edge
  # id : insert a directed edge
  # dd : delete a directed edge
  # mv : make a v-structure
  # rv : remove a v-structure
  # rd : reverse a directed edge
  
  G_und = G_undirected(A)
  
  A_dir = A
  A_dir[which(A == 1 & t(A) == 1,TRUE)] <- 0
  A_und = A + t(A)
  A_und[upper.tri(A_und,diag=TRUE)] = NA
  
  iu_set = c()
  du_set = c()
  id_set = c()
  dd_set = c()
  mv_set = c()
  rv_set = c()
  rd_set = c()
  
  if(n.edge(A) < m){
    
    # set of nodes for iu
    nodes = which(A_und == 0,TRUE)
    if(length(nodes)!=0){
      
    iu_set = cbind(1,nodes,NA)
    
    # set of nodes for id
    id_set = cbind(3,rbind(iu_set[,2:3],iu_set[,3:2]),NA)}
  
  }
  
  # set of nodes for du
  nodes = which(A_und == 2,TRUE)
  if(length(nodes)!=0)
    du_set = cbind(2,nodes,NA)
  
  # set of nodes for dd and rd
  nodes = which(A_dir==1,TRUE)
  if(length(nodes)!=0){
    dd_set = cbind(4,nodes,NA)
    rd_set = cbind(7,nodes,NA)
  }
  
  # set of nodes for mv
  nodes = v.undirect(A)
  if(length(nodes)!=0)
  mv_set = cbind(5,nodes)
  
  # set of nodes for rv
  nodes = v.direct(A)
  if(length(nodes)!=0)
  rv_set = cbind(6,nodes)
  
  O = rbind(iu_set,du_set,id_set,dd_set,mv_set,rv_set,rd_set)
  
  repeat {
    
    a = sample(dim(O)[1],1)
    
    act_to_eval = paste0(types[O[a,1]],"(A=A,c(",as.vector(O[a,2]),",",as.vector(O[a,3]),",",as.vector(O[a,4]),"),G_und=G_und)")
    val         = eval(parse(text=act_to_eval))
    
    if (val != 0){
      break
    }
  }
  
  action = paste0(actions[O[a,1]],"(A=A,c(",as.vector(O[a,2]),",",as.vector(O[a,3]),",",as.vector(O[a,4]),"))")
  A_new = eval(parse(text=action))
  
  A_new = as(A_new,"graphNEL")
  A_new = pcalg::pdag2dag(A_new,keepVstruct=TRUE)$graph
  A_new = dag2essgraph(A_new)
  A_new = as(A_new,"matrix")
  
  return(list(A_new = A_new, O_new = O[a,1], nodes = c(O[a,2],O[a,3],O[a,4])))
  
}
