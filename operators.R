# The following functions test the validity of each type of operator according to He et al. (2013)

# dd_1 condition (validity of DeleteD)

dd_1 <- function(A,x,y){
  neigh_y = find.neib(y,A)
  l.neigh = length(neigh_y)
  if(l.neigh <= 1){
    return(TRUE)} else{
      # is a clique?
      if(sum(A[neigh_y,neigh_y])<l.neigh*(l.neigh-1))
        return(FALSE)}
          return(TRUE)
}

# Perfect operator dd

perfect.dd <- function(A,nodes,G_und){
  x = nodes[1]
  y = nodes[2]
  if(dd_1(A,x,y)) TRUE else FALSE
}


# id_2 condition (validity of InsertD)

id_2_1 <- function(A,x,y){
  
  par_x = find.parents(x,A)
  par_y = find.parents(y,A)
  
  if(setequal(par_x,par_y)) FALSE else TRUE
  
}

id_2_2 <- function(A,x,y){
  
  neigh_y = find.neib(y,A)
  par_x = find.parents(x,A)
  Omega_xy = intersect(par_x,neigh_y)
  # is a clique?
  length.omega = length(Omega_xy)
  if(length.omega>1){
    if(sum(A[Omega_xy,Omega_xy])<length.omega*(length.omega-1)) FALSE else TRUE
  } else {return(TRUE)}
  
}

id_2_3 <- function(A,x,y){
  
  A_dir = A
  A_dir[which(A == 1 & t(A) == 1,TRUE)] <- 0
  
  G_dir = as(A_dir,"graphNEL")
  G_dir = as(G_dir,"igraph")
  
  allpath_xy = all_simple_paths(G_dir,y,x)
  
  # I have to exclude the "all" directed paths?
  # dirpath_xy = all_simple_paths(G_dir,y,x)
  # pdirpath_xy = setdiff(allpath_xy,dirpath_xy)
  
  pdirpath_xy = allpath_xy
  
  if(length(pdirpath_xy) == 0){
    return(TRUE)
  } else{
    
    neigh_y = find.neib(y,A)
    par_x = find.parents(x,A)
    Omega_xy = intersect(par_x,neigh_y)
    
    len = sapply(pdirpath_xy,function(i){
      length(intersect(i,Omega_xy))>0})
    if(all(len)) TRUE else FALSE
  }
} 

id_2 <- function(A,x,y){
    if(id_2_1(A,x,y) & id_2_2(A,x,y) & id_2_3(A,x,y)) TRUE else FALSE
}

# Perfect operator id

perfect.id <- function(A,nodes,G_und){
  x = nodes[1]
  y = nodes[2]
  
  if(id_2(A,x,y)) TRUE else FALSE
  
}


G_undirected = function(A){
  A_und = matrix(0,nrow(A),ncol(A))
  A_und[which(A == 1 & t(A) == 1,TRUE)]<-1
  G_und = as(A_und,"graphNEL")
  G_und = as(G_und,"igraph")  
}


# iu_2 condition (validity of InsertU)

# use the function all_simple_paths from 'igraph'
# require igraph

iu_2 <- function(A,x,y,G_und){
  
  par_x = find.parents(x,A)
  par_y = find.parents(y,A)
  
  if(setequal(par_x,par_y) == FALSE){
    return(FALSE)
  } else{
    
    neigh_x = find.neib(x,A)
    neigh_y = find.neib(y,A)
    neigh_xy = intersect(neigh_x,neigh_y)
    
    if(length(neigh_xy) == 0){
      return(TRUE)} else{
        
        undpath_xy = all_simple_paths(G_und,x,y)
        
        len = lapply(undpath_xy,function(i)
          length(intersect(i,neigh_xy)))
        
        if(all(len != 0)) TRUE else FALSE
      }
  }
}

# Perfect operator iu

perfect.iu <- function(A,nodes,G_und){
  x = nodes[1]
  y = nodes[2]
  if(iu_2(A,x,y,G_und)) TRUE else FALSE
}


# du_1 condition

# A   the adjacency matrix of the CPDAG in input

du_1 <- function(A,x,y){   # to be tested
  neigh_x = find.neib(x,A)
  neigh_y = find.neib(y,A)
  neigh_xy = intersect(neigh_x,neigh_y)
  length.xy = length(neigh_xy)
  if(length.xy <= 1){
    return(TRUE)} else{
      # is a clique?
      if(sum(A[neigh_xy,neigh_xy])<length.xy*(length.xy-1))
        return(FALSE)
    }
  return(TRUE)
}

# Perfect operator du

perfect.du <- function(A,nodes,G_und){
  x = nodes[1]
  y = nodes[2]
  if(du_1(A,x,y)) TRUE else FALSE
}


# mv_1 condition

# A   the adjacency matrix of the CPDAG in input
# (x,z,y) belongs to the set MakeV, so that  x - z - y
# we apply x -> z <- y

mv_1 <- function(A,x,z,y,G_und){
  
  neigh_x  = find.neib(x,A)
  neigh_y  = find.neib(y,A)
  neigh_xy = intersect(neigh_x,neigh_y)
  if(length(neigh_xy) <= 2){
    return(TRUE)
  } else{
    
    undpath_xy = all_simple_paths(G_und,x,y)
    len = sapply(undpath_xy,function(i){
      length(intersect(i,neigh_xy))>0})
    if(all(len)) TRUE else FALSE
  }
}

# Perfect operator mv

perfect.mv <- function(A,nodes,G_und){
  x = nodes[1]
  y = nodes[3]
  z = nodes[2]
  if(mv_1(A,x,z,y,G_und)) TRUE else FALSE
}


# rd_2 condition (validity of ReverseD) # See Chickering 2002 (a)

# Perfect operator rd

perfect.rd <- function(A,nodes,G_und){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 0 # RemoveD
  if(id_2(A,y,x)) TRUE else FALSE # InsertD x <- y
}


# rv_1

# A   the adjacency matrix of the CPDAG in input
# (x,z,y) belongs to the set RemoveV, so that  x -> z <- y
# we apply x - z - y

rv_1 <- function(A,x,z,y){
  
  par_x = find.parents(x,A)
  par_y = find.parents(y,A)
  
  if(setequal(par_x,par_y)) TRUE else FALSE
  
}


# rv_2

rv_2 <- function(A,x,z,y){
  
  par_x = find.parents(x,A)
  par_z = find.parents(z,A)
  neigh_x = find.neib(x,A)
  neigh_y = find.neib(x,A)
  neigh_xy = intersect(neigh_x,neigh_y)
  
  if(setequal(union(par_x,setdiff(neigh_xy,z)),setdiff(par_z,c(x,y)))) TRUE else FALSE
}


# rv_3

rv_3 <- function(A,x,z,y,G_und){
  neigh_x = find.neib(x,A)
  neigh_y = find.neib(x,A)
  neigh_xy = intersect(neigh_x,neigh_y)
  
  if(length(neigh_xy) == 0){
    return(TRUE)} else{
      
      undpath_xy = all_simple_paths(G_und,x,y)
      len = lapply(undpath_xy,function(i)
        length(intersect(i,neigh_xy)))
      
      if(all(len != 0)) TRUE else FALSE
      
    }
}


# Perfect operator rv

perfect.rv <- function(A,nodes,G_und){
  x = nodes[1]
  y = nodes[3]
  z = nodes[2]
  if(rv_1(A,x,z,y) == TRUE & rv_2(A,x,z,y) == TRUE & rv_3(A,x,z,y,G_und) == TRUE) TRUE else FALSE
}


# 6 actions

# A       adjacency matrix of CPDAG C
# x,y,z   nodes involved in the action

iu <- function(A,nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = A[y,x] = 1
  return(A)
}

du <- function(A,nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = A[y,x] = 0
  return(A)
}

id <- function(A,nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 1
  return(A)
}

dd <- function(A,nodes){
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 0
  return(A)
}

mv <- function(A,nodes){
  x = nodes[1]
  y = nodes[3]
  z = nodes[2]
  A[z,x] = A[z,y] = 0
  return(A)
}

rv <- function(A,nodes){ # check!
  x = nodes[1]
  y = nodes[3]
  z = nodes[2]
  A[z,x] = A[z,y] = 1
  return(A)
}

rd <- function(A,nodes){ # reverse D x -> y
  x = nodes[1]
  y = nodes[2]
  A[x,y] = 0 
  A[y,x] = 1
  return(A)
}


# Sets (each operator can be applied to)

# A adjacency matrix of CPDAG C

# A_low<-A
# A_upp<-A

# A_low[upper.tri(A,diag=TRUE)]<-NA
# A_upp[lower.tri(A,diag=TRUE)]<-NA

# diag(A_upp)<-diag(A_low)<-NA
# A_und<-A_low+t(A_upp)

# iu_set <- which(A_und == 0,TRUE)
# dd_set <- which(A_und == 2,TRUE)
# id_set <- which(A_und == 0,TRUE)
# dd_set <- which(A == 1 & t(A) == 0,TRUE)
# mv_set <- v.undirect(A)
# rv_set <- v.direct(A)
# rd_set <- dd_set

# To find v-structures (directed)

# A adjacency matrix of a graph

v.direct<-function(A){
  Av = matrix(c(0,1,0,0,0,0,0,1,0),byrow=T,3,3)
  v.str = c()
  for(i in 1:nrow(A)){
    for(j in 1:i){
      common = intersect(which(as.integer(A[i,]) != 0),which(as.integer(A[j,]) != 0))
      for(k in common){
        if(i != j & all(A[c(i,k,j),c(i,k,j)] == Av)){
          v.str = list(v.str,c(i,k,j))
        }
      }
    }
  }
  if(is.null(v.str)){
    NULL
  } else{
    matrix(unlist(v.str),byrow=T,ncol=3)}
}

# To find undirected v-structures

# A adjacency matrix of a graph

v.undirect<-function(A){
  Avu = matrix(c(0,1,0,1,0,1,0,1,0),byrow=T,3,3)
  v.und = c()
  for(i in 1:nrow(A)){
    for(j in 1:i){
      common = intersect(which(as.integer(A[i,]) != 0),which(as.integer(A[j,]) != 0))
      for(k in common){
        if(i != j & all(A[c(i,k,j),c(i,k,j)] == Avu)){
          v.und = list(v.und,c(i,k,j))
        }
      }
    }
  }
  if(is.null(v.und)){
    NULL
  } else{
    matrix(unlist(v.und),byrow=T,ncol=3)}
}


# To test strong protection

# 4 configurations

# An arrow is strongly potected if it occurs in at least one of the following four configurations (Andersson 1997)
# The test is on x -> y

# auxiliary functions
# to find parents
find.parents = function(x,A){
  which(as.integer(A[,x]) == 1 & as.integer(A[x,]) == 0)
} 
find.sons = function(x,A){
  which(as.integer(A[x,]) == 1 & as.integer(A[,x]) == 0)
} 
find.neib = function(x,A){
  which(as.integer(A[x,]) == 1 & as.integer(A[,x]) == 1)
} 

# Conf. a

protect.a <- function(x,y,A){
  par_x = which(as.integer(A[,x]) == 1 & as.integer(A[x,]) == 0)
  if(length(par_x) == 0){
    return(FALSE)} else{
      par_y = which(as.integer(A[,y]) == 1 & as.integer(A[y,]) == 0)
      if(length(intersect(par_x,par_y)) != length(setdiff(union(par_x,par_y),x))){
        return(TRUE)} else{
          return(FALSE)}
    }
}

# Conf. b

protect.b <- function(x,y,A){
  par_y = setdiff(which(as.integer(A[,y]) == 1 & as.integer(A[y,]) == 0),x)
  if(length(par_y) == 0){
    return(FALSE)} else{
      if(length(which((A[par_y,x]+A[x,par_y]) == 0)) != 0){
        return(TRUE)} else{
          return(FALSE)}
    }
}

# Conf. c

protect.c <- function(x,y,A){
  par_y = setdiff(which(as.integer(A[,y]) == 1 & as.integer(A[y,]) == 0),x)
  chi_x = setdiff(which(as.integer(A[x,]) == 1 & as.integer(A[,x]) == 0),y)
  if(length(par_y) == 0 & length(chi_x) == 0){
    return(FALSE)} else{
      if(length(intersect(par_y,chi_x)) != 0){
        return(TRUE)} else{
          return(FALSE)}
    }
}

# Conf. d

protect.d <- function(x,y,A){
  par_y = setdiff(which(as.integer(A[,y]) == 1 & A[y,] == 0),x)
  chi_x = intersect(which(as.integer(A[x,]) == 1),which(as.integer(A[,x]) == 1))
  if(length(par_y) == 0|length(chi_x) == 0|length(intersect(par_y,chi_x)) <= 1){
    return(FALSE)} else{
      int = intersect(par_y,chi_x)
      C = A[int,int]
      if(length(which(C+t(C) != 0)) != length(int)*(length(int)-1)){ # at least a non zero entry
        return(TRUE)} else{
          return(FALSE)
        }
    }
}


# Strong protection function

strong.protect <- function(x,y,A){
  if(protect.a(x,y,A) == TRUE | protect.b(x,y,A) == TRUE |
     protect.c(x,y,A) == TRUE | protect.d(x,y,A) == TRUE){
    return(TRUE)} else{
      return(FALSE)
    }
}
