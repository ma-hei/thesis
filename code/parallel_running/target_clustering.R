get_nn = function(c, adj_mat, max_n){
  
  temp = which(adj_mat[c,]>0)
  if (length(temp)>max_n){
    neighbors = order(adj_mat[c,], decreasing = T)[1:max_n]
    neighbors = union(neighbors,c)
    return (neighbors)
  } else {
    if (length(temp)==0){
      return(c)
    } else{
      neighbors = union(temp,c)
      return(neighbors)
    }
  }
  
}

get_cluster = function(c, adj_mat, max_n){
  
  target_cluster = c
  next_c = c
  already_queried = c()
  while (length(target_cluster)<max_n && next_c!=-1){
    
    #cat('querying ',next_c,'\n')
    
    temp = get_nn(next_c, adj_mat, 10)
    target_cluster = union(target_cluster, temp)
    already_queried = c(already_queried, next_c)
    
    temp = setdiff(target_cluster, already_queried)
    next_c = temp[1]
    
    #cat('cluster : ', target_cluster,'\n')
    #cat('length :', length(target_cluster),'\n')
    
  }
  return(target_cluster)
}