load('kiba_data.rda')

load('../data/mix_dataset.Rda')
target_sim = read.table('../data/mix_dataset_target_target_sim.txt')
target_sim = as.matrix(target_sim)


find_max_connected = function(set, sim_mat, thresh){
  
  max_n = 0
  d = -1
  for (i in 1:length(set)){
    nn = length(which(sim_mat[set[i],]>thresh & sim_mat[set[i],]<1))
    if (nn>max_n){
      max_n = nn
      d = set[i]
    }
  }
  
  return(list(d,max_n))
  
}

expand_set = function(set, sim_mat, thresh){
  
  neighbors = c()
  for (i in 1:length(set)){
    neighbors = union(neighbors,which(simi_mat[set[i],]>thresh & simi_mat[set[i],]<1))
  }
  
  new_neighbors = setdiff(neighbors,set)
  
  return (new_neighbors)
  
}

get_n_steps_neighbors = function(temp, sim_mat, it){
  
  for (i in 1:it){
    new_neighbors = expand_set(temp, sim_mat)
    new_neighbors = setdiff(new_neighbors, temp)
    temp = union(temp,new_neighbors)
  }
  
  return(temp)
  
}

get_n_steps_neighbors_ = function(temp, sim_mat, it, dont_touch, thresh){
  
  for (i in 1:it){
    new_neighbors = expand_set(temp, sim_mat, thresh)
    new_neighbors = setdiff(new_neighbors, temp)
    new_neighbors = setdiff(new_neighbors, dont_touch)
    temp = union(temp,new_neighbors)
  }
  
  return(temp)
  
}

#set = unique(dt_triplet[,1])
set = unique(mix_dataset[,2])

simi_mat = target_sim
steps = 2
already_clustered = c()
cluster_list = list()
n_clusters = 0
thresh = 0.4
while(length(set)>0){
  
  ans = find_max_connected(set, sim_mat, thresh)
  d = ans[[1]]
  nn = ans[[2]]
  cat(length(set),' drugs left to cluster\n')
  cat('max connected drug is: ',d,' with ',nn,' neighbors\n')
  if (nn>0){
    
    temp = 0
    cluster = get_n_steps_neighbors_(d, sim_mat, 1, already_clustered, thresh)
    cluster_next = get_n_steps_neighbors_(d, sim_mat, 2, already_clustered, thresh)
    cat('found ',length(cluster_next),' neighbors in ',2,' steps\n')
    steps = 3
    while(length(cluster)!=length(cluster_next)){
      cluster = cluster_next
      cluster_next = get_n_steps_neighbors_(d, sim_mat, steps, already_clustered,thresh)
      steps = steps+1
      cat('found ',length(cluster_next),' neighbors in ',steps,' steps\n')
    }
    
    already_clustered = union(already_clustered, cluster)
    cat('training data of cluster: ',length(which(mix_dataset[,2]%in%cluster)),'\n')
    n_clusters = n_clusters+1
    cluster_list[[n_clusters]] = cluster
    cat(cluster,'\n\n')
    set = setdiff(set,already_clustered)
  }
  else{
    break;
  }
  
}


drugs_a = mix_dataset[which(mix_dataset[,2]==205),1]
drugs_b = mix_dataset[which(mix_dataset[,2]==208),1]
intersect(drugs_a, drugs_b)


