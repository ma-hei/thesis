make_training_adj_mat = function(train_mat, drug_adj_mat, target_adj_mat, drug_set, target_set){
  
  inds = which(t(train_mat>=0))
  
  adj_mats = list()
  for (i in inds){
  #i = inds[1]
    adj_mats[[match(i,inds)]] = matrix(0, nrow = length(inds), ncol = length(inds))
    for (k in setdiff(inds,i)){
      
      drug_i = drug_set[floor(i/length(target_set))+1]
      drug_k = drug_set[floor(i/length(target_set))+1]
      temp = i%%length(target_set)
      if (temp == 0){
        target_i = target_set[length(target_set)]
      } else{
        target_i = target_set[temp]
      }
      temp = k%%length(target_set)
      if (temp == 0){
        target_k = target_set[length(target_set)]
      } else{
        target_k = target_set[temp]
      }
      
      if (drug_i==drug_k){
        similarity = target_adj_mat[target_i, target_k]
      } else if (target_i==target_k){
        similarity = drug_adj_mat[target_i, target_k]
      } else{
        similarity = target_adj_mat[target_i, target_k] * drug_adj_mat[target_i, target_k]
      }
      
      adj_mats[[match(i,inds)]][match(i,inds), match(k,inds)] = similarity
      adj_mats[[match(i,inds)]][match(k,inds), match(i,inds)] = similarity
    }
  }
  
  return(adj_mats)
  
}

make_prediction_adj_mats = function(train_mat, train_mask, drug_adj_mat, target_adj_mat, drug_set, target_set){
  
  inds = which(train_mask>=0)
  kron = kronecker(drug_adj_mat[drug_set,drug_set], target_adj_mat[target_set, target_set])
  
  pred_adj_mats = list()
  for (i in inds){
    drug_i = floor(i/length(target_set))+1
    temp = i%%length(target_set)
    if (temp == 0){
      target_i = length(target_set)
    } else{
      target_i = temp
    }
    
    pred_adj_mats[[match(i,inds)]] = matrix(0, nrow = length(train_mat), ncol = length(train_mat))
    for (d in 1:nrow(train_mat)){
      temp_i = (d-1)*length(target_set)+1
      temp_j = (drug_i-1)*length(target_set)+1
      temp_vec = kron[temp_i:(temp_i+length(target_set)-1), temp_j:(temp_j+length(target_set)-1)]
      pred_adj_mats[[match(i,inds)]][temp_i:(temp_i+length(target_set)-1), temp_j:(temp_j+length(target_set)-1)] = temp_vec
      pred_adj_mats[[match(i,inds)]][temp_j:(temp_j+length(target_set)-1), temp_i:(temp_i+length(target_set)-1)] = temp_vec
    }
  }
  
  return (pred_adj_mats)
}
