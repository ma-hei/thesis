load('kiba_data.rda')

dim(simi_mat)

n_drugs = length(unique(dt_triplet[,1]))
n_targets = length(unique(dt_triplet[,2]))

dt_triplet[,3] = dt_triplet[,3]*-1
dt_triplet[,3] = dt_triplet[,3] - min(dt_triplet[,3]) 
hist(dt_triplet[,3])

### run the code in mf_tests.R to get the libmf predictions of 5 fold cv on the training data
X = preds
###

### create a matrix out of the triplets
dt_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
dt_mat[dt_triplet[,c(1,2)]] = dt_triplet[,3]
###

### create the adjacency matrix.. this doesn't work
### adj_mat = matrix(0, nrow = n_drugs*n_targets, ncol = n_drugs*n_targets)

create_adj_mat_for_column = function(dt_mat, simi_mat, c){
  n_train = length(which(dt_mat[,c]>0))
  adj_temp = matrix(0, nrow = n_train, ncol = n_train)
  temp = which(dt_mat[,c]>0)
  for (k in temp){
    #cat('got observations for drug ',k,'\n')
    neighbors = which(simi_mat[k,]>0.9 & simi_mat[k,]<1)
    #cat('which is adj to ',neighbors,'\n')
    for (n in neighbors){
      if (n %in% temp){
        #cat('neigbor ',n,' was also observed for this target\n')
        #cat(match(n, temp),'\n')
        adj_temp[match(k, temp), match(n, temp)] = 1
        adj_temp[match(n, temp), match(k, temp)] = 1
      }
    }
  }
  return (adj_temp)
}

adj_mats = list()
for (i in 1:n_targets){
  adj_mats[[i]] = create_adj_mat_for_column(dt_mat, simi_mat, i)
  cat('column ',i,' ',length(which(adj_mats[[i]]>0))/2,' number of edges to train\n')
}

train_crf_by_row_2(train = dt_mat, adj_list = adj_mats, X = X, crf_iters = 10000, eta = 0.001)