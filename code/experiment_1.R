load('kiba_data.rda')

n_drugs = length(unique(dt_triplet[,1]))
n_targets = length(unique(dt_triplet[,2]))

dt_triplet[,3] = dt_triplet[,3]*-1
dt_triplet[,3] = dt_triplet[,3] - min(dt_triplet[,3]) 
hist(dt_triplet[,3])

n_folds = 5
n = nrow(dt_triplet)
folds = vector(n_folds, mode='list')
shuf = sample(n)
m = n %/% n_folds
for (i in 1:(n_folds-1)) {
  folds[[i]] = shuf[1:m]
  shuf = setdiff(shuf, folds[[i]])
}
folds[[n_folds]] = shuf
preds = rep(0,n)

adj_mat_pred = matrix(0, nrow = n_drugs, ncol = n_drugs)
for (i in 1:n_drugs){
  inds = which(simi_mat[i,]>0.9 & simi_mat[i,]<1)
  adj_mat_pred[inds,] = 1
  adj_mat_pred[,inds] = 1
}

for (i in 1:n_folds){
  cat('fold ',i,'\n')
  teind = folds[[i]]
  trind = setdiff(1:n,teind)
  
  X = get_mf_cv(dt_triplet[trind,])
  
  dt_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
  dt_mat[dt_triplet[trind,c(1,2)]] = dt_triplet[trind,3]
  
  adj_mats_train = list()
  for (t in 1:n_targets){
    adj_mats_train[[t]] = create_adj_mat_for_column(dt_mat, simi_mat, t)
    cat('column ',t,' ',length(which(adj_mats[[t]]>0))/2,' number of edges to train\n')
  }
  
  params = train_crf_by_row_2(train = dt_mat, adj_list = adj_mats, X = X, crf_iters = 200, eta = 0.001)
  
  X = get_libmf_prediction(dt_mat, 300)
  
  rmse.val = sqrt(mean((dt_triplet[teind,3]-X[dt_triplet[teind,c(1:2)]])^2))
  
  pred_mat = make_crf_predictions_by_row(params[[1]], params[[2]], dt_mat, adj_mat_pred, X)
  
  preds[teind] = pred_mat[dt_triplet[teind,c(1:2)]]
  
  rmse.val = sqrt(mean((dt_triplet[teind,3]-preds[teind])^2))
  
  cat('rmse: ', rmse.val,'\n')
  
}