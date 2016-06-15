#source('/vol/cluster-data/mheideme/parallel_running/crf_utils.R')
source('crf_utils.R')
#source('/vol/cluster-data/mheideme/parallel_running/recosystem.R')
source('recosystem.R')
#source('/vol/cluster-data/mheideme/parallel_running/target_clustering.R')
source('target_clustering.R')
require('Matrix')
require(foreach)
require(doParallel)
require(parallel)

args = commandArgs(trailingOnly = T)
filename = args[1]
from = as.numeric(as.numeric(args[2]))
to = as.numeric(as.numeric(args[3]))

filename = '5fold_cv_data_kiba.Rda'

load(file = filename)

n_folds = 5
n_drugs = length(unique(dataset_triplets[,1]))
n_targets = length(unique(dataset_triplets[,2]))

crf_predictions = rep(NA, nrow(dataset_triplets))
crf_predictions_combined = rep(NA, nrow(dataset_triplets))

for (fold in 1:n_folds){
  
  test_ind = test_folds[[fold]]
  train_ind = setdiff(1:nrow(dataset_triplets),test_ind)
  
  train_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
  train_mat[dataset_triplets[train_ind,c(1,2)]] = dataset_triplets[train_ind,3]
  
  #cat('getting MF predictions for train data.. \n')
  mf_preds_train = get_mf_cv(train_mat, 400)
  mf_preds_train_mat = mf_preds_train[[2]]
  
  #cat('getting MF predictions for complete matrix..\n')
  mf_preds_all = get_libmf_prediction(train_mat, 400)
#   
#   cl = makeCluster(4) # 4 means 4 physical core, does not work with hyper-threads.
#   registerDoParallel(cl) 
#   results = foreach (d = from:to, .packages='Matrix') %dopar% {
    for (d = from:to){
    
    if (length(which(train_mat[d,]>=0))>1){
      
      training_data = train_mat[d, which(train_mat[d,]>=0)]
      inds = which(train_mat[d,]>=0)
      training_adj_mat = target_adj_mat[inds,inds]
      inds = which(training_adj_mat>0, arr.ind = T)
      training_adj_mat = sparseMatrix(i = inds[,1], j = inds[,2], x = training_adj_mat[inds], dims = dim(training_adj_mat))
      X = mf_preds_train_mat[d, which(!is.na(mf_preds_train_mat[d,]))]
      
      adj_mats = list()
      adj_mats[[1]] = training_adj_mat
      params = train_crf(training_data = training_data, adj_mats = adj_mats, X = X, crf_iters = 600, eta = 0.002, 10)

      adj_mats = list()
      adj_mats[[1]] = target_adj_mat

      prediction = crf_predict(train_data = train_mat[d,], X = mf_preds_all[d,], alpha = params[[1]], adj_mats = adj_mats, beta = params[[2]])

      inds = which(dataset_triplets[test_ind,1] == d)
      crf_predictions[test_ind[inds]] = prediction[dataset_triplets[test_ind[inds],2]]
      
    } else {
      inds = which(dataset_triplets[test_ind,1] == d)
      crf_predictions[test_ind[inds]] =  mf_preds_all[d, dataset_triplets[test_ind[inds],2]]
    }
    
    return(crf_predictions)
  }
 
  stopCluster(cl)
  for (i in 1:length(results)){
  	crf_predictions_combined[which(!is.na(results[[i]]))] = results[[i]][which(!is.na(results[[i]]))]
  }

  
}

save(crf_predictions_combined, file = args[4])
