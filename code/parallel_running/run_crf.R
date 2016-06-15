source('crf_utils.R')
source('recosystem.R')
source('target_clustering.R')
require('Matrix')
require(foreach)
require(doParallel)
require(parallel)

load(file = "./5cv_data.Rda")

args = commandArgs(trailingOnly = T)
from = as.numeric(as.numeric(args[1]))
to = as.numeric(as.numeric(args[2]))
from = 1
to = 4

n_folds = 5
n_drugs = length(unique(dataset_triplets[,1]))
n_targets = length(unique(dataset_triplets[,2]))

for (fold in 1:n_folds){
  
  test_ind = test_folds[[fold]]
  train_ind = setdiff(1:nrow(dataset_triplets),test_ind)
  
  train_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
  train_mat[dataset_triplets[train_ind,c(1,2)]] = dataset_triplets[train_ind,3]
  
  cat('getting MF predictions for train data.. \n')
  mf_preds_train = get_mf_cv(train_mat, 400)
  mf_preds_train_mat = mf_preds_train[[2]]
  
  cat('getting MF predictions for complete matrix..\n')
  mf_preds_all = get_libmf_prediction(train_mat, 400)
  
  ## train crf for column 26
  
  for (t in 1:n_targets){
    
    cat('target ',t,'... ',length(which(train_mat[,t]>=0)),' observations\n')
    
    if (length(which(train_mat[,t]>=0))>1){
      
      training_data = train_mat[which(train_mat[,t]>0), t]
      inds = which(train_mat[,t]>=0)
      training_adj_mat = drug_adj_mat[inds,inds]
      inds = which(training_adj_mat>0, arr.ind = T)
      training_adj_mat = sparseMatrix(i = inds[,1], j = inds[,2], x = training_adj_mat[inds], dims = dim(training_adj_mat))
      X = mf_preds_train_mat[which(!is.na(mf_preds_train_mat[,t])),t]
      
      adj_mats = list()
      adj_mats[[1]] = training_adj_mat
      params = train_crf(training_data = training_data, adj_mats = adj_mats, X = X, crf_iters = 600, eta = 0.004)
      
      adj_mats = list()
      adj_mats[[1]] = drug_adj_mat
      prediction = crf_predict(train_data = train_mat[,t], X = mf_preds_all[,t], alpha = params[[1]], adj_mats = adj_mats, beta = params[[2]])
      
      inds = which(dataset_triplets[test_ind,2] == t)
      crf_predictions[test_ind[inds]] = prediction[dataset_triplets[test_ind[inds],1]]
      
      
    } else {
      inds = which(dataset_triplets[test_ind,2] == t)
      crf_predictions[test_ind[inds]] = mf_preds_all[dataset_triplets[test_ind[inds],1],t]
    }
    
    inds = which(!is.na(crf_predictions))
    crf_metrics = get_metrics(crf_predictions[inds], dataset_triplets[inds,3], 7.6)
    
    cat('auc so far, crf: ',round(crf_metrics[[2]], digits = 5),'\n')
    cat('aupr so far, crf: ',round(crf_metrics[[3]], digits = 5),'\n\n')
    
  }
  
}