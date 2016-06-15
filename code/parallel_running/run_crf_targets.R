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
filename = '5fold_cv_data_metz.Rda'
load(file = filename)

from = 81
to = 85

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
  
  cat('getting MF predictions for train data.. \n')
  mf_preds_train = get_mf_cv(train_mat, 400)
  mf_preds_train_mat = mf_preds_train[[2]]
  
  #cat('getting MF predictions for complete matrix..\n')
  mf_preds_all = get_libmf_prediction(train_mat, 800)
  
  ## train crf for column 26
#    cl = makeCluster(4) # 4 means 4 physical core, does not work with hyper-threads.
#    registerDoParallel(cl) 
#    results = foreach (t = from:to, .packages='Matrix') %dopar% {
#      crf_predictions = rep(NA, nrow(dataset_triplets))
    
 for (t in from:to){
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
      params = train_crf(training_data = training_data, adj_mats = adj_mats, X = X, crf_iters = 600, eta = 0.004,20)
      
      train_data = train_mat[,t]
#       known = which(train_mat[,t]>=0)
#       m = mean(train_mat[which(train_mat>=0)])
#       inds = which(dataset_triplets[test_ind,2]==t)
#       val_drugs = dataset_triplets[test_ind[inds],1]
#       nan_drugs = setdiff(c(1:n_drugs),c(known,val_drugs))
#       train_data[nan_drugs] = m
      
      adj_mats = list()
      adj_mats[[1]] = drug_adj_mat
      #prediction = crf_predict(train_data = train_mat[,t], X = mf_preds_all[,t], alpha = params[[1]], adj_mats = adj_mats, beta = params[[2]])
      prediction = crf_predict(train_data = train_data, X = mf_preds_all[,t], alpha = params[[1]], adj_mats = adj_mats, beta = params[[2]])
      
      inds = which(dataset_triplets[test_ind,2] == t)
      crf_predictions[test_ind[inds]] = prediction[dataset_triplets[test_ind[inds],1]]
      
      
    } else {
      inds = which(dataset_triplets[test_ind,2] == t)
      crf_predictions[test_ind[inds]] = mf_preds_all[dataset_triplets[test_ind[inds],1],t]
    }
    
#    inds = which(!is.na(crf_predictions))
#    crf_metrics = get_metrics(crf_predictions[inds], dataset_triplets[inds,3], 7.6)
    
#    cat('auc so far, crf: ',round(crf_metrics[[2]], digits = 5),'\n')
#    cat('aupr so far, crf: ',round(crf_metrics[[3]], digits = 5),'\n\n')

#      return(crf_predictions)
  }
#    stopCluster(cl)
# 	for (i in 1:length(results)){
#      		crf_predictions_combined[which(!is.na(results[[i]]))] = results[[i]][which(!is.na(results[[i]]))]
#    	}
}


save(crf_predictions_combined, file = args[4])
