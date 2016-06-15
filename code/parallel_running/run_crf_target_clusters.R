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
n_neighbors = as.numeric(as.numeric(args[4]))

filename = '5fold_cv_data_metz.Rda'
load(file = filename)
n_neighbors = 3
from = 19
to = 20


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
  mf_preds_train = get_mf_cv(train_mat, 800)
  mf_preds_train_mat = mf_preds_train[[2]]
  
  #cat('getting MF predictions for complete matrix..\n')
  mf_preds_all = get_libmf_prediction(train_mat, 800)
  
  # parallization by foreach and doparallel
#   cl = makeCluster(4) # 4 means 4 physical core, does not work with hyper-threads.
#   registerDoParallel(cl)
#   results = foreach (t = from:to, .packages='Matrix') %dopar% {
# 	crf_predictions = rep(NA,nrow(dataset_triplets))
    for (t in 1:n_targets){
    target_set = get_nn(c = t, adj_mat = target_adj_mat, max_n = n_neighbors)
    ## if target t has at least one neighbor, train and predict with a crf that integrates both similarities
    if (length(target_set)>1){
      
      # cat('training crf for targets ',target_set,'..\n')
      ## train_data is the available data for the target set
      ## X is the prediction of MF on this data
      ## training_adj_mats contains two matrices of dimension n_train x n_train
      ## the first matrix connects all cells of the training data that are drug neighbors
      ## the second matrix connects all cells of the training data that are target neighbors
      
	    train_data = t(train_mat[,target_set])[which(t(train_mat[,target_set])>=0)]
      X = t(mf_preds_train_mat[,target_set])[which(!is.na(t(mf_preds_train_mat[,target_set])))]
      training_adj_mats = make_training_adj_mat(train_mat = train_mat[,target_set], drug_adj_mat = drug_adj_mat, target_adj_mat = target_adj_mat, drug_set = c(1:n_drugs), target_set = target_set)
      
	    if (length(train_data)>2000){
		    eta = 0.001
	    } else {
		    eta = 0.004
	    }

	    params = train_crf(training_data = train_data, adj_mats = training_adj_mats, X = X, crf_iters = 600, eta = eta, 1)
      
      ## X is the vector of all predictions on the target set
      ## prediction_adj_mats contains two matrices of dimension (n_drugs*length(target_set)^2)
      ## the first matrix connects cells that are drug-neighbors
      ## the second matrix connects cells that are target-neighbors
      X = as.vector(t(mf_preds_all[,target_set]))
      prediction_adj_mats = make_prediction_adj_mat(train_mat = train_mat[,target_set], drug_adj_mat = drug_adj_mat, target_adj_mat = target_adj_mat, drug_set = c(1:n_drugs), target_set = target_set)

	    train_data = train_mat[,target_set]
#        for (t in target_set){
#          known = which(train_mat[,t]>=0)
#          m = mean(train_mat[which(train_mat>=0)])
#          inds = which(dataset_triplets[test_ind,2]==t)
#          val_drugs = dataset_triplets[test_ind[inds],1]
#          nan_drugs = setdiff(c(1:n_drugs),c(known,val_drugs))
#          train_data[nan_drugs, match(t, target_set)] = m
#        }
       train_data = as.vector(t(train_data))

      #prediction = crf_predict(train_data = as.vector(t(train_mat[,target_set])), X = X, alpha = params[[1]], adj_mats = prediction_adj_mats, beta = params[[2]])

	    prediction = crf_predict(train_data = train_data, X = X, alpha = params[[1]], adj_mats = prediction_adj_mats, beta = params[[2]])
      
      prediction_mat = matrix(prediction, nrow = n_drugs, byrow = T)
      
      ## get the prediction of the targets in target_set that are in the test set
      inds = which(dataset_triplets[test_ind,2]==t)
      drugs = dataset_triplets[test_ind[inds],1]
      labels = dataset_triplets[test_ind[inds],3]
      crf_predictions_test_col = prediction_mat[drugs, match(t,target_set)]
      crf_predictions[test_ind[inds]] = crf_predictions_test_col
      
    } else {
      ## else integrate only the drug similarity
      
      # cat('target ',t,'... ',length(which(train_mat[,t]>=0)),' observations\n')
      if (length(which(train_mat[,t]>=0))>1){
        
        training_data = train_mat[which(train_mat[,t]>=0), t]
        inds = which(train_mat[,t]>=0)
        training_adj_mat = drug_adj_mat[inds,inds]
        inds = which(training_adj_mat>0, arr.ind = T)
        training_adj_mat = sparseMatrix(i = inds[,1], j = inds[,2], x = training_adj_mat[inds], dims = dim(training_adj_mat))
        X = mf_preds_train_mat[which(!is.na(mf_preds_train_mat[,t])),t]
        
        adj_mats = list()
        adj_mats[[1]] = training_adj_mat
        params = train_crf(training_data = training_data, adj_mats = adj_mats, X = X, crf_iters = 600, eta = 0.004, 100)
        
        adj_mats = list()
        adj_mats[[1]] = drug_adj_mat

	

        prediction = crf_predict(train_data = train_mat[,t], X = mf_preds_all[,t], alpha = params[[1]], adj_mats = adj_mats, beta = params[[2]])
        
        inds = which(dataset_triplets[test_ind,2] == t)
        crf_predictions[test_ind[inds]] = prediction[dataset_triplets[test_ind[inds],1]]
        
      } else {
        inds = which(dataset_triplets[test_ind,2] == t)
        crf_predictions[test_ind[inds]] = mf_preds_all[dataset_triplets[test_ind[inds],1],t]
      }
    }
    
        inds = which(!is.na(crf_predictions))
        crf_metrics = get_metrics(crf_predictions[inds], dataset_triplets[inds,3], 7.6)
        
         cat('fold ',fold,'\n')
         cat('auc so far, crf: ',round(crf_metrics[[2]], digits = 5),'\n')
         cat('aupr so far, crf: ',round(crf_metrics[[3]], digits = 5),'\n\n')
    
#     return(crf_predictions)
   }
#   stopCluster(cl)
#   
#   for (i in 1:length(results)){
#     crf_predictions_combined[which(!is.na(results[[i]]))] = results[[i]][which(!is.na(results[[i]]))]
#   }
  
}

#inds = which(!is.na(crf_predictions_combined))
#crf_metrics = get_metrics(crf_predictions_combined[inds], dataset_triplets[inds,3], 7.6)

save(crf_predictions_combined, file = args[5])
