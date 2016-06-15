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
from_d = as.numeric(as.numeric(args[2]))
to_d = as.numeric(as.numeric(args[3]))

n_folds = 5
n_drugs = length(unique(dataset_triplets[,1]))
n_targets = length(unique(dataset_triplets[,2]))

crf_predictions = rep(NA, nrow(dataset_triplets))
crf_predictions_combined = rep(NA, nrow(dataset_triplets))

load('../kiba_drug_sim.rda')
load('../kiba_target_sim.rda')

drug_adj_mat = drug_sim_mat
target_adj_mat = sim_mat
diag(drug_adj_mat) = 0
diag(target_adj_mat) = 0
for (fold in 1:n_folds){
  
  test_ind = test_folds[[fold]]
  train_ind = setdiff(1:nrow(dataset_triplets),test_ind)
  
  train_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
  train_mat[dataset_triplets[train_ind,c(1,2)]] = dataset_triplets[train_ind,3]
  
  cat('getting MF predictions for train data.. \n')
  mf_preds_train = get_mf_cv(train_mat, 400)
  mf_preds_train_mat = mf_preds_train[[2]]
  
  cat('getting MF predictions for complete matrix..\n')
  mf_preds_all = get_libmf_prediction(train_mat, 800)
  
  cl = makeCluster(4) # 4 means 4 physical core, does not work with hyper-threads.
  registerDoParallel(cl)
  results = foreach (d = from_d:to_d, .packages='Matrix') %dopar% {
    crf_predictions = rep(NA,nrow(dataset_triplets))
#   for (d in from_d : to_d){
    inds = which(dataset_triplets[test_ind,1]==d)
    targets = dataset_triplets[test_ind[inds],2]
    for (t in targets){
      
      cat('drug ',d,' target ',t,'\n')
      
      #drug_inds = c(which(drug_adj_mat[d,]>0),d)
      drug_inds = c(d,order(drug_adj_mat[d,], decreasing = T)[1:40])
      #target_inds = c(which(target_adj_mat[t,]>0),t)
      target_inds = c(t,order(target_adj_mat[t,], decreasing = T)[1:40])
      if (length(which(train_mat[drug_inds,target_inds]>=0))>1){
        
        training_adj_mats = make_training_adj_mat(train_mat = train_mat[drug_inds,target_inds], drug_adj_mat = drug_adj_mat, target_adj_mat = target_adj_mat, drug_set = drug_inds, target_set = target_inds)
        
        inds_ = which(t(train_mat[drug_inds,target_inds]>=0))
        training_data = t(train_mat[drug_inds, target_inds])[inds_]
        X = as.vector(t(mf_preds_train_mat[drug_inds,target_inds])[inds_])
        eta = 0.004
        params = train_crf(training_data = training_data, adj_mats = training_adj_mats, X = X, crf_iters = 1000, eta = eta, 50)
        
        train_data = train_mat[drug_inds, target_inds]
        train_data = as.vector(t(train_data))
        X = as.vector(t(mf_preds_all[drug_inds,target_inds]))
        prediction_adj_mats = make_prediction_adj_mat(train_mat = train_mat[drug_inds, target_inds], drug_adj_mat = drug_adj_mat, target_adj_mat = target_adj_mat, drug_set = drug_inds, target_set = target_inds)
        prediction = crf_predict(train_data = train_data, X = X, alpha = params[[1]], adj_mats = prediction_adj_mats, beta = params[[2]])
        prediction_mat = matrix(prediction, ncol = length(target_inds), byrow = T)
      
        crf_predictions[test_ind[inds[match(t, targets)]]] = prediction_mat[match(d, drug_inds), match(t, target_inds)]
        
      } else {
        crf_predictions[test_ind[inds[match(t, targets)]]] = mf_preds_all[d, t]
      }
      
      inds__ = which(!is.na(crf_predictions))
      crf_metrics = get_metrics(crf_predictions[inds__], dataset_triplets[inds__,3], 12.1)
      #     
      cat('fold ',fold,'\n')
      cat('auc so far, crf: ',round(crf_metrics[[2]], digits = 5),'\n')
      cat('aupr so far, crf: ',round(crf_metrics[[3]], digits = 5),'\n\n')
      
      return(crf_predictions)
    }
  }
  stopCluster(cl)
  for (i in 1:length(results)){
    crf_predictions_combined[which(!is.na(results[[i]]))] = results[[i]][which(!is.na(results[[i]]))]
  }
}

save(crf_predictions_combined, file = args[5])
