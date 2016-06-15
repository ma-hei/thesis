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

filename = '5fold_cv_data_kiba.Rda'
load(file = filename)

drug_clusters = hclust(dist(drug_sim_mat))
target_clusters = hclust(dist(sim_mat))

clusterCut_drugs = cutree(drug_clusters,200)
clusterCut_targets = cutree(target_clusters,10)

n_folds = 5
n_drugs = length(unique(dataset_triplets[,1]))
n_targets = length(unique(dataset_triplets[,2]))

crf_predictions = rep(NA, nrow(dataset_triplets))
crf_predictions_combined = rep(NA, nrow(dataset_triplets))

target_adj_mat_below = target_adj_mat
target_adj_mat_below[which(target_adj_mat>=0.5)] = 0
target_adj_mat_above = target_adj_mat
target_adj_mat_above[which(target_adj_mat<0.5)] = 0

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
  
  # parallization by foreach and doparallel
  #    cl = makeCluster(4) # 4 means 4 physical core, does not work with hyper-threads.
  #    registerDoParallel(cl)
  #    results = foreach (d = from:to, .packages='Matrix') %dopar% {
  #      crf_predictions = rep(NA,nrow(dataset_triplets))
  for (drug_cluster in 1:200){
    for (target_cluster in 1:10){
      
      drug_set = which(clusterCut_drugs==drug_cluster)
      target_set = which(clusterCut_targets==target_cluster)
      
      cat('drug set: ',drug_set,'\n')
      cat('target set: ',target_set,'\n')
      
      train_data = t(train_mat[drug_set,target_set])[which(t(train_mat[drug_set,target_set])>=0)]
      X = t(mf_preds_train_mat[drug_set,target_set])[which(!is.na(t(mf_preds_train_mat[drug_set,target_set])))]
      
      training_adj_mats_1 = make_training_adj_mat(train_mat = train_mat[drug_set,target_set], drug_adj_mat = drug_adj_mat, target_adj_mat = target_adj_mat_below, drug_set = drug_set, target_set = target_set)
      training_adj_mats_2 = make_training_adj_mat(train_mat = train_mat[drug_set,target_set], drug_adj_mat = drug_adj_mat, target_adj_mat = target_adj_mat_above, drug_set = drug_set, target_set = target_set)
      
      adj_mats = list()
      adj_mats[[1]] = training_adj_mats_1[[1]]
      adj_mats[[2]] = training_adj_mats_1[[2]]
      adj_mats[[3]] = training_adj_mats_2[[2]]
      
      eta = 0.001
      
      params = train_crf(training_data = train_data, adj_mats = adj_mats, X = X, crf_iters = 5000, eta = eta, 100)
      
      ## X is the vector of all predictions on the target set
      ## prediction_adj_mats contains two matrices of dimension (n_drugs*length(target_set)^2)
      ## the first matrix connects cells that are drug-neighbors
      ## the second matrix connects cells that are target-neighbors
      X = as.vector(t(mf_preds_all[drug_set,target_set]))
      prediction_adj_mats_1 = make_prediction_adj_mat(train_mat = train_mat[drug_set,target_set], drug_adj_mat = drug_adj_mat, target_adj_mat = target_adj_mat_below, drug_set = drug_set, target_set = target_set)
      prediction_adj_mats_2 = make_prediction_adj_mat(train_mat = train_mat[drug_set,target_set], drug_adj_mat = drug_adj_mat, target_adj_mat = target_adj_mat_above, drug_set = drug_set, target_set = target_set)
      
      adj_mats = list()
      adj_mats[[1]] = prediction_adj_mats_1[[1]]
      adj_mats[[2]] = prediction_adj_mats_1[[2]]
      adj_mats[[3]] = prediction_adj_mats_2[[2]]
      
      train_data = train_mat[drug_set,target_set]
      train_data = as.vector(t(train_data))
      
      prediction = crf_predict(train_data = train_data, X = X, alpha = params[[1]], adj_mats = adj_mats, beta = params[[2]])
      
      prediction_mat = matrix(prediction, ncol = length(target_set), byrow = T)
      
      inds = which(dataset_triplets[test_ind,1]%in%drug_set & dataset_triplets[test_ind,2]%in%target_set)
      for (i in 1:length(inds)){
        temp_d = match(dataset_triplets[test_ind[inds[i]],1], drug_set)
        temp_t = match(dataset_triplets[test_ind[inds[i]],2], target_set)
        prediction = prediction_mat[temp_d, temp_t]
        crf_predictions[test_ind[inds[i]]] = prediction
      }
      
      inds = which(!is.na(crf_predictions))
      crf_metrics = get_metrics(crf_predictions[inds], dataset_triplets[inds,3], 12.1)
      
      cat('fold ',fold,'\n')
      cat('auc so far, crf: ',round(crf_metrics[[2]], digits = 5),'\n')
      cat('aupr so far, crf: ',round(crf_metrics[[3]], digits = 5),'\n\n')
      
      }
  }
  
}

