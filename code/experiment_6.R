# prepare the 5 folds

source('pipeline.R')
source('recosystem.R')

load('kiba_data.rda')

drug_sim = simi_mat
load('kiba_target_sim.rda')
target_sim = sim_mat

n_drugs = length(unique(dt_triplet[,1]))
n_targets = length(unique(dt_triplet[,2]))

hist(drug_sim)
hist(target_sim)

dt_triplet[,3] = dt_triplet[,3]*-1
dt_triplet[,3] = dt_triplet[,3] - min(dt_triplet[,3]) 
hist(dt_triplet[,3])

n_folds = 5
test_folds = get_folds(dt_triplet, n_folds)

crf_predictions = rep(NA, nrow(dt_triplet))
mf_predictions = rep(NA, nrow(dt_triplet))

drug_adj_mat = make_adjacency_mat(drug_sim)

mf_preds_train_mat_folds = list()
mf_preds_folds_all = list()


for (i in 1:n_folds){
  
  cat('fold ',i,'\n')
  
  test_ind = test_folds[[i]]
  train_ind = setdiff(1:nrow(dt_triplet),test_ind)
  
  dt_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
  dt_mat[dt_triplet[train_ind,c(1,2)]] = dt_triplet[train_ind,3]
  
  cat('getting MF cv-prediction on training data..\n')
  mf_preds_train = get_mf_cv(dt_mat, 400)
  mf_preds_train_mat = mf_preds_train[[2]]
  
  mf_preds_train_mat_folds[[i]] = mf_preds_train_mat
  
  cat('getting MF predictions for complete matrix..\n')
  mf_preds_all = get_libmf_prediction(dt_mat, 400)
  
  mf_preds_folds_all[[i]] = mf_preds_all
  
}

save(test_folds, drug_adj_mat, drug_sim, dt_triplet, mf_preds_train_mat_folds, mf_preds_folds_all, file="kiba_drugcrfs_data.rda")
