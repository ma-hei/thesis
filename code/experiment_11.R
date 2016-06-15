dataset = read.table('../data/known_drug-target_interaction_affinities_pKi__Metz_et_al.2011.txt')
dataset = as.matrix(dataset)

drug_sim = read.table('../data/drug-drug_similarities_2D__Metz_et_al.2011.txt')
drug_sim = as.matrix(drug_sim)
drug_sim = drug_sim/100

target_sim = read.table('../data/target-target_similarities_WS_normalized__Metz_et_al.2011.txt')
target_sim = as.matrix(target_sim)
target_sim = target_sim/100
hist(target_sim)


inds = which(!is.na(dataset), arr.ind = T)
dataset_triplets = matrix(0, nrow = nrow(inds), ncol = 3)
dataset_triplets[,c(1,2)] = inds
dataset_triplets[,3] = dataset[inds]

n_drugs = length(unique(dataset_triplets[,1]))
n_targets = length(unique(dataset_triplets[,2]))

n_folds = 5
test_folds = get_folds(dataset_triplets, n_folds)

mf_predictions = rep(NA, nrow(dataset_triplets))
target_crf_predictions = rep(NA, nrow(dataset_triplets))

drug_adj_mat = make_adjacency_mat_4n(drug_sim)
for (i in 1:n_folds){
  
  test_ind = test_folds[[i]]
  train_ind = setdiff(1:nrow(dataset_triplets),test_ind)
  
  dt_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
  dt_mat[dataset_triplets[train_ind,c(1,2)]] = dataset_triplets[train_ind,3]
  
  cat('getting MF cv-prediction on training data..\n')
  mf_preds_train = get_mf_cv(dt_mat, 400)
  mf_preds_train_mat = mf_preds_train[[2]]
  
  cat('getting MF predictions for complete matrix..\n')
  mf_preds_all = get_libmf_prediction(dt_mat, 400)
  
  temp_inds = dataset_triplets[,c(1,2)]
  preds_single = train_and_predict_single(mf_preds_train_mat = mf_preds_train_mat, dt_mat = dt_mat, sim_mat = drug_sim, temp_inds = temp_inds, mf_preds_all = mf_preds_all, target_set = c(1:n_targets), test_ind = test_ind, adj_mat = drug_adj_mat, crf_iters = 400, dataset_triplets)
  target_crf_predictions[which(!is.na(preds_single))] = preds_single[which(!is.na(preds_single))]
  
  inds = which(!is.na(target_crf_predictions))
  crf_metrics = get_metrics(target_crf_predictions[inds], dataset_triplets[inds,3], 7.6)
  
  cat('fold ',i,'\n')
  cat('auc so far, crf: ',round(crf_metrics[[2]], digits = 5),'\n')
  cat('aupr so far, crf: ',round(crf_metrics[[3]], digits = 5),'\n')
  
  
}

mf_metrics = get_metrics(mf_predictions, dataset_triplets[,3], 7.6)

