load('../data/mix_dataset.Rda')

n_drugs = length(unique(mix_dataset[,1]))
n_targets = length(unique(mix_dataset[,2]))

drug_sim = read.table('../data/mix_dataset_drug_drug_sim.txt')
drug_sim = as.matrix(drug_sim)
target_sim = read.table('../data/mix_dataset_target_target_sim.txt')

target_adj_mat = make_adjacency_mat_targets(target_sim, 0.5)
drug_adj_mat = make_adjacency_mat(drug_sim)

test_folds = get_folds(mix_dataset, 5)

crf_predictions = rep(NA, nrow(mix_dataset))
mf_predictions = rep(NA, nrow(mix_dataset))

i = 1
test_ind = test_folds[[i]]
train_ind = setdiff(1:nrow(mix_dataset),test_ind)

dt_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
dt_mat[mix_dataset[train_ind,c(1,2)]] = mix_dataset[train_ind,3]

cat('getting MF cv-prediction on training data..\n')
mf_preds_train = get_mf_cv(dt_mat, 400)
mf_preds_train_mat = mf_preds_train[[2]]

cat('getting MF predictions for complete matrix..\n')
mf_preds_all = get_libmf_prediction(dt_mat, 400)

sim_mat = drug_sim
adj_mat = drug_adj_mat

for (t in 1:n_targets){
  
#for (t in c(65,171)){
  
  cat('fold ',i,', target ',t,'... ',length(which(dt_mat[,t]>=0)),' observations\n')
  
  mf_pred_train_col_t = mf_preds_train_mat[which(!is.na(mf_preds_train_mat[,t])),t]
  adj_mat_train_col_t = make_training_adj_mat_for_column(dt_mat, sim_mat, t)
  training_vals_col_t = dt_mat[which(dt_mat[,t]>=0),t]
  
  
  cat(length(which(adj_mat_train_col_t>0)),'\n')
  if (length(which(dt_mat[,t]>=0))>500){
    eta = 0.001
  } else{
    eta = 0.01
  }
  
  params = train_crf_row(y = training_vals_col_t, X = mf_pred_train_col_t, adj_mat = adj_mat_train_col_t, crf_iters = 1000, eta = eta)  
  cat('learned parameters: ', params[[1]], params[[2]],'\n')
  
  inds = which(mix_dataset[test_ind,2] == t)
  labels_test_col = mix_dataset[test_ind[inds], 3]
  mf_prediction_col = mf_preds_all[,t]
  mf_prediction_test_col = mf_preds_all[cbind(mix_dataset[test_ind[inds],1],mix_dataset[test_ind[inds],2])]
  mf_predictions[test_ind[inds]] = mf_prediction_test_col
  
  cat('making crf predictions..\n')
  
  crf_prediction_col = make_crf_predictions_row(params[[1]], params[[2]], column = dt_mat[,t], adj_mat = adj_mat, X = mf_prediction_col)
  crf_prediction_test_col = crf_prediction_col[mix_dataset[test_ind[inds],1]]
  crf_predictions[test_ind[inds]] = crf_prediction_test_col
  
  mf_metrics = get_metrics(mf_prediction_test_col, labels_test_col, 7)
  crf_metrics = get_metrics(crf_prediction_test_col, labels_test_col, 7)
  
  cat('target rmse (mf, crf): ',mf_metrics[[1]],', ',crf_metrics[[1]],'\n')
  cat('target auc (mf, crf): ',mf_metrics[[2]],', ',crf_metrics[[2]],'\n')
  cat('target aupr (mf, crf): ',mf_metrics[[3]],', ',crf_metrics[[3]],'\n')
  
  inds = which(!is.na(mf_predictions))
  mf_metrics = get_metrics(mf_predictions[inds], mix_dataset[inds,3], 7)
  crf_metrics = get_metrics(crf_predictions[inds], mix_dataset[inds,3], 7)
  
  cat('all test rmse (mf, crf) so far: ',round(mf_metrics[[1]], digits = 3),', ',round(crf_metrics[[1]], digits = 3),'\n')
  cat('all test auc (mf, crf) so far: ',round(mf_metrics[[2]], digits = 3),', ',round(crf_metrics[[2]], digits = 3),'\n')
  cat('all test aupr (mf, crf) so far: ',round(mf_metrics[[3]], digits = 3),', ',round(crf_metrics[[3]], digits = 3),'\n\n')
  
}