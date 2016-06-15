## the dual crf with targets 147, 183, 415, 416, 479, 484

load('../data/mix_dataset.Rda')
target_sim = read.table('../data/mix_dataset_target_target_sim.txt')
target_sim = as.matrix(target_sim)
drug_sim = read.table('../data/mix_dataset_drug_drug_sim.txt')
drug_sim = as.matrix(drug_sim)

n_drugs = length(unique(mix_dataset[,1]))
n_targets = length(unique(mix_dataset[,2]))

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

target_adj_mat = make_adjacency_mat_targets(target_sim, 0.4)
drug_adj_mat = make_adjacency_mat(drug_sim)

target_set = c(32,33,203,490)
target_set = c(147,183,415,416,479,484)
#target_set = c(12, 180, 390 )
#target_set = c(11, 143, 287)
#target_set = c(8, 139, 516)
target_set = c(43, 28, 329, 406)
target_set = c(65,171)
target_set = c(134,137, 298)
target_set = c(253,472,473)
target_set = c(26,648,649,650)
target_set = c(113,325,336)
target_set = c(52,507)
target_set = c(19,256)
target_set = c(246,247)
target_set = c(21, 457)
target_set = c(11,14,143,170,287)
target_set = cluster_list[[1]]
target_set = cluster_list[[2]]
adj_mat_cols = target_adj_mat[target_set, target_set]
adj_mat_rows = drug_adj_mat
dt_mat_temp = dt_mat[,target_set]
mf_preds_train_mat_temp = mf_preds_train_mat[,target_set]
mf_preds_all_mat_temp = mf_preds_all[,target_set]
params = train_crf_dual_simple(training_data_mat = dt_mat_temp, row_adj_mat = adj_mat_rows, col_adj_mat = adj_mat_cols, mf_preds_train_mat = mf_preds_train_mat_temp, eta = 0.01, crf_iters = 600)
crf_preds_temp = simple_crf_predict(training_data_mat = dt_mat_temp, row_adj_mat = adj_mat_rows, col_adj_mat = adj_mat_cols, alpha = params[[1]], row_beta = params[[2]], col_beta = params[[3]], X_mat = mf_preds_all_mat_temp)

print_prediction_metrics_set(target_set = target_set, mix_dataset = mix_dataset, test_ind = test_ind, crf_preds = crf_preds_temp, mf_preds_all_mat = mf_preds_all_mat_temp)

train_and_predict_single(mf_preds_train_mat = mf_preds_train_mat, dt_mat = dt_mat, sim_mat = drug_sim, mix_dataset = mix_dataset, mf_preds_all = mf_preds_all, target_set = target_set, test_ind = test_ind, adj_mat = drug_adj_mat, crf_iters = 800)

crf_predictions = rep(NA, nrow(mix_dataset))
mf_predictions = rep(NA, nrow(mix_dataset))
for (i in 1:141){
  target_set = cluster_list[[i]]
  adj_mat_cols = target_adj_mat[target_set, target_set]
  adj_mat_rows = drug_adj_mat
  dt_mat_temp = dt_mat[,target_set]
  mf_preds_train_mat_temp = mf_preds_train_mat[,target_set]
  mf_preds_all_mat_temp = mf_preds_all[,target_set]
  params = train_crf_dual_simple(training_data_mat = dt_mat_temp, row_adj_mat = adj_mat_rows, col_adj_mat = adj_mat_cols, mf_preds_train_mat = mf_preds_train_mat_temp, eta = 0.01, crf_iters = 600)
  crf_preds_temp = simple_crf_predict(training_data_mat = dt_mat_temp, row_adj_mat = adj_mat_rows, col_adj_mat = adj_mat_cols, alpha = params[[1]], row_beta = params[[2]], col_beta = params[[3]], X_mat = mf_preds_all_mat_temp)
  for (t in target_set){
    inds = which(mix_dataset[test_ind,2]==t)
    drugs = mix_dataset[test_ind[inds],1]
    labels = mix_dataset[test_ind[inds],3]
    crf_predictions_test_col = crf_preds_temp[drugs, match(t,target_set)]
    mf_predictions_test_col =  mf_preds_all_mat_temp[drugs, match(t,target_set)]
    crf_predictions[test_ind[inds]] = crf_predictions_test_col
    mf_predictions[test_ind[inds]] = mf_predictions_test_col
  }
  
  inds = which(!is.na(mf_predictions))
  mf_metrics = get_metrics(mf_predictions[inds], mix_dataset[inds,3], 7)
  crf_metrics = get_metrics(crf_predictions[inds], mix_dataset[inds,3], 7)
  
  cat('all rmse (mf, crf) so far: ',round(mf_metrics[[1]], digits = 5),', ',round(crf_metrics[[1]], digits = 5),'\n')
  cat('all test auc (mf, crf) so far: ',round(mf_metrics[[2]], digits = 5),', ',round(crf_metrics[[2]], digits = 5),'\n')
  cat('all test aupr (mf, crf) so far: ',round(mf_metrics[[3]], digits = 5),', ',round(crf_metrics[[3]], digits = 5),'\n\n')
  
}

remaining_preds = train_and_predict_single(mf_preds_train_mat = mf_preds_train_mat, dt_mat = dt_mat, sim_mat = drug_sim, mix_dataset = mix_dataset, mf_preds_all = mf_preds_all, target_set = target_set, test_ind = test_ind, adj_mat = drug_adj_mat, crf_iters = 800)

all_crf_preds = rep(NA, nrow(mix_dataset))
all_crf_preds[which(!is.na(crf_predictions))] = crf_predictions[which(!is.na(crf_predictions))]
all_crf_preds[which(!is.na(remaining_preds[[1]]))] = remaining_preds[[1]][which(!is.na(remaining_preds[[1]]))]

all_mf_preds = rep(NA, nrow(mix_dataset))
all_mf_preds[which(!is.na(mf_predictions))] = mf_predictions[which(!is.na(mf_predictions))]
all_mf_preds[which(!is.na(remaining_preds[[2]]))] = remaining_preds[[2]][which(!is.na(remaining_preds[[2]]))]

inds = which(!is.na(all_crf_preds))
crf_metrics = get_metrics(all_crf_preds[inds], mix_dataset[inds,3], 7)
mf_metrics = get_metrics(all_mf_preds[inds], mix_dataset[inds,3], 7)

cat('all rmse (mf, crf) so far: ',round(mf_metrics[[1]], digits = 5),', ',round(crf_metrics[[1]], digits = 5),'\n')
cat('all test auc (mf, crf) so far: ',round(mf_metrics[[2]], digits = 5),', ',round(crf_metrics[[2]], digits = 5),'\n')
cat('all test aupr (mf, crf) so far: ',round(mf_metrics[[3]], digits = 5),', ',round(crf_metrics[[3]], digits = 5),'\n\n')

crf_predictions = rep(NA, nrow(mix_dataset))
mf_predictions = rep(NA, nrow(mix_dataset))
for (i in 1:50){
  target_set = cluster_list[[i]]
  target_set = which(!c(1:698%in%target_set_1))
  adj_mat_cols = target_adj_mat[target_set, target_set]
  adj_mat_rows = drug_adj_mat
  dt_mat_temp = dt_mat[,target_set]
  mf_preds_train_mat_temp = mf_preds_train_mat[,target_set]
  mf_preds_all_mat_temp = mf_preds_all[,target_set]
  params = train_crf_dual_simple(training_data_mat = dt_mat_temp, row_adj_mat = adj_mat_rows, col_adj_mat = adj_mat_cols, mf_preds_train_mat = mf_preds_train_mat_temp, eta = 0.01, crf_iters = 600)
  crf_preds_temp = simple_crf_predict(training_data_mat = dt_mat_temp, row_adj_mat = adj_mat_rows, col_adj_mat = adj_mat_cols, alpha = params[[1]], row_beta = params[[2]], col_beta = params[[3]], X_mat = mf_preds_all_mat_temp)
  for (t in target_set){
    inds = which(mix_dataset[test_ind,2]==t)
    drugs = mix_dataset[test_ind[inds],1]
    labels = mix_dataset[test_ind[inds],3]
    crf_predictions_test_col = crf_preds_temp[drugs, match(t,target_set)]
    mf_predictions_test_col =  mf_preds_all_mat_temp[drugs, match(t,target_set)]
    crf_predictions[test_ind[inds]] = crf_predictions_test_col
    mf_predictions[test_ind[inds]] = mf_predictions_test_col
  }
  
  inds = which(!is.na(mf_predictions))
  mf_metrics = get_metrics(mf_predictions[inds], mix_dataset[inds,3], 7)
  crf_metrics = get_metrics(crf_predictions[inds], mix_dataset[inds,3], 7)
  
  cat('all rmse (mf, crf) so far: ',round(mf_metrics[[1]], digits = 5),', ',round(crf_metrics[[1]], digits = 5),'\n')
  cat('all test auc (mf, crf) so far: ',round(mf_metrics[[2]], digits = 5),', ',round(crf_metrics[[2]], digits = 5),'\n')
  cat('all test aupr (mf, crf) so far: ',round(mf_metrics[[3]], digits = 5),', ',round(crf_metrics[[3]], digits = 5),'\n\n')
  
}


crf_predictions = rep(NA, nrow(mix_dataset))
mf_predictions = rep(NA, nrow(mix_dataset))
