load('kiba_data.rda')

n_drugs = length(unique(dt_triplet[,1]))
n_targets = length(unique(dt_triplet[,2]))

dim(simi_mat)

dt_triplet[,3] = dt_triplet[,3]*-1
dt_triplet[,3] = dt_triplet[,3] - min(dt_triplet[,3]) 
hist(dt_triplet[,3])

test_folds = get_folds(dt_triplet, 5)

crf_predictions = rep(NA, nrow(dt_triplet))
mf_predictions = rep(NA, nrow(dt_triplet))

adj_mat = make_adjacency_mat(simi_mat)

for (i in 1:n_folds){
  cat('fold ',i,'\n')
  
  test_ind = test_folds[[i]]
  train_ind = setdiff(1:nrow(dt_triplet),test_ind)
  
  dt_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
  dt_mat[dt_triplet[train_ind,c(1,2)]] = dt_triplet[train_ind,3]
  
  cat('getting MF cv-prediction on training data..\n')
  mf_preds_train = get_mf_cv(dt_mat, 400)
  mf_preds_train_mat = mf_preds_train[[2]]
  
  cat('getting MF predictions for complete matrix..\n')
  mf_preds_all = get_libmf_prediction(dt_mat, 400)
  
  for (t in 1:n_targets){
    
    cat('fold ',i,', target ',t,'... ',length(which(dt_mat[,t]>=0)),' observations\n')
    
    mf_pred_train_col_t = mf_preds_train_mat[which(!is.na(mf_preds_train_mat[,t])),t]
    adj_mat_train_col_t = make_training_adj_mat_for_column(dt_mat, simi_mat, t)
    training_vals_col_t = dt_mat[which(dt_mat[,t]>=0),t]
    
    if (length(which(dt_mat[,t]>=0))>500){
      eta = 0.001
    } else{
      eta = 0.01
    }
    
    params = train_crf_row(y = training_vals_col_t, X = mf_pred_train_col_t, adj_mat = adj_mat_train_col_t, crf_iters = 200, eta = eta)  
    cat('learned parameters: ', params[[1]], params[[2]],'\n')
    
    inds = which(dt_triplet[test_ind,2] == t)
    labels_test_col = dt_triplet[test_ind[inds], 3]
    mf_prediction_col = mf_preds_all[,t]
    mf_prediction_test_col = mf_preds_all[cbind(dt_triplet[test_ind[inds],1],dt_triplet[test_ind[inds],2])]
    mf_predictions[test_ind[inds]] = mf_prediction_test_col
    
    cat('making crf predictions..\n')
    
    crf_prediction_col = make_crf_predictions_row(params[[1]], params[[2]], column = dt_mat[,t], adj_mat = adj_mat, X = mf_prediction_col)
    crf_prediction_test_col = crf_prediction_col[dt_triplet[test_ind[inds],1]]
    crf_predictions[test_ind[inds]] = crf_prediction_test_col
    
    mf_metrics = get_metrics(mf_prediction_test_col, labels_test_col, 12.1)
    crf_metrics = get_metrics(crf_prediction_test_col, labels_test_col, 12.1)
 
    cat('target rmse (mf, crf): ',mf_metrics[[1]],', ',crf_metrics[[1]],'\n')
    cat('target auc (mf, crf): ',mf_metrics[[2]],', ',crf_metrics[[2]],'\n')
    cat('target aupr (mf, crf): ',mf_metrics[[3]],', ',crf_metrics[[3]],'\n')
    
    inds = which(!is.na(mf_predictions))
    mf_metrics = get_metrics(mf_predictions[inds], dt_triplet[inds,3], 12.1)
    crf_metrics = get_metrics(crf_predictions[inds], dt_triplet[inds,3], 12.1)
    
    cat('all test rmse (mf, crf) so far: ',mf_metrics[[1]],', ',crf_metrics[[1]],'\n')
    cat('all test auc (mf, crf) so far: ',mf_metrics[[2]],', ',crf_metrics[[2]],'\n')
    cat('all test aupr (mf, crf) so far: ',mf_metrics[[3]],', ',crf_metrics[[3]],'\n\n')
    
  }
}

save(crf_predictions, file="crf_predictions_eta01.rda")