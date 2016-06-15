load("experiment_3.rda")
length(test_folds)

load('kiba_target_sim.rda')

adj_mat = make_adjacency_mat_targets(sim_mat, 0.5)

crf_predictions = rep(NA, nrow(dt_triplet))
mf_predictions = rep(NA, nrow(dt_triplet))


#for (i in 1:5){
i = 1
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
  
  for (d in 1:n_drugs){
    
    cat('fold ',i,', drug ',d,'... ',length(which(dt_mat[d,]>=0)),' observations\n')
    
    mf_pred_train_row_d = mf_preds_train_mat[d,which(!is.na(mf_preds_train_mat[d,]))]
    adj_mat_train_row_d = make_training_adj_mat_for_row(dt_mat, sim_mat, d, 0.5)
    if (length(which(adj_mat_train_row_d>0))>0){
      
      training_vals_row_d = dt_mat[d,which(dt_mat[d,]>=0)]
      
      eta = 0.01
      
      params = train_crf_row(y = training_vals_row_d, X = mf_pred_train_row_d, adj_mat = adj_mat_train_row_d, crf_iters = 2000, eta = eta)  
      cat('learned parameters: ', params[[1]], params[[2]],'\n')
      
      inds = which(dt_triplet[test_ind,1] == d)
      labels_test_row = dt_triplet[test_ind[inds], 3]
      mf_prediction_row = mf_preds_all[d,]
      mf_prediction_test_row = mf_preds_all[cbind(dt_triplet[test_ind[inds],1],dt_triplet[test_ind[inds],2])]
      mf_predictions[test_ind[inds]] = mf_prediction_test_row
      
      cat('making crf predictions..\n')
      
      crf_prediction_row = make_crf_predictions_row(params[[1]], params[[2]], column = dt_mat[d,], adj_mat = adj_mat, X = mf_prediction_row)
      crf_prediction_test_row = crf_prediction_row[dt_triplet[test_ind[inds],2]]
      crf_predictions[test_ind[inds]] = crf_prediction_test_row
      
    } else {
      
      inds = which(dt_triplet[test_ind,1] == d)
      labels_test_row = dt_triplet[test_ind[inds], 3]
      mf_prediction_row = mf_preds_all[d,]
      mf_prediction_test_row = mf_preds_all[cbind(dt_triplet[test_ind[inds],1],dt_triplet[test_ind[inds],2])]
      mf_predictions[test_ind[inds]] = mf_prediction_test_row
      
      crf_prediction_test_row = mf_prediction_test_row
      crf_predictions[test_ind[inds]] = mf_prediction_test_row
      
    }
    
    
    mf_metrics = get_metrics(mf_prediction_test_row, labels_test_row, 12.1)
    crf_metrics = get_metrics(crf_prediction_test_row, labels_test_row, 12.1)
    
    cat('drug rmse (mf, crf): ',mf_metrics[[1]],', ',crf_metrics[[1]],'\n')
    cat('drug auc (mf, crf): ',mf_metrics[[2]],', ',crf_metrics[[2]],'\n')
    cat('drug aupr (mf, crf): ',mf_metrics[[3]],', ',crf_metrics[[3]],'\n')
    
    inds = which(!is.na(mf_predictions))
    mf_metrics = get_metrics(mf_predictions[inds], dt_triplet[inds,3], 12.1)
    crf_metrics = get_metrics(crf_predictions[inds], dt_triplet[inds,3], 12.1)
    
    cat('all test rmse (mf, crf) so far: ',round(mf_metrics[[1]], digits = 2),', ',round(crf_metrics[[1]], digits = 2),'\n')
    cat('all test auc (mf, crf) so far: ',round(mf_metrics[[2]], digits = 2),', ',round(crf_metrics[[2]], digits = 2),'\n')
    cat('all test aupr (mf, crf) so far: ',round(mf_metrics[[3]], digits = 2),', ',round(crf_metrics[[3]], digits = 2),'\n\n')
    
  }
#}

1456
