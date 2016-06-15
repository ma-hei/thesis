
load('../data/mix_dataset.Rda')
target_sim = read.table('../data/mix_dataset_target_target_sim.txt')
target_sim = as.matrix(target_sim)
drug_sim = read.table('../data/mix_dataset_drug_drug_sim.txt')
drug_sim = as.matrix(drug_sim)

n_drugs = length(unique(mix_dataset[,1]))
n_targets = length(unique(mix_dataset[,2]))

n_folds = 5
test_folds = get_folds(mix_dataset, n_folds)

crf_predictions = rep(NA, nrow(mix_dataset))
crf_predictions_single = rep(NA, nrow(mix_dataset))
mf_predictions = rep(NA, nrow(mix_dataset))
for (i in 1:n_folds){
  
  test_ind = test_folds[[i]]
  train_ind = setdiff(1:nrow(mix_dataset),test_ind)
  
  ## make the training data into matrix format
  dt_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
  dt_mat[mix_dataset[train_ind,c(1,2)]] = mix_dataset[train_ind,3]
  
  cat('getting MF cv-prediction on training data..\n')
  mf_preds_train = get_mf_cv(dt_mat, 400)
  mf_preds_train_mat = mf_preds_train[[2]]
  
  cat('getting MF predictions for complete matrix..\n')
  mf_preds_all = get_libmf_prediction(dt_mat, 400)
  
  target_adj_mat = make_adjacency_mat_targets(target_sim, 0.5)
  drug_adj_mat = make_adjacency_mat(drug_sim, 0.9)
  drug_adj_mat2 = make_adjacency_mat_4n(drug_sim)
  
  ans = get_target_clusters(target_sim, 0.5)
  cluster_list = ans[[1]]
  remaining_targets = ans[[2]]
  
  ## first go over all clusters and make predictions with the crf that
  ## integrates both similarities
  
   for (k in 1:length(cluster_list)){
     
    cat('\ntarget cluster ',k,'/',length(cluster_list),'.. targets of cluster: ',cluster_list[[k]],'\n')
    n_train = length(which(dt_mat[,cluster_list[[k]]]>=0))
    cat('number of observations in cluster:', n_train,'\n')
    ## take only the targets in the cluster..
    target_set = cluster_list[[k]]
    ## build the adjacencey mat for that subset..
    adj_mat_cols = target_adj_mat[target_set, target_set]
    adj_mat_rows = drug_adj_mat
    ## take only the targets in the cluster from the training data
    dt_mat_temp = dt_mat[,target_set]
    ## and the mf predictions for the targets in the cluster..
    mf_preds_train_mat_temp = mf_preds_train_mat[,target_set]
    mf_preds_all_mat_temp = mf_preds_all[,target_set]
    
    ## train the crf and make predictions
    if (n_train>1000){
      eta = 0.001
    } else{
      eta = 0.01
    }
    cat('learning params for crf..\n')
    #params = train_crf_dual_simple(training_data_mat = dt_mat_temp, row_adj_mat = adj_mat_rows, col_adj_mat = adj_mat_cols, mf_preds_train_mat = mf_preds_train_mat_temp, eta = eta, crf_iters = 600)
    #cat('making predictions for cluster\n')
    #crf_preds_temp = simple_crf_predict(training_data_mat = dt_mat_temp, row_adj_mat = adj_mat_rows, col_adj_mat = adj_mat_cols, alpha = params[[1]], row_beta = params[[2]], col_beta = params[[3]], X_mat = mf_preds_all_mat_temp)
    
    ## put the predicted values into the prediction vector
    ## now go over each target in the set, find the drugs which need to be predicted 
    ## for the target and put the prediction into the prediction vector
#     for (t in target_set){
#       inds = which(mix_dataset[test_ind,2]==t)
#       drugs = mix_dataset[test_ind[inds],1]
#       labels = mix_dataset[test_ind[inds],3]
#       crf_predictions_test_col = crf_preds_temp[drugs, match(t,target_set)]
#       crf_predictions[test_ind[inds]] = crf_predictions_test_col
#       mf_predictions_test_col =  mf_preds_all_mat_temp[drugs, match(t,target_set)]
#       mf_predictions[test_ind[inds]] = mf_predictions_test_col
#     }
    
    temp_inds = mix_dataset[,c(1,2)]
    preds_single = train_and_predict_single(mf_preds_train_mat = mf_preds_train_mat, dt_mat = dt_mat, sim_mat = drug_sim, temp_inds = temp_inds, mf_preds_all = mf_preds_all, target_set = target_set, test_ind = test_ind, adj_mat = drug_adj_mat2, crf_iters = 600)
    crf_predictions_single[which(!is.na(preds_single))] = preds_single[which(!is.na(preds_single))]
    
    inds = which(!is.na(crf_predictions_single))
    crf_metrics_single = get_metrics(crf_predictions_single[inds], mix_dataset[inds,3], 7)
    
    inds = which(!is.na(crf_predictions))
    crf_metrics = get_metrics(crf_predictions[inds], mix_dataset[inds,3], 7)
    mf_metrics = get_metrics(mf_predictions[inds], mix_dataset[inds,3], 7)
    
    cat('fold ',i,'\n')
    cat('rmse so far, crf: ',round(crf_metrics[[1]], digits = 5),' crf single: ',crf_metrics_single[[1]],' mf: ', round(mf_metrics[[1]], digits = 5),'\n')
    cat('auc so far, crf: ',round(crf_metrics[[2]], digits = 5),' crf single :',crf_metrics_single[[2]],' mf: ',round(mf_metrics[[2]], digits = 5),'\n')
    cat('aupr so far, crf: ',round(crf_metrics[[3]], digits = 5),' crf single: ',crf_metrics_single[[3]],' mf: ',round(mf_metrics[[3]], digits = 5),'\n')
  
#      cat('fold ',i,'\n')
#      cat('rmse so far, crf: ',round(crf_metrics[[1]], digits = 5),' mf: ', round(mf_metrics[[1]], digits = 5),'\n')
#      cat('auc so far, crf: ',round(crf_metrics[[2]], digits = 5),' mf: ',round(mf_metrics[[2]], digits = 5),'\n')
#      cat('aupr so far, crf: ',round(crf_metrics[[3]], digits = 5),' mf: ',round(mf_metrics[[3]], digits = 5),'\n')


  }
  
  
  #remaining_targets = c(1:n_targets)
  ## and then make predictions for the ramining targets, using only the drug similarity
  ## this function needs the indices of the triplets..
#   temp_inds = mix_dataset[,c(1,2)]
#   remaining_preds = train_and_predict_single(mf_preds_train_mat = mf_preds_train_mat, dt_mat = dt_mat, sim_mat = drug_sim, temp_inds = temp_inds, mf_preds_all = mf_preds_all, target_set = remaining_targets, test_ind = test_ind, adj_mat = drug_adj_mat, crf_iters = 600)
#   crf_predictions[which(!is.na(remaining_preds))] = remaining_preds[which(!is.na(remaining_preds))]
#   
#   inds = which(!is.na(crf_predictions))
#   crf_metrics = get_metrics(crf_predictions[inds], mix_dataset[inds,3], 7)
#   
#   cat('fold ',i,'\n')
#   cat('rmse: ',round(crf_metrics[[1]], digits = 5),'\n')
#   cat('auc: ',round(crf_metrics[[2]], digits = 5),'\n')
#   cat('aupr: ',round(crf_metrics[[3]], digits = 5),'\n\n')
  
}