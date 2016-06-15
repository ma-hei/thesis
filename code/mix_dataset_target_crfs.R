load('../data/mix_dataset.Rda')
drug_sim = read.table('../data/mix_dataset_drug_drug_sim.txt')
drug_sim = as.matrix(drug_sim)
n_drugs = length(unique(mix_dataset[,1]))
n_targets = length(unique(mix_dataset[,2]))

## create an adjacencey mat where all drugs with similarity > 0.9 are connected
## and each drug for which no such drug exists is connected to its nearest neighbor
drug_adj_mat = make_adj_mat_thresh(sim_mat = drug_sim, thresh = 0.9, TRUE)

n_folds = 5
test_folds = get_folds(mix_dataset, n_folds)

crf_predictions = rep(NA, nrow(mix_dataset))

for (i in 1:n_folds){
  
  test_ind = test_folds[[i]]
  train_ind = setdiff(1:nrow(mix_dataset),test_ind)
  
  ## make the training data into matrix format
  ## missing values are -1, everything else is >0
  train_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
  train_mat[mix_dataset[train_ind,c(1,2)]] = mix_dataset[train_ind,3]
  
  cat('getting MF cv-prediction on training data..\n')
  mf_preds_train = get_mf_cv(dt_mat, 400)
  mf_preds_train_mat = mf_preds_train[[2]]
  
  cat('getting MF predictions for complete matrix..\n')
  mf_preds_all = get_libmf_prediction(dt_mat, 400)
  
  for (t in 1:n_targets){
    
    cat('target ',t,'... ',length(which(train_mat[,t]>=0)),' observations\n')
    ## it only makes sense to train a crf, when more than 1 observation are given for the target
    if (length(which(train_mat[,t]>=0))>1){
      
      ## get the training data for the target
      training_data = train_mat[which(train_mat[,t]>0), t]
      inds = which(train_mat[,t]>=0)
      ## make an adjacency mat for the training data in sparse matrix format
      training_adj_mat = drug_adj_mat[inds,inds]
      inds = which(training_adj_mat>0, arr.ind = T)
      training_adj_mat = sparseMatrix(i = inds[,1], j = inds[,2], x = training_adj_mat[inds], dims = dim(training_adj_mat))
      
      ## get the prediction of MF on the training data
      X = mf_preds_train_mat[which(!is.na(mf_preds_train_mat[,t])),t]
      
      ## train the crf with the true training data, the adjacency mat of the training data and the MF prediction
      adj_mats = list()
      adj_mats[[1]] = training_adj_mat
      params = train_crf(training_data = training_data, adj_mats = adj_mats, X = X, crf_iters = 600, eta = 0.004)
      
      ## now make a prediction for the target
      adj_mats = list()
      adj_mats[[1]] = drug_adj_mat
      prediction = crf_predict(train_data = train_mat[,t], X = mf_preds_all[,t], alpha = params[[1]], adj_mats = adj_mats, beta = params[[2]])
      
      ## put the prediction into the prediction vector
      inds = which(mix_dataset[test_ind,2] == t)
      crf_predictions[test_ind[inds]] = prediction[mix_dataset[test_ind[inds],1]]
      
    } else {
      ## when only one observation is given for the target, just take the MF prediction
      inds = which(mix_dataset[test_ind,2] == t)
      crf_predictions[test_ind[inds]] = mf_preds_all[mix_dataset[test_ind[inds],1],t]
    }
    
    inds = which(!is.na(crf_predictions))
    crf_metrics = get_metrics(crf_predictions[inds], mix_dataset[inds,3], 7.0)
    
    cat('auc so far, crf: ',round(crf_metrics[[2]], digits = 5),'\n')
    cat('aupr so far, crf: ',round(crf_metrics[[3]], digits = 5),'\n\n')
    
  }
  
}

