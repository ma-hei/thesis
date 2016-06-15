get_folds = function(triplets, n_folds){
  
  n = nrow(triplets)
  folds = vector(n_folds, mode='list')
  shuf = sample(n)
  m = n %/% n_folds
  for (i in 1:(n_folds-1)) {
    folds[[i]] = shuf[1:m]
    shuf = setdiff(shuf, folds[[i]])
  }
  folds[[n_folds]] = shuf
  
  return(folds)
  
}

make_Bk_mats = function(adj_mats){
  
  n_train = nrow(similarity_mats[[1]])
  
  Bk_mats = list()
  for (i in 1:length(similarity_mats)){
    Bk_mats[[i]] = matrix(0, nrow = n_train, ncol = n_train)
    Bk_mats[[i]] = similarity_mats[[i]]*(-1)
    for (k in 1:n_train){
      Bk_mats[[i]][k,k] = sum(similarity_mats[[i]][k,])-similarity_mats[[i]][k,k]
    }
  }
  
  return(Bk_mats)
  
}

make_triplets = function(dataset_mat){
  
  inds = which(!is.na(dataset_mat), arr.ind = T)
  dataset_triplets = matrix(0, nrow = nrow(inds), ncol = 3)
  dataset_triplets[,c(1,2)] = inds
  dataset_triplets[,3] = dataset_mat[inds]
  
  return (dataset_triplets)
}

train_crf = function(training_data, adj_mats, X, crf_iters, eta, list_parameters_every_its){
  
  n = length(training_data)
  ## make Bk mats
  Bk_mats = list()
  for (i in 1:length(adj_mats)){
    Bk_mats[[i]] = -1*adj_mats[[i]]
    for (k in 1:n){
      Bk_mats[[i]][k,k] = sum(adj_mats[[i]][k,])-adj_mats[[i]][k,k]
    }
  }
  
  ## get A and B mats to compute gradient
  
  alpha = 0.01
  beta = list()
  beta_old = list()
  for (i in 1:length(adj_mats)){
    beta[[i]] = 0.01
  }
  y_vec = training_data
  for (it in 1:crf_iters){
    if ((it%%list_parameters_every_its)==0){
      cat('training iteration ',it,'\n')
      cat('alpha: ',alpha,'\n')
      for (k in 1:length(beta)){
        cat('beta ',k,': ',beta[[k]],'\n')
      }
    }
    
    A_train = sparseMatrix(i = c(1:n), j = c(1:n), x = rep(alpha, n))
    B_train = beta[[1]]*Bk_mats[[1]]
    if (length(Bk_mats)>1){
      for (i in 2:length(Bk_mats)){
        B_train = B_train+beta[[i]]*Bk_mats[[i]]
      }
    }
    
    Sigma1 = 2*(A_train+B_train)
    temp = chol(Sigma1) 
    Sigma = chol2inv(temp)
    
    b = 2*X*alpha
    mu = Sigma%*%b
    
    log_alpha = log(alpha)
    grad_alpha = -t(y_vec)%*%y_vec+2*t(y_vec)%*%X-2*t(X)%*%mu+t(mu)%*%mu+sum(diag(Sigma))
    grad_log_alpha = alpha*grad_alpha
    log_alpha = log_alpha + eta * grad_log_alpha
    alpha_old = alpha
    alpha = as.numeric(exp(log_alpha))
    
    for (i in 1:length(Bk_mats)){
      log_beta = log(beta[[i]])
      grad_beta = -t(y_vec)%*%Bk_mats[[i]]%*%y_vec+t(mu)%*%Bk_mats[[i]]%*%mu+t(as.vector(Sigma)%*%as.vector(Bk_mats[[i]]))
      grad_log_beta = beta[[i]]*grad_beta
      log_beta = log_beta + eta*grad_log_beta
      beta_old[[i]] = beta[[i]]
      beta[[i]] = as.numeric(exp(log_beta))
    }
    
    
    if(alpha>10){
      alpha=10
    }
    
    counter = 0
    for (i in 1:length(Bk_mats)){
      if (beta[[i]]>10){
        beta[[i]]=10
      }
      if(abs(beta_old[[i]]-beta[[i]])<0.00001){
        counter = counter+1
      }
    }
    
    if (counter==length(Bk_mats) && abs(alpha_old-alpha)){
      cat('stopping after training iteration ',it,'\n')
      cat('alpha: ',alpha,'\n')
      for (k in 1:length(beta)){
        cat('beta ',k,': ',beta[[k]],'\n')
      }
      return(list(alpha,beta))
    }
  }
  return(list(alpha,beta))
}

train_crf_ = function(training_data, adj_mats, X, crf_iters, eta, list_parameters_every_its, inds_){
  
  n = length(training_data)
  ## make Bk mats
  Bk_mats = list()
  for (i in 1:length(adj_mats)){
    Bk_mats[[i]] = -1*adj_mats[[i]]
    for (k in 1:n){
      Bk_mats[[i]][k,k] = sum(adj_mats[[i]][k,])-adj_mats[[i]][k,k]
    }
  }
  
  ## get A and B mats to compute gradient
  
  alpha = 0.01
  beta = list()
  beta_old = list()
  for (i in 1:length(adj_mats)){
    beta[[i]] = 0.01
  }
  y_vec = training_data
  for (it in 1:crf_iters){
    if ((it%%list_parameters_every_its)==0){
      cat('training iteration ',it,'\n')
      cat('alpha: ',alpha,'\n')
      for (k in 1:length(beta)){
        cat('beta ',k,': ',beta[[k]],'\n')
      }
    }
    
    temp = rep(alpha, n)
    #temp[inds_] = 0.000000001
    A_train = sparseMatrix(i = c(1:n), j = c(1:n), x = temp)
    B_train = beta[[1]]*Bk_mats[[1]]
    if (length(Bk_mats)>1){
      for (i in 2:length(Bk_mats)){
        B_train = B_train+beta[[i]]*Bk_mats[[i]]
      }
    }
    
    Sigma1 = 2*(A_train+B_train)
    temp = chol(Sigma1) 
    Sigma = chol2inv(temp)
    
    b = 2*X*alpha
    mu = Sigma%*%b
    
    log_alpha = log(alpha)
    grad_alpha = -t(y_vec)%*%y_vec+2*t(y_vec)%*%X-2*t(X)%*%mu+t(mu)%*%mu+sum(diag(Sigma))
    grad_log_alpha = alpha*grad_alpha
    log_alpha = log_alpha + eta * grad_log_alpha
    alpha_old = alpha
    alpha = as.numeric(exp(log_alpha))
    
    for (i in 1:length(Bk_mats)){
      log_beta = log(beta[[i]])
      grad_beta = -t(y_vec)%*%Bk_mats[[i]]%*%y_vec+t(mu)%*%Bk_mats[[i]]%*%mu+t(as.vector(Sigma)%*%as.vector(Bk_mats[[i]]))
      grad_log_beta = beta[[i]]*grad_beta
      log_beta = log_beta + eta*grad_log_beta
      beta_old[[i]] = beta[[i]]
      beta[[i]] = as.numeric(exp(log_beta))
    }
    
    
    if(alpha>10){
      alpha=10
    }
    
    counter = 0
    for (i in 1:length(Bk_mats)){
      if (beta[[i]]>10){
        beta[[i]]=10
      }
      if(abs(beta_old[[i]]-beta[[i]])<0.00001){
        counter = counter+1
      }
    }
    
    if (counter==length(Bk_mats) && abs(alpha_old-alpha)){
      cat('stopping after training iteration ',it,'\n')
      cat('alpha: ',alpha,'\n')
      for (k in 1:length(beta)){
        cat('beta ',k,': ',beta[[k]],'\n')
      }
      return(list(alpha,beta))
    }
  }
  return(list(alpha,beta))
}

## train_data is a vector, where the values to be predicted are NA
## the values in train_data which are not NA are in the training data
crf_predict = function(train_data, X, alpha, adj_mats, beta){
  
  n = length(train_data)
  A = sparseMatrix(i = c(1:n), j = c(1:n), x = rep(alpha,n))
  B = -1*beta[[1]]*adj_mats[[1]]
  if (length(adj_mats)>1){
    for (i in 2:length(adj_mats)){
      B = B-beta[[i]]*adj_mats[[i]]
    }
  }
  for (i in 1:n){
    sum = 0
    for (k in 1:length(adj_mats)){
      sum = sum+beta[[k]]*sum(adj_mats[[k]][i,])-beta[[k]]*adj_mats[[k]][i,i]
    }
    B[i,i] = sum
  }
  
  Sigma1 = 2*(A+B)
  Sigma = chol2inv(chol(Sigma1))
  
  b = 2*X*alpha
  
  mu = Sigma%*%b
  
  unknown = which(train_data<0)
  known = which(train_data>=0)
  
  Sigma12 = Sigma[unknown, known]
  Sigma22 = Sigma[known, known]
  Sigma221 = chol2inv(chol(Sigma22))
  mu_ = mu[unknown] + Sigma12%*%Sigma221%*%(train_data[known] - mu[known])
  mu_all = rep(0, length(mu))
  mu_all[unknown] = as.vector(mu_)
  mu_all[known] = train_data[known]
  
  return (mu_all)
  
}

crf_predict_ = function(train_data, X, alpha, adj_mats, beta, inds_){
  
  n = length(train_data)
  temp = rep(alpha,n)
  
  temp = rep(alpha, n)
  A = sparseMatrix(i = c(1:n), j = c(1:n), x = temp)
  B = -1*beta[[1]]*adj_mats[[1]]
  if (length(adj_mats)>1){
    for (i in 2:length(adj_mats)){
      B = B-beta[[i]]*adj_mats[[i]]
    }
  }
  for (i in 1:n){
    sum = 0
    for (k in 1:length(adj_mats)){
      sum = sum+beta[[k]]*sum(adj_mats[[k]][i,])-beta[[k]]*adj_mats[[k]][i,i]
    }
    B[i,i] = sum
  }
  
  Sigma1 = 2*(A+B)
  Sigma = chol2inv(chol(Sigma1))
  
  b = 2*X*alpha
  
  mu = Sigma%*%b
  
  unknown = which(train_data<0)
  known = which(train_data>=0)
  
  Sigma12 = Sigma[unknown, known]
  Sigma22 = Sigma[known, known]
  Sigma221 = chol2inv(chol(Sigma22))
  mu_ = mu[unknown] + Sigma12%*%Sigma221%*%(train_data[known] - mu[known])
  mu_all = rep(0, length(mu))
  mu_all[unknown] = as.vector(mu_)
  mu_all[known] = train_data[known]
  
  return (mu_all)
  
}

make_adj_mat_thresh = function(sim_mat, thresh, add_nn, nn){
  
  n = nrow(sim_mat)
  adj_mat = matrix(0, nrow = n, ncol = n)
  for (i in 1:n){
    inds = which(sim_mat[i,]>thresh)
    inds = setdiff(inds, i)
    if (add_nn==TRUE){
      if (length(inds)<1){
        inds = order(sim_mat[i,], decreasing = T)[1:nn]
        inds = setdiff(inds, i)
      }
    }
    adj_mat[i, inds] = sim_mat[i, inds]
    adj_mat[inds, i] = sim_mat[inds, i]
  }
  
  return(adj_mat)
}

get_mf_cv = function(train_mat, iters){
  
  n_drugs = nrow(train_mat)
  n_targets = ncol(train_mat)
  
  triplets = matrix(0, nrow = length(which(train_mat>=0)), ncol = 3)
  triplets[,1] = which(train_mat>=0, arr.ind = T)[,1]
  triplets[,2] = which(train_mat>=0, arr.ind = T)[,2]
  triplets[,3] = train_mat[which(train_mat>=0, arr.ind = T)]
  
  train = triplets
  k = 5
  n = nrow(train)
  folds = vector(k, mode='list')
  shuf = sample(n)
  m = n %/% k
  for (i in 1:(k-1)) {
    folds[[i]] = shuf[1:m]
    shuf = setdiff(shuf, folds[[i]])
  }
  folds[[k]] = shuf
  preds = rep(0,n)
  
  train.triplet = train[,1:3]
  train.triplet[,1] = train.triplet[,1] - 1
  train.triplet[,2] = train.triplet[,2] - 1
  
  for (i in 1:k) {
    teind = folds[[i]]
    trind = setdiff(1:n,teind)
    
    res = libmf(train.triplet[trind,], m = n_drugs, n = n_targets, k = 20, cost = 0.01, lrate = 0.01,
                
                niter = iters, nthread = 1, nmf = FALSE, verbose = FALSE)
    
    P = res[[2]][[1]]
    Q = res[[2]][[2]]
    
    estM = P %*% t(Q)
   
    tmp.pred = estM[train[teind,1:2]]
    preds[teind] = tmp.pred
  }
  
  preds_mat = matrix(NA, nrow = n_drugs, ncol = n_targets)
  for (i in 1:nrow(triplets)){
    preds_mat[triplets[i,1], triplets[i,2]] = preds[i]
  }
  
  return(list(preds, preds_mat))
  
}

get_libmf_prediction = function(train, iter){
  
  n_drugs = nrow(train)
  n_targets = ncol(train)
  
  train_mat_triplet = matrix(0, nrow = length(which(train>=0)), ncol = 3)
  train_mat_triplet[,1] = which(train>=0, arr.ind = T)[,1]
  train_mat_triplet[,2] = which(train>=0, arr.ind = T)[,2]
  train_mat_triplet[,3] = train[which(train>=0, arr.ind = T)]
  
  train_mat_triplet[,1] = train_mat_triplet[,1] - 1
  train_mat_triplet[,2] = train_mat_triplet[,2] - 1
  
  res = libmf(train_mat_triplet, m = n_drugs, n = n_targets, k = 20, cost = 0.01, lrate = 0.01,
              
              niter = iter, nthread = 1, nmf = FALSE, verbose = FALSE)
  
  P = res[[2]][[1]]
  Q = res[[2]][[2]]
  
  estM = P %*% t(Q)
  
  return(estM)
  
}

## this function creates an adjacency mat for the training data in training_mat
## training_mat is a matrix with the targets in target_set and the drugs in drug_set
## training_mat contains NA for missing values
make_training_adj_mat = function(train_mat, drug_adj_mat, target_adj_mat, drug_set, target_set){
  
  temp_adj_mat_drugs = matrix(0, nrow = nrow(train_mat)*ncol(train_mat), ncol = nrow(train_mat)*ncol(train_mat))
  temp_adj_mat_targets = matrix(0, nrow = nrow(train_mat)*ncol(train_mat), ncol = nrow(train_mat)*ncol(train_mat))
  
  ##integrate the target similarities
  for (i in 1:nrow(train_mat)){
    from = (i-1)*ncol(train_mat)+1
    to = i*ncol(train_mat)
    temp_adj_mat_targets[from:to,from:to] = target_adj_mat[target_set, target_set]
  }
  ## integrate the drug similarities
  for (i in 1:ncol(train_mat)){
    from = i
    to = (nrow(train_mat)-1)*ncol(train_mat)+i
    temp = seq(from = from, to = to, by = ncol(train_mat))
    temp_adj_mat_drugs[temp,temp] = drug_adj_mat[drug_set, drug_set]
  }
  
  inds = which(t(train_mat)>=0)
  temp_adj_mat_drugs = temp_adj_mat_drugs[inds,inds]
  inds = which(temp_adj_mat_drugs>0, arr.ind = T)
  sparse_adj_mat_drugs = sparseMatrix(i = inds[,1], j = inds[,2], x = temp_adj_mat_drugs[inds], dims = dim(temp_adj_mat_drugs))
  
  inds = which(t(train_mat)>=0)
  temp_adj_mat_targets = temp_adj_mat_targets[inds,inds]
  inds = which(temp_adj_mat_targets>0, arr.ind = T)
  sparse_adj_mat_targets = sparseMatrix(i = inds[,1], j = inds[,2], x = temp_adj_mat_targets[inds], dims = dim(temp_adj_mat_targets))
  
  return (list(sparse_adj_mat_drugs, sparse_adj_mat_targets))
  
}

## this function creates an adjacency mat for the training data in training_mat
## training_mat is a matrix with the targets in target_set and the drugs in drug_set
## training_mat contains NA for missing values
make_training_adj_mat_ = function(train_mat, drug_adj_mat, target_adj_mat, drug_set, target_set, inds_){
  
  temp_adj_mat_drugs = matrix(0, nrow = nrow(train_mat)*ncol(train_mat), ncol = nrow(train_mat)*ncol(train_mat))
  temp_adj_mat_targets = matrix(0, nrow = nrow(train_mat)*ncol(train_mat), ncol = nrow(train_mat)*ncol(train_mat))
  
  ##integrate the target similarities
  for (i in 1:nrow(train_mat)){
    from = (i-1)*ncol(train_mat)+1
    to = i*ncol(train_mat)
    temp_adj_mat_targets[from:to,from:to] = target_adj_mat[target_set, target_set]
  }
  ## integrate the drug similarities
  for (i in 1:ncol(train_mat)){
    from = i
    to = (nrow(train_mat)-1)*ncol(train_mat)+i
    temp = seq(from = from, to = to, by = ncol(train_mat))
    temp_adj_mat_drugs[temp,temp] = drug_adj_mat[drug_set, drug_set]
  }
  
  inds = inds_
  temp_adj_mat_drugs = temp_adj_mat_drugs[inds,inds]
  inds = which(temp_adj_mat_drugs>0, arr.ind = T)
  sparse_adj_mat_drugs = sparseMatrix(i = inds[,1], j = inds[,2], x = temp_adj_mat_drugs[inds], dims = dim(temp_adj_mat_drugs))
  
  inds = inds_
  temp_adj_mat_targets = temp_adj_mat_targets[inds,inds]
  inds = which(temp_adj_mat_targets>0, arr.ind = T)
  sparse_adj_mat_targets = sparseMatrix(i = inds[,1], j = inds[,2], x = temp_adj_mat_targets[inds], dims = dim(temp_adj_mat_targets))
  
  return (list(sparse_adj_mat_drugs, sparse_adj_mat_targets))
  
}

make_training_adj_mat__ = function(train_mat, drug_adj_mat, target_adj_mat, drug_set, target_set, inds_, d){
  
  temp_adj_mat_drugs = matrix(0, nrow = nrow(train_mat)*ncol(train_mat), ncol = nrow(train_mat)*ncol(train_mat))
  temp_adj_mat_targets = matrix(0, nrow = nrow(train_mat)*ncol(train_mat), ncol = nrow(train_mat)*ncol(train_mat))
  
  ##integrate the target similarities
  #for (i in 1:nrow(train_mat)){
    i = match(d, drug_set)
    from = (i-1)*ncol(train_mat)+1
    to = i*ncol(train_mat)
    temp_adj_mat_targets[from:to,from:to] = target_adj_mat[target_set, target_set]
  #}
  ## integrate the drug similarities
  for (i in 1:ncol(train_mat)){
    from = i
    to = (nrow(train_mat)-1)*ncol(train_mat)+i
    temp = seq(from = from, to = to, by = ncol(train_mat))
    temp_adj_mat_drugs[temp,temp] = drug_adj_mat[drug_set, drug_set]
  }
  
  inds = inds_
  temp_adj_mat_drugs = temp_adj_mat_drugs[inds,inds]
  inds = which(temp_adj_mat_drugs>0, arr.ind = T)
  sparse_adj_mat_drugs = sparseMatrix(i = inds[,1], j = inds[,2], x = temp_adj_mat_drugs[inds], dims = dim(temp_adj_mat_drugs))
  
  inds = inds_
  temp_adj_mat_targets = temp_adj_mat_targets[inds,inds]
  inds = which(temp_adj_mat_targets>0, arr.ind = T)
  sparse_adj_mat_targets = sparseMatrix(i = inds[,1], j = inds[,2], x = temp_adj_mat_targets[inds], dims = dim(temp_adj_mat_targets))
  
  return (list(sparse_adj_mat_drugs, sparse_adj_mat_targets))
  
}

make_prediction_adj_mat = function(train_mat, drug_adj_mat, target_adj_mat, drug_set, target_set){
  
  temp_adj_mat_drugs = matrix(0, nrow = nrow(train_mat)*ncol(train_mat), ncol = nrow(train_mat)*ncol(train_mat))
  temp_adj_mat_targets = matrix(0, nrow = nrow(train_mat)*ncol(train_mat), ncol = nrow(train_mat)*ncol(train_mat))
  
  ##integrate the target similarities
  for (i in 1:nrow(train_mat)){
    from = (i-1)*ncol(train_mat)+1
    to = i*ncol(train_mat)
    temp_adj_mat_targets[from:to,from:to] = target_adj_mat[target_set, target_set]
  }
  ## integrate the drug similarities
  for (i in 1:ncol(train_mat)){
    from = i
    to = (nrow(train_mat)-1)*ncol(train_mat)+i
    temp = seq(from = from, to = to, by = ncol(train_mat))
    temp_adj_mat_drugs[temp,temp] = drug_adj_mat[drug_set, drug_set]
  }
  
  inds = which(temp_adj_mat_drugs>0, arr.ind = T)
  sparse_adj_mat_drugs = sparseMatrix(i = inds[,1], j = inds[,2], x = temp_adj_mat_drugs[inds], dims = dim(temp_adj_mat_drugs))
  
  inds = which(temp_adj_mat_targets>0, arr.ind = T)
  sparse_adj_mat_targets = sparseMatrix(i = inds[,1], j = inds[,2], x = temp_adj_mat_targets[inds], dims = dim(temp_adj_mat_targets))
  
  return (list(sparse_adj_mat_drugs, sparse_adj_mat_targets))
  
}

make_prediction_adj_mat_ = function(train_mat, drug_adj_mat, target_adj_mat, drug_set, target_set, d){
  
  temp_adj_mat_drugs = matrix(0, nrow = nrow(train_mat)*ncol(train_mat), ncol = nrow(train_mat)*ncol(train_mat))
  temp_adj_mat_targets = matrix(0, nrow = nrow(train_mat)*ncol(train_mat), ncol = nrow(train_mat)*ncol(train_mat))
  
  ##integrate the target similarities
  #for (i in 1:nrow(train_mat)){
    i = match(d, drug_set)
    from = (i-1)*ncol(train_mat)+1
    to = i*ncol(train_mat)
    temp_adj_mat_targets[from:to,from:to] = target_adj_mat[target_set, target_set]
  #}
  ## integrate the drug similarities
  for (i in 1:ncol(train_mat)){
    from = i
    to = (nrow(train_mat)-1)*ncol(train_mat)+i
    temp = seq(from = from, to = to, by = ncol(train_mat))
    temp_adj_mat_drugs[temp,temp] = drug_adj_mat[drug_set, drug_set]
  }
  
  inds = which(temp_adj_mat_drugs>0, arr.ind = T)
  sparse_adj_mat_drugs = sparseMatrix(i = inds[,1], j = inds[,2], x = temp_adj_mat_drugs[inds], dims = dim(temp_adj_mat_drugs))
  
  inds = which(temp_adj_mat_targets>0, arr.ind = T)
  sparse_adj_mat_targets = sparseMatrix(i = inds[,1], j = inds[,2], x = temp_adj_mat_targets[inds], dims = dim(temp_adj_mat_targets))
  
  return (list(sparse_adj_mat_drugs, sparse_adj_mat_targets))
  
}

get_metrics = function(preds, labels, cutoff){
  rmse = sqrt(mean((preds - labels)^2))
  aupr.val = NA
  if (length(which(labels>=cutoff))>0 && length(preds)>1 && length(which(labels<cutoff))>0){
    
    pred.obj = ROCR::prediction(preds, as.numeric(labels>=cutoff))
    auc.obj = ROCR::performance(pred.obj, measure = 'auc')
    prec.obj = ROCR::performance(pred.obj, measure = 'prec')
    rec.obj = ROCR::performance(pred.obj, measure = 'rec')
    prec.val = prec.obj@y.values[[1]]
    rec.val = rec.obj@y.values[[1]]
    auc.val = auc.obj@y.values[[1]]
    
    if (is.na(prec.val[1])){
      prec.val[1] = 1
    }
    
    ## comment this out if aupr computation crashes
    func = approxfun(cbind(rec.val,prec.val), yleft = 1)
    #aupr.val = integrate(func, 0, 1, subdivisions = 1000L)$value
    #tryCatch(integrate(func, 0, 1, subdivisions = 1000L)$value, finally = (aupr.val = NA))
    try((aupr.val = integrate(func, 0, 1, subdivisions = 1000L)$value), silent = TRUE)
    ##
  } else {
    auc.val=NA
    aupr.val=NA
  }
  
  #return(list(rmse,auc.val,NA))
  return(list(rmse,auc.val,aupr.val))
  
}
