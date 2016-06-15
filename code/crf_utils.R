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

make_adjacency_mat_targets = function(sim_mat, thresh){
  
  n_targets = nrow(sim_mat)
  
  adj_mat = matrix(0, nrow = n_targets, ncol = n_targets)
  for (i in 1:n_targets){
    inds = which(sim_mat[i,]>thresh & sim_mat[i,]<1)
    adj_mat[i, inds] = 1
    adj_mat[inds,i ] = 1
    #adj_mat[i, inds] = sim_mat[i, inds]
    #adj_mat[inds, i] = sim_mat[inds, i]
  }
  
  return(adj_mat)
  
}

make_adjacency_mat = function(sim_mat, thresh){
  
  n_drugs = nrow(sim_mat)
  
  adj_mat = matrix(0, nrow = n_drugs, ncol = n_drugs)
  for (i in 1:n_drugs){
    inds = which(sim_mat[i,]>thresh & sim_mat[i,]<1)
    #inds = order(sim_mat[i,], decreasing = T)[2:5]
    #if (length(inds)<4){
    #  inds = which(sim_mat[i,]>0.8 & sim_mat[i,]<1)
    #}
    #if (length(inds)<4){
    #  inds = order(sim_mat[i,], decreasing = T)[1:5]
    #  inds = setdiff(inds, i)
    #}
    adj_mat[i, inds] = 1
    adj_mat[inds,i ] = 1
    #adj_mat[i, inds] = simi_mat[i, inds]
    #adj_mat[inds, i] = simi_mat[inds, i]
  }
  
  return(adj_mat)
  
}

make_adjacency_mat_4n = function(sim_mat){
  
  n_drugs = nrow(sim_mat)
  
  adj_mat = matrix(0, nrow = n_drugs, ncol = n_drugs)
  for (i in 1:n_drugs){
    inds = which(sim_mat[i,]>0.9 & sim_mat[i,]<1)
    #inds = order(sim_mat[i,], decreasing = T)[2:5]
    if (length(inds)<1){
      inds = order(sim_mat[i,], decreasing = T)[2]
      inds = setdiff(inds, i)
    }
    adj_mat[i, inds] = 1
    adj_mat[inds,i ] = 1
    #adj_mat[i, inds] = simi_mat[i, inds]
    #adj_mat[inds, i] = simi_mat[inds, i]
  }
  
  return(adj_mat)
  
}

make_B_mats_simple = function(training_data_mat, row_adj_mat, col_adj_mat){
  
  n_drugs = nrow(training_data_mat)
  n_targets = ncol(training_data_mat)
  
  n_train = length(which(training_data_mat>=0))
  
  inds = which(t(training_data_mat)>=0, arr.ind = T)
  
  B_row = matrix(0, nrow = n_train, ncol = n_train)
  for (d in 1:n_drugs){
    for (t in 1:n_targets){
      if (training_data_mat[d,t]>=0){
        given_targets = which(training_data_mat[d,]>=0)
        given_targets = setdiff(given_targets,t)
        neighbors = which(col_adj_mat[t,]>0)
        given_neighbors = given_targets[which(given_targets %in% neighbors)]
        cell_a = which(inds[,1]==t & inds[,2]==d)
        for (gn in given_neighbors){
          cell_b = which(inds[,1]==gn & inds[,2]==d)
          B_row[cell_a, cell_b] = -1
        }
        B_row[cell_a, cell_a] = length(given_neighbors)
      }
    }
  }
  
  B_col = matrix(0, nrow = n_train, ncol = n_train)
  for (d in 1:n_drugs){
    for (t in 1:n_targets){
      if (training_data_mat[d,t]>=0){
        given_drugs = which(training_data_mat[,t]>=0)
        given_drugs = setdiff(given_drugs,d)
        neighbors = which(row_adj_mat[d,]>0)
        given_neighbors = given_drugs[which(given_drugs %in% neighbors)]
        cell_a = which(inds[,1]==t & inds[,2]==d)
        for (gn in given_neighbors){
          cell_b = which(inds[,1]==t & inds[,2]==gn)
          B_col[cell_a, cell_b] = -1
        }
        B_col[cell_a, cell_a] = length(given_neighbors)
      }
    }
  }
  
  return(list(B_row, B_col))
  
}

get_gaussian_train_simple = function(training_data_mat, X_vec, alpha_train, row_beta_train, col_beta_train, B_row, B_col){
  
  n_train = length(which(training_data_mat>0))
  A_train = matrix(0, nrow = n_train, ncol = n_train)
  B_train = matrix(0, nrow = n_train, ncol = n_train)
  
  B_train = row_beta_train*B_row + col_beta_train*B_col
  for (i in 1:n_train){
    A_train[i,i] = alpha_train
  }
  
  Sigma1 = 2*(A_train+B_train)
  require("MASS")
  temp = chol(Sigma1) 
  Sigma = chol2inv(temp)
  
  b = 2*X_vec*alpha_train
  mu = Sigma%*%b
  
  return (list(Sigma, mu))
  
}

get_gradient_alpha = function(y_vec, X_vec, mu, Sigma){
  
  gradient_alpha = -t(y_vec)%*%y_vec+2*t(y_vec)%*%X_vec-2*t(X_vec)%*%mu+t(mu)%*%mu+sum(diag(Sigma))
  
  return(as.numeric(gradient_alpha))
  
}

get_B_gradient = function(y_vec, B, mu, Sigma){
  
  gradient_beta = -t(y_vec)%*%B%*%y_vec+t(mu)%*%B%*%mu+t(as.vector(Sigma)%*%as.vector(B))
  
  return(as.numeric(gradient_beta))
}

simple_crf_predict = function(training_data_mat, row_adj_mat, col_adj_mat, alpha, row_beta, col_beta, X_mat){
  
  n_drugs = nrow(training_data_mat)
  n_targets = ncol(training_data_mat)
  
  n_triplets_B = length(which(row_adj_mat>0))*n_targets + length(which(col_adj_mat>0))*n_drugs + n_drugs*n_targets
  n_triplets_A = n_drugs*n_targets
  triplets_A = matrix(0, nrow = n_triplets_A, ncol = 3)
  triplets_B = matrix(0, nrow = n_triplets_B, ncol = 3)
  counter_A = 1
  counter_B = 1
  for (d in 1:n_drugs){
    for (t in 1:n_targets){
      
      row_neighbors = which(row_adj_mat[d,]==1)
      col_neighbors = which(col_adj_mat[t,]==1)
      temp_a = (d-1)*n_targets+t
      for (rn in row_neighbors){
        temp_b = (rn-1)*n_targets+t
        triplets_B[counter_B,] = c(temp_a, temp_b, -col_beta)
        counter_B = counter_B+1
      }
      for (cn in col_neighbors){
        temp_b = (d-1)*n_targets+cn
        triplets_B[counter_B,] = c(temp_a, temp_b, -row_beta)
        counter_B = counter_B+1
      }
      nrow_neighbors = length(which(row_adj_mat[d,]>0))
      ncol_neighbors = length(which(col_adj_mat[t,]>0))
      sum = col_beta*nrow_neighbors+row_beta*ncol_neighbors
      triplets_A[counter_A,] = c(temp_a, temp_a, alpha)
      triplets_B[counter_B,] = c(temp_a, temp_a, sum)
      counter_B = counter_B+1
      counter_A = counter_A+1
    }
  }
  A = sparseMatrix(i = triplets_A[,1], j = triplets_A[,2], x = triplets_A[,3])
  B = sparseMatrix(i = triplets_B[,1], j = triplets_B[,2], x = triplets_B[,3])
  
  Sigma1 = 2*(A+B)
  inds = which(Sigma1!=0, arr.ind = T)
  sm = sparseMatrix(i = inds[,1], j = inds[,2], x = Sigma1[inds])
  sparsity = nrow(inds)/(nrow(Sigma1)*ncol(Sigma1))
  cat('sparsity of matrix to invert: ', sparsity,'\n')
  Sigma = chol2inv(chol(sm))
  
  X = as.vector(t(X_mat))
  b = 2*X*alpha
  
  mu = Sigma%*%b
  
  unknown = which(as.vector(t(training_data_mat))<0)
  known = which(as.vector(t(training_data_mat))>=0)
  
  Sigma12 = Sigma[unknown, known]
  Sigma22 = Sigma[known, known]
  Sigma221 = chol2inv(chol(Sigma22))
  mu_ = mu[unknown] + Sigma12%*%Sigma221%*%(as.vector(t(training_data_mat))[known] - mu[known])
  mu_all = rep(0, length(mu))
  mu_all[unknown] = as.vector(mu_)
  mu_all[known] = as.vector(t(dt_mat_temp))[known]
  
  pred_mat = matrix(mu_all, nrow = n_drugs, byrow = T)
  
  return (pred_mat)
}

train_crf_dual_simple = function(training_data_mat, row_adj_mat, col_adj_mat, mf_preds_train_mat, crf_iters, eta){
  
  ans = make_B_mats_simple(training_data_mat, row_adj_mat, col_adj_mat)
  B_row = ans[[1]]
  B_col = ans[[2]]
  
  y_vec = t(training_data_mat)[which(t(training_data_mat)>=0)]
  X_train = t(mf_preds_train_mat)[which(!is.na(t(mf_preds_train_mat)))]
  
  eta = eta
  row_beta_train = 0.001
  col_beta_train = 0.001
  alpha_train = 0.001
  
  for (it in 1:crf_iters){
    
    if ((it%%100)==0){
      cat('iteration ',it,'.. alpha ',alpha_train,' target beta ',row_beta_train,' drug beta ',col_beta_train,'\n')
    }
    #cat('iteration ',it,'.. alpha ',alpha_train,' row beta ',row_beta_train,' col beta ',col_beta_train,'\n')
    
    ans = get_gaussian_train_simple(training_data_mat = training_data_mat, X_vec = X_train, alpha_train = alpha_train, row_beta_train = row_beta_train, col_beta_train = col_beta_train, B_row = B_row, B_col = B_col)
    Sigma = ans[[1]]
    mu = ans[[2]]
    
    log_alpha_train = log(alpha_train)
    grad_alpha_train = get_gradient_alpha(y_vec = y_vec, X_vec = X_train, mu = mu, Sigma = Sigma)
    grad_log_alpha_train = alpha_train*grad_alpha_train
    log_alpha_train = log_alpha_train + eta * grad_log_alpha_train
    alpha_old = alpha_train
    alpha_train = exp(log_alpha_train)
    
    log_row_beta_train = log(row_beta_train)
    grad_row_beta_train = get_B_gradient(y_vec = y_vec, B = B_row, mu = mu, Sigma = Sigma)
    grad_log_row_beta_train = row_beta_train*grad_row_beta_train
    log_row_beta_train = log_row_beta_train + eta*grad_log_row_beta_train
    row_beta_old = row_beta_train
    row_beta_train = exp(log_row_beta_train)
    
    log_col_beta_train = log(col_beta_train)
    grad_col_beta_train = get_B_gradient(y_vec = y_vec, B = B_col, mu = mu, Sigma = Sigma)
    grad_log_col_beta_train = col_beta_train*grad_col_beta_train
    log_col_beta_train = log_col_beta_train + eta*grad_log_col_beta_train
    col_beta_old = col_beta_train
    col_beta_train = exp(log_col_beta_train)
    
    if (row_beta_train<1e-10){
      row_beta_train = 1e-10
    } else if (row_beta_train>2){
      row_beta_train = 2
    }
    
    if (col_beta_train<1e-10){
      col_beta_train = 1e-10
    } else if (col_beta_train>2){
      col_beta_train = 2
    }
    
    if (alpha_train<1e-10){
      alpha_train = 1e-10
    } else if (alpha_train>2){
      alpha_train = 2
    }
    
    if (abs(alpha_old-alpha_train)<0.00001 && abs(row_beta_old-row_beta_train)<0.00001 && abs(col_beta_old-col_beta_train)<0.00001){
      cat('learned params after iteration ',it,'.. alpha ',alpha_train,' target beta ',row_beta_train,' drug beta ',col_beta_train,'\n')
      
      return (list(alpha_train, row_beta_train, col_beta_train))
    }
    
  }
  cat('learned params after iteration ',it,'.. alpha ',alpha_train,' target beta ',row_beta_train,' drug beta ',col_beta_train,'\n')
  
  return (list(alpha_train, row_beta_train, col_beta_train))
  
}

make_training_adj_mat_for_column = function(dt_mat, sim_mat, col){
  
  n_train = length(which(dt_mat[,col]>=0))
  adj_temp = matrix(0, nrow = n_train, ncol = n_train)
  temp = which(dt_mat[,col]>0)
  for (k in temp){
    #cat('got observation for drug ',k,'\n')
    neighbors = which(sim_mat[k,]>0.9 & sim_mat[k,]<1)
    #neighbors = order(simi_mat[k,], decreasing = T)[2:5]
    if (length(neighbors)<4){
      inds = order(sim_mat[i,], decreasing = T)[1:5]
      inds = setdiff(inds, i)
      #neighbors = order(simi_mat[k,], decreasing = T)[2:5]
      neighbors = inds
    }
    #cat('which is adj to ',neighbors,'\n')
   
    for (n in neighbors){
      if (n %in% temp){
        adj_temp[match(k, temp), match(n, temp)] = 1
        adj_temp[match(n, temp), match(k, temp)] = 1
      }
    }
  }
  return (adj_temp)
}

train_crf_column = function(y, X, adj_mat, crf_iters, eta){
  
  n_train = length(y)
  
  B1 = matrix(0, nrow = n_train, ncol = n_train)
  for (i in 1:n_train){
    for (j in 1:n_train){
      if (i==j){
        B1[i,j] = sum(adj_mat[i,])
      } else {
        B1[i,j] = -adj_mat[i,j]
      }
    }
  }
  
  alpha = 0.01
  beta = 0.01
  
  log_alpha = log(alpha)
  log_beta = log(beta)
  
  for (it in 1:crf_iters){
    
    #cat('training iteration ',it,':',alpha,' ',beta,'\n')
    if ((it%%100)==0){
      cat('training iteration ',it,':',alpha,' ',beta,'\n')
    }
    
    gradients = compute_gradient(y,alpha, beta, adj_mat, B1, X, n_train)
    grad_alpha = gradients[[1]]
    grad_beta = gradients[[2]]
    
    grad_log_alpha = alpha*grad_alpha
    grad_log_beta = beta*grad_beta
    
    log_alpha = log_alpha + eta * grad_log_alpha
    log_beta = log_beta + eta * grad_log_beta
    
    alpha_old = alpha
    beta_old = beta
    
    alpha = exp(log_alpha)
    beta = exp(log_beta)
    
    if (beta>10){
      beta = 10
    }
    
    if (abs(alpha_old-alpha)<0.00001 && abs(beta_old-beta)<0.00001){
      cat('training iteration ',it,':',alpha,' ',beta,'\n')
      return(list(alpha,beta))
    }
    
  }
  return(list(alpha,beta))
}

make_crf_predictions_col = function(alpha, beta, column, adj_mat, X){
  
  n_drugs = nrow(adj_mat)
  
  A = matrix(0, nrow = n_drugs, ncol = n_drugs)
  B = matrix(0, nrow = n_drugs, ncol = n_drugs)
  
  #for (i in 1:n_drugs){
  #  for (j in 1:n_drugs){
  #    if (i==j){
  #      A[i,j] = sum(alpha)
  #      B[i,j] = beta*sum(adj_mat[i,])
  #    } else {
  #      B[i,j] = -beta*adj_mat[i,j]
  #    }
  #  }
  #}
  
  for (i in 1:n_drugs){
    A[i,i] = sum(alpha)
    B[i,i] = beta*sum(adj_mat[i,])
    neighbors = which(adj_mat[i,]>0)
    for (n in neighbors){
      B[i,n] = -beta
    }
  }
  
  #cat('inverting matrix..\n')
  Sigma1 = 2*(A+B)
  require("MASS")
  inds = which(Sigma1!=0, arr.ind = T)
  sm = sparseMatrix(i = inds[,1], j = inds[,2], x = Sigma1[inds])
  #temp = chol(Sigma1) 
  Sigma = chol2inv(chol(sm))
  #Sigma = chol2inv(temp)
  
  b = 2*X*alpha
  mu = Sigma%*%b
  
  train = column
  
  n_drugs = length(train)
  
  unknown = which(as.vector(t(train))<0)
  known = which(as.vector(t(train))>=0)
  
  Sigma12 = Sigma[unknown, known]
  Sigma22 = Sigma[known, known]
  Sigma221 = chol2inv(chol(Sigma22))
  
  mu_ = mu[unknown] + Sigma12%*%Sigma221%*%(as.vector(t(train))[known] - mu[known])
  
  mu_all = rep(0, length(mu))
  mu_all[unknown] = as.vector(mu_)
  mu_all[known] = as.vector(t(train))[known]
  
  preds = matrix(mu_all, nrow = n_drugs, byrow=T)
  
  return(preds)
  
}

compute_gradient = function(y, alpha, beta, adj_mat_train, B_, X, n_train){
  
  A = matrix(0, nrow = n_train, ncol = n_train)
  B = matrix(0, nrow = n_train, ncol = n_train)
  
  for (i in 1:(n_train)){
    for (j in 1:(n_train)){
      if (i==j){
        A[i,j] = sum(alpha)
        B[i,j] = beta*sum(adj_mat_train[i,]) 
      } else {
        B[i,j] = -(beta*adj_mat_train[i,j])
      }
    }
  }
  
  b = 2*X*alpha
  
  Sigma1 = 2*(A+B)
  require("MASS")
  #Sigma = ginv(Sigma1)
  inds = which(Sigma1!=0, arr.ind = T)
  sm = sparseMatrix(i = inds[,1], j = inds[,2], x = Sigma1[inds])
  Sigma = chol2inv(chol(sm))
  Sigma = as.matrix(Sigma)
  #temp = chol(Sigma1) 
  #Sigma = chol2inv(temp)
  
  mu = Sigma%*%b
  
  grad_alpha = -t(y)%*%y + 2*t(y)%*%X - 2*t(X)%*%mu + t(mu)%*%mu + sum(diag(Sigma))
  
  grad_beta = -t(y)%*%B_%*%y + t(mu)%*%B_%*%mu + t(as.vector(Sigma))%*%as.vector(B_)
  
  return(list(grad_alpha, grad_beta))
  
}

train_and_predict_single = function(mf_preds_train_mat, dt_mat, sim_mat, temp_inds, mf_preds_all, target_set, test_ind, adj_mat, crf_iters, dataset_triplets){
  
  crf_predictions = rep(NA, nrow(temp_inds))
  
  for (t in target_set){
    
    cat('target ',t,'... ',length(which(dt_mat[,t]>=0)),' observations\n')
    
    if (length(which(dt_mat[,t]>=0))>0){
      
      mf_pred_train_col_t = mf_preds_train_mat[which(!is.na(mf_preds_train_mat[,t])),t]
      adj_mat_train_col_t = make_training_adj_mat_for_column(dt_mat, sim_mat, t)
      training_vals_col_t = dt_mat[which(dt_mat[,t]>=0),t]
      
      if (length(which(dt_mat[,t]>=0))>400){
        eta = 0.008
      } else{
        eta = 0.01
      }
      
      mf_prediction_col = mf_preds_all[,t]
      params = train_crf_column(y = training_vals_col_t, X = mf_pred_train_col_t, adj_mat = adj_mat_train_col_t, crf_iters = crf_iters, eta = eta)  
      cat('learned parameters: ', params[[1]], params[[2]],'\n')
      cat('making predictions..\n')
      crf_prediction_col = make_crf_predictions_col(params[[1]], params[[2]], column = dt_mat[,t], adj_mat = adj_mat, X = mf_prediction_col)
      
      ## put the prediction into the vector
      inds = which(temp_inds[test_ind,2] == t)
      crf_prediction_test_col = crf_prediction_col[temp_inds[test_ind[inds],1]]
      crf_predictions[test_ind[inds]] = crf_prediction_test_col
      
      inds = which(!is.na(crf_predictions))
      crf_metrics = get_metrics(crf_predictions[inds], dataset_triplets[inds,3], 7.6)
      
    } else{
      inds = which(temp_inds[test_ind,2] == t)
      crf_predictions[test_ind[inds]] = mean(mf_preds_train_mat)
    }
    
    
    cat('auc so far, crf: ',round(crf_metrics[[2]], digits = 5),'\n')
    cat('aupr so far, crf: ',round(crf_metrics[[3]], digits = 5),'\n')
    
  }
  
  return (crf_predictions)
  
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