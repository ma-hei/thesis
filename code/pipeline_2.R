print_prediction_metrics_set = function(target_set, mix_dataset, test_ind, crf_preds, mf_preds_all_mat){
  
  crf_predictions = rep(NA, nrow(mix_dataset))
  mf_predictions = rep(NA, nrow(mix_dataset))
  for (t in target_set){
    inds = which(mix_dataset[test_ind,2]==t)
    drugs = mix_dataset[test_ind[inds],1]
    labels = mix_dataset[test_ind[inds],3]
    crf_predictions_test_col = crf_preds[drugs, match(t,target_set)]
    ans1 = get_metrics(preds = crf_predictions_test_col, labels = labels, cutoff = 7)
    mf_predictions_test_col = mf_preds_all_mat[drugs, match(t,target_set)]
    ans2 = get_metrics(preds = mf_predictions_test_col, labels = labels, cutoff = 7)
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

make_B_matrices = function(training_data_mat, row_adj_mat, col_adj_mat){
  
  n_drugs = nrow(training_data_mat)
  n_targets = ncol(training_data_mat)
  
  ## make Bk matrices.. one for each column
  Bk_matrices = list()
  for (t in 1:n_targets){
    obs = which(training_data_mat[,t]>=0)
    n_obs = length(obs)
    Bk_matrices[[t]] = matrix(0, nrow = n_obs, ncol = n_obs)
    for (o in obs){
      n_neighbors = 0
      for (o_ in setdiff(obs,o)){
        if (row_adj_mat[o,o_]==1){
          n_neighbors = n_neighbors + 1
          Bk_matrices[[t]][match(o,obs), match(o_,obs)] = -1 
        }
      }
      Bk_matrices[[t]][match(o,obs), match(o,obs)] = n_neighbors
    }
  }
  
  n_train = length(which(training_data_mat>0))
  
  B_row = matrix(0, nrow = n_train, ncol = n_train)
  
  inds = which(t(training_data_mat)>=0, arr.ind = T)
  for (t in 1:n_targets){
    neighbors = which(col_adj_mat[t,]>0)
    for (d in 1:n_drugs){
      if (training_data_mat[d,t]>=0){        
        cell = which(inds[,1]==t & inds[,2]==d)
        avail_targets = which(training_data_mat[d,]>=0)
        n_neighbors = length(which(neighbors%in%avail_targets))
        B_row[cell,cell] = n_neighbors
      }  
    }
  }
  
  for (t in 1:n_targets){
    neighbors = which(col_adj_mat[t,]>0)
    for (n in neighbors){
      for (d in 1:n_drugs){
        avail_targets = which(training_data_mat[d,]>0)
        if (t %in% avail_targets && n %in% avail_targets){
          if (d==1){
            if (t==1){
              cell_a = 1
            } else {
              cell_a =length(which(training_data_mat[1,1:(t-1)]>0)) + 1
            }
            if (n==1){
              cell_b = 1
            } else{
              cell_b =length(which(training_data_mat[1,1:(n-1)]>0)) + 1
            }
          } else {
            temp = length(which(training_data_mat[1:(d-1),]>0))
            if (t==1){
              cell_a = temp + 1
            } else {
              cell_a = temp + length(which(training_data_mat[d,1:(t-1)]>0)) + 1
            }
            if (n==1){
              cell_b = temp + 1
            } else {
              cell_b = temp + length(which(training_data_mat[d,1:(n-1)]>0)) + 1
            }
          }
          B_row[cell_a, cell_b] = -1
          B_row[cell_b, cell_a] = -1
        }
      }
    }
  }
  
  
  return (list(Bk_matrices, B_row))
  
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

get_Sigma_mu_train = function(training_data_mat, X_vec, target_beta_train_list, row_beta_train, alpha_train, row_adj_mat, col_adj_mat){
  
  n_train = length(which(training_data_mat>0))
  A_train = matrix(0, nrow = n_train, ncol = n_train)
  B_train = matrix(0, nrow = n_train, ncol = n_train)
  
  inds = which(t(training_data_mat)>=0, arr.ind = T)
  for (i in 1:nrow(inds)){
    drug = inds[i,2]
    target = inds[i,1]
    
    drug_neighbors = which(row_adj_mat[drug,]>0)
    n_drug_neighbors = 0
    for (d in drug_neighbors){
      if (training_data_mat[d, target]>=0){
        n_drug_neighbors = n_drug_neighbors+1
        cell_a = i
        cell_b = which(inds[,1]==target & inds[,2]==d)
        B_train[cell_a, cell_b] = -target_beta_train_list[[target]]
      }
    }
    
    target_neighbors = which(col_adj_mat[target,]>0)
    n_target_neighbors = 0
    for (t in target_neighbors){
      if (training_data_mat[drug, t]>=0){
        n_target_neighbors = n_target_neighbors+1
        cell_a = i
        cell_b = which(inds[,1]==t & inds[,2]==drug)
        B_train[cell_a, cell_b] = -row_beta_train
      }
    }
    B_train[i,i] = n_target_neighbors * row_beta_train + n_drug_neighbors*target_beta_train_list[[target]]
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

get_gaussian_predict_simple = function(training_data_mat, X_mat, alpha, beta_row, beta_train, B_row, B_col){
  
  n_drugs = nrow(training_data_mat)
  n_targets = ncol(training_data_mat)
  A = matrix(0, nrow = n_drugs*n_targets, ncol = n_drugs*n_targets)
  B = matrix(0, nrow = n_drugs*n_targets, ncol = n_drugs*n_targets)
  
  B = row_beta_train*B_row + col_beta_train*B_col
  for (i in 1:(n_drugs*n_targets)){
    A[i,i] = alpha
  }
  
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
    
    if ((it%%50)==0){
      cat('iteration ',it,'.. alpha ',alpha_train,' row beta ',row_beta_train,' col beta ',col_beta_train,'\n')
    }
    #cat('iteration ',it,'.. alpha ',alpha_train,' row beta ',row_beta_train,' col beta ',col_beta_train,'\n')
    
    ans = get_gaussian_train_simple(training_data_mat = training_data_mat, X_vec = X_train, alpha_train = alpha_train, row_beta_train = row_beta_train, col_beta_train = col_beta_train, B_row = B_row, B_col = B_col)
    Sigma = ans[[1]]
    mu = ans[[2]]
    
    log_alpha_train = log(alpha_train)
    grad_alpha_train = get_gradient_alpha(y_vec = y_vec, X_vec = X_train, mu = mu, Sigma = Sigma)
    grad_log_alpha_train = alpha_train*grad_alpha_train
    log_alpha_train = log_alpha_train + eta * grad_log_alpha_train
    alpha_train = exp(log_alpha_train)
   
    log_row_beta_train = log(row_beta_train)
    grad_row_beta_train = get_B_gradient(y_vec = y_vec, B = B_row, mu = mu, Sigma = Sigma)
    grad_log_row_beta_train = row_beta_train*grad_row_beta_train
    log_row_beta_train = log_row_beta_train + eta*grad_log_row_beta_train
    row_beta_train = exp(log_row_beta_train)
    
    log_col_beta_train = log(col_beta_train)
    grad_col_beta_train = get_B_gradient(y_vec = y_vec, B = B_col, mu = mu, Sigma = Sigma)
    grad_log_col_beta_train = col_beta_train*grad_col_beta_train
    log_col_beta_train = log_col_beta_train + eta*grad_log_col_beta_train
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
    
  }
  
  return (list(alpha_train, row_beta_train, col_beta_train))
  
}

simple_crf_predict = function(training_data_mat, row_adj_mat, col_adj_mat, alpha, row_beta, col_beta, X_mat){
  
  n_drugs = nrow(training_data_mat)
  n_targets = ncol(training_data_mat)
  
  #A = matrix(0, nrow = n_drugs*n_targets, ncol = n_drugs*n_targets)
  #B = matrix(0, nrow = n_drugs*n_targets, ncol = n_drugs*n_targets)
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
      #  B[temp_a,temp_b] = B[temp_a,temp_b] - col_beta
      #  triplets_B = rbind(triplets_B, c(temp_a, temp_b, -col_beta))
        triplets_B[counter_B,] = c(temp_a, temp_b, -col_beta)
        counter_B = counter_B+1
      }
      for (cn in col_neighbors){
        temp_b = (d-1)*n_targets+cn
      #  B[temp_a,temp_b] = B[temp_a,temp_b] - row_beta
      #  triplets_B = rbind(triplets_B, c(temp_a, temp_b, -row_beta))
        triplets_B[counter_B,] = c(temp_a, temp_b, -row_beta)
        counter_B = counter_B+1
      }
      #A[temp_a, temp_a] = alpha
      nrow_neighbors = length(which(row_adj_mat[d,]>0))
      ncol_neighbors = length(which(col_adj_mat[t,]>0))
      sum = col_beta*nrow_neighbors+row_beta*ncol_neighbors
      #B[temp_a, temp_a] = sum
      #triplets_B = rbind(triplets_B, c(temp_a, temp_a, sum))
      #triplets_A = rbind(triplets_A, c(temp_a, temp_a, alpha))
      triplets_A[counter_A,] = c(temp_a, temp_a, alpha)
      triplets_B[counter_B,] = c(temp_a, temp_a, sum)
      counter_B = counter_B+1
      counter_A = counter_A+1
    }
  }
  A = sparseMatrix(i = triplets_A[,1], j = triplets_A[,2], x = triplets_A[,3])
  B = sparseMatrix(i = triplets_B[,1], j = triplets_B[,2], x = triplets_B[,3])
  
  Sigma1 = 2*(A+B)
  #require("MASS")
  #temp = chol(Sigma1) 
  #Sigma = chol2inv(temp)
  inds = which(Sigma1!=0, arr.ind = T)
  sm = sparseMatrix(i = inds[,1], j = inds[,2], x = Sigma1[inds])
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

get_B_gradient = function(y_vec, B, mu, Sigma){
  
  gradient_beta = -t(y_vec)%*%B%*%y_vec+t(mu)%*%B%*%mu+t(as.vector(Sigma)%*%as.vector(B))
  
  return(as.numeric(gradient_beta))
}
  
get_gradient_alpha = function(y_vec, X_vec, mu, Sigma){
  
  gradient_alpha = -t(y_vec)%*%y_vec+2*t(y_vec)%*%X_vec-2*t(X_vec)%*%mu+t(mu)%*%mu+sum(diag(Sigma))
  
  return(as.numeric(gradient_alpha))
  
}

get_gradient_target_t = function(y_vec_target, B_mat_t, Sigma_t, mu_t){
  
  gradient_beta_t = -t(y_vec_target)%*%B_mat_t%*%y_vec_target+t(mu_t)%*%B_mat_t%*%mu_t+t(as.vector(Sigma_t))%*%as.vector(B_mat_t)
  
  return (gradient_beta_t)
  
}

get_gradient_rows = function(y_vec, B_row, mu, Sigma){
  
  gradient_beta_rows = -t(y_vec)%*%B_row%*%y_vec+t(mu)%*%B_row%*%mu+t(as.vector(Sigma)%*%as.vector(B_row))
  
  return(gradient_beta_rows)
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

train_crf_2 = function(training_data_mat, row_adj_mat, col_adj_mat, eta, crf_iters, mf_preds){
  
  n_targets = ncol(training_data_mat)
  
  ans = make_B_matrices(training_data_mat = training_data_mat, row_adj_mat = row_adj_mat, col_adj_mat = col_adj_mat)
  Bcol_mats = ans[[1]]
  Brow_mat = ans[[2]]
  
  #mf_train = get_mf_cv(training_data_mat, iters = 400)
  #X_train = t(mf_train[[2]])[which(!is.na(t(mf_train[[2]])))]
  y_vec = t(training_data_mat)[which(t(training_data_mat)>=0)]
  X_train = t(mf_preds)[which(!is.na(t(mf_preds)))]
  
  target_beta_train_list = list()
  for (i in 1:n_targets){
    #target_beta_train_list[[i]] = target_beta[[i]]
    target_beta_train_list[[i]] = 0.001
  }
  #row_beta_train = row_beta
  row_beta_train = 0.001
  #alpha_train = alpha
  alpha_train = 0.001
  
  eta = eta
  
  for (i in 1:crf_iters){
    #cat('iter ',i,'.. alpha: ',alpha_train,' row beta: ',row_beta_train,'\n')
    #cat('iter ',i,'.. col betas: ' , target_beta_train_list[[1]],' ',target_beta_train_list[[2]],' ',target_beta_train_list[[3]],' ',target_beta_train_list[[4]],'\n')
    #cat('iter ',i,'.. col betas: ' , target_beta_train_list[[1]],' ',target_beta_train_list[[2]],'\n')
    ans = get_Sigma_mu_train(training_data_mat = training_data_mat, X_vec = X_train, target_beta_train_list = target_beta_train_list, row_beta_train = row_beta_train, alpha_train = alpha_train, row_adj_mat = row_adj_mat, col_adj_mat = col_adj_mat)
    Sigma = ans[[1]]
    mu = ans[[2]]
    
    log_alpha_train = log(alpha_train)
    grad_alpha_train = get_gradient_alpha(y_vec = y_vec, X_vec = X_train, mu = mu, Sigma = Sigma)
    grad_log_alpha_train = alpha_train*grad_alpha_train
    log_alpha_train = log_alpha_train + eta * grad_log_alpha_train
    alpha_train = exp(log_alpha_train)
    if (alpha_train<1e-10){
      alpha_train = 1e-10
    } else if (alpha_train>2){
      alpha_train = 2
    }
    
    for (t in 1:n_targets){
      inds = which(t(training_data_mat)>=0, arr.ind = T)
      inds_ = which(inds[,1]==t)
      Sigma_t = Sigma[inds_, inds_]
      mu_t = mu[inds_]
      y_vec_t = training_data_mat[which(training_data_mat[,t]>=0),t]
      log_target_beta_train = log(target_beta_train_list[[t]])
      grad_target_beta_train = get_gradient_target_t(y_vec_target = y_vec_t, B_mat_t = Bcol_mats[[t]], Sigma_t = Sigma_t, mu_t = mu_t)
      grad_log_target_beta_train = target_beta_train_list[[t]]*grad_target_beta_train
      if (length(y_vec_t)<100){
        eta_ = 0.07
      } else{
        eta_ = 0.01
      }
      log_target_beta_train = log_target_beta_train + eta_ *grad_log_target_beta_train
      target_beta_train_list[[t]] = exp(log_target_beta_train)
      if (target_beta_train_list[[t]]<1e-10){
        target_beta_train_list[[t]] = 1e-10
      } else if (target_beta_train_list[[t]]>2){
        target_beta_train_list[[t]] = 2
      }
    }
    
    log_row_beta_train = log(row_beta_train)
    grad_row_beta_train = get_gradient_rows(y_vec = y_vec, B_row = Brow_mat, mu = mu, Sigma = Sigma)
    grad_log_row_beta_train = row_beta_train*grad_row_beta_train
    log_row_beta_train = log_row_beta_train + eta*grad_log_row_beta_train
    row_beta_train = exp(log_row_beta_train)
    if (row_beta_train<1e-10){
      row_beta_train = 1e-10
    } else if (row_beta_train>2){
      row_beta_train = 2
    }
  }
  
  return(list(alpha_train, row_beta_train, target_beta_train_list))
}

make_prediction_2 = function(training_data_mat, target_beta, row_beta, alpha, mf_preds_columns, adj_mat_cols, adj_mat_rows){
  
  n_targets = ncol(training_data_mat)
  A = matrix(0, nrow = n_drugs*n_targets, ncol = n_drugs*n_targets)
  B = matrix(0, nrow = n_drugs*n_targets, ncol = n_drugs*n_targets)
  
  for (d in 1:n_drugs){
    for (t in 1:n_targets){
      row_neighbors = which(adj_mat_rows[d,]==1)
      col_neighbors = which(adj_mat_cols[t,]==1)
      temp_a = (d-1)*n_targets+t
      for (rn in row_neighbors){
        temp_b = (rn-1)*n_targets+t
        B[temp_a,temp_b] = B[temp_a,temp_b] - target_beta[[t]]
      }
      for (cn in col_neighbors){
        temp_b = (d-1)*n_targets+cn
        B[temp_a,temp_b] = B[temp_a,temp_b] - row_beta
      }
      A[temp_a, temp_a] = sum(alpha)
      nrow_neighbors = length(which(adj_mat_rows[d,]>0))
      ncol_neighbors = length(which(adj_mat_cols[t,]>0))
      sum = target_beta[[t]]*nrow_neighbors+row_beta*ncol_neighbors
      B[temp_a, temp_a] = sum
    }
  }
  
  Sigma1 = 2*(A+B)
  require("MASS")
  temp = chol(Sigma1) 
  Sigma = chol2inv(temp)
  
  X = as.vector(t(mf_preds_columns))
  b = 2*X*alpha
  
  mu = Sigma%*%b
  
  unknown = which(as.vector(t(dt_mat_temp))<0)
  known = which(as.vector(t(dt_mat_temp))>=0)
  
  Sigma12 = Sigma[unknown, known]
  Sigma22 = Sigma[known, known]
  Sigma221 = chol2inv(chol(Sigma22))
  mu_ = mu[unknown] + Sigma12%*%Sigma221%*%(as.vector(t(dt_mat_temp))[known] - mu[known])
  mu_all = rep(0, length(mu))
  mu_all[unknown] = mu_
  mu_all[known] = as.vector(t(dt_mat_temp))[known]
  
  pred_mat = matrix(mu_all, nrow = n_drugs, byrow = T)
  
  return (pred_mat)
  
}

train_and_predict_2 = function(target_adj_mat, drug_adj_mat, dt_mat, target_set, mf_preds_train_mat, mf_preds_all_mat, test_ind, mix_dataset, crf_iters){
  
  crf_predictions = rep(NA, nrow(mix_dataset))
  mf_predictions = rep(NA, nrow(mix_dataset))
  
  adj_mat_cols = target_adj_mat[target_set, target_set]
  adj_mat_rows = drug_adj_mat
  
  dt_mat_temp = dt_mat[,target_set]
  
  mf_preds = mf_preds_train_mat[,target_set]
  
  cat('learning params..\n')
  
  params = train_crf_2(training_data_mat = dt_mat_temp, row_adj_mat = adj_mat_rows, col_adj_mat = adj_mat_cols, eta = 0.01, crf_iters = crf_iters, mf_preds = mf_preds)
  
  cat('learned params.. ',params[[1]],'.. ',params[[2]],'..\n')
  cat('making prediction..\n')
  
  pred_mat = make_prediction_2(training_data_mat = dt_mat_temp, target_beta = params[[3]], row_beta = params[[2]], alpha = params[[1]], mf_preds_columns = mf_preds_all_mat[,target_set], adj_mat_cols = adj_mat_cols, adj_mat_rows = adj_mat_rows)
  
  for (t in target_set){
    inds = which(mix_dataset[test_ind,2]==t)
    drugs = mix_dataset[test_ind[inds],1]
    labels = mix_dataset[test_ind[inds],3]
    crf_predictions_test_col = pred_mat[drugs, match(t,target_set)]
    ans1 = get_metrics(preds = crf_predictions_test_col, labels = labels, cutoff = 7)
    mf_predictions_test_col = mf_preds_all_mat[drugs, t]
    ans2 = get_metrics(preds = mf_predictions_test_col, labels = labels, cutoff = 7)
    cat('target ',t,' ','RMSE mf ',ans2[[1]],' crf ',ans1[[1]],'\n')
    cat('target ',t,' ','AUC mf ',ans2[[2]],' crf ',ans1[[2]],'\n')
    cat('target ',t,' ','AUPR mf ',ans2[[3]],' crf ',ans1[[3]],'\n\n')
    crf_predictions[test_ind[inds]] = crf_predictions_test_col
    mf_predictions[test_ind[inds]] = mf_predictions_test_col
  }
  inds = which(!is.na(mf_predictions))
  ans2 = get_metrics(mf_predictions[inds], mix_dataset[inds,3], 7)
  ans1 = get_metrics(crf_predictions[inds], mix_dataset[inds,3], 7)
  cat('all RMSE mf ',ans2[[1]],' crf ',ans1[[1]],'\n')
  cat('all AUC mf ',ans2[[2]],' crf ',ans1[[2]],'\n')
  cat('all AUPR mf ',ans2[[3]],' crf ',ans1[[3]],'\n\n')
  
}

train_and_predict_single = function(mf_preds_train_mat, dt_mat, sim_mat, mix_dataset, mf_preds_all, target_set, test_ind, adj_mat, crf_iters){
  
  
  crf_predictions = rep(NA, nrow(mix_dataset))
  mf_predictions = rep(NA, nrow(mix_dataset))
  
  for (t in target_set){
    
    cat('fold ',i,', target ',t,'... ',length(which(dt_mat[,t]>=0)),' observations\n')
    
    mf_pred_train_col_t = mf_preds_train_mat[which(!is.na(mf_preds_train_mat[,t])),t]
    adj_mat_train_col_t = make_training_adj_mat_for_column(dt_mat, sim_mat, t)
    training_vals_col_t = dt_mat[which(dt_mat[,t]>=0),t]
    
    if (length(which(dt_mat[,t]>=0))>500){
      eta = 0.001
    } else{
      eta = 0.01
    }
    
    params = train_crf_row(y = training_vals_col_t, X = mf_pred_train_col_t, adj_mat = adj_mat_train_col_t, crf_iters = crf_iters, eta = eta)  
    cat('learned parameters: ', params[[1]], params[[2]],'\n')
    
    inds = which(mix_dataset[test_ind,2] == t)
    labels_test_col = mix_dataset[test_ind[inds], 3]
    mf_prediction_col = mf_preds_all[,t]
    mf_prediction_test_col = mf_preds_all[cbind(mix_dataset[test_ind[inds],1],mix_dataset[test_ind[inds],2])]
    mf_predictions[test_ind[inds]] = mf_prediction_test_col
    
    #cat('making crf predictions..\n')
    
    crf_prediction_col = make_crf_predictions_row(params[[1]], params[[2]], column = dt_mat[,t], adj_mat = adj_mat, X = mf_prediction_col)
    crf_prediction_test_col = crf_prediction_col[mix_dataset[test_ind[inds],1]]
    crf_predictions[test_ind[inds]] = crf_prediction_test_col
    
    mf_metrics = get_metrics(mf_prediction_test_col, labels_test_col, 7)
    crf_metrics = get_metrics(crf_prediction_test_col, labels_test_col, 7)
    
    #cat('target rmse (mf, crf): ',mf_metrics[[1]],', ',crf_metrics[[1]],'\n')
    #cat('target auc (mf, crf): ',mf_metrics[[2]],', ',crf_metrics[[2]],'\n')
    #cat('target aupr (mf, crf): ',mf_metrics[[3]],', ',crf_metrics[[3]],'\n')
    
  }
  
  inds = which(!is.na(mf_predictions))
  mf_metrics = get_metrics(mf_predictions[inds], mix_dataset[inds,3], 7)
  crf_metrics = get_metrics(crf_predictions[inds], mix_dataset[inds,3], 7)
  
  cat('all test rmse (mf, crf) so far: ',round(mf_metrics[[1]], digits = 3),', ',round(crf_metrics[[1]], digits = 3),'\n')
  cat('all test auc (mf, crf) so far: ',round(mf_metrics[[2]], digits = 3),', ',round(crf_metrics[[2]], digits = 3),'\n')
  cat('all test aupr (mf, crf) so far: ',round(mf_metrics[[3]], digits = 3),', ',round(crf_metrics[[3]], digits = 3),'\n\n')
  
  return (list(crf_predictions, mf_predictions))
  
}