
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

get_Sigma_mu_train = function(training_data_mat, X_vec, target_beta_train_list, row_beta_train, alpha_train){
  
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

get_gradient_alpha = function(y_vec, X_vec, mu, Sigma){
  
  gradient_alpha = -t(y_vec)%*%y_vec+2*t(y_vec)%*%X_vec-2*t(X_vec)%*%mu+t(mu)%*%mu+sum(diag(Sigma))
  
  return(gradient_alpha)
  
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

train_crf = function(training_data_mat, row_adj_mat, col_adj_mat, eta, crf_iters, mf_preds){
  
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
    cat('iter ',i,'.. alpha: ',alpha_train,' row beta: ',row_beta_train,'\n')
    cat('iter ',i,'.. col betas: ' , target_beta_train_list[[1]],' ',target_beta_train_list[[2]],' ',target_beta_train_list[[3]],' ',target_beta_train_list[[4]],'\n')
    ans = get_Sigma_mu_train(training_data_mat = training_data_mat, X_vec = X_train, target_beta_train_list = target_beta_train_list, row_beta_train = row_beta_train, alpha_train = alpha_train)
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

make_prediction = function(params, X){
  
}