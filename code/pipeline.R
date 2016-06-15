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
  temp = chol(Sigma1) 
  Sigma = chol2inv(temp)
  
  mu = Sigma%*%b
  
  grad_alpha = -t(y)%*%y + 2*t(y)%*%X - 2*t(X)%*%mu + t(mu)%*%mu + sum(diag(Sigma))
  
  grad_beta = -t(y)%*%B_%*%y + t(mu)%*%B_%*%mu + t(as.vector(Sigma))%*%as.vector(B_)
  
  return(list(grad_alpha, grad_beta))
  
}

train_crf = function(train, adj, mf_iters, crf_iters){
  
  n_drugs = nrow(train)
  n_targets = ncol(train)
  
  n_train = length(which(train>0))
  
  #ind = which(!is.na(train), arr.ind = T)
  ind = which(train<0,arr.ind=T)
  cut_out = rep(0, length(which(train<0)))
  for (i in 1:length(which(train<0))){
    row = ind[i,1]
    col = ind[i,2]
    mapped = (row-1)*n_targets+col
    cut_out[i] = mapped
  }
  
  adj_train = adj[-cut_out, -cut_out]
  
  #train_triplet = matrix(0, nrow = length(which(!is.na(training_data))), ncol = 3)
  #train_triplet[,1] = which(!is.na(training_data), arr.ind = T)[,1]
  #train_triplet[,2] = which(!is.na(training_data), arr.ind = T)[,2]
  #train_triplet[,3] = which(!is.na(training_data), arr.ind = T)[,3]
  train_triplet = matrix(0, nrow = length(which(train>0)), ncol = 3)
  train_triplet[,1] = which(train>0, arr.ind = T)[,1]
  train_triplet[,2] = which(train>0, arr.ind = T)[,2]
  train_triplet[,3] = train[which(train>0, arr.ind = T)]
  
  res = nmf.cv(train_triplet[,1:3], nfolds = 5, m = n_drugs, n = n_targets, k = 10,
            lambda = 1, gamma = 1, tol = 0.001, maxiter = mf_iters, threshold = 2, bias=TRUE,
            interaction = TRUE, prediction = TRUE)
  
  temp = matrix(NA, nrow = nrow(train), ncol = ncol(train))
  for (i in 1:nrow(train_triplet)){
    temp[train_triplet[i,1], train_triplet[i,2]] = res[[3]][i]
  }
  
  y = t(train)[which(t(train>0), arr.ind = T)]
  X = t(temp)[which(t(!is.na(temp)), arr.ind = T)]
  
  B = matrix(0, nrow = n_train, ncol = n_train)
  for (i in 1:n_train){
    for (j in 1:n_train){
      if (i==j){
        B[i,j] = sum(adj_train[i,])
      } else {
        B[i,j] = -adj_train[i,j]
      }
    }
  }
  
  alpha = 0.001
  beta = 0.001
  
  log_alpha = log(alpha)
  log_beta = log(beta)
  
  eta = 0.01
  
  for (i in 1:crf_iters){
    
    cat('training iteration ',i,':',alpha,' ',beta,'\n')
    gradients = compute_gradient(y,alpha, beta, adj_train, B, X, n_train)
    grad_alpha = gradients[[1]]
    grad_beta = gradients[[2]]
    
    grad_log_alpha = alpha*grad_alpha
    grad_log_beta = beta*grad_beta
    
    log_alpha = log_alpha + eta * grad_log_alpha
    log_beta = log_beta + eta * grad_log_beta
    
    alpha = exp(log_alpha)
    beta = exp(log_beta)
    
  }
  
  return(list(alpha, beta))
  
}

train_crf_by_row = function(train, adj, mf_iters, crf_iters){
  
  train_triplet = matrix(0, nrow = length(which(train>0)), ncol = 3)
  train_triplet[,1] = which(train>0, arr.ind = T)[,1]
  train_triplet[,2] = which(train>0, arr.ind = T)[,2]
  train_triplet[,3] = train[which(train>0, arr.ind = T)]
  
  res = nmf.cv(train_triplet[,1:3], nfolds = 5, m = n_drugs, n = n_targets, k = 10,
               lambda = 1, gamma = 1, tol = 0.001, maxiter = 100, threshold = 2, bias=TRUE,
               interaction = TRUE, prediction = TRUE)
  
  temp = matrix(NA, nrow = nrow(train), ncol = ncol(train))
  for (i in 1:nrow(train_triplet)){
    temp[train_triplet[i,1], train_triplet[i,2]] = res[[3]][i]
  }
  
  alpha = 0.01
  beta = 0.01
  
  log_alpha = log(alpha)
  log_beta = log(beta)
  
  eta = 0.001
  
  for (it in 1:crf_iters){
    
    cat('training iteration ',it,':',alpha,' ',beta,'\n')
    
    ## randomly choose a row
    row = sample(1:n_drugs, 1)
    ## train on row
    adj_train = adj_mat[(row-1)*n_targets + which(train[row,]>0),(row-1)*n_targets + which(train[row,]>0)]
    y = train[row, which(train[row,]>0)]
    X = temp[row, which(!is.na(temp[row,]))]
    
    n_train = length(which(train[row,]>0))
    
    B1 = matrix(0, nrow = n_train, ncol = n_train)
    for (i in 1:n_train){
      for (j in 1:n_train){
        if (i==j){
          B1[i,j] = sum(adj_train[i,])
        } else {
          B1[i,j] = -adj_train[i,j]
        }
      }
    }
    
    gradients = compute_gradient(y,alpha, beta, adj_train, B1, X, n_train)
    grad_alpha = gradients[[1]]
    grad_beta = gradients[[2]]
    
    grad_log_alpha = alpha*grad_alpha
    grad_log_beta = beta*grad_beta
    
    log_alpha = log_alpha + eta * grad_log_alpha
    log_beta = log_beta + eta * grad_log_beta
    
    alpha = exp(log_alpha)
    beta = exp(log_beta)
    
  }
  
  return(list(alpha,beta))
}

## in case I want to compute the X beforehand
train_crf_by_row_1 = function(train, adj, X, crf_iters){
  
  train_triplet = matrix(0, nrow = length(which(train>0)), ncol = 3)
  train_triplet[,1] = which(train>0, arr.ind = T)[,1]
  train_triplet[,2] = which(train>0, arr.ind = T)[,2]
  train_triplet[,3] = train[which(train>0, arr.ind = T)]
  
  temp = matrix(NA, nrow = nrow(train), ncol = ncol(train))
  for (i in 1:nrow(train_triplet)){
    temp[train_triplet[i,1], train_triplet[i,2]] = X[i]
  }
  
  alpha = 0.01
  beta = 0.01
  
  log_alpha = log(alpha)
  log_beta = log(beta)
  
  eta = 0.001
  
  for (it in 1:crf_iters){
    
    cat('training iteration ',it,':',alpha,' ',beta,'\n')
    
    ## randomly choose a row
    row = sample(1:n_drugs, 1)
    ## train on row
    adj_train = adj_mat[(row-1)*n_targets + which(train[row,]>0),(row-1)*n_targets + which(train[row,]>0)]
    y = train[row, which(train[row,]>0)]
    X = temp[row, which(!is.na(temp[row,]))]
    
    n_train = length(which(train[row,]>0))
    
    B1 = matrix(0, nrow = n_train, ncol = n_train)
    for (i in 1:n_train){
      for (j in 1:n_train){
        if (i==j){
          B1[i,j] = sum(adj_train[i,])
        } else {
          B1[i,j] = -adj_train[i,j]
        }
      }
    }
    
    gradients = compute_gradient(y,alpha, beta, adj_train, B1, X, n_train)
    grad_alpha = gradients[[1]]
    grad_beta = gradients[[2]]
    
    grad_log_alpha = alpha*grad_alpha
    grad_log_beta = beta*grad_beta
    
    log_alpha = log_alpha + eta * grad_log_alpha
    log_beta = log_beta + eta * grad_log_beta
    
    alpha = exp(log_alpha)
    beta = exp(log_beta)
    
  }
  
  return(list(alpha,beta))
}

## in case I can't store the adjacency mat
train_crf_by_row_2 = function(train, adj_list, X, crf_iters, eta){
  
  n_drugs = nrow(train)
  n_targets = ncol(train)
  
  train_triplet = matrix(0, nrow = length(which(train>0)), ncol = 3)
  train_triplet[,1] = which(train>0, arr.ind = T)[,1]
  train_triplet[,2] = which(train>0, arr.ind = T)[,2]
  train_triplet[,3] = train[which(train>0, arr.ind = T)]
  
  temp = matrix(NA, nrow = nrow(train), ncol = ncol(train))
  for (i in 1:nrow(train_triplet)){
    temp[train_triplet[i,1], train_triplet[i,2]] = X[i]
  }
  
  alpha = 0.01
  beta = 0.01
  
  log_alpha = log(alpha)
  log_beta = log(beta)
  
  eta = eta
  
  for (it in 1:crf_iters){
    
    cat('training iteration ',it,':',alpha,' ',beta,'\n')
    
    ## randomly choose a row
    col = sample(1:n_targets, 1)
    #col = 31
    ## train on row
    ## adj_train = adj_mat[(row-1)*n_targets + which(train[row,]>0),(row-1)*n_targets + which(train[row,]>0)]
    adj_train = adj_list[[col]]
    # y = train[row, which(train[row,]>0)]
    y = train[which(train[,col]>0),col]
    # X = temp[row, which(!is.na(temp[row,]))]
    X = temp[which(!is.na(temp[,col])),col]
    
    # n_train = length(which(train[row,]>0))
    n_train = length(which(train[,col]>0))
    
      cat('training on column ', col,'.. ',n_train,' observations\n')
      B1 = matrix(0, nrow = n_train, ncol = n_train)
      for (i in 1:n_train){
        for (j in 1:n_train){
          if (i==j){
            B1[i,j] = sum(adj_train[i,])
          } else {
            B1[i,j] = -adj_train[i,j]
          }
        }
      }
      
      gradients = compute_gradient(y,alpha, beta, adj_train, B1, X, n_train)
      grad_alpha = gradients[[1]]
      grad_beta = gradients[[2]]
      
      grad_log_alpha = alpha*grad_alpha
      grad_log_beta = beta*grad_beta
      
      log_alpha = log_alpha + eta * grad_log_alpha
      log_beta = log_beta + eta * grad_log_beta
      
      alpha = exp(log_alpha)
      beta = exp(log_beta)
      
  }
  
  return(list(alpha,beta))
  
}

train_crf_by_row_3 = function(train, adj_list, X, crf_iters, eta, col){
  
  n_drugs = nrow(train)
  n_targets = ncol(train)
  
  train_triplet = matrix(0, nrow = length(which(train>0)), ncol = 3)
  train_triplet[,1] = which(train>0, arr.ind = T)[,1]
  train_triplet[,2] = which(train>0, arr.ind = T)[,2]
  train_triplet[,3] = train[which(train>0, arr.ind = T)]
  
  temp = matrix(NA, nrow = nrow(train), ncol = ncol(train))
  for (i in 1:nrow(train_triplet)){
    temp[train_triplet[i,1], train_triplet[i,2]] = X[i]
  }
  
  alpha = 0.01
  beta = 0.01
  
  log_alpha = log(alpha)
  log_beta = log(beta)
  
  eta = eta
  
  for (it in 1:crf_iters){
    
    cat('training iteration ',it,':',alpha,' ',beta,'\n')
    
    ## randomly choose a row
    #col = sample(1:n_targets, 1)
    col = col
    ## train on row
    ## adj_train = adj_mat[(row-1)*n_targets + which(train[row,]>0),(row-1)*n_targets + which(train[row,]>0)]
    adj_train = adj_list[[col]]
    # y = train[row, which(train[row,]>0)]
    y = train[which(train[,col]>0),col]
    # X = temp[row, which(!is.na(temp[row,]))]
    X = temp[which(!is.na(temp[,col])),col]
    
    # n_train = length(which(train[row,]>0))
    n_train = length(which(train[,col]>0))
    
    #cat('training on column ', col,'.. ',n_train,' observations\n')
    B1 = matrix(0, nrow = n_train, ncol = n_train)
    for (i in 1:n_train){
      for (j in 1:n_train){
        if (i==j){
          B1[i,j] = sum(adj_train[i,])
        } else {
          B1[i,j] = -adj_train[i,j]
        }
      }
    }
    
    gradients = compute_gradient(y,alpha, beta, adj_train, B1, X, n_train)
    grad_alpha = gradients[[1]]
    grad_beta = gradients[[2]]
    
    grad_log_alpha = alpha*grad_alpha
    grad_log_beta = beta*grad_beta
    
    log_alpha = log_alpha + eta * grad_log_alpha
    log_beta = log_beta + eta * grad_log_beta
    
    alpha = exp(log_alpha)
    beta = exp(log_beta)
    
  }
  
  return(list(alpha,beta))
}

train_crf_row = function(y, X, adj_mat, crf_iters, eta){
  
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
    
    cat('training iteration ',it,':',alpha,' ',beta,'\n')
    if ((it%%40)==0){
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
      return(list(alpha,beta))
    }
    
  }
  return(list(alpha,beta))
}

generate_latent_factor_mat = function(n_drugs, n_targets){
  
  k = 20
  lf_drugs = matrix(0, nrow = n_drugs, k)
  for (i in 1:n_drugs){
    latent_factors_temp = rnorm(k, 0, 1)
    ## make sure all latent factors > 0 
    if (min(latent_factors_temp)<0){
      latent_factors_temp = latent_factors_temp - min(latent_factors_temp)
    }
    lf_drugs[i,] = latent_factors_temp
  }
  
  lf_targets = matrix(0, nrow = n_targets, k)
  for (i in 1:n_targets){
    latent_factors_temp = rnorm(k, 0, 1)
    ## make sure all latent factors > 0 
    if (min(latent_factors_temp)<0){
      latent_factors_temp = latent_factors_temp - min(latent_factors_temp)
    }
    lf_targets[i,] = latent_factors_temp
  }
  
  ## now that I have latent factors for the drugs and targets, multiply them to generate
  ## a matrix with values that have underlying latent factors
  training_data_lf = lf_drugs %*% t(lf_targets)
  
  return (training_data_lf)
  
}

create_example_adj_mat = function(n_drugs, n_targets){
  
  clique_list  = list()
  for (i in 1:(n_targets/10)){
    clique_list[[i]] = c(((i-1)*10+1):(i*10))
  }
  
  adj_mat = matrix(0, nrow = n_drugs*n_targets, ncol = n_drugs*n_targets)
  ## now go over each clique and for each clique member a, connect it to all other clique members
  for (i in 1:(length(clique_list))){
    for (k in 1:length(clique_list[[i]])){
      for (j in 1:length(clique_list[[i]])){
        if (k!=j){
          ##for all rows, connect rows clique_list[[i]][k] and clique_list[[i]][j]
          col_a = clique_list[[i]][k]
          col_b = clique_list[[i]][j]
          for (r in 1:n_drugs){
            ## connect cell at position [row r, col a] with cell at [row r, col b]
            adj_mat[(r-1)*n_targets+col_a, (r-1)*n_targets+col_b] = 1
            adj_mat[(r-1)*n_targets+col_b, (r-1)*n_targets+col_a] = 1
          }
        }
      }
    }
  }
  
  return(adj_mat)
}

create_crf = function(alpha, beta, X, adj_mat){
  
  n_cell = nrow(adj_mat)
  A = matrix(0, nrow = n_cell, ncol = n_cell)
  B = matrix(0, nrow = n_cell, ncol = n_cell)
  
  for (i in 1:n_cell){
    for (j in 1:n_cell){
      if (i==j){
        A[i,j] = sum(alpha)
        B[i,j] = beta*sum(adj_mat[i,])
      } else {
        B[i,j] = -beta*adj_mat[i,j]
      }
    }
  }
  b = 2*X*alpha
  
  Sigma1 = 2*(A+B)
  require("MASS")
  #Sigma = ginv(Sigma1)
  temp = chol(Sigma1) 
  Sigma = chol2inv(temp)
  mu = Sigma%*%b
  
  return(list(Sigma,mu))
  
}

create_crf_2 = function(alpha, Sigma, X){
  
  b = 2*X*alpha
  
  mu = Sigma%*%b
  
  return(list(Sigma,mu))
  
}

sample_training_data = function(all_data){
  
  clique_list  = list()
  for (i in 1:(ncol(all_data)/10)){
    clique_list[[i]] = c(((i-1)*10+1):(i*10))
  }
  
  training_data = matrix(-1, nrow = nrow(all_data), ncol = ncol(all_data))
  
  ## go over each row, in each row select one group for which no training data will
  ## be availble, for all other groups, select 2 - 5 cells that go into training data
  
  for (i in 1:nrow(all_data)){
    empty_group = sample(1:length(clique_list),3)
    for (k in 1:length(clique_list)){
      if (!k%in%empty_group){
        select_size = sample(1:2,1)
        selection = sample(1:length(clique_list[[k]]),select_size)
        for (x in 1:select_size){
          training_data[i,clique_list[[k]][selection[x]]] = all_data[i,clique_list[[k]][selection[x]]]
        }
      }
    }
  }
  return (training_data)
}

sample_training_data_sparser = function(all_data){
  
  
  clique_list  = list()
  for (i in 1:(ncol(all_data)/10)){
    clique_list[[i]] = c(((i-1)*10+1):(i*10))
  }
  
  training_data = matrix(-1, nrow = nrow(all_data), ncol = ncol(all_data))
  
  ## go over each row, in each row select one group for which no training data will
  ## be availble, for all other groups, select 2 - 5 cells that go into training data
  
  for (i in 1:nrow(all_data)){
    empty_group = sample(1:length(clique_list),6)
    for (k in 1:length(clique_list)){
      if (!k%in%empty_group){
        select_size = sample(1:2,1)
        selection = sample(1:length(clique_list[[k]]),select_size)
        for (x in 1:select_size){
          training_data[i,clique_list[[k]][selection[x]]] = all_data[i,clique_list[[k]][selection[x]]]
        }
      }
    }
  }
  return (training_data)
  
}

get_mf_prediction = function(train, iter){
  
  n_drugs = nrow(train)
  n_targets = ncol(train)
  
  train_mat_triplet = matrix(0, nrow = length(which(train>0)), ncol = 3)
  train_mat_triplet[,1] = which(train>0, arr.ind = T)[,1]
  train_mat_triplet[,2] = which(train>0, arr.ind = T)[,2]
  train_mat_triplet[,3] = train[which(train>0, arr.ind = T)]
  
  res = nmf(train_mat_triplet[,1:3], m = n_drugs, n = n_targets, k = 10,
            lambda = 1, gamma = 1, tol = 0.001, maxiter = iter, threshold = 2, bias=TRUE, interaction = TRUE)
  
  pred_mf = nmf.predict(res,cbind(rep(1:n_drugs, 1, each = n_targets), rep(1:n_targets, n_drugs)),
                        bias = TRUE,interaction = TRUE)
  
  pred_mf_mat = matrix(pred_mf, nrow = n_drugs, byrow = T)
  
  X = as.vector(t(pred_mf_mat))
  
  return(X)
  
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

get_crf_prediction = function(Sigma, mu, train){
  
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
  
  mat = matrix(mu_all, nrow = n_drugs, byrow=T)
  
  return(mat)
  
}

generate_dataset = function(n_drugs, n_targets){
  
## generate a matrix of values with underlying latent factors
  lf_mat = generate_latent_factor_mat(n_drugs, n_targets)
  myImagePlot(lf_mat)
  
  ## create an example adjacency matrix 
  adj_mat = create_example_adj_mat(n_drugs, n_targets)
  
  ## get mu and Sigma of probability distr. thats defined by a crf
  ## where X is the latent factor matrix
  res = create_crf(1,2,as.vector(t(lf_mat)), adj_mat)
  
  ## take a sample of this matrix (make sure all values above 0)
  m =  mvrnorm(n = 1, res[[2]], res[[1]], tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  m = (m-min(m))+1
  mat = matrix(m, nrow = n_drugs, byrow=T)
  ##myImagePlot(mat)
  
  train = sample_training_data(mat)
  
  return(list(mat, train, lf_mat))
  
}

generate_dataset_scenario_1 = function(n_drugs, n_targets){
  
  ## generate a matrix of values with underlying latent factors
  lf_mat = generate_latent_factor_mat(n_drugs, n_targets)
  ##myImagePlot(lf_mat)
  
  ## create an example adjacency matrix 
  adj_mat = create_example_adj_mat(n_drugs, n_targets)
  
  ## and remove edges of last clique
  for (i in 1:n_drugs){
    for (k in 61:70){
      mapped = (i-1)*n_targets+k
      adj_mat[mapped,] = 0
      adj_mat[,mapped] = 0
    }
  }
  
  ## get mu and Sigma of probability distr. thats defined by a crf
  ## where X is the latent factor matrix
  res = create_crf(1,2,as.vector(t(lf_mat)), adj_mat)
  
  ## take a sample of this matrix (make sure all values above 0)
  m =  mvrnorm(n = 1, res[[2]], res[[1]], tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  m = (m-min(m))+1
  mat = matrix(m, nrow = n_drugs, byrow=T)
  ##myImagePlot(mat)
  
  train = sample_training_data(mat)
  
  return(list(mat, train))
  
}

generate_dataset_scenario_2 = function(n_drugs, n_targets){
  
  ## generate a matrix of values with underlying latent factors
  lf_mat = generate_latent_factor_mat(n_drugs, n_targets)
  ##myImagePlot(lf_mat)
  
  ## create an example adjacency matrix 
  adj_mat = create_example_adj_mat(n_drugs, n_targets)
  
  ## and remove edges of last clique
  for (i in 1:n_drugs){
    for (k in 51:70){
      mapped = (i-1)*n_targets+k
      adj_mat[mapped,] = 0
      adj_mat[,mapped] = 0
    }
  }
  
  ## get mu and Sigma of probability distr. thats defined by a crf
  ## where X is the latent factor matrix
  res = create_crf(1,2,as.vector(t(lf_mat)), adj_mat)
  
  ## take a sample of this matrix (make sure all values above 0)
  m =  mvrnorm(n = 1, res[[2]], res[[1]], tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  m = (m-min(m))+1
  mat = matrix(m, nrow = n_drugs, byrow=T)
  ##myImagePlot(mat)
  
  train = sample_training_data(mat)
  
  return(list(mat, train))
  
}

generate_dataset_scenario_3 = function(n_drugs, n_targets){
  
  ## generate a matrix of values with underlying latent factors
  lf_mat = generate_latent_factor_mat(n_drugs, n_targets)
  ##myImagePlot(lf_mat)
  
  ## create an example adjacency matrix 
  adj_mat = create_example_adj_mat(n_drugs, n_targets)
  
  ## and remove edges of last clique
  for (i in 1:n_drugs){
    for (k in 51:70){
      mapped = (i-1)*n_targets+k
      adj_mat[mapped,] = 0
      adj_mat[,mapped] = 0
    }
  }
  
  ## get mu and Sigma of probability distr. thats defined by a crf
  ## where X is the latent factor matrix
  res = create_crf(1,2,as.vector(t(lf_mat)), adj_mat)
  
  ## take a sample of this matrix (make sure all values above 0)
  m =  mvrnorm(n = 1, res[[2]], res[[1]], tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  m = (m-min(m))+1
  mat = matrix(m, nrow = n_drugs, byrow=T)
  ##myImagePlot(mat)
  
  train = sample_training_data_sparser(mat)
  
  return(list(mat, train))
  
}

generate_dataset_scenario_4 = function(n_drugs, n_targets){
  
  ## generate a matrix of values with underlying latent factors
  lf_mat = generate_latent_factor_mat(n_drugs, n_targets)
  ##myImagePlot(lf_mat)
  
  ## create an example adjacency matrix 
  adj_mat = create_example_adj_mat(n_drugs, n_targets)
  
  ## and remove edges of last clique
  for (i in 1:n_drugs){
    for (k in 41:70){
      mapped = (i-1)*n_targets+k
      adj_mat[mapped,] = 0
      adj_mat[,mapped] = 0
    }
  }
  
  ## get mu and Sigma of probability distr. thats defined by a crf
  ## where X is the latent factor matrix
  res = create_crf(1,2,as.vector(t(lf_mat)), adj_mat)
  
  ## take a sample of this matrix (make sure all values above 0)
  m =  mvrnorm(n = 1, res[[2]], res[[1]], tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  m = (m-min(m))+1
  mat = matrix(m, nrow = n_drugs, byrow=T)
  ##myImagePlot(mat)
  
  train = sample_training_data(mat)
  
  return(list(mat, train))
  
}

generate_dataset_scenario_9 = function(n_drugs, n_targets){
  
  ## generate a matrix of values with underlying latent factors
  lf_mat = generate_latent_factor_mat(n_drugs, n_targets)
  myImagePlot(lf_mat)
  
  ## create an example adjacency matrix 
  adj_mat = create_example_adj_mat(n_drugs, n_targets)
  
  ## get mu and Sigma of probability distr. thats defined by a crf
  ## where X is the latent factor matrix
  res = create_crf(0.4,0.27,as.vector(t(lf_mat)), adj_mat)
  
  ## take a sample of this matrix (make sure all values above 0)
  m =  mvrnorm(n = 1, res[[2]], res[[1]], tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  m = (m-min(m))+1
  mat = matrix(m, nrow = n_drugs, byrow=T)
  ##myImagePlot(mat)
  
  train = sample_training_data(mat)
  
  return(list(mat, train))
  
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

make_crf_predictions_by_row = function(alpha, beta, train, adj_mat, X){
  
  n_drugs = nrow(train)
  n_targets = ncol(train)
  
  preds = matrix(NA, nrow = n_drugs, ncol = n_targets)
  
  A = matrix(0, nrow = n_drugs, ncol = n_drugs)
  B = matrix(0, nrow = n_drugs, ncol = n_drugs)
  
  for (i in 1:n_drugs){
    for (j in 1:n_drugs){
      if (i==j){
        A[i,j] = sum(alpha)
        B[i,j] = beta*sum(adj_mat[i,])
      } else {
        B[i,j] = -beta*adj_mat[i,j]
      }
    }
  }
  cat('inverting matrix..\n')
  Sigma1 = 2*(A+B)
  require("MASS")
  temp = chol(Sigma1) 
  Sigma = chol2inv(temp)
  
  for (t in 1:n_targets){
    cat('making predictions for column ',t,'\n')
    b = 2*X[,t]*alpha
    mu = Sigma%*%b
    
    preds[,t] = get_crf_prediction(Sigma, mu, train[,t])
  }
  
  return(preds)
  
}

make_crf_predictions_for_row = function(alpha, beta, train, adj_mat, X, col){
  
  n_drugs = nrow(train)
  n_targets = ncol(train)
  
  A = matrix(0, nrow = n_drugs, ncol = n_drugs)
  B = matrix(0, nrow = n_drugs, ncol = n_drugs)
  
  for (i in 1:n_drugs){
    for (j in 1:n_drugs){
      if (i==j){
        A[i,j] = sum(alpha)
        B[i,j] = beta*sum(adj_mat[i,])
      } else {
        B[i,j] = -beta*adj_mat[i,j]
      }
    }
  }
  cat('inverting matrix..\n')
  Sigma1 = 2*(A+B)
  require("MASS")
  temp = chol(Sigma1) 
  Sigma = chol2inv(temp)
  
  b = 2*X[,col]*alpha
  mu = Sigma%*%b
  
  preds = get_crf_prediction(Sigma, mu, train[,col])
  
  return(preds)
  
}

make_crf_predictions_row = function(alpha, beta, column, adj_mat, X){
  
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
  
  cat('inverting matrix..\n')
  Sigma1 = 2*(A+B)
  require("MASS")
  inds = which(Sigma1!=0, arr.ind = T)
  sm = sparseMatrix(i = inds[,1], j = inds[,2], x = Sigma1[inds])
  #temp = chol(Sigma1) 
  Sigma = chol2inv(chol(sm))
  #Sigma = chol2inv(temp)
  
  b = 2*X*alpha
  mu = Sigma%*%b
  
  preds = get_crf_prediction(Sigma, mu, column)
  
  return(preds)
  
}

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

make_adjacency_mat = function(sim_mat){
  
  n_drugs = nrow(sim_mat)
  
  adj_mat = matrix(0, nrow = n_drugs, ncol = n_drugs)
  for (i in 1:n_drugs){
    inds = which(sim_mat[i,]>0.9 & sim_mat[i,]<1)
    #inds = order(simi_mat[i,], decreasing = T)[2:5]
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

make_adjacency_mat_targets = function(simi_mat, thresh){
  
  n_targets = nrow(simi_mat)
  
  adj_mat = matrix(0, nrow = n_targets, ncol = n_targets)
  for (i in 1:n_targets){
    inds = which(simi_mat[i,]>thresh & simi_mat[i,]<1)
    #inds = order(simi_mat[i,], decreasing = T)[2:5]
    #if (length(inds)<4){
    #  inds = order(simi_mat[i,], decreasing = T)[1:5]
    #  inds = setdiff(inds, i)
    #}
    adj_mat[i, inds] = 1
    adj_mat[inds,i ] = 1
    #adj_mat[i, inds] = simi_mat[i, inds]
    #adj_mat[inds, i] = simi_mat[inds, i]
  }
  
  return(adj_mat)
  
}

make_training_adj_mat_for_row = function(dt_mat, simi_mat, row, thresh){
  
  n_train = length(which(dt_mat[row,]>=0))
  adj_temp = matrix(0, nrow = n_train, ncol = n_train)
  temp = which(dt_mat[row,]>0)
  for (k in temp){
    neighbors = which(simi_mat[k,]>thresh & simi_mat[k,]<1)
    #neighbors = order(simi_mat[k,], decreasing = T)[2:5]
    #if (length(neighbors)<4){
    #  neighbors = order(simi_mat[k,], decreasing = T)[2:5]
    #}
    for (n in neighbors){
      if (n %in% temp){
        adj_temp[match(k, temp), match(n, temp)] = 1
        adj_temp[match(n, temp), match(k, temp)] = 1
      }
    }
  }
  
  return(adj_temp)
}

make_training_adj_mat_for_column = function(dt_mat, simi_mat, col){
  
  n_train = length(which(dt_mat[,col]>=0))
  adj_temp = matrix(0, nrow = n_train, ncol = n_train)
  temp = which(dt_mat[,col]>0)
  for (k in temp){
    #cat('got observation for drug ',k,'\n')
    neighbors = which(simi_mat[k,]>0.9 & simi_mat[k,]<1)
    neighbors = order(simi_mat[k,], decreasing = T)[2:5]
    #if (length(neighbors)<4){
    #  inds = order(sim_mat[i,], decreasing = T)[1:5]
    #  inds = setdiff(inds, i)
      #neighbors = order(simi_mat[k,], decreasing = T)[2:5]
    #  neighbors = inds
    #}
    #cat('which is adj to ',neighbors,'\n')
    #noise = sample(temp, 5)
    #for (n in noise){
    #  if (n!=k){
    #    adj_temp[match(k, temp), match(n, temp)] = 1
    #    adj_temp[match(n, temp), match(k, temp)] = 1
    #  }
    #}
    for (n in neighbors){
      if (n %in% temp){
        #cat('neigbor ',n,' was also observed for this target\n')
        #cat(match(n, temp),'\n')
        adj_temp[match(k, temp), match(n, temp)] = 1
        adj_temp[match(n, temp), match(k, temp)] = 1
        #adj_temp[match(k, temp), match(n, temp)] = simi_mat[k, n]
        #adj_temp[match(n, temp), match(k, temp)] = simi_mat[n, k]
      }
    }
  }
  return (adj_temp)
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