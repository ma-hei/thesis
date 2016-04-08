
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
    ## train on row
    ## adj_train = adj_mat[(row-1)*n_targets + which(train[row,]>0),(row-1)*n_targets + which(train[row,]>0)]
    adj_train = adj_list[[col]]
    # y = train[row, which(train[row,]>0)]
    y = train[which(train[,col]>0),col]
    # X = temp[row, which(!is.na(temp[row,]))]
    X = temp[which(!is.na(temp[,col])),col]
    
    # n_train = length(which(train[row,]>0))
    n_train = length(which(train[,col]>0))
    
    if (n_train>0){
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

get_crf_prediction = function(Sigma, mu, train){
  
  unknown = which(as.vector(t(train))<0)
  known = which(as.vector(t(train))>0)
  
  Sigma12 = Sigma[unknown, known]
  Sigma22 = Sigma[known, known]
  Sigma221 = chol2inv(chol(Sigma22))
  
  mu_ = mu[unknown] + Sigma12%*%Sigma221%*%(as.vector(t(train))[known] - mu[known])
  
  mu_all = rep(0, length(mu))
  mu_all[unknown] = mu_
  mu_all[known] = as.vector(t(train))[known]
  
  mat = matrix(mu_all, nrow = n_drugs, byrow=T)
  
  return(mat)
  
}

generate_dataset = function(n_drugs, n_targets){
  
## generate a matrix of values with underlying latent factors
  lf_mat = generate_latent_factor_mat(n_drugs, n_targets)
  ##myImagePlot(lf_mat)
  
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
  
  return(list(mat, train))
  
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
  ##myImagePlot(lf_mat)
  
  ## create an example adjacency matrix 
  adj_mat = create_example_adj_mat(n_drugs, n_targets)
  
  ## get mu and Sigma of probability distr. thats defined by a crf
  ## where X is the latent factor matrix
  res = create_crf(0.5,0.07,as.vector(t(lf_mat)), adj_mat)
  
  ## take a sample of this matrix (make sure all values above 0)
  m =  mvrnorm(n = 1, res[[2]], res[[1]], tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  m = (m-min(m))+1
  mat = matrix(m, nrow = n_drugs, byrow=T)
  ##myImagePlot(mat)
  
  train = sample_training_data(mat)
  
  return(list(mat, train))
  
}



