compute_gradient = function(y, alpha, beta1, adj_mat_train, B1, X, n_train){
  
  A = matrix(0, nrow = n_train, ncol = n_train)
  B = matrix(0, nrow = n_train, ncol = n_train)
  
  for (i in 1:(n_train)){
    for (j in 1:(n_train)){
      if (i==j){
        A[i,j] = sum(alpha)
        B[i,j] = beta1*sum(adj_mat_train[i,]) 
      } else {
        B[i,j] = -(beta1*adj_mat_train[i,j])
      }
    }
  }
  
  b = 2*X*alpha
  
  Sigma1 = 2*(A+B)
  Sigma = ginv(Sigma1)
  
  mu = Sigma%*%b
  
  grad_alpha = -t(y)%*%y + 2*t(y)%*%X - 2*t(X)%*%mu + t(mu)%*%mu + sum(diag(Sigma))
  
  grad_beta1 = -t(y)%*%B1%*%y + t(mu)%*%B1%*%mu + t(as.vector(Sigma))%*%as.vector(B1)
  
  return(list(grad_alpha, grad_beta1))
  
}

train_crf = function(train, adj, mf_iters, crf_iters){
  
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
  
  alpha = 0.001
  beta = 0.001
  
  log_alpha = log(alpha)
  log_beta = log(beta)
  
  eta = 0.01
  
  for (i in 1:crf_iters){
    
    cat(alpha,' ',beta,'\n')
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
  
  return(list(X,alpha, beta))
  
}

test = train_crf(train = training_data, adj = adj_mat, 100, 100)
