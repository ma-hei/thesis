require("MASS")

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

train_crf = function(training_data, adj_mats, X, crf_iters, eta){
  
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
    if ((it%%100)==0){
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
    
    counter = 0
    for (i in 1:length(Bk_mats)){
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

make_adj_mat_thresh = function(sim_mat, thresh, add_nn){
  
  n = nrow(sim_mat)
  adj_mat = matrix(0, nrow = n, ncol = n)
  for (i in 1:n){
    inds = which(sim_mat[i,]>thresh)
    inds = setdiff(inds, i)
    if (add_nn==TRUE){
      if (length(inds)<1){
          inds = order(sim_mat[i,], decreasing = T)[1:2]
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

dataset = read.table('../data/known_drug-target_interaction_affinities_pKi__Metz_et_al.2011.txt')
dataset = as.matrix(dataset)

dataset_triplets = make_triplets(dataset)

n_drugs = nrow(dataset)
n_targets = ncol(dataset)

drug_sim = read.table('../data/drug-drug_similarities_2D__Metz_et_al.2011.txt')
drug_sim = as.matrix(drug_sim)
drug_sim = drug_sim/100

target_sim = read.table('../data/target-target_similarities_WS_normalized__Metz_et_al.2011.txt')
target_sim = as.matrix(target_sim)
target_sim = target_sim/100

target_adj_mat = make_adj_mat_thresh(sim_mat = target_sim, thresh = 0.66, FALSE)
drug_adj_mat = make_adj_mat_thresh(sim_mat = drug_sim, thresh = 0.75, TRUE)

n_folds = 5
test_folds = get_folds(dataset_triplets, n_folds)

fold = 1
test_ind = test_folds[[fold]]
train_ind = setdiff(1:nrow(dataset_triplets),test_ind)

train_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
train_mat[dataset_triplets[train_ind,c(1,2)]] = dataset_triplets[train_ind,3]

c('getting MF predictions for train data.. \n')
mf_preds_train = get_mf_cv(train_mat, 400)
mf_preds_train_mat = mf_preds_train[[2]]

cat('getting MF predictions for complete matrix..\n')
mf_preds_all = get_libmf_prediction(train_mat, 400)

## train crf for column 26

crf_predictions = rep(NA, nrow(dataset_triplets))
for (t in 1:n_targets){
  
  cat('target ',t,'... ',length(which(train_mat[,t]>=0)),' observations\n')
  
  if (length(which(train_mat[,t]>=0))>1){
    
    training_data = train_mat[which(train_mat[,t]>0), t]
    inds = which(train_mat[,t]>=0)
    training_adj_mat = drug_adj_mat[inds,inds]
    inds = which(training_adj_mat>0, arr.ind = T)
    training_adj_mat = sparseMatrix(i = inds[,1], j = inds[,2], x = training_adj_mat[inds], dims = dim(training_adj_mat))
    X = mf_preds_train_mat[which(!is.na(mf_preds_train_mat[,t])),t]
    
    adj_mats = list()
    adj_mats[[1]] = training_adj_mat
    params = train_crf(training_data = training_data, adj_mats = adj_mats, X = X, crf_iters = 400, eta = 0.008)
    
    adj_mats = list()
    adj_mats[[1]] = drug_adj_mat
    prediction = crf_predict(train_data = train_mat[,t], X = mf_preds_all[,t], alpha = params[[1]], adj_mats = adj_mats, beta = params[[2]])
    
    inds = which(dataset_triplets[test_ind,2] == t)
    crf_predictions[test_ind[inds]] = prediction[dataset_triplets[test_ind[inds],1]]

    
  } else {
    inds = which(dataset_triplets[test_ind,2] == t)
    crf_predictions[test_ind[inds]] = mf_preds_all[dataset_triplets[test_ind[inds],1],t]
  }
  
  inds = which(!is.na(crf_predictions))
  crf_metrics = get_metrics(crf_predictions[inds], dataset_triplets[inds,3], 7.6)
  
  cat('auc so far, crf: ',round(crf_metrics[[2]], digits = 5),'\n')
  cat('aupr so far, crf: ',round(crf_metrics[[3]], digits = 5),'\n\n')
  
}

## train crf for row 142

## train crf for columns 10, 53, 10, 19

## train crf for rows 19, 294, 1200, 1





