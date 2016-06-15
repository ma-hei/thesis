#trying to build a CRF that learns different parameters for rows and columns at once

n_drugs = 40
n_targets = 60

target_beta = list()
# for each column, set a different beta
for (i in 1:10){
  target_beta[[i]] = 5
}
for (i in 11:20){
  target_beta[[i]] = 5
}
for (i in 21:30){
  target_beta[[i]] = 0.001
}
for (i in 31:40){
  target_beta[[i]] = 3
}
for (i in 41:50){
  target_beta[[i]] = 3
}
for (i in 51:60){
  target_beta[[i]] = 0.005
}

# I have 50 different adjacency mats.. one for each column
# make it into 1

row_adj_mat = matrix(0, nrow = n_drugs, ncol = n_drugs)
for (i in 2:10){
  row_adj_mat[1,i] = 1
  row_adj_mat[i,1] = 1
}
for (i in 22:30){
  row_adj_mat[21,i] = 1
  row_adj_mat[i,21] = 1
}

col_adj_mat = matrix(0, nrow = n_targets, ncol = n_targets)
col_adj_mat[1,2] = 1
col_adj_mat[2,1] = 1
col_adj_mat[1,3] = 1
col_adj_mat[3,1] = 1
col_adj_mat[1,4] = 1
col_adj_mat[4,1] = 1

# set one row beta
row_beta = 1

# one alpha
alpha = 1

# create matrices A and B
A = matrix(0, nrow = n_drugs*n_targets, ncol = n_drugs*n_targets)
B = matrix(0, nrow = n_drugs*n_targets, ncol = n_drugs*n_targets)

for (d in 1:n_drugs){
  for (t in 1:n_targets){
    row_neighbors = which(row_adj_mat[d,]==1)
    col_neighbors = which(col_adj_mat[t,]==1)
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
    nrow_neighbors = length(which(row_adj_mat[d,]>0))
    ncol_neighbors = length(which(col_adj_mat[t,]>0))
    sum = target_beta[[t]]*nrow_neighbors+row_beta*ncol_neighbors
    B[temp_a, temp_a] = sum
  }
}


##
lf_mat = generate_latent_factor_mat(n_drugs, n_targets)
myImagePlot(lf_mat)
X = as.vector(t(lf_mat))
b = 2*X*alpha

Sigma1 = 2*(A+B)
require("MASS")
temp = chol(Sigma1) 
Sigma = chol2inv(temp)
mu = Sigma%*%b
mu = (mu-min(mu))+1
mat = matrix(mu, nrow = n_drugs, byrow=T)
myImagePlot(mat)

## sample some training data
training_data = matrix(-1, nrow = n_drugs, ncol = n_targets)
for (d in 1:n_drugs){
  temp = sample(1:n_targets, 10)
  for (t in temp){
    training_data[d,t] = mat[d,t]
  }
}
myImagePlot(training_data)
## 

#training_data_mat = mat
training_data_mat = training_data

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

##
#mf_train = get_mf_cv(training_data_mat, iters = 400)
#X_train = mf_train[[1]]

ans = make_B_matrices(training_data_mat = training_data_mat, row_adj_mat = row_adj_mat, col_adj_mat = col_adj_mat)
Bcol_mats = ans[[1]]
Brow_mat = ans[[2]]

mf_train = get_mf_cv(training_data_mat, iters = 400)
X_train = t(mf_train[[2]])[which(!is.na(t(mf_train[[2]])))]
y_vec = t(training_data_mat)[which(t(training_data_mat)>=0)]
sqrt(mean((X_train - y_vec)^2))

target_beta_train_list = list()
for (i in 1:n_targets){
  #target_beta_train_list[[i]] = target_beta[[i]]
  target_beta_train_list[[i]] = 0.001
}
#row_beta_train = row_beta
row_beta_train = 0.001
#alpha_train = alpha
alpha_train = 0.001

eta = 0.01

for (i in 1:1000){
  cat('iter ',i,'.. alpha: ',alpha_train,' row beta: ',row_beta_train,'\n')
  cat('iter ',i,'.. col 1 beta: ' , target_beta_train_list[[1]],'\n')
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
  } else if (alpha_train>10){
    alpha_train = 10
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
    log_target_beta_train = log_target_beta_train + eta *grad_log_target_beta_train
    target_beta_train_list[[t]] = exp(log_target_beta_train)
    if (target_beta_train_list[[t]]<1e-10){
      target_beta_train_list[[t]] = 1e-10
    } else if (target_beta_train_list[[t]]>10){
      target_beta_train_list[[t]] = 10
    }
  }
  
  log_row_beta_train = log(row_beta_train)
  grad_row_beta_train = get_gradient_rows(y_vec = y_vec, B_row = Brow_mat, mu = mu, Sigma = Sigma)
  grad_log_row_beta_train = row_beta_train*grad_row_beta_train
  log_row_beta_train = log_row_beta_train + eta*grad_log_row_beta_train
  row_beta_train = exp(log_row_beta_train)
  if (row_beta_train<1e-10){
    row_beta_train = 1e-10
  } else if (row_beta_train>10){
    row_beta_train = 10
  }
}

train_crf(training_data_mat = training_data, row_adj_mat = row_adj_mat, col_adj_mat = col_adj_mat, eta = 0.01, crf_iters = 200)
