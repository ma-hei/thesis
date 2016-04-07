## first I generate training data for a 50 x 70 matrix (50 imaginary drugs and 70 imaginary targets)
## I generate the training data such that it has both underlying latent factors and pairs of targets
## which have similar binding values (so a mix of MF + CRF should work well to model this kind of data)

n_drugs = 50
n_targets = 70

## first, randomly generate latent factors for the drugs
## I will make all latent factors > 0, s.th. I get only positive numbers in the matrix
## assume dimension of latent factors k = 20
k = 20
lf_drugs = matrix(0, 50, k)
for (i in 1:n_drugs){
  latent_factors_temp = rnorm(k, 0, 1)
  ## make sure all latent factors > 0 
  if (min(latent_factors_temp)<0){
    latent_factors_temp = latent_factors_temp - min(latent_factors_temp)
  }
  lf_drugs[i,] = latent_factors_temp
}

lf_targets = matrix(0, 70, k)
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

myImagePlot(training_data_lf)

## now I'm building a CRF. First I build the adjacency matrix of the CRF nodes.
## every matrix cell is a CRF node, so the size of the adjacency matrix is 
## (n_drugs*n_targets)^2 (for n_drugs = 50, n_targets = 70 -> 12.250.000 nodes.. very large..)
## I will have to come up with ways to split up the CRF for the practical application

## I let the first clique be the columns 1-10, second clique be columns 11-20, ...
## inside each clique, I connect each column with all other columns in the clique
## I get 7 cliques a 10 columns

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

## now I construct a multivariate gaussian from the probability distribution
## that is defined by the CRF with above structure. This step is explained here 
## http://www.cl.cam.ac.uk/~pr10/publications/fg13.pdf

alpha = 1
beta = 2

n_cell = n_drugs*n_targets

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

#X = rep(0, n_drugs*n_targets)
X = as.vector(t(training_data_lf))
b = 2*X*alpha

Sigma1 = 2*(A+B)
require("MASS")
#Sigma = ginv(Sigma1)
Sigma = chol2inv(chol(Sigma1))
mu = Sigma%*%b

m =  mvrnorm(n = 1, mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
m = (m-min(m))+1
mat = matrix(m, nrow = n_drugs, byrow=T)

myImagePlot(mat)

## now I have a fully populated value matrix with underlying latent factors and 
## in which 5 groups of rows exists that have similar values 
## I will remove most of the values to generate a training matrix

training_data = matrix(-1, nrow = n_drugs, ncol = n_targets)

## go over each row, in each row select one group for which no training data will
## be availble, for all other groups, select 2 - 5 cells that go into training data

for (i in 1:n_drugs){
  empty_group = sample(1:length(clique_list),3)
  for (k in 1:length(clique_list)){
    if (!k%in%empty_group){
      select_size = sample(1:2,1)
      selection = sample(1:length(clique_list[[k]]),select_size)
      for (x in 1:select_size){
        training_data[i,clique_list[[k]][selection[x]]] = mat[i,clique_list[[k]][selection[x]]]
      }
    }
  }
}

myImagePlot(training_data)

## now we check how good matrix factorization can recover the data

train_mat_triplet = matrix(0, nrow = length(which(training_data>0)), ncol = 3)
train_mat_triplet[,1] = which(training_data>0, arr.ind = T)[,1]
train_mat_triplet[,2] = which(training_data>0, arr.ind = T)[,2]
train_mat_triplet[,3] = training_data[which(training_data>0, arr.ind = T)]

res = nmf(train_mat_triplet[,1:3], m = n_drugs, n = n_targets, k = 10,
          lambda = 1, gamma = 1, tol = 0.001, maxiter = 500, threshold = 2, bias=TRUE, interaction = TRUE)


pred_mf = nmf.predict(res,cbind(rep(1:n_drugs, 1, each = 70), rep(1:70, n_drugs)),
                      bias = TRUE,interaction = TRUE)

pred_mf_mat = matrix(pred_mf, nrow = n_drugs, byrow = T)

myImagePlot(pred_mf_mat)

sqrt(mean((pred_mf - m)^2))

## now I need to create training data for the CRF
## I need a new adjacency matrix that connects only the cells that are known in the training data

n_train = length(which(training_data>0))

ind = which(training_data<0,arr.ind=T)
cut_out = rep(0, length(which(training_data<0)))
for (i in 1:length(which(training_data<0))){
  row = ind[i,1]
  col = ind[i,2]
  mapped = (row-1)*n_targets+col
  cut_out[i] = mapped
}

##add noise to adj mat
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


adj_mat_train = adj_mat[-cut_out, -cut_out]

## y is the real values in the training data
y = t(training_data)[which(t(training_data>0), arr.ind = T)]
X = t(pred_mf_mat)[which(t(training_data>0), arr.ind = T)]

B1 = matrix(0, nrow = n_train, ncol = n_train)
for (i in 1:n_train){
  for (j in 1:n_train){
    if (i==j){
      B1[i,j] = sum(adj_mat_train[i,])
    } else {
      B1[i,j] = -adj_mat_train[i,j]
    }
  }
}

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


alpha = 0.001
beta1 = 0.001

log_alpha = log(alpha)
log_beta1 = log(beta1)

eta = 0.01

for (i in 1:200){
  
  cat(alpha,' ',beta1,'\n')
  gradients = compute_gradient(y,alpha, beta1, adj_mat_train, B1, X, n_train)
  grad_alpha = gradients[[1]]
  grad_beta1 = gradients[[2]]
  
  grad_log_alpha = alpha*grad_alpha
  grad_log_beta1 = beta1*grad_beta1
  
  log_alpha = log_alpha + eta * grad_log_alpha
  log_beta1 = log_beta1 + eta * grad_log_beta1
  
  alpha = exp(log_alpha)
  beta1 = exp(log_beta1)
  
}

## now lets see what the CRF predicts

alpha = alpha
beta = beta1

A = matrix(0, nrow = n_drugs*n_targets, ncol = n_drugs*n_targets)
B = matrix(0, nrow = n_drugs*n_targets, ncol = n_drugs*n_targets)

for (i in 1:(n_drugs*n_targets)){
  for (j in 1:(n_drugs*n_targets)){
    if (i==j){
      A[i,j] = sum(alpha)
      B[i,j] = beta*sum(adj_mat[i,])-beta*adj_mat[i,j]
    } else {
      B[i,j] = -beta*adj_mat[i,j]
    }
  }
}
X = as.vector(t(pred_mf_mat))
b = 2*X*alpha

Sigma1 = 2*(A+B)
require("MASS")
Sigma = chol2inv(chol(Sigma1))

mu = Sigma%*%b

pred_crf =  mu

mat = matrix(pred_crf, nrow = n_drugs, byrow=T)

myImagePlot(mat)

sqrt(mean((pred_crf - m)^2))

condition_crf = function(mu, Sigma, training_data){
  
  unknown = which(as.vector(t(training_data))<0)
  known = which(as.vector(t(training_data))>0)
  
  Sigma12 = Sigma[unknown, known]
  Sigma22 = Sigma[known, known]
  Sigma221 = chol2inv(chol(Sigma22))
  
  mu_ = mu[unknown] + Sigma12%*%Sigma221%*%(as.vector(t(training_data))[known] - mu[known])
  
  mu_all = rep(0, length(mu))
  mu_all[unknown] = mu_
  mu_all[known] = as.vector(t(training_data))[known]
  
  mat = matrix(mu_all, nrow = n_drugs, byrow=T)
  
  myImagePlot(mat)
  
  sqrt(mean((mu_all - m)^2))
  
}


