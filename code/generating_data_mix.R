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
for (i in 1:7){
  clique_list[[i]] = c(((i-1)*10+1):(i*10))
}

adj_mat = matrix(0, nrow = n_drugs*n_targets, ncol = n_drugs*n_targets)
## now go over each clique and for each clique member a, connect it to all other clique members
for (i in 1:length(clique_list)){
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

alpha = 0.5
beta = 1

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

X = rep(0, n_col*n_row)
b = 2*X*alpha

Sigma1 = 2*(A+B)
require("MASS")
Sigma = ginv(Sigma1)

mu = Sigma%*%b

m =  mvrnorm(n = 1, mu, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

m = (m-min(m))+1

mat = matrix(m, nrow = n_row, byrow=T)

myImagePlot(mat)


