## same as scenario 4 but with less training data

n_drugs = 50
n_targets = 70

data = generate_dataset_scenario_4(n_drugs, n_targets)
complete = data[[1]]
train = data[[2]]
## randomly remove an other 150 training values
temp = sample(length(which(train>0)), 150)
ind = which(train>0,arr.ind=T)
for (i in 1:150){
  row = ind[temp[i],1]
  col = ind[temp[i],2]
  train[row,col] = -1
}

myImagePlot(complete)
myImagePlot(train)

## now create an adjacency mat that connects all cliques (and thus contains some noise)
adj_mat = create_example_adj_mat(n_drugs, n_targets)
params = train_crf(train = train, adj = adj_mat, 100, 50)

## get the vector X of the mf predictions on the complete data
X = get_mf_prediction(train, 500)
## the rmse of the mf prediction is:
sqrt(mean((X - as.vector(t(complete)))^2))
##and looks as follows:
mf_pred_mat = matrix(X, nrow = n_drugs, byrow = T)
myImagePlot(mf_pred_mat)
## get the Sigma and mu of the CRF with learned params and observations X
res = create_crf(params[[1]], params[[2]], X, adj_mat)

## condition the resulting probability dist on the observed values in the training data
mat = get_crf_prediction(res[[1]],res[[2]],train)
myImagePlot(mat)
## the rmse of the crf predictions is:
sqrt(mean((as.vector(t(mat)) - as.vector(t(complete)))^2))