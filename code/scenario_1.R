## in this scenario I generate the data as before, but for the last group
## the values in the group are not similar. The cells in the last two groups
## are still connected when I train the CRF, this means there is now noise in the 
## CRF structure.

n_drugs = 50
n_targets = 70

data = generate_dataset_scenario_1(n_drugs, n_targets)
complete_1 = data[[1]]
train_1 = data[[2]]
myImagePlot(complete_1)
myImagePlot(train_1)

complete_1 = complete
train_1 = train

## now create an adjacency mat that connects all cliques (and thus contains some noise)
adj_mat = create_example_adj_mat(n_drugs, n_targets)
params = train_crf(train = train_1, adj = adj_mat, 500, 50)

## get the vector X of the mf predictions on the complete data
X = get_mf_prediction(train_1, 500)
## the rmse of the mf prediction is:
sqrt(mean((X - as.vector(t(complete_1)))^2))
##and looks as follows:
mf_pred_mat = matrix(X, nrow = n_drugs, byrow = T)
myImagePlot(mf_pred_mat)
## get the Sigma and mu of the CRF with learned params and observations X
res = create_crf(params[[1]], params[[2]], X, adj_mat)

## condition the resulting probability dist on the observed values in the training data
mat = get_crf_prediction(res[[1]],res[[2]],train_1)
myImagePlot(mat)
## the rmse of the crf predictions is:
sqrt(mean((as.vector(t(mat)) - as.vector(t(complete_1)))^2))