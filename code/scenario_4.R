## same as scenario 1 but the last three groups are noise

n_drugs = 50
n_targets = 70

data = generate_dataset_scenario_4(n_drugs, n_targets)
complete = data[[1]]
train = data[[2]]
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