
n_drugs = 50
n_targets = 70

data = generate_dataset_scenario_1(n_drugs, n_targets)
complete = data[[1]]
train = data[[2]]
myImagePlot(complete)
myImagePlot(train)

adj_mat = create_example_adj_mat(n_drugs, n_targets)
params_all = train_crf(train = train, adj = adj_mat, 100, 50)

params_byrow = train_crf_by_row(train = train, adj = adj_mat, 100, 1000)

X = get_mf_prediction(train, 500)
## the rmse of the mf prediction is:
sqrt(mean((X - as.vector(t(complete)))^2))
mf_pred_mat = matrix(X, nrow = n_drugs, byrow = T)
myImagePlot(mf_pred_mat)
## get the Sigma and mu of the CRF with learned params and observations X
res = create_crf(params_all[[1]], params_all[[2]], X, adj_mat)
