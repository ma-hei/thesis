## train the crf row wise
## input: a training matrix and an adjacency matrix
## split up the matrices and train the crf row wise


n_drugs = 50
n_targets = 70

data = generate_dataset(n_drugs, n_targets)
complete = data[[1]]
train = data[[2]]
myImagePlot(complete)
myImagePlot(train)

adj_mat = create_example_adj_mat(n_drugs, n_targets)

params = train_crf_by_row(train = train, adj = adj_mat, 100, 1000)

