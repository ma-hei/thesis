## the dual crf with targets 147, 183, 415, 416, 479, 484

load('../data/mix_dataset.Rda')
target_sim = read.table('../data/mix_dataset_target_target_sim.txt')
target_sim = as.matrix(target_sim)
drug_sim = read.table('../data/mix_dataset_drug_drug_sim.txt')
drug_sim = as.matrix(drug_sim)

n_drugs = length(unique(mix_dataset[,1]))
n_targets = length(unique(mix_dataset[,2]))

test_folds = get_folds(mix_dataset, 5)

crf_predictions = rep(NA, nrow(mix_dataset))
mf_predictions = rep(NA, nrow(mix_dataset))

i = 5
test_ind = test_folds[[i]]
train_ind = setdiff(1:nrow(mix_dataset),test_ind)

dt_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
dt_mat[mix_dataset[train_ind,c(1,2)]] = mix_dataset[train_ind,3]

cat('getting MF cv-prediction on training data..\n')
mf_preds_train = get_mf_cv(dt_mat, 400)
mf_preds_train_mat = mf_preds_train[[2]]


target_adj_mat = make_adjacency_mat_targets(target_sim, 0.6)
drug_adj_mat = make_adjacency_mat(drug_sim)

adj_mat_cols = target_adj_mat[c(147, 183, 415, 416, 479, 484), c(147, 183, 415, 416, 479, 484)]
adj_mat_rows = drug_adj_mat

n_drugs = length(unique(mix_dataset[,1]))
n_targets = length(unique(mix_dataset[,2]))

dt_mat_temp = dt_mat[,c(147, 183, 415, 416, 479, 484)]
mf_preds = mf_preds_train_mat[,c(147, 183, 415, 416, 479, 484)]

params = train_crf(training_data_mat = dt_mat_temp, row_adj_mat = adj_mat_rows, col_adj_mat = adj_mat_cols, eta = 0.01, crf_iters = 1000, mf_preds = mf_preds)




