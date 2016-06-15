
load('../data/mix_dataset.Rda')
target_sim = read.table('../data/mix_dataset_target_target_sim.txt')
target_sim = as.matrix(target_sim)
drug_sim = read.table('../data/mix_dataset_drug_drug_sim.txt')
drug_sim = as.matrix(drug_sim)

n_drugs = length(unique(mix_dataset[,1]))
n_targets = length(unique(mix_dataset[,2]))

n_folds = 5
test_folds = get_folds(mix_dataset, n_folds)

crf_predictions = rep(NA, nrow(mix_dataset))
crf_predictions_single = rep(NA, nrow(mix_dataset))
mf_predictions = rep(NA, nrow(mix_dataset))

ans = get_target_clusters(target_sim, 0.6)
target_cluster_list = ans[[1]]
remaining_targets = ans[[2]]

drug_adj_mat = make_adjacency_mat_4n(drug_sim)
drug_adj_mat = make_adjacency_mat(sim_mat = drug_sim, thresh = 0.9)
cluster_1_adj_mat = drug_adj_mat[drugs_1, drugs_1]
cluster_2_adj_mat = drug_adj_mat[drugs_2, drugs_2]

target_adj_mat = make_adjacency_mat_targets(target_sim, 0.6)

i = 1
test_ind = test_folds[[i]]
train_ind = setdiff(1:nrow(mix_dataset),test_ind)

## make the training data into matrix format
dt_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
dt_mat[mix_dataset[train_ind,c(1,2)]] = mix_dataset[train_ind,3]

cat('getting MF cv-prediction on training data..\n')
mf_preds_train = get_mf_cv(dt_mat, 400)
mf_preds_train_mat = mf_preds_train[[2]]

cat('getting MF predictions for complete matrix..\n')
mf_preds_all = get_libmf_prediction(dt_mat, 400)

k = 1
cat('\ntarget cluster ',k,'/',length(target_cluster_list),'.. targets of cluster: ',target_cluster_list[[k]],'\n')
n_train = length(which(dt_mat[,target_cluster_list[[k]]]>=0))
cat('number of observations in cluster:', n_train,'\n')

target_set = target_cluster_list[[k]]
adj_mat_cols = target_adj_mat[target_set, target_set]
adj_mat_rows = drug_adj_mat

dt_mat_temp = dt_mat[,target_set]
mf_preds_train_mat_temp = mf_preds_train_mat[,target_set]
mf_preds_all_mat_temp = mf_preds_all[,target_set]

eta = 0.01
params = train_crf_dual_simple(training_data_mat = dt_mat_temp, row_adj_mat = adj_mat_rows, col_adj_mat = adj_mat_cols, mf_preds_train_mat = mf_preds_train_mat_temp, eta = eta, crf_iters = 600)
crf_preds_temp = simple_crf_predict(training_data_mat = dt_mat_temp, row_adj_mat = adj_mat_rows, col_adj_mat = adj_mat_cols, alpha = params[[1]], row_beta = params[[2]], col_beta = params[[3]], X_mat = mf_preds_all_mat_temp)



