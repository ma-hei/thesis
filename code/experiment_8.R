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

i = 1
test_ind = test_folds[[i]]
train_ind = setdiff(1:nrow(mix_dataset),test_ind)

dt_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
dt_mat[mix_dataset[train_ind,c(1,2)]] = mix_dataset[train_ind,3]

cat('getting MF cv-prediction on training data..\n')
mf_preds_train = get_mf_cv(dt_mat, 400)
mf_preds_train_mat = mf_preds_train[[2]]

cat('getting MF predictions for complete matrix..\n')
mf_preds_all = get_libmf_prediction(dt_mat, 400)

target_adj_mat = make_adjacency_mat_targets(target_sim, 0.6)
drug_adj_mat = make_adjacency_mat(drug_sim)

#target_set = c(32,33,203,490)
target_set = c(147,183,415,416,479,484)
#target_set = c(12, 180, 390 )
#target_set = c(11, 143, 287)
target_set = c(8, 139, 516)
adj_mat_cols = target_adj_mat[target_set, target_set]
adj_mat_rows = drug_adj_mat
dt_mat_temp = dt_mat[,target_set]
mf_preds_train_mat_temp = mf_preds_train_mat[,target_set]
mf_preds_all_mat_temp = mf_preds_all[,target_set]
params = train_crf_dual_simple(training_data_mat = dt_mat_temp, row_adj_mat = adj_mat_rows, col_adj_mat = adj_mat_cols, mf_preds_train_mat = mf_preds_train_mat_temp, eta = 0.01, crf_iters = 400)
crf_preds_temp = simple_crf_predict(training_data_mat = dt_mat_temp, row_adj_mat = adj_mat_rows, col_adj_mat = adj_mat_cols, alpha = params[[1]], row_beta = params[[2]], col_beta = params[[3]], X_mat = mf_preds_all_mat_temp)

print_prediction_metrics_set(target_set = target_set, mix_dataset = mix_dataset, test_ind = test_ind, crf_preds = crf_preds_temp, mf_preds_all_mat = mf_preds_all_mat_temp)

train_and_predict_single(mf_preds_train_mat = mf_preds_train_mat, dt_mat = dt_mat, sim_mat = drug_sim, mix_dataset = mix_dataset, mf_preds_all = mf_preds_all, target_set = target_set, test_ind = test_ind, adj_mat = drug_adj_mat)
