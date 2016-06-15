source('crf_utils.R')
source('recosystem.R')
source('target_clustering.R')
require('Matrix')

dataset = read.table('../../data/drug-target_interaction_affinities_Kd__Davis_et_al.2011.txt')
dataset = as.matrix(dataset)
dataset = -log(dataset/(1e9),10)

dataset_triplets = make_triplets(dataset)

n_drugs = nrow(dataset)
n_targets = ncol(dataset)

drug_sim = read.table('../../data/drug-drug_similarities_2D.txt')
drug_sim = as.matrix(drug_sim)

target_sim = read.table('../../data/target-target_similarities_WS_normalized.txt')
target_sim = as.matrix(target_sim)
target_sim = target_sim/100

args = commandArgs(trailingOnly = T)
drug_thresh = as.numeric(as.numeric(args[1]))
target_thresh = as.numeric(as.numeric(args[2]))

target_adj_mat = make_adj_mat_thresh(sim_mat = target_sim, thresh = target_thresh, FALSE)
drug_adj_mat = make_adj_mat_thresh(sim_mat = drug_sim, thresh = drug_thresh, TRUE, 2)

n_folds = 5
test_folds = get_folds(dataset_triplets, n_folds)

save(test_folds, dataset_triplets, drug_adj_mat, target_adj_mat, file = "5fold_cv_data_davis.Rda")

