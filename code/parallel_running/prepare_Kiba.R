source('crf_utils.R')
source('recosystem.R')
require('Matrix')
source('target_clustering.R')
require(foreach)
require(doParallel)
require(parallel)

load('../kiba_2_triplet.rda')
load('../kiba_2_target_sim.rda')
load('../kiba_2_drug_sim.rda')

dt_triplet[,3] = dt_triplet[,3]*-1
dt_triplet[,3] = dt_triplet[,3] - min(dt_triplet[,3]) 

dataset_triplets = dt_triplet


target_sim = sim_mat
drug_sim = drug_sim_mat

n_drugs = length(unique(dt_triplet[,1]))
n_targets = length(unique(dt_triplet[,2]))

args = commandArgs(trailingOnly = T)
drug_thresh = as.numeric(as.numeric(args[1]))
target_thresh = as.numeric(as.numeric(args[2]))

target_adj_mat = make_adj_mat_thresh(sim_mat = target_sim, thresh = target_thresh, FALSE)
drug_adj_mat = make_adj_mat_thresh(sim_mat = drug_sim, thresh = drug_thresh, FALSE)

n_folds = 5
test_folds = get_folds(dataset_triplets, n_folds)

save(test_folds, dataset_triplets, drug_adj_mat, target_adj_mat, file = "5fold_cv_data_kiba.Rda")
