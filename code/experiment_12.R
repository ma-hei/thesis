## metz data set drug crfs

dataset = read.table('../data/known_drug-target_interaction_affinities_pKi__Metz_et_al.2011.txt')
dataset = as.matrix(test)

drug_sim = read.table('../data/drug-drug_similarities_2D__Metz_et_al.2011.txt')
drug_sim = as.matrix(drug_sim)
drug_sim = drug_sim/100

target_sim = read.table('../data/target-target_similarities_WS_normalized__Metz_et_al.2011.txt')
target_sim = as.matrix(target_sim)
target_sim = target_sim/100
target_clusters = get_target_clusters(dataset_triplets,target_sim, 0.69)
hist(target_sim)


inds = which(!is.na(dataset), arr.ind = T)
dataset_triplets = matrix(0, nrow = nrow(inds), ncol = 3)
dataset_triplets[,c(1,2)] = inds
dataset_triplets[,3] = dataset[inds]

n_drugs = length(unique(dataset_triplets[,1]))
n_targets = length(unique(dataset_triplets[,2]))

n_folds = 5
test_folds = get_folds(dataset_triplets, n_folds)

mf_predictions = rep(NA, nrow(dataset_triplets))
target_crf_predictions = rep(NA, nrow(dataset_triplets))

for (d in 1:n_drugs){
    
  
  
}