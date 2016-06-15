
dt_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)

for (i in 1:length(train_ind)){
  dt_mat[dt_triplet[train_ind[i],1], dt_triplet[train_ind[i],2]] = dt_triplet[train_ind[i],3]
}

for (i in 1:length(train_ind)){
  if (mf_preds_train_mat[dt_triplet[train_ind[i],1], dt_triplet[train_ind[i],2]]<0){
    cat(i)
  }
}

inds = which(dt_triplet[test_ind,2] == 1)
mf_preds_all[dt_triplet[test_ind[inds],c(1:2)]]


for (i in 1:n_folds){
  cat('fold ',i,'\n')
  
  test_ind = test_folds[[i]]
  train_ind = setdiff(1:nrow(dt_triplet),test_ind)
  
  dt_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
  dt_mat[dt_triplet[train_ind,c(1,2)]] = dt_triplet[train_ind,3]
  
  cat('getting MF predictions..\n')
  mf_preds_all = get_libmf_prediction(dt_mat, 400)

  rmse_fold = sqrt(mean((dt_triplet[test_ind,3]-mf_preds_all[dt_triplet[test_ind,c(1:2)]])^2))
  
  cat('rmse: ',rmse_fold,'\n')
  
  mf_preds_ = matrix(NA, nrow = n_drugs, ncol = n_targets)
  
  for (t in 1:n_targets){
    
    inds = which(dt_triplet[test_ind,2] == t)
    labels_test_col = dt_triplet[test_ind[inds], 3]
    mf_prediction_col = mf_preds_all[cbind(dt_triplet[test_ind[inds],1],dt_triplet[test_ind[inds],2])]
    
    mf_preds_[dt_triplet[test_ind[inds],1],t] = mf_prediction_col
    
  }
  
  rmse_fold_ = sqrt(mean((dt_triplet[test_ind,3]-mf_preds_[dt_triplet[test_ind,c(1:2)]])^2))
  
  cat('compare: ',rmse_fold_,'\n')
   
}

counter = 0
for (i in 1:n_drugs){
  if (length(which(drug_sim[i,]>0.8 & drug_sim[i,]<1))>0){
    counter = counter+1
  }
}

counter = 0
for (i in 1:n_targets){
  if (length(which(sim_mat[i,]>0.5 & sim_mat[i,]<1))>0){
    counter = counter+1
  }
}


counter = 0
for (i in 1:n_targets){
  if (length(which(adj_mat[i,]>0))>0){
    counter = counter+1
  }
}

counter = 0
for (i in 1:n_targets){
  if (length(which(target_adj_mat[i,]>0))>1){
    counter = counter+1
  }
}

inds = which(adj_mat_train_col_t>0,arr.ind=T)
for (i in 1:nrow(inds)){
  
  cat(training_vals_col_t[inds[i,1]],' ',training_vals_col_t[inds[i,2]],'\n')
}

load('kiba_triplet.rda')
load('kiba_drug_sim.rda')

n_drugs = length(unique(dt_triplet[,1]))
n_targets = length(unique(dt_triplet[,2]))

dt_triplet[,3] = dt_triplet[,3]*-1
dt_triplet[,3] = dt_triplet[,3] - min(dt_triplet[,3]) 
hist(dt_triplet[,3])

n_folds = 5
test_folds = get_folds(dt_triplet, n_folds)

adj_mat = make_adjacency_mat(drug_sim_mat)

i = 1
test_ind = test_folds[[i]]
train_ind = setdiff(1:nrow(dt_triplet),test_ind)

dt_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
dt_mat[dt_triplet[train_ind,c(1,2)]] = dt_triplet[train_ind,3]

cat('getting MF cv-prediction on training data..\n')
mf_preds_train = get_mf_cv(dt_mat, 400)
mf_preds_train_mat = mf_preds_train[[2]]

cat('getting MF predictions for complete matrix..\n')
mf_preds_all = get_libmf_prediction(dt_mat, 400)
t = 1


crf_predictions = rep(NA, nrow(dt_triplet))
mf_predictions = rep(NA, nrow(dt_triplet))


for (t in 1:n_targets){
  
  cat('fold ',i,', target ',t,'... ',length(which(dt_mat[,t]>=0)),' observations\n')
  
mf_pred_train_col_t = mf_preds_train_mat[which(!is.na(mf_preds_train_mat[,t])),t]
adj_mat_train_col_t = make_training_adj_mat_for_column(dt_mat, drug_sim_mat, t)
training_vals_col_t = dt_mat[which(dt_mat[,t]>=0),t]

if (length(which(dt_mat[,t]>=0))>500){
  eta = 0.001
} else{
  eta = 0.01
}

if(length(which(dt_mat[,t]>=0))>100){
  niter = 200
} else{
  niter = 1000
}

counter = 0
for (t in 1:n_targets){
  cat(length(which(target_sim[t,]>0.4 & target_sim[t,]<1)),'\n')
  if(length(which(target_sim[t,]>0.4 & target_sim[t,]<1))>0){
    counter = counter+1
  }
}

sum = 0
for (d in 1:n_drugs){
  sum = sum+length(which(row_adj_mat[d,]>0))*n_targets
}
for (t in 1:n_targets){
  sum = sum+length(which(col_adj_mat[t,]>0))*n_drugs
}
sum = 0
for (i in 1:50){
  sum = sum + length(cluster_list[[i]])
}

target_set_1 = cluster_list[[1]]
for (i in 2:141){
  target_set_1 = union(target_set_1, cluster_list[[i]])
}

n_drugs = length(unique(dataset_triplets[,1]))

counter = 0
drug_list = list()
for (i in 1:n_drugs){
  if (length(which(drug_adj_mat[i,]>0.9))>10){
    counter = counter+1
    drug_list[[counter]] = i
  }
}

counter = 0
for (i in 1:n_targets){
  if (length(which(target_adj_mat[i,]>0))>0){
    counter = counter+1
  }
}

counter = 0
for (i in 1:592){
  if (length(which(training_adj_mat[i,]>0))>0){
    counter = counter+1
  }
}

inds = which(dataset_triplets[test_ind,1]==d)
dataset_triplets[test_ind[inds[which(dataset_triplets[test_ind[inds],3]>=12.1)]],c(1:3)]

temp = 123
train_mat[drug_set,c(which(target_adj_mat[temp,]>0),temp)]
prediction_mat[1,temp]
mf_preds_all[1,temp]

which(train_mat[drug_set,c(which(target_adj_mat[temp,]>0),temp)]>=12.1)
which(train_mat[drug_set,c(which(target_adj_mat[temp,]>0),temp)]<12.1 & train_mat[drug_set,c(which(target_adj_mat[temp,]>0),temp)]>-1)
