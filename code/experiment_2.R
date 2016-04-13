load('kiba_data.rda')

n_drugs = length(unique(dt_triplet[,1]))
n_targets = length(unique(dt_triplet[,2]))

dt_triplet[,3] = dt_triplet[,3]*-1
dt_triplet[,3] = dt_triplet[,3] - min(dt_triplet[,3]) 
hist(dt_triplet[,3])

n_folds = 5
n = nrow(dt_triplet)
folds = vector(n_folds, mode='list')
shuf = sample(n)
m = n %/% n_folds
for (i in 1:(n_folds-1)) {
  folds[[i]] = shuf[1:m]
  shuf = setdiff(shuf, folds[[i]])
}
folds[[n_folds]] = shuf
preds = rep(0,n)

adj_mat_pred = matrix(0, nrow = n_drugs, ncol = n_drugs)
for (i in 1:n_drugs){
  inds = which(simi_mat[i,]>0.9 & simi_mat[i,]<1)
  adj_mat_pred[i, inds] = 1
  adj_mat_pred[inds,i ] = 1
}

for (i in 1:n_folds){
  cat('fold ',i,'\n')
  
  dt_mat = matrix(-1, nrow = n_drugs, ncol = n_targets)
  dt_mat[dt_triplet[trind,c(1,2)]] = dt_triplet[trind,3]
  
  adj_mats_train = list()
  for (t1 in 1:n_targets){
    adj_mats_train[[t1]] = make_training_adj_mat_for_column(dt_mat, simi_mat, t1)
    #cat('column ',t,' ',length(which(adj_mats[[t]]>0))/2,' number of edges to train\n')
  }
  
  teind = folds[[i]]
  trind = setdiff(1:n,teind)
  
  for (t in 1:n_targets){
    cat('target ',t,'... ',length(which(dt_mat[,t]>0)),' observations\n')
    
    X = get_mf_cv(dt_triplet[trind,])
    
    cat('training crf.. \n')
    
    params = train_crf_by_row_3(train = dt_mat, adj_list = adj_mats_train, X = X, crf_iters = 200, eta = 0.01, t)
    
    cat('parameters: ', params[[1]], params[[2]],'\n')
    
    cat('getting MF predictions..\n')
    
    X = get_libmf_prediction(dt_mat, 400)
    rmse_all = sqrt(mean((dt_triplet[teind,3]-X[dt_triplet[teind,c(1:2)]])^2))
    inds = which(dt_triplet[teind,2] == t)
    preds_col = X[dt_triplet[teind[inds],c(1:2)]]
    label_col = dt_triplet[teind[inds], 3]
    mf_rmse = sqrt(mean((preds_col-label_col)^2))
    if (length(which(label_col>=12.1))>0){
      pred.obj = ROCR::prediction(preds_col, as.numeric(label_col>=12.1))
      auc.obj = ROCR::performance(pred.obj, measure = 'auc')
      prec.obj = ROCR::performance(pred.obj, measure = 'prec')
      rec.obj = ROCR::performance(pred.obj, measure = 'rec')
      prec.val = prec.obj@y.values[[1]]
      rec.val = rec.obj@y.values[[1]]
      
      if (is.na(prec.val[1])){
        prec.val[1]=1
      }
      func = approxfun(cbind(rec.val,prec.val), yleft = 1)
      #func = tryCatch(approxfun(cbind(rec.val,prec.val), yleft = 1), error=function(e) NA)
      mf_aupr.val = integrate(func, 0, 1, subdivisions = 1000L)$value
      mf_auc.val = auc.obj@y.values[[1]]
    } else{
      mf_aupr.val = 0
      mf_auc.val = 0
    }
    
    cat('making crf predictions..\n')
    
    crf_preds = make_crf_predictions_for_row(params[[1]], params[[2]], train = dt_mat, adj_mat_pred, X, t)
    
    preds_col = crf_preds[dt_triplet[teind[inds],1]]
    label_col = dt_triplet[teind[inds], 3]
    crf_rmse = sqrt(mean((dt_triplet[teind[inds], 3] - crf_preds[dt_triplet[teind[inds],1]])^2))
   
    if (length(which(label_col>=12.1))>0){
      pred.obj = ROCR::prediction(preds_col, as.numeric(label_col>=12.1))
      auc.obj = ROCR::performance(pred.obj, measure = 'auc')
      prec.obj = ROCR::performance(pred.obj, measure = 'prec')
      rec.obj = ROCR::performance(pred.obj, measure = 'rec')
      prec.val = prec.obj@y.values[[1]]
      rec.val = rec.obj@y.values[[1]]
      crf_auc.val = auc.obj@y.values[[1]]
      func = approxfun(cbind(rec.val,prec.val), yleft = 1)
      crf_aupr.val = integrate(func, 0, 1, subdivisions = 1000L)$value
    } else {
      crf_auc.val = 0
      crf_aupr.val = 0
    }
    
    cat('rmse (mf, crf): ',mf_rmse,', ',crf_rmse,'\n')
    cat('auc (mf, crf): ',mf_auc.val,', ',crf_auc.val,'\n')
    cat('aupr (mf, crf): ',mf_aupr.val,', ',crf_aupr.val,'\n') 
  }
}