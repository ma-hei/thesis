
# train = mix_dataset
train = dt_triplet
k = 5
n = nrow(train)
folds = vector(k, mode='list')
shuf = sample(n)
m = n %/% k
for (i in 1:(k-1)) {
  folds[[i]] = shuf[1:m]
  shuf = setdiff(shuf, folds[[i]])
}
folds[[k]] = shuf
preds = rep(0,n)

train.triplet = train[,1:3]
train.triplet[,1] = train.triplet[,1] - 1
train.triplet[,2] = train.triplet[,2] - 1

for (i in 1:k) {
  teind = folds[[i]]
  trind = setdiff(1:n,teind)
  res = libmf(train.triplet[trind,], m = n_drugs, n = n_targets, k = 20, cost = 0.01, lrate = 0.01,
              
              niter = 600, nthread = 1, nmf = FALSE, verbose = FALSE)
  
  P = res[[2]][[1]]
  Q = res[[2]][[2]]
  
  estM = P %*% t(Q)
  tmp.pred = estM[train[teind,1:2]]
  preds[teind] = tmp.pred
}

rmse.val = sqrt(mean((train[,3]-preds)^2))
pred.obj = ROCR::prediction(preds, as.numeric(train[,3]>cutoff))
auc.obj = ROCR::performance(pred.obj, measure = 'auc')
prec.obj = ROCR::performance(pred.obj, measure = 'prec')
rec.obj = ROCR::performance(pred.obj, measure = 'rec')
f.obj = ROCR::performance(pred.obj, measure = 'f')

# AUC
auc.val = auc.obj@y.values[[1]]

# Precision, Recall and F-score
prec.val = prec.obj@y.values[[1]]
rec.val = rec.obj@y.values[[1]]
f.val = f.obj@y.values[[1]]

# AUPR
func = approxfun(cbind(rec.val,prec.val), yleft = 1)
aupr.val = integrate(func, 0, 1, subdivisions = 1000L)$value

# plot
plot(rec.val, prec.val, type='l')
plot(preds, train[,3], pch = 20)
abline(h = cutoff, col = 2)
abline(v = cutoff, col = 2)
auc.val
# 0.9718129
aupr.val
# 0.9324747
rmse.val

