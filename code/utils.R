CalcPrecision = function(pred,truth,k){ # Length smaller than k!
    k = min(k,length(pred))
    retrieved_ind = order(pred,decreasing = TRUE)[1:k]
    true_positives = sum(truth[retrieved_ind])
    return(true_positives/k)
}

CalcRecall = function(pred,truth,k){
    k = min(k,length(pred))
    relevant = sum(truth==1)
    retrieved_ind = order(pred,decreasing = TRUE)[1:k]
    true_positives = sum(truth[retrieved_ind])
    if (relevant>0)
        return(true_positives/relevant)
    else
        return(0)
}

CalcPR = function(pred,triplet,k,threshold = NULL,verbose = FALSE)
{
    if (!is.null(threshold))
    {
        if (length(threshold)==1)
            triplet[,3] = realToBinary(triplet[,3],threshold)
        else
            triplet[,3] = realToBinary(triplet[,3],threshold[triplet[,1]])
    }
    total_precision = 0
    total_recall = 0
    truth_and_prediction_ind = split(1:nrow(triplet),triplet[,1])
    number_rows = length(truth_and_prediction_ind)
    for (i in 1:number_rows){
        ind = truth_and_prediction_ind[[i]]
        ligand_i_truth = triplet[ind,3]
        ligand_i_prediction = pred[ind]
        
        temp_precision = CalcPrecision(ligand_i_prediction,ligand_i_truth,k)
        if (is.na(temp_precision)) browser()
        temp_recall = CalcRecall(ligand_i_prediction, ligand_i_truth,k)
        total_precision = total_precision + temp_precision
        total_recall = total_recall + temp_recall
        if (verbose)
        {
            cat("ligand : ",i," precision: ")
            cat(temp_precision," recall: ")
            cat(temp_recall,"\n")
        }
    }
    precision_avg = total_precision/number_rows
    recall_avg = total_recall/number_rows
    
    return(c(precision_avg,recall_avg))
}

CalcAUC = function(pred,label)
{
    pred = (pred-min(pred))/(max(pred)-min(pred))
    return(auc(roc(pred,as.factor(label))))
}

realToBinary = function(val,threshold)
{
    res = as.numeric(val>threshold)
    return(res)
}

hinge = function(x,threshold)
{
    res = max(0,x-threshold)
    return(res)
}

rmse = function(preds, labels) sqrt(mean((preds-labels)^2))

svds.cv = function(triplet,nfolds,m,n,k)
{
    # split indices
    N = nrow(triplet)
    ind = sample(1:N,N)
    cv.ind = vector(nfolds,mode='list')
    rmd = (1:N)%%nfolds
    rmd[which(rmd==0)] = nfolds
    
    for (i in 1:nfolds)
    {
        tmp = which(rmd==i)
        cv.ind[[i]] = ind[tmp]
    }
    
    # Start to seperately do it
    Result = vector(nfolds,mode='list')
    Rmse.Train = rep(0,nfolds)
    Rmse.Test = rep(0,nfolds)
    for (i in 1:nfolds)
    {
        cat('Starting Round',i,'\n')
        tr.ind = setdiff(1:N,cv.ind[[i]])
        te.ind = cv.ind[[i]]
        subtrip = triplet[tr.ind,]
        Mat = spMatrix(i = subtrip[,1],
                       j = subtrip[,2],
                       x = subtrip[,3],
                       nrow = m, ncol = n)
        Mat = as(Mat,'dgCMatrix')
        res = svds(Mat,k)
        Result[[i]] = res
        
        res = res$u %*% diag(res$d) %*% t(res$v)
        rmse.train = rmse(res[as.matrix(triplet[tr.ind,1:2])],
                          triplet[tr.ind,3])
        rmse.test = rmse(res[as.matrix(triplet[te.ind,1:2])],
                         triplet[te.ind,3])
        Rmse.Train[i] = rmse.train
        Rmse.Test[i] = rmse.test
        cat('Train RMSE',rmse.train,'Test RMSE',rmse.test,'\n\n')
    }
    cat('Total Evaluation:\n')
    cat('Train RMSE',mean(Rmse.Train),'+',sd(Rmse.Train),
        ', Test RMSE',mean(Rmse.Test),'+',sd(Rmse.Test))
    return(list(Result,cv.ind))
}

PlotPR = function(pred,truth)
{
    n = length(pred)
    retrieved = order(pred,decreasing = TRUE)
    relevant = sum(truth==1)
    prec.res = rep(0,n)
    rec.res = rep(0,n)
    for (k in 1:n)
    {
        retrieved_ind = retrieved[1:k]
        true_positives = sum(truth[retrieved_ind])
        rec.res[k] = true_positives/relevant
        prec.res[k] = true_positives/k
    }
    return(list(prec.res,rec.res))
}
