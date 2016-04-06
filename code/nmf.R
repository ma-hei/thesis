# source('utils.R')

# Home brew nmf 
calcObj = function(triplet,p,q,bu,bi,mu,lambda,bias = TRUE, interaction = TRUE)
{
    N = nrow(triplet)
    r = triplet[,3]
    ans = 0
    for (ui in 1:N)
    {
        u = triplet[ui,1]
        i = triplet[ui,2]
        if (bias && interaction)
        {
            ans = ans+(r[ui]-sum(p[u,]*q[i,])-bu[u]-bi[i]-mu)^2+
                lambda*sum(p[u,]^2+q[i,]^2+bu[u]^2+bi[i]^2)
        } else if (interaction) {
            ans = ans+(r[ui]-sum(p[u,]*q[i,]))^2+lambda*sum(p[u,]^2+q[i,]^2)
        } else if (bias) {
            ans = ans+(r[ui]-bu[u]-bi[i]-mu)^2+lambda*sum(bu[u]^2+bi[i]^2)
        } else {
            ans = ans+(r[ui]-mu)^2
        }
    }
    return(ans)
}

CalcRMSE = function(res,triplet,bias = TRUE, interaction = TRUE)
{
    estr = nmf.predict(res,triplet[,1:2],
                       bias = bias,interaction = interaction)
    r = triplet[,3]
    return(rmse(estr,r))
}

CalcPrecRec = function(res,triplet,k = 5, threshold = 2, bias = TRUE, interaction = TRUE)
{
    estr = nmf.predict(res,triplet[,1:2],
                       bias = bias,interaction = interaction)
    return(CalcPR(estr, triplet, k, threshold = threshold))
}

# triplet: the data
# m: the number of rows
# n: the number of cols
# k,lambda,gamma: parameters for the model
# tol,maxiter: controlling the iteration
nmf = function(triplet,m,n,k,lambda,gamma,tol,maxiter,valid = NULL, threshold = 2,
               bias = TRUE, interaction = TRUE, alpha = 0, eps = 1, regular = NULL,
               matInput = NULL, gradInput = NULL, verbose = TRUE)
{
    r = triplet[,3]
    if (is.null(matInput)) {
        p = matrix(rnorm(m*k,0,0.01),m,k)
        q = matrix(rnorm(n*k,0,0.01),n,k)
        bu = rnorm(m,0,0.01)
        bi = rnorm(n,0,0.01)
        mu = rnorm(1,0,0.01)
    } else {
        p = matInput$p
        q = matInput$q
        bu = matInput$bu
        bi = matInput$bi
        mu = matInput$mu
    }
    
    # Ada Gradient
    if (is.null(gradInput)) {
        pG = matrix(1e-8,m,k)
        qG = matrix(1e-8,n,k)
        buG = rep(1e-8,m)
        biG = rep(1e-8,n)
        muG = 1e-8
    } else {
        pG = gradInput$pG
        qG = gradInput$qG
        buG = gradInput$buG
        biG = gradInput$biG
        muG = gradInput$muG
    }
    
    # Initialization
    err = Inf
    lastObj = Inf
    nStep = 1
    N = nrow(triplet)
    if (length(threshold)!=m)
    {
        if (length(threshold)==1)
            threshold = rep(threshold,m)
        else
            stop('use threshold of length 1 or m')
    }
    
    #Starting Training
    while (err>tol && nStep<=maxiter)
    {
        shuf = sample(1:N,N)
        # try to find the way to balance
#         ind1 = which(triplet[,3] == 1)
#         ind0 = which(triplet[,3] == 0)
#         nn = length(ind1)
#         ind0 = sample(ind0, nn)
#         shuf = sample(c(ind0,ind1), 2*nn)
        # shuf = 1:N
#         one.weight = sum(triplet[,3]==0)/nrow(triplet)
#         zero.weight = 1 - one.weight
        pgrad = matrix(0,m,k)
        qgrad = matrix(0,n,k)
        bugrad = rep(0,m)
        bigrad = rep(0,n)
        mugrad = 0
        for (ui in shuf)
        {
            # get index
            u = triplet[ui,1]
            i = triplet[ui,2]
            
            # calculate error
            if (bias && interaction) {
                eui = sum(p[u,]*q[i,])+bu[u]+bi[i]+mu
            } else if (interaction) {
                eui = sum(p[u,]*q[i,])
            } else if (bias) {
                eui = bu[u]+bi[i]+mu
            } else {
                eui = mu
            }
            if (is.null(regular))
                eui = r[ui]-eui
            else if (regular=='linear')
                eui = (realToBinary(r[ui],threshold[u])-eui)*(1+alpha*hinge(r[ui],threshold[u]))#eui = (r[ui]-eui)*(1+alpha*r[ui])
            else if (regular=='log')
                eui = (realToBinary(r[ui],threshold[u])-eui)*(1+alpha*log(1+abs(hinge(r[ui],threshold[u])/eps)))#eui = (r[ui]-eui)*(1+alpha*log(1+abs(r[ui])/eps))
            else
                stop('No regularization found')
#            eui = 1/(1+exp(-eui)) - r[ui]
#             if (r[ui] == 0) {
#                 eui = eui*zero.weight
#             } else {
#                 eui = eui*one.weight
#             }
            
            # calculate grad
            pgrad[u,] = pgrad[u,]+eui*q[i,]-lambda*p[u,]
            qgrad[i,] = qgrad[i,]+eui*p[u,]-lambda*q[i,]
            bugrad[u] = bugrad[u]+eui-lambda*bu[u]
            bigrad[i] = bigrad[i]+eui-lambda*bi[i]
            mugrad = mugrad+eui
            
            # Update the Grad matrix
#             pG[u,] = pG[u,] + pgrad^2
#             qG[i,] = qG[i,] + qgrad^2
#             buG[u] = buG[u] + bugrad^2
#             biG[i] = biG[i] + bigrad^2
#             muG = muG + mugrad^2
            
            # update
#             p[u,] = p[u,]+gamma*pgrad/sqrt(pG[u,])
#             q[i,] = q[i,]+gamma*qgrad/sqrt(qG[i,])
#             bu[u] = bu[u]+gamma*bugrad/sqrt(buG[u])
#             bi[i] = bi[i]+gamma*bigrad/sqrt(biG[i])
#             mu = mu+gamma*mugrad/sqrt(muG)
            
            # cat(calcObj(triplet,p,q,lambda),'\n')
        }
        # update the grad matrix
        pG = pG + pgrad^2
        qG = qG + qgrad^2
        buG = buG + bugrad^2
        biG = biG + bigrad^2
        muG = muG + mugrad^2
        # update
        p = p+gamma*pgrad/sqrt(pG)
        q = q+gamma*qgrad/sqrt(qG)
        bu = bu+gamma*bugrad/sqrt(buG)
        bi = bi+gamma*bigrad/sqrt(biG)
        mu = mu+gamma*mugrad/sqrt(muG)


        obj = calcObj(triplet,p,q,bu,bi,mu,lambda,
                      bias = bias, interaction = interaction)
        if (verbose)
        {
            cat(nStep,obj,'\n')
            if (is.null(valid))
            {
                if (is.null(regular))
                {
                    rmse.step = CalcRMSE(list(p = p,q = q,bu = bu,bi = bi,mu = mu),
                                         triplet, bias = bias, interaction = interaction)
                    cat('RMSE',rmse.step,'\t\t')
                }
                precrec.step = CalcPrecRec(list(p = p,q = q,bu = bu,bi = bi,mu = mu),
                                           triplet, 5, threshold,
                                           bias = bias, interaction = interaction)
                cat('Precision',precrec.step[1],'Recall',precrec.step[2],'\n')
            }
            else
            {
                if (is.null(regular))
                {
                    rmse.train = CalcRMSE(list(p = p,q = q,bu = bu,bi = bi,mu = mu),
                                          triplet,bias = bias, interaction = interaction)
                    rmse.valid = CalcRMSE(list(p = p,q = q,bu = bu,bi = bi,mu = mu),
                                          valid,bias = bias, interaction = interaction)
                }
                precrec.train = CalcPrecRec(list(p = p,q = q,bu = bu,bi = bi,mu = mu),
                                            triplet, 5, threshold,
                                            bias = bias, interaction = interaction)
                precrec.valid = CalcPrecRec(list(p = p,q = q,bu = bu,bi = bi,mu = mu),
                                            valid, 5, threshold,
                                            bias = bias, interaction = interaction)
                cat('Train ')
                if (is.null(regular))
                {
                    cat('RMSE',rmse.train,'\t\t')
                }
                cat('Precision',precrec.train[1],'Recall',precrec.train[2],'\n')
                cat('Test ')
                if (is.null(regular))
                {
                    cat('Valid RMSE',rmse.valid,'\t\t')
                }
                cat('Precision',precrec.valid[1],'Recall',precrec.valid[2],'\n')
            }
        }
        err = abs(lastObj-obj)
        nStep = nStep+1
        lastObj = obj
    }
    return(list(p = p,q = q,bu = bu,bi = bi,mu = mu,
                pG = pG, qG = qG, buG = buG, biG = biG, muG = muG))
}

nmf.predict = function(res,index,bias = TRUE,interaction = TRUE, returnMat = FALSE)
{
    p = res$p
    q = res$q
    bu = res$bu
    bi = res$bi
    mu = res$mu
    m = length(bu)
    n = length(bi)
    N = nrow(index)
    
    mat = matrix(0,m,n)
    
    if (bias && interaction) {
        mat = p %*% t(q)
        offset = rep(bu,length(bi)) + rep(bi,each = length(bu))
        offset = matrix(offset,length(bu),length(bi))
        mat = mat+offset+mu
    } else if (interaction) {
        mat = p %*% t(q)
    } else if (bias) {
        offset = rep(bu,length(bi)) + rep(bi,each = length(bu))
        offset = matrix(offset,length(bu),length(bi))
        mat = mat+offset+mu
    } else {
        mat = mat+mu
    }
    if (returnMat)
        return(mat)
    pred = mat[as.matrix(index[,1:2])]
    return(pred)
}

# The cross validation function
# The only added parameter is nfolds
nmf.cv = function(triplet,nfolds,m,n,k,lambda,gamma,tol,maxiter, threshold = 2,
                  bias = TRUE, interaction = TRUE, alpha = 0, eps = 1, 
                  regular = NULL, prediction = FALSE)
{
    timestp = proc.time()
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
    
    preds = rep(0,nrow(triplet))
    matInput = vector(nfolds,mode='list')
    gradInput = vector(nfolds,mode='list')
    for (iter in 1:maxiter)
    {
        cat('\nStarting Iteration',iter,'\n')
        Rmse.Train = rep(0,nfolds)
        Rmse.Test = rep(0,nfolds)
        PrecRec.Train = matrix(0,nfolds,2)
        PrecRec.Test = matrix(0,nfolds,2)
        cat('Fold  ')
        for (i in 1:nfolds)
        {
            cat(i,' ')
            tr.ind = setdiff(1:N,cv.ind[[i]])
            te.ind = cv.ind[[i]]
            subtrip = triplet[tr.ind,]
            res = nmf(subtrip, m = m, n = n, k = k,
                      lambda = lambda, gamma = gamma, tol = tol, maxiter = 1,
                      threshold = threshold, valid = triplet[te.ind,],bias = bias, 
                      interaction = interaction, alpha = alpha, regular = regular,
                      matInput = matInput[[i]], gradInput = gradInput[[i]], verbose = FALSE)
            matInput[[i]] = res[1:5]
            gradInput[[i]] = res[6:10]
            Result[[i]] = res
            if (prediction)
                preds[te.ind] = nmf.predict(res,triplet[te.ind,1:2],
                                            bias = bias,interaction = interaction)
            rmse.train = CalcRMSE(res,triplet[tr.ind,],bias = bias,
                                  interaction = interaction)
            rmse.test = CalcRMSE(res,triplet[te.ind,],bias = bias,
                                 interaction = interaction)
            precrec.train = CalcPrecRec(res,triplet[tr.ind,],5,threshold,bias = bias,
                                        interaction = interaction)
            precrec.test = CalcPrecRec(res,triplet[te.ind,],5,threshold,bias = bias,
                                       interaction = interaction)
            Rmse.Train[i] = rmse.train
            Rmse.Test[i] = rmse.test
            PrecRec.Train[i,] = precrec.train
            PrecRec.Test[i,] = precrec.test
        }
        cat('\n')
        if (is.null(regular))
        {
            cat('Training RMSE',mean(Rmse.Train),'+',sd(Rmse.Train),'\t\t')
        }
        cat('Training Precision',mean(PrecRec.Train[,1]),'+',sd(PrecRec.Train[,1]),
            'Recall',mean(PrecRec.Train[,2]),'+',sd(PrecRec.Train[,2]),'\n')
        if (is.null(regular))
        {
            cat('Test RMSE',mean(Rmse.Test),'+',sd(Rmse.Test),'\t\t')
        }
        cat('Test Precision',mean(PrecRec.Test[,1]),'+',sd(PrecRec.Test[,1]),
            'Recall',mean(PrecRec.Test[,2]),'+',sd(PrecRec.Test[,2]),'\n')
    }
    cat('Time',(proc.time()-timestp)[3],'\n')

    if (prediction)
        return(list(Result,cv.ind,preds))
    else
        return(list(Result,cv.ind))
}















