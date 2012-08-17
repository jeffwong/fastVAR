#x is a multivariate time series where each col is a time series
#p is the order of the AR component
fastVAR = function(y, p=1, getdiag=T) {
  if(p < 1) {
    stop("p must be a positive integer")
  }
  varz = VARZ(y, p)
  model = lm(varz$y.p ~ varz$Z)
  if(any(is.na(model$coefficients))) {
    stop("Multivariate lm has invalid coefficients.  
      Check the rank of the design matrix")
  }
  if(getdiag) {
    return ( list (
      model = model,
      diag = VARdiag(varz$y.p, varz$Z, model$coefficients,
        varz$n, varz$T, (varz$k+1), p, varz$dof)
    ))
  } else {
    return (model)
  }
}

#j is the response vector from y.p
.VARlassocore = function(j, y, p, lambda, Z, y.spec) {
  colIndex = j[1]
  y.to.remove = which(!y.spec[colIndex,])
  z.to.remove = c()
  if(length(y.to.remove) > 0) {
    z.to.remove = as.vector(sapply(1:p, function(i) {
      y.to.remove + ncol(y)*(i-1)
    }))
  }
  if(length(z.to.remove) > 0) {
    Z.reduced = Z[,-z.to.remove]
  } else {
    Z.reduced = Z
  }
       
  if(is.null(lambda)) {
    j.lasso = cv.glmnet(Z.reduced, j[-1])
    j.lasso.s = j.lasso$lambda.min
  } else {
    j.lasso = glmnet(Z.reduced, j[-1])
    j.lasso.s =lambda[colIndex]
  }

  B = rep(0, ncol(Z)+1)  #make vector of coefficients full size
  B.reduced = rep(0, ncol(Z.reduced) + 1)

  j.coef = coef(j.lasso, j.lasso.s)
  B.reduced[(j.coef@i)+1] = j.coef@x
  if(length(z.to.remove) > 0) {
    B[-(z.to.remove+1)] = B.reduced
    return (B)
  }
  else {
    return (B.reduced)
  }
}

#Note that mclapply is not available for windows
VARlasso = function(y, p, lambda=NULL, y.spec=matrix(1,nrow=ncol(y),ncol=ncol(y)),
  getdiag=T, numcore=1, ...) {
  if(p < 1) stop("p must be a positive integer")
  if(!is.matrix(y)) {
    stop("y must be a matrix (not a data frame).  Consider using as.matrix(y)")
  }
  varz = VARZ(y,p)
  Z = varz$Z
  y.augmented = rbind(1:ncol(y),varz$y.p )

  if(numcore == 1 | !require(multicore)) {
    var.lasso = apply(y.augmented, 2, .VARlassocore, y=y,p=p, lambda=lambda,
      Z=Z, y.spec=y.spec)
  } else {
    y.augmented.list = c(unname(as.data.frame(y.augmented)))
    var.lasso = mclapply(y.augmented.list, .VARlassocore, y=y,p=p, lambda=lambda,
      Z=Z, y.spec=y.spec, mc.cores=numcore, ...)
    var.lasso = matrix(unlist(var.lasso), nrow=(varz$k+1), ncol=ncol(y))
    colnames(var.lasso) = colnames(y)
  }

  rownames(var.lasso) = c('intercept', colnames(Z))

  if(getdiag) {
    return( list (
      B = var.lasso,
      diag = VARdiag(varz$y.p, Z, var.lasso,
        varz$n, varz$T, (varz$k+1), p, varz$dof)
    ))
  } else {
    return (var.lasso)
  }
}

VARpredict = function(y, p, B, n.ahead=1, threshold) {
  if(p < 1) {
    stop("p must be a positive integer")
  }
  y.pred = matrix(nrow=n.ahead, ncol=ncol(y))
  colnames(y.pred) = colnames(y)
  for(i in 1:n.ahead) {
    Z.ahead = c(1,as.vector(t(y[
      ((nrow(y)):
      (nrow(y)-p+1))
    ,])))
    y.ahead = Z.ahead %*% B
    if(!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if(length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if(i == n.ahead) break
    y = rbind(y, y.ahead)
  }
  return (y.pred)
}

#Construct the design matrix for a VAR(p) model
VARZ = function(y, p) {
  if(p < 1) stop("p must be a positive integer")
      
  if(is.null(colnames(y))) {
    colnames(y) = sapply(1:ncol(y), function(j) {
      paste('y',j,sep='')
    })
  }

  n = ncol(y)
  T = nrow(y)
  k = n*p
  dof = T - p - k
        
  y.p = y[(1+p):T,]
  Z = do.call('cbind', (sapply(0:(p-1), function(k) {
    y[(p-k):(T-1-k),]
  }, simplify=F)))
  Z.names = colnames(Z)
  colnames(Z) = sapply(1:p, function(i) {
    startIndex = (i-1)*n + 1
    endIndex = i*n
    paste(Z.names[startIndex:endIndex], '.l', i, sep='')
  })

  return ( list (
    n=n,
    T=nrow(y.p),
    k=k,
    dof=dof,
    y.p=y.p,
    Z=Z
  ))
}

#Slightly buggy
#VARdiag should return a vector for BIC, not a scalar
VARstep = function(y, max.p) {
  bic = sapply(1:max.p, function(i) {
    BIC = fastVAR(y, i, getdiag=T)$diag$BIC
  })
}

VARdiag = function(y.p, Z, B, n, T, k, p, dof) {
  Z = cbind(1, Z)
  colnames(Z)[1] = 'intercept'
  ZTZ.inv = solve(t(Z) %*% Z)
  residual = y.p - Z %*% B
  residual.cov = t(residual) %*% residual  #residual covariance matrix
  sigma.hat = (1 / (T - k)) * 
    (residual.cov)
  cov.hat = kronecker(sigma.hat, ZTZ.inv)
  se = matrix(sqrt(diag(cov.hat)), nrow=nrow(B), ncol=ncol(B))
  tvalue = B / se
  pvalue = 2*(1 - pt(abs(tvalue),dof)) 
  rownames(se) = colnames(Z)
  rownames(tvalue) = colnames(Z)
  rownames(pvalue) = colnames(Z)
  colnames(se) = colnames(y.p)
  colnames(tvalue) = colnames(y.p)
  colnames(pvalue) = colnames(y.p)

  #Get model selection statistics
  sigma.tilde = 1 / T * residual.cov
  sigma.tilde.det = log(det(sigma.tilde))
  
  AIC = sigma.tilde.det + 2/T * p * n^2
  BIC = sigma.tilde.det + log(T)/T * p * n^2
  HQ = sigma.tilde.det + 2*log(log(T))/T * p * n^2

  return (list(
    se=se,
    tvalue=tvalue,
    pvalue=pvalue,
    residual=residual,
    cov=cov.hat,
    AIC=AIC,
    BIC=BIC,
    HQ=HQ
  ))
}
