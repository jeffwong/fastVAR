.sparseVARX = function(j, y, x, p, b, lambda, Z, y.spec, x.spec) {
  colIndex = j[1]
  x.to.remove = which(!x.spec[colIndex,])
  y.to.remove = which(!y.spec[colIndex,])
  np = ncol(y) * p
  z.x.to.remove = c()
  z.y.to.remove = c()
  if(length(x.to.remove) > 0) {
    z.x.to.remove = as.vector(sapply(1:b, function(ii) {
      x.to.remove + np + ncol(x)*(ii-1)
    }))
  }
  if(length(y.to.remove) > 0) {
    z.y.to.remove = as.vector(sapply(1:p, function(ii) {
      y.to.remove + ncol(y)*(ii-1)
    }))
  }
  z.to.remove = unlist(c(z.x.to.remove, z.y.to.remove))
  if(length(z.to.remove > 0)) {
    Z.reduced = Z[,-z.to.remove]
  }
  else {
    Z.reduced = Z
  }
       
  if(is.null(lambda)) {
    j.lasso = cv.glmnet(Z.reduced, j[-1])
    j.lasso.s = j.lasso$lambda.min
  } else {
    j.lasso = glmnet(Z.reduced, j[-1])
    j.lasso.s = lambda[colIndex]
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

SparseVARX = function(y, x, p, b, lambda=NULL, 
  y.spec=matrix(1,nrow=ncol(y),ncol=ncol(y)), 
  x.spec=matrix(1,nrow=ncol(y),ncol=ncol(x)),
  numcore=1, ...) {
 
  if(p < 1) stop("p must be a positive integer")
  if(!is.matrix(y)) {
    stop("y must be a matrix (not a data frame).  Consider using as.matrix(y)")
  }
  varx.z = VARX.Z(y,x,p,b)
  Z = varx.z$Z
  y.augmented = rbind(1:ncol(y),varx.z$y.p)
  
  if(numcore==1 | !require(multicore)) {
    var.lasso = apply(y.augmented, 2, .sparseVARX, y=y,x=x,p=p,b=b,
      lambda=lambda, Z=Z, y.spec=y.spec, x.spec=x.spec)
  } else {
    y.augmented.list = c(unname(as.data.frame(y.augmented)))
    var.lasso = mclapply(y.augmented.list, .sparseVARX, y=y,x=x,p=p,b=b,
      lambda=lambda, Z=Z, y.spec=y.spec, x.spec=x.spec, mc.cores=numcore, ...)
    var.lasso = matrix(unlist(var.lasso), nrow=(varx.z$k+1), ncol=ncol(y))
    colnames(var.lasso) = colnames(y)
  }
  rownames(var.lasso) = c('intercept', colnames(Z))

  return(structure (list(
              coef = var.lasso,
              varx.z = varx.z 
  ), class="fastVAR.SparseVARX"))
}

predict.fastVAR.SparseVARX = function(SparseVARX, xnew, n.ahead=1, threshold) {
  y.pred = matrix(nrow=n.ahead, ncol=ncol(SparseVARX$varx.z$y.orig))
  colnames(y.pred) = colnames(SparseVARX$varx.z$y.orig)
  for(i in 1:n.ahead) {
    Z.ahead.y = as.vector(t(SparseVARX$varx.z$y.orig[
      ((nrow(SparseVARX$varx.z$y.orig)):
      (nrow(SparseVARX$varx.z$y.orig)-SparseVARX$varx.z$p+1))
    ,]))
    if(SparseVARX$varx.z$b == 0) {
      Z.ahead.x = xnew[1,]
    } else {
      Z.ahead.x = as.vector(t(SparseVARX$varx.z$x.orig[
        ((nrow(SparseVARX$varx.z$x.orig)):
        (nrow(SparseVARX$varx.z$x.orig)-SparseVARX$varx.z$b+1))
      ,]))
    }
    Z.ahead = c(1, Z.ahead.y, Z.ahead.x)
    y.ahead = Z.ahead %*% SparseVARX$coef
    if(!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if(length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
      y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if(i == n.ahead) break
    SparseVARX$varx.z$y.orig = rbind(SparseVARX$varx.z$y.orig, y.ahead)
    SparseVARX$varx.z$x.orig = rbind(SparseVARX$varx.z$x.orig, xnew[1,])
  }
  return (y.pred)
}
