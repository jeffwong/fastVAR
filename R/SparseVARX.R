.sparseVARX = function(j, y, x, p, b, intercept, y.spec, x.spec) {
  colIndex = j[1]
  x.to.remove = which(!x.spec[colIndex,])
  y.to.remove = which(!y.spec[colIndex,])
  np = ncol(y) * p
  z.x.to.remove = c()
  z.y.to.remove = c()
  if(length(x.to.remove) > 0) {
    z.x.to.remove = as.vector(sapply(1:b, function(ii) {
      if(intercept) (x.to.remove + 1) + np + ncol(x)*(ii-1)
      else x.to.remove + np + ncol(x)*(ii-1)
    }))
  }
  if(length(y.to.remove) > 0) {
    z.y.to.remove = as.vector(sapply(1:p, function(ii) {
      if(intercept) (y.to.remove + 1) + ncol(y)*(ii-1)
      else y.to.remove + ncol(y)*(ii-1)
    }))
  }
  z.to.remove = unlist(c(z.x.to.remove, z.y.to.remove))
  if(length(z.to.remove > 0)) {
    Z.reduced = Z[,-z.to.remove]
  }
  else {
    Z.reduced = Z
  }
       
  return (cv.glmnet(Z.reduced, j[-1]))
}

SparseVARX = function(y, x, p, b, intercept=F, 
  y.spec=matrix(1,nrow=ncol(y),ncol=ncol(y)), 
  x.spec=matrix(1,nrow=ncol(y),ncol=ncol(x)),
  numcore=1, ...) {
 
  if(p < 1) stop("p must be a positive integer")
  if(!is.matrix(y)) {
    stop("y must be a matrix (not a data frame).  Consider using as.matrix(y)")
  }
  var.z = VARX.Z(y,x,p,b,intercept)
  Z = var.z$Z
  y.augmented = rbind(1:ncol(y),var.z$y.p)
  
  if(numcore==1 | !require(multicore)) {
    var.lasso = apply(y.augmented, 2, .sparseVARX, y=y,x=x,p=p,b=b,
                      intercept=intercept, y.spec=y.spec, x.spec=x.spec)
  } else {
    y.augmented.list = c(unname(as.data.frame(y.augmented)))
    var.lasso = mclapply(y.augmented.list, .sparseVARX,
                         y=y,x=x,p=p,b=b,
                         intercept=intercept,
                         y.spec=y.spec, x.spec=x.spec,
                         mc.cores=numcore, ...)
  }

  return(structure (list(
              model = var.lasso,
              var.z = var.z 
  ), class="fastVAR.SparseVARX"))
}

coef.fastVAR.SparseVARX = function(SparseVARX, l1penalty) {
  if(missing(l1penalty)) {
    B = lapply(SparseVARX$model, function(model) {
          as.vector(coef(model, model$lambda.min))
    })
  }
  else if(length(l1penalty) == 1) {
    B = lapply(SparseVARX$model, function(model) {
          as.vector(coef(model, l1penalty))
    })
  } else {
    B = matrix(0, nrow(SparseVARX$var.z$Z), ncol(SparseVARX$var.z$Z))
    for (i in 1:length(l1penalty)) {
      B[,i] = as.vector(coef(SparseVARX$model[[i]], l1penalty[i]))
    }
  }
  colnames(B) = colnames(SparseVARX$var.z$y.orig)
  rownames(B) = colnames(SparseVARX$var.z$Z)

  return (B)
}

predict.fastVAR.SparseVARX = function(SparseVARX, xnew, n.ahead=1, threshold, ...) {
  y.pred = matrix(nrow=n.ahead, ncol=ncol(SparseVARX$var.z$y.orig))
  colnames(y.pred) = colnames(SparseVARX$var.z$y.orig)
  for(i in 1:n.ahead) {
    Z.ahead.y = as.vector(t(SparseVARX$var.z$y.orig[
      ((nrow(SparseVARX$var.z$y.orig)):
      (nrow(SparseVARX$var.z$y.orig)-SparseVARX$var.z$p+1))
    ,]))
    if(SparseVARX$var.z$b == 0) {
      Z.ahead.x = xnew[1,]
    } else {
      Z.ahead.x = as.vector(t(SparseVARX$var.z$x.orig[
        ((nrow(SparseVARX$var.z$x.orig)):
        (nrow(SparseVARX$var.z$x.orig)-SparseVARX$var.z$b+1))
      ,]))
    }
    if(SparseVARX$var.z$intercept) Z.ahead = c(1, Z.ahead.y, Z.ahead.x)
    else Z.ahead = c(Z.ahead.y, Z.ahead.x)
    y.ahead = Z.ahead %*% coef(SparseVARX, ...)
    if(!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if(length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if(i == n.ahead) break
    SparseVARX$var.z$y.orig = rbind(SparseVARX$var.z$y.orig, y.ahead)
    SparseVARX$var.z$x.orig = rbind(SparseVARX$var.z$x.orig, xnew[1,])
  }
  return (y.pred)
}
