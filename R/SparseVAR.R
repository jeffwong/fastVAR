.sparseVAR = function(j, y, p, lambda, Z, y.spec) {
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

#' Sparse Vector Autoregression
#' Fit a vector autoregressive model with lasso penalty.
#' The VAR model is estimated using a multiresponse linear regression.
#' The sparse VAR fits multiple uniresponse linear regressions with lasso penalty
#' mclapply from multicore can be used to fit the individual uniresponse
#' linear regressions in parallel.  Note that mclapply is not available for windows
#' @param y
#' @param p
#' @param lambda
#' @param y.spec
#' @param numcore
SparseVAR = function(y, p, lambda=NULL, y.spec=matrix(1,nrow=ncol(y),ncol=ncol(y)),
                     numcore=1, ...) {
  if(p < 1) stop("p must be a positive integer")
  if(!is.matrix(y)) {
    stop("y must be a matrix (not a data frame).  Consider using as.matrix(y)")
  }
  var.z = VAR.Z(y,p)
  Z = var.z$Z
  y.augmented = rbind(1:ncol(y),var.z$y.p )

  if(numcore == 1 | !require(multicore)) {
    var.lasso = apply(y.augmented, 2, .sparseVAR, y=y,p=p, lambda=lambda,
      Z=Z, y.spec=y.spec)
  } else {
    y.augmented.list = c(unname(as.data.frame(y.augmented)))
    var.lasso = mclapply(y.augmented.list, .sparseVAR, y=y,p=p, lambda=lambda,
      Z=Z, y.spec=y.spec, mc.cores=numcore, ...)
    var.lasso = matrix(unlist(var.lasso), nrow=(var.z$k+1), ncol=ncol(y))
    colnames(var.lasso) = colnames(y)
  }

  rownames(var.lasso) = c('intercept', colnames(Z))
  return(structure (list(
              coef = var.lasso,
              var.z = var.z 
  ), class="fastVAR.SparseVAR"))
}

predict.fastVAR.SparseVAR = function(SparseVAR, n.ahead=1, threshold) {
  y.pred = matrix(nrow=n.ahead, ncol=ncol(SparseVAR$var.z$y.orig))
  colnames(y.pred) = colnames(SparseVAR$var.z$y.orig)
  for(i in 1:n.ahead) {
    Z.ahead = c(1,as.vector(t(SparseVAR$var.z$y.orig[
      ((nrow(SparseVAR$var.z$y.orig)):
      (nrow(SparseVAR$var.z$y.orig)-SparseVAR$var.z$p+1))
    ,])))
    y.ahead = Z.ahead %*% SparseVAR$coef
    if(!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if(length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if(i == n.ahead) break
    SparseVAR$var.z$y.orig = rbind(SparseVAR$var.z$y.orig, y.ahead)
  }
  return (y.pred)
}
