.sparseVAR = function(j, y, p, intercept, y.spec) {
  colIndex = j[1]
  y.to.remove = which(!y.spec[colIndex,])
  z.to.remove = c()
  if(length(y.to.remove) > 0) {
    z.to.remove = as.vector(sapply(1:p, function(i) {
      if(intercept) (y.to.remove + 1) + ncol(y)*(i-1)
      else y.to.remove + ncol(y)*(i-1)
    }))
  }
  if(length(z.to.remove) > 0) {
    Z.reduced = Z[,-z.to.remove]
  } else {
    Z.reduced = Z
  }
       
  return (cv.glmnet(Z.reduced, j[-1]))
}

#' Sparse Vector Autoregression
#' Fit a vector autoregressive model with lasso penalty.
#' The VAR model is estimated using a multiresponse linear regression.
#' The sparse VAR fits multiple uniresponse linear regressions with lasso penalty
#' mclapply from multicore can be used to fit the individual uniresponse
#' linear regressions in parallel.  Note that mclapply is not available for windows
#' @param y
#' @param p
#' @param intercept
#' @param y.spec
#' @param numcore
SparseVAR = function(y, p, intercept = F,
                     y.spec=matrix(1,nrow=ncol(y),ncol=ncol(y)),
                     numcore=1, ...) {
  if(p < 1) stop("p must be a positive integer")
  if(!is.matrix(y)) {
    stop("y must be a matrix (not a data frame).  Consider using as.matrix(y)")
  }
  var.z = VAR.Z(y,p,intercept)
  Z = var.z$Z
  y.augmented = rbind(1:ncol(y),var.z$y.p)

  if(numcore == 1 | !require(multicore)) {
    var.lasso = apply(y.augmented, 2, .sparseVAR, y=y,p=p,
                      intercept=intercept, y.spec=y.spec)
  } else {
    y.augmented.list = c(unname(as.data.frame(y.augmented)))
    var.lasso = mclapply(y.augmented.list, .sparseVAR, y=y,p=p,
                         Z=Z, y.spec=y.spec, mc.cores=numcore, ...)
  }

  return(structure(list(
              model = var.lasso,
              var.z = var.z 
  ), class="fastVAR.SparseVAR"))
}

#' Coefficients of a SparseVAR model
#'
#' @param l1penalty The l1 penalty to be applied to the SparseVAR
#'   model.  
coef.fastVAR.SparseVAR = function(SparseVAR, l1penalty) {
  if(missing(l1penalty)) {
    B = lapply(SparseVAR$model, function(model) {
          as.vector(coef(model, model$lambda.min))
    })
  }
  else if(length(l1penalty) == 1) {
    B = lapply(SparseVAR$model, function(model) {
          as.vector(coef(model, l1penalty))
    })
  } else {
    B = rep(0, ncol(SparseVAR$var.z$Z))
    for (i in 1:length(l1penalty)) {
      B[i] = as.vector(coef(SparseVAR$model[[i]], l1penalty[i]))
    }
  }
  colnames(B) = colnames(SparseVAR$var.z$y.orig)
  rownames(B) = colnames(SparseVAR$var.z$Z)

  return (B)
}

predict.fastVAR.SparseVAR = function(SparseVAR, n.ahead=1, threshold, ...) {
  y.pred = matrix(nrow=n.ahead, ncol=ncol(SparseVAR$var.z$y.orig))
  colnames(y.pred) = colnames(SparseVAR$var.z$y.orig)
  for(i in 1:n.ahead) {
    if(SparseVAR$var.z$intercept) {
      Z.ahead = c(1,as.vector(t(SparseVAR$var.z$y.orig[
        ((nrow(SparseVAR$var.z$y.orig)):
        (nrow(SparseVAR$var.z$y.orig)-SparseVAR$var.z$p+1))
      ,])))
    } else {
      Z.ahead = as.vector(t(SparseVAR$var.z$y.orig[
        ((nrow(SparseVAR$var.z$y.orig)):
        (nrow(SparseVAR$var.z$y.orig)-SparseVAR$var.z$p+1))
      ,]))
    }
    y.ahead = Z.ahead %*% coef(SparseVAR, ...)
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
