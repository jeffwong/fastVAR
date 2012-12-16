.sparseVARX = function(j, y, x, p, b, Z, y.spec, x.spec) {
  colIndex = j[1]
  x.to.remove = which(!x.spec[colIndex,])
  y.to.remove = which(!y.spec[colIndex,])
  np = ncol(y) * p
  z.x.to.remove = c()
  z.y.to.remove = c()
  if(length(x.to.remove) > 0) {
    z.x.to.remove = as.vector(sapply(1:b, function(ii) {
      (x.to.remove + 1) + np + ncol(x)*(ii-1)
    }))
  }
  if(length(y.to.remove) > 0) {
    z.y.to.remove = as.vector(sapply(1:p, function(ii) {
      (y.to.remove + 1) + ncol(y)*(ii-1)
    }))
  }
  z.to.remove = unlist(c(z.x.to.remove, z.y.to.remove))
  if(length(z.to.remove > 0)) {
    Z.reduced = Z[,-z.to.remove]
  }
  else {
    Z.reduced = Z
  }
       
  return (cv.glmnet(Z.reduced[,-1], j[-1]))
}

#' Sparse Vector Autoregression with Exogenous Inputs
#'
#' Fit a VARX model with lasso penalty.
#' The VAR model is estimated using a multiresponse linear regression.
#' The sparse VAR fits multiple uniresponse linear regressions with lasso penalty.
#' mclapply from multicore can be used to fit the individual uniresponse
#' linear regressions in parallel.  Note that mclapply is not available for windows
#' @param y A matrix of endogenous variables where each column represents an individual time series
#' @param x A matrix of exogenous variables where each column represents an individual time series
#' @param p the number of lags of Y to include in the design matrix
#' @param b the number of lags to X include in the design matrix
#' @param y.spec A binary matrix that can constrain the number of lagged predictor variables.  
#'   If y.spec[i][j] = 0, the ith time series in y will not be regressed on the jth
#'   time series of y, or any of its lags.
#' @param x.spec A binary matrix that can constrain the number of lagged exogenous variables.  
#'   If x.spec[i][j] = 0, the ith time series in y will not be regressed on the jth
#'   time series of x, or any of its lags.
#' @param numcore number of cpu cores to use to parallelize this function
#' @examples
#'   data(Canada)
#'   x = matrix(rnorm(84*4), 84, 4)
#'   SparseVARX(Canada, x = x, p = 3, b = 2)
#' @export
SparseVARX = function(y, freq=rep(NA,ncol(y)), x, p, b, 
  y.spec=matrix(1,nrow=ncol(y),ncol=ncol(y)), 
  x.spec=matrix(1,nrow=ncol(y),ncol=ncol(x)),
  numcore=1, ...) {
 
  if (p < 1) stop("p must be a positive integer")
  if (!is.matrix(y)) {
    stop("y must be a matrix (not a data frame).  Consider using as.matrix(y)")
  }
  if (missing(x)) {
    return (SparseVAR(y, freq, p, y.spec, numcore))
  }

  y.seasons = deseason(y, freq)
  var.z = VARX.Z(y.seasons$remaining,x,p,b,intercept=T)
  Z = var.z$Z
  y.augmented = rbind(1:ncol(y),var.z$y.p)
  
  if (numcore==1) {
    var.lasso = apply(y.augmented, 2, .sparseVARX, y=y,x=x,p=p,b=b,
                      Z=Z, y.spec=y.spec, x.spec=x.spec)
  } else {
    y.augmented.list = c(unname(as.data.frame(y.augmented)))
    var.lasso = mclapply(y.augmented.list, .sparseVARX,
                         y=y,x=x,p=p,b=b,
                         Z=Z,
                         y.spec=y.spec, x.spec=x.spec,
                         mc.cores=numcore, ...)
  }

  return(structure (list(
                         model = var.lasso,
                         var.z = var.z,
                         seasons = y.seasons),
                    class="fastVAR.SparseVARX"))
}

#' Coefficients of a SparseVARX model
#'
#' The underlying library, glmnet, computes the full path to the lasso.
#' This means it is computationally easy to compute the lasso solution
#' for any penalty term.  This function allows you to pass in the desired
#' l1 penalty and will return the coefficients
#' @param sparseVARX an object of class fastVAR.SparseVARX
#' @param l1penalty The l1 penalty to be applied to the SparseVARX
#'   model.  
coef.fastVAR.SparseVARX = function(sparseVARX, l1penalty) {
  if (missing(l1penalty)) {
    B = data.frame(lapply(sparseVARX$model, function(model) {
          as.vector(coef(model, model$lambda.min))
    }))
  }
  else if (length(l1penalty) == 1) {
    B = data.frame(lapply(sparseVARX$model, function(model) {
          as.vector(coef(model, l1penalty))
    }))
  } else {
    B = matrix(0, nrow=ncol(sparseVARX$var.z$Z), ncol=ncol(sparseVARX$var.z$y.p))
    for (i in 1:length(l1penalty)) {
      B[,i] = as.vector(coef(sparseVARX$model[[i]], l1penalty[i]))
    }
  }
  colnames(B) = colnames(sparseVARX$var.z$y.orig)
  rownames(B) = colnames(sparseVARX$var.z$Z)

  return (as.matrix(B))
}

#' SparseVARX Predict
#'
#' Predict n steps ahead from a fastVAR.SparseVARX object
#' @param sparseVARX an object of class fastVAR.SparseVARX returned from SparseVARX
#' @param xnew a matrix of future values for the exogenous inputs.  Should contain
#'   n.ahead rows
#' @param n.ahead number of steps to predict
#' @param threshold threshold prediction values to be greater than this value
#' @param ... extra parameters to pass into the coefficients method
#'   for objects of type fastVAR.VAR
#' @examples
#'   data(Canada)
#'   x = matrix(rnorm(84*4), 84, 4)
#'   predict(SparseVARX(Canada, x = x, p = 3, b = 2), xnew=matrix(rnorm(2*4),2,4), n.ahead=2)
#' @export
predict.fastVAR.SparseVARX = function(sparseVARX, xnew, n.ahead=1, threshold, ...) {
  if (nrow(xnew) != n.ahead) stop("xnew should have n.ahead rows")
  y.pred = matrix(nrow=n.ahead, ncol=ncol(sparseVARX$var.z$y.orig))
  colnames(y.pred) = colnames(sparseVARX$var.z$y.orig)
  for (i in 1:n.ahead) {
    Z.ahead.y = as.vector(t(sparseVARX$var.z$y.orig[
      ((nrow(sparseVARX$var.z$y.orig)):
      (nrow(sparseVARX$var.z$y.orig)-sparseVARX$var.z$p+1))
    ,]))
    if(sparseVARX$var.z$b == 0) {
      Z.ahead.x = xnew[1,]
    } else {
      Z.ahead.x = as.vector(t(sparseVARX$var.z$x.orig[
        ((nrow(sparseVARX$var.z$x.orig)):
        (nrow(sparseVARX$var.z$x.orig)-sparseVARX$var.z$b+1))
      ,]))
    }
    Z.ahead = c(1, Z.ahead.y, Z.ahead.x)
    y.ahead = Z.ahead %*% coef(sparseVARX, ...)
    if(!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if(length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if (i == n.ahead) break
    sparseVARX$var.z$y.orig = rbind(sparseVARX$var.z$y.orig, y.ahead)
    sparseVARX$var.z$x.orig = rbind(sparseVARX$var.z$x.orig, xnew[1,])
  }
  freq = sparseVARX$seasons$freq
  freq.indices = which(!is.na(sparseVARX$seasons$freq))
  if (length(freq.indices) > 0) {
    lastSeason = lastPeriod(sparseVARX$seasons) #returns a list
    y.pred.seasonal = sapply(freq.indices, function(i) {
      season.start = periodIndex(freq[i], nrow(sparseVARX$var.z$y.orig + 1))
      season.end = season.start + n.ahead - 1
      rep(lastSeason[[i]], ceiling(n.ahead / freq[i]))[season.start : season.end]
    })
    y.pred[,freq.indices] = y.pred[,freq.indices] + y.pred.seasonal
    return (y.pred)
  }
  else return (y.pred)
}
