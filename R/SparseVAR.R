.sparseVAR = function(j, y, p, Z, y.spec) {
  colIndex = j[1]
  y.to.remove = which(!y.spec[colIndex,])
  z.to.remove = c()
  if(length(y.to.remove) > 0) {
    z.to.remove = as.vector(sapply(1:p, function(i) {
      (y.to.remove + 1) + ncol(y)*(i-1)
    }))
  }
  if(length(z.to.remove) > 0) {
    Z.reduced = Z[,-z.to.remove]
  } else {
    Z.reduced = Z
  }
       
  return (cv.glmnet(Z.reduced[,-1], j[-1]))
}

#' Sparse Vector Autoregression
#'
#' Fit a vector autoregressive model with lasso penalty.
#' The VAR model is estimated using a multiresponse linear regression.
#' The sparse VAR fits multiple uniresponse linear regressions with lasso penalty.
#' mclapply from multicore can be used to fit the individual uniresponse
#' linear regressions in parallel.  Note that mclapply is not available for windows
#' @param y A matrix where each column represents an individual time series
#' @param p the number of lags to include in the design matrix
#' @param y.spec A binary matrix that can constrain the number of lagged predictor variables.  
#'   If y.spec[i][j] = 0, the ith time series in y will not be regressed on the jth
#'   time series of y, or any of its lags.
#' @param numcore number of cpu cores to use to parallelize this function
#' @examples 
#'   data(Canada)
#'   SparseVAR(Canada, 3)
#' @export
SparseVAR = function(y, freq=NULL, p
                     y.spec=matrix(1,nrow=ncol(y),ncol=ncol(y)),
                     numcore=1, ...) {
  if(p < 1) stop("p must be a positive integer")
  if(!is.matrix(y)) {
    stop("y must be a matrix (not a data frame).  Consider using as.matrix(y)")
  }
  y.seasons = deseason(y, freq)
  var.z = VAR.Z(y$remaining,p,intercept = T)
  Z = var.z$Z
  y.augmented = rbind(1:ncol(y),var.z$y.p)

  if(numcore == 1 | !require(multicore)) {
    var.lasso = apply(y.augmented, 2, .sparseVAR, y=y,p=p,
                      Z=Z, y.spec=y.spec)
  } else {
    y.augmented.list = c(unname(as.data.frame(y.augmented)))
    var.lasso = mclapply(y.augmented.list, .sparseVAR, y=y,p=p,
                         Z=Z, y.spec=y.spec, mc.cores=numcore, ...)
  }

  return(structure(list(
                        model = var.lasso,
                        var.z = var.z,
                        seasons = y.seasons),
                   class="fastVAR.SparseVAR"))
}

#' Coefficients of a SparseVAR model
#'
#' The underlying library, glmnet, computes the full path to the lasso.
#' This means it is computationally easy to compute the lasso solution
#' for any penalty term.  This function allows you to pass in the desired
#' l1 penalty and will return the coefficients
#' @param sparseVAR an object of class fastVAR.SparseVAR
#' @param l1penalty The l1 penalty to be applied to the SparseVAR
#'   model.
#' @export  
coef.fastVAR.SparseVAR = function(sparseVAR, l1penalty) {
  if (missing(l1penalty)) {
    B = data.frame(lapply(sparseVAR$model, function(model) {
          as.vector(coef(model, model$lambda.min))
    }))
  }
  else if (length(l1penalty) == 1) {
    B = data.frame(lapply(sparseVAR$model, function(model) {
          as.vector(coef(model, l1penalty))
    }))
  } else {
    B = matrix(0, nrow=ncol(sparseVARX$var.z$Z), ncol=ncol(sparseVARX$var.z$y.p))
    for (i in 1:length(l1penalty)) {
      B[,i] = as.vector(coef(sparseVAR$model[[i]], l1penalty[i]))
    }
  }
  colnames(B) = colnames(sparseVAR$var.z$y.orig)
  rownames(B) = colnames(sparseVAR$var.z$Z)

  return (as.matrix(B))
}

#' SparseVAR Predict
#'
#' Predict n steps ahead from a fastVAR.SparseVAR object
#' @param sparseVAR an object of class fastVAR.SparseVAR returned from SparseVAR
#' @param n.ahead number of steps to predict
#' @param threshold threshold prediction values to be greater than this value
#' @param ... extra parameters to pass into the coefficients method
#'   for objects of type fastVAR.VAR
#' @examples
#'   data(Canada)
#'   predict(SparseVAR(Canada, 3), 1)
#' @export
predict.fastVAR.SparseVAR = function(sparseVAR, n.ahead=1, threshold, ...) {
  y.pred = matrix(nrow=n.ahead, ncol=ncol(sparseVAR$var.z$y.orig))
  colnames(y.pred) = colnames(sparseVAR$var.z$y.orig)
  for (i in 1:n.ahead) {
    Z.ahead = c(1,as.vector(t(sparseVAR$var.z$y.orig[
        ((nrow(sparseVAR$var.z$y.orig)):
        (nrow(sparseVAR$var.z$y.orig)-sparseVAR$var.z$p+1))
      ,])))
    y.ahead = Z.ahead %*% coef(sparseVAR, ...)
    if (!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if(length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if (i == n.ahead) break
    sparseVAR$var.z$y.orig = rbind(sparseVAR$var.z$y.orig, y.ahead)
  }
  season = lastPeriod(sparseVAR$seasons)
  return (y.pred + season)
}
