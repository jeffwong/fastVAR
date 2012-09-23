#' Vector Autoregression
#'
#' Fit a VAR model by creating the lagged design matrix
#' and fitting a multivariate response matrix to it.  Note that
#' the VAR response matrix omits the first p responses from the input
#' matrix Y.  This is because observations in Y are modeled by the
#' p previous values, so the first p observations cannot be modeled.
#' @param y A matrix where each column represents an individual time series
#' @param p the number of lags to include in the design matrix
#' @param intercept logical.  If true, include an intercept term in the model
#' @param weights weights applied to the multiresponse linear regression.
#'   Better predictions might come from weighting observations far in the past
#'   less so they impact the objective value less.  Either NULL,
#'   "exponential", or "linear"
#' @param l2penalty a ridge regression penalty, useful when the design matrix is 
#'   very wide, which may happen if p is large.
#' @param getdiag logical.  If true, return diagnostics
#' @export
VAR = function(y, p=1, intercept = T, weights=NULL, l2penalty=NULL, getdiag=T) {
  if(p < 1) {
    stop("p must be a positive integer")
  }
  var.z = VAR.Z(y, p, intercept)
  if(is.null(l2penalty) {
    if(!is.null(weights) & !is.vector(weights)) {
      weights = switch(weights,
                       exponential = exponentialWeights(var.z$Z, var.z$y.p)                 
                       linear = linearWeights(var.z$Z, var.z$y.p)
                      )
    }
    model = lm(var.z$y.p ~ var.z$Z, weights = weights)
    if(any(is.na(model$coefficients))) {
      stop("Multivariate lm has invalid coefficients.  
            Check the rank of the design matrix")
    }
    result = structure(list(
                            model = model,
                            var.z = var.z
                            ), class="fastVAR.VAR")
  } else {
    #Compute full path ridge solution
    z.svd = svd(var.z$Z)
    ridgePath = function(l2penalty) {
        z.svd$v %*% diag(1 / (z.svd$d^2 + l2penalty)) %*% diag(z.svd$d) %*% t(z.svd$u) %*% var.z$y.p
    }
    ridge.coef = ridgePath(l2penalty)
    result = structure(list(
                            model = structure(ridgePath, class="fastVAR.RidgePath"),
                            var.z = var.z), class="fastVAR.VAR")
  }

  if(getdiag) result$diag = VAR.diag(result)

  return (result)
}

coef.fastVAR.VAR = function(VAR, ...) {
  coef(VAR$model, ...)
}

#' VAR Predict
#'
#' Predict n steps ahead from a fastVAR.VAR object
#' @param VAR an object of class fastVAR.VAR returned from from VAR
#' @param n.ahead number of steps to predict
#' @param threshold threshold prediction values to be greater than this value
#' @param ... extra parameters to pass into the coefficients method
#'   for objects of type fastVAR.VAR
predict.fastVAR.VAR = function(VAR, n.ahead=1, threshold, ...) {
  y.pred = matrix(nrow=n.ahead, ncol=ncol(VAR$var.z$y.orig))
  colnames(y.pred) = colnames(VAR$var.z$y.orig)
  for(i in 1:n.ahead) {
    if(VAR$var.z$intercept) {
      Z.ahead = c(1,as.vector(t(VAR$var.z$y.orig[
        ((nrow(VAR$var.z$y.orig)):
        (nrow(VAR$var.z$y.orig)-VAR$var.z$p+1))
      ,])))
    } else {
      Z.ahead = as.vector(t(VAR$var.z$y.orig[
        ((nrow(VAR$var.z$y.orig)):
        (nrow(VAR$var.z$y.orig)-VAR$var.z$p+1))
      ,]))
    }
    y.ahead = Z.ahead %*% coef(VAR, ...)
    if(!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if(length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if(i == n.ahead) break
    VAR$var.z$y.orig = rbind(VAR$var.z$y.orig, y.ahead)
  }
  return (y.pred)
}
