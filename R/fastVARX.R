#' Vector Autoregression with Exogenous Inputs
#'
#' Fit a VAR model by creating the lagged design matrix
#' and fitting a multivariate response matrix to it.  Note that
#' the VARX response matrix omits the first max(p,b) responses from the input
#' matrix Y.  This is because observations in Y are modeled by the
#' max(p,b) previous values, so the first max(p,b) observations cannot be modeled.
#' @param y A matrix of endogenous variables where each column represents an individual time series
#' @param x A matrix of exogenous variables where each column represents an individual time series
#' @param p the number of lags of Y to include in the design matrix
#' @param b the number of lags to X include in the design matrix
#' @param intercept logical.  If true include the intercept term in the model
#' @param weights weights applied to the multiresponse linear regression.
#'   Better predictions might come from weighting observations far in the past
#'   less so they impact the objective value less.  Either NULL,
#'   "exponential", or "linear"
#' @param l2penalty a ridge regression penalty, useful when the design matrix is 
#'   very wide, which may happen if p is large.
#' @param getdiag logical.  If true, return diagnostics
#' @export
VARX = function(y, x, p=1, b=1, intercept=T, weights=NULL, l2penalty=NULL, getdiag=T) {
  if(p < 1) stop("p must be a positive integer")
  if(missing(x)) {
    return (VAR(y, p, weights, intercept, l2penalty, getdiag))
  }
  var.z = VARX.Z(y, x, p, b, intercept)
  if(is.null(l2penalty)) {
    if(!is.null(weights) & !is.vector(weights)) {
      weights = switch(weights,
                       exponential = exponentialWeights(var.z$Z, var.z$y.p),
                       linear = linearWeights(var.z$Z, var.z$y.p))
    }
    model = lm(var.z$y.p ~ var.z$Z)
    if(sum(is.na(model$coefficients)) > 0) {
      stop("Multivariate lm has invalid coefficients.
            Check the rank of the design matrix")
    }
    result = structure(list(
                            model = model,
                            var.z = var.z
                            ), class="fastVAR.VARX")
 
  } else {
    #Compute full path ridge solution
    z.svd = svd(var.z$Z)
    ridgePath = function(l2penalty) {
        z.svd$v %*% diag(1 / (z.svd$d^2 + l2penalty)) %*% diag(z.svd$d) %*% t(z.svd$u) %*% var.z$y.p
    }
    result = structure(list(
                            model = structure(ridgePath, class="fastVAR.RidgePath"),
                            var.z = var.z), class="fastVAR.VARX")
  }
  if(getdiag) result$diag = VAR.diag(result)

  return (result) 
}

coef.fastVAR.VARX = function(VARX, ...) {
  coef(VARX$model, ...)
}

#' VARX Predict
#'
#' Predict n steps ahead from a fastVAR.VARX object
#' @param VARX an object of class fastVAR.VARX returned from VARX
#' @param xnew new values for the exogenous variables
#' @param n.ahead number of steps to predict
#' @param threshold threshold prediction values to be greater than this value
#' @param ... extra parameters to pass into the coefficients method
#'   for objects of type fastVAR.VARX
predict.fastVAR.VARX = function(VARX, xnew, n.ahead=1, threshold, ...) {
  y.pred = matrix(nrow=n.ahead, ncol=ncol(VARX$var.z$y.orig))
  colnames(y.pred) = colnames(VARX$var.z$y.orig)
  for(i in 1:n.ahead) {
    Z.ahead.y = as.vector(t(VARX$var.z$y.orig[
      ((nrow(VARX$var.z$y.orig)):
      (nrow(VARX$var.z$y.orig)-VARX$var.z$p+1))
    ,]))
    if(VARX$var.z$b == 0) {
      Z.ahead.x = xnew[1,]
    } else {
      Z.ahead.x = as.vector(t(VARX$var.z$x.orig[
        ((nrow(VARX$var.z$x.orig)):
        (nrow(VARX$var.z$x.orig)-VARX$var.z$b+1))
      ,]))
    }
    Z.ahead = c(1, Z.ahead.y, Z.ahead.x)
    y.ahead = Z.ahead %*% coef(VARX)
    if(!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if(length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
      y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if(i == n.ahead) break
    VARX$var.z$y.orig = rbind(VARX$var.z$y.orig, y.ahead)
    VARX$var.z$x.orig = rbind(VARX$var.z$x.orig, xnew[1,])
  }
  return (y.pred)
}
