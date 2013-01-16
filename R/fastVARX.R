#' Vector Autoregression with Exogenous Inputs
#'
#' Fit a VAR model by creating the lagged design matrix
#' and fitting a multivariate response matrix to it.  Note that
#' the VARX response matrix omits the first max(p,b) responses from the input
#' matrix Y.  This is because observations in Y are modeled by the
#' max(p,b) previous values, so the first max(p,b) observations cannot be modeled.
#' @param y A matrix of endogenous variables where each column represents an individual time series
#' @param freq only used if the time series are periodic.  freq is a vector of
#'   frequencies for each of the time series, as in 'ts(y, freq = ...)'.
#'   If the time series are not periodic, then this vector can be a vector of NA
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
#' @examples
#'   data(Canada)
#'   x = matrix(rnorm(84*4), 84, 4)
#'   VARX(Canada, x = x, p = 3, b = 2, intercept=F)
#' @export
VARX = function(y, freq = rep(NA,ncol(y)), x, p=1, b=1, intercept=T, weights=NULL, l2penalty=NULL, getdiag=T) {
  if (p < 1) stop("p must be a positive integer")
  if (missing(x)) {
    return (VAR(y, freq, p, intercept, weights, l2penalty, getdiag))
  }
  y.seasons = deseason(y, freq)
  var.z = VARX.Z(y.seasons$remaining, x, p, b, intercept)
  if (!is.null(weights) & !is.vector(weights)) {
      weights = switch(weights,
                       exponential = exponentialWeights(var.z$Z, var.z$y.p),
                       linear = linearWeights(var.z$Z, var.z$y.p))
  }
  if (is.null(l2penalty)) {
    #model = lm(var.z$y.p ~ -1 + var.z$Z, weights = weights)
    model = fastMlm(var.z$Z, var.z$y.p, weights)
    if (sum(is.na(coef(model))) > 0) {
      warning("Multivariate lm has invalid coefficients.
               Check the rank of the design matrix")
    }
    result = structure(list(
                            model = model,
                            var.z = var.z,
                            seasons = y.seasons
                           ), class="fastVAR.VARX")
  } else {
    #Compute full path ridge solution
    result = structure(list(
                            model = structure(list(ridgePath = ridgePath(var.z$y.p, var.z$Z, weights),
                                                   l2penalty = l2penalty),
                                              class="fastVAR.RidgePath"),
                            var.z = var.z,
                            seasons = y.seasons),
                       class="fastVAR.VARX")
  }
  if (getdiag) result$diag = VAR.diag(result)

  return (result) 
}

#' VARX Coefficients
#'
#' If the VARX object was fit using a l2 penalty, then the full ridge path was
#' calculated and stored in the object.  This means the user can adjust the ridge penalty
#' term here and recompute the coefficients of the VARX
#' @param VARX an object of class fastVAR.VARX
#' @param ... if VAR was fit using a l2 penalty, the user can specify a different
#'   l2 penalty here and have the coefficients recomputed
#' @return The coefficients for the VARX model
#' @method coef fastVAR.VARX
#' @S3method coef fastVAR.VARX
coef.fastVAR.VARX = function(VARX, ...) {
  coef(VARX$model, ...)
}

#' VARX Predict
#'
#' Predict n steps ahead from a fastVAR.VARX object
#' @param VARX an object of class fastVAR.VARX returned from VARX
#' @param xnew a matrix of future values for the exogenous inputs.  Should contain
#'   n.ahead rows
#' @param n.ahead number of steps to predict
#' @param threshold threshold prediction values to be greater than this value
#' @param ... extra parameters to pass into the coefficients method
#'   for objects of type fastVAR.VARX
#' @examples
#'   data(Canada)
#'   x = matrix(rnorm(84*4), 84, 4)
#'   predict(VARX(Canada, x = x, p = 3, b = 2, intercept = F), xnew = matrix(rnorm(2*4),2,4), n.ahead = 2)
#' @method predict fastVAR.VARX
#' @S3method predict fastVAR.VARX
predict.fastVAR.VARX = function(VARX, xnew, n.ahead=1, threshold, ...) {
  freq = VARX$seasons$freq
  freq.indices = which(!is.na(VARX$seasons$freq))
  if (missing(xnew)) {
    if (length(freq.indices) > 0)
      return (VARX$var.z$Z %*% coef(VARX) +
              VARX$seasons$seasonal[-(1:VARX$var.z$p),])
    else
      return (VARX$var.z$Z %*% coef(VARX))
  }
  if (nrow(xnew) != n.ahead) stop("xnew should have n.ahead rows")
  y.pred = matrix(nrow=n.ahead, ncol=ncol(VARX$var.z$y.orig))
  colnames(y.pred) = colnames(VARX$var.z$y.orig)
  y.orig = VARX$var.z$y.orig
  for (i in 1:n.ahead) {
    Z.ahead.y = as.vector(t(y.orig[
      ((nrow(y.orig)):
      (nrow(y.orig)-VARX$var.z$p+1))
    ,]))
    if (VARX$var.z$b == 0) {
      Z.ahead.x = xnew[i,]
    } else {
      Z.ahead.x = as.vector(t(VARX$var.z$x.orig[
        ((nrow(VARX$var.z$x.orig)):
        (nrow(VARX$var.z$x.orig)-VARX$var.z$b+1))
      ,]))
    }
    if(VARX$var.z$intercept) Z.ahead = c(1, Z.ahead.y, Z.ahead.x)
    else Z.ahead = c(Z.ahead.y, Z.ahead.x)
    y.ahead = Z.ahead %*% coef(VARX)
    if (!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if (length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
      y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if (i == n.ahead) break
    y.orig = rbind(y.orig, y.ahead)
    VARX$var.z$x.orig = rbind(VARX$var.z$x.orig, xnew[i,])
  }
  if (length(freq.indices) > 0) {
    lastSeason = lastPeriod(VARX$seasons) #returns a list
    y.pred.seasonal = sapply(freq.indices, function(i) {
      season.start = periodIndex(freq[i], nrow(VARX$var.z$y.orig) + 1)
      season.end = season.start + n.ahead - 1
      rep(lastSeason[[i]], ceiling(n.ahead / freq[i]))[season.start : season.end]
    })
    y.pred[,freq.indices] = y.pred[,freq.indices] + y.pred.seasonal
    return (y.pred)
  }
  else return (y.pred)
}
