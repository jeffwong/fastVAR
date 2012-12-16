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
#' @examples
#'   data(Canada)
#'   VAR(Canada, p = 3, intercept = F)
#' @export
VAR = function(y, freq = rep(NA,ncol(y)), p=1, intercept = T, weights=NULL, l2penalty=NULL, getdiag=T) {
  if (p < 1) {
    stop("p must be a positive integer")
  }
  y.seasons = deseason(y, freq)
  var.z = VAR.Z(y.seasons$remaining, p, intercept)
  if (!is.null(weights) & !is.vector(weights)) {
      weights = switch(weights,
                       exponential = exponentialWeights(var.z$Z, var.z$y.p),
                       linear = linearWeights(var.z$Z, var.z$y.p)
                      )
  }
  if (is.null(l2penalty)) {
    #model = lm(var.z$y.p ~ -1 + var.z$Z, weights = weights)
    model = fastMlm(var.z$Z, var.z$y.p, weights)
    if (any(is.na(coef(model)))) {
      warning("Multivariate lm has invalid coefficients.  
               Check the rank of the design matrix")
    }
    result = structure(list(
                            model = model,
                            var.z = var.z,
                            seasons = y.seasons
                            ),
                       class="fastVAR.VAR")
  } else {
    #Compute full path ridge solution
    ridge.coef = ridgePath(var.z$y.p, var.z$Z, weights)
    result = structure(list(
                            model = structure(list(ridgePath=ridge.coef, l2penalty=l2penalty), class="fastVAR.RidgePath"),
                            var.z = var.z,
                            seasons = y.seasons),
                       class="fastVAR.VAR")
  }

  if (getdiag) result$diag = VAR.diag(result)

  return (result)
}

#' VAR Coefficients
#'
#' If the VAR object was fit using a l2 penalty, then the full ridge path was
#' calculated and stored in the object.  This means the user can adjust the ridge penalty
#' term here and recompute the coefficients of the VAR
#' @param VAR an object of class fastVAR.VAR
#' @param ... if VAR was fit using a l2 penalty, the user can specify a different
#'   l2 penalty here and have the coefficients recomputed
#' @return The coefficients for the VAR model
#' @export
coef.fastVAR.VAR = function(VAR, ...) {
  coef(VAR$model, ...)
}

#' VAR Predict
#'
#' Predict n steps ahead from a fastVAR.VAR object
#' @param VAR an object of class fastVAR.VAR returned from VAR
#' @param n.ahead number of steps to predict
#' @param threshold threshold prediction values to be greater than this value
#' @param ... extra parameters to pass into the coefficients method
#'   for objects of type fastVAR.VAR
#' @examples
#'   data(Canada)
#'   predict(VAR(Canada, p = 3, intercept = F), 1)
#' @export
predict.fastVAR.VAR = function(VAR, n.ahead=1, threshold, ...) {
  y.pred = matrix(nrow=n.ahead, ncol=ncol(VAR$var.z$y.orig))
  colnames(y.pred) = colnames(VAR$var.z$y.orig)
  for (i in 1:n.ahead) {
    if (VAR$var.z$intercept) {
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
    if (!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if (length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if (i == n.ahead) break
    VAR$var.z$y.orig = rbind(VAR$var.z$y.orig, y.ahead)
  }
  freq = VAR$seasons$freq
  freq.indices = which(!is.na(VAR$seasons$freq))
  if (length(freq.indices) > 0) {
    lastSeason = lastPeriod(VAR$seasons) #returns a list
    y.pred.seasonal = sapply(freq.indices, function(i) {
      season.start = periodIndex(freq[i], nrow(VAR$var.z$y.orig + 1))
      season.end = season.start + n.ahead - 1
      rep(lastSeason[[i]], ceiling(n.ahead / freq[i]))[season.start : season.end]
    })
    y.pred[,freq.indices] = y.pred[,freq.indices] + y.pred.seasonal
    return (y.pred)
  }
  else return (y.pred)
}
