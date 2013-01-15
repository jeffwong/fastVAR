#' Group Vector Autoregression with Exogenous Inputs via Group Lasso
#'
#' Fit a VAR model by creating the lagged design matrix
#' and fitting a multivariate response matrix to it.  Note that
#' the GroupVARX response matrix omits the first max(p,b) responses from the input
#' matrix Y.  This is because observations in Y are modeled by the
#' max(p,b) previous values, so the first max(p,b) observations cannot be modeled.

#' While multivariate response regressions can be solved as multiple
#' univariate response regressions, this multivariate response problem
#' can better be solved by using Group Lasso.  Instead of seeking sparsity
#' in the coefficients for each univariate response, Group Lasso attempts to
#' find sparsity in the coefficient matrix as a whole.
#' @param y A matrix of endogenous variables where each column represents an individual time series
#' @param freq only used if the time series are periodic.  freq is a vector of
#'   frequencies for each of the time series, as in 'ts(y, freq = ...)'.
#'   If the time series are not periodic, then this vector can be a vector of NA
#' @param x A matrix of exogenous variables where each column represents an individual time series
#' @param p the number of lags of Y to include in the design matrix
#' @param b the number of lags to X include in the design matrix
#' @param weights weights applied to the multiresponse linear regression.
#'   Better predictions might come from weighting observations far in the past
#'   less so they impact the objective value less.  Either NULL,
#'   "exponential", or "linear"
#' @param alpha the mixing parameter between group lasso and quadratic, as in 'glmnet'
#' @param getdiag logical.  If true, return diagnostics
#' @examples
#'   data(Canada)
#'   x = matrix(rnorm(84*4), 84, 4)
#'   GroupVARX(y, x=x, p=3, b=2)
#' @export
GroupVARX = function(y, freq = rep(NA,ncol(y)), x, p=1, b=1, weights=NULL, alpha = 0.4, getdiag=T) {
  if (p < 1) stop("p must be a positive integer")
  if (missing(x)) {
    return (GroupVAR(y, freq, p, weights, getdiag))
  }
  y.seasons = deseason(y, freq)
  var.z = VARX.Z(y.seasons$remaining, x, p, b, intercept = T)
  if (is.null(weights)) weights = 1:nrow(var.z$y.p)
  if (is.character(weights)) {
      weights = switch(weights,
                       exponential = exponentialWeights(var.z$Z, var.z$y.p),
                       linear = linearWeights(var.z$Z, var.z$y.p))
  }
  model = cv.glmnet(var.z$Z[,-1], var.z$y.p, family = "mgaussian", weights = weights, alpha = alpha)
  result = structure(list(model = model,
                          var.z = var.z,
                          seasons = y.seasons),
                     class = "fastVAR.GroupVARX")

  if (getdiag) result$diag = VAR.diag(result)

  return (result) 
}

#' GroupVARX Coefficients
#'
#' @param GroupVARX an object of class fastVAR.GroupVARX
#' @param ... extra parameters to pass into the coefficients method
#'   for objects of type fastVAR.GroupVARX
#' @return The coefficients for the GroupVARX model
#' @method coef fastVAR.GroupVARX
#' @S3method coef fastVAR.GroupVARX
coef.fastVAR.GroupVARX = function(GroupVARX, ...) {
  do.call('cbind', lapply(coef(GroupVAR$model, ...), as.matrix))
  colnames(coefs) = colnames(GroupVAR$var.z$y.orig)
  return (coefs)
}

#' GroupVARX Predict
#'
#' Predict n steps ahead from a fastVAR.GroupVARX object
#' @param GroupVARX an object of class fastVAR.GroupVARX returned from GroupVARX
#' @param xnew a matrix of future values for the exogenous inputs.  Should contain
#'   n.ahead rows
#' @param n.ahead number of steps to predict
#' @param threshold threshold prediction values to be greater than this value
#' @param ... extra parameters to pass into the coefficients method
#'   for objects of type fastVAR.GroupVARX
#' @examples
#'   data(Canada)
#'   x = matrix(rnorm(84*4), 84, 4)
#'   predict(GroupVARX(Canada, x = x, p = 3, b = 2), xnew = matrix(rnorm(2*4),2,4), n.ahead = 2)
#' @method predict fastVAR.GroupVARX
#' @S3method predict fastVAR.GroupVARX
predict.fastVAR.GroupVARX = function(GroupVARX, xnew, n.ahead=1, threshold, ...) {
  freq = GroupVARX$seasons$freq
  freq.indices = which(!is.na(GroupVARX$seasons$freq))

  if (missing(xnew)) {
    if (length(freq.indices) > 0)
      return (GroupVARX$var.z$Z %*% coef(GroupVARX) +
              GroupVARX$seasons$seasonal[-(1:GroupVARX$var.z$p),])
    else
      return (GroupVARX$var.z$Z %*% coef(GroupVARX))
  }
  if (nrow(xnew) != n.ahead) stop("xnew should have n.ahead rows")
  y.pred = matrix(nrow=n.ahead, ncol=ncol(GroupVARX$var.z$y.orig))
  colnames(y.pred) = colnames(GroupVARX$var.z$y.orig)
  y.orig = GroupVARX$var.z$y.orig
  for (i in 1:n.ahead) {
    Z.ahead.y = as.vector(t(y.orig[
      ((nrow(y.orig)):
      (nrow(y.orig)-GroupVARX$var.z$p+1))
    ,]))
    if (GroupVARX$var.z$b == 0) {
      Z.ahead.x = xnew[i,]
    } else {
      Z.ahead.x = as.vector(t(GroupVARX$var.z$x.orig[
        ((nrow(GroupVARX$var.z$x.orig)):
        (nrow(GroupVARX$var.z$x.orig)-GroupVARX$var.z$b+1))
      ,]))
    }
    if(GroupVARX$var.z$intercept) Z.ahead = c(1, Z.ahead.y, Z.ahead.x)
    else Z.ahead = c(Z.ahead.y, Z.ahead.x)
    y.ahead = Z.ahead %*% coef(GroupVARX)
    if (!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if (length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
      y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if (i == n.ahead) break
    y.orig = rbind(y.orig, y.ahead)
    GroupVARX$var.z$x.orig = rbind(GroupVARX$var.z$x.orig, xnew[i,])
  }
  if (length(freq.indices) > 0) {
    lastSeason = lastPeriod(GroupVARX$seasons) #returns a list
    y.pred.seasonal = sapply(freq.indices, function(i) {
      season.start = periodIndex(freq[i], nrow(GroupVARX$var.z$y.orig) + 1)
      season.end = season.start + n.ahead - 1
      rep(lastSeason[[i]], ceiling(n.ahead / freq[i]))[season.start : season.end]
    })
    y.pred[,freq.indices] = y.pred[,freq.indices] + y.pred.seasonal
    return (y.pred)
  }
  else return (y.pred)
}
