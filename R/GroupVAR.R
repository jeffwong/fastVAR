#' Group Vector Autoregression via Group Lasso
#'
#' Fit a VAR model by creating the lagged design matrix
#' and fitting a multivariate response matrix to it.  Note that
#' the VAR response matrix omits the first p responses from the input
#' matrix Y.  This is because observations in Y are modeled by the
#' p previous values, so the first p observations cannot be modeled.
#'
#' While multivariate response regressions can be solved as multiple
#' univariate response regressions, this multivariate response problem
#' can better be solved by using Group Lasso.  Instead of seeking sparsity
#' in the coefficients for each univariate response, Group Lasso attempts to
#' find sparsity in the coefficient matrix as a whole.
#' @param y A matrix where each column represents an individual time series
#' @param freq only used if the time series are periodic.  freq is a vector of
#'   frequencies for each of the time series, as in 'ts(y, freq = ...)'.
#'   If the time series are not periodic, then this vector can be a vector of NA
#' @param p the number of lags to include in the design matrix
#' @param weights weights applied to the multiresponse linear regression.
#'   Better predictions might come from weighting observations far in the past
#'   less so they impact the objective value less.  Either NULL,
#'   "exponential", or "linear"
#' @param alpha the elastic net mixing parameter, as defined in 'glmnet'
#' @param getdiag logical.  If true, return diagnostics
#' @examples
#'   data(Canada)
#'   GroupVAR(Canada, p = 3)
#' @export
GroupVAR = function(y, freq = rep(NA,ncol(y)), p=1, weights=NULL, alpha = 0.4, getdiag=T) {
  if (p < 1) {
    stop("p must be a positive integer")
  }
  y.seasons = deseason(y, freq)
  var.z = VAR.Z(y.seasons$remaining, p, intercept = T)
  if (is.null(weights)) weights = 1:nrow(var.z$y.p)
  if (is.character(weights)) {
      weights = switch(weights,
                       exponential = exponentialWeights(var.z$Z, var.z$y.p),
                       linear = linearWeights(var.z$Z, var.z$y.p)
                      )
  }
  model = cv.glmnet(var.z$Z[,-1], var.z$y.p, family = "mgaussian", weights = weights, alpha = alpha)
  result = structure(list(model = model,
                          var.z = var.z,
                          seasons = y.seasons),
                     class = "fastVAR.GroupVAR")

  if (getdiag) result$diag = VAR.diag(result)

  return (result)
}

#' GroupVAR Coefficients
#'
#' @param GroupVAR an object of class fastVAR.GroupVAR
#' @param ... extra parameters to pass into the coefficients method
#'   for objects of type fastVAR.GroupVAR
#' @return The coefficients for the VAR model
#' @method coef fastVAR.GroupVAR
#' @S3method coef fastVAR.GroupVAR
coef.fastVAR.GroupVAR = function(GroupVAR, ...) {
  coefs = do.call('cbind', lapply(coef(GroupVAR$model, ...), as.matrix))
  colnames(coefs) = colnames(GroupVAR$var.z$y.orig)
  return (coefs)
}

#' GroupVAR Predict
#'
#' Predict n steps ahead from a fastVAR.GroupVAR object
#' @param GroupVAR an object of class fastVAR.GroupVAR returned from GroupVAR
#' @param n.ahead number of steps to predict
#' @param threshold threshold prediction values to be greater than this value
#' @param ... extra parameters to pass into the coefficients method
#'   for objects of type fastVAR.GroupVAR
#' @examples
#'   data(Canada)
#'   predict(GroupVAR(Canada, p = 3), 1)
#' @method predict fastVAR.GroupVAR
#' @S3method predict fastVAR.GroupVAR
predict.fastVAR.GroupVAR = function(GroupVAR, n.ahead, threshold, ...) {
  freq = GroupVAR$seasons$freq
  freq.indices = which(!is.na(GroupVAR$seasons$freq))
  if (missing(n.ahead)) {
    if (length(freq.indices) > 0) 
      return (GroupVAR$var.z$Z %*% coef(GroupVAR) +
              GroupVAR$seasons$seasonal[-(1:GroupVAR$var.z$p),])
    else
      return (GroupVAR$var.z$Z %*% coef(GroupVAR))
  }
  y.pred = matrix(nrow=n.ahead, ncol=ncol(GroupVAR$var.z$y.orig))
  colnames(y.pred) = colnames(GroupVAR$var.z$y.orig)
  for (i in 1:n.ahead) {
    if (GroupVAR$var.z$intercept) {
      Z.ahead = c(1,as.vector(t(GroupVAR$var.z$y.orig[
        ((nrow(GroupVAR$var.z$y.orig)):
        (nrow(GroupVAR$var.z$y.orig)-GroupVAR$var.z$p+1))
      ,])))
    } else {
      Z.ahead = as.vector(t(GroupVAR$var.z$y.orig[
        ((nrow(GroupVAR$var.z$y.orig)):
        (nrow(GroupVAR$var.z$y.orig)-GroupVAR$var.z$p+1))
      ,]))
    }
    y.ahead = Z.ahead %*% coef(GroupVAR, ...)
    if (!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if (length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if (i == n.ahead) break
    GroupVAR$var.z$y.orig = rbind(GroupVAR$var.z$y.orig, y.ahead)
  }
  if (length(freq.indices) > 0) {
    lastSeason = lastPeriod(GroupVAR$seasons) #returns a list
    y.pred.seasonal = sapply(freq.indices, function(i) {
      season.start = periodIndex(freq[i], nrow(GroupVAR$var.z$y.orig + 1))
      season.end = season.start + n.ahead - 1
      rep(lastSeason[[i]], ceiling(n.ahead / freq[i]))[season.start : season.end]
    })
    y.pred[,freq.indices] = y.pred[,freq.indices] + y.pred.seasonal
    return (y.pred)
  }
  else return (y.pred)
}
