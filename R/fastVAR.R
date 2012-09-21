#' Vector Autoregression
#'
#' Fit a VAR model by creating the lagged design matrix
#' and fitting a multivariate response matrix to it
#' @param y A matrix where each column represents an individual time series
#' @param p the number of lags to include in the design matrix
#' @param getdiag logical.  If true, return diagnostics
#' @param weights weights applied to the multiresponse linear regression.
#'   Better predictions might come from weighting observations far in the past
#'   less so they impact the objective value less.
#' @param l2penalty a ridge regression penalty, useful when the design matrix is 
#'   very wide, which may happen if p is large.
#' @export
VAR = function(y, p=1, weights=NULL, l2penalty=NULL, getdiag=T) {
  if(p < 1) {
    stop("p must be a positive integer")
  }
  var.z = VAR.Z(y, p)
  model = lm(var.z$y.p ~ var.z$Z)
  if(any(is.na(model$coefficients))) {
    stop("Multivariate lm has invalid coefficients.  
          Check the rank of the design matrix")
  }
  if(getdiag) {
    return ( structure ( list (
      model = model,
      var.z = var.z,
      diag = VAR.diag(var.z$y.p, var.z$Z, model$coefficients,
        var.z$n, var.z$T, (var.z$k+1), p, var.z$dof)
    ), class="fastVAR.VAR"))
  } else {
    return (structure( list (
                model = model,
                var.z = var.z
    ), class="fastVAR.VAR"))
  }
}

#' Exponential Weights
#' 
#' Weights that decay exponentially.  Values in the past receive smaller weights
#' @param x
#' @param y
#' @export
exponentialWeights = function(x, y) {

}

predict.fastVAR.VAR = function(VAR, n.ahead=1, threshold) {
  y.pred = matrix(nrow=n.ahead, ncol=ncol(VAR$var.z$y.orig))
  colnames(y.pred) = colnames(VAR$var.z$y.orig)
  for(i in 1:n.ahead) {
    Z.ahead = c(1,as.vector(t(VAR$var.z$y.orig[
      ((nrow(VAR$var.z$y.orig)):
      (nrow(VAR$var.z$y.orig)-VAR$var.z$p+1))
    ,])))
    y.ahead = Z.ahead %*% VAR$model$coef
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

#' Inverse Covariance Matrix
#'
#' Computes the inverse of the covariance matrix
#' using an svd
#' @param Z the design matrix
.inverseCovariance = function(Z) {
    Z.centered = scale(Z)
    Z.centered.svd = svd(Z.centered.svd)
    d2.inv.1 = 1/(Z.centered.svd$d[which(Z.centered.svd$d != 0)])^2
    d2.inv.2 = rep(0, length(which(Z.centered.svd$d == 0)))
    if(length(d2.inv.2) > 0) warning("Warning, covariance matrix not invertible.  Setting the inverse of the singular values to 0")
    Z.centered.svd$v %*% diag(c(d2.inv.1, d2.inv.2)) %*% t(Z.centered.svd$v)
}

VAR.diag = function(y.p, Z, B, n, T, k, p, dof) {
  Z = cbind(1, Z)
  colnames(Z)[1] = 'intercept'
  ZTZ.inv = .inverseCovariance(Z)
  residual = y.p - Z %*% B
  residual.cov = t(residual) %*% residual  #residual covariance matrix
  sigma.hat = (1 / (T - k)) * 
    (residual.cov)
  cov.hat = kronecker(sigma.hat, ZTZ.inv)
  se = matrix(sqrt(diag(cov.hat)), nrow=nrow(B), ncol=ncol(B))
  tvalue = B / se
  pvalue = 2*(1 - pt(abs(tvalue),dof)) 
  rownames(se) = colnames(Z)
  rownames(tvalue) = colnames(Z)
  rownames(pvalue) = colnames(Z)
  colnames(se) = colnames(y.p)
  colnames(tvalue) = colnames(y.p)
  colnames(pvalue) = colnames(y.p)

  #Get model selection statistics
  sigma.tilde = 1 / T * residual.cov
  sigma.tilde.det = log(det(sigma.tilde))
  
  AIC = sigma.tilde.det + 2/T * p * n^2
  BIC = sigma.tilde.det + log(T)/T * p * n^2
  HQ = sigma.tilde.det + 2*log(log(T))/T * p * n^2

  return (list(
    se=se,
    tvalue=tvalue,
    pvalue=pvalue,
    residual=residual,
    cov=cov.hat,
    AIC=AIC,
    BIC=BIC,
    HQ=HQ
  ))
}
