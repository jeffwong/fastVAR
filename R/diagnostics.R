#' Inverse Covariance Matrix
#'
#' Computes the inverse of the covariance matrix
#' using an svd
#' @param Z the design matrix
.inverseCovariance = function(Z) {
    Z.centered = scale(Z, scale=F)
    Z.centered.svd = svd(Z.centered)
    d2.inv.1 = 1/(Z.centered.svd$d[which(Z.centered.svd$d != 0)])^2
    d2.inv.2 = rep(0, length(which(Z.centered.svd$d == 0)))
    if(length(d2.inv.2) > 0) warning("Warning, covariance matrix not invertible.  Setting the inverse of the singular values to 0")
    Z.centered.svd$v %*% diag(c(d2.inv.1, d2.inv.2)) %*% t(Z.centered.svd$v)
}

#' VAR Diagnostics
#'
#' Compute diagnostics for VAR models returned from VAR
#' @param VAR a model of class fastVAR.VAR returned by function VAR
#' @export
VAR.diag = function(VAR) {
  ZTZ.inv = .inverseCovariance(VAR$var.z$Z)
  residual = VAR$var.z$y.p - VAR$var.z$Z %*% coef(VAR)
  residual.cov = t(residual) %*% residual  #residual covariance matrix
  sigma.hat = (1 / (VAR$var.z$T - VAR$var.z$k)) * 
    (residual.cov)
  cov.hat = kronecker(sigma.hat, ZTZ.inv)
  se = matrix(sqrt(diag(cov.hat)), nrow=nrow(coef(VAR)), ncol=ncol(coef(VAR)))
  tvalue = coef(VAR) / se
  pvalue = 2*(1 - pt(abs(tvalue),VAR$var.z$dof)) 
  rownames(se) = colnames(VAR$var.z$Z)
  rownames(tvalue) = colnames(VAR$var.z$Z)
  rownames(pvalue) = colnames(VAR$var.z$Z)
  colnames(se) = colnames(VAR$var.z$y.p)
  colnames(tvalue) = colnames(VAR$var.z$y.p)
  colnames(pvalue) = colnames(VAR$var.z$y.p)

  #Get model selection statistics
  sigma.tilde = 1 / VAR$var.z$T * residual.cov
  sigma.tilde.det = log(det(sigma.tilde))
  
  n = ncol(VAR$var.z$y.p)

  AIC = sigma.tilde.det + 2/VAR$var.z$T * VAR$var.z$p * n^2
  BIC = sigma.tilde.det + log(VAR$var.z$T)/VAR$var.z$T * VAR$var.z$p * n^2
  HQ = sigma.tilde.det + 2*log(log(VAR$var.z$T))/VAR$var.z$T * VAR$var.z$p * n^2

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
