VARX = function(y, x, p, b, getdiag=T) {
  if(p < 1) stop("p must be a positive integer")
  if(missing(x)) {
    return (VAR(y, p, getdiag))
  }
  varx.z = VARX.Z(y, x, p, b)
  model = lm(varx.z$y.p ~ varx.z$Z)
  if(sum(is.na(model$coefficients)) > 0) {
    stop("Multivariate lm has invalid coefficients.
          Check the rank of the design matrix")
  }
  if(getdiag) {
    return ( structure ( list (
      model = model,
      varx.z = varx.z,
      diag = VARdiag(varx.z$y.p, varx.z$Z, model$coefficients,
        varx.z$ny, varx.z$T, (varx.z$k+1), varx.z$p.max, varx.z$dof)
    ), class="fastVAR.VARX"))
  } else {
    return (structure( list( 
                model = model,
                varx.z = varx.z
    ), class="fastVAR.VARX"))
  }
}

predict.fastVAR.VARX = function(VARX, xnew, n.ahead=1, threshold) {
  y.pred = matrix(nrow=n.ahead, ncol=ncol(VARX$varx.z$y.orig))
  colnames(y.pred) = colnames(VARX$varx.z$y.orig)
  for(i in 1:n.ahead) {
    Z.ahead.y = as.vector(t(VARX$varx.z$y.orig[
      ((nrow(VARX$varx.z$y.orig)):
      (nrow(VARX$varx.z$y.orig)-VARX$varx.z$p+1))
    ,]))
    if(VARX$varx.z$b == 0) {
      Z.ahead.x = xnew[1,]
    } else {
      Z.ahead.x = as.vector(t(VARX$varx.z$x.orig[
        ((nrow(VARX$varx.z$x.orig)):
        (nrow(VARX$varx.z$x.orig)-VARX$varx.z$b+1))
      ,]))
    }
    Z.ahead = c(1, Z.ahead.y, Z.ahead.x)
    y.ahead = Z.ahead %*% VARX$model$coef
    if(!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if(length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
      y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if(i == n.ahead) break
    VARX$varx.z$y.orig = rbind(VARX$varx.z$y.orig, y.ahead)
    VARX$varx.z$x.orig = rbind(VARX$varx.z$x.orig, xnew[1,])
  }
  return (y.pred)
}
