GraphicalVAR = function(y, freq = rep(NA,ncol(y)), p=1, intercept = T, weights=NULL, rho=0.01, getdiag=T) {
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
  model = graphicalLm(var.z$Z, var.z$y.p, weights, rho)
  if (any(is.na(coef(model)))) {
      warning("Multivariate lm has invalid coefficients.  
               Check the rank of the design matrix")
  }
  result = structure(list(
                          model = model,
                          var.z = var.z,
                          seasons = y.seasons
                         ),
                     class="fastVAR.GVAR")

  if (getdiag) result$diag = VAR.diag(result)

  return (result)
}

coef.fastVAR.GVAR = function(GVAR, ...) {
  coef(GVAR$model, ...)
}

predict.fastVAR.GVAR = function(GVAR, n.ahead=1, threshold, ...) {
  y.pred = matrix(nrow=n.ahead, ncol=ncol(GVAR$var.z$y.orig))
  colnames(y.pred) = colnames(GVAR$var.z$y.orig)
  for (i in 1:n.ahead) {
    if (GVAR$var.z$intercept) {
      Z.ahead = c(1,as.vector(t(GVAR$var.z$y.orig[
        ((nrow(GVAR$var.z$y.orig)):
        (nrow(GVAR$var.z$y.orig)-GVAR$var.z$p+1))
      ,])))
    } else {
      Z.ahead = as.vector(t(GVAR$var.z$y.orig[
        ((nrow(GVAR$var.z$y.orig)):
        (nrow(GVAR$var.z$y.orig)-GVAR$var.z$p+1))
      ,]))
    }
    y.ahead = Z.ahead %*% coef(GVAR, ...)
    if (!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if (length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if (i == n.ahead) break
    GVAR$var.z$y.orig = rbind(GVAR$var.z$y.orig, y.ahead)
  }
  freq = GVAR$seasons$freq
  freq.indices = which(!is.na(GVAR$seasons$freq))
  if (length(freq.indices) > 0) {
    lastSeason = lastPeriod(GVAR$seasons) #returns a list
    y.pred.seasonal = sapply(freq.indices, function(i) {
      season.start = periodIndex(freq[i], nrow(GVAR$var.z$y.orig + 1))
      season.end = season.start + n.ahead - 1
      rep(lastSeason[[i]], ceiling(n.ahead / freq[i]))[season.start : season.end]
    })
    y.pred[,freq.indices] = y.pred[,freq.indices] + y.pred.seasonal
    return (y.pred)
  }
  else return (y.pred)
}
