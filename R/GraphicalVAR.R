#' @export
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
                     class="fastVAR.GraphicalVAR")

  if (getdiag) result$diag = VAR.diag(result)

  return (result)
}

#' @method coef fastVAR.GraphicalVAR
#' @S3method coef fastVAR.GraphicalVAR
coef.fastVAR.GraphicalVAR = function(GraphicalVAR, ...) {
  coef(GraphicalVAR$model, ...)
}

#' @method predict fastVAR.GraphicalVAR
#' @S3method predict fastVAR.GraphicalVAR
predict.fastVAR.GraphicalVAR = function(GraphicalVAR, n.ahead, threshold, ...) {
  freq = GraphicalVAR$seasons$freq
  freq.indices = which(!is.na(GraphicalVAR$seasons$freq))

  if (missing(n.ahead)) {
    if (length(freq.indices) > 0) 
      return (GraphicalVAR$var.z$Z %*% coef(GraphicalVAR) +
              GraphicalVAR$seasons$seasonal[-(1:GraphicalVAR$var.z$p),])
    else
      return (GraphicalVAR$var.z$Z %*% coef(GraphicalVAR))
  }

  y.pred = matrix(nrow=n.ahead, ncol=ncol(GraphicalVAR$var.z$y.orig))
  colnames(y.pred) = colnames(GraphicalVAR$var.z$y.orig)
  y.orig = GraphicalVAR$var.z$y.orig
  for (i in 1:n.ahead) {
    if (GraphicalVAR$var.z$intercept) {
      Z.ahead = c(1,as.vector(t(y.orig[
        ((nrow(y.orig)):
        (nrow(y.orig)-GraphicalVAR$var.z$p+1))
      ,])))
    } else {
      Z.ahead = as.vector(t(y.orig[
        ((nrow(y.orig)):
        (nrow(y.orig)-GraphicalVAR$var.z$p+1))
      ,]))
    }
    y.ahead = Z.ahead %*% coef(GraphicalVAR, ...)
    if (!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if (length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if (i == n.ahead) break
    y.orig = rbind(y.orig, y.ahead)
  }
  if (length(freq.indices) > 0) {
    lastSeason = lastPeriod(GraphicalVAR$seasons) #returns a list
    y.pred.seasonal = sapply(freq.indices, function(i) {
      season.start = periodIndex(freq[i], nrow(GraphicalVAR$var.z$y.orig) + 1)
      season.end = season.start + n.ahead - 1
      rep(lastSeason[[i]], ceiling(n.ahead / freq[i]))[season.start : season.end]
    })
    y.pred[,freq.indices] = y.pred[,freq.indices] + y.pred.seasonal
    return (y.pred)
  }
  else return (y.pred)
}
