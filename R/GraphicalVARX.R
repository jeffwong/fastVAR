#' @export
GraphicalVARX = function(y, freq = rep(NA,ncol(y)), x, p=1, b=1, intercept=T, weights=NULL, rho=0.01, getdiag=T) {
  if (p < 1) stop("p must be a positive integer")
  if (missing(x)) {
    return (GVAR(y, freq, p, intercept, weights, rho, getdiag))
  }
  y.seasons = deseason(y, freq)
  var.z = VARX.Z(y.seasons$remaining, x, p, b, intercept)
  if (!is.null(weights) & !is.vector(weights)) {
      weights = switch(weights,
                       exponential = exponentialWeights(var.z$Z, var.z$y.p),
                       linear = linearWeights(var.z$Z, var.z$y.p))
  }
  model = graphicalLm(var.z$Z, var.z$y.p, weights)
  if (sum(is.na(coef(model))) > 0) {
      warning("Multivariate lm has invalid coefficients.
               Check the rank of the design matrix")
  }
  result = structure(list(
                          model = model,
                          var.z = var.z,
                          seasons = y.seasons
                         ), class="fastVAR.GraphicalVARX")
  if (getdiag) result$diag = VAR.diag(result)

  return (result) 
}

#' GraphicalVARX Coefficients
#'
#' @param VARX an object of class fastVAR.VARX
#' @param ...
#' @return The coefficients for the VARX model
#' @export
coef.fastVAR.GraphicalVARX = function(GraphicalVARX, ...) {
  coef(GraphicalVARX$model, ...)
}

#' GraphicalVARX Predict
#'
#' Predict n steps ahead from a fastVAR.GraphicalVARX object
#' @param GraphicalVARX an object of class fastVAR.GraphicalVARX returned from GraphicalVARX
#' @param xnew a matrix of future values for the exogenous inputs.  Should contain
#'   n.ahead rows
#' @param n.ahead number of steps to predict
#' @param threshold threshold prediction values to be greater than this value
#' @param ... extra parameters to pass into the coefficients method
#'   for objects of type fastVAR.GraphicalVARX
#' @examples
#'   data(Canada)
#'   x = matrix(rnorm(84*4), 84, 4)
#'   predict(GraphicalVARX(Canada, x = x, p = 3, b = 2, intercept = F), xnew = matrix(rnorm(2*4),2,4), n.ahead = 2)
#' @export
predict.fastVAR.GraphicalVARX = function(GraphicalVARX, xnew, n.ahead=1, threshold, ...) {  
  freq = GraphicalVARX$seasons$freq
  freq.indices = which(!is.na(GraphicalVARX$seasons$freq))
  
  if (missing(xnew)) {
    if (length(freq.indices) > 0) 
      return (GraphicalVARX$var.z$Z %*% coef(GraphicalVARX) +
              GraphicalVARX$seasons$seasonal[-(1:GraphicalVARX$var.z$p),])
    else
      return (GraphicalVARX$var.z$Z %*% coef(GraphicalVARX))
  }

  if (nrow(xnew) != n.ahead) stop("xnew should have n.ahead rows")
  y.pred = matrix(nrow=n.ahead, ncol=ncol(GraphicalVARX$var.z$y.orig))
  colnames(y.pred) = colnames(GraphicalVARX$var.z$y.orig)
  for (i in 1:n.ahead) {
    Z.ahead.y = as.vector(t(GraphicalVARX$var.z$y.orig[
      ((nrow(GraphicalVARX$var.z$y.orig)):
      (nrow(GraphicalVARX$var.z$y.orig)-GraphicalVARX$var.z$p+1))
    ,]))
    if (GraphicalVARX$var.z$b == 0) {
      Z.ahead.x = xnew[i,]
    } else {
      Z.ahead.x = as.vector(t(GraphicalVARX$var.z$x.orig[
        ((nrow(GraphicalVARX$var.z$x.orig)):
        (nrow(GraphicalVARX$var.z$x.orig)-GraphicalVARX$var.z$b+1))
      ,]))
    }
    if(GraphicalVARX$var.z$intercept) Z.ahead = c(1, Z.ahead.y, Z.ahead.x)
    else Z.ahead = c(Z.ahead.y, Z.ahead.x)
    y.ahead = Z.ahead %*% coef(GraphicalVARX)
    if (!missing(threshold)) {
      threshold.indices = which(y.ahead < threshold)
      if (length(threshold.indices) > 0)
        y.ahead[threshold.indices] = threshold
      y.ahead[threshold.indices] = threshold
    }
    y.pred[i,] = y.ahead
    if (i == n.ahead) break
    GraphicalVARX$var.z$y.orig = rbind(GraphicalVARX$var.z$y.orig, y.ahead)
    GraphicalVARX$var.z$x.orig = rbind(GraphicalVARX$var.z$x.orig, xnew[i,])
  }
  if (length(freq.indices) > 0) {
    lastSeason = lastPeriod(GraphicalVARX$seasons) #returns a list
    y.pred.seasonal = sapply(freq.indices, function(i) {
      season.start = periodIndex(freq[i], nrow(GraphicalVARX$var.z$y.orig + 1))
      season.end = season.start + n.ahead - 1
      rep(lastSeason[[i]], ceiling(n.ahead / freq[i]))[season.start : season.end]
    })
    y.pred[,freq.indices] = y.pred[,freq.indices] + y.pred.seasonal
    return (y.pred)
  }
  else return (y.pred)
}
