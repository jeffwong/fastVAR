#' Ridge Path
#' 
#' Gives back a function to compute the ridge regression for any
#' l2 penalty
#' @param y response vector or matrix
#' @param design matrix
#' @export
ridgePath = function(Y, X, weights = NULL) {
  if (!is.null(weights)) {X = scale.rows(X, sqrt(weights)); Y = scale.rows(Y, sqrt(weights))}
  x.svd = svd(X)
  duty = diag(x.svd$d) %*% t(x.svd$u) %*% Y
  return (function(l2penalty) {
    x.svd$v %*% diag(1 / (x.svd$d^2 + l2penalty)) %*% duty
  })
}

#' Ridge Coefficients
#'
#' The full ridge path.  Allows a user to calculate
#' the coefficients of a linear model for any l2 penalty
#' @param ridgePath an object of class fastVAR.RidgePath
#' @param lambda the desired l2penalty
#' @export
coef.fastVAR.RidgePath = function(model, l2penalty) {
  if (missing(l2penalty)) model$ridgePath(model$l2penalty)
  else model$ridgePath(l2penalty)
}
