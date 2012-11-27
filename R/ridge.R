#' Ridge Path
#' 
#' Gives back a function to compute the ridge regression for any
#' l2 penalty
#' @param y response vector or matrix
#' @param design matrix
#' @export
ridgePath = function(y, x) {
  x.svd = svd(x)
  duty = diag(x.svd$d) %*% t(x.svd$u) %*% y
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
coef.fastVAR.RidgePath = function(ridgePath, lambda) {
  ridgePath(lambda)
}
