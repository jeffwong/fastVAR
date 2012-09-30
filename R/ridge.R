#' Ridge Path
#'
#' The full ridge path.  Allows a user to calculate
#' the coefficients of a linear model for any l2 penalty
#' @param ridgePath an object of class fastVAR.RidgePath
#' @param lambda the desired l2penalty
#' @export
coef.fastVAR.RidgePath = function(ridgePath, lambda) {
  ridgePath(lambda)
}
