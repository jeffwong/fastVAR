#' Exponential Weights
#' 
#' Weights that decay exponentially.  Values in the past receive smaller weights
#' @param x
#' @param y
#' @export
exponentialWeights = function(x, y) {
  if(is.vector(y)) l = length(y)
  else if(is.data.frame(y) | is.matrix(y)) l = nrow(y)
  cumsum(1:l) / sum(1:l)
}

#' Linear Weights
#' 
#' Weights that decay linearly.  Values in the past receive smaller weights
#' @param x
#' @param y
#' @export
linearWeights = function(x, y) {
  if(is.vector(y)) l = length(y)
  else if(is.data.frame(y) | is.matrix(y)) l = nrow(y)
  1:l
}
