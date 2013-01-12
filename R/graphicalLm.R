graphicalLm = function(X, Y, weights, rho) {
  if (is.null(weights)) weights = rep(1, nrow(Y))
  X = scale(scale.rows(X, sqrt(weights)), center = T, scale = F)
  Y = scale.rows(Y, sqrt(weights))
  X.glasso = glasso(var(X) * (nrow(X) - 1), rho)
  coefficients = X.glasso$wi %*% crossprod(X, Y)
  return (structure(list(coefficients = coefficients),
                    class = "fastVAR.graphicalLm"))
}

coef.fastVAR.graphicalLm = function(graphicalLm, ...) {
  return (graphicalLm$coefficients)
}
