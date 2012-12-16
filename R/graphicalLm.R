graphicalLm = function(X, Y, weights, rho) {
  X = scale.rows(X, sqrt(weights))
  Y = scale.rows(Y, sqrt(weights))
  cov.inv = glasso(var(X), rho) / (nrow(X) - 1)
  coefficients = cov.inv %*% crossprod(X, Y)
  return (structure(list(coefficients = coefficients),
                    class = "fastVAR.graphicalLm"))
}

coef.fastVAR.GVAR = function(GVAR, ...) {
  return (GVAR$coefficients)
}
