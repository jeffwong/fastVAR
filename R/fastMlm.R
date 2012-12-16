#' @export
fastMlm = function(X, Y, weights = NULL) {
    if (!is.null(weights)) {X = scale.rows(X, sqrt(weights)); Y = scale.rows(Y, sqrt(weights))}
    mlm = fast_mlm(X, Y)
    return (structure(mlm, class = "fastVAR.fastMlm"))
}

coef.fastVAR.fastMlm = function(fastMlm) {
    return (fastMlm$coefficients)
}
