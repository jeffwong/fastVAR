#' @export
fastMlm = function(X, Y, weights = NULL) {
    if (!is.null(weights)) {X = scale.rows(X, sqrt(weights)); Y = scale.rows(Y, sqrt(weights))}
    mlm = fast_mlm(X, Y)
    rownames(mlm$coefficients) = colnames(X)
    colnames(mlm$coefficients) = colnames(Y)
    return (structure(mlm, class = "fastVAR.fastMlm"))
}

#' @method coef fastVAR.fastMlm
#' @S3method coef fastVAR.fastMlm
coef.fastVAR.fastMlm = function(fastMlm) {
    return (fastMlm$coefficients)
}
