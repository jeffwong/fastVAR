\name{fastVAR-package}
\alias{fastVAR-package}
\alias{fastVAR-package}
\docType{package}
\title{
Compute large VAR and VARX models
}
\description{
This package is designed for time series data.  Uses fast implementations to 
estimate Vector Autoregressive models and Vector Autoregressive models with 
Exogenous Inputs.  For speedup, fastVAR can use multiple cpu cores to
calculate the estimates.  For very large systems, fastVAR uses Lasso penalty to 
return very sparse coefficient matrices.  Regression diagnostics can be used 
to compare models, and prediction functions can be used to calculate the 
n-step ahead prediction.

Includes Canada data set from package vars to validate results.
}
\details{
\tabular{ll}{
Package: \tab fastVAR\cr
Type: \tab Package\cr
Version: \tab 1.9.9\cr
Date: \tab 2012-09-30\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
Depends: glmnet
Suggests: multicore
}

Very few functions:
\itemize{
  \item VAR - estimate a standard VAR model
  \item VARX - estimate a standard VARX model with exogenous inputs
  \item SparseVAR - use lasso penalty to estimate a VAR model
  \item SparseVARX - use lasso penalty to estimate a VARX model
  \item VAR.diag - get regression diagnostics for a VAR or VARX model
}

VAR, VARX, SparseVAR, and SparseVARX have accompanying prediction methods, which
can be invoked by predict(one of the four VAR models, n.ahead)

}
\author{
Jeffrey Wong <jeff.ct.wong@gmail.com>
}
\references{
Robert Tibshirani <http://www-stat.stanford.edu/~tibs/lasso.html>
}
\keyword{ package }
\seealso{
glmnet, vars
}
