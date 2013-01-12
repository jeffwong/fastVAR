\name{SparseVARX}
\alias{SparseVARX}
\title{Sparse Vector Autoregression with Exogenous Inputs}
\usage{
  SparseVARX(y, freq = rep(NA, ncol(y)), x, p, b,
    y.spec = matrix(1, nrow = ncol(y), ncol = ncol(y)),
    x.spec = matrix(1, nrow = ncol(y), ncol = ncol(x)),
    numcore = 1, alpha = 0.8, ...)
}
\arguments{
  \item{y}{A matrix of endogenous variables where each
  column represents an individual time series}

  \item{freq}{only used if the time series are periodic.
  freq is a vector of frequencies for each of the time
  series, as in 'ts(y, freq = ...)'.  If the time series
  are not periodic, then this vector can be a vector of NA}

  \item{x}{A matrix of exogenous variables where each
  column represents an individual time series}

  \item{p}{the number of lags of Y to include in the design
  matrix}

  \item{b}{the number of lags to X include in the design
  matrix}

  \item{y.spec}{A binary matrix that can constrain the
  number of lagged predictor variables.  If y.spec[i][j] =
  0, the ith time series in y will not be regressed on the
  jth time series of y, or any of its lags.}

  \item{x.spec}{A binary matrix that can constrain the
  number of lagged exogenous variables.  If x.spec[i][j] =
  0, the ith time series in y will not be regressed on the
  jth time series of x, or any of its lags.}

  \item{numcore}{number of cpu cores to use to parallelize
  this function}

  \item{alpha}{the elastic net mixing parameter, as defined
  in 'glmnet'}
}
\description{
  Fit a VARX model with lasso penalty. The VAR model is
  estimated using a multiresponse linear regression. The
  sparse VAR fits multiple uniresponse linear regressions
  with lasso penalty. mclapply from multicore can be used
  to fit the individual uniresponse linear regressions in
  parallel.  Note that mclapply is not available for
  windows
}
\examples{
data(Canada)
  x = matrix(rnorm(84*4), 84, 4)
  SparseVARX(Canada, x = x, p = 3, b = 2)
}

