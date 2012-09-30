\name{SparseVAR}
\alias{SparseVAR}
\title{Sparse Vector Autoregression}
\usage{
  SparseVAR(y, p, intercept = F,
    y.spec = matrix(1, nrow = ncol(y), ncol = ncol(y)),
    numcore = 1, ...)
}
\arguments{
  \item{y}{A matrix where each column represents an
  individual time series}

  \item{p}{the number of lags to include in the design
  matrix}

  \item{intercept}{logical.  If true, include an intercept
  term in the model}

  \item{y.spec}{A binary matrix that can constrain the
  number of lagged predictor variables.  If y.spec[i][j] =
  0, the ith time series in y will not be regressed on the
  jth time series of y, or any of its lags.}

  \item{numcore}{number of cpu cores to use to parallelize
  this function}
}
\description{
  Fit a vector autoregressive model with lasso penalty. The
  VAR model is estimated using a multiresponse linear
  regression. The sparse VAR fits multiple uniresponse
  linear regressions with lasso penalty. mclapply from
  multicore can be used to fit the individual uniresponse
  linear regressions in parallel.  Note that mclapply is
  not available for windows
}

