\name{VAR}
\alias{VAR}
\title{Vector Autoregression}
\usage{
  VAR(y, freq = rep(NA, ncol(y)), p = 1, intercept = T,
    weights = NULL, l2penalty = NULL, getdiag = T)
}
\arguments{
  \item{y}{A matrix where each column represents an
  individual time series}

  \item{p}{the number of lags to include in the design
  matrix}

  \item{intercept}{logical.  If true, include an intercept
  term in the model}

  \item{weights}{weights applied to the multiresponse
  linear regression.  Better predictions might come from
  weighting observations far in the past less so they
  impact the objective value less.  Either NULL,
  "exponential", or "linear"}

  \item{l2penalty}{a ridge regression penalty, useful when
  the design matrix is very wide, which may happen if p is
  large.}

  \item{getdiag}{logical.  If true, return diagnostics}
}
\description{
  Fit a VAR model by creating the lagged design matrix and
  fitting a multivariate response matrix to it.  Note that
  the VAR response matrix omits the first p responses from
  the input matrix Y.  This is because observations in Y
  are modeled by the p previous values, so the first p
  observations cannot be modeled.
}
\examples{
data(Canada)
  VAR(Canada, p = 3, intercept = F)
}

