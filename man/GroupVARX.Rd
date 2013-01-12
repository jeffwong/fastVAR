\name{GroupVARX}
\alias{GroupVARX}
\title{Group Vector Autoregression with Exogenous Inputs via Group Lasso}
\usage{
  GroupVARX(y, freq = rep(NA, ncol(y)), x, p = 1, b = 1,
    weights = NULL, alpha = 0.4, getdiag = T)
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

  \item{weights}{weights applied to the multiresponse
  linear regression.  Better predictions might come from
  weighting observations far in the past less so they
  impact the objective value less.  Either NULL,
  "exponential", or "linear"}

  \item{alpha}{the mixing parameter between group lasso and
  quadratic, as in 'glmnet'}

  \item{getdiag}{logical.  If true, return diagnostics}
}
\description{
  Fit a VAR model by creating the lagged design matrix and
  fitting a multivariate response matrix to it.  Note that
  the GroupVARX response matrix omits the first max(p,b)
  responses from the input matrix Y.  This is because
  observations in Y are modeled by the max(p,b) previous
  values, so the first max(p,b) observations cannot be
  modeled. While multivariate response regressions can be
  solved as multiple univariate response regressions, this
  multivariate response problem can better be solved by
  using Group Lasso.  Instead of seeking sparsity in the
  coefficients for each univariate response, Group Lasso
  attempts to find sparsity in the coefficient matrix as a
  whole.
}
\examples{
data(Canada)
  x = matrix(rnorm(84*4), 84, 4)
  GroupVARX(y, x=x, p=3, b=2)
}

