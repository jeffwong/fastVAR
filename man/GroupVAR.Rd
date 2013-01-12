\name{GroupVAR}
\alias{GroupVAR}
\title{Group Vector Autoregression via Group Lasso}
\usage{
  GroupVAR(y, freq = rep(NA, ncol(y)), p = 1,
    weights = NULL, alpha = 0.4, getdiag = T)
}
\arguments{
  \item{y}{A matrix where each column represents an
  individual time series}

  \item{freq}{only used if the time series are periodic.
  freq is a vector of frequencies for each of the time
  series, as in 'ts(y, freq = ...)'.  If the time series
  are not periodic, then this vector can be a vector of NA}

  \item{p}{the number of lags to include in the design
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
  the VAR response matrix omits the first p responses from
  the input matrix Y.  This is because observations in Y
  are modeled by the p previous values, so the first p
  observations cannot be modeled.
}
\details{
  While multivariate response regressions can be solved as
  multiple univariate response regressions, this
  multivariate response problem can better be solved by
  using Group Lasso.  Instead of seeking sparsity in the
  coefficients for each univariate response, Group Lasso
  attempts to find sparsity in the coefficient matrix as a
  whole.
}
\examples{
data(Canada)
  GroupVAR(Canada, p = 3)
}

