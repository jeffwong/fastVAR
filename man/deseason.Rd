\name{deseason}
\alias{deseason}
\title{Deseason}
\usage{
  deseason(mts, frequency = rep(NA, ncol(mts)), auto = F)
}
\arguments{
  \item{mts}{either a vector representing 1 time series or
  a a data frame or matrix representing multiple time
  series}

  \item{frequency}{the number of observations per period -
  either a single numeric in the univariate case or a
  vector in the multivariate case}

  \item{auto}{if true try to use the spectral density to
  estimate the frequency}
}
\description{
  Deseasonalize a time series so that models can target the
  unexplainable components.
}

