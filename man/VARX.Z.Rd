\name{VARX.Z}
\alias{VARX.Z}
\title{VARX Design Matrix}
\usage{
  VARX.Z(y, x, p, b, intercept = F)
}
\arguments{
  \item{y}{a multivariate time series matrix where rows
  represent time and columns represent different time
  series}

  \item{y}{a multivariate time series matrix where rows
  represent time and columns represent different time
  series}

  \item{p}{the lag order of the y matrix}

  \item{b}{the lag order of the x matrix}

  \item{intercept}{logical.  If true, include a 1 vector
  for the intercept term}
}
\value{
  \item{n}{Number of endogenous time series that are being
  measured} \item{T }{The number of time points in the
  reduced response matrix} \item{k }{The total number of
  predictor variables used to model each endogenous time
  series} \item{dof }{The degrees of freedom of the
  residuals} \item{y.p }{The reduced response matrix}
  \item{Z }{The design matrix} \item{y.orig}{The original
  input matrix} \item{p}{The lag order of the y matrix}
  \item{x.orig}{The original input matrix} \item{b}{The lag
  order of the x matrix} \item{intercept}{logical.  If
  true, include a 1 vector for the intercept term}
}
\description{
  Compute the Design Matrix for a standard VAR model with
  exogenous inputs
}

