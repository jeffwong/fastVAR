\name{predict.fastVAR.VARX}
\alias{predict.fastVAR.VARX}
\title{VARX Predict}
\usage{
  \method{predict}{fastVAR.VARX} (VARX, xnew, n.ahead = 1,
    threshold, ...)
}
\arguments{
  \item{VARX}{an object of class fastVAR.VARX returned from
  VARX}

  \item{xnew}{a matrix of future values for the exogenous
  inputs.  Should contain n.ahead rows}

  \item{n.ahead}{number of steps to predict}

  \item{threshold}{threshold prediction values to be
  greater than this value}

  \item{...}{extra parameters to pass into the coefficients
  method for objects of type fastVAR.VARX}
}
\description{
  Predict n steps ahead from a fastVAR.VARX object
}
\examples{
data(Canada)
  x = matrix(rnorm(84*4), 84, 4)
  predict(VARX(Canada, x = x, p = 3, b = 2, intercept = F), xnew = matrix(rnorm(2*4),2,4), n.ahead = 2)
}

