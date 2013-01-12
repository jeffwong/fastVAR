\name{predict.fastVAR.GroupVARX}
\alias{predict.fastVAR.GroupVARX}
\title{GroupVARX Predict}
\usage{
  \method{predict}{fastVAR.GroupVARX} (GroupVARX, xnew,
    n.ahead = 1, threshold, ...)
}
\arguments{
  \item{GroupVARX}{an object of class fastVAR.GroupVARX
  returned from GroupVARX}

  \item{xnew}{a matrix of future values for the exogenous
  inputs.  Should contain n.ahead rows}

  \item{n.ahead}{number of steps to predict}

  \item{threshold}{threshold prediction values to be
  greater than this value}

  \item{...}{extra parameters to pass into the coefficients
  method for objects of type fastVAR.GroupVARX}
}
\description{
  Predict n steps ahead from a fastVAR.GroupVARX object
}
\examples{
data(Canada)
  x = matrix(rnorm(84*4), 84, 4)
  predict(GroupVARX(Canada, x = x, p = 3, b = 2), xnew = matrix(rnorm(2*4),2,4), n.ahead = 2)
}

