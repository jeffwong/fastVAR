\name{predict.fastVAR.SparseVARX}
\alias{predict.fastVAR.SparseVARX}
\title{SparseVARX Predict}
\usage{
  predict.fastVAR.SparseVARX(SparseVARX, xnew, n.ahead = 1,
    threshold, ...)
}
\arguments{
  \item{sparseVARX}{an object of class fastVAR.SparseVARX
  returned from SparseVARX}

  \item{xnew}{a matrix of future values for the exogenous
  inputs.  Should contain n.ahead rows}

  \item{n.ahead}{number of steps to predict}

  \item{threshold}{threshold prediction values to be
  greater than this value}

  \item{...}{extra parameters to pass into the coefficients
  method for objects of type fastVAR.VAR}
}
\description{
  Predict n steps ahead from a fastVAR.SparseVARX object
}
\examples{
data(Canada)
  x = matrix(rnorm(84*4), 84, 4)
  predict(SparseVARX(Canada, x, 3, 2), xnew=matrix(rnorm(2*4),2,4), n.ahead=2)
}

