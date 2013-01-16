\name{predict.fastVAR.GraphicalVARX}
\alias{predict.fastVAR.GraphicalVARX}
\title{GraphicalVARX Predict}
\usage{
  predict.fastVAR.GraphicalVARX(GraphicalVARX, xnew,
    n.ahead = 1, threshold, ...)
}
\arguments{
  \item{GraphicalVARX}{an object of class
  fastVAR.GraphicalVARX returned from GraphicalVARX}

  \item{xnew}{a matrix of future values for the exogenous
  inputs.  Should contain n.ahead rows}

  \item{n.ahead}{number of steps to predict}

  \item{threshold}{threshold prediction values to be
  greater than this value}

  \item{...}{extra parameters to pass into the coefficients
  method for objects of type fastVAR.GraphicalVARX}
}
\description{
  Predict n steps ahead from a fastVAR.GraphicalVARX object
}
\examples{
data(Canada)
  x = matrix(rnorm(84*4), 84, 4)
  predict(GraphicalVARX(Canada, x = x, p = 3, b = 2, intercept = F), xnew = matrix(rnorm(2*4),2,4), n.ahead = 2)
}

