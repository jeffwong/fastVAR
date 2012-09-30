\name{predict.fastVAR.VAR}
\alias{predict.fastVAR.VAR}
\title{VAR Predict}
\usage{
  predict.fastVAR.VAR(VAR, n.ahead = 1, threshold, ...)
}
\arguments{
  \item{VAR}{an object of class fastVAR.VAR returned from
  VAR}

  \item{n.ahead}{number of steps to predict}

  \item{threshold}{threshold prediction values to be
  greater than this value}

  \item{...}{extra parameters to pass into the coefficients
  method for objects of type fastVAR.VAR}
}
\description{
  Predict n steps ahead from a fastVAR.VAR object
}
\examples{
data(Canada)
  predict(VAR(Canada, 3, intercept=F), 1)
}

