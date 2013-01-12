\name{predict.fastVAR.GroupVAR}
\alias{predict.fastVAR.GroupVAR}
\title{GroupVAR Predict}
\usage{
  \method{predict}{fastVAR.GroupVAR} (GroupVAR, n.ahead,
    threshold, ...)
}
\arguments{
  \item{GroupVAR}{an object of class fastVAR.GroupVAR
  returned from GroupVAR}

  \item{n.ahead}{number of steps to predict}

  \item{threshold}{threshold prediction values to be
  greater than this value}

  \item{...}{extra parameters to pass into the coefficients
  method for objects of type fastVAR.GroupVAR}
}
\description{
  Predict n steps ahead from a fastVAR.GroupVAR object
}
\examples{
data(Canada)
  predict(GroupVAR(Canada, p = 3), 1)
}

