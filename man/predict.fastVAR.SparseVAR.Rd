\name{predict.fastVAR.SparseVAR}
\alias{predict.fastVAR.SparseVAR}
\title{SparseVAR Predict}
\usage{
  predict.fastVAR.SparseVAR(sparseVAR, n.ahead = 1,
    threshold, ...)
}
\arguments{
  \item{sparseVAR}{an object of class fastVAR.SparseVAR
  returned from SparseVAR}

  \item{n.ahead}{number of steps to predict}

  \item{threshold}{threshold prediction values to be
  greater than this value}

  \item{...}{extra parameters to pass into the coefficients
  method for objects of type fastVAR.VAR}
}
\description{
  Predict n steps ahead from a fastVAR.SparseVAR object
}
\examples{
data(Canada)
  predict(SparseVAR(Canada, 3), 1)
}

