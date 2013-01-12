\name{coef.fastVAR.SparseVAR}
\alias{coef.fastVAR.SparseVAR}
\title{Coefficients of a SparseVAR model}
\usage{
  \method{coef}{fastVAR.SparseVAR} (sparseVAR, l1penalty)
}
\arguments{
  \item{sparseVAR}{an object of class fastVAR.SparseVAR}

  \item{l1penalty}{The l1 penalty to be applied to the
  SparseVAR model.}
}
\description{
  The underlying library, glmnet, computes the full path to
  the lasso. This means it is computationally easy to
  compute the lasso solution for any penalty term.  This
  function allows you to pass in the desired l1 penalty and
  will return the coefficients
}

