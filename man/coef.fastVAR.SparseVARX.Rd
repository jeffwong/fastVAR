\name{coef.fastVAR.SparseVARX}
\alias{coef.fastVAR.SparseVARX}
\title{Coefficients of a SparseVARX model}
\usage{
  \method{coef}{fastVAR.SparseVARX} (sparseVARX, l1penalty)
}
\arguments{
  \item{sparseVARX}{an object of class fastVAR.SparseVARX}

  \item{l1penalty}{The l1 penalty to be applied to the
  SparseVARX model.}
}
\description{
  The underlying library, glmnet, computes the full path to
  the lasso. This means it is computationally easy to
  compute the lasso solution for any penalty term.  This
  function allows you to pass in the desired l1 penalty and
  will return the coefficients
}

