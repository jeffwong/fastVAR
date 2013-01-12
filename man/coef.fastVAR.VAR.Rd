\name{coef.fastVAR.VAR}
\alias{coef.fastVAR.VAR}
\title{VAR Coefficients}
\usage{
  \method{coef}{fastVAR.VAR} (VAR, ...)
}
\arguments{
  \item{VAR}{an object of class fastVAR.VAR}

  \item{...}{if VAR was fit using a l2 penalty, the user
  can specify a different l2 penalty here and have the
  coefficients recomputed}
}
\value{
  The coefficients for the VAR model
}
\description{
  If the VAR object was fit using a l2 penalty, then the
  full ridge path was calculated and stored in the object.
  This means the user can adjust the ridge penalty term
  here and recompute the coefficients of the VAR
}

