\name{coef.fastVAR.VARX}
\alias{coef.fastVAR.VARX}
\title{VARX Coefficients}
\usage{
  \method{coef}{fastVAR.VARX} (VARX, ...)
}
\arguments{
  \item{VARX}{an object of class fastVAR.VARX}

  \item{...}{if VAR was fit using a l2 penalty, the user
  can specify a different l2 penalty here and have the
  coefficients recomputed}
}
\value{
  The coefficients for the VARX model
}
\description{
  If the VARX object was fit using a l2 penalty, then the
  full ridge path was calculated and stored in the object.
  This means the user can adjust the ridge penalty term
  here and recompute the coefficients of the VARX
}

