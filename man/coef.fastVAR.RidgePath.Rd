\name{coef.fastVAR.RidgePath}
\alias{coef.fastVAR.RidgePath}
\title{Ridge Coefficients}
\usage{
  \method{coef}{fastVAR.RidgePath} (model, l2penalty)
}
\arguments{
  \item{ridgePath}{an object of class fastVAR.RidgePath}

  \item{lambda}{the desired l2penalty}
}
\description{
  The full ridge path.  Allows a user to calculate the
  coefficients of a linear model for any l2 penalty
}

