\name{VARXZ}
\alias{VARXZ}
\title{
Design Matrix for VARX model
}
\description{
Determine the design matrix of a VARX model
}
\usage{
VARXZ(y, x, p, b)
}
\arguments{
  \item{y}{
Response matrix where each column is a separate time series
}
  \item{x}{
Matrix of exogenous inputs, where each column is a separate exogenous variable
}
  \item{p}{
The autoregressive lag of the response variable
}
  \item{b}{
The autoregressive lag of the exogenous variables
}
}
\details{
Creates a design matrix using p lags of y and b lags of x.

Creates the reduced response matrix, which is a subset of the y matrix.
In the VAR model, each point in time is modeled after values from previous
points in time.  In a VAR(p) model, the first p time points cannot be modeled
as they require data from p periods ago.  Likewise, in a VARX(p,b) model,
the first max(p,b) points cannot be modeled.  The reduced response matrix
is the matrix with these points dropped.
}
\value{
\item{ny }{Number of endogenous time series that are being measured}
\item{nx }{Number of exogenous time series that are being used for prediction}
\item{T }{The number of time points in the reduced response matrix}
\item{k }{The total number of predictor variables used to model each endogenous time series}
\item{p.max }{The greater of p and b}
\item{dof }{The degrees of freedom of the residuals}
\item{y.p }{The reduced response matrix}
\item{Z }{The design matrix}
}
\author{
Jeffrey Wong
}
\seealso{
VARZ
}
\examples{
data(Canada)
x = matrix(rnorm(84*3),84,3)
VARXZ(Canada, x, 3, 2)
}
