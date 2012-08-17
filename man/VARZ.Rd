\name{VARZ}
\alias{VARZ}
\title{
Design Matrix for VAR model
}
\description{
Determine the design matrix of a VARX model
}
\usage{
VARZ(y, p)
}
\arguments{
  \item{y}{
Response matrix where each column is a separate time series
}
  \item{p}{
The autoregressive lag of the response variable
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
\item{n}{Number of endogenous time series that are being measured}
\item{T }{The number of time points in the reduced response matrix}
\item{k }{The total number of predictor variables used to model each endogenous time series}
\item{dof }{The degrees of freedom of the residuals}
\item{y.p }{The reduced response matrix}
\item{Z }{The design matrix}

}
\author{
Jeffrey Wong
}
\seealso{
VARXZ
}
\examples{
data(Canada)
VARZ(Canada,3)
}
