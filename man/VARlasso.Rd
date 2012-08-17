\name{VARlasso}
\alias{VARlasso}
\title{
Vector Autoregressive Model with Lasso Penalty
}
\description{
Fit a VAR model with lasso penalty to multivariate response data
}
\usage{
VARlasso(y, p, lambda=NULL, y.spec = matrix(1, nrow = ncol(y), ncol = ncol(y)), getdiag = T, numcore = 1, ...)
}
\arguments{
  \item{y}{
Response matrix where each column is a separate time series
}
  \item{p}{
The autoregressive lag of the VAR(p) model
}
  \item{lambda}{
A vector where the ith element represents the L1 penalty on the ith time series.  
If NULL, fastVAR uses glmnet's cv.glmnet function to choose lambda
}
  \item{y.spec}{
A binary matrix that can constrain the number of lagged predictor variables.  
If y.spec[i][j] = 0, the ith time series in y will not be regressed on the jth
time series of y, or any of its lags.
}
  \item{getdiag}{
Logical.  If true, report diagnostics
}
  \item{numcore}{
The number of cpu cores that can be used to evaluate this function
}
  \item{\dots}{
Extra arguments that are passed to mclapply from the 'multicore' package
}
}
\details{
Creates a design matrix of p lags, then applies the cv.glmnet function of the 
'glmnet' package to run a linear regression on each time series of y with 
lasso penalty.  The result is a very sparse coefficient matrix.
}
\value{
If getdiag is false, return a sparse coefficient matrix of the multivariate response regression.
If getdiag is true, return a list with components:
\item{B }{The sparse coefficient matrix of the multivariate regression}
\item{diag }{A list of diagnostics.  See VARdiag}
}
\references{
Robert Tibshirani, Stanford University <http://www-stat.stanford.edu/~tibs/lasso.html>
}
\author{
Jeffrey Wong
}
\section{errors}{
  It is possible that the regression returns a coefficient matrix with NAs in it.
  If so, the diagnostic calculations cannot be done, and an error will return.
  To explore this, one can rerun the function with getdiag=F.

  To run VARlasso with multiple cores requires the 'multicore' package to be
  installed.
}
\seealso{
fastVAR, fastVARX, VARXlasso
}
\examples{
data(Canada)
VARlasso(Canada,3, getdiag=FALSE)
}
\keyword{ lm }
\keyword{ VAR }
