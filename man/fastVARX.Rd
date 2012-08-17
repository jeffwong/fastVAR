\name{fastVARX}
\alias{fastVARX}
\title{
Vector Autoregressive Model with Exogenous Inputs
}
\description{
Fit a VARX to multivariate response data
}
\usage{
fastVARX(y, x, p, b, getdiag = T)
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
  \item{getdiag}{
Logical.  If true, report diagnostics
}
}
\details{
Creates a design matrix using p lags of y and b lags of x,
then runs a multivariate response regression.
}
\value{
If getdiag is false, return an object of lm class representing the multivariate response regression.
If getdiag is true, return a list with components:
\item{model }{The lm object representing the multivariate regression}
\item{diag }{A list of diagnostics.  See VARdiag}
}
\author{
Jeffrey Wong
}
\section{errors}{
  It is possible that the regression returns a coefficient matrix with NAs in it.
  If so, the diagnostic calculations cannot be done, and an error will return.
  To explore this, one can rerun the function with getdiag=F.
}

\seealso{
fastVAR, VARlasso, VARXlasso
}
\examples{
data(Canada)
x = matrix(rnorm(84*3),84,3)
fastVARX(Canada,x,3,2,getdiag=FALSE)
}
\keyword{ lm }
\keyword{ VAR }
