\name{fastVAR}
\alias{fastVAR}
\title{
Vector Autoregressive Model
}
\description{
Fit a VAR to multivariate response data
}
\usage{
fastVAR(y, p = 1, getdiag = T)
}
\arguments{
  \item{y}{
Response matrix where each column is a separate time series
}
  \item{p}{
The autoregressive lag of the VAR(p) model
}
  \item{getdiag}{
Logical.  If true, report diagnostics
}
}
\details{
Uses multivariate OLS to fit the response data to p lags of itself.
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
fastVARX, VARlasso, VARXlasso
}
\examples{
data(Canada)
fastVAR(Canada,3)$model$coefficients
}
\keyword{ lm }
\keyword{ VAR }
