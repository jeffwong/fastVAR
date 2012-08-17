\name{VARpredict}
\alias{VARpredict}
\title{
Prediction with VAR model
}
\description{
Return the n.ahead-step ahead prediction of the VAR model
}
\usage{
VARpredict(y, p, B, n.ahead = 1, threshold)
}
\arguments{
  \item{y}{
Response matrix where each column is a separate time series
}
  \item{p}{
The autoregressive lag of the VAR(p) model
}
  \item{B}{
The coefficient matrix typically calculated by fastVAR, fastVARX, VARlasso or VARXlasso
}
  \item{n.ahead}{
The number of steps ahead to predict
}
  \item{threshold}{
Optional.  If any predicted values fall below threshold, replace them with threshold
}
}
\details{
Use the coefficient matrix to simulate the next n.ahead values.  

Sometimes it is known that values cannot fall below 0.  Set threshold=0 to
get predictions that are always positive.  Any predicted value that was below 0
will be set to 0.
}
\value{
A matrix where the rows represent the n.ahead steps, and the columns represent 
the predicted values for each time series.
}
\author{
Jeffrey Wong
}
\seealso{
VARXpredict
}
\examples{
data(Canada)
B = fastVAR(Canada,3,getdiag=FALSE)$coefficients
VARpredict(Canada, 3, B, threshold=0)
}
\keyword{ predict }
