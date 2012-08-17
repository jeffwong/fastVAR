\name{VARXpredict}
\alias{VARXpredict}
\title{
Prediction with VARX model
}
\description{
Return the n.ahead-step ahead prediction of the VARX model
}
\usage{
VARXpredict(y, x, p, b, B, xnew, n.ahead = 1, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Response matrix where each column is a separate time series
}
  \item{x}{
Matrix of exogenous variables where each column is a separate time series
}
  \item{p}{
The autoregressive lag of the VAR(p) model
}
  \item{b}{
The autoregressive lag of the exogenous variables
}
  \item{B}{
The coefficient matrix typically calculated by fastVAR, fastVARX, VARlasso or VARXlasso
}
  \item{xnew}{
The new values for the exogenous variables
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
VARpredict
}
\examples{
data(Canada)
x = matrix(rnorm(84*3),84,3)
B = fastVARX(Canada,x,3,2,getdiag=FALSE)$coefficients
VARXpredict(Canada, x, 3, 2, B, matrix(rnorm(3),3,1), threshold=0)
}
\keyword{ predict }
