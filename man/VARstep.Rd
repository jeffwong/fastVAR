\name{VARstep}
\alias{VARstep}
\title{
Model Selection for VAR(p)
}
\description{
Step-like function that measures the BIC of a VAR model at various values of p
}
\usage{
VARstep(y, max.p)
}
\arguments{
  \item{y}{
Response matrix where each column is a separate time series
}
  \item{max.p}{
The maximum lag value to use
}
}
\details{
Similar to R's step function for determining what variables to use in a model.
This funciton uses BIC to determine what lag order, p, should be used.
}
\value{
The optimal value for p
}
\author{
Jeffrey Wong
}
\examples{
data(Canada)
VARstep(Canada, max.p=10)
}
