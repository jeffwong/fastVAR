\name{VARdiag}
\alias{VARdiag}
\title{
VAR and VARX Diagnostics
}
\description{
Report diagnostics such as standard errors, t values, p values, AIC and BIC.
This function is primarily used for internal purposes, but can be used as a 
public function.
}
\usage{
VARdiag(y.p, Z, B, n, T, k, p, dof)
}
\arguments{
  \item{y.p}{
The reduced response matrix.  Typically obtained by calling VARZ or VARXZ
}
  \item{Z}{
The design matrix.  Typically obtained by calling VARZ or VARXZ
}
  \item{B}{
The matrix of coefficients.  Typically obtained by fastVAR, fastVARX, VARlasso or VARXlasso
}
  \item{n}{
The number of endogenous time series that are being measured
}
  \item{T}{
The number of time periods of each time series
}
  \item{k}{
The number of variables of the design matrix
}
  \item{p}{
The autoregressive lag of the response variables
}
  \item{dof}{
The degrees of freedom for the residuals
}
}
\details{
Reports standard regression diagnostics, such as the standard errors of the coefficients,
their t values and p values, as well as the AIC and BIC of the model.  Warning, this function
is slightly buggy and does not measure the AIC and BIC accurately.
}
\value{
\item{se }{A matrix of the standard errors of the coefficients}
\item{tvalue }{A matrix of the t values of the coefficients}
\item{pvalue }{A matrix of the p values of the coefficients}
\item{residual }{A matrix of the residuals}
\item{cov }{The covariance matrix of the residuals}
\item{AIC }{The AIC of the model}
\item{BIC }{The BIC of the model}
\item{HQ }{The HQ of the model}
}
\author{
Jeffrey Wong
}
\seealso{
VARXdiag
}
\examples{
data(Canada)
canada.fit = fastVAR(Canada,3,getdiag=TRUE)
canada.fit$diag
}
