\name{confint}
\alias{confint}
\alias{confint,ipdlme-method}
\title{Compute confidence intervals for population effects of an \code{ipdlme} object}
\description{
  Computes Wald-type confidence intervals for all fixed effects of the fit of an IPD LME model.
}
\usage{
confint(object,parm,level=0.95,\dots)
}
\arguments{
 \item{object}{object of the \code{ipdlme} class}
 \item{parm}{term names, currently ignored}
 \item{level}{numeric value of confidence level to be used}
 \item{\dots}{additional arguments, currently not implemented}
}

\value{
Matrix of the estimate, lower and upper confidence intervals for the fixed effects of the IPD LME model.
}

\seealso{\code{\link{ipdlme-class}}}

\examples{
data(regress_chol)

confint(ipdlme(n, Y, S2))
}

\keyword{models}
