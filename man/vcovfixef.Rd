\name{VcovFixef}
\alias{VcovFixef}
\alias{VcovFixef,ipdlme-method}
\title{Extract the variance-covariance for the population (fixed effect) estimates for an \code{ipdlme} object}
\description{
  Extract the variance-covariance for the fixed effects estimates of the fit of an IPD LME model.
}
\usage{
VcovFixef(object,\dots)
}
\arguments{
 \item{object}{object of the \code{ipdlme} class}
 \item{\dots}{additional arguments, currently not implemented}
}

\value{
Matrix of the variance-covariance for the fitted fixed effects of an IPD LME model.
}

\seealso{\code{\link{ipdlme-class}}}

\examples{
data(regress_chol)

VcovFixef(ipdlme(n, Y, S2))
}

\keyword{models}
