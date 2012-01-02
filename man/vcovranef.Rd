\name{VcovRanef}
\alias{VcovRanef}
\alias{VcovRanef,ipdlme-method}
\title{Extract the variance-covariance matrix for study random effect estimates for an \code{ipdlme} object}
\description{
  Extract the variance-covariance for the estimated BLUPs for the study baseline and treatment terms for the fit of an IPD LME model.
}
\usage{
VcovRanef(object,\dots)
}
\arguments{
 \item{object}{object of the \code{ipdlme} class}
 \item{\dots}{additional arguments, currently not implemented}
}

\value{
Matrix of the variance-covariance for the study randoms effects of an IPD LME model.
}

\seealso{\code{\link{ipdlme-class}}}

\examples{
data(regress_chol)

VcovRanef(ipdlme(n, Y, S2))
}

\keyword{models}
