\name{sigma2}
\alias{sigma2}
\alias{sigma2,ipdlme-method}
\title{Extract the estimated residual variance for an \code{ipdlme} object}
\description{
  Extract the estimated residual variance for the fit of an IPD LME model.
}
\usage{
sigma2(object)
}
\arguments{
 \item{object}{object of the \code{ipdlme} class}
}

\value{
Numeric, maximum likelihood estimate for the residual variance of the fitted IPD LME model.
}

\seealso{\code{\link{ipdlme-class}}}

\examples{
data(regress_chol)

sigma2(ipdlme(n, Y, S2))
}

\keyword{models}
