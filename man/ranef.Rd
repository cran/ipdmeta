\name{ranef}
\alias{ranef}
\alias{ranef,ipdlme-method}
\title{Extract the study random effects for an \code{ipdlme} object}
\description{
  Extract random effects estimates of the fit of an IPD LME model.
}
\usage{
ranef(object,\dots)
}
\arguments{
 \item{object}{object of the \code{ipdlme} class}
 \item{\dots}{additional arguments, currently not implemented}
}

\value{
Numeric vector of the fitted random effects of an IPD LME model.
}

\seealso{\code{\link{ipdlme-class}}}

\examples{
data(regress_chol)

ranef(ipdlme(n, Y, S2))
}

\keyword{models}
