\name{fixef}
\alias{fixef}
\alias{fixef,ipdlme-method}
\title{Extract the population (fixed effect) estimates for an \code{ipdlme} object}
\description{
  Extract fixed effects estimates of the fit of an IPD LME model.
}
\usage{
fixef(object,\dots)
}
\arguments{
 \item{object}{object of the \code{ipdlme} class}
 \item{\dots}{additional arguments, currently not implemented}
}

\value{
Numeric vector of the fitted fixed effects of an IPD LME model.
}

\seealso{\code{\link{ipdlme-class}}}

\examples{
data(regress_chol)

fixef(ipdlme(n, Y, S2))
}

\keyword{models}
