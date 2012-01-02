\name{plot}
\alias{plot}
\alias{plot,ipdlme-method}
\title{Forest plot of study random effects for an \code{ipdlme} object}
\description{
Forest plot of the study random effects.
}
\usage{
plot(x,y,\dots)
}
\arguments{
 \item{x}{object of the \code{ipdlme} class}
 \item{y}{vector of label names for studies. This will be used for labels along the y-axis.}
 \item{\dots}{additional arguments passed to \code{graphics} \code{plot} function}
}

\value{
Side-by-side forest plots with intercept effects on the left and treatment effects on the right.
}

\seealso{\code{\link{ipdlme-class}}}

\examples{
data(regress_chol)

plot(ipdlme(n, Y, S2))
}

\keyword{models}
