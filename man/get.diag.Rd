\name{get.diag}
\alias{get.diag}

\title{Extract Diagonal Elements from Matrix}

\description{
Wrapper for \code{diag} that allows for scalar argument.
}

\usage{
get.diag(x,...)
}

\arguments{

  \item{x}{Matrix or scalar}
  \item{...}{Additional arguments to \code{diag}}

}


\value{
Diagonal as vector.
}

\details{
Same as \code{diag} but interprets scalar as scalar rather than the
elements of an identity matrix. Thus if x is scalar or vector it is
simply returned.
}

\author{S. Kovalchik \email{skoval@ucla.edu}}

\seealso{\code{\link{diag}}}

\examples{

get.diag(matrix(1:4,2,2))

get.diag(4)

}
