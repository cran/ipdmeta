\name{cluster.subvector.matrix}
\alias{cluster.subvector.matrix}

\title{Construct matrix from grouped vector}

\description{
Given a matrix organized by effect from some 1 to \code{n.groups} clusters, returns a matrix with each row corresponding to a cluster.
}

\usage{
cluster.subvector.matrix(b,n.groups)
}

\arguments{

  \item{b}{Vector of clustered values stacked by effect, within each effect ordered form 1 to \code{n.groups} groups.}
  \item{n.groups}{Number of clusters}

}

\value{

Matrix with row a cluster and column 1 to the number of effects.

}


\author{S. Kovalchik \email{skoval@ucla.edu}}

\examples{

#CONSIDER 10 GROUPS WITH BIVARIATE FRAILTY 

frailties <- 1:20
cluster.subvector.matrix(frailties,10)

}
