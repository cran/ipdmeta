\name{which.formula}
\alias{which.formula}

\title{
Formula Comparison
}

\description{
Give index of shared terms of formulae. Returns the index of the
shared components in terms of the vector of terms for \code{f1}.
}

\usage{
which.formula(f1,f2)
}

\arguments{

\item{f1}{reference formula}
\item{f2}{formula}

} 

\value{

Vector of shared index.

}

\author{
S. A. Kovalchik <skoval@ucla.edu>
}

\examples{

f1 <- ~x*y
f2 <- ~y

which.formula(f1,f2)

}

