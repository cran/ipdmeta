\name{frailty.model.matrix}
\alias{frailty.model.matrix}

\title{Design matrix for random effects in Cox mixed effects model}

\description{
Returns matrix of n by number of cluster effects design matrix.
}

\usage{
frailty.model.matrix(random,data)
}

\arguments{

  \item{random}{formula for the random components as would be stated in \code{coxme}}
  \item{data}{Data frame containing the observed outcomes}
}


\value{

Matrix of cluster effect design matrices.

}


\author{S. Kovalchik \email{skoval@ucla.edu}}

\examples{

data(cancer.rats)

#RANDOM EFFECTS "DESIGN MATRIX"
#BASELINE FRAILTY, 50 LITTERS 3 RATS IN EACH LITTER

frailty.model.matrix(~(1|litter),cancer.rats)

}
