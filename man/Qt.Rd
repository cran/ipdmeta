\name{Qt}
\alias{Qt}

\title{Measures of covariate heterogeneity}

\description{
Measures of covariate heterogeneity proposed by Simmonds and Higgins (2007) for assessing the power of a meta-regression
}

\usage{
Qt(m,n,sigma2)
}

\arguments{

  \item{m}{vector of study-level covariate means}
  \item{n}{vector of study sample sizes}
  \item{sigma2}{vector of covariate sample variances}
}

\value{

A list with the following named components: \code{t}, \code{Qd}, \code{Qe}, \code{bar.Qd},\code{bar.Qe}, \code{tilde.Qd}, \code{tilde.Qe}

}

\references{

 Simmonds, M. C., Higgins, J. P. T., (2007), Covariate heterogeneity in meta-analysis:
  criteria for deciding between meta-regression and individual patient data, \emph{Statistics in Medicine},  26 (15): 2982-99.

}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(poynard)

#COVARIATE HETEROGENEITY FOR AGE

m <- (poynard$n0*poynard$age0+poynard$n1*poynard$age1)/(poynard$n0+poynard$n1)
n <- poynard$n0+poynard$n1
sigma2 <- ((poynard$n0-1)*poynard$age.s20+(poynard$n1-1)*poynard$age.s21)/(poynard$n0+poynard$n1-2)


Q <- Qt(m,n,sigma2)

lapply(Q,function(x){x/Q$t})

}
