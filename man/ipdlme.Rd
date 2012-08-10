\name{ipdlme}
\alias{ipdlme}
\title{IPD Linear Mixed-Effects Models from Aggregate Data}
\description{
Fits the IPD linear mixed-effects model for a treatment and single categorical covariate and their interaction using study-level summary statistics.
}
\usage{
ipdlme(n, y, s2, max.iter = 100, tol = 1e-10,equal.residual.var=TRUE) 
}
\arguments{
  \item{n}{a list ordered by trial with each element a data frame having the columns \code{trt} and \code{ctrl} and the row the sample sizes in each group for the covariate of interest, e.g. males/females, never/former/current smokers.}
  \item{y}{a list ordered by trial with each element a data frame having the columns \code{trt} and \code{ctrl} and the row the average outcome (the response variable of the model) in each group for the covariate of interest.
}
\item{s2}{A vector of the study overall sample variances for the outcome or a list of the treatment-subgroup sample variances following the structure of \code{n} and \code{y}.}
\item{max.iter}{An integer indicating the maximum number of iterations to perform in the optimization algorithm}
\item{tol}{A numeric value that is used for the termination rule (see details).}
\item{equal.residual.var}{logical, indicates whether the residual variances for the outcomes are assumed constant between studies or are heterogeneous}
}

\details{


	The last group in the covariates, designated by the final row in the list arguments, will be treated as the reference variable.

	The order of the list components should agree for all arguments.

	The optimization algorithm performs maximum likelihood estimation for the IPD LME model. Given the treatment indicator \code{z} the factor \code{x} and the study variable \code{study}, the fitted model is equivalent to \code{y~x*z+(z|study)} fitted with \code{lmer} with the option \code{REML=FALSE}. 

	The termination rule looks at the percentage change for the fixed and study random effects after the first iteration. Convergence is declared at the first iteration for which the maximum percentage change is less than \code{tol}.
}

\value{
  An object of class \code{\link{ipdlme-class}}, for which a number of methods
  are available.
}

\seealso{
  The \code{\link{ipdlme-class}}
}

\examples{
data(regress_age)

# Homogeneous residual variance
metafit <- ipdlme(n,Y,S2)

summary(metafit)

fixef(metafit)
confint(metafit)

ranef(metafit)

labs <- paste("Center",c(1,10,11,2:9))
plot(metafit, y=labs)

#Heterogeneous residual variance
metafit <- ipdlme(n,Y,S2,equal=FALSE)

summary(metafit)

fixef(metafit)
confint(metafit)

ranef(metafit)

labs <- paste("Center",c(1,10,11,2:9))
plot(metafit, y=labs)
}

\author{Stephanie A. Kovalchik <kovalchiksa@nih.gov>}

\keyword{models}

\keyword{methods}
