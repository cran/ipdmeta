\name{coxmeta.initialize}
\alias{coxmeta.initialize}

\title{Initialization for Mixed-Level Meta-Analysis Estimation}

\description{
Returns list of components for MLMA estimation.
}

\usage{
coxmeta.initialize(

   ipd.formula,
   meta.formula,
   random.formula,
   ipd.data,
   meta.data,
   ipd.groups,
   meta.groups,
   sigma2,
   study.group.interaction,
   beta.index

)
}

\arguments{

  \item{ipd.formula}{Formula for patient-level Cox model as would be
      supplied to \code{coxph}, \code{Surv(time,event)~trt*x}}
  \item{meta.formula}{Formula for study-level data \code{surv~log(time)+trt*x}}
  \item{ipd.data}{Patient-level data frame}
  \item{meta.data}{Study-level data frame containing survival estimates by study and treatment group}
 \item{ipd.groups}{Scalar number of clusters of patient-level type}
   \item{meta.groups}{Scalar number of clusters of study-level type}
  \item{sigma2}{Estimates variances for aggregate survival estimates of \code{meta.data}}
  \item{study.group.interaction}{Factor that is the study and treatment group interaction for \code{meta.data}}
  \item{beta.index}{List return object of \code{beta.indices}}

}


\value{

List of \code{model} and \code{params}. The components of \code{model}
are the output of \code{coxmeta.data.objects}, which consists of
\code{cox} (IPD) and \code{surv} (Study) components. These are lists
having formula, data, random design matrix and group information.

The \code{params} are return list of
\code{paramter.initialization}. This also has the features of the T
proposal density for the start of the algorithm (location, scale,
degrees of freedom) by each evidence type.

}

\details{

Dependent function \code{paramter.initialize} makes use of
\code{coxmcem.mixed.initialize} to get starting values for the frailty
location and variance parameters.

}

\author{S. Kovalchik \email{skoval@ucla.edu}}

\seealso{\code{\link{beta.indices}}}
