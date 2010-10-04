\name{mixed.biased}
\alias{mixed.biased}

\title{Simulated Mixed-Level Meta-Analysis Survival Data with
  Ecological Bias}

\description{

List of \code{ipd}, patient-level data, and \code{meta}, study-level
data. 5 trials each with 200 observations for each trial. The
Weibull-based model had a -.2 log(hazard) main effect for the centered
covariate \code{x} and -.4 interaction effect between the centered
covariate and treatment factor. There was no additional contribution
for the mean of the covariate for main or interaction effect. This
mimics an ecological bias in which there is a patient-level main and
interaction effect but no effects at the study-level.

Intended for illustrating use of MLMA with separation of study-level
and patient-level estimation.

}

\format{

For \code{ipd}:

\tabular{ll}{
\code{group} \tab Cluster factor, 1 to 10\cr
\code{trt} \tab Treatment indicator, 1 or 0\cr
\code{x} \tab Continuous covariate\cr
\code{x.bar} \tab Mean of x by \code{group} \cr
\code{time} \tab Survival time\cr
\code{event} \tab Indicator if event was observed, 1 or 0\cr
\code{trt.group} \tab Factor of treatment and group
}

For \code{meta}:

\tabular{ll}{
\code{time} \tab Time corresponding to survival estimate\cr
\code{surv} \tab Survival estimate from treatment.group KM curve\cr
\code{sigma2} \tab Corresponding variance of KM survival estimate\cr
\code{group} \tab Cluster factor, 1 to 10\cr
\code{x} \tab Treatment within study sample mean of covariate\cr
\code{x.bar} \tab Same as x \cr
\code{sub.group} \tab Factor of treatment and group
}

}

\seealso{
\code{\link{make.clustered.ipd}}, \code{\link{mlma}}
}

\author{S. Kovalchik \email{skoval@ucla.edu}}

