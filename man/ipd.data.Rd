\name{ipd.data}
\alias{ipd.data}

\title{Simulated patient-level clustered survival data}

\description{

Patient-level survival data with clustering by group. Meant to be an example of a partial survival meta-analysis dataset. The complete data set has 10 studies (group) with 200 observations per study. Covariates include treatment and a continuous covariate. There is only a treatment effect (beta = -0.5) on the log-hazard scale. There is between-group heterogeneity for baseline and treatment group, a moderate between-group heterogeneity for the covariate.

The simulated mixed data meta-analysis supposes that 2 of the 10 studies did not provide patient-level data. The evidence is instead in the form of survival estimates by treatment group within study, which are given in \code{\link{meta.data}}.

Intended for illustrating use of  mixed data survival meta-analysis estimation.

}

\format{

\tabular{ll}{
\code{group} \tab Cluster factor, 1 to 10\cr
\code{trt} \tab Treatment indicator, 1 or 0\cr
\code{x} \tab Continuous covariate\cr
\code{time} \tab Survival time\cr
\code{event} \tab Indicator if event was observed, 1 or 0\cr
\code{group.trt} \tab Factor of treatment and group
}
}

\seealso{
\code{\link{make.clustered.ipd}}, \code{\link{meta.data}}
}

\author{S. Kovalchik \email{skoval@ucla.edu}}

