\name{meta.data}
\alias{meta.data}

\title{Simulated study-level clustered survival data}

\description{

Aggregate of clustered survival data, summarized by a set of survival estimates from the KM plot of each treatment group within study. This is the aggregate component of a simulated mixed data meta dataset which supposes that 2 of 10 eligible studies did not provide patient-level data.

Intended for illustrating use of  mixed data survival meta-analysis estimation.

}

\format{

\tabular{ll}{
\code{time} \tab Time corresponding to survival estimate\cr
\code{surv} \tab Survival estimate from treatment.group KM curve\cr
\code{sigma2} \tab Corresponding variance of KM survival estimate\cr
\code{group} \tab Cluster factor, 1 to 10\cr
\code{x} \tab Treatment within study sample mean of covariate\cr
\code{sub.group} \tab Factor of treatment and group
}

}

\seealso{
\code{\link{survfit}}, \code{\link{ipd.data}}
}

\author{S. Kovalchik \email{skoval@ucla.edu}}

