\name{meta.data}
\alias{meta.data}

\title{Simulated study-level clustered survival data}

\description{

Study-level trials of MLMA dataset. 10 eligible randomized controlled
parallel-group trials, 8 providing
data at the patient-level (\code{ipd.data}) and 2 providing KM plots
by treatment arm. 


 This is the aggregate component of a simulated mixed data meta dataset which supposes that 2 of 10 eligible studies did not provide patient-level data.

Intended for illustrating use of mixed-level meta-analysis estimation
for a time-to-event outcome.

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

\examples{

data(meta.data)

xyplot(surv~time|factor(trt,lab=c("Control","Treated")),
                        xlim=c(0,6),xlab="Time",ylab="S(t)",
                        group=factor(group,lab=paste("Trial",9:10)),
                        meta.data,type="b",aspect=1,
                        auto.key=TRUE)

}
