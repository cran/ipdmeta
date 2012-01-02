\name{ipd.sep}
\alias{ipd.sep}

\title{IPD meta-analysis Subgroup Effect Power Estimator}

\description{
The function estimates the power of an IPD meta-analysis to detect a specified subgroup effect (covariate-treatment interaction) based on summary statistics. 
}

\usage{
ipd.sep(
	      effect,
               event0=NULL,
               event1=NULL,
               mean0=NULL,
               mean1=NULL,
               var0=NULL,
               var1=NULL,
               x0=NULL,
               x1=NULL,
               s20=NULL,
               s21=NULL,
               n0=NULL,
               n1=NULL,
	       data,
               alpha=.05
	       )
}

\arguments{

All vector entries are grouped by study.

  \item{effect}{scalar, subgroup effect under alternative hypothesis}
  \item{event0}{vector, for binary outcome, events in group 0}
  \item{event1}{vector, for binary outcome, events in group 1}
  \item{mean0}{vector, for continuous outcome, mean in group 0}
  \item{mean1}{vector, for continuous outcome, mean in group 1}
  \item{var0}{vector, for contunuous outcome, sample variances for responses in group 0}
  \item{var1}{vector, for contunuous outcome, sample variances for responses in group 1}
  \item{x0}{vector of subgroup covariate means for group 0}
  \item{x1}{vector of subgroup covariate means for group 1}
  \item{s20}{vector of covariate sample variances for control group}
  \item{s21}{vector of covariate sample variances for treatment group}
  \item{n0}{vector of number of subjects for group 0}
  \item{n1}{vector of number of subjects for group 1}
  \item{data}{data frame containing the objects specified in response or covariate arguments}
 \item{alpha}{scalar significance level of Wald test (two-sided)}

}

\details{

If a data frame is supplied, then the object indicated in each vector argument is looked for in \code{data}.

For a patient-level binary outcome, \code{mean0}, \code{mean1}, \code{var0} and \code{var1} should not be specified. Zero event counts will be corrected with a 0.5 factor. For a continuous response, \code{event0} and \code{event1} should not be specified.

For a covariate that is a mean proportion, such as proportion male, no sample variances need to be specified. If no values are given for the sample variances \code{s20} and \code{s21} it will be assumed that the covariate is a mean proportion and the sample variances will be determined from the proportions.

The SEP for the IPD meta-analysis is based on a generalized linear mixed model for the patient-level analysis. The model has intercept, treatment, covariate and interaction fixed effects and independent random effects for the baseline and treatment by study. Under this model, an estimator for the subgroup effect variance, that is, the variance for the estimate of the covariate-treatment interaction, for either an identity or logistic GLMM, can be obtained from the study sample statistics. This variance is then used to estimate the power of the IPD meta-analysis for a specified subgroup effect based on a two-sided Wald test.

}

\value{

A list with the following named components:

\tabular{ll}{
 \code{esimated.power} \tab The estmated IPD meta-analysis interactive effect power \cr
 \code{power.lower} \tab Lower bound for \code{level} CI \cr
 \code{power.upper} \tab Upper bound for \code{level} CI \cr
 \code{estimated.se} \tab Estimated standard error of IPD meta-analysis interaction effect \cr
 \code{se.lower} \tab Lower bound for \code{level} CI \cr
 \code{se.upper} \tab Upper bound for \code{level} CI \cr
 \code{sigma} \tab The mean of the study residual variance \cr
\code{sigma0} \tab Estimate of intercept random effect variance from simple RE meta-analysis with DSL estimator \cr
\code{sigma1} \tab Estimate of treatment random effect variance simple RE meta-analysis with DSL estimator \cr
\code{level} \tab confidence level for Wald test
 }
}

\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(poynard)

#AGE SEP FOR IPD META-ANALYSIS OF BETA-ANTAGONISTS TO PREVENT GI BLEEDING EVENTS

#ALTERNATIVE HYPOTHESIS FOR AGE-TREATMENT EFFECT
#WITH 10 YEARS CHANGE TO OR TREATMENT EFFECT exp(beta*10)
#EFFECT MODIFIER CHANGES TREATMENT EFFECT BY 30%

beta = log(1.3)/10

age.sep <- 

ipd.sep(
 effect=beta,
 event0=bleed0,
 event1=bleed1,
 n0=n0,
 n1=n1,
 x0=age0,
 x1=age1,
 s20=age.s20,
 s21=age.s21,
 data=poynard
)

age.sep

#GENDER SUBGROUP EFFECT; 30% OR CHANGE BY GENDER

beta <- log(1.3)

gender.sep <- 

ipd.sep(
 effect=beta,
 event0=bleed0,
 event1=bleed1,
 n0=n0,
 n1=n1,
 x0=male0,
 x1=male1,
 data=poynard
)

gender.sep


}
