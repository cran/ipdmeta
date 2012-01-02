\name{regress_age}
\alias{Y}
\alias{S2}
\alias{n}
\alias{S2.subgroup}
\alias{regress_age}

\title{Aggregated data set of REGRESS trial outcomes for IPD meta-analysis of treatment covariate effect}

\description{

The data come from the Regression Growth Evaluation Statin Study (REGRESS) which was a double-blind, placebo-controlled multicenter study to assess the effects of 2 years of treatment with the 3-hydroxy-3-methylglutaryl coenzyme A reductase inhibitor pravastatin on progression and regression of coronary atherosclerosis in male patients with a serum cholesterol level between 4 and 8 mmol/L (155 and 310 mg/dL) by quantitative coronary arteriography. The primary end point of the trial was the average mean segment diameter per patient.

These data are used to demonstrate how a patient-level meta-analysis to estimate the effect modificaton of a categorical factor can be estimated from aggregate data. In this case centers of REGRESS are treated as trials and the categorical factor is baseline age (<50, 50-59, 60+).

Because the data are for demonstrative purposes only, observations have been permuted within treatment group and trial so that the summary statistics are not equivalent to the outcomes of the original study.
}

 
\format{


The data are formatted as needed for use with \code{\link{ipdlme}}. Note that the treatment group columns of each data frame must be named \code{trt} and \code{ctrl}.

\tabular{ll}{
\code{Y} \tab List of trial specific data frames containing the mean outcome (change is mean segment diameter) within each treatment-age subgroup. The treatment groups are by column and age groups by row. \cr
\code{n} \tab  List of trial specific data frames containing the sample sizes of each treatment-age subgroup.\cr
\code{S2} \tab Vector of trial sample variances in change in mean segment diameter.\cr
\code{S2.subgroup} \tab List of trial specific data frames containing the sample variance of each treatment-age subgroup.
 }
}

\references{
Jukema JW, Bruschke AVG, van Boven AJ, Reiber JHC, Bal ET, Zwinderman AH, et al. on behalf of the REGRESS Study Group Interuniversity Cardiology Institute Utrecht Netherlands. Effects of Lipid Lowering by Pravastatin on Progression and Regression of Coronary Artery Disease in Symptomatic Men With Normal to Moderately Elevated Serum Cholesterol Levels
\emph{Circulation.} 1995;91:2528-2540.
}


\author{S. Kovalchik \email{s.a.kovalchik@gmail.com}}

\examples{

data(regress_age)

n[1:3]

Y[1:3]

S2[1:3]

}