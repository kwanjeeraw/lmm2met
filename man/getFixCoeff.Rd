% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getFixCoeff.R
\name{getFixCoeff}
\alias{getFixCoeff}
\title{Get outputs of LMM fitting and chi-square test}
\usage{
getFixCoeff(x, save=TRUE)
}
\arguments{
\item{x}{\code{lmm2met} object from \code{fitLmm}}

\item{save}{logical; export to csv file, \code{defaut = TRUE}}
}
\value{
a dataframe
}
\description{
Get fixed-effects coefficients and p-values from chi-square test.
}
\details{
The function summarizes outputs from \code{summary} and \code{drop1} function.
These include overall means of measured variables (Intercept), coefficients of fixed effects and
significance levels (Pr(Chi)) from chi-square tests.
}
\examples{
#fitMet = fitLmm(fix=c('Sex','Age','BMI','Stage','Location','Tissue'), 
#random='(1|Id)', data=adipose, start=10, end=14)
#getFixCoeff(fitMet)
}
\references{
https://cran.r-project.org/web/packages/lme4/index.html
}
\seealso{
\code{\link{lmer}}, \code{\link{drop1}}, \code{\link{summary}}
}
\author{
Kwanjeera W \email{kwanjeera.wan@mahidol.ac.th}
}
