#' @title Get outputs of LMM fitting and chi-square test
#' @name getFixCoeff
#' @description Get fixed-effects coefficients and p-values from chi-square test.
#' @usage getFixCoeff(x, save=TRUE)
#' @param x \code{lmm2met} object from \code{fitLmm}
#' @param save logical; export to csv file, \code{defaut = TRUE}
#' @details The function summarizes outputs from \code{summary} and \code{drop1} function.
#' These include overall means of measured variables (Intercept), coefficients of fixed effects and
#' significance levels (Pr(Chi)) from chi-square tests.
#' @return a dataframe
#'
#' @author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#' @references https://cran.r-project.org/web/packages/lme4/index.html
#' @seealso \code{\link{lmer}}, \code{\link{drop1}}, \code{\link{summary}}
#' @examples
#' #fitMet = fitLmm(fix=c('Sex','Age','BMI','Stage','Location','Tissue'), 
#' #random='(1|Id)', data=adipose, start=10, end=14)
#' #getFixCoeff(fitMet)
#' @export
getFixCoeff <- function(x, save=TRUE) {
  cat('Exporting table of length', 1, '\n')
  tmp = fixeff.tab(x)
  value = tmp$table

  #Export result
  if(save){
    write.csv(value, file = paste0(getwd(),"/output.csv"))
    cat('Succcessfully saved to -> ', paste0(getwd(),"/output.csv"), '\n')
  }
  value
}
