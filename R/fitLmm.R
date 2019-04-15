#' @title Fit linear mixed-effects models (LMMs)
#' @name fitLmm
#' @description Fit a linear mixed-effects model (LMM) to each measured variable (e.g. metabolite) via \code{lmer} function of \pkg{lme4}.
#' @usage fitLmm(fix, random, data, start, end=NULL, auto=FALSE, major=NULL, pval=0.05, ...)
#' @param fix list or vector; list or vector of fixed effects
#' @param random list or vector; list or vector of random effects
#' @param data dataframe; a dataframe contains metadata and measured variables
#' @param start integer; index of the first measured variable
#' @param end integer; index of the last measured variable, If not given, all variables are fitted
#' @param auto logical; automatically include only significant fixed effects based on chi-square test, \code{defaut = FALSE}
#' @param major list or vector; list or vector of fixed effects always be included, if provided, \code{auto} is \code{= TRUE}
#' @param pval numeric; value of significance level for chi-square test, \code{defaut = 0.05}
#' @param ... other arguments of \code{lmer} function
#' @details The function processes a given data matrix (e.g. metabolomics data) by fitting a LMM to each measured variable via \code{lmer} function of \pkg{\link{lme4}}.
#' Significance level of the fixed-effect term is assessed by chi-square test as implemented in \code{drop1}.
#' If \code{auto = FALSE}, all fixed effects are included in models.
#' If \code{auto = TRUE}, only significant fixed effects with \code{Pr(Chi) < pval}, are included in models.
#' If \code{major} is provided, \code{major}-fixed effects and significant fixed effects with \code{Pr(Chi) < pval} are included in models.
#' Random effect must be provided to form LMMs. The random-effect term can be in several forms as given in Table 2 of Bates et al. (2015).
#' @return
#' \code{fitLmm} returns a LMM for each measured variable. Model information is returned as an object of class \code{lmm2met},
#' with the following components:
#' \code{completeMod} = list of LMMs containing all given fixed and random effects
#'
#' \code{updateMod} = list of LMMs containing fixed and random effects, if \code{auto = TRUE} or \code{major} is provided
#'
#' \code{testRes} = list of outputs from chi-square tests
#'
#' \code{fittedDat} = a dataframe of processed data
#'
#' \code{rawDat} = a dataframe of original data
#'
#' \code{dataIndex} = vector of indices of the first and last measured variable
#'
#' @author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#' @references https://cran.r-project.org/web/packages/lme4/index.html
#' @references Bates, D., et al., Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 2015. 67(1): p. 1-48.
#' @seealso \code{\link{lmer}}, \code{\link{drop1}}, \code{\link{update}}, \code{\link{summary}}
#' @examples
#' ##Ex 1. LMMs containing all given fixed and random effects
#' #fitMet = fitLmm(fix=c('Sex','Age','BMI','Stage','Location','Tissue'), 
#' #random='(1|Id)', data=adipose, start=10, end=14)
#' ##Ex 2. LMMs containing only significant fixed effects and random effects
#' #fitMet = fitLmm(fix=c('Sex','Age','BMI','Stage','Location','Tissue'), 
#' #random='(1|Id)', data=adipose, start=10, end=14, auto=TRUE)
#' ##Ex 3. LMMs containing only significant fixed effects, Tissue effect and random effects
#' #fitMet = fitLmm(fix=c('Sex','Age','BMI','Stage','Location','Tissue'), 
#' #random='(1|Id)', data=adipose, start=10, end=14, major='Tissue')
#' #structure(fitMet) #show lmm2met object
#' #summary(fitMet$updateMod[[1]]) #summarize outputs of fitting a LMM to variable X1
#' @importFrom stats as.formula drop1 predict update
#' @importFrom lme4 lmer
#' @export
fitLmm <- function(fix, random, data, start, end=NULL, auto=FALSE, major=NULL, pval=0.05, ...) {
  if(missing(fix) || missing(random) || missing(data) || missing(start)){
    stop("Error: argument is missing, with no default")
  }
  if(!is.character(fix) && !is.list(fix) && !is.vector(fix)){
    stop("Error: incorrect type of argument, 'fix' is not a character or list or vector")
  }

  #set default values
  end <- ifelse(is.null(end),ncol(data),end)
  completeMod <- list();  updateMod <- list();  testRes <- list();  fit <- data[,1:start-1]
  for(i in start:end){#loop through dependent variables
    val <- colnames(data)[i]
  cat("Fitting model to",val,"\n")
    #fit a linear mixed model (LMM)
    completeMod[[val]] <- lme4::lmer(as.formula(paste(val,"~", paste(do.call(c,list(fix,random)), collapse = "+"))), data=data, ...)
    #test all fixed effects
    testRes[[val]] <- drop1(completeMod[[val]], test = "Chisq", trace = F)
    if(auto){#select fixed effects to fit from Chisq test & pval cutoff
      fixlist <- fix[which(testRes[[val]]$`Pr(Chi)` < pval)-1]
      if(length(fixlist)>0){ dofix = fixlist }else{ dofix = 1 }
    }
    else if(!is.null(major)){#select fixed effects to fit from Chisq test & pval cutoff
      fixlist <- unique(c(fix[which(testRes[[val]]$`Pr(Chi)` < pval)-1], major))
      if(length(fixlist)>0){ dofix = fixlist }else{ dofix = 1 }
    }
    else{#use complete model
      dofix <- list()
    }
    if(length(dofix) > 0){#update with selected fixed effects
      updateMod[[val]] <- update(completeMod[[val]], as.formula(paste(val,"~", paste(do.call(c,list(dofix,random)), collapse = "+"))))
      fit <- cbind(fit, predict(updateMod[[val]]))
      names(fit)[ncol(fit)] <- val
    }
    else{#use complete model
      updateMod[[val]] <- completeMod[[val]]
      fit <- cbind(fit, predict(completeMod[[val]]))
      names(fit)[ncol(fit)] <- val
    }
  }
  value <- list(completeMod=completeMod, updateMod=updateMod, testRes=testRes, fittedDat=fit, rawDat=data, dataIndex=c(start,end))
  attr(value, "class") <- "lmm2met"
  value
}
