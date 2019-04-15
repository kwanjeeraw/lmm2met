#' @title Plot \code{lmm2met} object
#' @name plot
#' @description plot \code{lmm2met} object.
#' @usage plot.lmm2met(x, type, fix, ...)
#' @param x \code{lmm2met} object from \code{fitLmm}
#' @param type character; type of plot. See Details
#' @param fix character; name of a fixed effect, if \code{type ='fitted'}
#' @param ... other arguments of \code{plot} function
#' @details Types of plots include:
#'
#' \code{residual} - residual and a normal Q-Q plot
#'
#' \code{randeff} - random effects plot
#'
#' \code{fixeff} - heatmaps of coefficients and significance levels
#'
#' \code{fitted} - plots of measured variables before and after LMM fitting
#'
#' \code{chisq} - html table of p-values from chi-square test
#'
#' \code{coeff} - html table of fixed-effects coefficients
#'
#' @return
#' plot
#'
#' @author Kwanjeera W \email{kwanjeera.wan@@mahidol.ac.th}
#' @references https://cran.r-project.org/web/packages/lme4/index.html
#' @examples
#' #fitMet = fitLmm(fix=c('Sex','Age','BMI','Stage','Location','Tissue'), 
#' #random='(1|Id)', data=adipose, start=10, end=14)
#' #plot(fitMet, type='residual')
#' @import ggplot2
#' @import dplyr
#' @importFrom lme4 ranef
#' @importFrom gridExtra grid.arrange marrangeGrob
#' @importFrom lattice dotplot
#' @importFrom gplots heatmap.2
#' @importFrom knitr kable
#' @importFrom kableExtra kable_styling cell_spec
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom stats fitted qnorm quantile residuals
#' @importFrom utils write.csv
#' @export
plot.lmm2met <- function(x, type, fix=NULL, ...) {
  tmparg <- try(type <- match.arg(tolower(type), c("residual","randeff","fixeff","fitted","chisq","coeff"), several.ok = FALSE), silent = TRUE)
  if (class(tmparg) == "try-error") {
    stop("argument 'type' is not valid, choose one from the list: residual,randeff,fixeff,fitted,chisq,coeff")
  }
  #residual plot
  if(type=='residual'){
    #examine residual
    md = x$updateMod
    pl = list(); j = 1
    cat('Ploting graph of length', length(md),' ...\n')
    for(i in 1:length(md)){
      tmp = data.frame(fitt=fitted(md[[i]]),resi=residuals(md[[i]]))
      y = quantile(tmp$resi[!is.na(tmp$resi)], c(0.25, 0.75))
      x = qnorm(c(0.25, 0.75))
      slp = diff(y)/diff(x)
      intc = y[1L] - slp * x[1L]
      pl[[j]] = ggplot(tmp,aes(fitt,resi)) + ggtitle(names(md[i])) + geom_point() +
        stat_smooth(method="loess",colour="gray35",linetype="dashed") + geom_hline(yintercept=0, col="red") +
        xlab("Fitted values") + ylab("Residuals") + theme_bw() + theme(plot.margin = margin(20, 20, 20, 20, "pt"))
      pl[[j+1]] = ggplot(tmp, aes(sample=resi)) + ggtitle(names(md[i])) + stat_qq() +
        geom_abline(intercept=intc, slope=slp, colour="red") +
        xlab("Theoretical Quantiles") + ylab("Sample Quantiles") + theme_bw() + theme(plot.margin = margin(20, 20, 20, 20, "pt"))
      j=j+2
    }
    glist = lapply(pl, ggplotGrob)
    ggsave(paste0(getwd(),"/residual.pdf"), marrangeGrob(grobs=glist, nrow=2, ncol=2), width = 18, height=18)
    cat('Succcessfully saved to -> ', paste0(getwd(),"/residual.pdf"), '\n')
  }

  #random intercept & slope
  if(type=='randeff'){
    #examine random effects
    md = x$updateMod
    pl = list()
    cat('Ploting graph of length', length(md),' ...\n')
    for(i in 1:length(md)){
      key.settings <- list(title =names(md[i]), space="right", text = list(""))
      pl = c(pl,dotplot(ranef(md[[i]],condVar=T), scales = list(x = list(relation = "free")), key=key.settings,
                        par.settings=list(layout.widths=list(right.padding=5))))
    }
    ggsave(paste0(getwd(),"/randomEffect.pdf"), marrangeGrob(grobs=pl, nrow=2, ncol=2), width = 18, height=18)
    cat('Succcessfully saved to -> ', paste0(getwd(),"/randomEffect.pdf"), '\n')
  }

  #plot fixed effects heatmap
  if(type=='fixeff'){
    cat('Ploting graph of length', 1,' ...\n')
    #heatmap
    tmp = fixeff.tab(x)
    dt = tmp$table
    md = data.frame(dt[ , tmp$pindex:ncol(dt)]) #pchi
    colnames(md) = colnames(dt)[tmp$pindex:ncol(dt)]
    mc = data.frame(dt[ , 2:(tmp$pindex-1)]) #coeff
    colnames(mc) = colnames(dt)[2:(tmp$pindex-1)]
    h = 8 + (2*ceiling(nrow(md)/30))
    hei = 4 + ceiling(nrow(md)/30)
    breaks <- seq(0, 1, by=0.01)
    cols <- colorRampPalette(c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac"))(length(breaks)-1)
    pdf(file = paste0(getwd(),'/fixEffect.pdf'), width = 18, height = h)
    if(is.null(dim(md))){
      md = as.matrix(md)
      row.names(md)=row.names(dt)
      colnames(md)=colnames(dt)[-1:-ceiling(ncol(dt)/2)]
      heatmap.2(cbind(md,md),Rowv=NULL,Colv=NULL,
        col = cols, breaks = breaks,
        scale = "none", trace="none",  na.rm = TRUE, dendrogram="none",
        density.info="histogram", denscol="black", labCol = "",
        lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(1, hei), lwid=c(2, 5, 3),
        margins=c(10,0),key.par=list(mar=c(5,0,1,0)), key.title = "", key.xlab = colnames(md),
        keysize = 0.8, cexCol = 0.9, cexRow = 0.9)
      mc = as.matrix(mc)
      row.names(mc)=row.names(dt)
      colnames(mc)=colnames(dt)[ceiling(ncol(dt)/2)]
      heatmap.2(cbind(mc,mc),Rowv=NULL,Colv=NULL, col = colorRampPalette(c('blue','white', 'red'))(128),
                scale = "none", trace="none",  na.rm = TRUE, dendrogram="none",
                density.info="histogram", denscol="black", labCol = "",
                lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(1, hei), lwid=c(2, 5, 3),
                margins=c(10,0),key.par=list(mar=c(5,0,1,0)), key.title = "", key.xlab = colnames(mc),
                keysize = 0.8, cexCol = 0.9, cexRow = 0.9)
    }else{
      heatmap.2(as.matrix(md),Rowv=NULL,Colv=NULL,
        col = cols, breaks = breaks,
        scale = "none", trace="none",  na.rm = TRUE, dendrogram="none",
        density.info="histogram", denscol="black",
        lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(1, hei), lwid=c(2, 5, 3),
        margins=c(10,0),key.par=list(mar=c(5,0,1,0)), key.title = "", key.xlab = "Pr(Chi)",
        keysize = 1, cexCol = 1, cexRow = 1, srtCol=30)
      heatmap.2(as.matrix(mc),Rowv=NULL,Colv=NULL, col = colorRampPalette(c('blue','white', 'red'))(128),
                scale = "col", trace="none",  na.rm = TRUE, dendrogram="none",
                density.info="histogram", denscol="black",
                lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(1, hei), lwid=c(2, 5, 3),
                margins=c(10,0),key.par=list(mar=c(5,0,1,0)), key.title = "", key.xlab = "Coefficient",
                keysize = 1, cexCol = 1, cexRow = 1, srtCol=30)
    }
    dev.off()
    cat('Succcessfully saved to -> ', paste0(getwd(),"/fixEffect.pdf"), '\n')
  }

  #plot result
  if(type=='fitted'){
    if(missing(fix)){
      stop("Error: argument 'fix' is missing, with no default")
    }
    start = x$dataIndex[1]; end = x$dataIndex[2]; ind = c(start:end)
    data1 = x$rawDat
    data2 = data1[,1:end]
    df = rbind(data2, x$fittedDat)
    pl = list(); j = 1
    tmp = df[,1:start-1]
    tmp$Model = factor(rep(c('before','after'), each = nrow(data1)), levels = c('before','after'))
    cat('Ploting graph of length', length(ind),' ...\n')
    for(i in 1:length(ind)){
      coln = colnames(df)[ind[i]]
      tmp$val = df[,coln]
      pl[[j]] = ggplot(data = tmp, aes(x=as.factor(tmp[,fix]),y=val)) +
        geom_boxplot(aes(color=Model), outlier.size = 3.5, outlier.shape = 1) +
        geom_point(aes(color=Model), alpha=0.5, position = position_dodge(width=0.75)) + scale_color_brewer(palette="Set1") + labs(x=fix, y=coln) +
        theme_light() + theme(plot.margin = margin(20, 20, 20, 20, "pt"))

      pl[[j+1]] = ggplot(data = tmp, aes(x=as.factor(tmp[,fix]),y=val)) +
        geom_smooth(aes(group=Model,color=Model, linetype=Model, fill=Model), se = TRUE, method = 'loess', alpha=0.2) +
        geom_point(aes(color=Model), alpha=0.5, position = position_dodge(width=0.75)) + scale_color_brewer(palette="Set1") + labs(x=fix, y=coln) +
        theme_light() + theme(plot.margin = margin(20, 20, 20, 20, "pt"))
      j=j+2
    }
    glist = lapply(pl, ggplotGrob)
    ggsave(paste0(getwd(),"/fittedModel.pdf"), marrangeGrob(grobs=glist, nrow=2, ncol=2), width = 18, height=18)
    cat('Succcessfully saved to -> ', paste0(getwd(),"/fittedModel.pdf"), '\n')
  }

  #Pr.Chisq table
  if(type=='chisq'){
    pl=plotFix(x)
    return(pl)
  }

  #Coefficient table
  if(type=='coeff'){
    pl=plotCoeff(x)
    return(pl)
  }
}
