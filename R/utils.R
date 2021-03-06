fixeff.tab <- function(model){
  tb = data.frame()
  for(i in 1:length(model$completeMod)){
    rw = format.row(model$completeMod[[i]],model$testRes[[i]])
    tb = rbind(tb,rw)
  }
  row.names(tb) = names(model$completeMod)
  pchi = ncol(tb)-(nrow(model$testRes[[i]])-1)+1
  value = list(table=tb,pindex=pchi)
  return(value)
}

format.row <- function(x,y){
  coeff = t(data.frame(summary(x)$coefficients[,1]))
  pchi = t(data.frame(t(y)[4,-1]))
  colnames(pchi) = paste0('Pr(Chi).',rownames(y)[-1])
  cbind(coeff,pchi)
}

plotFix <- function(model){
  cat('Ploting table of length', 1, ' ...\n')
  tmp = fixeff.tab(model)
  dt = tmp$table
  md = dt[, tmp$pindex:ncol(dt)]
  md %>%
    mutate(
      ID = row.names(.)
    ) %>%
    mutate_if(is.numeric, function(y) {
      cell_spec(round(y,3), "html", color = ifelse(y < 0.05 & !is.na(y) , "red", "black"))
    }) %>%
    select(ID,colnames(md)) %>%
    kable("html", escape = F, align = "c") %>% kable_styling("striped", full_width = F)
}

plotCoeff <- function(model){
  cat('Ploting table of length', 1, ' ...\n')
  tmp = fixeff.tab(model)
  dt = tmp$table
  md = dt[, 1:tmp$pindex-1]
  md %>%
    mutate(
      ID = row.names(.)
    ) %>%
    mutate_if(is.numeric, function(y) {
      cell_spec(round(y,3), "html", color = "black")
    }) %>%
    select(ID,colnames(md)) %>%
    kable("html", escape = F, align = "c") %>% kable_styling("striped", full_width = F)
}

#rmarkdown::render("vignettes/lmm2met_vignette.Rmd", output_format="all")
#rmarkdown::render("vignettes/lmm2met_vignette.Rmd", output_format="word_document")
#rmarkdown::render("vignettes/lmm2met_vignette.Rmd", output_format="pdf_document")
