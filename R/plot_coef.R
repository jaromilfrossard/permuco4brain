#' Plot of coefficients by effect, for multiple channels through sample.
#'
#' @description Function to provide plot for exploratory analysis. Shows the effects (beta) through samples for several channels.
#'
#' @param object a formula or brainperm object
#' @param ... other arguments including: \code{data} a dataframe of the design, \code{signal} a 3D array of the signal, \code{effect} a vector of the effect to plot, of \code{coding_sum}, a logical indicating if the effect should be computed using a contr.sum coding.
#' @importFrom graphics par
#' @importFrom stats ts.plot
#' @export
plot_coef = function(object,...){UseMethod("plot_coef")}

#' @export
plot_coef.formula = function(object, data, signal, effect = NULL, coding_sum = T,...){
  formula <- object


  terms <- terms(formula,special="Error",data=data)
  ind_error <- attr(terms, "specials")$Error
  error_term <- attr(terms, "variables")[[1 + ind_error]]

  formula_f <- update(formula, paste("~ .-",deparse(error_term, width.cutoff = 500L, backtick = TRUE)))
  e_term <- deparse(error_term[[2L]], width.cutoff = 500L,backtick = TRUE)

  formula_allfixed_design <- as.formula(paste(c("~",formula_f[[2]],"+",e_term),collapse=""))

  mf_design <- model.frame(formula = formula_allfixed_design, data = data)
  if(coding_sum){mf_design <- permuco:::changeContrast(mf_design, contr = contr.sum)}

  mf_f <- model.frame(formula = formula_f, data = mf_design)

  mm_f <- model.matrix(attr(mf_f, "terms"), data = mf_f)

  if(is.null(effect)){effect = attr(mm_f,"assign")[attr(mm_f,"assign")>0]}

  effect_names = c("intercept",attr(terms,"term.labels"))[effect+1]

  lisftofxmat = lapply(effect,function(assigni){
    mm_f[,attr(mm_f,"assign")==assigni,drop=F]
  })

  names(lisftofxmat) = effect_names

  plot_coef_listofx(signal, lisftofxmat)


}


plot_coef_listofx <- function(y, listofx){
  dimy = dim(y)
  ymat = matrix(y,nrow = dimy[1],ncol=prod(dimy[-1]))

  listof_coef = lapply(listofx,function(xi){
    qr_xi = qr(xi)
    coefi = matrix(qr.coef(qr_xi,ymat),ncol= prod(dimy[-1]))
    coefi = array(coefi,dim=c(nrow(coefi),dimy[-1]))
    coefi
  })

  par(mfrow=c(length(listof_coef),1))
  for(i in 1:length(listof_coef)){
    dimci = dim(listof_coef[[i]])
    coefi = aperm(listof_coef[[i]],c(2,1,3))
    coefi = matrix(coefi,nrow=dimci[2],ncol=prod(dimci[-2]))
    ts.plot(coefi,main=names(listofx)[i])

  }
  return(listof_coef)


  }
