#' Summary method for \code{brainperm} object
#'
#' @description List per effect of the results of the \code{brainperm} object.
#'
#' @param object a \code{brainperm} object.
#' @param multcomp a character string indicating the multiple comparison procedure to display. Only relevant if several multiple comparison procedures are computed in the \code{brainperm} object.
#' @param table_type a character string indicating the type of table. The default value, \code{table_type = "cluster"}, displays results by clusters. If \code{table_type = "full"}, the table shows the results by individual tests.
#' @param alternative a character string indicating the alternative hypothesis for t-test. The default value if \code{"two.sided"}. \code{"less"} and \code{"greater"} are also available.
#' @param ... other arguments pass to print.
#'
#' @return a list of tables containing the results of the test for each effect.
#'
#' @family summary functions
#' @export
summary.brainperm = function(object, multcomp = NULL, table_type = "cluster",
                                alternative = "two.sided",...){

  if(is.null(multcomp)){multcomp = object$multcomp[1]}
  if(!(multcomp%in%object$multcomp)){
    stop(paste0("Choose multcomp argument in:", object$multcomp))
  }

  table_type = match.arg(table_type,c("cluster","full"))

  ##switch alternative
  alternative = match.arg(alternative,c("two.sided","less","greater"))
  switch(alternative,
          "two.sided" = {multiple_comparison <- object$multiple_comparison},
         "less" = {multiple_comparison <- object$multiple_comparison_less},
         "greater" = {multiple_comparison <- object$multiple_comparison_greater})


  out = list()

  for(i in 1:length(multiple_comparison)){

  switch(table_type,
         "cluster" = {

    if(multcomp == "clustermass"){
      out[[i]] <- brainperm_cluster_table(multiple_comparison[[i]],effect_name = names(multiple_comparison)[i])
      }else if(multcomp == "troendle"){
      out[[i]] <- brainperm_pseudocluster_table(multiple_comparison[[i]],effect_name = names(multiple_comparison)[i],
                                                    multcomp = "troendle",... = ...)
      }else if(multcomp == "tfce"){
             out[[i]] <- brainperm_pseudocluster_table(multiple_comparison[[i]],effect_name = names(multiple_comparison)[i],
                                                       multcomp = "tfce",... = ...)
      }else if(multcomp == "clusterdepth"){
        out[[i]] <- brainperm_pseudocluster_table(multiple_comparison[[i]],effect_name = names(multiple_comparison)[i],
                                                  multcomp = "clusterdepth",... = ...)
      }else if(multcomp == "stepdownmaxT"){
        out[[i]] <- brainperm_pseudocluster_table(multiple_comparison[[i]],effect_name = names(multiple_comparison)[i],
                                                  multcomp = "stepdownmaxT",... = ...)}

  },
  "full" = {
    out[[i]] <- brainperm_full_table(multiple_comparison[[i]],effect_name = names(multiple_comparison)[i],
                                     multcomp = multcomp,... = ...)
  })}
  names(out) <- names(multiple_comparison)



  class(out) <- append("listof_multcomp_table",class(out))
  return(out)

}
