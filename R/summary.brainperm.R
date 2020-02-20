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
      out[[i]] <- brainperm_cluster_table(multiple_comparison[[i]],effect_name = names(multiple_comparison)[i])}
    else if(multcomp == "troendle"){
      out[[i]] <- brainperm_pseudocluster_table(multiple_comparison[[i]],effect_name = names(multiple_comparison)[i],
                                                    multcomp = "troendle",... = ...)}
  },
  "full" = {
    out[[i]] <- brainperm_full_table(multiple_comparison[[i]],effect_name = names(multiple_comparison)[i],
                                     multcomp = multcomp,... = ...)
  })}



  class(out) <- append("listof_multcomp_table",class(out))
  return(out)

}
