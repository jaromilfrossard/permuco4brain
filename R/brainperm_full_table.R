brainperm_full_table = function(x,effect_name, multcomp,...){
  if(is.null(x)){ return(NULL)}
  test_info <- x$uncorrected$test_info

  tab <- x[[multcomp]]$data

  if(multcomp=="clustermass"){
  attr(tab,"threshold") <- x$clustermass$threshold
  attr(tab,"fun_name") <- test_info$fun_name}
  attr(tab,"effect_name") <- effect_name
  attr(tab,"multcomp") <- multcomp
  attr(tab,"nDV") <- test_info$nDV
  attr(tab,"method") <- test_info$method
  attr(tab,"test") <- test_info$test
  attr(tab,"alternative") <- test_info$alternative
  attr(tab,"type") = test_info$type
  attr(tab,"df") <- test_info$df
  attr(tab,"np") <- test_info$np
  attr(tab,"table_type") <- "cluster"
  tab




}
