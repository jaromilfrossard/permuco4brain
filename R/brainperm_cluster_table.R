brainperm_cluster_table = function(x,effect_name){

  if(is.null(x)){ return(NULL)}
  test_info <- x$uncorrected$test_info

  data <- x$clustermass$data

  tab <- unique(data[,c(4,5,6)])
  tab <- tab[tab$cluster_id!=0,,drop=F]


  if(nrow(tab)==0){
    tab <- data.frame()
    attr(tab,"nocluster") = T
  }else{

  range_out <- do.call("rbind",lapply(tab$cluster_id,function(id_i){
    dfi <- data[data$cluster_id ==id_i, ,drop=F]
    tabi <- table(dfi$channel)
    max_sample = max(tabi)
    maxchan <- names(tabi)[tabi==max_sample]
    if(length(maxchan)>1){
      maxchan = paste0(maxchan[1]," (",length(maxchan),")")
    }
    data.frame(`First sample` = min(dfi$sample),`Last sample`=max(dfi$sample),`N.chan.` = sum(tabi!=0),`Main chan.` = maxchan,
      `Main chan. length` = max_sample, `N. test` = nrow(dfi))}
    ))

  tab = cbind(tab,range_out)[, c(2,4,5,6,7,8,9,3,1)]
  row.names(tab)= NULL
  colnames(tab) <- c("Cluster id", "First sample", "Last sample", "N. chan.","Main chan.", "Main chan. length",
                     "N. test", "Clustermass", "P(>mass)")
  attr(tab,"nocluster") = F
  }


  attr(tab,"threshold") <- x$clustermass$threshold
  attr(tab,"fun_name") <- test_info$fun_name
  attr(tab,"effect_name") <- effect_name
  attr(tab,"multcomp") <- "clustermass"
  attr(tab,"nDV") <- test_info$nDV
  attr(tab,"method") <- test_info$method
  attr(tab,"test") <- test_info$test
  attr(tab,"alternative") <- test_info$alternative
  attr(tab,"type") = test_info$type
  attr(tab,"df") <- test_info$df
  attr(tab,"np") <- test_info$np
  attr(tab,"table_type") <- "cluster"

  class(tab) <- append("multcomp_table", class(tab))
  tab


  }




# brainperm_table = function(x){
#   ct = lapply(1:length(x),function(i){
#     effect = x[[i]]
#     if(is.null(x[[i]])){append_to_name = ": effect not tested"}else{append_to_name = ""}
#     tab= data.frame(size = effect[[2]]$cluster$csize,
#                     mass = effect[[2]]$cluster$mass_statistic,
#                     pvalue = effect[[2]]$cluster$pvalue)
#     attr(tab,"threshold") = effect[[2]]$threshold
#     attr(tab,"effect_name") = paste0(names(x)[i],append_to_name)
#     class(tab) = append("cluster_table",class(tab))
#     tab
#   })
#   class(ct) = append("listof_cluster_table",class(ct))
#   names(ct) = names(x)
#   ct
# }
