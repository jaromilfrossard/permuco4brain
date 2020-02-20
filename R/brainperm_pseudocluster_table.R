#' @importFrom igraph delete.vertices
brainperm_pseudocluster_table = function(x,effect_name, multcomp,...){
  if(is.null(x)){ return(NULL)}
  test_info <- x$uncorrected$test_info

  dotargs = list(...)
  if(is.null( dotargs$alpha)){ dotargs$alpha = 0.05}

  data <- x[[multcomp]]$data
  data$cluster_id = 0

  graph <- x[[multcomp]]$graph
  graph = set.vertex.attribute(graph,"cluster_id",value = 0)

  # pvi = get.vertex.attribute(graph, "pvalue")
  # graph <- set.vertex.attribute(graph, "pvalue",index = c(20,21,200,2,150),value = c(0.01,0.01,0.01,0.01,0.01))
  #

  gi <- delete.vertices(graph,V(graph)[vertex_attr(graph, "pvalue")>=dotargs$alpha])

  cc = clusters(gi,mode ="weak")
  if(length(cc$membership!=0)){
    graph = set.vertex.attribute(graph,"cluster_id", index = names(cc$membership),value = cc$membership)

  }
  data$cluster_id = vertex_attr(graph, "cluster_id")

  tab <- unique(data[,5,drop = F])
  tab$pvalue = paste0("<",round(dotargs$alpha,4))
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

    tab = cbind(tab,range_out)[, c(1,3:8,2)]
    row.names(tab)= NULL
    colnames(tab) <- c("Cluster id", "First sample", "Last sample", "N. chan.","Main chan.", "Main chan. length",
                       "N. test", "P(>mass)")
    attr(tab,"nocluster") = F
  }


  attr(tab,"threshold") <- x$clustermass$threshold
  attr(tab,"fun_name") <- test_info$fun_name
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

  class(tab) <- append("multcomp_table", class(tab))
  tab


}



