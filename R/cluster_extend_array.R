cluster_extend_array <- function(distribution, graph, threshold){
  if(length(threshold)==1L){
    return(cluster_extend_array_scaltau(distribution = distribution, graph = graph, threshold = threshold))
  }else if(length(threshold)>1L){
    return(cluster_extend_array_vectau(distribution = distribution, graph = graph, threshold = threshold))
  }
}


cluster_extend_array_vectau <- function(distribution, graph, threshold){

  #remove larger threshold
  th_max <- max(sum(threshold<max(distribution)),1)
  threshold <- threshold[seq_len(th_max)]


  df_n <- data.frame(name = vertex_attr(graph,"name"))
  df_n$id <-  1:length(df_n$name)

  graph <- set_vertex_attr(graph, name="statistic", value = as.numeric(aperm(distribution,perm = c(3,2,1))))

  cc_list <- list()
  gi <- graph
  for(thi in seq_along(threshold)){
    gi <-  delete_vertices(gi,V(gi)[vertex_attr(gi, "statistic")<threshold[thi]])
    cc_list[[thi]]  <- clusters(gi,mode ="weak")

  }


  df_list <- lapply(cc_list,function(cc){
    left_join(
      data.frame(
        name = names(cc$membership),
        extend =cc$csize[as.integer(cc$membership)]),
      df_n, by = "name")})

  v0 <- vector(mode =  "integer",length= length(df_n$name))

  lapply(df_list, function(dfi){
    v0[dfi$id] <- dfi$extend
    aperm(array(v0,dim(distribution)[c(3,2,1)]),perm = c(3,2,1))
  })


}

########################
cluster_extend_array_scaltau <- function(distribution, graph, threshold){


  df_n <- data.frame(name = vertex_attr(graph,"name"))

  graph <- set_vertex_attr(graph, name="statistic", value = as.numeric(aperm(distribution,perm = c(3,2,1))))
  gi <-  delete_vertices(graph,V(graph)[vertex_attr(graph, "statistic")<threshold])

  cc <- clusters(gi,mode ="weak")

  df_cc<-
      data.frame(
        name = names(cc$membership),
        extend =cc$csize[as.integer(cc$membership)])


  df_n <- left_join(df_n, df_cc, by = c("name" = "name"))
  df_n$extend[is.na(df_n$extend)] <- 0L

  aperm(array(df_n$extend,dim(distribution)[c(3,2,1)]),perm = c(3,2,1))

}


