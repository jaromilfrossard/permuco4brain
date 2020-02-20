# Compute a adjacency graph through time
#
# @description Extend a graph of a scalp to t time point
# @param graph A igraph object representing the adjacency of channels in a scalp.
# @param t A numeric representing the number of samples (time points).
# @return An igraph object with t \times V(graph) vertices.
#' @importFrom igraph as_adjacency_matrix graph_from_adjacency_matrix vertex.attributes
#' @importFrom Matrix bdiag Diagonal Matrix
#' @importFrom dplyr left_join
full_graph = function(graph, t){
  attr_df = data.frame(vertex.attributes(graph),stringsAsFactors = F)
  full_graph = as_adjacency_matrix(graph, type = c("upper"))
  nv = ncol(full_graph)
  names_df = expand.grid(channel = rownames(full_graph), sample = 1:t,stringsAsFactors = F)
  if(!is.null(attr_df$name)){
    names_df <- left_join(names_df,attr_df, by = c("channel" = "name"))
  }else{warning("The names channels should be the name attributes of the graph.")}
  names = paste(names_df[,1],"_",names_df[,2],sep="")
  full_graph = lapply(1:t,function(i)full_graph)
  full_graph = bdiag(full_graph)

  full_graph = full_graph +rbind(cbind(Matrix(0,nrow = nv*(t-1),ncol = nv),Diagonal(nv*(t-1))),
                                 cbind(Matrix(0,ncol = nv*(t),nrow = nv)))
  rownames(full_graph) = names
  colnames(full_graph) = names
  full_graph = graph_from_adjacency_matrix(full_graph,mode="upper")
  for(i in 1:ncol(names_df)){
    full_graph = set.vertex.attribute(full_graph,name = colnames(names_df)[i], value = names_df[,i])
  }
  return(full_graph)}
