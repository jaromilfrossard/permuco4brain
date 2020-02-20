#'Create igraph object from position of the channels
#'
#' Create an igraph object defining adjacency of the channels. Adjacency is defined when 2 channels have a (EUclidean) distance inferior to delta.
#'
#' @param data a data.frame containing the name and position of channels.
#' @param delta an double defining the maximal distance for adjacency of two channels.
#' @param name the colomn in the data containing the name of the channels.
#' @param x the colomn in the data containing the X position of the channels.
#' @param y the colomn in the data containing the Y position of the channels.
#' @param z the colomn in the data containing the Z position of the channels.
#' @importFrom stats dist
#' @export
position_to_graph <- function(data, delta = 4, name = "name", x = "x", y = "y", z = "z"){
  mc <- match.call()
  f <- formals()
  f <- f[!(names(f)%in%names(mc))]
  ## mc with default
  mc <- as.call(c(as.list(mc),f))

  layout <- list()
  data <- as.data.frame(data)

  for(i in c("name","x","y","z")){
    if(is.null(data[[mc[[i]]]])){stop(paste0("Error in argument ",i, ": the column ",as.character(mc[[i]])," is missing in the dataframe."))}
    layout[[i]] <- data[[mc[[i]]]]
  }
  layout <- as.data.frame(layout)

  distance_matrix <- dist(layout[, -1])

  adjacency_matrix <- as.matrix(distance_matrix) < delta
  diag(adjacency_matrix) <- FALSE

  dimnames(adjacency_matrix) <- list(layout$name, layout$name)


  graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")

  cl <- clusters(graph,mode = "weak")
  if(cl$no!=1){
    warning("This graph has more than 1 cluster. It is not valid for a clustermass test. Check the position of the channels or increase delta.")
  }


  graph <- delete_vertices(graph, V(graph)[!vertex_attr(graph, "name")%in%(layout[[1]])])


  graph <-set_vertex_attr(graph,"x", value = layout[match(vertex_attr(graph,"name"),layout[[1]]),][[2]])
  graph <-set_vertex_attr(graph,"y", value = layout[match(vertex_attr(graph,"name"),layout[[1]]),][[3]])
  graph <-set_vertex_attr(graph,"z", value = layout[match(vertex_attr(graph,"name"),layout[[1]]),][[4]])
  graph


}
