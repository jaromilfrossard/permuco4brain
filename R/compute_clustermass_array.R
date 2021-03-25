#' Cluster-mass test
#'
#' @description Compute the cluster-mass test with adjacency define by a graph.
#' @param distribution An 3d array representing the null distribution of multiple signal. The first dimension is the permutations, the second the samples, the third is the channels.
#' @param threshold The threshold used to compute the clusters.
#' @param aggr_FUN The function that aggregate the cluster into a scalar (cluster mass).
#' @param graph A igraph object representing the adjacency of the channels.
#' @param alternative a character string indicating the alternative hypothesis. Either \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#' @return graph A list containing a igraph object, a data frame containing the channels, time, statistic, cluster-mass and p-values,
#' @importFrom igraph set_vertex_attr delete_vertices clusters vertex_attr V set.vertex.attribute
#' @family MCP
#' @export
compute_clustermass_array = function(distribution, threshold, aggr_FUN, graph, alternative){
  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  switch(alternative,
         "greater" = {
    threshold <- abs(threshold)
    is_selected = function(distr,tau){distr>tau}
    extreme <- function(x) max(x, na.rm = T)
  }, "less" = {
    threshold <- -abs(threshold)
    is_selected = function(distr,tau){distr<tau}
    extreme <- function(x) max(x, na.rm = T)
  }, "two.sided" = {
    distribution <- abs(distribution)
    is_selected = function(distr,tau){distr>tau}
    extreme <- function(x) max(x, na.rm = T)
  })




  graph = full_graph(graph, t = dim(distribution)[2])
  ###null
  mass_distribution = apply(distribution,c(1),function(stat){
    gi = set_vertex_attr(graph, name="statistic", value = as.numeric(t(stat)))
    gi = delete_vertices(gi,V(gi)[!is_selected(vertex_attr(gi, "statistic"),threshold)])
    cc = clusters(gi,mode = "weak")
    if(length(cc$membership)==0){
      return(0)}else{
        return(extreme(sapply(1:max(1,max(cc$membership)),function(i){
          aggr_FUN(vertex_attr(gi,name = "statistic", index = names(cc$membership)[cc$membership == i]))
        })))}

  })

  ##observed
  graph = set_vertex_attr(graph, name="statistic", value = as.numeric(t(distribution[1,,])))
  g = delete_vertices(graph,V(graph)[!is_selected(vertex_attr(graph, "statistic"),threshold)])
  cc = clusters(g,mode ="weak")
  if(length(cc$membership)>0){
  clustermass = sapply(1:max(cc$membership), function(i) {
    aggr_FUN(vertex_attr(g,name = "statistic", index = names(cc$membership)[cc$membership == i]))
  })
  pvalue = sapply(clustermass, function(mi) permuco:::compute_pvalue(stat = mi,
                                                                        distribution = mass_distribution, alternative  = alternative))
  }else{
    clustermass = numeric()
    pvalue = numeric()
  }
  cc$clustermass = clustermass
  cc$pvalue = pvalue

  graph = set.vertex.attribute(graph, "pvalue",value= 1)
  graph = set.vertex.attribute(graph, "cluster_id",value=0)
  graph = set.vertex.attribute(graph, "clustermass",value=NA)


  graph = set.vertex.attribute(graph,"pvalue", index = names(cc$membership),value = cc$pvalue[cc$membership])
  graph = set.vertex.attribute(graph,"cluster_id", index = names(cc$membership),value = cc$membership)
  graph = set.vertex.attribute(graph,"clustermass", index = names(cc$membership),value = cc$clustermass[cc$membership])

  df = data.frame(channel = vertex_attr(graph,name=c("channel")),
                  sample =as.numeric(vertex_attr(graph,name=c("sample"))),
                  statistic = as.numeric(vertex_attr(graph,name=c("statistic"))),
                  pvalue = as.numeric(vertex_attr(graph,name=c("pvalue"))),
                  cluster_id = as.numeric(vertex_attr(graph,name=c("cluster_id"))),
                  clustermass = as.numeric(vertex_attr(graph,name=c("clustermass"))),
                  stringsAsFactors = F)

  df$pvalue[is.na(df$pvalue)]=1
  df$cluster_id[is.na(df$cluster_id)]=0
  df$clustermass[is.na(df$clustermass)]=0


  return(list(graph = graph, data = df,cluster = cc, distribution = mass_distribution,threshold = threshold))
}

