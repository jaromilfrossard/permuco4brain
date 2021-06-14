#' Step-down version of the max T multiple comparison procedure
#'
#' @description Compute the step-down version of the max T multiple comparison procedure
#' @param distribution An 3d array representing the null distribution of multiple signal. The first dimension is the permutations, the second the samples, the third is the channels.
#' @param graph A igraph object representing the adjacency of the channels.
#' @param alternative a character string indicating the alternative hypothesis. Either \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#' @return graph a list containing an igraph and a data.frame, with the results for each sample, channel.
#' @family MCP
#' @export
compute_stepdownmaxT_array = function(distribution,alternative, graph){
  alternative <- match.arg(alternative, c("greater","less","two.sided"))



  sdmaxt <- permuco::compute_stepdownmaxT(matrix(distribution,nrow=dim(distribution)[1],
                                         ncol=prod(dim(distribution)[-1])),alternative = alternative)

  pvalue = t(matrix(sdmaxt$main[,2], ncol = dim(distribution)[3]))


  graph <- full_graph(graph, t = dim(distribution)[2])

  graph <- set_vertex_attr(graph, name = "statistic",
                           value = as.numeric(t(distribution[1, , ])))

  graph <- set_vertex_attr(graph, name = "pvalue",
                           value = as.numeric(pvalue))


  df <- data.frame(channel = vertex_attr(graph, name = c("channel")),
                   sample = as.numeric(vertex_attr(graph, name = c("sample"))),
                   statistic = as.numeric(vertex_attr(graph, name = c("statistic"))),
                   pvalue = as.numeric(vertex_attr(graph, name = c("pvalue"))),
                   stringsAsFactors = F)



  return(list(graph = graph, data = df, cluster = NULL, distribution = NULL,
              threshold = NULL))



}




# distribution_mat = matrix(distribution,nrow=dim(distribution)[1],
#                           ncol=prod(dim(distribution)[-1]))

# distribution_rank = apply(-distribution_mat,2,function(coli){
#   rank(coli,ties.method = "max")
# })
# rank_uncorr <- rank(distribution_rank[1,],ties.method = "min")
# order_test = integer(0)
# p_corrected <- numeric(length = 0)
#
# urank = sort(unique(rank_uncorr))[1]
# ntest = length(sort(unique(rank_uncorr)))
# pb = txtProgressBar(min = 1, max = ntest, style = 3)
# testi = 0
# for(urank in sort(unique(rank_uncorr))){
#   testi=testi+1
#   setTxtProgressBar(pb, testi)
#   which_test <- which(urank==rank_uncorr)
#   order_test <- c(order_test,which_test)
#   pvali <- distribution_rank[,which(urank<=rank_uncorr),drop=F]
#   distr_min <- apply(pvali,1,min)
#   p_corrected <- c(p_corrected,
#                    permuco:::compute_pvalue(distribution = distr_min,
#                                             stat = matrix(distribution_rank[,which_test],nrow=1),
#                                             alternative = "less"))
# }
# p_corrected = cummax(p_corrected)[order(order_test)]
# p_corrected = matrix(p_corrected,nrow= dim(distribution)[2],ncol= dim(distribution)[3])
# ### in graph


# graph = full_graph(graph, t = dim(distribution)[2])
#
# graph = set_vertex_attr(graph, name = "statistic",
#                         value = as.numeric(t(distribution[1, , ])))
#
# graph = set_vertex_attr(graph, name = "pvalue",
#                         value = as.numeric(t(p_corrected)))
#
# g = delete_vertices(graph, V(graph)[get.vertex.attribute(graph,
#                                                          "pvalue") >= alpha])
# cc = clusters(g, mode = "weak")
# cc$mass_statistic = rep(NA,length(cc$csize))
# cc$pvalue = rep(paste0("<",alpha),length(cc$csize))
#
#
#
#
# graph = set.vertex.attribute(graph, "cluster_id",
#                              value = NA)
#
# graph = set.vertex.attribute(graph, "mass_statistic",
#                              value = NA)
#
# graph = set.vertex.attribute(graph, "cluster_id", index = names(cc$membership),
#                              value = cc$membership)
#
# df = data.frame(electrode = get.vertex.attribute(graph, name = c("electrode")),
#                 sample = as.numeric(get.vertex.attribute(graph, name = c("sample"))),
#                 statistic = as.numeric(get.vertex.attribute(graph, name = c("statistic"))),
#                 pvalue = as.numeric(get.vertex.attribute(graph, name = c("pvalue"))),
#                 cluster_id = as.numeric(get.vertex.attribute(graph, name = c("cluster_id"))),
#                 mass_statistic = as.numeric(get.vertex.attribute(graph, name = c("mass_statistic"))))
#
# df$cluster_id[is.na(df$cluster_id)] = 0
# df$mass_statistic[is.na(df$mass_statistic)] = 0
