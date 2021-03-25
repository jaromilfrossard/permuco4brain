#' Threshold-free cluster-enhancement
#'
#' @description Compute the TFCE statistics and p-value
#' @param distribution An 3d array representing the null distribution of multiple signal. The first dimension is the permutations, the second the samples, the third is the channels.
#' @param graph A igraph object representing the adjacency of the channels.
#' @param alternative a character string indicating the alternative hypothesis. Either \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#' @param E  a numeric indicating the parameter associated to the cluster extend. Suggestion \code{E = 0.5}.
#' @param H a numeric indicating the parameter associated to the cluster height. Suggestion for t statistic: \code{H = 1}, for F statistic: \code{H = 2}.
#' @param ndh a numeric indicating the number of step for the approximation of the integral. Suggestion \code{ndh = 500}.
#' @return graph a list containing an igraph and a data.frame, with the results for each sample, channel.
#' @importFrom future.apply future_lapply
#' @importFrom abind abind
#' @family MCP
#' @export
compute_tfce_array <- function(distribution, graph , alternative, E, H, ndh){
  alternative <- match.arg(alternative, c("greater","less","two.sided"))

  switch(alternative,
         "greater" = {

         }, "less" = {
           distribution <- -distribution
         }, "two.sided" = {
           distribution <- abs(distribution)
         })


  range_d <- range(abs(distribution))
  dhi <- seq(from = range_d[1], to =range_d[2], length.out = ndh)
  dh <- dhi[2]-dhi[1]

  graph <- full_graph(graph, t = dim(distribution)[2])

  #out <- array(0,dim= dim(distribution))

  HH <- dh*dhi^H



  cat("Computing the TFCE statistic\n")

  out <- future_lapply(seq_len(dim(distribution)[1]),function(permi){
    outi <- cluster_extend_array(distribution = distribution[permi,,,drop=F],graph = graph,threshold = dhi)
    outi <- mapply(function(exti,hhi){exti^E*hhi},exti= outi, hhi = HH[seq_along(outi)], SIMPLIFY = F)
    Reduce(`+`,outi)
  })

  out <- do.call("abind",list(out,along = 1))


  # for(permi in seq_len(dim(distribution)[1])){
  #   if(progress){
  #     cat("Permutation", permi, "out of", dim(distribution)[1],"        \r")}
  #
  #   outi <- cluster_extend_array(distribution = distribution[permi,,,drop=F],graph = graph,threshold = dhi)
  #   outi <- mapply(function(exti,hhi){exti^E*hhi},exti= outi, hhi = HH, SIMPLIFY = F)
  #   out[permi,,] <- Reduce(`+`,outi)
  #
  # }

  tfce_distri = apply(out,1,max)


  graph <- set_vertex_attr(graph, name = "statistic",
                           value = as.numeric(t(distribution[1, , ])))

  graph <- set_vertex_attr(graph, name = "tfce",
                           value = as.numeric(t(out[1,,])))

  graph <- set_vertex_attr(graph, name = "pvalue",
                           value = permuco:::compute_pvalue(tfce_distri,as.numeric(t(out[1,,]))))


  df <- data.frame(channel = vertex_attr(graph,name=c("channel")),
                   sample =as.numeric(vertex_attr(graph,name=c("sample"))),
                   statistic =as.numeric(vertex_attr(graph,name=c("statistic"))),
                   tfce =as.numeric(vertex_attr(graph,name=c("tfce"))),
                   pvalue =as.numeric(vertex_attr(graph,name=c("pvalue"))),
                   stringsAsFactors = F)




  return(list(graph = graph, data = df, cluster= NULL, threshold = NULL, distribution =  tfce_distri,E = E, H = H, ndh = ndh, dh = dh))

}
