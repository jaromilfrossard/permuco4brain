#' Cluster Depth Tests
#'
#' @description Compute the Cluster- multiple comparisons procedure.
#' @param distribution An 3d array representing the null distribution of multiple signal. The first dimension is the permutations, the second the samples, the third is the channels.
#' @param threshold A numeric specifying the height of the threshold to compute the clusters.
#' @param graph A igraph object representing the adjacency of the channels (used for the outputs only).
#' @param alternative A character string indicating the alternative hypothesis. Either \code{"greater"}, \code{"less"} or \code{"two.sided"}.
#' @param border A character string indicating the method to handles the border. Either \code{"reverse"}, reversing the computation of the cluster-depth for clusters at the border of the time frame, or \code{"ignore"}.
#' @return A list containing an igraph and a data.frame, with the results for each sample, channel.
#' @references Frossard, J., & Renaud, O. (2021). The Cluster Depth Tests: Toward Point-Wise Strong Control of the Family-Wise Error Rate in Massively Univariate Tests with Application to M/EEG. arXiv preprint arXiv:2105.07514.
#' @family MCP
#' @export
compute_clusterdepth_array <- function(distribution, threshold, graph, alternative, border){

  border <-  match.arg(border, c("reverse", "ignore"))

  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  switch(alternative,
         "greater" = {
           threshold <- abs(threshold)
         }, "less" = {
           threshold <- -abs(threshold)
           distribution <- - distribution
         }, "two.sided" = {
           distribution <- abs(distribution)
           threshold <- abs(threshold)
         })



  ## compute the cluster depth and the depth distribution
  cd_distr <-
    lapply(seq_len(dim(distribution)[3]),function(chani) {

      cluster <- permuco:::get_cluster_matrix(distribution[,,chani],threshold = threshold)
      depth_head <-  permuco:::get_clusterdepth_head(cluster, border = border)
      depth_tail <-  permuco:::get_clusterdepth_tail(cluster, border = border)

      distr_head<- permuco:::depth_distribution(distribution[,,chani], head_mat = depth_head)
      distr_tail<- permuco:::depth_distribution(distribution[,,chani], tail_mat = depth_tail)


      list(depth_head = depth_head, depth_tail = depth_tail,distr_head = distr_head,distr_tail = distr_tail)
    })


  depth_head <- do.call("abind",list(lapply(cd_distr,function(cd)cd$depth_head),along=3))
  depth_tail <- do.call("abind",list(lapply(cd_distr,function(cd)cd$depth_tail),along=3))


  ### maximum over all channels
  mxh <- max(sapply(cd_distr,function(x)ncol(x$distr_head)))
  distr_head <-lapply(cd_distr,function(cd){
    matrix(c(cd$distr_head,rep(0,mxh*nrow(cd$distr_head)-length(cd$distr_head))),ncol=mxh)})
  distr_head <-do.call("pmax",distr_head)

  mxt <- max(sapply(cd_distr,function(x)ncol(x$distr_tail)))
  distr_tail <-lapply(cd_distr,function(cd){
    matrix(c(rep(0,mxt*nrow(cd$distr_tail)-length(cd$distr_tail)),cd$distr_tail),ncol=mxt)})
  distr_tail <- do.call("pmax",distr_tail)



  #### troendle for each clusters
    cluster <- permuco:::get_cluster_matrix(t(distribution[1,,]),threshold = threshold)


    stats<-
      lapply(seq_len(nrow(cluster)),function(chani){
        lapply(seq_len(max(cluster[chani,])),function(cli){
          distribution[1,cluster[chani,]==cli,chani]
          } )
        })

    ## pvalue and maximum head tail
    pval_list <-
      lapply(stats, function(stats_chan){
        lapply(stats_chan,function(stats_cli){
          pmax(
            compute_troendle(rbind(stats_cli,distr_head[,seq_along(stats_cli),drop=F]),alternative = alternative)$main[,2],
            compute_troendle(rbind(stats_cli,distr_tail[,seq_along(stats_cli)+ncol(distr_tail)-length(stats_cli),drop=F]),alternative = alternative)$main[,2])
          })
        })

    pval_list <-lapply(pval_list,function(x){
      do.call("c",x)
    })


    ## reshaping outputs



    pval_mat = cluster
    pval_mat[pval_mat==0]<-NA



    pval_mat<-
      lapply(seq_len(nrow(cluster)),function(chani){
        out <-  pval_mat[chani,]
        if(sum(!is.na(out))>0){
          out[which(!is.na(out))] <- pval_list[[chani]]

        }
        out
      })

    fgraph <- full_graph(graph, t = dim(distribution)[2])
    fgraph = set_vertex_attr(fgraph, name="statistic", value = as.numeric(t(distribution[1,,])))

    df_g <- data.frame(channel = vertex_attr(fgraph,name = "channel"),
                       sample =   vertex_attr(fgraph,name = "sample"),
                       statistic =   vertex_attr(fgraph,name = "statistic"),
                       pvalue = as.numeric(do.call("rbind",pval_mat)),
                       stringsAsFactors = F)


    df_g$pvalue[is.na(df_g$pvalue)] = 1
    return(list(graph = fgraph, data = df_g, threshold = threshold))





}


