#' Plot a statistical map of one effect.
#'
#' @description plot the significant test for 1 effect, all samples and all channels. Channels are in the y-axis and the samples in the x-axis. Non-significant cluster are shown in grey and the significant one in color from yellow to red as a function of the individual statistics.
#' @param x a brainperm object.
#' @param effect an integer indicating which effect to plot.
#' @param ... other argument pass to image().,
#' @importFrom tidyr pivot_wider
#' @importFrom graphics image box axis
#' @export
image.brainperm = function(x, effect = NULL,...){


  multcomp = x$multcomp[1]




  if(is.null(effect)){effect = x$effect[1]}

  dotargs = list(...)


  if(is.null(dotargs$alpha)){dotargs$alpha = 0.05}


  if(is.null(dotargs$alternative)){dotargs$alternative = "two.sided"}
  switch(dotargs$alternative,
         "two.sided" = {multiple_comparison = x$multiple_comparison},
         "greater" = {multiple_comparison = x$multiple_comparison_greater},
         "less" = {multiple_comparison = x$multiple_comparison_less})

  if(is.null(dotargs$xlab)){dotargs$xlab = "sample"}

  if(is.null(dotargs$ylab)){dotargs$ylab = "channel"}

  ##if(is.null(dotargs$add_border)){dotargs$add_border = FALSE}


  if(is.null(dotargs$main)){
    dotargs$main = names(multiple_comparison)[effect]
  }


  ## order from channel position
  order_channel = order(-vertex_attr(x$graph,"y"),vertex_attr(x$graph,"x"))
  enames = vertex_attr(x$graph,"name")[order_channel]

  ## non corrected fvalue
  #fvalue=multiple_comparison[[effect]]$uncorrected$statistic[order_channel,]

  data=multiple_comparison[[effect]][[multcomp]]$data

  #data = data[data$cluster_id!=0,]

  ## transform pvalue into matrix channel-sample

  pval_mat = pivot_wider(data[,c("channel", "sample", "pvalue")],names_from = "sample",values_from = "pvalue")
  enames_pv = as.character(pval_mat$channel)
  pval_mat = pval_mat[order(match(as.character(pval_mat$channel),vertex_attr(x$graph,"name"))),]
  pval_mat = as.matrix(pval_mat[,-1])

  ## transform fvalue into matrix channel-sample

  f_mat = pivot_wider(data[,c("channel", "sample", "statistic")],names_from = "sample",values_from ="statistic")
  f_mat = f_mat[order(match(as.character(f_mat$channel),vertex_attr(x$graph,"name"))),]
  f_mat = as.matrix(f_mat[,-1])


  ## transform fvalue into matrix channel-sample
  if(multcomp=="clustermass"){
    cl_mat = pivot_wider(data[,c("channel", "sample", "cluster_id")],names_from = "sample",values_from ="cluster_id")
    cl_mat = cl_mat[order(match(as.character(cl_mat$channel),vertex_attr(x$graph,"name"))),]
    cl_mat = as.matrix(cl_mat[,-1])}else{
      cl_mat = matrix(as.numeric(pval_mat<dotargs$alpha),nrow = nrow(pval_mat))
    }

  # if(dotargs$add_border){
  #   b_mat <- data
  #   egos <- ego(multiple_comparison[[effect]][[multcomp]]$graph, V(multiple_comparison[[effect]][[multcomp]]$graph))
  #   b_mat$border <- sapply(egos,function(egoi){
  #     ids = vertex_attr(multiple_comparison[[effect]][[multcomp]]$graph,name = "cluster_id",index = egoi)
  #     if(ids[1]==0){
  #       return(F)}else{
  #         return(sum(ids[1]!= ids[-1])!=0)}})
  #
  #
  #   b_mat <- pivot_wider(b_mat[,c("channel", "sample", "border")],names_from = "sample",values_from ="border")
  #   b_mat = b_mat[order(match(as.character(b_mat$channel),vertex_attr(x$graph,"name"))),]
  #   b_mat = as.matrix(b_mat[,-1])
  #   b_mat[!b_mat] = NA
  #   b_mat = t(apply(b_mat, 2, rev))
  #
  # }



  ## non significant sample in grey
  ns_mat = f_mat
  ns_mat[pval_mat<dotargs$alpha]=NA
  ns_mat[pval_mat>=dotargs$alpha]=1
  ns_mat[cl_mat==0]=NA

  ## plot only significant in color
  sign_mat = f_mat
  sign_mat[pval_mat>=dotargs$alpha]=NA

  xlim= c(1,nrow(f_mat))
  ylim= c(1,ncol(f_mat))


  #rev matrix and plot
  sign_mat = t(apply(sign_mat, 2, rev))
  ns_mat = t(apply(ns_mat, 2, rev))

  image(1:nrow(ns_mat),1:ncol(ns_mat), ns_mat,yaxt="n",col="grey", xlab=dotargs$xlab, main=dotargs$main,ylab=dotargs$ylab)
  # if(dotargs$add_border){image(1:nrow(b_mat),1:ncol(b_mat), b_mat,add=T,col= "black")}
  if(sum(!is.na(sign_mat))!=0){image(1:nrow(sign_mat),1:ncol(sign_mat), sign_mat,add=T)}
  box()
  axis(side=2,at = c(xlim[1]:xlim[2])[c(xlim[1]:xlim[2])%% 2 == 0],
       labels = rev(enames[c(xlim[1]:xlim[2])%% 2 == 0]))
  axis(side=4,at = c(xlim[1]:xlim[2])[c(xlim[1]:xlim[2])%% 2 != 0],
       labels = rev(enames[c(xlim[1]:xlim[2])%% 2 != 0]))
}






