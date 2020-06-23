#' Plot a graph of the channels for 1 sample/time-point
#'
#' @description Plot the graphs of the channels for several samples. Red indicates significant clusters, grey non-significant and white represent below the threshold.
#'
#' @param x brainperm object with a clustermass test.
#' @param effect an integer specifying which effect to plot.
#' @param samples a vetor of integers specifying the samples to plots.
#' @param ... other arguments including: \code{effect}: an integer specifying the effect to plot, \code{alpha}: argument to pass in par() and argument to pass in plot().
#'
#' @return a plot of a graph of the channels.
#'
#' @importFrom graphics par plot text
#' @importFrom igraph get.edgelist
#' @family plotting functions
#' @export
plot.brainperm <- function(x, effect = 1, samples,...){
  # save parameters
  par0 <- par()

  dotargs = list(...)
  dotargs_par <- dotargs[names(dotargs) %in% names(par())]


  if(is.null(dotargs$alpha)){dotargs$alpha = 0.05}

  if (is.null(dotargs_par$mar)) {
    dotargs_par$mar = c(0, 0, 2, 0)
  }
  if (is.null(dotargs_par$oma)) {
    dotargs_par$oma = c(4, 4, 4, 1)
  }
  par(dotargs_par)

  #Vertex data
  graph <- x$graph
  df <- data.frame(vertex_attr(graph),stringsAsFactors = F)
  df <- with(df,xyz2polar(name,x,y,z))

  #radius of channels
  if(length(df$x)!=1){
    rad <- min(dist(cbind(df$x,df$y)),na.rm = T)*(1/4)}
  else{
    rad <-0.4
  }


  #Edge data
  df_edge <- data.frame(get.edgelist(graph),stringsAsFactors = F)
  colnames(df_edge)= c("from","to")
  df_edge <- left_join(df_edge,df,by = c("from"="name") )
  colnames(df_edge)[3:4] <- c("x.from","y.from")
  df_edge <- left_join(df_edge,df,by = c("to"="name") )
  colnames(df_edge)[5:6] <- c("x.to","y.to")

  #plot dimension
  if(length(df$x)==1){
    xlim <- c(df$x-0.5,df$x+0.5)
  }else{
    xlim <- range(df$x)
  }
  if(length(df$y)==1){
    ylim <- c(df$y-0.5,df$y+0.5)
  }else{
    ylim <- range(df$y)
  }

  xlim <- c(xlim[1]-rad,xlim[2]+rad)
  ylim <- c(ylim[1]-rad,ylim[2]+rad)


  #mfrow args
  mfrow <- floor(sqrt(length(samples)))
  mfrow <- c(mfrow,ceiling(length(samples)/mfrow))

  if (is.null(dotargs_par$mfrow)) {
    dotargs_par$mfrow <- mfrow
  }


  title <- names(x$multiple_comparison)[effect]


  par(dotargs_par)



  for(sample_id in 1:length(samples)){
    main <- paste0("sample: ",samples[sample_id])

  plot(0,type="n",xlim = xlim,ylim = ylim, bty ="n",
       xaxt = "n",yaxt = "n",xlab = "",ylab = "",main = main,...=...)
  with(df_edge,segments(x.from,y.from,x.to,y.to))
  df_sampi <- left_join(df,subset(x$multiple_comparison[[effect]]$clustermass$data,
                            sample==(samples[sample_id])), by = c("name" = "channel"))

  u_cl <- unique(df_sampi$cluster_id)
  for(cli in u_cl){
    df_cli = df_sampi[df_sampi$cluster_id == cli,,drop=F]
    if(cli==0){
      col = "white"
    }else if(df_cli$pvalue[1]>=dotargs$alpha){
      col = "grey"
    }else if(df_cli$pvalue[1]<dotargs$alpha){
      col = "red"}

    for(i in 1:nrow(df_cli)){
      draw_circle(df_cli$x[i],df_cli$y[i],radius = rad,col =col)
    }
  }
  text(labels = df$name,x=df$x,y=df$y)
  }

  title(title, outer = T, cex = 2)


  par0 <- par0[!names(par0) %in% c("cin", "cra", "csi", "cxy", "din", "page")]
  par(par0)

}
