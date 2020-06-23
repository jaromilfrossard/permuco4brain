#' @importFrom graphics polygon
draw_circle <- function(x = 0,y=0,radius=1,border = "black",col = NA,nseg = 360,...){
  xx <- x + radius*cos( seq(0,2*pi, length.out=nseg) )
  xx <- c(xx,xx[1])
  yy <- y + radius*sin( seq(0,2*pi, length.out=nseg) )
  yy <- c(yy,yy[1])
  polygon(xx,yy,border = border,col = col)}
