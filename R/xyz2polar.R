# copy paste from eeguana
xyz2polar <- function (name, x, y, z, scale = TRUE) {
  az <- atan2(y, x)
  el <- atan2(z, sqrt(x^2 + y^2))
  x <- (pi/2 - el) * cos(az)
  y <- (pi/2 - el) * sin(az)
  if (scale) {
    k <- pi/2
  }
  else {
    k <- 1
  }
  data.frame(name = name, x = x/k, y = y/k,stringsAsFactors = F)
}
