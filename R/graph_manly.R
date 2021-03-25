#' @importFrom future.apply future_apply
graph_manly <- function (args) {
  switch(args$test, fisher = {
    funT = function(qr_rdx, qr_mm, py) {
      colSums(qr.fitted(qr_rdx, py)^2)/colSums(qr.resid(qr_mm,
                                                        py)^2) * (NROW(py) - qr_mm$rank)/(qr_rdx$rank)
    }
  }, t = {
    funT = function(qr_rdx, qr_mm, py) {
      as.numeric(qr.coef(qr_rdx, py))/sqrt(colSums(qr.resid(qr_mm,
                                                            py)^2)/sum(rdx^2)) * sqrt(NROW(args$y) - qr_mm$rank)
    }
  })
  select_x <- c(1:length(attr(args$mm, "assign"))) %in%
    args$colx
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[, !select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  qr_rdx <- qr(rdx)
  qr_1 <- qr(rep(1, NROW(args$y)))
  r1y <- qr.resid(qr_1, args$y)
  h1y <- qr.fitted(qr_1, args$y)
  type = attr(args$P, "type")

  out = t(future_apply(permuco:::as.matrix.Pmat(args$P),2,function(pi){
    funT(qr_rdx = qr_rdx, qr_mm = qr_mm, py = permuco::Pmat_product(x = r1y,
                                                           P = pi, type = type) + h1y)}))
    return(out)
}
