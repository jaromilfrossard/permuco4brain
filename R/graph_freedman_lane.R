graph_freedman_lane <- function(args){
  #select x
  ##test selection
  switch(args$test,
         "fisher"= {funT = function(qr_rdx, qr_mm, prdy){
           colSums(qr.fitted(qr_rdx, prdy)^2)/colSums(qr.resid(qr_mm, prdy)^2)* (NROW(prdy)-qr_mm$rank)/(qr_rdx$rank)}
         },
         "t" = {funT = function(qr_rdx, qr_mm, prdy){
           as.numeric(qr.coef(qr_rdx, prdy))/sqrt(colSums(qr.resid(qr_mm, prdy)^2)/sum(rdx^2)) * sqrt(NROW(args$y)-qr_mm$rank)}
         })

  ##effect selection
  select_x <- c(1:length(attr(args$mm,"assign"))) %in% args$colx

  ##data reduction
  qr_mm <- qr(args$mm)
  qr_d <- qr(args$mm[,!select_x, drop = F])
  rdx <- qr.resid(qr_d, args$mm[, select_x, drop = F])
  qr_rdx <- qr(rdx)
  rdy <- qr.resid(qr_d, args$y)

  #####permutation
  cl <- makeCluster(args$ncores)
  out = t(parApply(cl = cl,permuco:::as.matrix.Pmat(args$P),2,function(pi){
    funT(qr_rdx = qr_rdx, qr_mm = qr_mm,
         prdy = permuco::Pmat_product(x = rdy, P =pi,type= attr(args$P,"type")))}))
  stopCluster(cl)
  return(out)
}