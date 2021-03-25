#' @importFrom future.apply future_apply
graph_Rd_kheradPajouh_renaud = function(args){
  #select x
  mm = args$mm
  assign = attr(mm,"assign")
  select_x = assign==args$i
  select_between = assign%in%args$link[2,]
  select_within = assign == (args$link[3,args$i])

  ##zmat
  z = permuco:::khatrirao(a = args$mm_id, b = mm[,select_within,drop=F])
  z = qr.resid(qr(mm),z)
  qr_z = qr(z)

  qr_d = qr(mm[,!select_x,drop=F])
  rdx = qr.resid(qr_d,mm[,select_x,drop=F])


  qr_rdx = qr(rdx)
  qr_rdz = qr(qr.resid(qr_d,z))
  ry = qr.resid(qr_d,args$y)



  #####permutation
  out = t(future_apply(permuco:::as.matrix.Pmat(args$P),2,function(pi){
    #colSums(qr.fitted(qr_rdx,ry[pi,,drop=F])^2)/colSums(qr.fitted(qr_rdz,ry[pi,,drop=F])^2)}))
    pry <- permuco::Pmat_product(x = ry, P =pi,type= attr(args$P,"type"))
    colSums(qr.fitted(qr_rdx,pry)^2)/colSums(qr.fitted(qr_rdz,pry)^2)}))


  out = out*(qr_rdz$rank/qr_rdx$rank)

  #

  return(out)}





# #'@importFrom parallel parApply makeCluster detectCores stopCluster
# graph_Rd_kheradPajouh_renaud = function(args){
#   #select x
#   mm = args$mm
#   assign = attr(mm,"assign")
#   select_x = assign==args$i
#   select_between = assign%in%args$link[2,]
#   select_within = assign == (args$link[3,args$i])
#
#   ##zmat
#   z = permuco:::khatrirao(a = args$mm_id, b = mm[,select_within,drop=F])
#   z = qr.resid(qr(mm),z)
#   qr_z = qr(z)
#
#   qr_d = qr(mm[,!select_x,drop=F])
#   rdx = qr.resid(qr_d,mm[,select_x,drop=F])
#
#
#   qr_rdx = qr(rdx)
#   qr_rdz = qr(qr.resid(qr_d,z))
#   ry = qr.resid(qr_d,args$y)
#
#
#
#   #####permutation
#   cl <- makeCluster(args$ncores)
#   out = t(parApply(cl = cl,permuco:::as.matrix.Pmat(args$P),2,function(pi){
#     #colSums(qr.fitted(qr_rdx,ry[pi,,drop=F])^2)/colSums(qr.fitted(qr_rdz,ry[pi,,drop=F])^2)}))
#     pry <- permuco::Pmat_product(x = ry, P =pi,type= attr(args$P,"type"))
#     colSums(qr.fitted(qr_rdx,pry)^2)/colSums(qr.fitted(qr_rdz,pry)^2)}))
#
#   stopCluster(cl)
#
#   out = out*(qr_rdz$rank/qr_rdx$rank)
#
#   #
#
#   return(out)}





