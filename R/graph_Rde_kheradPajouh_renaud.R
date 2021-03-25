#' @importFrom future.apply future_apply
graph_Rde_kheradPajouh_renaud = function(args){
  #select x
  mm = args$mm
  assign = attr(mm,"assign")
  select_x = assign==args$i
  select_between = assign%in%args$link[2,]
  select_within = assign == (args$link[3,args$i])

  within0 = unique(args$link[3, ])
  select_within_e = assign %in% (within0[!(within0 %in% args$link[3,
                                                                  args$i])])

  ##zmat
  z = permuco:::khatrirao(a = args$mm_id, b = mm[,select_within,drop=F])
  z = qr.resid(qr(mm),z)
  qr_z = qr(z)
  ee = permuco:::khatrirao(a = args$mm_id, b = mm[, select_within_e,
                                                  drop = F])

  ee = qr.resid(qr(mm), ee)
  ee = qr.resid(qr_z, ee)

  qr_d = qr(mm[,!select_x,drop=F])
  qr_de = qr(cbind(mm[, !select_x, drop = F], ee))

  rdex = qr.resid(qr_de, mm[, select_x, drop = F])
  qr_rdex = qr(rdex)

  qr_rdez = qr(qr.resid(qr_de, z))

  ry = qr.resid(qr_de,args$y)



  #####permutation
  out = t(future_apply(permuco:::as.matrix.Pmat(args$P),2,function(pi){
    #colSums(qr.fitted(qr_rdx,ry[pi,,drop=F])^2)/colSums(qr.fitted(qr_rdz,ry[pi,,drop=F])^2)}))
    pry <- permuco::Pmat_product(x = ry, P =pi,type= attr(args$P,"type"))
    colSums(qr.fitted(qr_rdex,pry)^2)/colSums(qr.fitted(qr_rdez,pry)^2)}))


  out = out*(qr_rdez$rank/qr_rdex$rank)

  #

  return(out)}




# #'@importFrom parallel parApply makeCluster detectCores stopCluster
# graph_Rde_kheradPajouh_renaud = function(args){
#   #select x
#   mm = args$mm
#   assign = attr(mm,"assign")
#   select_x = assign==args$i
#   select_between = assign%in%args$link[2,]
#   select_within = assign == (args$link[3,args$i])
#
#   within0 = unique(args$link[3, ])
#   select_within_e = assign %in% (within0[!(within0 %in% args$link[3,
#                                                                   args$i])])
#
#   ##zmat
#   z = permuco:::khatrirao(a = args$mm_id, b = mm[,select_within,drop=F])
#   z = qr.resid(qr(mm),z)
#   qr_z = qr(z)
#   ee = permuco:::khatrirao(a = args$mm_id, b = mm[, select_within_e,
#                                         drop = F])
#
#   ee = qr.resid(qr(mm), ee)
#   ee = qr.resid(qr_z, ee)
#
#   qr_d = qr(mm[,!select_x,drop=F])
#   qr_de = qr(cbind(mm[, !select_x, drop = F], ee))
#
#   rdex = qr.resid(qr_de, mm[, select_x, drop = F])
#   qr_rdex = qr(rdex)
#
#   qr_rdez = qr(qr.resid(qr_de, z))
#
#   ry = qr.resid(qr_de,args$y)
#
#
#
#   #####permutation
#   cl <- makeCluster(args$ncores)
#   out = t(parApply(cl = cl,permuco:::as.matrix.Pmat(args$P),2,function(pi){
#     #colSums(qr.fitted(qr_rdx,ry[pi,,drop=F])^2)/colSums(qr.fitted(qr_rdz,ry[pi,,drop=F])^2)}))
#     pry <- permuco::Pmat_product(x = ry, P =pi,type= attr(args$P,"type"))
#     colSums(qr.fitted(qr_rdex,pry)^2)/colSums(qr.fitted(qr_rdez,pry)^2)}))
#
#   stopCluster(cl)
#
#   out = out*(qr_rdez$rank/qr_rdex$rank)
#
#   #
#
#   return(out)}
