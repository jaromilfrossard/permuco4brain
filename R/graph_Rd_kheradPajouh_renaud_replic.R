#' @importFrom future.apply future_apply
graph_Rd_kheradPajouh_renaud_replic = function(args){
  #select x
  mm = args$mm
  assign = attr(mm,"assign")
  select_x = assign==args$i
  select_between = assign%in%args$link[2,]
  select_within = assign == (args$link[3,args$i])


  #####
  zid <- args$mm_id
  zid_0 <- as.integer(zid%*%t(contr.sum(ncol(zid)+1)))
  zid_0 = matrix( (zid_0==1L)+ (zid_0==ncol(zid)),ncol= ncol(zid)+1)
  factors = names(attr(mm,"contrasts"))
  eff = colnames(args$link)
  wfact_link <- sapply(strsplit(eff,":"),function(effi){sum(!effi%in%factors)==0})
  link_fact <- args$link[,wfact_link,drop=F]
  wfact_mm <- (attr(mm,"assign"))%in%c(0,link_fact[2,,drop=F]+link_fact[3,,drop=F])
  ZE = permuco:::khatrirao(zid_0,mm[,wfact_mm,drop=F])
  mm_gr=cbind(mm[,wfact_mm,drop=F],ZE)
  qr_gr = qr(mm_gr)
  #hii = diag(qr.fitted(qr_gr,diag(length(args$y))))
  duplic = duplicated.matrix(mm_gr)

  ######


  ##zmat
  z = permuco:::khatrirao(a = args$mm_id, b = mm[,select_within,drop=F])
  z = qr.resid(qr(mm),z)
  qr_z = qr(z)

  qr_d = qr(mm[,!select_x,drop=F])
  rdx = qr.resid(qr_d,mm[,select_x,drop=F])


  qr_rdx = qr(rdx)
  qr_rdz = qr(qr.resid(qr_d,z))
  ry = qr.resid(qr_d,args$y)



  ###
  #fpry_dup = fpry[!duplic,,drop=F]
  qr_d_dup = qr(mm[!duplic,!select_x,drop=F])
  rdx_dup = qr.resid(qr_d_dup,mm[!duplic,select_x,drop=F])
  qr_rdx_dup = qr(rdx_dup)


  z_dup = permuco:::khatrirao(a = args$mm_id[!duplic,,drop=F], b = mm[!duplic,select_within,drop=F])
  z_dup  = qr.resid(qr(mm[!duplic,,drop=F]),z_dup )
  qr_rdz_dup = qr(qr.resid(qr_d_dup,z_dup))
  ###



  #####permutation
  out = t(future_apply(permuco:::as.matrix.Pmat(args$P),2,function(pi){
    #colSums(qr.fitted(qr_rdx,ry[pi,,drop=F])^2)/colSums(qr.fitted(qr_rdz,ry[pi,,drop=F])^2)}))
    pry <- permuco::Pmat_product(x = ry, P =pi,type= attr(args$P,"type"))
    pry <-qr.fitted(qr_gr,pry)[!duplic,,drop=F]
    colSums(qr.fitted(qr_rdx_dup,pry)^2)/colSums(qr.fitted(qr_rdz_dup,pry)^2)}))


  out = out*(qr_rdz_dup$rank/qr_rdx_dup$rank)

  #

  return(out)}



# graph_Rd_kheradPajouh_renaud_replic = function(args){
#   #select x
#   mm = args$mm
#   assign = attr(mm,"assign")
#   select_x = assign==args$i
#   select_between = assign%in%args$link[2,]
#   select_within = assign == (args$link[3,args$i])
#
#
#   #####
#   zid <- args$mm_id
#   zid_0 <- as.integer(zid%*%t(contr.sum(ncol(zid)+1)))
#   zid_0 = matrix( (zid_0==1L)+ (zid_0==ncol(zid)),ncol= ncol(zid)+1)
#   factors = names(attr(mm,"contrasts"))
#   eff = colnames(args$link)
#   wfact_link <- sapply(strsplit(eff,":"),function(effi){sum(!effi%in%factors)==0})
#   link_fact <- args$link[,wfact_link,drop=F]
#   wfact_mm <- (attr(mm,"assign"))%in%c(0,link_fact[2,,drop=F]+link_fact[3,,drop=F])
#   ZE = permuco:::khatrirao(zid_0,mm[,wfact_mm,drop=F])
#   mm_gr=cbind(mm[,wfact_mm,drop=F],ZE)
#   qr_gr = qr(mm_gr)
#   #hii = diag(qr.fitted(qr_gr,diag(length(args$y))))
#   duplic = duplicated.matrix(mm_gr)
#
#   ######
#
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
#   ###
#   #fpry_dup = fpry[!duplic,,drop=F]
#   qr_d_dup = qr(mm[!duplic,!select_x,drop=F])
#   rdx_dup = qr.resid(qr_d_dup,mm[!duplic,select_x,drop=F])
#   qr_rdx_dup = qr(rdx_dup)
#
#
#   z_dup = permuco:::khatrirao(a = args$mm_id[!duplic,,drop=F], b = mm[!duplic,select_within,drop=F])
#   z_dup  = qr.resid(qr(mm[!duplic,,drop=F]),z_dup )
#   qr_rdz_dup = qr(qr.resid(qr_d_dup,z_dup))
#   ###
#
#
#
#   #####permutation
#   cl <- makeCluster(args$ncores)
#   out = t(parApply(cl = cl,permuco:::as.matrix.Pmat(args$P),2,function(pi){
#     #colSums(qr.fitted(qr_rdx,ry[pi,,drop=F])^2)/colSums(qr.fitted(qr_rdz,ry[pi,,drop=F])^2)}))
#     pry <- permuco::Pmat_product(x = ry, P =pi,type= attr(args$P,"type"))
#     pry <-qr.fitted(qr_gr,pry)[!duplic,,drop=F]
#     colSums(qr.fitted(qr_rdx_dup,pry)^2)/colSums(qr.fitted(qr_rdz_dup,pry)^2)}))
#
#   stopCluster(cl)
#
#   out = out*(qr_rdz_dup$rank/qr_rdx_dup$rank)
#
#   #
#
#   return(out)}
