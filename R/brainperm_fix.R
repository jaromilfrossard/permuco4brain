#' @importFrom stats rnorm qt pt
brainperm_fix <- function(formula, data, method, threshold, np, P, graph, effect, coding_sum, test,type,
                             aggr_FUN, multcomp, return_distribution,ncores,new_method,rnd_rotation){

  ##Method$
  if(is.null(method)){method = "freedman_lane"}
  if(!new_method){
    method <- match.arg(method,c("freedman_lane","manly"))}

  ###Aggr FUN
  if(is.null(aggr_FUN)){
    switch(test,
           "t"={
             fun_name <- "the sum of squares"
             aggr_FUN <- function(x)sum(x^2)},
           "fisher"={
             fun_name <- "the sum"
             aggr_FUN <- function(x)sum(x)})
  }else{
    fun_name <- "a user-defined function"
  }


  ###FUN Permut

  switch(paste(method,sep="_"),
         "freedman_lane" = {
           funP = function(...) {graph_freedman_lane(...)}
         },
         "manly" = {
           funP = function(...) {graph_manly(...)}
         },
         {
           funP=function(...){eval(parse(text=paste("graph_",method,"(...)",sep="",collpase="")))}
         })


  ###FUN multcomp
  switch(multcomp,
         "clustermass" = {
           funMultComp = function(distribution,threshold,aggr_FUN,graph,alternative){
             compute_clustermass_array(distribution = distribution,threshold = threshold,
                                       aggr_FUN = aggr_FUN,graph = graph, alternative = alternative)}},
         "troendle" = {funMultComp = function(distribution,threshold,aggr_FUN,graph,alternative){
           compute_troendle_array(distribution = distribution,graph = graph, alternative = alternative)
         }})

  #Formula transforamtion

  ff <-  formula
  ff[[2]] <- NULL

  #model frames


  mf <- model.frame(ff, data = data)
  if(coding_sum){mf <- permuco:::changeContrast(mf, contr = contr.sum)}
  mm <- model.matrix(attr(mf, "terms"), data = mf)

  # terms
  terms <- terms(ff, special = "Error", data = data)
  indError <- attr(terms, "specials")$Error


  ##response
  signal <- eval(formula[[2]], parent.frame(n=2))
  dim_y = dim(signal)
  dnames = dimnames(signal)
  dim(signal) = c(dim_y[1],dim_y[2]*dim_y[3])


  ####select test==============================
  switch(test,
         "fisher" = {
           col_ref <- attr(mm,"assign")
           colx <- 1:max(attr(mm,"assign"))
           names(colx) <- attr(terms,"term.labels")
         },
         "t" = {
           col_ref <- 1:length(attr(mm,"assign"))
           qr_mm <- qr(mm)
           colx <- qr_mm$pivot[1:qr_mm$rank]
           if(method != "huh_jhun"){
             colx <- colx[(attr(mm,"assign")[colx])!=0]}
           names(colx) <- colnames(mm)[colx]})
  if(is.null(effect)){effect = colx}else{colx = colx[colx==effect]}



  #permutation matrix====================================

  #check dim of P
  if (!is.null(P)) {
    check_P <- permuco:::check_P(
      P = P,
      method = method,
      test = test,
      n = dim_y[1],
      ncol_x2 = as.numeric(table(attr(mm,"assign")[attr(mm,"assign")!=0])),
      ncol_x = NCOL(mm)
    )
    if (!check_P) {
      np <- permuco:::np.Pmat(P)
      P <- NULL
      warnings("P argument is not valid and will be recomputed")
    }
  }

  #create permutation matrices ==============================
  if(is.null(P)){switch(method,
                        "huh_jhun" = {
                          switch (test,
                                  "t" = {P <- Pmat(np = np, n = dim_y[1] - NCOL(mm) + 1,type = type)},
                                  "fisher" = {{
                                    P <- lapply(as.numeric(table(col_ref))[-1],
                                                function(cx){Pmat(np = np, n = dim_y[1] - NCOL(mm) + cx,type = type)})}}
                          )},
                        {P <- Pmat(np = np, n = dim_y[1],type = type)})}


  if(sum(permuco:::np.Pmat(P) <= 3999)>0){
    warning("The number of permutations is below 4000, p-values might be unreliable.")
  }
  np <- permuco:::np.Pmat(P)

  #create rnd_rotation matrices==============================
  if(method=="huh_jhun" & is.null(rnd_rotation)){
    rnd_rotation <- matrix(rnorm(dim_y[1]^2),ncol=dim_y[1])
  }

  ###initialisze output
  multiple_comparison <- list()
  length(multiple_comparison) <- length(colx)
  names(multiple_comparison) <-names(colx)

  if(test == "t"){
    multiple_comparison_less <- multiple_comparison_greater <- list()
    length(multiple_comparison_less) <- length(multiple_comparison_greater) <- length(colx)
    names(multiple_comparison_less) <- names(multiple_comparison_greater) <-names(colx)
  }else{
    multiple_comparison_less <- multiple_comparison_greater <- NULL
  }

  ##adjust multiple threshold
  switch(test,
         "t" = {
           df <- permuco:::compute_degree_freedom_fix(test = test,mm = mm,assigni = colx)},
         "fisher" = {
           df <- permuco:::compute_degree_freedom_fix(test = test,mm = mm,assigni = attr(mm,"assign"))})


  if(is.null(threshold)){
    switch(test,
           "t" = {
             threshold <- qt(p = 0.975,df = df)},
           "fisher" = {
             threshold <- qf(p = 0.95, df1 = df[,1],df2 = df[,2])})
  }else if(length(threshold)==1){
    threshold <- rep(threshold,length(colx))
  }else if(length(threshold)>1){
    threshold <- as.numeric(matrix(threshold, nrow = length(colx)))
  }

  ##permute graph dimensions
  if(is.null(dnames[[3]])){
    warning("dimnames of the 3D matrix are missing. Make sure that the vertices of igraph match the 3rd dimension of the 3D matrix.")
  }else{
    perm = match(vertex_attr(graph,"name"),dnames[[3]])
    if( !isTRUE(all.equal(perm,1:length(perm)))){
      warning("Reordering the graph vertices to match the 3rd dimension of the 3D matrix.")
      if(sum(is.na(perm))>=1){
        stop("Names of the graph vertices and channels must match.")
      }
      graph = permute.vertices(graph,perm )
    }

  }



  args <- list(y = signal, mm = mm, P = P, rnd_rotation = rnd_rotation, test = test, ncores = ncores)

  cat("Computing Effect:\n")

  for(i in 1:length(colx)){
    cat(i," (",names(colx)[i], ") of ", length(colx), ". Start at ", as.character(Sys.time()),". ",sep = "")
    cat("\n")

    ###dfi
    switch(test,
           "t" = {dfi <- df[i]},
           "fisher" = {dfi <- df[i,]}
    )

    ###test info

    test_info = list(test = test, df = dfi, alternative = "two.sided", method = method, np = np,
                     nDV = dim_y[2]*dim_y[3], fun_name = fun_name,type = attr(args$P,"type"))

    ##compute distribution
    args$colx <- which(col_ref == colx[i])

  if(method == "huh_jhun"&test =="fisher"){args$P <- P[[i]]}

  distribution <- funP(args = args)

  pvalue <- apply(distribution,2,function(col){
    permuco:::compute_pvalue(distribution = col)})
  ##inarray
  pvalue = matrix(pvalue,nrow = dim_y[2],ncol = dim_y[3])
  dim(distribution) = c(np,dim_y[2],dim_y[3])
  #dim(signal) = dim_y

  multiple_comparison[[i]]=list()
  multiple_comparison[[i]]$uncorrected = list(statistic = t(distribution[1,,]),pvalue = pvalue,
                                              test_info = test_info)
  if(return_distribution){multiple_comparison[[i]]$uncorrected$distribution = distribution}

  multiple_comparison[[i]][[2]] =
    funMultComp(distribution = distribution,
                threshold = threshold[i], aggr_FUN = aggr_FUN, graph = graph, alternative = test_info$alternative)
  names(multiple_comparison[[i]])[2] = multcomp

  if(test=="t"){
    test_info$alternative <- "greater"
    pvalue_para <- pt(distribution[1,,],df = dfi,lower.tail = F)
    pvalue <- apply(as.matrix(distribution,nrow=np),2,function(col)permuco:::compute_pvalue(distribution = col,alternative = test_info$alternative))
    pvalue <- matrix(pvalue,nrow = dim_y[2],ncol = dim_y[3])

    multiple_comparison_greater[[i]]$uncorrected <- list(statistic = t(distribution[1,,]),pvalue = pvalue,
                                                         test_info = test_info)
    multiple_comparison_greater[[i]]$clustermass =
      funMultComp(distribution = distribution,
                  threshold = threshold[i], aggr_FUN = aggr_FUN, graph = graph, alternative = test_info$alternative)

    test_info$alternative <- "less"
    pvalue_para <- pt(distribution[1,,],df = dfi,lower.tail = T)
    pvalue <- apply(as.matrix(distribution,nrow=np),2,function(col)permuco:::compute_pvalue(distribution = col,alternative = test_info$alternative))
    multiple_comparison_less[[i]]$uncorrected <- list(statistic = t(distribution[1,,]),pvalue = pvalue,
                                                         test_info = test_info)
    multiple_comparison_less[[i]]$clustermass =
      funMultComp(distribution = distribution,
                  threshold = threshold[i], aggr_FUN = aggr_FUN, graph = graph, alternative = test_info$alternative)



  }



  }

  dim(signal) = dim_y
  dimnames(signal) = dnames
  attr(mf,"terms")=NULL

  out <- list()
  out$y <- signal
  out$model.matrix <- mm
  out$link <- link
  out$P <- P
  out$multiple_comparison <- multiple_comparison
  out$multiple_comparison_greater <- multiple_comparison_greater
  out$multiple_comparison_less <- multiple_comparison_less
  #out$table = table
  out$data<-mf
  out$method <- method
  out$multcomp <- multcomp
  out$threshold <- threshold
  out$effect <- effect
  out$graph <- graph
  class(out) <- "brainperm"
  return(out)

  return(out)




}
