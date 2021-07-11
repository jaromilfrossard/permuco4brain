# #'Computer clustermass test based on multiple signal.
# #'
# #'@description compute the maris oostenveld permutation test of on multiple signal. Used for full scalp EEG analysis.
# #'
# #'@param formula the formula of the model
# #'@param data a dataframe containing the design
# #'@param signal a 3 dimentional array. The the row are the observations, colomn are the samples and the third dimention are the nodes of the graph.
# #'@param method a charcter string specifying the method
# #'@param P a Pmat object from \code{permuco}
# #'@param graph a igraph object of an undirected graph specifing the neighborgoods relationship between the nodes
# #'@param multcomp the multiple comparison procedure only \code{"maris_oostenveld"} is available.
# #'@param return_distribution If set to true return the null distribution by permutation.
# #'@param coding_sum logical. If \code{TRUE}, set the coding of the design to sum, if \code{FALSE}, take the coding define in the dataframe.
# #'@param threshold see \code{clusterlm}.
# #'@param np a scalar indicating the number of permutations. Will be overwrite by \code{P} if specified.
# #'@param aggr_FUN the function to aggregate individual statistics into cluster mass.
# #'@param effect a number indicating the effect to test. Refer to the \code{assign} attribute of the \code{model.matrix} object. If \code{NULL} it will test all the effects.
# #'@param ... further arguments
#'@importFrom stats update as.formula contr.sum model.frame contrasts<- model.matrix qf
#'@importFrom igraph permute.vertices
# #'@export
brainperm_rnd <- function(formula, data, method, threshold, np, P, graph, effect, coding_sum, test,type,
                               aggr_FUN, multcomp, return_distribution, new_method, E, H, ndh){


##Method$
  if(is.null(method)){method = "Rd_kheradPajouh_renaud"}
  if(!new_method){
    method <- match.arg(method,c("Rd_kheradPajouh_renaud","Rd_kheradPajouh_renaud_replic","Rde_kheradPajouh_renaud"))}

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
      "Rd_kheradPajouh_renaud" = {
        funP = function(...) {graph_Rd_kheradPajouh_renaud(...)}
        },
      "Rde_kheradPajouh_renaud" = {
        funP = function(...) {graph_Rde_kheradPajouh_renaud(...)}
      },
      {
        funP=function(...){eval(parse(text=paste("graph_",method,"(...)",sep="",collpase="")))}
        })

    ###FUN multcomp
  switch(multcomp,
         "clustermass" = {
           funMultComp = function(distribution,threshold,aggr_FUN,graph,alternative, E, H, ndh){
             compute_clustermass_array(distribution = distribution,threshold = threshold,
                                       alternative = alternative, aggr_FUN = aggr_FUN,graph = graph)}},
         "troendle" = {funMultComp = function(distribution,threshold,aggr_FUN,graph,alternative, E, H, ndh){
           compute_troendle_array(distribution = distribution,graph = graph, alternative = alternative)
         }},
         "stepdownmaxT" = {funMultComp = function(distribution,threshold,aggr_FUN,graph,alternative, E, H, ndh){
           compute_stepdownmaxT_array(distribution = distribution,graph = graph, alternative = alternative)
         }},
         "tfce" = {funMultComp = function(distribution,threshold,aggr_FUN,graph,alternative, E, H, ndh){
           compute_tfce_array(distribution = distribution,graph = graph, alternative = alternative,E = E, H = H, ndh = ndh)
         }},
         "clusterdepth" = {
           funMultComp = function(distribution,threshold,aggr_FUN,graph,alternative, E, H, ndh){
             compute_clusterdepth_array(distribution = distribution,threshold = threshold,
                                       alternative = alternative, graph = graph)
             }})





  #Formula transforamtion
  ff = formula
  ff[[2]]=NULL
  terms<-terms(ff,special="Error",data=data)
  ind_error <- attr(terms, "specials")$Error
  error_term <- attr(terms, "variables")[[1 + ind_error]]
  formula_f <- update(ff, paste("~ .-",deparse(error_term, width.cutoff = 500L, backtick = TRUE)))
  e_term <- deparse(error_term[[2L]], width.cutoff = 500L,backtick = TRUE)




    #random/fix formula

    formula_allfixed <- as.formula(paste(c("~",formula_f[[2]],"+",e_term),collapse=""))
    formula_allfixed_design <- as.formula(paste(c("~",formula_f[[2]],"+",e_term),collapse=""))

    formula_within <- formula(paste("~", e_term, collapse = ""))
    formula_within<- formula(paste("~",deparse(error_term[[2]][[3]]),collapse=""))
    formula_id <- formula(paste("~",deparse(error_term[[2]][[2]]),collapse = ""))


    #Model frame
    mf <- model.frame(formula = formula_allfixed, data = data)
    mf_design <- model.frame(formula = formula_allfixed_design, data = data)
    if(coding_sum){mf_design <- permuco:::changeContrast(mf_design, contr = contr.sum)}

    mf_f <- model.frame(formula = formula_f, data = mf_design)
    terms_f <- attr(mf_f,"terms")
    mf_id <- model.frame(formula = formula_id, data = as.data.frame(lapply(mf_design,function(col){
      col = as.factor(col)
      contrasts(col) = contr.sum
      col})))

    signal <- eval(formula[[2]], parent.frame(n=2))

    if(class(signal) != "array"){
      stop("Convert signal into a 3 dimensional array (Design x Sample x Channel)")
    }


    ##response
    mf <- eval(mf, parent.frame(n=1))
    dim_y = dim(signal)
    dnames = dimnames(signal)
    dim(signal) = c(dim_y[1],dim_y[2]*dim_y[3])

    ##link fixed random
    link = link(formula_f=formula_f,formula_within=formula_within)

    ###model .matrix
    mm_f <- model.matrix(attr(mf_f, "terms"), data = mf_f)
    mm_id <- model.matrix(attr(mf_id, "terms"), data = mf_id)[,-1,drop=F]
    name <- colnames(mm_f)

    ##check data


    permuco:::checkBalancedData(fixed_formula = formula_f, data = cbind(mf))

    #compute permutation
    if (is.null(P)) {P = Pmat(np = np, n = dim(signal)[1],type = type)}
    np = permuco:::np.Pmat(P)

    ##distribution
    args <- list(mm = mm_f, mm_id = mm_id, link = link, P = P, y = signal)

    if(is.null(effect)){effect = 1:max(attr(mm_f,"assign"))}


    multiple_comparison <- list()
    length(multiple_comparison) <- max(attr(mm_f,"assign"))
    names(multiple_comparison) <- attr(attr(mf_f, "terms"), "term.labels")


    ##adjust multiple threshold
    if(is.null(threshold)){
      df = permuco:::compute_degree_freedom_rnd(test = "fisher",mm = mm_f,assigni = attr(mm_f,"assign"),mm_id = mm_id,link = link)
      threshold = qf(p = 0.95, df1 = df[,1],df2 =df[,2])
    }else if(length(threshold)==1){threshold = rep(threshold,length(multiple_comparison))
    } else if(length(threshold)>1){
      threshold = as.numeric(matrix(threshold,nrow=length(multiple_comparison)))
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

    cat("Computing Effect:\n")

    for(i in effect){
      cat(i," (",attr(terms_f,"term.labels")[i], ") of ", length(effect), ". Start at ", as.character(Sys.time()),". ",sep = "")
      cat("\n")
      args$i = i

      test_info = list(test = test, df = df[i,], alternative = "two.sided", method = method, np = np,
                       nDV = dim_y[2]*dim_y[3], fun_name = fun_name,type = attr(args$P,"type"))

      ###initialisze output
      #argi<<- args
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
                    threshold = threshold[i], aggr_FUN = aggr_FUN, graph = graph,
                    alternative = test_info$alternative, E = E, H = H, ndh = ndh)
      names(multiple_comparison[[i]])[2] = multcomp
    }

    dim(signal) = dim_y
    dimnames(signal) = dnames
    attr(mf,"terms")=NULL


    #table = brainperm_table(multiple_comparison)


    out=list()
    out$y = signal
    out$model.matrix = mm_f
    out$model.matrix_id = mm_id = mm_id
    out$link = link
    out$P = P
    out$multiple_comparison = multiple_comparison
    #out$table = table
    out$data=mf
    out$method = method
    out$multcomp = multcomp
    out$test <- test
    out$threshold = threshold
    out$effect = effect
    out$graph = graph
    class(out) <- "brainperm"
    return(out)

  }
