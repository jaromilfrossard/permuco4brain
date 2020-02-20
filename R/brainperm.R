#'Computer clustermass test based on multiple signal.
#'
#'@description Compute permutation test of on multiple signals spatially distributed like full-scalp EEG analysis.
#'
#'@param formula A formula object defining the design of the model. The left part should a 3 dimentional array. Its row are the observations (or design), its colomn are the samples (or time-point) and the third dimention are the nodes of the graph. Each row should correspond to the \code{data}.
#'@param data A dataframe containing the design.
#'@param graph An igraph object representing an undirected graph. It specifies the neighborgoods/spatial relationship between the nodes.
#'@param np A scalar indicating the number of permutations. Will be overwrite by \code{P} if specified.
#'@param test A character string to specify the name of the test. Default is \code{"fisher"}. \code{"t"} is available for the fixed effects model.
#'@param method A character string specifying the re-sampling method. See \code{permuco}.
#'@param threshold See \code{clusterlm} in \code{permuco}.
#'@param aggr_FUN A function used as mass function. It should aggregate the statistics of a cluster into one scalar. Default is the sum of squares fot t statistic and sum for F statistic.
#'@param multcomp The multiple comparison procedure only \code{"clustermass"} (default) and \code{"troendle"}is available.
#'@param effect An integer indicating the effect to test. Refer to the \code{assign} attribute of the \code{model.matrix} object. If \code{NULL} it will test all the effects.
#'@param ... further arguments
#' @details
#' The random effects model is only avaible with a F statistic.\cr
#'
#' Other arguments could be pass in \code{...} :\cr \cr
#' \code{P} : A matrix containing the permutation of class \code{matrix} or \code{Pmat}; which is used for the reproductibility of the results. The first column must be the identity. \code{P} overwrites \code{np} argument.\cr \cr
#' \code{return_distribution = FALSE} : return the permutation distribution of the statistics. Warnings : return one high dimentional matrices (number of test times number of permutation) for each test.\cr
#' \code{coding_sum} : a logical defining the coding of the design matrix to \code{contr.sum}: set by default to \code{TRUE} for ANOVA (when the argument \code{test} is \code{"fisher"} ) to tests main effects and is set to \code{FALSE} when \code{test} is \code{"t"}.  If \code{coding_sum} is set to \code{FALSE} the design matrix is computed with the coding defined in the dataframe and the tests of simple effets are possible with a coding of the dataframe set to \code{contr.treatment}. \cr
#' \code{ncores} : An integer specifiying the number of cores for parrallel computing. Default is \code{detectCores()-1}.
#'
#'@import permuco
#'@export
brainperm <- function(formula, data, graph, np = 5000, test = "fisher", aggr_FUN = NULL,
                        method = NULL, threshold = NULL, multcomp = "clustermass", effect = NULL,...){


  Terms <- terms(formula, special = "Error", data = data)
  indError <- attr(Terms, "specials")$Error
  dotargs = list(...)


  if (is.null(dotargs$return_distribution)) {
    dotargs$return_distribution = F
  }

  if (is.null(dotargs$new_method)) {
    dotargs$new_method = F
  }


  if(is.null(dotargs$ncores)){ncores = detectCores()-1}

  if (is.null(dotargs$coding_sum)) {
    switch(test, t = {
      dotargs$coding_sum = F
    }, fisher = {
      dotargs$coding_sum = T
    })
  }

  multcomp <- match.arg(multcomp, c("clustermass", "troendle"), several.ok = F)
  if (is.null(indError)) {
    result <- brainperm_fix(formula = formula, data = data, method = method, threshold = threshold, np = np, P = dotargs$P,
                               graph = graph, effect = effect, coding_sum = dotargs$coding_sum, test = test,
                               aggr_FUN = aggr_FUN, multcomp = multcomp, ncores = ncores,
                               return_distribution = dotargs$return_distribution,new_method = dotargs$new_method,
                               rnd_rotation = dotargs$rnd_rotation)
  }
  else if (!is.null(indError)) {
    if (test != "fisher") {
      warning("Random effects model only accept fisher statistics. Test statistic is set to fisher.")
      test = "fisher"
    }
    result <- brainperm_rnd(formula = formula, data = data, method = method, threshold = threshold, np = np, P = dotargs$P,
                            graph = graph, effect = effect, coding_sum = dotargs$coding_sum, test = test,
                            aggr_FUN = aggr_FUN, multcomp = multcomp, ncores = ncores,
                            return_distribution = dotargs$return_distribution,new_method = dotargs$new_method)
  }
  return(result)
}
