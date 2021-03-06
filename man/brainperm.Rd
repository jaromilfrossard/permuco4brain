% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brainperm.R
\name{brainperm}
\alias{brainperm}
\title{Multiple comparisons procedure for multiple signals}
\usage{
brainperm(
  formula,
  data,
  graph,
  np = 5000,
  method = NULL,
  type = "permutation",
  test = "fisher",
  aggr_FUN = NULL,
  threshold = NULL,
  multcomp = "clustermass",
  effect = NULL,
  ...
)
}
\arguments{
\item{formula}{A formula defining the design of the model. The left part should be a 3 dimensional array. Its rows are the observations (or design), its column are the samples (or time-point) and the third dimension are the space (or the nodes of the graph). Each row should correspond to the design in the \code{data} argument.}

\item{data}{A dataframe containing the design.}

\item{graph}{An \code{igraph} object. It specifies the neighborhoods/adjacency between the signals/nodes.}

\item{np}{A scalar indicating the number of permutations. It will be overwrite if \code{P} is manually specified.}

\item{method}{A character string to specify the re-sampling method. See \code{permuco} for details on permutation methods.}

\item{type}{A character string to specify the type of re-sampling transformation. Default is \code{"permutation"} and \code{"signflip"} is also available. Is overridden if \code{P} is specified. See help from \code{Pmat} in \code{permuco}.}

\item{test}{A character string to specify the name of the test. Default is \code{"fisher"}. \code{"t"} is available for the fixed effects model.}

\item{aggr_FUN}{A function used as mass function. It should aggregate the statistics of a cluster into one scalar. Default is the sum of squares for t statistic and sum for F statistic.}

\item{threshold}{See \code{clusterlm} in \code{permuco}.}

\item{multcomp}{The multiple comparison procedure only \code{"clustermass"} (default), \code{"tfce"}, \code{"troendle"}is available.}

\item{effect}{An integer indicating the effect to test. It refers to the \code{assign} attribute of the \code{model.matrix} object. The default (\code{effect = NULL}) compute all effects.}

\item{...}{further arguments}
}
\value{
a \code{brainperm} object.
}
\description{
Compute permutation test with with multiples comparison procedure (cluster-mass, TFCE or Troendle).
}
\details{
The random effects model is only available with a F statistic.\cr

Other arguments could be pass in \code{...} :\cr \cr
\code{P} : A matrix containing the permutation of class \code{matrix} or \code{Pmat}; which is used for the reproducibility of the results. The first column must be the identity. \code{P} overwrites \code{np} argument.\cr \cr
\code{return_distribution = FALSE} : return the permutation distribution of the statistics. Warnings : return one high dimensional matrices (number of test times number of permutation) for each test.\cr
\code{coding_sum = TRUE} : a logical defining the coding of the design matrix to \code{contr.sum}: set by default to \code{TRUE} for ANOVA (when the argument \code{test} is \code{"fisher"} ) to tests main effects and is set to \code{FALSE} when \code{test} is \code{"t"}.  If \code{coding_sum} is set to \code{FALSE} the design matrix is computed with the coding defined in the dataframe and the tests of simple effects are possible with a coding of the dataframe set to \code{contr.treatment}. \cr
\code{E = 0.5} : a numeric for extend parameter of the TFCE.
\code{H} : a numeric for height parameter of the TFCE. When \code{test = "t"}, the default is \code{H = 2} and when \code{test = "fisher"}, the default is \code{H = 1}.
\code{ndh = 500} : an integer defining the number of steps when estimating the integral of the TFCE.
}
\concept{main functions}
