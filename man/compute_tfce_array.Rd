% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_tfce_array.R
\name{compute_tfce_array}
\alias{compute_tfce_array}
\title{Threshold-free cluster-enhancement (TFCE)}
\usage{
compute_tfce_array(distribution, graph, alternative, E, H, ndh)
}
\arguments{
\item{distribution}{An 3d array representing the null distribution of multiple signal. The first dimension is the permutations, the second the samples, the third is the channels.}

\item{graph}{A igraph object representing the adjacency of the channels.}

\item{alternative}{A character string indicating the alternative hypothesis. Either \code{"greater"}, \code{"less"} or \code{"two.sided"}.}

\item{E}{A numeric indicating the parameter associated to the cluster extend. Suggestion \code{E = 0.5}.}

\item{H}{A numeric indicating the parameter associated to the cluster height. Suggestion for t statistic: \code{H = 1}, for F statistic: \code{H = 2}.}

\item{ndh}{A numeric indicating the number of step for the approximation of the integral. Suggestion \code{ndh = 500}.}
}
\value{
A list containing an igraph and a data.frame, with the results for each sample, channel.
}
\description{
Compute the TFCE statistics and the corresponding p-values.
}
\references{
Smith, S. M., & Nichols, T. E. (2009). Threshold-free cluster enhancement: addressing problems of smoothing, threshold dependence and localisation in cluster inference. Neuroimage, 44(1), 83-98, https://doi.org/10.1016/j.neuroimage.2008.03.061.
}
\seealso{
Other MCP: 
\code{\link{compute_clusterdepth_array}()},
\code{\link{compute_clustermass_array}()},
\code{\link{compute_stepdownmaxT_array}()},
\code{\link{compute_troendle_array}()}
}
\concept{MCP}
