#' Print method for \code{brainperm} object
#'
#' @description Display the results of the \code{brainperm} object.
#'
#' @param x a \code{brainperm} object.
#' @param ... others arguments pass to code{summary.brainperm}.
#'
#' @return print a summary object
#'
#' @seealso summary.brainperm
#' @family summary functions
#' @export
print.brainperm = function(x,...){
  print(summary.brainperm(x,...))
}
