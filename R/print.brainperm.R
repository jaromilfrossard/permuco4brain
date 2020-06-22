#' Print method for \code{brainperm} object
#'
#' @description Display the results of the \code{brainperm} object.
#'
#' @param x a \code{brainperm} object.
#' @param ... others arguments pass to code{summary.brainperm}.
#' @seealso summary.brainperm
#' @export
print.brainperm = function(x,...){
  print(summary.brainperm(x,...))
}
