#' @title Modify Default Parameters For Base Learner
#' @description Modify default parameters for methods provided in the package.
#' @param FUN Method.
#' @param ... Modified arguments.
#' @export modify_parameter
#' @return It returns a new function with modified parameters.
#' @examples
#' glmnet_lasso <- modify_parameter(glmnet1, alpha = 1)
#' glmnet_ridge <- modify_parameter(glmnet1, alpha = 0)

modify_parameter <- function(FUN, ...) {
  if (!is.function(FUN)) {
    stop("FUN must be a valid function")
  }
  
  .FUN <- FUN  # Keep the original function reference
  args <- list(...)
  invisible(lapply(seq_along(args), function(i) {
    formals(.FUN)[[names(args)[i]]] <<- args[[i]]
  }))
  
  return(.FUN)  # Return the modified function
}
