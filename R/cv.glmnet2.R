#' @export cv_glmnet2
cv_glmnet2 <- function(xx, yy, foldid, alpha = seq(0, 10, by = 2) / 10, 
                       lambda = exp(seq(log(10^-8), log(5), length.out = 100)), ...) {
  # Initialize a list to store fits for different alpha values
  fits <- vector("list", length(alpha))
  names(fits) <- alpha
  
  # Loop through each alpha value and fit the model
  for (ii in seq_along(alpha)) {
    fits[[ii]] <- cv.glmnet(x = xx, y = yy, foldid = foldid, alpha = alpha[ii], lambda = lambda, ...)
  }
  
  # Select the best model based on cross-validation error (cvm)
  if (length(alpha) == 1) {
    idx <- 1
  } else {
    # Find the index of the model with the minimum cvm at lambda.1se
    idx <- which.min(sapply(fits, function(model) model$cvm[which(model$lambda == model$lambda.1se)]))
  }
  
  # Add the selected alpha to the result
  fits[[idx]]$alpha <- alpha[idx]
  
  # Return the best fit model
  return(fits[[idx]])
}
