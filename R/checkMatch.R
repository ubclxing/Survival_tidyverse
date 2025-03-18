#' @title Check Model-Family Match
#' @description Check if the specified model function matches the given family type.
#' @param family A character string indicating the family (e.g., "gaussian", "binomial", "survival").
#' @param FUN The model function to check.
#' @return A logical value indicating whether the function matches the family type.
#' @export check_match
#' @examples
#' # Example usage of the check_match function
#' check_match("gaussian", lm)
#' check_match("binomial", glm)
#' check_match("survival", surv)

check_match <- function(family, FUN) {
  # Define model functions for different families
  continuous <- c(lm1, glm1, glmnet1, rpart1, glmnet_lasso, glmnet_ridge)
  binary <- c(glm1, glmnet1, rpart1, glmnet_lasso, glmnet_ridge, lda1, qda1, KNN1, svm1)
  survival_models <- c(ela1, surv, lasso1)
  
  # Length of each model type
  continuous_length <- length(continuous)
  binary_length <- length(binary)
  survival_length <- length(survival_models)
  
  # Initialize the user-defined flag
  user_defined <- TRUE
  
  # Check if FUN matches any continuous model
  for (ii in seq_len(continuous_length)) {
    if (identical(FUN, continuous[[ii]])) {
      user_defined <- FALSE
      if (family == "gaussian") {
        return(TRUE)
      }
    }
  }
  
  # Check if FUN matches any binary model
  for (ii in seq_len(binary_length)) {
    if (identical(FUN, binary[[ii]])) {
      user_defined <- FALSE
      if (family == "binomial") {
        return(TRUE)
      }
    }
  }
  
  # Check if FUN matches any survival model
  for (ii in seq_len(survival_length)) {
    if (identical(FUN, survival_models[[ii]])) {
      user_defined <- FALSE
      if (family == "survival") {
        return(TRUE)
      }
    }
  }
  
  # Uncomment if you want to provide a warning when the method is undefined
  # if (user_defined) {
  #   warning("Please make sure that the method is consistent with the outcome type.")
  #   return(TRUE)
  # }
  
  # Return FALSE if no match was found
  return(FALSE)
}
