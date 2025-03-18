#' @title Concordance Index (C-index)
#' @description Compute the C-index in survival analysis.
#' @details The most frequently used evaluation metric of survival models is the concordance index (C-index).
#' It is a measure of rank correlation between predicted risk scores. Only uncensored observations are considered.
#' @param pre A numeric vector of predicted survival times.
#' @param object A dataframe or matrix, where the first column is observed survival time and the second column is survival status.
#' @export c_index
#' @return The value of the C-index.
#' @examples
#' set.seed(1)
#' # Simulate predicted survival times
#' x1 <- rnorm(100)
#' # Simulate observed survival times
#' x2 <- rnorm(100)
#' # Simulate observed survival status
#' x3 <- rep(c(1, 0), 50)
#' # Calculate C-index
#' c_index(x1, cbind(x2, x3))

c_index <- function(pre, object) {
  total_pairs <- 0
  c <- 0
  test <- cbind(pre, object)
  
  # Sort the test data by the observed survival time (second column)
  test <- test[order(test[, 2]), ]
  
  for (i in 1:nrow(test)) {
    if (i == nrow(test)) break
    else if (test[i, 3] == 0) next  # Only consider uncensored observations
    else {
      for (j in (i + 1):nrow(test)) {
        total_pairs <- total_pairs + 1
        
        # If observed survival times are equal, consider the predicted survival times
        if (test[j, 2] == test[i, 2]) {
          if (test[j, 1] > test[i, 1] & test[j, 1] == test[i, 1]) {
            c <- c + 0.5
          }
        }
        
        # If the second observed survival time is greater than the first, and the predicted survival time is also greater, increment by 1
        if (test[j, 2] > test[i, 2]) {
          if (test[j, 1] > test[i, 1]) {
            c <- c + 1
          }
        }
        
        # If the second observed survival time is greater, and the predicted survival times are equal, increment by 0.5
        if (test[j, 2] > test[i, 2]) {
          if (test[j, 1] == test[i, 1]) {
            c <- c + 0.5
          }
        }
      }
    }
  }
  
  # Return the C-index as the ratio of concordant pairs
  return(c / total_pairs)
}
