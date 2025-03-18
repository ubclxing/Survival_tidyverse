#' @title Mean Square Error Loss
#' @description Compute the mean square error loss in survival analysis.
#' @details The mean square error is calculated on log-scale, only considering uncensored observations.
#' @param pre A numeric vector of predicted survival time.
#' @param object A dataframe or matrix, the first column is observed survival time, the second column is survival status.
#' @export mse
#' @return The value of MSE.
#' @examples
#' set.seed(1)
#' # simulate predicted survival time
#' x1 <- rnorm(100) + 5
#' # simulate observed survival time
#' x2 <- rnorm(100) + 10
#' # simulate observed survival status
#' x3 <- rep(c(1, 0), 50)
#' # calculate MSE
#' mse(x1, cbind(x2, x3))

mse <- function(pre, object) {
  # Find indices of uncensored observations (status = 1)
  uncensored_indices <- which(object[, 2] == 1)
  
  # Extract event times and predicted event times for uncensored observations
  event_time <- object[uncensored_indices, 1]
  predicted_event_time <- pre[uncensored_indices]
  
  # Ensure the event times and predicted times have the same names
  names(event_time) <- names(predicted_event_time)
  
  # Compute the mean squared error on the log scale
  mse_value <- mean((log(event_time) - log(predicted_event_time))^2)
  
  return(mse_value)
}
