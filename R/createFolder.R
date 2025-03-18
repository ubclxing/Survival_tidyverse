#' @export createFolds
createFolds <- function(y, k = 10, list = TRUE, return_train = FALSE) {
  # Check if y is a Surv object, and extract time if so
  if (class(y)[1] == "Surv") {
    y <- y[, "time"]
  }
  
  # For numeric outcome, split into quantiles
  if (is.numeric(y)) {
    cuts <- floor(length(y) / k)
    
    # Ensure that cuts are within reasonable bounds
    if (cuts < 2) cuts <- 2
    if (cuts > 5) cuts <- 5
    
    breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
    y <- cut(y, breaks, include.lowest = TRUE)
  }
  
  # Create folds for factors
  if (k < length(y)) {
    y <- factor(as.character(y))
    num_in_class <- table(y)
    fold_vector <- vector(mode = "integer", length = length(y))
    
    # Loop through each class and assign folds
    for (i in seq_along(num_in_class)) {
      min_reps <- num_in_class[i] %/% k
      if (min_reps > 0) {
        spares <- num_in_class[i] %% k
        seq_vector <- rep(1:k, min_reps)
        
        if (spares > 0) {
          seq_vector <- c(seq_vector, sample(1:k, spares))
        }
        
        fold_vector[which(y == names(num_in_class)[i])] <- sample(seq_vector)
      } else {
        fold_vector[which(y == names(num_in_class)[i])] <- sample(1:k, size = num_in_class[i])
      }
    }
  } else {
    fold_vector <- seq_along(y)
  }
  
  # Split into list or return fold vector
  if (list) {
    out <- split(seq_along(y), fold_vector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq_along(out))), sep = "")
    
    # Return training set if requested
    if (return_train) {
      out <- lapply(out, function(data, y) y[-data], y = seq_along(y))
    }
  } else {
    out <- fold_vector
  }
  
  return(out)
}
