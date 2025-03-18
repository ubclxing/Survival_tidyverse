max_finite <- function(xx) max(xx[is.finite(xx)])
min_finite <- function(xx) min(xx[is.finite(xx)])

# calculate residuals for binary predictions, based on observed binary values and predicted probabilities
resid_bin <- function(ymat, yhat, xmat = NULL,
                      type = c("deviance", "pearson", "raw"),
                      resid_std = FALSE) {
  
  # Check if standardized residuals are requested
  if (resid_std) {
    if (is.null(xmat)) stop("Standardized residuals require information about the hat matrix.")
    
    # Calculate the hat matrix
    QQ <- qr.Q(qr(xmat))  # Q matrix from QR decomposition
    hats <- rowSums(QQ^2)  # Diagonal of the hat matrix, X(X'X)^(-1)X'
  } else {
    hats <- 0
  }
  
  # Calculate raw residuals (for continuous outcome models)
  resi0 <- ymat - yhat
  
  # Calculate residuals based on the specified type
  if (type == "deviance") {
    resi <- ifelse(ymat == 1, sqrt(-2 * log(yhat)), -sqrt(-2 * (log(1 - yhat))))
  }
  if (type == "pearson") {
    resi <- ifelse(ymat == 1, exp(log(1 - yhat) / 2 - log(yhat) / 2), -exp(log(yhat) / 2 - log(1 - yhat) / 2))
  }
  if (type == "raw") {
    resi <- ymat - yhat
  }
  
  # Standardize residuals if needed
  if (resid_std) {
    resi <- resi / sqrt(1 - hats)
  }
  
  # Avoid infinity residuals by capping the values
  m1 <- max_finite(resi)
  m2 <- min_finite(resi)
  resi[resi > m1] <- m1 + 100
  resi[resi < m2] <- m2 - 100
  
  return(resi)
}