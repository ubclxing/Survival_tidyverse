#' @title Area Under the Curve (AUC)
#' @description The AUC function calculates the numeric value of the area under the ROC curve (AUC) using the trapezoidal rule and optionally plots the ROC curve.
#' @details The ROC curve is created by plotting the true positive rate (TPR) against the false positive rate (FPR) at different threshold settings.
#' By default, the total area under the curve is computed. A truncated AUC statistic can be specified with the `cutoff` argument.
#' The `cutoff` argument specifies the bounds of FPR. Common cutoff values are 1 (i.e., no truncation) or 0.2 (i.e., specificity > 0.8).
#' @param prob A numeric vector of predicted probabilities.
#' @param outcome A numeric vector of observed binary outcomes.
#' @param cutoff A number between 0 and 1 specifying where the threshold of the ROC curve should be truncated. The default value is 1 (no truncation).
#' @param roc_plot Logical. If `TRUE`, the ROC curve will be plotted.
#' @export AUC
#' @return The value of the area under the curve (AUC).
#' @importFrom graphics abline
#' @examples
#' set.seed(1)
#' # Simulate predictors
#' x1 <- rnorm(200)
#' x2 <- rnorm(200)
#' # Simulate outcome
#' pr <- 1 / (1 + exp(-(3 * x1 + 2 * x2 + 1)))
#' y <- rbinom(200, 1, pr)
#' df <- data.frame(y = y, x1 = x1, x2 = x2)
#' # Fit logistic regression model on the first 100 observations
#' lg_model <- glm(y ~ x1 + x2, data = df[1:100, ], family = "binomial")
#' # Predict outcome for the last 100 observations
#' prob <- predict(lg_model, df[101:200, c("x1", "x2")], type = "response")
#' # Calculate AUC and plot the ROC Curve
#' AUC(prob, y[101:200], roc_plot = TRUE)
#' # Calculate AUC and plot the ROC Curve with cutoff
#' AUC(prob, y[101:200], cutoff = 0.2, roc_plot = TRUE)


library(ggplot2)

# AUC function to calculate the Area Under the Curve
AUC <- function(prob, outcome, cutoff = 1, roc_plot = FALSE) {
  nn <- length(outcome)
  
  # Check if the length of prob and outcome are the same
  if (length(prob) != nn) {
    stop("prob and binary should be the same length")
  }
  
  # Check if prob values are within [0, 1]
  if ((max(prob) > 1) | (min(prob) < 0)) {
    stop("prob values should be in [0,1]")
  }
  
  # Remove missing values from both prob and outcome
  data <- tibble(prob = prob, outcome = outcome) %>%
    filter(!is.na(prob), !is.na(outcome))
  
  # Sort by probability in descending order
  data <- data %>%
    arrange(desc(prob))
  
  # Calculate AUC components
  cp <- sum(data$outcome)  # Condition positive
  cn <- nn - cp            # Condition negative
  tp <- cumsum(data$outcome)  # True positive
  fp <- (1:nn) - tp            # False positive
  tpr <- tp / cp                # True positive rate (sensitivity)
  fpr <- fp / cn                # False positive rate (1 - specificity)
  
  # If ROC plot is requested, create the plot
  if (roc_plot) {
    ggplot(data, aes(x = fpr, y = tpr)) +
      geom_line() +
      geom_vline(xintercept = cutoff) +
      labs(x = "False Positive Rate", y = "True Positive Rate")
  }
  
  # Calculate AUC using the trapezoid method
  auc_value <- trapezoid(fpr, tpr, cutoff)
  
  return(auc_value)
}

# Helper function to calculate the AUC using the trapezoid method
trapezoid <- function(fpr, tpr, cc) {
  ord_fpr <- order(fpr)
  fpr <- fpr[ord_fpr]
  tpr <- tpr[ord_fpr]
  
  # If FPR is less than the cutoff, use the whole vector
  if (max(fpr) < cc) {
    fpr_cut <- fpr
    tpr_cut <- tpr
  } else if (min(fpr) > cc) {
    stop("threshold smaller than smallest allowed value")
  } else {
    # Identify the point where FPR is equal to the cutoff
    fpr_lcc <- sum(fpr < cc)
    fpr1 <- fpr[fpr_lcc]
    fpr2 <- fpr[fpr_lcc + 1]
    tpr1 <- tpr[fpr_lcc]
    tpr2 <- tpr[fpr_lcc + 1]
    
    # Interpolate to find TPR at the cutoff point
    tpr_cc <- (tpr1 * (cc - fpr1) + tpr2 * (fpr2 - cc)) / (fpr2 - fpr1)
    
    # Update FPR and TPR vectors to include the cutoff point
    fpr_cut <- c(fpr[1:fpr_lcc], cc)
    tpr_cut <- c(tpr[1:fpr_lcc], tpr_cc)
  }
  
  # Calculate the AUC using the trapezoidal rule
  return(sum(diff(fpr_cut) * (tpr_cut[-1] + head(tpr_cut, -1))) / 2)
}
