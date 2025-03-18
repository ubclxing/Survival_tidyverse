#' @title Evaluation using Cross-Validation
#' @description Use cross-validation to evaluate model performance.
#' @details The most frequently used evaluation metric of survival models is the concordance index (C-index).
#' It is a measure of rank correlation between predicted risk scores. Only uncensored observations are used.
#' @param xmat Matrix of predictors, each row is an observation vector.
#' @param ymat Matrix of outcomes.
#' Quantitative for family = \strong{"gaussian"}.
#' A factor of two levels for family = \strong{"binomial"}.
#' A survival object for family = \strong{"survival"}.
#' @param family Response type for each response.
#' If all response variables are within the same family, it can be \strong{"gaussian"}, \strong{"binomial"} or \strong{"survival"},
#' otherwise it is a vector with elements \strong{"gaussian"}, \strong{"binomial"} and \strong{"survival"} to indicate each response family.
#' @param nfolds Integer, number of folds for Cross-Validation to evaluate the performance of stacking algorithms.
#' @param cv Logical, indicate if use Cross-Validation Stacking algorithm.
#' @param residual Logical, indicate if use Residual Stacking algorithm.
#' @param cv_stacking_nfold Integer, number of folds for Cross-Validation Stacking algorithm. The default value is 5.
#' @param method_step1 Base Learners for fitting models in Step 1 of Stacking Algorithm. It can be one base learner function for all outcomes or a list of base learner functions for each outcome. The list of all base learners can be obtained by \code{list.learners()}.
#' @param method_step2 Base Learners for fitting models in Step 2 of Stacking Algorithm.
#' @param resid_type The residual type for Residual Stacking.
#' @param resid_std Logical, whether or not to use standardized residuals.
#' @param dist1 Assumed distribution for survival response.
#' If the argument is a character string, it is assumed to name an element from survreg.distributions.
#' These include \strong{"weibull"}, \strong{"exponential"}, \strong{"gaussian"}, \strong{"logistic"}, \strong{"lognormal"} and \strong{"loglogistic"}.
#' Otherwise, it is assumed to be a user-defined list conforming to the format described in \code{"survreg.distributions"}.
#' @param weights Weight for survival response.
#' @export cv_MTPS
#' @return It returns the mean squared error of continuous outcomes, and AUC, accuracy, recall, and precision for binary outcomes of predictions using cross-validation.
#' @examples
#' data("HIV")
#' cv_MTPS(xmat = XX, ymat = YY, family = "gaussian", nfolds = 2, method_step1 = rpart1, method_step2 = lm1)

cv_MTPS <- function(xmat, ymat, family, nfolds = 5, cv = FALSE, residual = TRUE,
                    cv_stacking_nfold = 5, method_step1, method_step2,
                    resid_type = c("deviance", "pearson", "raw"), resid_std = FALSE,
                    dist1 = NULL, weights = NULL) {
  
  resid_type <- match.arg(resid_type)
  
  ny <- ncol(ymat)
  
  # Check family input consistency
  if (length(family) == 1) {
    if (!family %in% c("gaussian", "binomial", "survival")) {
      stop("family must be gaussian, binomial, or survival")
    }
    family <- rep(family, ny)
  }
  
  # Validate family length consistency with response variables
  if (length(family) != ny) {
    stop("length of family must be consistent with response")
  }
  if (sum(family %in% c("gaussian", "binomial", "survival")) != ny) {
    stop("family must be gaussian, binomial, survival, or their combination")
  }
  
  # Ensure method consistency for step1 and step2
  method_step1 <- rep(list(method_step1), ny)
  method_step2 <- rep(list(method_step2), ny)
  
  # Validate method lengths
  if (length(method_step1) != ny) stop("length of method_step1 must be 1 or the same as response columns")
  if (length(method_step2) != ny) stop("length of method_step2 must be 1 or the same as response columns")
  
  # Metrics initialization
  cindex <- which(family == "gaussian")
  bindex <- which(family == "binomial")
  
  metrics_ctn <- matrix(NA, nrow = 1, ncol = length(cindex))
  colnames(metrics_ctn) <- colnames(ymat[, cindex])
  rownames(metrics_ctn) <- "MSE"
  metrics_ctn <- as.data.frame(metrics_ctn)
  
  metrics_bin <- matrix(NA, nrow = 4, ncol = length(bindex))
  colnames(metrics_bin) <- colnames(ymat[, bindex])
  rownames(metrics_bin) <- c("AUC", "Accuracy", "Recall", "Precision")
  metrics_bin <- as.data.frame(metrics_bin)
  
  # Cross-validation setup
  idx_cv <- createFolds(rowMeans(xmat), k = nfolds, list = FALSE)
  pred <- ymat
  pred[!is.na(pred)] <- NA
  
  for (i_fold in 1:nfolds) {
    # Train and test data for the current fold
    y_train <- ymat[idx_cv != i_fold, ]
    y_test <- ymat[idx_cv == i_fold, ]
    x_train <- xmat[idx_cv != i_fold, ]
    x_test <- xmat[idx_cv == i_fold, ]
    
    # Fit the model using MTPS
    fit <- MTPS(xmat = x_train, ymat = y_train, family = family,
                cv = cv, residual = residual,
                nfold = cv_stacking_nfold,
                method_step1 = method_step1,
                method_step2 = method_step2,
                resid_type = resid_type, resid_std = resid_std, dist1 = dist1, weights = weights)
    
    # Make predictions
    pred[idx_cv == i_fold, ] <- predict(fit, x_test)
  }
  
  # Metrics computation for continuous outcomes
  if (length(cindex) > 0) {
    metrics_ctn["MSE", ] <- apply((pred[, cindex] - ymat[, cindex])^2, 2, mean)
  }
  
  # Metrics computation for binary outcomes
  for (jj in bindex) {
    metrics_bin["AUC", which(jj == bindex)] <- AUC(pred[, jj], outcome = ymat[, jj])
    table <- table((pred[, jj] > 0.5) * 1, ymat[, jj])
    metrics_bin["Accuracy", which(jj == bindex)] <- (table[1, 1] + table[2, 2]) / sum(table)
    metrics_bin["Recall", which(jj == bindex)] <- table[2, 2] / (table[2, 2] + table[1, 2])
    metrics_bin["Precision", which(jj == bindex)] <- table[2, 2] / (table[2, 2] + table[2, 1])
  }
  
  metrics <- list(continuous = metrics_ctn, binary = metrics_bin)
  return(metrics)
}
