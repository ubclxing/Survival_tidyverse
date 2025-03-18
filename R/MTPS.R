#' @title Fit Models using Revised Stacking Algorithm
#' @aliases MTPS
#' @description Fit a model using standard stacking algorithm or revised stacking algorithms to simultaneously predict multiple outcomes.
#' @param xmat Predictor matrix, each row is an observation vector.
#' @param ymat Responses matrix. Quantitative for family = "gaussian". A factor of two levels for family = "binomial". A survival object for family = "survival".
#' @param xmat_list The user defines a list to specify the number of subset columns of the xmat for each outcome.
#' @param family Response type for each response. If all response variable are within the same family, it can be "gaussian", "binomial" or "survival", otherwise it is a vector with elements "gaussian", "binomial", and "survival" to indicate each response family.
#' @param cv Logical, indicate if use Cross-Validation Stacking algorithm.
#' @param residual Logical, indicate if use Residual Stacking algorithm.
#' @param nfold Integer, number of folds for Cross-Validation Stacking algorithm. The default value is 5.
#' @param method_step1 Base Learners for fitting models in Step 1 of Stacking Algorithm. It can be one base learner function for all outcomes or a list of base learner functions for each outcome. The list of all base learners can be obtained by \code{list.learners()}.
#' @param method_step2 Base Learners for fitting models in Step 2 of Stacking Algorithm.
#' @param resid_type The residual type for Residual Stacking.
#' @param resid_std Logical, whether or not use standardized residual.
#' @param dist1 Assumed distribution for response variable. If the argument is a character string, then it is assumed to name an element from \code{"survreg.distributions"}. These include "weibull", "exponential", "gaussian", "logistic", "lognormal", and "loglogistic".
#' @param weights A logical indicating whether to use weight for response in step1.
#' @param weight2 A logical indicating whether to use weight for response in step2. The default is FALSE.
#' @return It returns an MTPS object. It is a list of 4 parameters containing information about step 1 and step 2 models and the revised stacking algorithm method.
#' @export MTPS
#' @examples
#' data("HIV")
#' set.seed(1)
#' xmat <- as.matrix(XX)
#' ymat <- as.matrix(YY)
#' id <- createFolds(rowMeans(XX), k=5, list=FALSE)
#' training_id <- id != 1
#' y_train <- ymat[training_id, ]
#' y_test  <- ymat[!training_id, ]
#' x_train <- xmat[training_id, ]
#' x_test  <- xmat[!training_id, ]
#'
#' # Residual Stacking
#' fit_rs <- MTPS(xmat = x_train, ymat = y_train,
#'                family = "gaussian", cv = FALSE, residual = TRUE,
#'                method_step1 = rpart1, method_step2 = lm1)
#' pre1 <- predict(fit_rs, x_test)
#'
#' # Using different base learners for different outcomes
#' fit_mix_out <- MTPS(xmat = x_train, ymat = y_train,
#'                     family = "gaussian", cv = FALSE, residual = TRUE,
#'                     method_step1 = c(rpart1, lm1, rpart1, lm1, lm1),
#'                     method_step2 = c(rpart1, lm1, lm1, lm1, lm1))
#' pre2 <- predict(fit_mix_out, x_test)
#'
#' # Residual Stacking for Survival Analysis
#' set.seed(1)
#' data("simdat_mtps")
#' id_train <- sample(1:100, 80)
#' xmat_train <- xmat[id_train, ]
#' xmat_test <- xmat[-id_train, ]
#' ymat_train <- cbind(list(survival::Surv(ymat[id_train, "time01"], ymat[id_train, "status01"])),
#'                     list(survival::Surv(ymat[id_train, "time02"], ymat[id_train, "status02"])))
#' weights <- find_km_weights_mat(ymat[id_train, ], num_outcome = 2)
#'
#' xmat_list <- list(c(1), c(2))
#' fit <- MTPS(xmat_train, ymat_train, xmat_list = xmat_list, family = 'survival',
#'             cv = FALSE, residual = TRUE, nfold = 5, method_step1 = surv,
#'             method_step2 = lm1, dist1 = "lognormal", weights = weights)
#' pre3 <- predict.MTPS(fit, xmat_test)

MTPS <- function(xmat, ymat, xmat_list = NULL, family,
                 cv = FALSE, residual = TRUE,
                 nfold = 5, method_step1, method_step2,
                 resid_type = c("deviance", "pearson", "raw"),
                 resid_std = FALSE, dist1 = NULL, weights = NULL, weight2 = FALSE) {

  ny <- ncol(ymat)  # Number of outcomes

  # Check family input
  if (length(family) == 1) {
    if (!family %in% c("gaussian", "binomial", "survival")) {
      stop("family must be gaussian, binomial or survival")
    }
    family <- rep(family, ny)  # Replicate family type for each outcome
  }

  if (length(family) != ny) {
    stop("length of family must be consistent with response")
  }

  # Check consistency of method_step1 and method_step2
  if (length(method_step1) == 1) {
    method_step1 <- rep(list(method_step1), ny)
  }
  if (length(method_step2) == 1) {
    method_step2 <- rep(list(method_step2), ny)
  }

  if (length(method_step1) != ny || length(method_step2) != ny) {
    stop("length of method_step1 and method_step2 must be 1 or the same as the response column")
  }

  # Step 1: Fit base learners using method_step1
  if (cv) {
    fit1 <- cv_multiFit(xmat = xmat, ymat = ymat, nfold = nfold,
                        method = method_step1, xmat_list = xmat_list,
                        family = family, dist1 = dist1, weights = weights)
  } else {
    fit1 <- multiFit(xmat = xmat, ymat = ymat,
                     method = method_step1, xmat_list = xmat_list,
                     family = family, dist1 = dist1, weights = weights)
  }

  pred1 <- fit1$y_fitted  # Predictions from step 1

  # Step 2: Prepare response matrix for step 2
  if (sum(family %in% c('survival')) == ny) {
    ymat_step2 <- matrix(NA, nrow = nrow(xmat), ncol = ny)
    for (tt in 1:ny) {
      ymat_step2[, tt] <- ymat[[tt]][, 1]
    }
    family <- rep("gaussian", ny)
    pred1 <- log(pred1)  # Log-transform predictions for step 2
    ymat_step2 <- log(ymat_step2)  # Log-transform response for step 2
  } else {
    ymat_step2 <- ymat
  }

  # Step 2: Fit base learners using method_step2
  if (weight2) {
    weight2_matrix <- lapply(1:ny, function(i) {
      cbind(ymat[[i]][, 2], (ymat_step2 - pred1)[, i])
    })
    w2 <- do.call(cbind, weight2_matrix)
    w2 <- as.data.frame(w2)

    # Assigning column names for weight2
    fields <- c("statusph", "ageph")
    numbers <- sprintf("%02d", 1:ny)
    names(w2) <- as.vector(outer(fields, numbers, paste0))

    weight1 <- find_km_weights_mat(w2)
  } else {
    weight1 <- NULL
  }

  # Step 2: Residual fitting
  if (residual) {
    fit2 <- rs_multiFit(yhat = pred1, ymat = ymat_step2, xmat = xmat,
                        family = family, resid_type = resid_type,
                        resid_std = resid_std,
                        method = method_step2, dist1 = NULL,
                        weights = weight1)
  } else {
    fit2 <- multiFit(xmat = pred1, ymat = ymat_step2,
                     method = method_step2, family = family,
                     dist1 = NULL, weights = weight1)
  }

  fit <- list(fit1 = fit1, fit2 = fit2,
              cv = cv, residual = residual)
  class(fit) <- "MTPS"
  return(fit)
}

#' @method predict MTPS
#' @export predict.MTPS
#' @export
predict.MTPS <- function(object, newdata, ...) {

  if (object$cv) {
    pred1 <- predict(object$fit1, newdata, ynew = NULL)
  } else {
    pred1 <- predict(object$fit1, newdata, ynew = NULL)
  }

  # reverse the log-transformation
  if (sum(object$fit1$family %in% c('survival')) == length(object$fit1$family)) pred1 <- log(pred1)

  if (object$residual) {
    pred2 <- predict(object$fit2, pred1, object$mybest, ynew = NULL)
  } else {
    pred2 <- predict(object$fit2, pred1, object$mybest, ynew = NULL)
  }
  # reverse the log-transformation
  if (sum(object$fit$family %in% c('survival')) == length(object$fit$family))
    pred2 <- exp(pred2)

  return(pred2)
}
