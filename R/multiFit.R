#' @title Fit models on multiple outcomes.
#' @param xmat Matrix of predictors, each row is an observation vector.
#' @param ymat Matrix of outcomes.
#' Quantitative for family = "gaussian".
#' A factor of two levels for family = "binomial".
#' A survival object for family = "survival".
#' @param xmat_list The user defines a numeric vector to specify the number of subset columns of the xmat.
#' @param method Method for fitting models.
#' It can be one base learner function for all outcomes or a list of base learner functions for each outcome.
#' The list of all base learners can be obtained by \code{list.learners()}.
#' @param family Response type for each response.
#' If all response variable are within the same family it can be "gaussian", "binomial" or "survival",
#' otherwise it is a vector with elements "gaussian", "binomial" and "survival" to indicate each response family.
#' @param dist1 Assumed distribution for response variable.
#' If the argument is a character string, then it is assumed to name an element from \code{"survreg.distributions"}.
#' These include "weibull", "exponential", "gaussian", "logistic", "lognormal" and "loglogistic".
#' @param weights Weight for response.
#' @return It returns a multiFit object. It is a list of 6 parameters containing information about the fitted models and fitted values for each outcome.
#' @description This function fits individual models to predict each outcome separately.
#' @export multiFit
#' @examples
#' data("HIV")
#' set.seed(1)
#' xmat <- as.matrix(XX)
#' ymat <- as.matrix(YY)
#' id <- createFolds(rowMeans(XX), k = 5, list = FALSE)
#' training_id <- id != 1
#' y_train <- ymat[training_id, ]
#' y_test  <- ymat[!training_id, ]
#' x_train <- xmat[training_id, ]
#' x_test  <- xmat[!training_id, ]
#' fit <- multiFit(xmat = x_train, ymat = y_train,
#'                 method = rpart1, family = "gaussian")
#' pre1 <- predict(fit, x_test)
#' 
#' # Using different base learners for different outcomes
#' fit_mix_out <- multiFit(xmat = x_train, ymat = y_train,
#'                        method = c(rpart1, rpart1, lm1, lm1, lm1),
#'                        family = "gaussian")
#' pre2 <- predict(fit_mix_out, x_test)

multiFit <- function(xmat, ymat, xmat_list = NULL,
                     method, family, dist1 = NULL, weights = NULL) {
  
  nx <- ncol(xmat)
  ny <- ncol(ymat)

  # Check family input
  if (length(family) == 1) {
    if (!family %in% c("gaussian", "binomial", "survival")) {
      stop("family must be gaussian, binomial, or survival")
    }
    family <- rep(family, ny)
  }
  
  # Check consistency of method for each outcome
  if (length(method) == 1) {
    method <- rep(list(method), ny)
  }
  
  # Initialize matrices to save fitted results
  y_fitted <- matrix(NA, nrow = nrow(xmat), ncol = ny)
  models <- vector("list", ny)
  colnames(y_fitted) <- names(models) <- colnames(ymat)
  fit <- vector("list", ny)
  mybest <- vector("list", ny)
  
  # Ensure that column names of predictor matrix are correct
  colnames(xmat) <- paste0("X", 1:nx)
  xmat0 <- xmat
  
  # Check xmat_list
  if (is.null(xmat_list)) {
    xmat_list <- vector("list", ny)
    for (kk in 1:ny) {
      xmat_list[[kk]] <- 1:nx
    }
  }
  
  # Fit models for each outcome
  for (kk in 1:ny) {
   
    xmat <- xmat0
    xmat <- xmat[, xmat_list[[kk]]]
    xmat <- as.data.frame(xmat)
    names(xmat) <- colnames(xmat0)[xmat_list[[kk]]]
    
    # Fit model based on family type and distribution
    if (is.null(dist1)) {
      if (family[kk] == "survival") {
        fit[[kk]] <- method[[kk]](xmat, ymat[[kk]], family = family[kk], weights = weights[, kk])
      } else {
        fit[[kk]] <- method[[kk]](xmat, ymat[, kk], family = family[kk], weights = weights[, kk])
      }
      models[[kk]] <- fit[[kk]]$model
      y_fitted[, kk] <- fit[[kk]]$y_fitted
      mybest[[kk]] <- fit[[kk]]$mybest
    } else {
      if (family[kk] == "survival") {
        fit[[kk]] <- method[[kk]](xmat, ymat[[kk]], family = family[kk], weights = weights[, kk], dist1 = dist1)
      } else {
        fit[[kk]] <- method[[kk]](xmat, ymat[, kk], family = family[kk], weights = weights[, kk], dist1 = dist1)
      }
      models[[kk]] <- fit[[kk]]$model
      y_fitted[, kk] <- fit[[kk]]$y_fitted
      mybest[[kk]] <- fit[[kk]]$mybest
    }
  }
  
  # Return multiFit object
  multiFit_fits <- list(fit = fit, xmat_list = xmat_list,
                        y_fitted = y_fitted,
                        model = models,
                        mybest = mybest,
                        method = method,
                        family = family)
  class(multiFit_fits) <- "multiFit"
  return(multiFit_fits)
}

#' @method predict multiFit
#' @export predict.multiFit
#' @export
predict.multiFit <- function(object, newdata, ...) {
  ny <- length(object$model)
  newdata <- as.data.frame(newdata)
  pred <- matrix(NA, nrow = nrow(newdata), ncol = ny)
  
  colnames(newdata) <- paste0("X", 1:ncol(newdata))
  newdata0 <- newdata
  
  for (ii in 1:ny) {
    newdata <- newdata0
    if (is.null(object$xmat_list)) {
      xmat_list[[ii]] <- 1:ncol(newdata)
    } else {
      xmat_list <- object$xmat_list
    }
    newdata <- as.data.frame(newdata[, object$xmat_list[[ii]]])
    names(newdata) <- names(newdata0)[xmat_list[[ii]]]
    model <- object$model[[ii]]
    pred[, ii] <- object$fit[[ii]]$predFun(model, newdata)
    
    if (object$family[ii] == "binomial") {
      # Predicted probabilities should be within [0,1]
      pred[, ii][pred[, ii] > 1] <- 1
      pred[, ii][pred[, ii] < 0] <- 0
    }
  }
  return(pred)
}
