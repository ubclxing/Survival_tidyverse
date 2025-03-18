#' @title Cross-validation for Multi-Fit Model
#' @description Use cross-validation to evaluate a multi-fit model's performance.
#' @details The function performs cross-validation using different methods for fitting models and evaluates them based on the family type (gaussian, binomial, survival).
#' @param xmat Matrix of predictors, each row is an observation vector.
#' @param ymat Matrix of outcomes (quantitative for family = "gaussian", factor of two levels for family = "binomial", survival object for family = "survival").
#' @param xmat_list A list of predictor matrices (optional).
#' @param nfold Integer, number of folds for cross-validation.
#' @param method A model fitting method or a list of methods for each outcome.
#' @param family The response type for each response (can be "gaussian", "binomial", or "survival").
#' @param dist1 Assumed distribution for survival response (optional).
#' @param weights Weights for the survival response (optional).
#' @export cv_multiFit
#' @return A list containing the fitted models, predictions, and other model details.
#' @examples
#' cv_multiFit(xmat = XX, ymat = YY, family = "gaussian", nfold = 2, method = rpart1)

cv_multiFit <- function(xmat, ymat, xmat_list = NULL, nfold = 5, method, family = family, dist1 = NULL, weights = NULL) {
  
  ny <- ncol(ymat)
  nx <- ncol(xmat)
  
  # Check family input
  if (length(family) == 1) {
    if (!family %in% c("gaussian", "binomial", "survival")) {
      stop("family must be gaussian, binomial, or survival")
    }
    family <- rep(family, ny)
  }
  
  if (length(family) != ny) {
    stop("length of family must be consistent with response")
  }
  
  if (sum(family %in% c("gaussian", "binomial", "survival")) != ny) {
    stop("family must be gaussian, binomial, survival, or their combination")
  }
  
  # Check method consistency
  if (length(method) == 1) {
    method <- rep(list(method), ny)
  }
  
  if (length(method) != ny) {
    stop("length of method must be 1 or the same as the number of response columns")
  }
  
  # Initialize variables for survival vs other models
  if (sum(family %in% c("survival")) == 0) {
    ymat_stand <- scale(ymat)
    km1 <- kmeans(ymat_stand, 10, nstart = 100)
    idx_cv <- createFolds(factor(km1$cluster), k = nfold, list = FALSE)
    y_fitted <- ymat
    y_fitted[!is.na(y_fitted)] <- NA
    models <- vector("list", nfold)
    names(models) <- paste0("fold", 1:nfold)
    fits <- vector("list", nfold)
  } else {
    ymat_time <- matrix(NA, ncol = ny, nrow = nrow(xmat))
    ymat_status <- matrix(NA, ncol = ny, nrow = nrow(xmat))
    for (i in 1:ny) {
      ymat_time[, i] <- ymat[[i]][, 1]
      ymat_status[, i] <- ymat[[i]][, 2]
    }
    km1 <- kmeans(ymat_status, nfold, nstart = 100)
    idx_cv <- createFolds(factor(km1$cluster), k = nfold, list = FALSE)
    y_fitted <- ymat_time
    y_fitted[!is.na(y_fitted)] <- NA
    models <- vector("list", nfold)
    names(models) <- paste0("fold", 1:nfold)
    fits <- vector("list", nfold)
  }
  
  # Cross-validation for non-survival responses
  if (sum(family %in% c("survival")) == 0) {
    for (ii in 1:nfold) {
      # Make train and test data for the ii-th fold
      y_train <- ymat[idx_cv != ii, ]
      y_test <- ymat[idx_cv == ii, ]
      x_train <- xmat[idx_cv != ii, ]
      x_test <- xmat[idx_cv == ii, ]
      colnames(x_test) <- paste0("X", 1:ncol(x_test))
      
      # Fit the model for the ii-th fold
      fits[[ii]] <- multiFit(xmat = x_train, ymat = y_train, xmat_list = xmat_list, method, family)
      y_fitted[idx_cv == ii, ] <- predict.multiFit(fits[[ii]], x_test)
      models[[ii]] <- fits[[ii]]$model
    }
  } else {
    for (ii in 1:nfold) {
      # Make train and test data for the ii-th fold (survival)
      y_train_time <- ymat_time[idx_cv != ii, ]
      y_train_status <- ymat_status[idx_cv != ii, ]
      y_train <- vector("list", ny)
      for (j in 1:ny) {
        y_train[[j]] <- Surv(y_train_time[, j], y_train_status[, j])
      }
      
      y_test_time <- ymat_time[idx_cv == ii, ]
      y_test_status <- ymat_status[idx_cv == ii, ]
      y_test <- vector("list", ny)
      for (j in 1:ny) {
        y_test[[j]] <- Surv(y_test_time[, j], y_test_status[, j])
      }
      
      x_train <- xmat[idx_cv != ii, ]
      x_test <- xmat[idx_cv == ii, ]
      colnames(x_test) <- paste0("X", 1:ncol(x_test))
      
      cv_train <- cbind(y_train_time, y_train_status, x_train)
      names(cv_train)[1:(2 * ny)] <- c(paste0("ageph0", 1:ny), paste0("statusph0", 1:ny))
      
      if (is.null(weights)) {
        cv_weights <- NULL
      } else {
        cv_weights <- find_km_weights_mat(cv_train)
      }
      
      # Fit model for survival outcomes
      if (is.null(dist1)) {
        fits[[ii]] <- multiFit(xmat = x_train, ymat = y_train, xmat_list = xmat_list, method,
                               family, dist1 = dist1, weights = cv_weights)
      } else {
        fits[[ii]] <- multiFit(xmat = x_train, ymat = y_train, xmat_list = xmat_list,
                               method, family, dist1 = dist1, weights = cv_weights)
      }
      
      y_fitted[idx_cv == ii, ] <- predict.multiFit(fits[[ii]], x_test)
      models[[ii]] <- fits[[ii]]$model
    }
  }
  
  # Return the results as a list
  cv_multiFit_fits <- list(fit = fits, y_fitted = y_fitted, model = models, method = method, family = family)
  class(cv_multiFit_fits) <- "cv.multiFit"
  return(cv_multiFit_fits)
}

#' @method predict cv.multiFit
#' @export predict.cv.multiFit
#' @export
predict.cv.multiFit <- function(object, newdata, ...) {
  ny <- ncol(object$y_fitted)
  nfold <- length(object$model)
  temp0 <- array(NA, c(nrow(newdata), ny, nfold))
  clas <- object$method
  
  colnames(newdata) <- paste0("X", 1:ncol(newdata))
  for (ii in 1:ny) {
    for (jj in 1:nfold) {
      model <- object$model[[jj]][[ii]]
      temp0[, ii, jj] <- object$fit[[jj]][["fit"]][[ii]]$predFun(model, newdata)
    }
  }
  
  # Calculate the mean prediction
  pred <- apply(temp0, c(1, 2), mean)
  
  # Correct predicted probabilities for binary outcomes
  bindex <- object$family == "binomial"
  pred[, bindex][pred[, bindex] > 1] <- 1
  pred[, bindex][pred[, bindex] < 0] <- 0
  
  return(pred)
}
