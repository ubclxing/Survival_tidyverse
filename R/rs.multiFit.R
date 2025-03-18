#' @export rs_multiFit
rs_multiFit <- function(yhat, ymat, xmat = NULL,
                        family,
                        resid_type = c("deviance", "pearson", "raw"), resid_std = FALSE,
                        method, dist1 = NULL, weights = NULL) {
  
  resid_type <- match.arg(resid_type)
  
  ny <- ncol(ymat)
  
  # Check family input
  if (length(family) == 1) {
    if (!family %in% c("gaussian", "binomial", "survival")) {
      stop("family must be gaussian, binomial or survival")
    }
    family <- rep(family, ny)
  }
  
  if (length(family) != ny) {
    stop("length of family must be consistent with response")
  }
  if (sum(family %in% c("gaussian", "binomial", "survival")) != ny) {
    stop("family must be gaussian, binomial or survival or their combination")
  }
  
  if (length(method) == 1) {
    method <- rep(list(method), ny)
  }
  
  y_fitted <- ymat
  y_fitted[!is.na(y_fitted)] <- NA
  models <- vector("list", ny)
  colnames(y_fitted) <- names(models) <- colnames(ymat)
  fit <- vector("list", ny)
  
  # Ensure column names of predictor matrix are correct
  colnames(yhat) <- paste0("X", 1:ny)
  
  # Derive residuals to be predicted by other outcome variables
  for (kk in 1:ny) {
    if (family[kk] == "gaussian") {
      resi <- ymat - yhat
    }
    if (family[kk] == "binomial") {
      resi <- resid_bin(ymat, yhat, xmat, type = resid_type, resid_std = resid_std)
    }
    
    # Prepare xmat for fitting the second stage model
    xmat <- yhat[, -kk]
    if (ny == 2) {
      xmat <- as.data.frame(xmat)
      names(xmat) <- paste0("X", c(1:ny)[-kk])
    }
    
    fit[[kk]] <- method[[kk]](xmat = xmat, ymat = resi[, kk], family = "gaussian", weights = weights[, kk])
    models[[kk]] <- fit[[kk]]$model
    y_fitted[, kk] <- fit[[kk]]$y_fitted
  }
  
  rs_multi_fit_fits <- list(fit = fit,
                            y_fitted = y_fitted,
                            family = family,
                            models = models,
                            method = method,
                            resid_type = resid_type,
                            resid_std = resid_std)
  
  class(rs_multi_fit_fits) <- "rs.multiFit"
  return(rs_multi_fit_fits)
}

#' @method predict rs.multiFit
#' @export predict.rs.multiFit
#' @export
predict.rs.multiFit <- function(object, newdata, ...) {
  ny <- ncol(object$y_fitted)
  resid_mat <- newdata
  resid_mat[!is.na(resid_mat)] <- NA
  mtd <- object$method
  pred <- newdata
  pred[!is.na(pred)] <- NA
  
  colnames(newdata) <- paste0("X", 1:ny)
  
  for (ii in 1:ny) {
    xx <- newdata[, -ii]
    
    if (ny == 2) {
      xx <- as.data.frame(xx)
      names(xx) <- paste0("X", c(1:ny)[-ii])
    }
    
    model <- object$model[[ii]]
    resid_mat[, ii] <- object$fit[[ii]]$predFun(model, xx)
  }
  
  cindex <- object$family == "gaussian"
  bindex <- object$family == "binomial"
  
  # Continuous outcome
  pred[, cindex] <- newdata[, cindex] + resid_mat[, cindex]
  
  # Binary outcome
  if (sum(bindex) != 0) {
    if (object$resid_std) {
      QQ <- qr.Q(qr(newdata[, bindex]))  # Q matrix of QR decomposition
      hats <- rowSums(QQ^2)  # Diagonal of hat matrix X(X'X)^(-1)X', used to standardize residuals
    } else {
      hats <- rep(0, nrow(matrix(newdata[, bindex])))  # Dr. Xing revised the code
    }
    
    if (object$resid_type == "pearson") {
      pred[, bindex] <- newdata[, bindex] + resid_mat[, bindex] * sqrt(newdata[, bindex] * (1 - newdata[, bindex])) * sqrt(1 - hats)
    }
    if (object$resid_type == "raw") {
      pred[, bindex] <- newdata[, bindex] + resid_mat[, bindex]
    }
    if (object$resid_type == "deviance") {
      dev0 <- -sqrt(-2 * log(1 - newdata[, bindex]))
      dev1 <- sqrt(-2 * log(newdata[, bindex]))
      
      # Avoid infinity residuals
      m1 <- max_finite(dev1)
      m2 <- min_finite(dev0)
      dev1[dev1 > m1] <- m1 + 100
      dev0[dev0 < m2] <- m2 - 100
      
      res0 <- abs(resid_mat[, bindex] - dev0)
      res1 <- abs(resid_mat[, bindex] - dev1)
      pred[, bindex] <- res0 / (res0 + res1)
    }
    
    pred[, bindex][pred[, bindex] > 1] <- 1
    pred[, bindex][pred[, bindex] < 0] <- 0
  }
  
  return(pred)
}
