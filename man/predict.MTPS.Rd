\name{predict.MTPS}
\alias{predict.MTPS}
\title{Make predictions from a "MTPS" model}
\description{
This function makes predictions from a revised stacking model.}
\usage{
\method{predict}{MTPS}(object, newdata, ...)
}
\arguments{
  \item{object}{A fitted object from \code{"MTPS"}}

  \item{newdata}{Matrix of new predictors at which predictions are to be made}

  \item{...}{additional arguments affecting the predictions produced}
}
\value{
  The predicted value from new predictors.
}

\examples{
data("HIV")
set.seed(1)
xmat <- as.matrix(XX)
ymat <- as.matrix(YY)
id <- createFolds(rowMeans(XX), k=5, list=FALSE)
training_id <- id != 1
y_train <- ymat[training_id, ]
y_test  <- ymat[!training_id, ]
x_train <- xmat[training_id, ]
x_test  <- xmat[!training_id, ]
# Cross-Validation Residual Stacking
fit_rs <- MTPS(xmat = x_train, ymat = y_train,
  family = "gaussian",cv = FALSE, residual = TRUE,
  method_step1 = rpart1, method_step2 = lm1)
pred_rs <- predict(fit_rs, x_test)

# Residual Stacking for Survival Analysis
  set.seed(1)
  data("simdat_mtps")

 # prepare training and test set
 id_train <- sample(1:100, 80)
 xmat_train <- xmat[id_train,]
 xmat_test <- xmat[-id_train,]
 ymat_train <- cbind(list(survival::Surv(ymat[id_train,"time01"],ymat[id_train,"status01"])),
 list(survival::Surv(ymat[id_train,"time02"],ymat[id_train,"status02"])))
 # Produce the Kaplan-Meier estimator
 weights <- find_km_weights_mat(ymat[id_train,],num_outcome = 2)
 # fit Residual Stacking Model for Survival Data
 fit <- MTPS(xmat_train, ymat_train, family = 'survival', cv=FALSE,
 residual = TRUE, nfold=5, method_step1 = surv,
 method_step2 = lm1, dist1 = "lognormal", weights = weights)
 # predict the survival time on test set
 pre <- predict.MTPS(fit, xmat_test)



}
