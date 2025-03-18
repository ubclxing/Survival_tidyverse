\name{predict.multiFit}
\alias{predict.multiFit}
\title{Make predictions for multiple outcomes}
\description{
This function makes predictions from a multiFit object.}
\usage{
\method{predict}{multiFit}(object, newdata, ...)
}
\arguments{
  \item{object}{A fitted object from \code{"multiFit"}}

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
fit <- multiFit(xmat = x_train, ymat = y_train,
                method = rpart1, family = "gaussian")
predict(fit, x_test)
}
