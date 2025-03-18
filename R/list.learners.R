#' @title List Available Base Learners
#' @name list.learners
#' @description This function lists all base learners provided in the package.
#' @details 
#' - lm1: Linear regression
#' - glm1: Generalized linear models
#' - glmnet1: Performs k-fold cross-validation to choose the best alpha and lambda for generalized linear models via penalized maximum likelihood.
#' - glmnet_lasso: LASSO, lambda is chosen by k-fold cross-validation for glmnet.
#' - glmnet_ridge: Ridge regression, lambda is chosen by k-fold cross-validation for glmnet.
#' - rpart1: Regression tree
#' - lda1: Linear discriminant analysis
#' - qda1: Quadratic discriminant analysis
#' - KNN1: K-nearest neighbor classification, k is chosen by cross-validation.
#' - svm1: Support vector machine
#' - surv: Parametric AFT model
#' - ela1: Elastic Net AFT model, weights are the weight of survival time.
#' @export list_learners
#' @return The name of all base learners provided in the package.
#' @examples
#' list_learners()

list_learners <- function() {
  # Define the available base learners
  learners <- c("lm1", "glm1", "glmnet1", "glmnet_lasso", "glmnet_ridge", "rpart1",
                "lda1", "qda1", "KNN1", "svm1", "ela1", "surv", "lasso1")
  
  # Print the models that can be chosen
  message("Models that can be chosen are:")
  print(learners)
}
