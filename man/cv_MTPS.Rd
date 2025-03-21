% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.MTPS.R
\name{cv_MTPS}
\alias{cv_MTPS}
\title{Evaluation using Cross-Validation}
\usage{
cv_MTPS(
  xmat,
  ymat,
  family,
  nfolds = 5,
  cv = FALSE,
  residual = TRUE,
  cv_stacking_nfold = 5,
  method_step1,
  method_step2,
  resid_type = c("deviance", "pearson", "raw"),
  resid_std = FALSE,
  dist1 = NULL,
  weights = NULL
)
}
\arguments{
\item{xmat}{Matrix of predictors, each row is an observation vector.}

\item{ymat}{Matrix of outcomes.
Quantitative for family = \strong{"gaussian"}.
A factor of two levels for family = \strong{"binomial"}.
A survival object for family = \strong{"survival"}.}

\item{family}{Response type for each response.
If all response variables are within the same family, it can be \strong{"gaussian"}, \strong{"binomial"} or \strong{"survival"},
otherwise it is a vector with elements \strong{"gaussian"}, \strong{"binomial"} and \strong{"survival"} to indicate each response family.}

\item{nfolds}{Integer, number of folds for Cross-Validation to evaluate the performance of stacking algorithms.}

\item{cv}{Logical, indicate if use Cross-Validation Stacking algorithm.}

\item{residual}{Logical, indicate if use Residual Stacking algorithm.}

\item{cv_stacking_nfold}{Integer, number of folds for Cross-Validation Stacking algorithm. The default value is 5.}

\item{method_step1}{Base Learners for fitting models in Step 1 of Stacking Algorithm. It can be one base learner function for all outcomes or a list of base learner functions for each outcome. The list of all base learners can be obtained by \code{list.learners()}.}

\item{method_step2}{Base Learners for fitting models in Step 2 of Stacking Algorithm.}

\item{resid_type}{The residual type for Residual Stacking.}

\item{resid_std}{Logical, whether or not to use standardized residuals.}

\item{dist1}{Assumed distribution for survival response.
If the argument is a character string, it is assumed to name an element from survreg.distributions.
These include \strong{"weibull"}, \strong{"exponential"}, \strong{"gaussian"}, \strong{"logistic"}, \strong{"lognormal"} and \strong{"loglogistic"}.
Otherwise, it is assumed to be a user-defined list conforming to the format described in \code{"survreg.distributions"}.}

\item{weights}{Weight for survival response.}
}
\value{
It returns the mean squared error of continuous outcomes, and AUC, accuracy, recall, and precision for binary outcomes of predictions using cross-validation.
}
\description{
Use cross-validation to evaluate model performance.
}
\details{
The most frequently used evaluation metric of survival models is the concordance index (C-index).
It is a measure of rank correlation between predicted risk scores. Only uncensored observations are used.
}
\examples{
data("HIV")
cv_MTPS(xmat = XX, ymat = YY, family = "gaussian", nfolds = 2, method_step1 = rpart1, method_step2 = lm1)
}
