\name{simdat_mtps}
\alias{simdat_mtps}
\alias{xmat}
\alias{ymat}
\title{A Simulated Database for Survival analysis}
\description{
The data is simulated for predicting multiple survival outcomes. \code{ymat} contains 2 survival response variables with 10\% censoring for 100 observations and \code{xmat} contains 2 predictors of those observations.
}
\details{
In the simdat_mtps database, the simulated survival time and status are used as multivariate outcomes, including \strong{time01},\strong{status01},\strong{time02},\strong{status02}. The predictors both follow a standard normal distribution. In this simulated data, "ymat" refers multiple survival outcomes and "xmat" refers the predictor data.
}
\examples{
data(simdat_mtps)
}
\keyword{datasets}
