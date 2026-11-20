\name{predict.curecurveboost}
\alias{predict.curecurveboost}
\title{Predict a finite-horizon refined survival curve}
\usage{
\method{predict}{curecurveboost}(object, newdata,
  type = c("survival", "susceptible", "cure", "all"), times = NULL, ...)
}
\arguments{
\item{object}{An object returned by \code{curecurveboost}.}
\item{newdata}{A data frame with the model predictors.}
\item{type}{Population survival, susceptible survival, fixed tail anchor, or all outputs.}
\item{times}{Times between zero and the fitted horizon, inclusive. Defaults to the nonzero interval boundaries.}
\item{...}{Reserved for compatibility.}
}
\description{
Returns nonincreasing survival curves on the fitted horizon. Population
survival is the fixed tail anchor plus its complement times susceptible survival.
Requests outside that horizon fail rather than imply unsupported extrapolation.
}
\value{
A matrix with subjects in rows and times in columns for curve predictions,
a numeric vector for \code{type = "cure"}, or a named list for \code{type = "all"}.
The cure output is exactly the original forest tail estimate, not a newly
identified or corrected lifetime probability.
}
\seealso{\code{curecurveboost}}
