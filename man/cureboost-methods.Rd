\name{cureboost-methods}
\alias{predict.cureboost}
\alias{print.cureboost}
\alias{summary.cureboost}
\alias{print.summary.cureboost}
\title{Methods for Tail-Targeted Cure Boosting}
\usage{
\method{predict}{cureboost}(object, newdata = NULL,
  type = c("cure", "response", "raw", "terminal", "contributions", "all"),
  nrounds = NULL, ...)

\method{print}{cureboost}(x, ...)

\method{summary}{cureboost}(object, ...)

\method{print}{summary.cureboost}(x, ...)
}
\arguments{
  \item{object, x}{A fitted \code{cureboost} object, or its summary.}
  \item{newdata}{Optional data frame of predictors.}
  \item{type}{Clipped finite-tail prediction, raw additive prediction,
  terminal-node identifiers, per-tree contributions, or all components.}
  \item{nrounds}{Number of selected boosting trees used for prediction.}
  \item{...}{Additional arguments reserved for compatibility.}
}
\value{
The predict method returns a numeric vector, matrix, or component list
according to \code{type}. The print methods return their input invisibly, and
the summary method returns an object of class \code{summary.cureboost}.
}
\description{
Prediction, printing, and summary methods for tail-targeted cure boosting.
}
\details{
\code{type = "cure"} is a convenient name retained for compatibility with
\code{cureforest}; formally it is the clipped finite-tail prediction. For
\code{type = "contributions"}, row sums equal the raw additive prediction.
}
\seealso{\code{\link{cureboost}}, \code{\link{cureboost_importance}}}
