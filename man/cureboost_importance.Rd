\name{cureboost_importance}
\alias{cureboost_importance}
\title{Variable Importance for Tail-Targeted Cure Boosting}
\usage{
cureboost_importance(object, aggregate = TRUE)
}
\arguments{
  \item{object}{A fitted \code{cureboost} object.}
  \item{aggregate}{Whether to aggregate encoded model-matrix columns to
  formula terms.}
}
\value{
A data frame with split counts, weighted total and mean gains, and root-split
frequencies.
}
\description{
Summarizes residual split gains over the selected boosting trees. Split scores
are weighted by the squared coefficient of the corresponding stagewise
update.
}
\seealso{\code{\link{cureboost}}}
