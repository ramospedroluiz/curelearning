\name{curelearning-package}
\alias{curelearning-package}
\alias{curelearning}
\docType{package}
\title{Statistical Learning for Cure Models and Long-Term Survival}
\description{
Fits nearest-neighbor, tree, forest and boosting models for right-censored
survival data. Nearest-neighbor predictions target survival at specified late
times and may be withheld when local follow-up is insufficient.
}
\details{
Use \code{\link{cureknn}} for late-survival nearest neighbors. Tree-based methods include \code{curetree()} for a single tree,
\code{cureforest()} for a randomized ensemble, and \code{cureboost()} for a
sequence of residual tail corrections. The original functions
\code{survcure()} and \code{randomforestcure()} remain available for backward
compatibility. Node estimates use either a
tail-averaged Kaplan--Meier estimator or an inverse-probability-of-censoring
weighted estimator. Explicit maturity constraints control the number at risk
at the beginning and end of the tail window.

The cure probability is a population-level statistical estimand. It should not
be interpreted as an observed biological-cure label for an individual.
}
\author{
Pedro Luiz Ramos and Lucas Osses.
}
\seealso{
\code{\link{survcure}}, \code{\link{randomforestcure}},
\code{\link{cureboost}}, \code{\link{predict.survcure}},
\code{\link{cure_importance}}, and \code{\link{cureboost_importance}}.
}
\keyword{survival}
\keyword{tree}
