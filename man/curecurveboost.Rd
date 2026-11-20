\name{curecurveboost}
\alias{curecurveboost}
\title{Flexible finite-horizon refinement of a CureForest survival curve}
\usage{
curecurveboost(formula, data, breaks, nrounds = 150L,
  learning.rate = 0.1, maxdepth = 2L, minbucket = 30L,
  validation.fraction = 0.2, anchor.fraction = 0.5,
  forest.control = list(), seed = 1L, verbose = FALSE,
  selection = c("gated", "minimum"), gate.fraction = 0.5,
  gate.alpha = 0.05, minimum.gain = 0)
}
\arguments{
\item{formula}{A formula with a right-censored \code{survival::Surv} response.}
\item{data}{A data frame with complete response and predictor values.}
\item{breaks}{Strictly increasing interval boundaries, starting at zero. The last boundary is the prediction horizon. Choose these without using test outcomes.}
\item{nrounds}{Maximum number of correction trees. Zero is allowed.}
\item{learning.rate}{Maximum shrinkage step, strictly between zero and one.}
\item{maxdepth}{Maximum depth of each correction tree.}
\item{minbucket}{Minimum number of person-interval records per correction leaf.}
\item{validation.fraction}{Total fraction reserved for selection and, by default, independent acceptance.}
\item{anchor.fraction}{Fraction used only to fit the CureForest anchor. The remainder after validation is used to fit corrections.}
\item{forest.control}{Named arguments passed to \code{randomforestcure}, excluding \code{formula} and \code{data}. Includes tail-window and maturity controls.}
\item{seed}{Seed for the disjoint split and anchor fit.}
\item{verbose}{Whether to report correction progress at quarter increments.}
\item{selection}{The default \code{"gated"} chooses a candidate on one validation subset and accepts it only after an independent paired comparison. \code{"minimum"} reproduces the earlier minimum-validation-score rule.}
\item{gate.fraction}{Fraction of the validation pool reserved for the independent acceptance sample.}
\item{gate.alpha}{One-sided significance level for the approximate paired acceptance bound; default 0.05.}
\item{minimum.gain}{Minimum absolute mean Brier reduction required by the lower bound; default zero. Prespecify this without inspecting test outcomes.}
}
\description{
Refines the susceptible survival curve implied by a cure forest, while
preserving its estimated tail probability. Time and predictors enter the
correction trees jointly; proportional hazards are not imposed.
}
\details{
By default four disjoint samples are used for the anchor, correction,
candidate selection, and independent acceptance. The \code{"minimum"} rule
uses three samples and no independent acceptance step.
The anchor curve is converted to piecewise constant susceptible hazards.
Small positive hazard regularization allows correction in initially flat
intervals. The log hazard is updated by shallow regression trees fitted to
posterior-susceptibility-weighted working responses. A line search accepts
only decreases in the observed mixture negative log likelihood, with exact
partial exposure for censoring inside an interval. These are EM-surrogate
working weights, not an assertion that the observed Hessian is positive.

Validation selects the number of rounds using a mean IPCW Brier score at the
nonzero interval boundaries. The regularized grid representation with zero
corrections is a candidate. This representation can differ slightly from the
raw forest. Default validation weights use a censoring Kaplan--Meier fit on
the correction sample and require censoring independent of the event time and
predictors; conditional censoring requires suitable conditional weights in a
separately specified evaluation pipeline.

Under \code{"gated"}, only the single selected candidate is evaluated on the
acceptance sample. Subject-level differences in mean IPCW Brier loss between
round zero and the candidate are computed, keeping all horizons of a subject
together. The candidate is accepted only if
\eqn{\bar d-t_{1-\alpha,n-1}s_d/\sqrt{n}>\text{minimum.gain}}.
Otherwise the regularized round-zero curve is returned. No further model is
tried on that acceptance sample. This paired t bound is approximate, conditional
on the fitted models and censoring weights; it is not an exact finite-sample
error-control theorem. The acceptance sample must not be reused for a
hyperparameter search, censoring-model tuning, or repeated redesign.

The incidence anchor is fixed, so correction does not improve its cure
probabilities. A tail estimate represents biological cure only under the
identification assumptions and adequate follow-up supporting that interpretation.
No extrapolation beyond the specified horizon is supplied. Neither consistency
nor an inferential guarantee for this new correction follows from the forest
theorems. Preprocessing must also exclude validation and test outcomes.
}
\value{
An object of class \code{curecurveboost} containing the anchor, selected trees,
intervals, fitting history, selection rule, acceptance diagnostics, and indices
of the disjoint samples. With gated selection, \code{candidate_nrounds} records
the candidate length and \code{nrounds} records the accepted length (zero when
rejected). Acceptance diagnostics are stored in \code{gate}.
}
\seealso{\code{randomforestcure}, \code{predict.curecurveboost}}
