\name{cureboost}
\alias{cureboost}
\title{Fit Tail-Targeted Boosting for Cure Probability}
\usage{
cureboost(formula, data, nrounds = 100L, learning.rate = 0.05,
  maxdepth = 1L, nodesize = NULL, mtry = NULL, nsplit = 10L,
  tail.start = NULL, tail.quantile = 0.8, tail.end = NULL,
  tail.times = NULL, tail.grid.size = 20L,
  min.tail.at.risk = NULL, min.tail.end.at.risk = NULL,
  min.events.before.tail = 5L, lambda = 1e-6, min.gain = 0,
  ipcw.model = c("cox", "marginal", "user"), g.fun = NULL,
  ipcw.folds = 5L, g.lower = 0.02, validation.fraction = 0.2,
  early.stopping.rounds = 15L, stopping.tolerance = 1e-8,
  min.correlation = 1e-8, seed = 1L,
  na.action = stats::na.omit, keep.data = FALSE, verbose = FALSE)
}
\arguments{
  \item{formula}{A two-sided formula with a
  \code{survival::Surv(time, status)} response.}
  \item{data}{A data frame containing the response and predictors.}
  \item{nrounds}{Maximum number of shallow boosting trees.}
  \item{learning.rate}{Shrinkage in \eqn{(0,1]}.}
  \item{maxdepth}{Maximum depth of each weak tree. Depth one gives stumps.}
  \item{nodesize}{Minimum observations in each daughter node.}
  \item{mtry}{Number of encoded predictors considered at each node.}
  \item{nsplit}{Number of quantile cutpoints considered per predictor.}
  \item{tail.start, tail.end}{Start and end of the finite-tail window.}
  \item{tail.quantile}{Event-time quantile used when \code{tail.start} is
  omitted.}
  \item{tail.times}{Optional prespecified tail grid.}
  \item{tail.grid.size}{Maximum number of tail-grid times.}
  \item{min.tail.at.risk, min.tail.end.at.risk}{Tail maturity minima imposed
  on every candidate daughter node.}
  \item{min.events.before.tail}{Minimum failures before the tail begins.}
  \item{lambda}{Positive residual-score stabilizer.}
  \item{min.gain}{Minimum residual split score.}
  \item{ipcw.model}{Cross-fitted Cox, marginal KM, or user-supplied censoring
  survival model.}
  \item{g.fun}{Function of time and data used when
  \code{ipcw.model = "user"}.}
  \item{ipcw.folds}{Number of censoring-model folds.}
  \item{g.lower}{Lower truncation bound for censoring survival.}
  \item{validation.fraction}{Fraction reserved for early stopping. Set to zero
  to use the full sample for tree fitting.}
  \item{early.stopping.rounds}{Number of consecutive validation iterations
  without improvement before stopping.}
  \item{stopping.tolerance}{Minimum empirical or validation loss improvement.}
  \item{min.correlation}{Minimum normalized gradient correlation required to
  accept a weak tree.}
  \item{seed}{Random seed.}
  \item{na.action}{Training-data missing-value action.}
  \item{keep.data}{Whether to retain the analysis data.}
  \item{verbose}{Whether to report progress at four milestones.}
}
\value{An object of class \code{cureboost}.}
\description{
Fits an experimental stagewise ensemble for the observable finite-tail
survival target. Each shallow tree corrects the residual tail structure left
by the current ensemble while retaining the cure forest's follow-up maturity
restrictions.
}
\details{
The fitting target is a cross-fitted IPCW pseudo-outcome averaged over the
tail window. Candidate splits compare daughter-node mean residuals after
checking node size, late risk-set support, and events before the tail.
Accepted corrections must reduce empirical squared-error loss.

The response estimates finite-tail survival. It represents exact cure only
when residual survival among susceptible subjects is negligible over the
chosen tail. This initial implementation does not fit a separate latency
component and therefore does not return complete survival curves.
}
\examples{
data(cure_latency_demo)
train <- subset(cure_latency_demo, cohort == "train")
fit <- cureboost(
  survival::Surv(time, event) ~ cure_marker + latency_marker + age + stage,
  data = train,
  nrounds = 3,
  learning.rate = 0.1,
  tail.start = 8,
  tail.end = 10,
  nodesize = 60,
  min.tail.at.risk = 20,
  min.tail.end.at.risk = 5,
  validation.fraction = 0,
  ipcw.model = "marginal",
  seed = 7
)
predict(fit, train[1:5, ], type = "cure")
}
\seealso{\code{\link{cureforest}},
  \code{\link{predict.cureboost}},
  \code{\link{cureboost_importance}}}
\keyword{survival}
\keyword{tree}
