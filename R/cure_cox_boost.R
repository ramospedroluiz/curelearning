#' Fit cure-anchored Cox boosting
#'
#' Fits a semiparametric mixture-cure model in two connected stages. A
#' cure-fraction forest supplies the nonparametric incidence anchor unless an
#' existing cure model is provided. Conditional on that fixed cure
#' probability, shallow trees boost the susceptible Cox relative-risk
#' predictor using posterior susceptibility weights and a Breslow baseline.
#' The resulting marginal survival curve is coherent by construction.
#'
#' @param formula A formula with a `survival::Surv(time, status)` response.
#' @param data A data frame.
#' @param incidence.fit Optional fitted `randomforestcure`, `survcure`, or
#'   `cureboost` model. When omitted, a cure-fraction forest is fitted.
#' @param incidence.trees Number of trees in the automatically fitted incidence
#'   forest.
#' @param latency.rounds Maximum number of Cox boosting iterations.
#' @param latency.steps Number of latency trees per iteration.
#' @param learning.rate Shrinkage for latency trees.
#' @param incidence.maxdepth Maximum depth of incidence-anchor trees.
#' @param latency.maxdepth Maximum depth of Cox boosting trees.
#' @param nodesize Minimum node size; model defaults are used when `NULL`.
#' @param mtry Predictors considered at each split; defaults are used when
#'   `NULL`.
#' @param nsplit Candidate cut points per numeric predictor.
#' @param tail.start Start of the mature follow-up window.
#' @param tail.quantile Quantile used when `tail.start` is omitted.
#' @param tail.end End of the mature follow-up window.
#' @param tail.times Optional tail evaluation grid.
#' @param tail.grid.size Size of an automatically generated tail grid.
#' @param min.tail.at.risk Minimum number at risk at `tail.start`.
#' @param min.tail.end.at.risk Minimum number at risk at `tail.end`.
#' @param min.events.before.tail Minimum number of events before `tail.start`.
#' @param lambda Split regularization constant.
#' @param min.gain Minimum split gain.
#' @param validation.fraction Fraction held out for latency early stopping.
#' @param early.stopping.rounds Non-improving iterations before stopping.
#' @param stopping.tolerance Minimum validation improvement.
#' @param posterior.epsilon Numerical clipping constant.
#' @param honesty Whether the automatically fitted incidence forest is honest.
#' @param honesty.fraction Structure fraction for honest incidence trees.
#' @param sampsize Incidence-forest subsampling fraction or size.
#' @param replace Whether incidence trees sample with replacement.
#' @param n.cores Parallel workers for the incidence forest.
#' @param seed Reproducibility seed.
#' @param na.action Missing-data action.
#' @param keep.data Whether to retain training data.
#' @param verbose Whether to report progress at 25-percent increments.
#' @return An object inheriting from `curecoxboost` and `cureboost_joint`.
#' @export
curecoxboost <- function(formula,
                         data,
                         incidence.fit = NULL,
                         incidence.trees = 500L,
                         latency.rounds = 200L,
                         latency.steps = 1L,
                         learning.rate = 0.10,
                         incidence.maxdepth = 5L,
                         latency.maxdepth = 1L,
                         nodesize = NULL,
                         mtry = NULL,
                         nsplit = 10L,
                         tail.start = NULL,
                         tail.quantile = 0.80,
                         tail.end = NULL,
                         tail.times = NULL,
                         tail.grid.size = 20L,
                         min.tail.at.risk = NULL,
                         min.tail.end.at.risk = NULL,
                         min.events.before.tail = 5L,
                         lambda = 1e-6,
                         min.gain = 0,
                         validation.fraction = 0.20,
                         early.stopping.rounds = 20L,
                         stopping.tolerance = 1e-6,
                         posterior.epsilon = 1e-5,
                         honesty = TRUE,
                         honesty.fraction = 0.50,
                         sampsize = 0.80,
                         replace = FALSE,
                         n.cores = 1L,
                         seed = 1L,
                         na.action = stats::na.omit,
                         keep.data = FALSE,
                         verbose = FALSE) {
  call <- match.call()
  if (is.null(incidence.fit)) {
    incidence.fit <- randomforestcure(
      formula = formula,
      data = data,
      ntree = incidence.trees,
      splitrule = "cure",
      estimator = "km",
      tail.start = tail.start,
      tail.quantile = tail.quantile,
      tail.end = tail.end,
      tail.times = tail.times,
      tail.grid.size = tail.grid.size,
      curve.grid.size = max(40L, as.integer(tail.grid.size)),
      nodesize = nodesize,
      min.tail.at.risk = min.tail.at.risk,
      min.tail.end.at.risk = min.tail.end.at.risk,
      min.events.before.tail = min.events.before.tail,
      maxdepth = incidence.maxdepth,
      mtry = mtry,
      nsplit = nsplit,
      lambda = lambda,
      min.gain = min.gain,
      honesty = honesty,
      honesty.fraction = honesty.fraction,
      sampsize = sampsize,
      replace = replace,
      oob = FALSE,
      importance = TRUE,
      n.cores = n.cores,
      seed = seed,
      na.action = na.action,
      keep.data = FALSE
    )
  }

  object <- cureboost_joint(
    formula = formula,
    data = data,
    em.iterations = latency.rounds,
    incidence.steps = 0L,
    latency.steps = latency.steps,
    initialization.rounds = 0L,
    learning.rate.incidence = 0.05,
    learning.rate.latency = learning.rate,
    latency.baseline = "breslow",
    incidence.posterior.weight = 1,
    incidence.mode = "anchored",
    incidence.fit = incidence.fit,
    maxdepth = latency.maxdepth,
    nodesize = nodesize,
    mtry = mtry,
    nsplit = nsplit,
    tail.start = tail.start,
    tail.quantile = tail.quantile,
    tail.end = tail.end,
    tail.times = tail.times,
    tail.grid.size = tail.grid.size,
    min.tail.at.risk = min.tail.at.risk,
    min.tail.end.at.risk = min.tail.end.at.risk,
    min.events.before.tail = min.events.before.tail,
    lambda = lambda,
    min.gain = min.gain,
    ipcw.model = "marginal",
    validation.fraction = validation.fraction,
    early.stopping.rounds = early.stopping.rounds,
    stopping.tolerance = stopping.tolerance,
    posterior.epsilon = posterior.epsilon,
    seed = seed,
    na.action = na.action,
    keep.data = keep.data,
    verbose = verbose
  )
  object$call <- call
  object$target <- "cure_anchored_cox_latency"
  class(object) <- c("curecoxboost", "cureboost_joint")
  object
}
