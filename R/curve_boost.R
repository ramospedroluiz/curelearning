# Flexible finite-horizon refinement of a cure forest. The anchor and the
# correction sample are disjoint; validation is outside both fitted stages.
.ccrb_response <- function(formula, data) {
  mf <- stats::model.frame(formula, data, na.action = stats::na.fail)
  y <- stats::model.response(mf)
  if (!inherits(y, "Surv") || attr(y, "type") != "right")
    stop("A right-censored Surv response is required.")
  list(time = as.numeric(y[, 1]), status = as.integer(y[, 2]))
}

.ccrb_exposure <- function(time, breaks) {
  matrix(vapply(seq_len(length(breaks) - 1L), function(k)
    pmax(0, pmin(time, breaks[k + 1L]) - breaks[k]), numeric(length(time))),
    nrow = length(time), ncol = length(breaks) - 1L)
}

.ccrb_anchor <- function(anchor, data, breaks, hazard.floor) {
  p <- as.numeric(predict(anchor, data, type = "cure"))
  if (any(!is.finite(p)) || any(p >= 1))
    stop("The anchor must give finite probabilities below one; latency is undefined at one.")
  s <- predict(anchor, data, type = "survival", times = breaks[-1L])
  u <- pmin(pmax((s - p) / (1 - p), 1e-8), 1)
  u <- t(apply(cbind(1, u), 1L, cummin))
  cumulative <- -log(u)
  rate <- sweep(cumulative[, -1L, drop = FALSE] -
                  cumulative[, -ncol(cumulative), drop = FALSE],
                2, diff(breaks), "/")
  list(cure = p, rate = pmax(rate, hazard.floor), raw_survival = s)
}

.ccrb_design <- function(anchor, data, breaks) {
  x <- as.data.frame(.cfr_new_matrix(anchor, data))
  names(x) <- paste0("x", seq_len(ncol(x)))
  k <- length(breaks) - 1L
  z <- x[rep(seq_len(nrow(x)), times = k), , drop = FALSE]
  z$.interval <- rep((breaks[-1L] + breaks[-length(breaks)]) / 2,
                     each = nrow(x))
  rownames(z) <- NULL
  z
}

.ccrb_state <- function(eta, exposure, cure, status, event_interval) {
  rate <- exp(eta)
  H <- rowSums(rate * exposure)
  logu <- log1p(-cure) - H
  logp <- log(pmax(cure, .Machine$double.xmin))
  a <- pmax(logp, logu)
  logS <- a + log(exp(logp - a) + exp(logu - a))
  w <- exp(logu - logS)
  event <- status == 1L
  w[event] <- 1
  ll <- logS
  if (any(event)) ll[event] <- log1p(-cure[event]) +
    eta[cbind(which(event), event_interval[event])] - H[event]
  list(loss = -mean(ll), w = w, rate = rate)
}

.ccrb_curves <- function(rate, cure, breaks, times) {
  H <- rate %*% t(.ccrb_exposure(times, breaks))
  u <- exp(-H)
  list(survival = cure + (1 - cure) * u, susceptible = u)
}

.ccrb_km_weights <- function(train_y, valid_y, times, lower = 0.02) {
  km <- survival::survfit(survival::Surv(train_y$time, 1L-train_y$status) ~ 1)
  at <- function(t, before = FALSE) {
    if (before) t <- t - 1e-8 * pmax(1, abs(t))
    p <- findInterval(t, km$time)
    pmax(lower, c(1, km$surv)[p + 1L])
  }
  sapply(times, function(t) {
    ifelse(valid_y$time > t, 1/at(t),
           ifelse(valid_y$status == 1 & valid_y$time <= t,
                  1/at(valid_y$time, TRUE), 0))
  })
}

# One independent, subject-level paired comparison, after choosing ONE model.
# This is an approximate one-sided bound, not an exact finite-sample guarantee.
.ccrb_gate_decision <- function(improvement, alpha=0.05, minimum.gain=0) {
  .cfr_assert_scalar(alpha, "gate.alpha", 0, 1, inclusive=FALSE)
  .cfr_assert_scalar(minimum.gain, "minimum.gain", 0, Inf)
  n <- length(improvement)
  if (n < 2L || any(!is.finite(improvement))) stop("Invalid paired gate losses.")
  gain <- mean(improvement)
  se <- stats::sd(improvement)/sqrt(n)
  lower <- gain - stats::qt(1-alpha, df=n-1L)*se
  list(accepted=isTRUE(lower > minimum.gain), n=n, alpha=alpha,
       minimum_gain=minimum.gain, mean_improvement=gain,
       standard_error=se, lower_bound=lower)
}

.ccrb_gate <- function(object, gate.data, gate.weights, alpha, minimum.gain) {
  gy <- .ccrb_response(object$formula, gate.data)
  target <- outer(gy$time, object$score_times, ">")*1
  if (nrow(gate.data) < 30L || !identical(dim(gate.weights),dim(target)) ||
      any(!is.finite(gate.weights)) || any(gate.weights < 0))
    stop("Invalid independent gate sample or IPCW weights.")
  zero <- object
  zero$trees <- list()
  zero$nrounds <- 0L
  before <- predict(zero,gate.data,times=object$score_times)
  after <- predict(object,gate.data,times=object$score_times)
  # Average within subject before estimating variability across subjects.
  improvement <- rowMeans(gate.weights*((target-before)^2-(target-after)^2))
  decision <- .ccrb_gate_decision(improvement,alpha,minimum.gain)
  decision$candidate_rounds <- object$nrounds
  decision$reason <- if(object$nrounds==0L) "no_candidate" else
    if(decision$accepted) "positive_independent_lower_bound" else "insufficient_evidence"
  object$candidate_nrounds <- object$nrounds
  object$candidate_validation_brier <- object$validation_brier
  if (!decision$accepted) {
    object$trees <- list()
    object$nrounds <- 0L
    object$validation_brier <- object$history$validation_brier[1L]
  }
  object$gate <- decision
  object$gate_n <- nrow(gate.data)
  object$selection <- "gated"
  object
}

# Internal engine also used by the benchmark to share exactly the same anchor.
.ccrb_fit <- function(anchor, formula, data, validation.data, breaks,
                      score.times, score.weights = NULL, nrounds = 150L,
                      learning.rate = 0.1, maxdepth = 2L, minbucket = 30L,
                      patience = 20L, hazard.floor = 1e-5, verbose = FALSE,
                      selection=c("minimum","gated"), gate.data=NULL,
                      gate.weights=NULL, gate.alpha=0.05, minimum.gain=0) {
  selection <- match.arg(selection)
  if (selection=="gated" && is.null(gate.data))
    stop("Gated selection requires a separate, untouched gate.data sample.")
  if (!requireNamespace("rpart", quietly = TRUE)) stop("Install rpart.")
  if (length(breaks) < 3L || breaks[1] != 0 ||
      any(!is.finite(breaks)) || any(diff(breaks) <= 0)) stop("Invalid breaks.")
  if (any(score.times <= 0 | score.times > max(breaks))) stop("Invalid scoring times.")
  .cfr_assert_scalar(nrounds, "nrounds", 0, 10000, integer = TRUE)
  .cfr_assert_scalar(learning.rate, "learning.rate", 0, 1, inclusive = FALSE)
  .cfr_assert_scalar(hazard.floor, "hazard.floor", 0, Inf, inclusive = FALSE)
  y <- .ccrb_response(formula, data)
  vy <- .ccrb_response(formula, validation.data)
  tau <- max(breaks)
  if (any(y$time <= 0 & y$status == 1)) stop("Event times must be positive.")
  status <- y$status * (y$time <= tau)
  event_interval <- pmin(length(breaks)-1L,
                        pmax(1L, findInterval(y$time, breaks, left.open=TRUE)))
  exposure <- .ccrb_exposure(y$time, breaks)
  base <- .ccrb_anchor(anchor, data, breaks, hazard.floor)
  vb <- .ccrb_anchor(anchor, validation.data, breaks, hazard.floor)
  design <- .ccrb_design(anchor, data, breaks)
  vd <- .ccrb_design(anchor, validation.data, breaks)
  eta <- log(base$rate)
  veta <- log(vb$rate)
  observed <- exposure > 0
  counts <- matrix(0, nrow(data), ncol(exposure))
  counts[cbind(which(status == 1), event_interval[status == 1])] <- 1
  if (is.null(score.weights)) score.weights <- .ccrb_km_weights(y, vy, score.times)
  target <- outer(vy$time, score.times, ">") * 1
  if (!identical(dim(score.weights), dim(target)) || any(!is.finite(score.weights)))
    stop("Invalid validation IPCW weights.")
  score <- function(e) mean(score.weights *
    (target - .ccrb_curves(exp(e), vb$cure, breaks, score.times)$survival)^2)
  state <- .ccrb_state(eta, exposure, base$cure, status, event_interval)
  best <- score(veta)
  history <- data.frame(iteration=0L, train_nll=state$loss, validation_brier=best,
                        step=0)
  trees <- list()
  selected <- 0L
  stale <- 0L
  for (m in seq_len(nrounds)) {
    expected <- state$w * state$rate * exposure
    gradient <- counts - expected
    # Positive EM-surrogate working weights on observed interval exposure.
    # No tail-maturity filters are applied to latency splits.
    hessian <- pmax(expected[observed], 1e-5)
    dd <- design[as.vector(observed), , drop=FALSE]
    dd$.working <- gradient[observed] / hessian
    dd$.weight <- hessian
    tree <- rpart::rpart(stats::reformulate(names(design), response=".working"),
                        data=dd, weights=dd$.weight,
                        method="anova", control=rpart::rpart.control(
                          cp=0, maxdepth=maxdepth, minbucket=minbucket,
                          minsplit=2L*minbucket, xval=0, maxsurrogate=0))
    # Prediction needs only design columns; do not serialize training closures.
    attr(tree$terms, ".Environment") <- baseenv()
    direction <- matrix(pmax(-3, pmin(3, as.numeric(predict(tree, design)))),
                        nrow(data))
    step <- learning.rate
    accepted <- FALSE
    for (ls in seq_len(16L)) {
      proposed <- eta + step * direction
      next_state <- .ccrb_state(proposed, exposure, base$cure, status, event_interval)
      if (is.finite(next_state$loss) && next_state$loss < state$loss - 1e-10) {
        accepted <- TRUE
        break
      }
      step <- step / 2
    }
    if (!accepted) break
    eta <- proposed
    state <- next_state
    veta <- veta + step * matrix(pmax(-3, pmin(3, as.numeric(predict(tree, vd)))),
                                  nrow(validation.data))
    trees[[m]] <- list(tree=tree, step=step)
    vscore <- score(veta)
    history <- rbind(history, data.frame(iteration=m, train_nll=state$loss,
                                        validation_brier=vscore, step=step))
    if (vscore < best - 1e-7) {
      best <- vscore
      selected <- m
      stale <- 0L
    } else stale <- stale + 1L
    if (verbose && m %in% unique(round(nrounds*c(.25,.5,.75,1))))
      message(round(100*m/nrounds), "% of flexible correction rounds")
    if (stale >= patience) break
  }
  object <- list(anchor=anchor, formula=formula, breaks=breaks, horizon=tau,
                 hazard_floor=hazard.floor, trees=trees[seq_len(selected)],
                 nrounds=selected, history=history, score_times=score.times,
                 validation_brier=best, correction_n=nrow(data),
                 validation_n=nrow(validation.data),
                 selection="minimum", gate_n=0L,
                 target="finite_horizon_curve_with_fixed_tail_anchor")
  class(object) <- "curecurveboost"
  if (selection=="gated") {
    if (is.null(gate.weights)) gate.weights <- .ccrb_km_weights(
      y,.ccrb_response(formula,gate.data),score.times)
    object <- .ccrb_gate(object,gate.data,gate.weights,gate.alpha,minimum.gain)
  }
  object
}

#' Refine a CureForest curve using flexible time-dependent hazard boosting
#' @export
curecurveboost <- function(formula, data, breaks, nrounds=150L,
                           learning.rate=0.1, maxdepth=2L, minbucket=30L,
                           validation.fraction=0.2, anchor.fraction=0.5,
                           forest.control=list(), seed=1L, verbose=FALSE,
                           selection=c("gated","minimum"), gate.fraction=0.5,
                           gate.alpha=0.05, minimum.gain=0) {
  selection <- match.arg(selection)
  .cfr_assert_scalar(gate.fraction,"gate.fraction",0,1,inclusive=FALSE)
  .cfr_assert_scalar(gate.alpha,"gate.alpha",0,1,inclusive=FALSE)
  .cfr_assert_scalar(minimum.gain,"minimum.gain",0,Inf)
  .cfr_assert_scalar(validation.fraction, "validation.fraction", 0, .5, FALSE, FALSE)
  .cfr_assert_scalar(anchor.fraction, "anchor.fraction", 0, 1-validation.fraction,
                     FALSE, FALSE)
  n <- nrow(data)
  set.seed(seed)
  perm <- sample.int(n)
  ia <- perm[seq_len(floor(n*anchor.fraction))]
  iv <- utils::tail(perm, floor(n*validation.fraction))
  ib <- setdiff(perm, c(ia, iv))
  ig <- integer(0)
  if (selection=="gated") {
    ng <- floor(length(iv)*gate.fraction)
    ig <- utils::tail(iv,ng)
    iv <- utils::head(iv,length(iv)-ng)
    if(length(ig)<30L) stop("Insufficient independent gate sample.")
  }
  if (min(length(ia),length(iv),length(ib)) < 30L) stop("Insufficient split samples.")
  args <- utils::modifyList(list(ntree=300L, n.cores=1L, seed=seed,
                                 honesty=TRUE, oob=FALSE, estimator="km"), forest.control)
  if (any(c("formula","data") %in% names(args))) stop("Do not override formula/data.")
  anchor <- do.call(randomforestcure, c(list(formula=formula, data=data[ia,,drop=FALSE]), args))
  object <- .ccrb_fit(anchor, formula, data[ib,,drop=FALSE], data[iv,,drop=FALSE],
                      breaks, breaks[-1L], nrounds=nrounds,
                      learning.rate=learning.rate, maxdepth=maxdepth,
                      minbucket=minbucket, verbose=verbose, selection=selection,
                      gate.data=if(selection=="gated")data[ig,,drop=FALSE] else NULL,
                      gate.alpha=gate.alpha,minimum.gain=minimum.gain)
  object$split <- list(anchor=ia, correction=ib, validation=iv, gate=ig)
  object$seed <- seed
  object
}

#' Predict a finite-horizon corrected survival curve
#' @export
predict.curecurveboost <- function(object, newdata,
                                   type=c("survival","susceptible","cure","all"),
                                   times=NULL, ...) {
  type <- match.arg(type)
  if (type == "cure") return(as.numeric(predict(object$anchor,newdata,type="cure")))
  if (is.null(times)) times <- object$breaks[-1L]
  if (!length(times) || any(!is.finite(times)) || any(times < 0 | times > object$horizon))
    stop("Prediction times must lie in [0, horizon]; extrapolation is not identified.")
  base <- .ccrb_anchor(object$anchor, newdata, object$breaks, object$hazard_floor)
  eta <- log(base$rate)
  design <- .ccrb_design(object$anchor, newdata, object$breaks)
  for (b in object$trees) eta <- eta + b$step * matrix(
    pmax(-3,pmin(3,as.numeric(predict(b$tree, design)))), nrow(newdata))
  curves <- .ccrb_curves(exp(eta), base$cure, object$breaks, times)
  colnames(curves$survival) <- colnames(curves$susceptible) <- as.character(times)
  if (type != "all") return(curves[[type]])
  c(curves, list(cure=base$cure, times=times, horizon=object$horizon))
}
