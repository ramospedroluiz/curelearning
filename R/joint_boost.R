.cbj_baseline <- function(time, status, posterior_uncured, latency_raw,
                          index, model = c("weibull", "breslow")) {
  model <- match.arg(model)
  time <- time[index]
  status <- status[index]
  posterior_uncured <- posterior_uncured[index]
  latency_raw <- latency_raw[index]
  if (model == "weibull") {
    event_count <- sum(status)
    if (event_count <= 0) {
      stop("The Weibull latency baseline requires at least one event.",
           call. = FALSE)
    }
    safe_time <- pmax(time, .Machine$double.eps)
    clipped_latency <- pmin(latency_raw, 20)
    objective <- function(log_shape) {
      shape <- exp(log_shape)
      log_terms <- log(pmax(posterior_uncured, .Machine$double.eps)) +
        clipped_latency + shape * log(safe_time)
      maximum <- max(log_terms)
      log_denominator <- maximum + log(sum(exp(log_terms - maximum)))
      log_scale <- log(event_count) - log_denominator
      event_part <- sum(status * (
        log_shape + log_scale + (shape - 1) * log(safe_time) +
          clipped_latency
      ))
      -(event_part - event_count)
    }
    optimum <- stats::optimize(
      objective, interval = log(c(0.10, 10)), tol = 1e-7
    )
    shape <- exp(optimum$minimum)
    log_terms <- log(pmax(posterior_uncured, .Machine$double.eps)) +
      clipped_latency + shape * log(safe_time)
    maximum <- max(log_terms)
    log_denominator <- maximum + log(sum(exp(log_terms - maximum)))
    return(list(
      type = "weibull",
      shape = shape,
      log_scale = log(event_count) - log_denominator
    ))
  }
  relative_risk <- posterior_uncured * exp(pmin(latency_raw, 20))
  tab <- data.frame(
    time = time,
    risk = relative_risk,
    events = status,
    stringsAsFactors = FALSE
  )
  aggregate <- stats::aggregate(
    cbind(risk, events) ~ time, data = tab, FUN = sum
  )
  aggregate <- aggregate[order(aggregate$time), , drop = FALSE]
  aggregate$risk_set <- rev(cumsum(rev(aggregate$risk)))
  event_rows <- aggregate$events > 0
  event_table <- aggregate[event_rows, , drop = FALSE]
  event_table$increment <- event_table$events /
    pmax(event_table$risk_set, .Machine$double.eps)
  event_table$cumulative_hazard <- cumsum(event_table$increment)
  list(
    type = "breslow",
    table = event_table[, c(
      "time", "increment", "cumulative_hazard"
    ), drop = FALSE]
  )
}

.cbj_baseline_at <- function(baseline, query, increment = FALSE) {
  # Data frames are accepted for backward compatibility with early
  # experimental cureboost_joint objects.
  if (is.data.frame(baseline)) {
    table <- baseline
    type <- "breslow"
  } else {
    type <- baseline$type
    table <- baseline$table
  }
  if (identical(type, "weibull")) {
    query <- pmax(query, 0)
    if (increment) {
      return(exp(baseline$log_scale) * baseline$shape *
               pmax(query, .Machine$double.eps)^(baseline$shape - 1))
    }
    return(exp(baseline$log_scale) * query^baseline$shape)
  }
  if (is.null(table) || !nrow(table)) return(numeric(length(query)))
  position <- findInterval(query, table$time)
  value <- numeric(length(query))
  use <- position > 0L
  column <- if (increment) "increment" else "cumulative_hazard"
  value[use] <- table[[column]][position[use]]
  value
}

.cbj_log_hazard_at <- function(baseline, query, epsilon) {
  if (!is.data.frame(baseline) && identical(baseline$type, "weibull")) {
    return(
      baseline$log_scale + log(baseline$shape) +
        (baseline$shape - 1) * log(pmax(query, .Machine$double.eps))
    )
  }
  log(pmax(.cbj_baseline_at(baseline, query, increment = TRUE), epsilon))
}

.cbj_posterior <- function(time, status, cure_probability, latency_raw,
                           baseline, epsilon) {
  cumulative_hazard <- .cbj_baseline_at(baseline, time)
  susceptible_survival <- exp(
    -cumulative_hazard * exp(pmin(latency_raw, 20))
  )
  susceptible_probability <- 1 - cure_probability
  denominator <- cure_probability +
    susceptible_probability * susceptible_survival
  posterior <- susceptible_probability * susceptible_survival /
    pmax(denominator, epsilon)
  posterior[status == 1L] <- 1
  .cfr_clip(posterior, epsilon, 1 - epsilon)
}

.cbj_observed_loss <- function(time, status, cure_probability, latency_raw,
                               baseline, index, epsilon) {
  time <- time[index]
  status <- status[index]
  cure_probability <- cure_probability[index]
  latency_raw <- latency_raw[index]
  cumulative_hazard <- .cbj_baseline_at(baseline, time)
  susceptible_survival <- exp(
    -cumulative_hazard * exp(pmin(latency_raw, 20))
  )
  log_likelihood <- numeric(length(index))
  event <- status == 1L
  if (any(event)) {
    log_likelihood[event] <-
      log(pmax(1 - cure_probability[event], epsilon)) +
      .cbj_log_hazard_at(baseline, time[event], epsilon) +
      latency_raw[event] -
      cumulative_hazard[event] * exp(pmin(latency_raw[event], 20))
  }
  if (any(!event)) {
    marginal_survival <- cure_probability[!event] +
      (1 - cure_probability[!event]) * susceptible_survival[!event]
    log_likelihood[!event] <- log(pmax(marginal_survival, epsilon))
  }
  -mean(log_likelihood)
}

.cbj_add_tree <- function(index, residual, context, learning_rate,
                          component, iteration) {
  tree <- .cb_grow_tree(index, residual, context)
  direction <- .cb_predict_tree(
    tree, context$x[index, , drop = FALSE]
  )
  if (any(!is.finite(direction)) ||
      max(abs(direction)) <= .Machine$double.eps) {
    return(NULL)
  }
  tree$coefficient <- learning_rate
  tree$component <- component
  tree$iteration <- iteration
  tree$gradient_correlation <- mean(residual * direction)
  tree
}

.cbj_refine_posterior <- function(time, status, cure_probability, latency_raw,
                                  train, posterior, epsilon,
                                  baseline_model,
                                  iterations = 20L) {
  for (iteration in seq_len(iterations)) {
    baseline <- .cbj_baseline(
      time, status, posterior, latency_raw, train, baseline_model
    )
    updated <- .cbj_posterior(
      time, status, cure_probability, latency_raw, baseline, epsilon
    )
    if (max(abs(updated[train] - posterior[train])) < 1e-7) {
      posterior <- updated
      break
    }
    posterior <- updated
  }
  baseline <- .cbj_baseline(
    time, status, posterior, latency_raw, train, baseline_model
  )
  list(posterior = posterior, baseline = baseline)
}

#' Fit joint incidence-latency cure boosting
#'
#' Fits an experimental mixture-cure boosting model with separate additive
#' predictors for cure incidence and susceptible latency. An EM iteration
#' alternates posterior susceptibility weights, logistic incidence updates,
#' and proportional-hazards latency updates. The resulting marginal survival
#' curve is coherent by construction.
#'
#' @param formula A formula with a `survival::Surv(time, status)` response.
#' @param data A data frame.
#' @param em.iterations Maximum number of EM boosting iterations.
#' @param incidence.steps Number of incidence trees per EM iteration.
#' @param latency.steps Number of latency trees per EM iteration.
#' @param initialization.rounds Tail-pseudo-outcome incidence trees used only
#'   to initialize the EM algorithm.
#' @param learning.rate.incidence Shrinkage for incidence trees.
#' @param learning.rate.latency Shrinkage for latency trees.
#' @param latency.baseline Baseline susceptible latency model. `"weibull"`
#'   supplies the structural tail constraint needed to distinguish cure from
#'   very slow failures; `"breslow"` provides a semiparametric sensitivity
#'   analysis.
#' @param incidence.posterior.weight Weight assigned to the EM posterior cure
#'   target. The default of one gives the likelihood-based EM update. Values
#'   below one are retained only for experimental tail-target sensitivity.
#' @param incidence.mode Whether incidence is updated jointly by EM or held at
#'   an externally estimated nonparametric cure anchor.
#' @param incidence.fit A fitted `cureboost`, `survcure`, or
#'   `randomforestcure` object supplying the cure anchor when
#'   `incidence.mode = "anchored"`.
#' @param validation.fraction Fraction held out for early stopping.
#' @param early.stopping.rounds Non-improving EM iterations before stopping.
#' @param maxdepth Maximum depth of component trees.
#' @param nodesize Minimum node size.
#' @param mtry Predictors considered at each split.
#' @param nsplit Candidate cut points per numeric predictor.
#' @param tail.start Start of the mature follow-up window.
#' @param tail.quantile Quantile used when `tail.start` is omitted.
#' @param tail.end End of the mature follow-up window.
#' @param tail.times Optional tail evaluation grid.
#' @param tail.grid.size Size of an automatically generated tail grid.
#' @param min.tail.at.risk Minimum number at risk at `tail.start`.
#' @param min.tail.end.at.risk Minimum number at risk at `tail.end`.
#' @param min.events.before.tail Minimum events before `tail.start`.
#' @param lambda Split regularization constant.
#' @param min.gain Minimum split gain.
#' @param ipcw.model Censoring model used for the optional tail target.
#' @param g.fun User censoring-survival function when requested.
#' @param ipcw.folds Cross-fitting folds for the censoring model.
#' @param g.lower Lower clipping bound for censoring survival.
#' @param stopping.tolerance Minimum validation improvement.
#' @param posterior.epsilon Numerical clipping constant.
#' @param seed Reproducibility seed.
#' @param na.action Missing-data action.
#' @param keep.data Whether to retain training data.
#' @param verbose Whether to report progress at 25-percent increments.
#' @return An object of class `cureboost_joint`.
#' @export
cureboost_joint <- function(formula,
                            data,
                            em.iterations = 50L,
                            incidence.steps = 2L,
                            latency.steps = 2L,
                            initialization.rounds = 20L,
                            learning.rate.incidence = 0.05,
                            learning.rate.latency = 0.05,
                            latency.baseline = c("weibull", "breslow"),
                            incidence.posterior.weight = 1,
                            incidence.mode = c("joint", "anchored"),
                            incidence.fit = NULL,
                            maxdepth = 1L,
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
                            ipcw.model = c("cox", "marginal", "user"),
                            g.fun = NULL,
                            ipcw.folds = 5L,
                            g.lower = 0.02,
                            validation.fraction = 0.20,
                            early.stopping.rounds = 8L,
                            stopping.tolerance = 1e-6,
                            posterior.epsilon = 1e-5,
                            seed = 1L,
                            na.action = stats::na.omit,
                            keep.data = FALSE,
                            verbose = FALSE) {
  call <- match.call()
  latency.baseline <- match.arg(latency.baseline)
  incidence.mode <- match.arg(incidence.mode)
  .cfr_assert_scalar(em.iterations, "em.iterations", 1, 10000, integer = TRUE)
  .cfr_assert_scalar(incidence.steps, "incidence.steps", 0, 100, integer = TRUE)
  if (incidence.mode == "joint" && incidence.steps < 1L) {
    stop("incidence.steps must be positive for joint incidence updating.",
         call. = FALSE)
  }
  anchor_classes <- c("cureboost", "survcure", "randomforestcure")
  if (incidence.mode == "anchored" &&
      !inherits(incidence.fit, anchor_classes)) {
    stop(paste0(
      "incidence.fit must be a fitted cureboost, survcure, or ",
      "randomforestcure object in anchored mode."
    ), call. = FALSE)
  }
  .cfr_assert_scalar(latency.steps, "latency.steps", 1, 100, integer = TRUE)
  .cfr_assert_scalar(initialization.rounds, "initialization.rounds", 0, 10000,
                     integer = TRUE)
  .cfr_assert_scalar(learning.rate.incidence, "learning.rate.incidence", 0, 1)
  .cfr_assert_scalar(learning.rate.latency, "learning.rate.latency", 0, 1)
  if (learning.rate.incidence <= 0 || learning.rate.latency <= 0) {
    stop("Both learning rates must be strictly positive.", call. = FALSE)
  }
  .cfr_assert_scalar(
    incidence.posterior.weight, "incidence.posterior.weight", 0, 1
  )
  .cfr_assert_scalar(validation.fraction, "validation.fraction", 0, 0.5)
  .cfr_assert_scalar(early.stopping.rounds, "early.stopping.rounds", 1,
                     10000, integer = TRUE)
  .cfr_assert_scalar(stopping.tolerance, "stopping.tolerance", 0, Inf)
  .cfr_assert_scalar(posterior.epsilon, "posterior.epsilon", 1e-12, 0.05)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }

  fit <- .cfr_prepare_fit(
    formula = formula,
    data = data,
    na.action = na.action,
    splitrule = "cure",
    estimator = "ipcw",
    tail.start = tail.start,
    tail.quantile = tail.quantile,
    tail.end = tail.end,
    tail.times = tail.times,
    tail.grid.size = tail.grid.size,
    curve.grid.size = max(20L, as.integer(tail.grid.size)),
    nodesize = nodesize,
    min.tail.at.risk = min.tail.at.risk,
    min.tail.end.at.risk = min.tail.end.at.risk,
    min.events.before.tail = min.events.before.tail,
    maxdepth = maxdepth,
    mtry = mtry,
    nsplit = nsplit,
    epsilon = 1e-4,
    lambda = lambda,
    min.gain = min.gain,
    ipcw.model = ipcw.model,
    g.fun = g.fun,
    ipcw.folds = ipcw.folds,
    g.lower = g.lower,
    seed = seed
  )
  context <- fit$context
  n <- fit$prepared$n
  split <- .cb_validation_split(n, validation.fraction, seed, context)
  train <- split$train
  validation <- split$validation
  epsilon <- posterior.epsilon
  pseudo_cure <- .cfr_clip(context$ipcw_tail_mean, epsilon, 1 - epsilon)

  anchor_cure <- NULL
  if (incidence.mode == "anchored") {
    analysis_data <- fit$prepared$data
    anchor_cure <- .cfr_clip(
      as.numeric(predict(incidence.fit, analysis_data, type = "cure")),
      epsilon, 1 - epsilon
    )
    if (length(anchor_cure) != n || any(!is.finite(anchor_cure))) {
      stop("The incidence anchor did not produce one finite prediction per row.",
           call. = FALSE)
    }
  }
  initial_cure <- if (is.null(anchor_cure)) {
    .cfr_clip(mean(pseudo_cure[train]), 0.02, 0.98)
  } else {
    mean(anchor_cure[train])
  }
  incidence_initial <- if (is.null(anchor_cure)) {
    stats::qlogis(initial_cure)
  } else 0
  latency_initial <- 0
  incidence_raw <- if (is.null(anchor_cure)) {
    rep(incidence_initial, n)
  } else stats::qlogis(anchor_cure)
  latency_raw <- rep(latency_initial, n)
  incidence_trees <- list()
  latency_trees <- list()

  set.seed(seed)
  if (is.null(anchor_cure) && initialization.rounds > 0L) {
    for (round in seq_len(initialization.rounds)) {
      residual <- pseudo_cure[train] - stats::plogis(incidence_raw[train])
      tree <- .cbj_add_tree(
        train, residual, context, learning.rate.incidence,
        "incidence_initialization", 0L
      )
      if (is.null(tree)) break
      incidence_trees[[length(incidence_trees) + 1L]] <- tree
      incidence_raw <- incidence_raw + tree$coefficient *
        .cb_predict_tree(tree, context$x)
    }
  }

  cure_probability <- stats::plogis(incidence_raw)
  posterior <- ifelse(context$status == 1L, 1, 1 - cure_probability)
  refined <- .cbj_refine_posterior(
    context$time, context$status, cure_probability, latency_raw,
    train, posterior, epsilon, latency.baseline
  )
  posterior <- refined$posterior
  baseline <- refined$baseline

  history <- vector("list", em.iterations)
  best_validation <- Inf
  best_iteration <- 0L
  best_incidence_count <- length(incidence_trees)
  best_latency_count <- 0L
  patience <- 0L
  stop_reason <- "maximum EM iterations reached"
  progress_points <- unique(pmax(
    1L, round(em.iterations * c(0.25, 0.50, 0.75, 1))
  ))
  start <- proc.time()[[3L]]

  for (iteration in seq_len(em.iterations)) {
    cure_probability <- stats::plogis(incidence_raw)
    refined <- .cbj_refine_posterior(
      context$time, context$status, cure_probability, latency_raw,
      train, posterior, epsilon, latency.baseline, iterations = 5L
    )
    posterior <- refined$posterior
    baseline <- refined$baseline

    tail_hazard <- .cbj_baseline_at(baseline, context$tail_grid)
    susceptible_tail <- rowMeans(exp(-outer(
      exp(pmin(latency_raw, 20)), tail_hazard
    )))
    decontaminated_cure <-
      (pseudo_cure - susceptible_tail) /
      pmax(1 - susceptible_tail, epsilon)
    decontaminated_cure <- .cfr_clip(
      decontaminated_cure, epsilon, 1 - epsilon
    )
    incidence_target <- incidence.posterior.weight * (1 - posterior) +
      (1 - incidence.posterior.weight) * decontaminated_cure

    if (incidence.mode == "joint") {
      for (step in seq_len(incidence.steps)) {
        residual <- incidence_target[train] -
          stats::plogis(incidence_raw[train])
        tree <- .cbj_add_tree(
          train, residual, context, learning.rate.incidence,
          "incidence", iteration
        )
        if (is.null(tree)) break
        incidence_trees[[length(incidence_trees) + 1L]] <- tree
        incidence_raw <- incidence_raw + tree$coefficient *
          .cb_predict_tree(tree, context$x)
      }
    }

    for (step in seq_len(latency.steps)) {
      baseline <- .cbj_baseline(
        context$time, context$status, posterior, latency_raw, train,
        latency.baseline
      )
      cumulative_hazard <- .cbj_baseline_at(
        baseline, context$time[train]
      )
      residual <- context$status[train] - posterior[train] *
        exp(pmin(latency_raw[train], 20)) * cumulative_hazard
      residual <- pmin(2, pmax(-2, residual))
      tree <- .cbj_add_tree(
        train, residual, context, learning.rate.latency,
        "latency", iteration
      )
      if (is.null(tree)) break
      latency_trees[[length(latency_trees) + 1L]] <- tree
      latency_raw <- latency_raw + tree$coefficient *
        .cb_predict_tree(tree, context$x)
      latency_raw <- pmin(10, pmax(-10, latency_raw))
    }

    cure_probability <- stats::plogis(incidence_raw)
    refined <- .cbj_refine_posterior(
      context$time, context$status, cure_probability, latency_raw,
      train, posterior, epsilon, latency.baseline, iterations = 5L
    )
    posterior <- refined$posterior
    baseline <- refined$baseline
    training_loss <- .cbj_observed_loss(
      context$time, context$status, cure_probability, latency_raw,
      baseline, train, epsilon
    )
    validation_loss <- if (length(validation)) {
      .cbj_observed_loss(
        context$time, context$status, cure_probability, latency_raw,
        baseline, validation, epsilon
      )
    } else training_loss

    history[[iteration]] <- data.frame(
      iteration = iteration,
      training_loss = training_loss,
      validation_loss = validation_loss,
      mean_cure = mean(cure_probability[train]),
      mean_posterior_uncured = mean(posterior[train]),
      mean_incidence_target = mean(incidence_target[train]),
      mean_decontaminated_cure = mean(decontaminated_cure[train]),
      incidence_trees = length(incidence_trees),
      latency_trees = length(latency_trees),
      row.names = NULL
    )
    if (validation_loss < best_validation - stopping.tolerance) {
      best_validation <- validation_loss
      best_iteration <- iteration
      best_incidence_count <- length(incidence_trees)
      best_latency_count <- length(latency_trees)
      patience <- 0L
    } else {
      patience <- patience + 1L
    }
    if (verbose && iteration %in% progress_points) {
      message(round(100 * iteration / em.iterations),
              "% of requested joint boosting iterations completed")
    }
    if (length(validation) && patience >= early.stopping.rounds) {
      stop_reason <- "validation early stopping"
      break
    }
  }

  completed_iterations <- sum(vapply(history, Negate(is.null), logical(1)))
  history <- do.call(rbind, history[seq_len(completed_iterations)])
  if (!length(validation)) {
    best_iteration <- completed_iterations
    best_incidence_count <- length(incidence_trees)
    best_latency_count <- length(latency_trees)
  }
  incidence_trees <- if (best_incidence_count) {
    incidence_trees[seq_len(best_incidence_count)]
  } else list()
  latency_trees <- if (best_latency_count) {
    latency_trees[seq_len(best_latency_count)]
  } else list()
  history$selected <- history$iteration <= best_iteration

  incidence_raw <- if (is.null(anchor_cure)) {
    .cb_fitted_raw(incidence_initial, incidence_trees, context$x)
  } else stats::qlogis(anchor_cure)
  latency_raw <- .cb_fitted_raw(
    latency_initial, latency_trees, context$x
  )
  cure_probability <- stats::plogis(incidence_raw)
  posterior <- ifelse(context$status == 1L, 1, 1 - cure_probability)
  refined <- .cbj_refine_posterior(
    context$time, context$status, cure_probability, latency_raw,
    train, posterior, epsilon, latency.baseline
  )
  posterior <- refined$posterior
  baseline <- refined$baseline

  metadata <- .cfr_fit_metadata(
    fit$prepared, context, call, keep.data
  )
  incidence_importance <- if (is.null(incidence.fit)) {
    .cb_build_importance(incidence_trees, fit$prepared$variable_names)
  } else if (inherits(incidence.fit, "cureboost")) {
    cureboost_importance(incidence.fit, aggregate = FALSE)
  } else if (inherits(incidence.fit, "randomforestcure")) {
    cure_importance(incidence.fit, aggregate = FALSE)
  } else {
    .cfr_build_importance(list(incidence.fit), fit$prepared$variable_names)
  }
  latency_importance <- .cb_build_importance(
    latency_trees, fit$prepared$variable_names
  )
  object <- c(list(
    incidence_trees = incidence_trees,
    latency_trees = latency_trees,
    incidence_initial = incidence_initial,
    latency_initial = latency_initial,
    cure_probability = cure_probability,
    incidence_raw = incidence_raw,
    latency_raw = latency_raw,
    posterior_uncured = posterior,
    baseline = baseline,
    train_index = train,
    validation_index = validation,
    em_iterations = best_iteration,
    em_iterations_completed = completed_iterations,
    em_iterations_requested = as.integer(em.iterations),
    initialization_rounds = as.integer(initialization.rounds),
    learning_rate_incidence = learning.rate.incidence,
    learning_rate_latency = learning.rate.latency,
    incidence_posterior_weight = incidence.posterior.weight,
    incidence_mode = incidence.mode,
    incidence_anchor = incidence.fit,
    latency_baseline = latency.baseline,
    history = history,
    stop_reason = stop_reason,
    incidence_importance = incidence_importance,
    latency_importance = latency_importance,
    target = "mixture_cure_incidence_latency",
    elapsed = proc.time()[[3L]] - start,
    posterior_epsilon = epsilon
  ), metadata)
  class(object) <- "cureboost_joint"
  object
}

#' Predict from joint incidence-latency cure boosting
#'
#' @param object A fitted `cureboost_joint` object.
#' @param newdata New predictor data.
#' @param type Cure probability, susceptible survival, marginal survival,
#'   latency risk, or all components.
#' @param times Times for survival predictions.
#' @param ... Unused additional arguments.
#' @export
predict.cureboost_joint <- function(object,
                                    newdata = NULL,
                                    type = c("cure", "survival", "susceptible",
                                             "latency", "all"),
                                    times = NULL,
                                    ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    if (is.null(object$training_data)) {
      if (type == "cure") return(object$cure_probability)
      if (type == "latency") return(object$latency_raw)
      stop("newdata is required because training data were not retained.",
           call. = FALSE)
    }
    newdata <- object$training_data
  }
  x <- .cfr_new_matrix(object, newdata)
  incidence_raw <- if (!is.null(object$incidence_anchor)) {
    stats::qlogis(.cfr_clip(as.numeric(predict(
      object$incidence_anchor, newdata, type = "cure"
    )), object$posterior_epsilon, 1 - object$posterior_epsilon))
  } else {
    .cb_fitted_raw(object$incidence_initial, object$incidence_trees, x)
  }
  latency_raw <- .cb_fitted_raw(
    object$latency_initial, object$latency_trees, x
  )
  cure <- stats::plogis(incidence_raw)
  names(cure) <- names(latency_raw) <- rownames(newdata)
  if (type == "cure") return(cure)
  if (type == "latency") return(latency_raw)
  if (is.null(times) || !length(times) || any(!is.finite(times)) ||
      any(times < 0)) {
    stop("times must contain finite nonnegative values.", call. = FALSE)
  }
  times <- as.numeric(times)
  cumulative_hazard <- .cbj_baseline_at(object$baseline, times)
  susceptible <- exp(-outer(
    exp(pmin(latency_raw, 20)), cumulative_hazard
  ))
  survival <- cure + (1 - cure) * susceptible
  colnames(susceptible) <- colnames(survival) <- format(
    times, trim = TRUE, scientific = FALSE
  )
  rownames(susceptible) <- rownames(survival) <- rownames(newdata)
  if (type == "susceptible") return(susceptible)
  if (type == "survival") return(survival)
  list(
    cure = cure,
    susceptible = susceptible,
    survival = survival,
    latency = latency_raw,
    incidence_raw = incidence_raw,
    times = times
  )
}

#' Variable importance for joint cure boosting
#'
#' @param object A fitted `cureboost_joint` object.
#' @param component Incidence, latency, or both components.
#' @param aggregate Aggregate encoded columns to formula terms.
#' @export
cureboost_joint_importance <- function(
    object, component = c("both", "incidence", "latency"), aggregate = TRUE) {
  if (!inherits(object, "cureboost_joint")) {
    stop("object must inherit from cureboost_joint.", call. = FALSE)
  }
  component <- match.arg(component)
  prepare <- function(value, label) {
    if (isTRUE(aggregate)) {
      value <- .cfr_aggregate_importance(value, object$variable_map)
    }
    value$component <- label
    value
  }
  out <- switch(
    component,
    incidence = prepare(object$incidence_importance, "incidence"),
    latency = prepare(object$latency_importance, "latency"),
    both = rbind(
      prepare(object$incidence_importance, "incidence"),
      prepare(object$latency_importance, "latency")
    )
  )
  out[order(out$component, -out$total_gain), , drop = FALSE]
}

#' @export
print.cureboost_joint <- function(x, ...) {
  cat("Joint incidence-latency cure boosting\n")
  cat("  Observations / encoded predictors:", x$n, "/", x$p, "\n")
  cat("  Selected / completed / requested EM iterations:",
      x$em_iterations, "/", x$em_iterations_completed, "/",
      x$em_iterations_requested, "\n")
  cat("  Incidence / latency trees:", length(x$incidence_trees), "/",
      length(x$latency_trees), "\n")
  if (!is.null(x$incidence_anchor)) {
    anchor_trees <- if (inherits(x$incidence_anchor, "cureboost")) {
      x$incidence_anchor$nrounds
    } else if (inherits(x$incidence_anchor, "randomforestcure")) {
      length(x$incidence_anchor$trees)
    } else 1L
    cat("  Anchored incidence trees:", anchor_trees, "\n")
  }
  cat("  Incidence mode:", x$incidence_mode, "\n")
  cat("  Susceptible latency baseline:", x$latency_baseline, "\n")
  cat("  Mean fitted cure probability:",
      format(mean(x$cure_probability), digits = 5), "\n")
  cat("  Stop reason:", x$stop_reason, "\n")
  cat("  Elapsed seconds:", format(x$elapsed, digits = 5), "\n")
  invisible(x)
}

#' @export
summary.cureboost_joint <- function(object, ...) {
  selected <- object$history[object$history$selected, , drop = FALSE]
  out <- list(
    call = object$call,
    n = object$n,
    p = object$p,
    em_iterations = object$em_iterations,
    incidence_trees = length(object$incidence_trees),
    latency_trees = length(object$latency_trees),
    latency_baseline = object$latency_baseline,
    mean_cure = mean(object$cure_probability),
    mean_posterior_uncured = mean(object$posterior_uncured),
    final_training_loss = utils::tail(selected$training_loss, 1L),
    best_validation_loss = utils::tail(selected$validation_loss, 1L),
    stop_reason = object$stop_reason,
    importance = cureboost_joint_importance(object),
    elapsed = object$elapsed
  )
  class(out) <- "summary.cureboost_joint"
  out
}

#' @export
print.summary.cureboost_joint <- function(x, ...) {
  cat("Joint incidence-latency cure boosting summary\n")
  cat("  n =", x$n, "; p =", x$p, "; EM iterations =",
      x$em_iterations, "\n")
  cat("  incidence / latency trees =", x$incidence_trees, "/",
      x$latency_trees, "\n")
  cat("  susceptible latency baseline =", x$latency_baseline, "\n")
  cat("  mean cure =", format(x$mean_cure, digits = 5),
      "; posterior uncured =",
      format(x$mean_posterior_uncured, digits = 5), "\n")
  cat("  training / validation loss =",
      format(x$final_training_loss, digits = 5), "/",
      format(x$best_validation_loss, digits = 5), "\n")
  cat("  stop reason:", x$stop_reason, "\n")
  invisible(x)
}
