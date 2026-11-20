.cb_node_maturity <- function(index, context) {
  time <- context$time[index]
  status <- context$status[index]
  n_tail <- sum(time >= context$tail_start)
  n_tail_end <- sum(time >= context$tail_end)
  n_events <- sum(status == 1L & time <= context$tail_start)
  list(
    valid = length(index) >= context$nodesize &&
      n_tail >= context$min_tail_at_risk &&
      n_tail_end >= context$min_tail_end_at_risk &&
      n_events >= context$min_events_before_tail,
    n_tail = n_tail,
    n_tail_end = n_tail_end,
    n_events = n_events
  )
}

.cb_residual_variance <- function(residual) {
  if (length(residual) <= 1L) return(0)
  value <- stats::var(residual) / length(residual)
  if (is.finite(value)) value else 0
}

.cb_find_best_split <- function(index, residual, context) {
  p <- ncol(context$x)
  variables <- if (context$mtry >= p) {
    seq_len(p)
  } else {
    sample.int(p, context$mtry)
  }
  x_node <- context$x[index, , drop = FALSE]
  best <- NULL
  best_score <- -Inf

  for (variable in variables) {
    cutpoints <- .cfr_cutpoints(x_node[, variable], context$nsplit)
    if (!length(cutpoints)) next
    for (cutpoint in cutpoints) {
      left <- x_node[, variable] <= cutpoint
      if (anyNA(left)) next
      n_left <- sum(left)
      n_right <- length(left) - n_left
      if (n_left < context$nodesize || n_right < context$nodesize) next

      left_index <- index[left]
      right_index <- index[!left]
      left_maturity <- .cb_node_maturity(left_index, context)
      right_maturity <- .cb_node_maturity(right_index, context)
      if (!left_maturity$valid || !right_maturity$valid) next

      left_residual <- residual[left]
      right_residual <- residual[!left]
      difference <- mean(left_residual) - mean(right_residual)
      balance <- (n_left * n_right) / length(index)^2
      denominator <- .cb_residual_variance(left_residual) +
        .cb_residual_variance(right_residual) + context$lambda
      score <- balance * difference^2 / denominator

      if (is.finite(score) && score > best_score) {
        best_score <- score
        best <- list(
          variable = variable,
          cutpoint = cutpoint,
          score = score,
          left = left,
          default_left = n_left >= n_right
        )
      }
    }
  }
  best
}

.cb_grow_tree <- function(index, residual, context) {
  max_nodes <- 2^(context$maxdepth + 1L) - 1L
  node_store <- vector("list", max_nodes)
  node_count <- 0L

  grow <- function(node_index, node_residual, depth) {
    node_count <<- node_count + 1L
    node_id <- node_count
    maturity <- .cb_node_maturity(node_index, context)
    node_store[[node_id]] <<- data.frame(
      node = node_id,
      depth = depth,
      split_variable = NA_integer_,
      split_value = NA_real_,
      split_score = NA_real_,
      left = NA_integer_,
      right = NA_integer_,
      default_left = TRUE,
      value = mean(node_residual),
      n = length(node_index),
      n_tail = maturity$n_tail,
      n_tail_end = maturity$n_tail_end,
      n_events = maturity$n_events,
      terminal = TRUE
    )

    can_split <- depth < context$maxdepth &&
      length(node_index) >= 2L * context$nodesize
    split <- if (can_split) {
      .cb_find_best_split(node_index, node_residual, context)
    } else NULL
    if (is.null(split) || split$score <= context$min_gain) return(node_id)

    left_id <- grow(
      node_index[split$left], node_residual[split$left], depth + 1L
    )
    right_id <- grow(
      node_index[!split$left], node_residual[!split$left], depth + 1L
    )
    node_store[[node_id]]$split_variable <<- split$variable
    node_store[[node_id]]$split_value <<- split$cutpoint
    node_store[[node_id]]$split_score <<- split$score
    node_store[[node_id]]$left <<- left_id
    node_store[[node_id]]$right <<- right_id
    node_store[[node_id]]$default_left <<- split$default_left
    node_store[[node_id]]$value <<- NA_real_
    node_store[[node_id]]$terminal <<- FALSE
    node_id
  }

  grow(index, residual, 0L)
  nodes <- do.call(rbind, node_store[seq_len(node_count)])
  rownames(nodes) <- NULL
  list(nodes = nodes, variable_names = context$variable_names)
}

.cb_route_tree <- function(tree, x) {
  .cfr_route_tree(tree, x)
}

.cb_predict_tree <- function(tree, x, terminal = FALSE) {
  node <- .cb_route_tree(tree, x)
  if (terminal) return(node)
  tree$nodes$value[node]
}

.cb_validation_split <- function(n, fraction, seed, context) {
  if (fraction <= 0) {
    return(list(train = seq_len(n), validation = integer(0)))
  }
  n_validation <- max(1L, floor(n * fraction))
  if (n - n_validation < 2L * context$nodesize) {
    stop(
      "The requested validation fraction leaves too few observations for boosting.",
      call. = FALSE
    )
  }
  for (attempt in 0:49) {
    split_seed <- as.integer(
      (as.double(seed) + 104729 * attempt) %% .Machine$integer.max
    )
    set.seed(split_seed)
    validation <- sample.int(n, n_validation, replace = FALSE)
    train <- setdiff(seq_len(n), validation)
    if (.cb_node_maturity(train, context)$valid) {
      return(list(train = train, validation = validation))
    }
  }
  stop(
    "Could not construct a mature boosting sample after 50 validation splits.",
    call. = FALSE
  )
}

.cb_build_importance <- function(trees, variable_names) {
  count <- numeric(length(variable_names))
  gain <- numeric(length(variable_names))
  root <- numeric(length(variable_names))
  if (!length(trees)) {
    return(data.frame(
      variable = variable_names,
      split_count = count,
      total_gain = gain,
      mean_gain = gain,
      root_frequency = root,
      row.names = NULL
    ))
  }
  for (tree in trees) {
    internal <- which(!tree$nodes$terminal)
    if (!length(internal)) next
    variable <- tree$nodes$split_variable[internal]
    weighted_gain <- tree$nodes$split_score[internal] * tree$coefficient^2
    count <- count + tabulate(variable, nbins = length(variable_names))
    for (j in seq_along(variable)) {
      gain[variable[j]] <- gain[variable[j]] + weighted_gain[j]
    }
    root_variable <- tree$nodes$split_variable[1L]
    if (is.finite(root_variable)) root[root_variable] <- root[root_variable] + 1
  }
  data.frame(
    variable = variable_names,
    split_count = count,
    total_gain = gain,
    mean_gain = ifelse(count > 0, gain / count, 0),
    root_frequency = root / length(trees),
    row.names = NULL
  )
}

.cb_fitted_raw <- function(initial, trees, x) {
  value <- rep(initial, nrow(x))
  if (!length(trees)) return(value)
  for (tree in trees) {
    value <- value + tree$coefficient * .cb_predict_tree(tree, x)
  }
  value
}

#' Fit tail-targeted boosting for cure probability
#'
#' Fits an experimental gradient-boosting ensemble for the observable
#' finite-tail survival target. Each shallow tree corrects residual
#' tail-survival structure using the cure-forest maturity restrictions.
#'
#' @param formula A formula with a `survival::Surv(time, status)` response.
#' @param data A data frame.
#' @param nrounds Maximum number of boosting trees.
#' @param learning.rate Shrinkage applied to each correction.
#' @param validation.fraction Fraction held out for early stopping; use zero to
#'   fit all observations.
#' @param early.stopping.rounds Number of non-improving validation iterations
#'   before stopping.
#' @return An object of class `cureboost`.
#' @export
cureboost <- function(formula,
                      data,
                      nrounds = 100L,
                      learning.rate = 0.05,
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
                      early.stopping.rounds = 15L,
                      stopping.tolerance = 1e-8,
                      min.correlation = 1e-8,
                      seed = 1L,
                      na.action = stats::na.omit,
                      keep.data = FALSE,
                      verbose = FALSE) {
  call <- match.call()
  .cfr_assert_scalar(nrounds, "nrounds", 1, 100000, integer = TRUE)
  .cfr_assert_scalar(learning.rate, "learning.rate", 0, 1)
  if (learning.rate <= 0) {
    stop("learning.rate must be strictly positive.", call. = FALSE)
  }
  .cfr_assert_scalar(lambda, "lambda", 0, Inf)
  if (lambda <= 0) stop("lambda must be strictly positive.", call. = FALSE)
  .cfr_assert_scalar(validation.fraction, "validation.fraction", 0, 0.5)
  .cfr_assert_scalar(early.stopping.rounds, "early.stopping.rounds", 1,
                     100000, integer = TRUE)
  .cfr_assert_scalar(stopping.tolerance, "stopping.tolerance", 0, Inf)
  .cfr_assert_scalar(min.correlation, "min.correlation", 0, 1)
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
    curve.grid.size = max(2L, as.integer(tail.grid.size)),
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
  pseudo <- context$ipcw_tail_mean
  if (any(!is.finite(pseudo))) {
    stop("The IPCW tail pseudo-outcome contains non-finite values.",
         call. = FALSE)
  }
  split <- .cb_validation_split(
    fit$prepared$n, validation.fraction, seed, context
  )
  train <- split$train
  validation <- split$validation
  initial <- mean(pseudo[train])
  fitted_train <- rep(initial, length(train))
  fitted_validation <- rep(initial, length(validation))
  best_validation <- if (length(validation)) {
    mean((pseudo[validation] - fitted_validation)^2) / 2
  } else Inf
  best_iteration <- 0L
  patience <- 0L
  trees <- vector("list", nrounds)
  history <- vector("list", nrounds)
  accepted <- 0L
  stop_reason <- "maximum rounds reached"
  progress_points <- unique(pmax(1L, round(nrounds * c(0.25, 0.50, 0.75, 1))))

  start <- proc.time()[[3L]]
  set.seed(seed)
  for (m in seq_len(nrounds)) {
    residual <- pseudo[train] - fitted_train
    risk_before <- mean(residual^2) / 2
    tree <- .cb_grow_tree(train, residual, context)
    raw_direction <- .cb_predict_tree(
      tree, context$x[train, , drop = FALSE]
    )
    direction_norm <- sqrt(mean(raw_direction^2))
    if (!is.finite(direction_norm) || direction_norm <= .Machine$double.eps) {
      stop_reason <- "no nonzero mature-tail correction remained"
      break
    }

    tree$nodes$value <- tree$nodes$value / direction_norm
    direction_train <- raw_direction / direction_norm
    correlation <- mean(residual * direction_train)
    if (correlation < 0) {
      tree$nodes$value <- -tree$nodes$value
      direction_train <- -direction_train
      correlation <- -correlation
    }
    residual_norm <- sqrt(mean(residual^2))
    correlation_ratio <- correlation / max(
      residual_norm, .Machine$double.eps
    )
    if (!is.finite(correlation_ratio) ||
        correlation_ratio < min.correlation) {
      stop_reason <- "the accepted tree failed the numerical weak-learning check"
      break
    }

    coefficient <- learning.rate * correlation
    proposed_train <- fitted_train + coefficient * direction_train
    risk_after <- mean((pseudo[train] - proposed_train)^2) / 2
    reduction <- risk_before - risk_after
    if (!is.finite(reduction) || reduction <= stopping.tolerance) {
      stop_reason <- "empirical loss reduction was below stopping.tolerance"
      break
    }

    tree$coefficient <- coefficient
    tree$gradient_correlation <- correlation
    tree$gradient_correlation_ratio <- correlation_ratio
    tree$training_reduction <- reduction
    accepted <- accepted + 1L
    trees[[accepted]] <- tree
    fitted_train <- proposed_train

    validation_risk <- NA_real_
    if (length(validation)) {
      direction_validation <- .cb_predict_tree(
        tree, context$x[validation, , drop = FALSE]
      )
      fitted_validation <- fitted_validation +
        coefficient * direction_validation
      validation_risk <- mean(
        (pseudo[validation] - fitted_validation)^2
      ) / 2
      if (validation_risk < best_validation - stopping.tolerance) {
        best_validation <- validation_risk
        best_iteration <- accepted
        patience <- 0L
      } else {
        patience <- patience + 1L
      }
    } else {
      best_iteration <- accepted
    }

    internal <- which(!tree$nodes$terminal)
    history[[accepted]] <- data.frame(
      iteration = accepted,
      training_risk = risk_after,
      validation_risk = validation_risk,
      reduction = reduction,
      theoretical_reduction =
        (learning.rate - learning.rate^2 / 2) * correlation^2,
      correlation = correlation,
      correlation_ratio = correlation_ratio,
      terminal_nodes = sum(tree$nodes$terminal),
      root_variable = if (length(internal)) {
        context$variable_names[tree$nodes$split_variable[1L]]
      } else NA_character_,
      row.names = NULL
    )

    if (verbose && m %in% progress_points) {
      message(round(100 * m / nrounds), "% of requested boosting rounds completed")
    }
    if (length(validation) && patience >= early.stopping.rounds) {
      stop_reason <- "validation early stopping"
      break
    }
  }

  trees <- trees[seq_len(accepted)]
  history <- if (accepted) do.call(rbind, history[seq_len(accepted)]) else {
    data.frame(
      iteration = integer(0), training_risk = numeric(0),
      validation_risk = numeric(0), reduction = numeric(0),
      theoretical_reduction = numeric(0), correlation = numeric(0),
      correlation_ratio = numeric(0), terminal_nodes = integer(0),
      root_variable = character(0)
    )
  }
  if (!length(validation)) best_iteration <- accepted
  if (best_iteration < accepted) {
    trees <- if (best_iteration) trees[seq_len(best_iteration)] else list()
  }
  if (nrow(history)) history$selected <- history$iteration <= best_iteration

  fitted_raw <- .cb_fitted_raw(initial, trees, context$x)
  metadata <- .cfr_fit_metadata(
    fit$prepared, context, call, keep.data
  )
  importance <- .cb_build_importance(trees, fit$prepared$variable_names)
  object <- c(
    list(
      trees = trees,
      nrounds = length(trees),
      nrounds_requested = as.integer(nrounds),
      accepted_rounds = accepted,
      best_iteration = best_iteration,
      learning_rate = learning.rate,
      initial = initial,
      fitted = .cfr_clip(fitted_raw),
      fitted_raw = fitted_raw,
      pseudo_outcome = pseudo,
      train_index = train,
      validation_index = validation,
      validation_fraction = validation.fraction,
      history = history,
      stop_reason = stop_reason,
      importance = importance,
      target = "finite_tail",
      elapsed = proc.time()[[3L]] - start
    ),
    metadata
  )
  class(object) <- "cureboost"
  object
}

#' Variable importance for tail-targeted cure boosting
#'
#' @param object A fitted `cureboost` object.
#' @param aggregate If `TRUE`, aggregate encoded columns to formula terms.
#' @return A data frame ordered by total weighted split gain.
#' @export
cureboost_importance <- function(object, aggregate = TRUE) {
  if (!inherits(object, "cureboost")) {
    stop("object must inherit from cureboost.", call. = FALSE)
  }
  importance <- object$importance
  if (isTRUE(aggregate)) {
    importance <- .cfr_aggregate_importance(
      importance, object$variable_map
    )
  }
  importance[order(importance$total_gain, decreasing = TRUE), , drop = FALSE]
}

#' Predict from tail-targeted cure boosting
#'
#' @param object A fitted `cureboost` object.
#' @param newdata New predictor data. If omitted, stored fitted values are used
#'   for ordinary response prediction.
#' @param type Prediction type: clipped finite-tail response, raw additive
#'   predictor, terminal nodes, per-tree contributions, or all outputs.
#' @param nrounds Number of fitted boosting trees to use.
#' @export
predict.cureboost <- function(object,
                              newdata = NULL,
                              type = c("cure", "response", "raw", "terminal",
                                       "contributions", "all"),
                              nrounds = NULL,
                              ...) {
  type <- match.arg(type)
  if (type == "response") type <- "cure"
  if (is.null(nrounds)) nrounds <- object$nrounds
  .cfr_assert_scalar(nrounds, "nrounds", 0, object$nrounds, integer = TRUE)

  if (is.null(newdata) && nrounds == object$nrounds &&
      type %in% c("cure", "raw")) {
    value <- if (type == "cure") object$fitted else object$fitted_raw
    names(value) <- rownames(object$training_data)
    return(value)
  }
  if (is.null(newdata)) {
    if (is.null(object$training_data)) {
      stop(
        "newdata is required for this prediction because the training data were not retained.",
        call. = FALSE
      )
    }
    newdata <- object$training_data
  }
  x <- .cfr_new_matrix(object, newdata)
  use_trees <- if (nrounds) object$trees[seq_len(nrounds)] else list()
  contributions <- matrix(
    object$initial, nrow = nrow(x), ncol = nrounds + 1L
  )
  colnames(contributions) <- c(
    "initial", if (nrounds) paste0("tree", seq_len(nrounds)) else character(0)
  )
  terminal <- matrix(NA_integer_, nrow(x), nrounds)
  if (nrounds) {
    for (m in seq_len(nrounds)) {
      tree <- use_trees[[m]]
      contributions[, m + 1L] <- tree$coefficient *
        .cb_predict_tree(tree, x)
      terminal[, m] <- .cb_predict_tree(tree, x, terminal = TRUE)
    }
    colnames(terminal) <- paste0("tree", seq_len(nrounds))
  }
  raw <- rowSums(contributions)
  cure <- .cfr_clip(raw)
  names(raw) <- names(cure) <- rownames(newdata)
  rownames(contributions) <- rownames(terminal) <- rownames(newdata)

  if (type == "cure") return(cure)
  if (type == "raw") return(raw)
  if (type == "terminal") return(terminal)
  if (type == "contributions") return(contributions)
  list(
    cure = cure,
    raw = raw,
    contributions = contributions,
    terminal = terminal,
    nrounds = nrounds
  )
}

#' @export
print.cureboost <- function(x, ...) {
  cat("Tail-targeted cure boosting\n")
  cat("  Observations / encoded predictors:", x$n, "/", x$p, "\n")
  cat("  Selected / accepted / requested rounds:", x$nrounds, "/",
      x$accepted_rounds, "/", x$nrounds_requested, "\n")
  cat("  Learning rate / maximum depth:", x$learning_rate, "/",
      x$control$maxdepth, "\n")
  cat("  Tail window:", format(x$tail_start, digits = 5), "to",
      format(x$tail_end, digits = 5), "\n")
  cat("  Target / censoring model:", x$target, "/", x$ipcw_model, "\n")
  cat("  Stop reason:", x$stop_reason, "\n")
  cat("  Elapsed seconds:", format(x$elapsed, digits = 5), "\n")
  invisible(x)
}

#' @export
summary.cureboost <- function(object, ...) {
  selected_history <- if (nrow(object$history)) {
    object$history[object$history$selected, , drop = FALSE]
  } else object$history
  out <- list(
    call = object$call,
    n = object$n,
    p = object$p,
    nrounds = object$nrounds,
    accepted_rounds = object$accepted_rounds,
    learning_rate = object$learning_rate,
    tail_window = c(object$tail_start, object$tail_end),
    target = object$target,
    ipcw_model = object$ipcw_model,
    initial = object$initial,
    final_training_risk = if (nrow(selected_history)) {
      utils::tail(selected_history$training_risk, 1L)
    } else mean((object$pseudo_outcome[object$train_index] - object$initial)^2) / 2,
    best_validation_risk = if (length(object$validation_index)) {
      if (nrow(selected_history)) {
        utils::tail(selected_history$validation_risk, 1L)
      } else mean((object$pseudo_outcome[object$validation_index] -
                   object$initial)^2) / 2
    } else NA_real_,
    stop_reason = object$stop_reason,
    importance = cureboost_importance(object),
    elapsed = object$elapsed
  )
  class(out) <- "summary.cureboost"
  out
}

#' @export
print.summary.cureboost <- function(x, ...) {
  cat("Tail-targeted cure boosting summary\n")
  cat("  n =", x$n, "; p =", x$p, "; selected rounds =", x$nrounds, "\n")
  cat("  finite-tail window =", format(x$tail_window[1L], digits = 5),
      "to", format(x$tail_window[2L], digits = 5), "\n")
  cat("  final training risk =", format(x$final_training_risk, digits = 5),
      "; validation risk =", format(x$best_validation_risk, digits = 5), "\n")
  cat("  stop reason:", x$stop_reason, "\n")
  cat("\nLeading variables by weighted gain:\n")
  print(utils::head(x$importance, 10L), row.names = FALSE)
  invisible(x)
}
