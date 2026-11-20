if (file.exists(file.path("R", "internals.R"))) {
  source(file.path("R", "internals.R"))
  source(file.path("R", "fit.R"))
  source(file.path("R", "predict.R"))
  source(file.path("R", "boost.R"))
} else {
  library(curelearning)
}

simulate_boost_data <- function(n = 800, seed = 1) {
  set.seed(seed)
  x1 <- runif(n)
  x2 <- factor(rbinom(n, 1, 0.5), labels = c("low", "high"))
  x3 <- rnorm(n)
  pi0 <- plogis(-1 + 4 * (x1 - 0.5))
  cured <- rbinom(n, 1, pi0)
  susceptible_time <- pmin(rexp(n, ifelse(x2 == "high", 0.35, 0.14)), 9)
  event_time <- ifelse(cured == 1, Inf, susceptible_time)
  censor_time <- runif(n, 6, 12)
  data.frame(
    time = pmin(event_time, censor_time),
    status = as.integer(event_time <= censor_time),
    x1 = x1,
    x2 = x2,
    x3 = x3,
    pi0 = pi0
  )
}

data <- simulate_boost_data()
fit <- cureboost(
  survival::Surv(time, status) ~ x1 + x2 + x3,
  data = data,
  nrounds = 40,
  learning.rate = 0.10,
  maxdepth = 1,
  nodesize = 35,
  mtry = 3,
  nsplit = 10,
  tail.start = 7,
  tail.end = 9,
  tail.times = seq(7, 9, length.out = 5),
  min.tail.at.risk = 10,
  min.tail.end.at.risk = 3,
  min.events.before.tail = 4,
  ipcw.model = "marginal",
  validation.fraction = 0.20,
  early.stopping.rounds = 8,
  stopping.tolerance = 1e-12,
  seed = 2026,
  keep.data = TRUE
)

selected <- fit$history[fit$history$selected, , drop = FALSE]
stopifnot(
  inherits(fit, "cureboost"),
  fit$nrounds >= 1L,
  fit$nrounds <= fit$accepted_rounds,
  nrow(selected) == fit$nrounds,
  all(selected$reduction > 0),
  all(diff(selected$training_risk) <= 1e-12),
  max(abs(selected$reduction - selected$theoretical_reduction)) < 1e-10
)

prediction <- predict(fit, data[1:30, ], type = "cure")
raw <- predict(fit, data[1:30, ], type = "raw")
contributions <- predict(fit, data[1:30, ], type = "contributions")
all_prediction <- predict(fit, data[1:30, ], type = "all")
baseline <- predict(fit, data[1:30, ], type = "raw", nrounds = 0)
stopifnot(
  length(prediction) == 30L,
  all(is.finite(prediction)),
  all(prediction >= 0 & prediction <= 1),
  isTRUE(all.equal(raw, rowSums(contributions), tolerance = 1e-12)),
  isTRUE(all.equal(all_prediction$cure, prediction, tolerance = 1e-12)),
  all(baseline == fit$initial),
  identical(dim(all_prediction$terminal), c(30L, fit$nrounds)),
  isTRUE(all.equal(unname(predict(fit)), unname(fit$fitted), tolerance = 0))
)

importance <- cureboost_importance(fit)
stopifnot(
  nrow(importance) == 3L,
  sum(importance$split_count) > 0,
  all(importance$total_gain >= 0),
  "x1" %in% importance$variable
)

summary_fit <- summary(fit)
stopifnot(
  inherits(summary_fit, "summary.cureboost"),
  is.finite(summary_fit$final_training_risk),
  summary_fit$final_training_risk <
    mean((fit$pseudo_outcome[fit$train_index] - fit$initial)^2) / 2
)

bad_lambda <- try(cureboost(
  survival::Surv(time, status) ~ x1 + x2,
  data = data,
  nrounds = 2,
  lambda = 0,
  tail.start = 7,
  tail.end = 9
), silent = TRUE)
stopifnot(inherits(bad_lambda, "try-error"))

cat("cureboost tests passed\n")
