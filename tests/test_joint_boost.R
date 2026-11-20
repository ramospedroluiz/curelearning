if (file.exists(file.path("R", "internals.R"))) {
  source(file.path("R", "internals.R"))
  source(file.path("R", "fit.R"))
  source(file.path("R", "predict.R"))
  source(file.path("R", "boost.R"))
  source(file.path("R", "joint_boost.R"))
  source(file.path("R", "cure_cox_boost.R"))
} else {
  library(curelearning)
}

simulate_joint_data <- function(n = 700L, seed = 20260816L) {
  set.seed(seed)
  incidence <- runif(n)
  latency <- rnorm(n)
  noise <- rnorm(n)
  cure_probability <- plogis(-0.5 + 3 * (incidence - 0.5))
  cured <- rbinom(n, 1L, cure_probability)
  susceptible_time <- rexp(n, rate = 0.20 * exp(1.2 * latency))
  event_time <- ifelse(cured == 1L, Inf, susceptible_time)
  censor_time <- runif(n, 6, 13)
  data.frame(
    time = pmin(event_time, censor_time),
    status = as.integer(event_time <= censor_time),
    incidence = incidence,
    latency = latency,
    noise = noise,
    cure_probability = cure_probability,
    cured = cured
  )
}

data <- simulate_joint_data()
fit <- cureboost_joint(
  survival::Surv(time, status) ~ incidence + latency + noise,
  data = data,
  em.iterations = 15L,
  incidence.steps = 2L,
  latency.steps = 2L,
  initialization.rounds = 8L,
  learning.rate.incidence = 0.10,
  learning.rate.latency = 0.10,
  maxdepth = 1L,
  nodesize = 30L,
  mtry = 3L,
  nsplit = 8L,
  tail.start = 7,
  tail.end = 10,
  tail.times = seq(7, 10, length.out = 5L),
  min.tail.at.risk = 10L,
  min.tail.end.at.risk = 3L,
  min.events.before.tail = 4L,
  ipcw.model = "marginal",
  validation.fraction = 0.20,
  early.stopping.rounds = 5L,
  seed = 20260816L,
  keep.data = TRUE
)

prediction <- predict(
  fit, data[1:50, ], type = "all", times = c(2, 5, 8, 10)
)
stopifnot(
  inherits(fit, "cureboost_joint"),
  identical(fit$latency_baseline, "weibull"),
  identical(fit$baseline$type, "weibull"),
  is.finite(fit$baseline$shape),
  fit$baseline$shape > 0,
  fit$em_iterations >= 1L,
  length(fit$incidence_trees) >= fit$initialization_rounds,
  length(fit$latency_trees) >= 1L,
  all(is.finite(prediction$cure)),
  all(prediction$cure > 0 & prediction$cure < 1),
  all(prediction$susceptible >= 0 & prediction$susceptible <= 1),
  all(prediction$survival >= 0 & prediction$survival <= 1),
  all(apply(prediction$susceptible, 1L, function(x) all(diff(x) <= 1e-12))),
  all(apply(prediction$survival, 1L, function(x) all(diff(x) <= 1e-12))),
  max(abs(
    prediction$survival -
      (prediction$cure + (1 - prediction$cure) * prediction$susceptible)
  )) < 1e-12,
  all(fit$posterior_uncured[data$status == 1L] > 0.99),
  all(is.finite(fit$history$training_loss)),
  all(is.finite(fit$history$validation_loss))
)

importance <- cureboost_joint_importance(fit)
stopifnot(
  all(c("incidence", "latency") %in% importance$component),
  "incidence" %in% importance$variable,
  "latency" %in% importance$variable
)

summary_fit <- summary(fit)
stopifnot(
  inherits(summary_fit, "summary.cureboost_joint"),
  is.finite(summary_fit$final_training_loss),
  is.finite(summary_fit$best_validation_loss)
)

# The earlier semiparametric baseline remains available as a sensitivity
# analysis and old saved objects remain predictable after the baseline wrapper
# introduced for the Weibull implementation.
fit_breslow <- cureboost_joint(
  survival::Surv(time, status) ~ incidence + latency + noise,
  data = data,
  em.iterations = 3L,
  incidence.steps = 1L,
  latency.steps = 1L,
  initialization.rounds = 0L,
  latency.baseline = "breslow",
  maxdepth = 1L,
  nodesize = 30L,
  mtry = 3L,
  nsplit = 5L,
  tail.start = 7,
  tail.end = 10,
  tail.times = seq(7, 10, length.out = 5L),
  min.tail.at.risk = 10L,
  min.tail.end.at.risk = 3L,
  min.events.before.tail = 4L,
  ipcw.model = "marginal",
  validation.fraction = 0,
  seed = 20260817L
)
breslow_prediction <- predict(
  fit_breslow, data[1:10, ], type = "survival", times = c(2, 5, 8)
)
stopifnot(
  identical(fit_breslow$latency_baseline, "breslow"),
  identical(fit_breslow$baseline$type, "breslow"),
  all(is.finite(breslow_prediction)),
  all(apply(breslow_prediction, 1L, function(x) all(diff(x) <= 1e-12)))
)

anchor <- cureboost(
  survival::Surv(time, status) ~ incidence + latency + noise,
  data = data,
  nrounds = 12L,
  learning.rate = 0.10,
  maxdepth = 1L,
  nodesize = 30L,
  mtry = 3L,
  nsplit = 5L,
  tail.start = 7,
  tail.end = 10,
  tail.times = seq(7, 10, length.out = 5L),
  min.tail.at.risk = 10L,
  min.tail.end.at.risk = 3L,
  min.events.before.tail = 4L,
  ipcw.model = "marginal",
  validation.fraction = 0,
  seed = 20260818L
)
fit_anchored <- cureboost_joint(
  survival::Surv(time, status) ~ incidence + latency + noise,
  data = data,
  em.iterations = 5L,
  incidence.steps = 0L,
  latency.steps = 2L,
  initialization.rounds = 0L,
  latency.baseline = "breslow",
  incidence.mode = "anchored",
  incidence.fit = anchor,
  maxdepth = 1L,
  nodesize = 30L,
  mtry = 3L,
  nsplit = 5L,
  tail.start = 7,
  tail.end = 10,
  tail.times = seq(7, 10, length.out = 5L),
  min.tail.at.risk = 10L,
  min.tail.end.at.risk = 3L,
  min.events.before.tail = 4L,
  ipcw.model = "marginal",
  validation.fraction = 0,
  seed = 20260818L
)
anchored_prediction <- predict(
  fit_anchored, data[1:20, ], type = "all", times = c(2, 5, 8)
)
stopifnot(
  identical(fit_anchored$incidence_mode, "anchored"),
  identical(fit_anchored$latency_baseline, "breslow"),
  max(abs(
    anchored_prediction$cure -
      predict(anchor, data[1:20, ], type = "cure")
  )) < 1e-12,
  all(apply(anchored_prediction$survival, 1L,
            function(x) all(diff(x) <= 1e-12)))
)

fit_cox <- curecoxboost(
  survival::Surv(time, status) ~ incidence + latency + noise,
  data = data,
  incidence.fit = anchor,
  latency.rounds = 3L,
  latency.steps = 1L,
  learning.rate = 0.10,
  latency.maxdepth = 1L,
  nodesize = 30L,
  mtry = 3L,
  nsplit = 5L,
  tail.start = 7,
  tail.end = 10,
  tail.times = seq(7, 10, length.out = 5L),
  min.tail.at.risk = 10L,
  min.tail.end.at.risk = 3L,
  min.events.before.tail = 4L,
  validation.fraction = 0,
  seed = 20260818L
)
stopifnot(
  inherits(fit_cox, "curecoxboost"),
  inherits(fit_cox, "cureboost_joint"),
  identical(fit_cox$target, "cure_anchored_cox_latency"),
  identical(fit_cox$latency_baseline, "breslow"),
  identical(fit_cox$incidence_mode, "anchored")
)

cat("joint cureboost tests passed\n")
