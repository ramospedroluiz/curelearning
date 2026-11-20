# curelearning

Statistical learning for cure models and long-term survival with right-censored data. Fit nearest neighbors (`cureknn()`), trees (`curetree()`), forests (`cureforest()`) and boosting models (`cureboost()` and `cureboost_joint()`) through a formula interface. The nearest-neighbor method estimates survival at specified times, not an individual cure probability.

## Installation from GitHub

```r
# install.packages("remotes")  # if needed
remotes::install_github("ramospedroluiz/curelearning", upgrade = "never")
library(curelearning)
```

To reproduce a particular analysis, install the commit recorded with that
analysis using the `ref` argument to `install_github()`.

## Installation from a source archive

Install the required dependencies (`survival`, `rpart` and `FNN`) and run:

```r
install.packages("curelearning_0.1.0.tar.gz", repos = NULL, type = "source")
library(curelearning)
```

For new analyses use `library(curelearning)` or a qualified call such as `curelearning::cureforest()`. The original `cureforest` package is a separate historical distribution. The model functions and their prediction interfaces keep their names.

## Nearest neighbors for late survival

A runnable example with the included synthetic cohort:

```r
data(cure_latency_demo, package = "curelearning")
train <- subset(cure_latency_demo, cohort == "train")
patients <- head(subset(cure_latency_demo, cohort == "test"), 3)
fit <- cureknn(survival::Surv(time, event) ~ age + stage + cure_marker,
               data = train, times = c(7, 8, 9), seed = 1)
predict(fit, patients)
head(predict(fit, patients, type = "neighbors"))
```

The synthetic cohort has an administrative cutoff at year 10, so the example predicts at years 7--9. Extending the window to the cutoff can withhold a prediction. Unsupported estimates are missing, not zero. See `?cureknn` for the follow-up check, censoring model and fixed sample splits for temporal designs.

## Two-minute comparison with randomForestSRC

The package includes a prespecified train/test data set in which one marker
changes cure probability and another changes failure timing among susceptible
subjects. The two forests are fitted explicitly, just as in a standard
random-forest example:

```r
install.packages("randomForestSRC")
library(curelearning)
library(survival)

data(cure_latency_demo)
train <- subset(cure_latency_demo, cohort == "train")
test <- subset(cure_latency_demo, cohort == "test")
form <- Surv(time, event) ~
  cure_marker + latency_marker + age + stage + biomarker + noise

fit.cure <- cureforest(
  form, train, ntree = 200, mtry = 6, nodesize = 80, maxdepth = 1,
  tail.start = 7, tail.end = 9.5, n.cores = 6, seed = 20260725
)

fit.rsf <- randomForestSRC::rfsrc(
  form, train, ntree = 200, mtry = 6, nodesize = 80, nodedepth = 1,
  splitrule = "logrank", nthread = 6, seed = 20260725
)

pred.cure <- predict(fit.cure, test, type = "cure")
head(pred.cure)

comparison <- compare_forest_fits(fit.cure, fit.rsf, test, n.cores = 6)
comparison$performance
```

The returned table reports the held-out C-index for the proposed cure forest
and for the ordinary log-rank random survival forest, together with Brier
scores and IBS. The complete script is installed with the package:

```r
system.file("examples", "cureforest-vs-rsf.R", package = "curelearning")
```

The example is a controlled illustration of the difference between cure and
log-rank splitting, not an empirical superiority claim.

## Cure-fraction survival tree

```r
library(survival)
library(curelearning)

fit_tree <- curetree(
  Surv(time, status) ~ age + stage + biomarker,
  data = cancer_data,
  tail.start = 8, tail.end = 10,
  nodesize = 50, maxdepth = 3, nsplit = 10
)

print(fit_tree)
summary(fit_tree)
plot(fit_tree)
predict(fit_tree, newdata = cancer_data, type = "cure")
```

## Cure-fraction random survival forest

```r
fit_forest <- cureforest(
  Surv(time, status) ~ age + stage + biomarker,
  data = cancer_data,
  ntree = 500,
  mtry = 2,
  nodesize = 50,
  maxdepth = 4,
  nsplit = 10,
  sampsize = 0.8,
  honesty = TRUE,
  n.cores = 4,
  seed = 2026
)

fit_forest$predicted                  # out-of-bag cure predictions
cure_importance(fit_forest)
predict(fit_forest, cancer_data, type = "cure")
predict(fit_forest, cancer_data, type = "survival", times = c(5, 8, 10))
predict(fit_forest, cancer_data, type = "susceptible", times = c(5, 8, 10))
predict(fit_forest, cancer_data, type = "dynamic", times = c(5, 8, 10))
```

## Tail-targeted cure boosting

`cureboost()` adds shallow trees sequentially. Each tree corrects finite-tail
structure not explained by the preceding trees, and every candidate daughter
must satisfy the same follow-up maturity restrictions used by the cure forest.

```r
fit_boost <- cureboost(
  form, train,
  nrounds = 100,
  learning.rate = 0.05,
  maxdepth = 1,
  nodesize = 80,
  tail.start = 7,
  tail.end = 9.5,
  ipcw.model = "cox",
  validation.fraction = 0.2,
  early.stopping.rounds = 15,
  seed = 2026
)

pred.boost <- predict(fit_boost, test, type = "cure")
cureboost_importance(fit_boost)
fit_boost$history
```

The complete console example is installed at:

```r
system.file("examples", "cureboost-example.R", package = "curelearning")
```

The current implementation estimates the observable finite-tail target using
cross-fitted IPCW pseudo-outcomes. It does not yet fit a separate susceptible-
latency component or return a complete survival curve.

## Covariate-dependent censoring

For censoring that depends on observed covariates, a cross-fitted Cox model can
estimate the censoring distribution:

```r
fit_ipcw <- randomforestcure(
  Surv(time, status) ~ age + stage + biomarker,
  data = cancer_data,
  estimator = "ipcw",
  ipcw.model = "cox",
  ipcw.folds = 5,
  ntree = 500
)
```

Alternatively, `ipcw.model = "user"` accepts a function `g.fun(t, newdata)`
that returns `P(C > t | X)` for every row of `newdata`.

## Tree and forest predictions

For fitted trees and forests, `predict()` supports:

- `"cure"`: long-term survival or cure probability;
- `"survival"`: population survival;
- `"susceptible"`: survival among susceptible individuals;
- `"dynamic"`: probability of cure conditional on survival to time `t`;
- `"terminal"`: terminal-node identifiers;
- `"all"`: all available functionals.

Within every tree, population survival is constrained to be no smaller than
that tree's estimated tail probability; the coherent tree curves are then
averaged. This order matters in finite forests. The returned population and
susceptible curves satisfy the mixture identity numerically.

## Practical guidance

- Choose the tail window from scientific considerations and follow-up support,
  not by maximizing test-set performance.
- Inspect `summary()` for terminal-node fallback and maturity diagnostics.
- Use `splitrule = "logrank"` with the same controls for a matched comparison.
- Treat IJ standard errors as experimental diagnostics. Requesting them
  requires a fit made with `inference = TRUE`, which enforces one-shot
  fixed-size subsampling without replacement, a sublinear integer subsample
  schedule, prespecified tail times, positive score stabilization, and
  external user-supplied IPCW weights when IPCW is used. Calibration around
  the forest expectation does not remove localization or residual-tail bias,
  so target-centered intervals require additional bias rates;
  `se.fit = FALSE` remains the default.
- When IJ output is requested, inspect `ij.unstable` and
  `finite.B.correction.fraction`; a warning means that the finite-tree
  subtraction exhausted the raw IJ variance and more trees are needed before
  using the corrected standard error.
- `nsplit` values between 8 and 20 are often a useful starting range.
- `n.cores > 1` uses PSOCK workers and is supported on Windows, macOS, and
  Unix-like systems.

The original functions `survcure()` and `randomforestcure()` remain available
for backward compatibility.

See `vignette("curelearning-introduction", package = "curelearning")` for a complete
reproducible example.
