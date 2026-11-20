# curelearning 0.1.0

- First release candidate under the new package name. The statistical methods,
  exported fitting functions, and prediction interfaces are unchanged from
  cureforest 0.2.0.9002. Use `library(curelearning)` in new work.
- The entries below record the history of the predecessor package.

# cureforest 0.2.0.9002 (development)

- Added `cureknn()` for late survival, with independent tuning, a local
  follow-up check and access to the neighbors of each component.

# cureforest 0.2.0.9001 (development)

- Added independent acceptance to `curecurveboost()`. The default `selection
  = "gated"` splits the validation pool into candidate selection and acceptance
  samples. A candidate is retained only when a one-sided paired lower bound
  for its mean Brier improvement is positive; otherwise round zero is returned.
- Kept `selection = "minimum"` for the earlier rule. The internal fitting
  engine retains its earlier default so existing benchmark scripts do not
  silently change their statistical procedure.
- Added gate diagnostics, disjoint split indices, a minimum-gain option, and
  tests that gate outcomes cannot affect training or candidate selection.
- The acceptance bound is approximate, not a new inferential guarantee for
  cure forests. The cure/tail anchor remains unchanged.

# cureforest 0.1.3

- Corrected forest survival aggregation: each tree curve is first constrained
  by that tree's tail-supported probability, and the coherent tree curves are
  then averaged.
- Restricted `inference = TRUE` to the kernel covered by the accompanying
  theory: one fixed-size subsample without replacement per tree, honest
  splitting, a sublinear integer subsample schedule, prespecified tail times,
  positive score stabilization, and external user-supplied IPCW weights when
  IPCW is used.
- Added a same-subsample root fallback for inference trees, eliminating
  reject-and-resample behavior while preserving a bounded prediction.
- IJ standard errors now require a fit created explicitly with
  `inference = TRUE`.
- Clarified throughout that a fixed tail window targets a finite-tail
  survival functional; interpreting it as exact cure additionally requires
  negligible residual susceptible survival.

# cureforest 0.1.2

- Added `compare_forest_fits()` so that a fitted cure forest and a fitted
  ordinary random survival forest can be compared directly on held-out data.
- Reworked the installed demonstration to show the two model fits and their
  predictions explicitly before calculating held-out C-index and Brier scores.

# cureforest 0.1.1

- Added the built-in `cure_latency_demo` data set for a self-contained,
  train/test illustration of cure-latency separation.
- Added `compare_cureforest()` to fit `cureforest` and `randomForestSRC` on the
  same split and report held-out Harrell C-index and marginal-IPCW Brier scores.
- Added `cureforest_demo()` as a one-command reproducible demonstration.

# cureforest 0.1.0

- Initial CRAN candidate release.
- Added `survcure()` for honest cure-fraction survival trees.
- Added `randomforestcure()` with OOB prediction and Windows parallelism.
- Added KM tail averages with linear-time Greenwood variance calculations.
- Added marginal, cross-fitted Cox, and user-supplied IPCW censoring weights.
- Added coherent predictions for population survival, susceptible survival,
  and dynamic cure probabilities.
- Added split-count, total-gain, and root-frequency variable importance.
- Added the fixed-size, without-replacement calibration factor to
  infinitesimal-jackknife variances and influence contributions. Simulations
  support variance calibration around the forest expectation; intervals for a
  cure target still require negligible forest and tail bias.
- Added a reproducible introductory vignette, package citation, CRAN submission
  notes, and portable API tests.
# cureforest 0.2.0

* Added the experimental `cureboost()` interface for stagewise estimation of
  the observable finite-tail target using mature-tail residual trees.
* Added validation-based early stopping, monotone training-loss diagnostics,
  `predict.cureboost()`, and `cureboost_importance()`.
* The first boosting implementation deliberately targets finite-tail survival
  through cross-fitted IPCW pseudo-outcomes; it does not yet fit a separate
  susceptible-latency component.
