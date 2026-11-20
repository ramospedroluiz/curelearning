\name{cureknn}
\alias{cureknn}
\alias{predict.cureknn}
\alias{print.cureknn}
\title{Nearest neighbors for survival at prespecified late times}
\usage{
cureknn(formula, data, times, k = NULL, n_splits = 5L, seed = 1L,
        censoring = c("cox", "marginal"), censor_vars = NULL,
        active = 2L, rho = 0.20, delta = 0.05, g_floor = 0.02,
        omega_max = 1/g_floor, min_coverage = 0.80,
        splits = NULL, id = NULL, distance_weights = NULL)
\method{predict}{cureknn}(object, newdata,
        type = c("survival", "scalar", "neighbors"), ...)
\method{print}{cureknn}(x, ...)
}
\arguments{
\item{formula}{Right-censored \code{survival::Surv(time, event) ~ x1 + x2}
with untransformed baseline columns. Use 1 for the event, 0 for censoring.}
\item{data,newdata}{Training data or new patients. New patients require only
the distance predictors. Use factors for categorical numeric codes. Numeric
missing values are imputed from construction medians. Rare and new categories
are pooled. Do not supply identifiers or post-baseline variables as predictors.}
\item{times}{Distinct increasing positive horizons in the units of the response.}
\item{k}{Candidate neighborhood sizes. Default: 25, 50, 100, 200, 400,
capped at the smallest estimation sample and deduplicated.}
\item{n_splits,seed}{Number of prespecified random splits and their seed.
Ignored if explicit splits are supplied. The caller's RNG state is restored.}
\item{censoring,censor_vars}{Conditional Cox model (default), or marginal
Kaplan--Meier under marginal independent censoring. Optionally specify baseline
column names for the Cox fit separately from the distance predictors. Constant
covariates are excluded; an unobserved censor-model factor level uses the
reference category. With no censoring the censoring survival is one.}
\item{active}{Maximum number of nonzero distance weights per split.}
\item{distance_weights}{Optional named externally fixed nonnegative weights
summing to one. Overrides weight learning (including \code{active}); for
replaying a frozen distance, not for selecting weights on evaluation data.}
\item{rho,delta}{Minimum known-status proportion after subtracting
\eqn{\sqrt{\log(4J/\delta)/(2k)}} for a grid of size \eqn{J}.}
\item{g_floor,omega_max}{Floor used for inverse weighting and maximum allowed
\emph{untruncated} inverse censoring weight used to assess support.}
\item{min_coverage}{Minimum fraction supported in the independent coverage role.}
\item{splits}{Optional list of splits. Each is a named list of four disjoint
integer row-index vectors: \code{construction}, \code{estimation},
\code{coverage}, \code{tuning}. Each has at least two rows. Use explicit
roles for chronological or grouped designs; automatic splitting assumes one
independent row per patient and uses approximately 22.5/52.5/7.5/17.5 percent.
Keep any final test cohort outside training.}
\item{id}{Optional unique baseline tie-breaking ID vector or name of its
column. Defaults to training row order.}
\item{object,x}{A fitted \code{cureknn} object.}
\item{type}{Horizon-specific survival, scalar late survival, or neighbor membership.}
\item{...}{No additional prediction options; reserved for S3 compatibility.}
}
\description{Selects nearest neighbors using characteristics associated with
late survival and estimates survival using inverse censoring weights. Predicts
survival over the specified horizons, not a biological cure probability.}
\details{Construction data determine preprocessing, sparse weights and the
censoring fit for estimation. Neighborhood size minimizes mean IPCW Brier loss
on the same supported tuning patients among candidates meeting the coverage
requirement. Validation weights use a separate censoring fit on the coverage
role. The final horizon determines known status. Each curve is clipped and
projected to be nonincreasing; the scalar clips the average of the raw
estimates. All components must report before their predictions are averaged.
The censoring model, positivity and follow-up assumptions matter for interpretation.}
\value{The fit stores splits, settings, component weights, selected sizes and
tuning diagnostics. Default prediction returns one row per patient, with
\code{reported}, \code{reason}, and one survival column per time. Withholding
produces \code{NA}, not zero. Scalar output uses \code{late_survival}. Neighbor
output contains one row per query, component and neighbor, with the original
training row and ID, rank and reporting flags. Components without an eligible
size have no neighbors. The fitted object and neighbor output are potentially
identifying and should not be distributed as anonymous summaries.}
\examples{
set.seed(32)
n <- 320
d <- data.frame(age = rnorm(n, 60, 8), stage = factor(sample(1:3, n, TRUE)))
t <- ifelse(runif(n) < 0.35, Inf, rexp(n, as.numeric(d$stage)/10))
cens <- runif(n, 10, 20)
d$time <- pmin(t, cens)
d$event <- as.integer(t <= cens)
fit <- cureknn(survival::Surv(time, event) ~ age + stage, d,
               times = 8:10, k = c(25, 50), n_splits = 2, censoring = "marginal")
patient <- data.frame(age = 60, stage = factor(2))
predict(fit, patient)
head(predict(fit, patient, type = "neighbors"))
}
