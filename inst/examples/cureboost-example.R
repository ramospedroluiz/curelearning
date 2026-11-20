library(curelearning)
library(survival)

data(cure_latency_demo)
train <- subset(cure_latency_demo, cohort == "train")
test <- subset(cure_latency_demo, cohort == "test")

form <- Surv(time, event) ~
  cure_marker + latency_marker + age + stage + biomarker + noise

fit.boost <- cureboost(
  form,
  train,
  nrounds = 80,
  learning.rate = 0.05,
  maxdepth = 1,
  mtry = 6,
  nodesize = 80,
  nsplit = 10,
  tail.start = 7,
  tail.end = 9.5,
  min.tail.at.risk = 20,
  min.tail.end.at.risk = 5,
  ipcw.model = "marginal",
  validation.fraction = 0.20,
  early.stopping.rounds = 15,
  seed = 2026
)

pred.boost <- predict(fit.boost, test, type = "cure")

print(fit.boost)
print(summary(fit.boost))

cat("\nFirst finite-tail predictions:\n")
print(head(pred.boost))

cat("\nMean held-out prediction by the true cure marker:\n")
print(tapply(pred.boost, test$cure_marker, mean))

cat("\nBoosting variable importance:\n")
print(cureboost_importance(fit.boost))

selected <- subset(fit.boost$history, selected)
cat("\nTraining loss was monotone:",
    all(diff(selected$training_risk) <= 1e-12), "\n")
