if (file.exists('R/curve_boost.R')) {
  for (p in list.files('R', pattern='[.]R$', full.names=TRUE)) source(p)
} else {
  library(curelearning)
  .ccrb_state <- getFromNamespace('.ccrb_state', 'curelearning')
  .ccrb_gate_decision <- getFromNamespace('.ccrb_gate_decision', 'curelearning')
  .ccrb_fit <- getFromNamespace('.ccrb_fit', 'curelearning')
}
set.seed(918)
n <- 1000L
x <- runif(n,-1,1)
z <- rnorm(n)
p <- plogis(-.5+2*x)
t <- ifelse(runif(n)<p, Inf, rexp(n,exp(z)/2))
c <- runif(n,5,10)
d <- data.frame(time=pmin(t,c),event=as.integer(t<=c),x=x,z=z)
fit <- curecurveboost(survival::Surv(time,event)~x+z,d,
  breaks=c(0,.5,1,2,4,6,8),nrounds=20L,maxdepth=2L,minbucket=15L,
  forest.control=list(ntree=15L,maxdepth=3L,nodesize=15L,
    min.tail.at.risk=4L,min.tail.end.at.risk=2L,min.events.before.tail=2L,
    tail.start=6,tail.end=8,nsplit=5L),seed=91L)
stopifnot(length(unique(unlist(fit$split)))==n,
          length(unlist(fit$split))==n,
          all(diff(fit$history$train_nll)<=1e-9),
          fit$validation_brier<=fit$history$validation_brier[1]+1e-10)
pred <- predict(fit,d[1:20,],type='all',times=c(0,.5,2,8))
stopifnot(all(is.finite(pred$survival)),all(pred$survival>=0 & pred$survival<=1),
  all(pred$survival[,1]==1),all(apply(pred$survival,1,function(s)all(diff(s)<=1e-12))),
  max(abs(pred$survival-(pred$cure+(1-pred$cure)*pred$susceptible)))<1e-12,
  max(abs(pred$cure-predict(fit$anchor,d[1:20,],type='cure')))<1e-12,
  inherits(try(predict(fit,d[1:2,],times=9),silent=TRUE),'try-error'))
# Analytic mixture-likelihood gradient, including a censored subject.
eta <- matrix(log(c(.3,.4,.5,.6)),2,2)
E <- matrix(c(1,1,1,.4),2,2)
pi <- c(.3,.4)
event <- c(1L,0L)
idx <- c(2L,2L)
s <- .ccrb_state(eta,E,pi,event,idx)
obs <- matrix(0,2,2); obs[1,2] <- 1
g <- (s$w*s$rate*E-obs)/2
for(i in seq_along(eta)) {
  a <- b <- eta; a[i]<-a[i]+1e-6; b[i]<-b[i]-1e-6
  num <- (.ccrb_state(a,E,pi,event,idx)$loss-
            .ccrb_state(b,E,pi,event,idx)$loss)/2e-6
  stopifnot(abs(num-g[i])<1e-6)
}
restored <- unserialize(serialize(fit,NULL))
stopifnot(identical(predict(fit,d[1:10,],times=4),
                    predict(restored,d[1:10,],times=4)))
stopifnot(identical(dim(predict(fit,d[1,,drop=FALSE],times=4)),c(1L,1L)))
# A single held-out paired gate: no-change and weak evidence must be rejected.
stopifnot(!.ccrb_gate_decision(rep(0,100))$accepted,
          !.ccrb_gate_decision(rep(-.01,100))$accepted,
          !.ccrb_gate_decision(rep(c(-.02,.022),50))$accepted,
          .ccrb_gate_decision(rep(c(.02,.04),50))$accepted,
          !.ccrb_gate_decision(rep(.01,100),minimum.gain=.02)$accepted)
dd <- rep(c(-.01,.03),50)
stopifnot(abs(.ccrb_gate_decision(dd)$lower_bound-
  stats::t.test(dd,alternative='greater')$conf.int[1])<1e-12)
# Gate outcomes cannot change training, early stopping, or candidate selection.
args <- list(anchor=fit$anchor,formula=fit$formula,
  data=d[fit$split$correction,],validation.data=d[fit$split$validation,],
  breaks=fit$breaks,score.times=fit$score_times,nrounds=15L,
  maxdepth=2L,minbucket=15L,selection='gated',gate.data=d[fit$split$gate,])
g1 <- do.call(.ccrb_fit,args)
args$gate.data$time <- rev(args$gate.data$time)
args$gate.data$event <- rev(args$gate.data$event)
g2 <- do.call(.ccrb_fit,args)
stopifnot(identical(g1$history,g2$history),
          identical(g1$candidate_nrounds,g2$candidate_nrounds),
          identical(fit$selection,'gated'),length(fit$split$gate)>0L)
if(!fit$gate$accepted)stopifnot(fit$nrounds==0L,length(fit$trees)==0L)
cat('flexible curve correction tests passed\n')
