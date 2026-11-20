library(curelearning)
set.seed(41)
n <- 440
d <- data.frame(x=rnorm(n), site=factor(sample(letters[1:3],n,TRUE)),
                key=sample.int(n))
d$time <- sample(1:14,n,TRUE); d$event <- 1L
f <- survival::Surv(time,event) ~ x + site
s <- list(construction=1:100,estimation=101:260,coverage=261:320,tuning=321:400)
fit <- function(dat=d, ts=8:10, ...) cureknn(f,dat,ts,k=c(25,50),splits=list(s),
                                              id="key",...)
m <- fit(); q <- d[401:410,c("x","site")]
p <- predict(m,q); nn <- predict(m,q,"neighbors"); sc <- predict(m,q,"scalar")
stopifnot(all(p$reported),all(p$reason=="reported"),nrow(p)==10,
          all(sc$late_survival>=0 & sc$late_survival<=1))
for (i in seq_len(nrow(q))) {
  rows <- nn$training_row[nn$query==i]
  vals <- vapply(8:10,function(t) mean(d$time[rows]>t),numeric(1))
  stopifnot(all(rows %in% s$estimation),all(nn$id[nn$query==i]==d$key[rows]),
            max(abs(vals-unlist(p[i,3:5])))<1e-12,abs(mean(vals)-sc$late_survival[i])<1e-12)
}
for (ts in list(8,8:9,6:10)) {
  pp <- predict(fit(ts=ts),q[1,,drop=FALSE])
  stopifnot(ncol(pp)==length(ts)+2L, all(diff(unlist(pp[1,-(1:2)]))<=1e-14))
}
old <- .Random.seed
a <- cureknn(f,d,8:10,k=25,n_splits=2,seed=511)
stopifnot(identical(.Random.seed,old),length(a$splits)==2L)
b <- cureknn(f,d,8:10,k=25,n_splits=2,seed=511)
stopifnot(identical(predict(a,q),predict(b,q)))
for (s0 in a$splits) stopifnot(!anyDuplicated(unlist(s0)))
newq <- q[1,,drop=FALSE];newq$x <- NA_real_;newq$site <- "not-in-training"
stopifnot(nrow(predict(a,newq))==1L)
# No eligible candidate: keep an auditable fit and never emit false zeros.
dc <- d; dc$time <- 1; dc$event <- 0
no <- fit(dc,censoring="marginal")
pp <- predict(no,q);stopifnot(all(!pp$reported),all(pp$reason=="no_eligible_k"),
  all(is.na(pp$S_8)),nrow(predict(no,q,"neighbors"))==0L)
# Clipping before shape projection is not the same as clipping the scalar once.
z <- curelearning:::.ck_local(list(pseudo=matrix(c(1.8,.2,.4),1),known=TRUE,omega=2),
                            matrix(1L,1,1),1,1,0,.95,50)
stopifnot(max(abs(z$curve-c(1,.3,.3)))<1e-12,abs(z$scalar-.8)<1e-12)
# An event-time tie uses the left limit only for an observed event.
g <- list(method="km",time=c(2,3),surv=c(.8,.4))
stopifnot(identical(curelearning:::.ck_g(g,q[1:2,],c(2,3),TRUE),c(1,.8)))
# One unsupported component vetoes the whole ensemble.
veto <- a; veto$components[[2]]$k <- NA_integer_; veto$components[[2]]$reason <- "no_common_tuning_support"
stopifnot(all(is.na(predict(veto,q)$S_8)),all(!predict(veto,q,"neighbors")$reported))
# Baseline tie ID, not event status, determines neighbors at equal distances.
dt <- d;dt$x<-1;dt$site<-"all"
tfit <- fit(dt); nnt <- predict(tfit,dt[1,,drop=FALSE],"neighbors")
stopifnot(identical(nnt$id,head(sort(dt$key[s$estimation]),tfit$components[[1]]$k)))
fails <- function(expr) inherits(try(force(expr),silent=TRUE),"try-error")
stopifnot(fails(fit(ts=c(9,8))),fails(fit(ts=c(8,8))),fails(predict(a,q["x"])),
          fails(cureknn(f,d,8,k=2000)),fails(cureknn(f,d,8,id=rep(1,n))),
          fails(cureknn(f,d,8,splits=list(within(s,{tuning<-coverage})))),
          fails(cureknn(f,d,8,censor_vars="event")))
cat("cureknn public API and edge cases passed.\n")
shift <- d
shift$site <- rep("A",n);shift$site[101:140] <- "B"
shift$site[1:40] <- "B"; shift$time[101:140] <- 1;shift$event[101:140] <- 0
local <- cureknn(f,shift,8:10,k=25,splits=list(s),id="key",censoring="marginal",
                 distance_weights=c(x=0,site=1))
st <- predict(local,data.frame(x=0,site=c("A","B")))
stopifnot(identical(st$reported,c(TRUE,FALSE)),is.na(st$S_10[2]),
          st$reason[2]=="insufficient_followup")
# The raw (not floored) censoring weights trigger an independent withholding rule.
cp <- m$components[[1]]; cp$omega[] <- 60
x0 <- m;x0$components[[1]] <- cp
stopifnot(all(predict(x0,q)$reason=="unstable_censoring_weights"))
stopifnot(fails(cureknn(f,d,8,distance_weights=c(site=2,x=-1))),
          fails(cureknn(f,d,8,distance_weights=c(0,1))))
