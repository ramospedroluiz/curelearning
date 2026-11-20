# Public fit/predict interface. No registry variables or filesystem dependencies.
cureknn <- function(formula, data, times, k = NULL, n_splits = 5L, seed = 1L,
                    censoring = c("cox", "marginal"), censor_vars = NULL,
                    active = 2L, rho = 0.20, delta = 0.05, g_floor = 0.02,
                    omega_max = 1/g_floor, min_coverage = 0.80,
                    splits = NULL, id = NULL, distance_weights = NULL) {
  call <- match.call(); censoring <- match.arg(censoring)
  if (!is.data.frame(data) || anyDuplicated(names(data))) stop("Supply a data frame with unique column names.")
  data <- as.data.frame(data); n <- nrow(data)
  if (!inherits(formula,"formula") || length(formula)!=3L) stop("Use Surv(time, event) ~ predictors.")
  if (!is.numeric(times) || !length(times) || any(!is.finite(times) | times<=0) ||
      is.unsorted(times, strictly=TRUE)) stop("times must be positive, finite and strictly increasing.")
  times <- as.numeric(times)
  prob <- function(v, zero, one) length(v)==1L && is.numeric(v) && is.finite(v) &&
    (if (zero) v>=0 else v>0) && (if (one) v<=1 else v<1)
  if (!prob(rho,TRUE,TRUE) || !prob(delta,FALSE,FALSE) || !prob(g_floor,FALSE,TRUE) ||
      !prob(min_coverage,FALSE,TRUE) || length(omega_max)!=1L || !is.numeric(omega_max) ||
      !is.finite(omega_max) || omega_max<1) stop("Invalid coverage or censoring thresholds.")
  if (!.ck_int(active) || length(active)!=1L || !.ck_int(n_splits) || length(n_splits)!=1L ||
      !.ck_int(seed,0L) || length(seed)!=1L) stop("Invalid active, n_splits or seed.")
  if (!is.null(k) && (!.ck_int(k) || anyDuplicated(k))) stop("k must contain distinct positive integers.")
  tt <- stats::terms(formula,data=data)
  vars <- attr(tt,"term.labels")
  if (!length(vars) || any(!vars %in% names(data)) || any(attr(tt,"order")!=1L) ||
      length(attr(tt,"offset"))) stop("Use additive, untransformed predictor columns.")
  outcome <- all.vars(formula[[2L]])
  if (!is.null(distance_weights)) {
    if (!is.numeric(distance_weights) || !setequal(names(distance_weights),vars) ||
        anyDuplicated(names(distance_weights)) || any(!is.finite(distance_weights)) ||
        any(distance_weights<0) || abs(sum(distance_weights)-1)>1e-10)
      stop("External distance_weights must be named nonnegative finite weights summing to one.")
    distance_weights <- distance_weights[vars]
  }
  if (!all(outcome %in% names(data))) stop("Outcome variables must be columns of data.")
  if (is.null(censor_vars)) censor_vars <- vars
  if (!is.character(censor_vars) || anyDuplicated(censor_vars) || any(!censor_vars %in% names(data)) ||
      any(c(vars,censor_vars) %in% outcome)) stop("Invalid predictors or outcome leakage.")
  if (length(id)==1L && is.character(id) && id %in% names(data)) {
    if (id %in% c(vars,censor_vars,outcome)) stop("Do not use the ID as a predictor or response.")
    id <- data[[id]]
  }
  if (is.null(id)) id <- seq_len(n)
  if ((!is.numeric(id) && !is.character(id)) || length(id)!=n || anyNA(id) ||
      anyDuplicated(id) || (is.numeric(id) && any(!is.finite(id)))) stop("IDs must be unique, complete numeric or character values.")
  y <- stats::model.response(stats::model.frame(tt,data,na.action=stats::na.pass))
  if (!inherits(y,"Surv") || attr(y,"type")!="right" || anyNA(y) ||
      any(!is.finite(y)) || any(y[,1L]<0)) stop("Complete, finite right-censored Surv responses are required.")
  time <- y[,1L]; event <- y[,2L]
  roles <- c("construction","estimation","coverage","tuning")
  if (is.null(splits)) {
    if (n<40L) stop("Automatic splitting needs at least 40 rows; supply more training data.")
    old <- if (exists(".Random.seed", envir=.GlobalEnv,inherits=FALSE)) get(".Random.seed",envir=.GlobalEnv)
    on.exit(if (is.null(old)) {if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv)} else assign(".Random.seed",old,envir=.GlobalEnv),add=TRUE)
    set.seed(seed)
    sizes <- pmax(2L,round(n*c(.225,.525,.075))); sizes <- c(sizes,n-sum(sizes))
    splits <- lapply(seq_len(n_splits),function(i)
      stats::setNames(lapply(split(sample.int(n),rep(1:4,sizes)),sort),roles))
  }
  if (!is.list(splits) || !length(splits)) stop("splits must be a nonempty list of four-role splits.")
  for (s in splits) {
    if (!is.list(s) || length(s)!=4L || !setequal(names(s), roles) ||
        !all(vapply(s,function(i) .ck_int(i) && all(i<=n) && length(i)>=2L,logical(1))) ||
        anyDuplicated(unlist(s,use.names=FALSE))) stop("Split roles must have disjoint, valid row indices (at least two per role).")
  }
  cap <- min(vapply(splits,function(s) length(s$estimation),integer(1)))
  if (is.null(k)) k <- unique(pmin(cap,c(25L,50L,100L,200L,400L)))
  if (max(k)>cap) stop("k exceeds the smallest estimation sample.")
  k <- sort(as.integer(k))
  fit_one <- function(s) {
    m <- s$construction; e <- s$estimation; c0 <- s$coverage; v <- s$tuning
    allvars <- unique(c(vars,if(censoring=="cox") censor_vars))
    specs <- lapply(data[m,allvars,drop=FALSE],.ck_spec)
    prep <- function(ix) .ck_frame(data[ix,,drop=FALSE],specs)
    M <- prep(m); E <- prep(e); C <- prep(c0); V <- prep(v)
    cv <- if (censoring=="cox") censor_vars else character()
    gc <- .ck_gfit(M[cv],time[m],event[m],censoring)
    gv <- .ck_gfit(C[cv],time[c0],event[c0],censoring)
    pseudo <- function(X,i,g) outer(time[i],times,">") / pmax(g,g_floor)
    gm <- .ck_gmat(gc,M[cv],times); ge <- .ck_gmat(gc,E[cv],times)
    yy <- rowMeans(pseudo(M,m,gm))
    w <- vapply(M[vars],.ck_eta,numeric(1),y=yy); w <- pmax(w,1e-4)
    # Match normalization at learning and at sparsification in the reference.
    w <- w/sum(w); w <- w/sum(w)
    keep <- order(w,decreasing=TRUE)[seq_len(min(active,length(w)))]
    w[-keep] <- 0; w <- w/sum(w)
    if (!is.null(distance_weights)) w <- distance_weights
    design <- function(x) .ck_distance_x(x,specs,w)
    component <- list(specs=specs[vars],weights=w,train=design(E),id=id[e],row=e,
      pseudo=pseudo(E,e,ge), omega=apply(1/ge,1L,max),
      known=event[e]==1 & time[e]<=max(times) | time[e]>max(times))
    query <- rbind(design(C),design(V))
    ix <- .ck_neighbors(component$train,query,k,component$id)
    calc <- lapply(k,function(kk) .ck_local(component,ix,kk,length(k),rho,delta,omega_max))
    ci <- seq_along(c0); vi <- length(c0)+seq_along(v)
    coverage <- vapply(calc,function(z) mean(z$ok[ci]),numeric(1))
    eligible <- coverage>=min_coverage
    common <- if (any(eligible)) Reduce("&",lapply(calc[eligible],function(z) z$ok[vi])) else rep(FALSE,length(v))
    risk <- rep(NA_real_,length(k)); reason <- "reported"
    if (!any(eligible)) reason <- "no_eligible_k"
    else if (!any(common)) reason <- "no_common_tuning_support"
    else {
      vv <- v[common]; X <- V[common,cv,drop=FALSE]
      alive <- outer(time[vv],times,">")
      died <- (!alive) & (event[vv]==1L)
      gtime <- pmax(.ck_gmat(gv,X,times),g_floor)
      gdeath <- pmax(.ck_g(gv,X,time[vv],left=TRUE),g_floor)
      for (j in which(eligible)) {
        pred <- calc[[j]]$curve[vi[common],,drop=FALSE]
        risk[j] <- mean(alive*(1-pred)^2/gtime + died*pred^2/gdeath)
      }
    }
    component$k <- if (any(is.finite(risk))) k[which.min(risk)] else NA_integer_
    component$reason <- reason
    component$tuning <- data.frame(k=k,coverage=coverage,eligible=eligible,
                                   n_common=sum(common),brier=risk)
    component
  }
  structure(list(call=call,times=times,vars=vars,k_grid=k,splits=splits,
    settings=list(rho=rho,delta=delta,omega_max=omega_max,g_floor=g_floor,
                  min_coverage=min_coverage,censoring=censoring,censor_vars=censor_vars,
                  active=active,seed=seed,distance_weights=distance_weights),
    components=lapply(splits,fit_one)),class="cureknn")
}

predict.cureknn <- function(object,newdata,type=c("survival","scalar","neighbors"),...) {
  type <- match.arg(type)
  if (length(list(...))) stop("Unexpected predict argument. Predictions are restricted to the fitted times.")
  if (missing(newdata) || !is.data.frame(newdata) || !nrow(newdata) ||
      anyDuplicated(names(newdata)) || !all(object$vars %in% names(newdata)))
    stop("Supply nonempty newdata with all baseline predictors.")
  newdata <- as.data.frame(newdata); n <- nrow(newdata); h <- length(object$times)
  S <- object$settings; total <- matrix(0,n,h); tail <- numeric(n); reason <- rep("reported",n)
  all_neighbors <- list()
  for (b in seq_along(object$components)) {
    cp <- object$components[[b]]
    if (is.na(cp$k)) {
      reason[reason=="reported"] <- cp$reason
      next
    }
    q <- .ck_distance_x(.ck_frame(newdata,cp$specs),cp$specs,cp$weights)
    ix <- .ck_neighbors(cp$train,q,cp$k,cp$id)
    z <- .ck_local(cp,ix,cp$k,length(object$k_grid),S$rho,S$delta,S$omega_max)
    fail <- !z$ok & reason=="reported"
    reason[fail] <- ifelse(z$lower[fail]<S$rho,"insufficient_followup","unstable_censoring_weights")
    if (type=="neighbors") {
      all_neighbors[[length(all_neighbors)+1L]] <- data.frame(query=rep(seq_len(n),each=cp$k),
        split=b,rank=rep(seq_len(cp$k),n),training_row=cp$row[as.vector(t(ix))],
        id=cp$id[as.vector(t(ix))],component_supported=rep(z$ok,each=cp$k))
    } else {total <- total+z$curve; tail <- tail+z$scalar}
  }
  ok <- reason=="reported"; B <- length(object$components)
  if (type=="neighbors") {
    ans <- if (length(all_neighbors)) do.call(rbind,all_neighbors) else
      data.frame(query=integer(),split=integer(),rank=integer(),training_row=integer(),
                 id=object$components[[1]]$id[FALSE],component_supported=logical())
    ans$reported <- ok[ans$query]; return(ans)
  }
  ans <- data.frame(reported=ok,reason=reason,row.names=rownames(newdata))
  if (type=="scalar") {ans$late_survival <- ifelse(ok,tail/B,NA_real_); return(ans)}
  total <- total/B; total[!ok,] <- NA_real_
  colnames(total) <- paste0("S_",as.character(object$times))
  cbind(ans,as.data.frame(total,check.names=FALSE))
}

print.cureknn <- function(x,...) {
  cat("Late-survival nearest neighbors (cureknn)\nHorizons:",x$times,
      "\nSelected k by split:",vapply(x$components,function(z) z$k,integer(1)),"\n")
  invisible(x)
}
