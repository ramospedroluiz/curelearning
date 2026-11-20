# Internal primitives for cureknn; all preprocessing objects are fit on one role.
.ck_int <- function(x, lo = 1L) {
  is.numeric(x) && length(x) > 0L && all(is.finite(x) & x == floor(x) &
    x >= lo & x <= .Machine$integer.max)
}

.ck_spec <- function(x) {
  if (is.numeric(x) && !is.object(x) && is.null(dim(x))) {
    good <- x[is.finite(x)]
    if (!length(good)) stop("A numeric predictor has no finite construction values.")
    m <- stats::median(good); s <- stats::IQR(good) / 1.349
    if (!is.finite(s) || s <= 0) s <- stats::sd(good)
    return(list(num = TRUE, center = m, scale = if (is.finite(s) && s > 0) s else 1))
  }
  if (!(is.factor(x) || is.character(x) || is.logical(x)) || !is.null(dim(x)))
    stop("Predictors must be numeric vectors, factors, logicals or characters.")
  v <- as.character(x); v[is.na(v) | !nzchar(v)] <- "unknown"
  n <- table(v); keep <- names(n)[n >= max(100, ceiling(.0025 * length(x)))]
  if (length(keep) < min(2L, length(n)))
    keep <- names(sort(n, decreasing = TRUE))[seq_len(min(2L, length(n)))]
  list(num = FALSE, levels = unique(c(sort(keep), "other")))
}

.ck_apply <- function(x, s) {
  if (s$num) {
    if (!is.numeric(x) || is.object(x) || !is.null(dim(x))) stop("Numeric predictor type changed.")
    x[!is.finite(x)] <- s$center
    return(x)
  }
  if (!(is.factor(x) || is.character(x) || is.logical(x))) stop("Categorical predictor type changed.")
  v <- as.character(x); v[is.na(v) | !nzchar(v)] <- "unknown"
  v[!v %in% s$levels] <- "other"
  factor(v, levels = s$levels)
}

.ck_frame <- function(x, specs) {
  if (!all(names(specs) %in% names(x))) stop("Missing predictor(s): ",
    paste(setdiff(names(specs), names(x)), collapse = ", "))
  as.data.frame(Map(.ck_apply, x[names(specs)], specs), optional = TRUE,
                stringsAsFactors = FALSE, check.names = FALSE)
}

.ck_eta <- function(x, y) {
  if (length(y) < 20L || stats::var(y) <= 0) return(0)
  if (is.numeric(x) && length(unique(x)) > 10L) {
    b <- unique(stats::quantile(x, 0:10/10, names = FALSE))
    if (length(b) < 3L) return(0)
    group <- cut(x, b, include.lowest = TRUE)
  } else group <- factor(x)
  eta <- stats::var(stats::ave(y, group, FUN = mean)) / stats::var(y)
  null <- (nlevels(group) - 1) / (length(y) - 1)
  max(0, (eta - null) / max(1e-12, 1 - null))
}

.ck_distance_x <- function(x, specs, weights) {
  # Retain all columns (including zero-weight columns) as in the reference
  # implementation, for reproducible floating point boundary ties.
  do.call(cbind, lapply(names(weights), function(n) {
    v <- x[[n]]; s <- specs[[n]]; w <- weights[[n]]
    if (s$num) return(matrix(sqrt(w) * (v - s$center) / s$scale, ncol = 1L))
    m <- matrix(0, length(v), length(s$levels))
    m[cbind(seq_along(v), as.integer(v))] <- sqrt(w/2)
    m
  }))
}

.ck_neighbors <- function(train, query, grid, id) {
  k <- max(grid); h <- min(k+1L, nrow(train))
  z <- FNN::get.knnx(train, query, k = h, algorithm = "kd_tree")
  bounds <- grid[grid < h]; sq <- rowSums(train^2)
  ix <- matrix(0L, nrow(query), k)
  for (i in seq_len(nrow(query))) {
    d <- z$nn.dist[i, ]
    if (any(abs(d[bounds]-d[bounds+1L]) <= 100*.Machine$double.eps *
        pmax(1, d[bounds], d[bounds+1L]))) {
      ds <- pmax(0, sq + sum(query[i, ]^2) - 2*as.numeric(train %*% query[i, ]))
      ix[i, ] <- order(ds, id, method = "radix")[seq_len(k)]
    } else {
      cand <- z$nn.index[i, seq_len(k)]
      ix[i, ] <- cand[order(d[seq_len(k)], id[cand], method = "radix")]
    }
  }
  ix
}

.ck_gfit <- function(x, time, event, method) {
  if (!any(event == 0L)) return(list(method = "none"))
  keep <- names(x)[vapply(x, function(v) length(unique(v)) > 1L, logical(1))]
  if (method == "marginal" || !length(keep)) {
    fit <- survival::survfit(survival::Surv(time, 1-event) ~ 1)
    return(list(method = "km", time = fit$time, surv = fit$surv))
  }
  d <- x[keep]; levels_fit <- lapply(d, function(z) if (is.factor(z)) levels(droplevels(z)))
  for (n in keep) if (is.factor(d[[n]])) d[[n]] <- droplevels(d[[n]])
  names(d) <- paste0("v", seq_along(d))
  f <- stats::reformulate(names(d), response="survival::Surv(.t, .d)")
  d$.t <- time; d$.d <- 1-event
  fit <- withCallingHandlers(survival::coxph(f, data=d, ties="efron", model=TRUE),
    warning=function(w) {
      if (grepl("did not converge|infinite|ran out", conditionMessage(w), ignore.case=TRUE))
        stop("Unstable Cox censoring fit; revise its predictors or explicitly use marginal censoring.")
    })
  bh <- survival::basehaz(fit, centered=TRUE)
  list(method="cox", fit=fit, keep=keep, levels=levels_fit, time=bh$time, hazard=bh$hazard)
}

.ck_g <- function(g, x, t, left=FALSE) {
  if (length(t)==1L) t <- rep(t, nrow(x))
  if (g$method=="none") return(rep(1, length(t)))
  # Include censoring at t for survivors and exclude censoring at t for deaths.
  idx <- findInterval(t, g$time, left.open=left)
  if (g$method=="km") return(c(1,g$surv)[idx+1L])
  d <- x[g$keep]
  for (n in g$keep) if (length(g$levels[[n]])) {
    v <- as.character(d[[n]]); lev <- g$levels[[n]]
    v[!v %in% lev] <- lev[1L]; d[[n]] <- factor(v,levels=lev)
  }
  names(d) <- paste0("v", seq_along(d))
  lp <- as.numeric(stats::predict(g$fit, newdata=d,type="lp",reference="sample"))
  h <- c(0,g$hazard)[idx+1L]
  result <- rep(1,length(t)); pos <- h>0
  result[pos] <- exp(-h[pos]*exp(lp[pos]))
  if (anyNA(result)) stop("Nonfinite censoring-model prediction.")
  result
}

.ck_gmat <- function(g,x,times) {
  matrix(unlist(lapply(times,function(t) .ck_g(g,x,t)),use.names=FALSE),
         nrow=nrow(x),ncol=length(times))
}

.ck_curve <- function(raw) {
  z <- pmin(pmax(raw,0),1)
  if (ncol(z)>1L) for (i in seq_len(nrow(z))) {
    if (is.unsorted(-z[i,])) z[i,] <- -stats::isoreg(-z[i,])$yf
  }
  z
}

.ck_local <- function(component, ix, k, J, rho, delta, omega.max) {
  ix <- ix[,seq_len(k),drop=FALSE]; n <- nrow(ix)
  colmean <- function(v) rowMeans(matrix(v[ix],n,k))
  raw <- matrix(unlist(lapply(seq_len(ncol(component$pseudo)),function(j)
    colmean(component$pseudo[,j])),use.names=FALSE),nrow=n)
  known <- colmean(component$known)
  omega <- apply(matrix(component$omega[ix],n,k),1L,max)
  lower <- known - sqrt(log(4*J/delta)/(2*k))
  list(curve=.ck_curve(raw), scalar=pmin(1,pmax(0,rowMeans(raw))), known=known,
       lower=lower,omega=omega,ok=lower>=rho & omega<=omega.max)
}
