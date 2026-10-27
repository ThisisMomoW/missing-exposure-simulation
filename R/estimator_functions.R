# AIPW, imputation, package drcmd and TR-AIPW estimators for families 1–5.
# Propensity/missingness bounds are 0.025/0.975 for AIPW and drcmd.
# Interaction outcome models use A:X terms; treatment models use main effects.
# TR-AIPW uses inverse-weight shrinkage and percentile bootstrap intervals.

suppressPackageStartupMessages({
  library(AIPW)
  library(SuperLearner)
  library(drcmd)
  library(mice)
  library(rootSolve)
  library(parallel)
  library(doParallel)
  library(foreach)
})


DEFAULT_TRIM <- 0.025
NUMERICAL_EPS <- 1e-6

bound_probability <- function(p, eps = NUMERICAL_EPS) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

empty_result <- function() {
  data.frame(estimate=NA_real_, se=NA_real_, ci_lower=NA_real_, ci_upper=NA_real_)
}

# TR-AIPW successful fits return additional bootstrap and weight diagnostics.
# Failed fits must return the same schema so that the 500 Monte Carlo rows can
# always be combined with rbind(), even when some replicates fail validity or
# bootstrap-success checks.
empty_trme_result <- function() {
  data.frame(
    estimate=NA_real_,
    se=NA_real_,
    ci_lower=NA_real_,
    ci_upper=NA_real_,
    n_boot_valid=NA_integer_,
    n_boot_invalid=NA_integer_,
    boot_median=NA_real_,
    boot_q025=NA_real_,
    boot_q975=NA_real_,
    boot_min=NA_real_,
    boot_max=NA_real_,
    q01_treated_ps=NA_real_,
    q99_treated_inv_weight=NA_real_,
    q99_control_inv_weight=NA_real_,
    q99_missingness_weight=NA_real_,
    max_missingness_weight=NA_real_
  )
}

safe_result <- function(expr, empty=empty_result()) {
  tryCatch(expr, error=function(e) {
    warning(conditionMessage(e)); empty
  })
}

# DGP-informed parametric SuperLearner wrapper.
#
# For the outcome regression Q(A, X), the true interaction DGP contains only
# treatment-by-covariate interactions:
#   A:X1, A:X2, and A:X3. with the same coefficient.
#
# The wrapper therefore fits Y ~ A + X1 + X2 + X3 + A:I(X1 + X2 + X3) when A is present in the
# learner design matrix. For nuisance models whose design matrix does not
# contain A (for example, the treatment propensity model), it falls back to a
# main-effects GLM rather than adding unnecessary X:X interactions.
SL.glm.AX <- function(
    Y,
    X,
    newX,
    family,
    obsWeights,
    id,
    ...
) {
  X <- as.data.frame(X)
  newX <- as.data.frame(newX)

  dat <- data.frame(
    Y = Y,
    X,
    check.names = FALSE
  )

  if (all(c("A", "X1", "X2", "X3") %in% names(X))) {
    fit_formula <- Y ~ A + X1 + X2 + X3 + A:I(X1 + X2 + X3)
  } else {
    fit_formula <- Y ~ .
  }

  fit_glm <- glm(
    fit_formula,
    data = dat,
    family = family,
    weights = obsWeights
  )

  pred <- predict(
    fit_glm,
    newdata = newX,
    type = "response"
  )

  fit_object <- list(
    object = fit_glm
  )
  class(fit_object) <- "SL.glm.AX"

  list(
    pred = as.numeric(pred),
    fit = fit_object
  )
}


predict.SL.glm.AX <- function(
    object,
    newdata,
    ...
) {
  newdata <- as.data.frame(newdata)

  pred <- predict(
    object$object,
    newdata = newdata,
    type = "response"
  )

  as.numeric(pred)
}


q_sl_library <- function(nuisance_spec = c("main", "interaction")) {
  nuisance_spec <- match.arg(nuisance_spec)

  if (nuisance_spec == "main") {
    "SL.glm"
  } else {
    "SL.glm.AX"
  }
}


g_sl_library <- function(nuisance_spec = c("main", "interaction")) {
  nuisance_spec <- match.arg(nuisance_spec)

  # The exposure DGP contains only main effects, so g(X) remains a main-effects
  # logistic model in both specifications. This isolates the comparison to the
  # outcome nuisance model.
  "SL.glm"
}


run_aipw_result <- function(df, nuisance_spec="main", trim=DEFAULT_TRIM) {
  q_learner <- q_sl_library(nuisance_spec)
  g_learner <- g_sl_library(nuisance_spec)

  invisible(capture.output(
    fit <- AIPW::AIPW$new(
      Y=df$Y, A=df$A, W=df[,c("X1","X2","X3")],
      Q.SL.library=q_learner,
      g.SL.library=g_learner
    )$fit()
  ))
  # g.bound is applied when the estimand and influence-curve SE are computed.
  # Capture the printed table once; repeated summary() calls create very large
  # Slurm logs and materially slow the Monte Carlo jobs.
  summary_text <- capture.output(
    sm <- fit$summary(g.bound=trim)
  )

  # Prefer structured output. Fall back to parsing printed output for package
  # versions whose summary method invisibly returns the R6 object.
  candidates <- list(sm, fit$result, fit$results)
  for (obj in candidates) {
    if (is.data.frame(obj) || is.matrix(obj)) {
      rn <- rownames(obj)
      hit <- grep("Risk Difference", rn, fixed=TRUE)
      if (length(hit)) {
        z <- as.numeric(obj[hit[1],])
        z <- z[is.finite(z)]
        if (length(z) >= 2) {
          return(data.frame(estimate=z[1], se=z[2],
                            ci_lower=if(length(z)>=3) z[3] else z[1]-1.96*z[2],
                            ci_upper=if(length(z)>=4) z[4] else z[1]+1.96*z[2]))
        }
      }
    }
  }
  line <- summary_text[grepl("Risk Difference",summary_text,fixed=TRUE)][1]
  if (is.na(line)) stop("Risk Difference row not found in AIPW summary")
  nums <- as.numeric(regmatches(line,gregexpr("[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?",line,perl=TRUE))[[1]])
  out <- empty_result(); out$estimate <- nums[1]
  if(length(nums)>=2) out$se <- nums[2]
  if(length(nums)>=4) {out$ci_lower<-nums[3]; out$ci_upper<-nums[4]}
  else if(is.finite(out$se)) {out$ci_lower<-out$estimate-1.96*out$se; out$ci_upper<-out$estimate+1.96*out$se}
  out
}

run_single_imp_result <- function(df_missing, nuisance_spec="main", trim=DEFAULT_TRIM) {
  cc <- df_missing[!is.na(df_missing$A),]

  # Under the interaction DGP, the exposure-imputation model is allowed to
  # depend on Y:X terms. We do not add unrelated X:X interactions.
  imp_formula <- if (nuisance_spec == "main") {
    A ~ X1 + X2 + X3 + Y
  } else {
    A ~ X1 + X2 + X3 + Y * (X1 + X2 + X3)
  }

  fit <- glm(imp_formula,data=cc,family=binomial())
  ph <- bound_probability(predict(fit,newdata=df_missing,type="response"))
  d <- df_missing; mis <- is.na(d$A)
  d$A[mis] <- rbinom(sum(mis),1,ph[mis])
  run_aipw_result(d,nuisance_spec,trim)
}

pool_rubin <- function(q, u, conf_level = 0.95) {
  keep <- is.finite(q) & is.finite(u) & u >= 0
  q <- q[keep]
  u <- u[keep]
  m <- length(q)

  if (m < 2L) {
    return(data.frame(estimate=if (m == 1L) q else NA_real_,
                      se=NA_real_, ci_lower=NA_real_, ci_upper=NA_real_))
  }

  qbar <- mean(q)
  ubar <- mean(u)
  b <- stats::var(q)
  tvar <- ubar + (1 + 1 / m) * b
  se <- sqrt(tvar)

  # Rubin large-sample degrees of freedom. When between-imputation variance is
  # essentially zero, the normal approximation is recovered.
  if (!is.finite(b) || b <= .Machine$double.eps) {
    df <- Inf
  } else {
    r <- ((1 + 1 / m) * b) / ubar
    df <- (m - 1) * (1 + 1 / r)^2
  }
  alpha <- 1 - conf_level
  crit <- if (is.finite(df)) stats::qt(1 - alpha / 2, df=df) else stats::qnorm(1 - alpha / 2)

  data.frame(estimate=qbar, se=se,
             ci_lower=qbar-crit*se, ci_upper=qbar+crit*se)
}

run_mice_result <- function(df_missing,m=10L,maxit=10L,seed=1L,nuisance_spec="main",trim=DEFAULT_TRIM) {
  d <- df_missing[,c("A","Y","X1","X2","X3")]

  if (nuisance_spec == "interaction") {
    # DGP-informed predictors for P(A | X, Y). These terms allow the exposure
    # imputation model to represent Y-by-covariate interaction without adding
    # unrelated X-by-X interactions.
    d$YX1 <- d$Y * d$X1
    d$YX2 <- d$Y * d$X2
    d$YX3 <- d$Y * d$X3
  }

  d$A <- factor(d$A,levels=c(0,1))

  meth <- mice::make.method(d)
  meth[] <- ""
  meth["A"] <- "logreg"

  pm <- mice::make.predictorMatrix(d)
  pm[,] <- 0

  predictors_A <- c("Y","X1","X2","X3")
  if (nuisance_spec == "interaction") {
    predictors_A <- c(predictors_A,"YX1","YX2","YX3")
  }
  pm["A",predictors_A] <- 1

  imp <- mice::mice(
    d,
    m=m,
    maxit=maxit,
    method=meth,
    predictorMatrix=pm,
    seed=seed,
    printFlag=FALSE
  )

  rr <- lapply(seq_len(m),function(k){
    z <- mice::complete(imp,k)
    z$A <- as.integer(as.character(z$A))
    run_aipw_result(z,nuisance_spec,trim)
  })

  q <- vapply(rr,`[[`,numeric(1),"estimate")
  se <- vapply(rr,`[[`,numeric(1),"se")

  if(any(!is.finite(se))) {
    return(data.frame(
      estimate=mean(q),
      se=NA,
      ci_lower=NA,
      ci_upper=NA
    ))
  }

  pool_rubin(q,se^2)
}

run_drcmd_result <- function(df_missing,nuisance_spec="main",trim=DEFAULT_TRIM) {
  learner <- if (nuisance_spec == "main") "SL.glm" else "SL.glm.AX"
  fit <- drcmd::drcmd(Y=df_missing$Y,A=df_missing$A,
                     X=df_missing[,c("X1","X2","X3")],
                     default_learners=learner, cutoff=trim, k=10L)
  est_obj <- fit$results$estimates
  se_obj  <- fit$results$ses
  est <- if(is.data.frame(est_obj)) est_obj$psi_hat_ate[1] else unlist(est_obj)["psi_hat_ate"]
  se  <- if(is.data.frame(se_obj))  se_obj$psi_hat_ate[1]  else unlist(se_obj)["psi_hat_ate"]
  if(!is.finite(est)) stop("drcmd ATE estimate was not found")
  if(!is.finite(se)) stop("drcmd ATE SE was not found in fit$results$ses")
  data.frame(estimate=as.numeric(est),se=as.numeric(se),
             ci_lower=as.numeric(est)-1.96*as.numeric(se),
             ci_upper=as.numeric(est)+1.96*as.numeric(se))
}

run_trme_rd_result <- function(df_missing, nuisance_spec="main",
                               missing_mechanism="MAR", shrink_rate=0.99,
                               B=100L, seed=1L) {
  # trme defines r=1 as exposure missing and sets it using is.na(A).
  # The simulation uses R=1 for observed exposure, but A already contains NA
  # exactly where exposure is missing, so no manual R recoding is needed.
  if (!is.numeric(shrink_rate) || length(shrink_rate) != 1L ||
      !is.finite(shrink_rate) || shrink_rate <= 0 || shrink_rate > 1) {
    stop("shrink_rate must be a finite scalar in (0, 1].")
  }
  if (identical(shrink_rate, 1)) {
    warning("TR-AIPW shrink_rate=1 applies no inverse-weight shrinkage.",
            call.=FALSE)
  }

  fit <- trme_rd(
    covs=c("X1","X2","X3"), Y="Y", A="A", data=df_missing,
    imp_model=TRUE, shrink_rate=shrink_rate, ci_alpha=0.95,
    method="TR-AIPW", bootstrap=TRUE, B=B, seed=seed,
    nuisance_spec=nuisance_spec,
    missing_mechanism=missing_mechanism
  )

  # The risk-difference point estimate is exactly tau1 - tau0. Prefer the
  # package-level Estimate when available, but verify that it agrees with the
  # two intervention-risk components.
  tau1 <- as.numeric(fit$tau1)
  tau0 <- as.numeric(fit$tau0)
  estimate_from_tau <- tau1 - tau0
  estimate_returned <- suppressWarnings(as.numeric(fit$Estimate))
  if (length(estimate_returned) != 1L || !is.finite(estimate_returned)) {
    estimate_returned <- estimate_from_tau
  }
  if (!is.finite(estimate_from_tau)) {
    stop("trme TR-AIPW intervention risks were not finite.")
  }
  if (abs(estimate_returned - estimate_from_tau) > 1e-8) {
    stop("trme returned Estimate does not equal tau1 - tau0 on the RD scale.")
  }
  estimate <- estimate_from_tau

  ## trme_rd() has already excluded failed or theoretically invalid
  ## bootstrap roots (tau0/tau1 outside [0,1], RD outside [-1,1], or
  ## inconsistent RD != tau1-tau0).
  boot <- as.numeric(fit$boot_est)
  boot <- boot[is.finite(boot)]
  min_boot <- max(20L, ceiling(0.8 * B))

  if (!is.finite(estimate)) stop("trme TR-AIPW risk difference was not finite")
  if (length(boot) < min_boot) {
    stop("Too few valid TR-AIPW bootstrap estimates: ",
         length(boot), " of ", B, ".")
  }

  ## Keep bootstrap SD as the reported SE diagnostic, but use percentile
  ## bootstrap limits for the confidence interval.
  se <- stats::sd(boot)
  ci <- stats::quantile(boot, probs=c(0.025,0.975),
                        na.rm=TRUE, names=FALSE, type=1)
  ci_lower <- as.numeric(ci[1])
  ci_upper <- as.numeric(ci[2])

  ps <- bound_probability(as.numeric(fit$fit_ps))
  mw <- as.numeric(fit$miss_weights)
  mw <- mw[is.finite(mw)]
  diag <- data.frame(
    n_boot_valid=length(boot),
    n_boot_invalid=B-length(boot),
    boot_median=stats::median(boot),
    boot_q025=unname(stats::quantile(boot,0.025,na.rm=TRUE,type=1)),
    boot_q975=unname(stats::quantile(boot,0.975,na.rm=TRUE,type=1)),
    boot_min=min(boot),
    boot_max=max(boot),
    q01_treated_ps=unname(stats::quantile(ps, 0.01, na.rm=TRUE)),
    q99_treated_inv_weight=unname(stats::quantile(1/ps, 0.99, na.rm=TRUE)),
    q99_control_inv_weight=unname(stats::quantile(1/(1-ps), 0.99, na.rm=TRUE)),
    q99_missingness_weight=if(length(mw)) unname(stats::quantile(mw, 0.99, na.rm=TRUE)) else NA_real_,
    max_missingness_weight=if(length(mw)) max(mw, na.rm=TRUE) else NA_real_
  )
  cbind(data.frame(estimate=estimate,se=se,ci_lower=ci_lower,ci_upper=ci_upper),diag)
}

fit_one_method <- function(method,dat,seed=1L,mice_m=10L,mice_maxit=10L,
                           nuisance_spec="main",trim=DEFAULT_TRIM,trme_boot_B=100L,trme_shrink_rate=0.99) {
  dd<-make_analysis_data(dat)
  # Generated datasets do not store the mechanism column, so fit_engine passes it
  # as an attribute before dispatch.
  mechanism<-attr(dat,"missing_mechanism") %||% "MAR"
  switch(method,
    full_aipw=safe_result(run_aipw_result(dd$full,nuisance_spec,trim)),
    cc_aipw=safe_result(run_aipw_result(dd$complete_cases,nuisance_spec,trim)),
    single_imp_aipw=safe_result(run_single_imp_result(dd$missing,nuisance_spec,trim)),
    mice_aipw=safe_result(run_mice_result(dd$missing,mice_m,mice_maxit,seed,nuisance_spec,trim)),
    drcmd=safe_result(run_drcmd_result(dd$missing,nuisance_spec,trim)),
    tr_aipw=safe_result(
      run_trme_rd_result(dd$missing,nuisance_spec,mechanism,
                         trme_shrink_rate,trme_boot_B,seed),
      empty=empty_trme_result()
    ),
    stop("Unknown method: ",method))
}

`%||%` <- function(x,y) if(is.null(x)||length(x)==0) y else x
