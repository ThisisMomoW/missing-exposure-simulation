# AIPW, imputation and package drcmd estimators for the nonlinear binary DGP.

bound_probability <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

empty_result <- function(status = "failed", error = NA_character_) {
  data.frame(
    estimate = NA_real_, se = NA_real_,
    ci_lower = NA_real_, ci_upper = NA_real_,
    status = status, error = error,
    stringsAsFactors = FALSE
  )
}

safe_result <- function(expr) {
  tryCatch({
    ans <- expr
    ans$status <- "ok"
    ans$error <- NA_character_
    ans
  }, error = function(e) empty_result(error = conditionMessage(e)))
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

run_aipw <- function(df, nuisance_spec, trim) {
  q_library <- if (identical(nuisance_spec, "rf")) "SL.ranger" else "SL.glm"

  invisible(capture.output(
    fit <- AIPW::AIPW$new(
      Y = df$Y,
      A = df$A,
      W = df[, c("X1", "X2", "X3"), drop = FALSE],
      Q.SL.library = q_library,
      g.SL.library = "SL.glm"
    )$fit()
  ))

  extract_aipw_result(fit, trim)
}

make_analysis_data <- function(dat) {
  full <- transform(dat, A = as.integer(A))
  missing <- transform(dat, A = as.integer(A_obs))
  cc <- missing[missing$R == 1L, , drop = FALSE]
  list(full = full, missing = missing, complete_cases = cc)
}

run_single_imp <- function(df_missing, nuisance_spec, trim) {
  cc <- df_missing[!is.na(df_missing$A), , drop = FALSE]
  if (nuisance_spec == "glm") {
    imp_fit <- glm(A ~ X1 + X2 + X3 + Y, data=cc, family=binomial())
    ph <- bound_probability(predict(imp_fit, newdata=df_missing, type="response"))
  } else {
    Xcc <- cc[, c("X1", "X2", "X3", "Y"), drop = FALSE]
    Xall <- df_missing[, c("X1", "X2", "X3", "Y"), drop = FALSE]
    imp_fit <- SuperLearner::SuperLearner(
      Y = cc$A, X = Xcc, family = binomial(), SL.library = "SL.ranger"
    )
    ph <- bound_probability(as.numeric(
      predict(imp_fit, newdata = Xall, onlySL = TRUE)$pred
    ))
  }
  d <- df_missing
  mis <- is.na(d$A)
  d$A[mis] <- rbinom(sum(mis), 1, ph[mis])
  run_aipw(d, nuisance_spec, trim)
}

run_mice <- function(df_missing, nuisance_spec, m, maxit, seed, trim, ntree=500L) {
  d <- df_missing[, c("A", "Y", "X1", "X2", "X3")]
  d$A <- factor(d$A, levels = c(0, 1))

  meth <- mice::make.method(d)
  meth[] <- ""
  meth["A"] <- if (nuisance_spec == "rf") "rf" else "logreg"

  pm <- mice::make.predictorMatrix(d)
  pm[,] <- 0
  pm["A", c("Y", "X1", "X2", "X3")] <- 1

  imp <- if (nuisance_spec == "rf") {
    mice::mice(d, m=m, maxit=maxit, method=meth, predictorMatrix=pm,
      seed=seed, printFlag=FALSE, rfPackage="ranger", ntree=ntree)
  } else {
    mice::mice(d, m=m, maxit=maxit, method=meth, predictorMatrix=pm,
      seed=seed, printFlag=FALSE)
  }

  rr <- lapply(seq_len(m), function(k) {
    z <- mice::complete(imp, k)
    z$A <- as.integer(as.character(z$A))
    run_aipw(z, nuisance_spec, trim)
  })

  q <- vapply(rr, `[[`, numeric(1), "estimate")
  se <- vapply(rr, `[[`, numeric(1), "se")
  retain_mice_inputs(pool_rubin(q, se^2), q, se^2)
}

run_drcmd <- function(df_missing, nuisance_spec, trim) {
  # Compare m/po both GLM versus m/po both RF; g and r remain GLM.
  m_library <- if (identical(nuisance_spec, "rf")) "SL.ranger" else "SL.glm"

  fit <- drcmd::drcmd(
    Y = df_missing$Y,
    A = df_missing$A,
    X = df_missing[, c("X1", "X2", "X3"), drop = FALSE],
    m_learners = m_library,
    g_learners = "SL.glm",
    r_learners = "SL.glm",
    po_learners = m_library,
    cutoff = trim,
    k = 10L
  )

  est_obj <- fit$results$estimates
  se_obj <- fit$results$ses
  est <- if (is.data.frame(est_obj)) est_obj$psi_hat_ate[1] else
    unlist(est_obj)["psi_hat_ate"]
  se <- if (is.data.frame(se_obj)) se_obj$psi_hat_ate[1] else
    unlist(se_obj)["psi_hat_ate"]

  if (!is.finite(est) || !is.finite(se)) {
    stop("drcmd ATE estimate or SE was not found.")
  }

  data.frame(
    estimate = as.numeric(est),
    se = as.numeric(se),
    ci_lower = as.numeric(est) - 1.96 * as.numeric(se),
    ci_upper = as.numeric(est) + 1.96 * as.numeric(se)
  )
}

fit_one <- function(method, nuisance_spec, dat, seed, mice_rf_trees=500L) {
  dd <- make_analysis_data(dat)

  switch(
    method,
    full_aipw = safe_result(run_aipw(dd$full, nuisance_spec, trim_cutoff)),
    cc_aipw = safe_result(run_aipw(dd$complete_cases, nuisance_spec,
                                   trim_cutoff)),
    single_imp_aipw = safe_result(run_single_imp(dd$missing, nuisance_spec,
                                                 trim_cutoff)),
    mice_aipw = safe_result(run_mice(dd$missing, nuisance_spec, mice_m,
                                     mice_maxit, seed, trim_cutoff, ntree=mice_rf_trees)),
    drcmd = safe_result(run_drcmd(dd$missing, nuisance_spec, trim_cutoff)),
    stop("Unknown method: ", method)
  )
}
