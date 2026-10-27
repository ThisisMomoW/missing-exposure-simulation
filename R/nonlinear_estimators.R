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

pool_rubin <- function(q, u) {
  keep <- is.finite(q) & is.finite(u)
  q <- q[keep]
  u <- u[keep]
  m <- length(q)
  if (m < 2L) stop("Fewer than two valid imputations.")
  qbar <- mean(q)
  ubar <- mean(u)
  b <- stats::var(q)
  total <- ubar + (1 + 1 / m) * b
  se <- sqrt(total)
  data.frame(
    estimate = qbar,
    se = se,
    ci_lower = qbar - 1.96 * se,
    ci_upper = qbar + 1.96 * se
  )
}

extract_aipw_result <- function(fit, trim) {
  invisible(capture.output(fit$summary(g.bound = trim)))
  
  res <- fit$result
  
  if (!is.matrix(res) && !is.data.frame(res)) {
    stop("AIPW fit$result is not a matrix or data frame.")
  }
  
  # Binary outcome only: do not silently accept continuous-outcome data.
  hit <- which(
    trimws(rownames(res)) == "Risk Difference"
  )
  
  if (length(hit) != 1L) {
    stop(
      "Cannot identify the ATE row. Available rows: ",
      paste(rownames(res), collapse = ", ")
    )
  }
  
  required <- c("Estimate", "SE", "95% LCL", "95% UCL")
  
  if (!all(required %in% colnames(res))) {
    stop(
      "Unexpected AIPW columns: ",
      paste(colnames(res), collapse = ", ")
    )
  }
  
  values <- as.numeric(unlist(
    res[hit, required, drop = FALSE],
    use.names = FALSE
  ))
  
  if (any(!is.finite(values)) || values[2] <= 0) {
    stop("AIPW returned nonfinite values or a nonpositive SE.")
  }
  
  data.frame(
    estimate = values[1],
    se = values[2],
    ci_lower = values[3],
    ci_upper = values[4]
  )
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
  # RF exposure imputation in both outcome-nuisance comparison arms.
  Xcc <- cc[, c("X1", "X2", "X3", "Y"), drop = FALSE]
  Xall <- df_missing[, c("X1", "X2", "X3", "Y"), drop = FALSE]
  imp_fit <- SuperLearner::SuperLearner(
    Y = cc$A, X = Xcc, family = binomial(), SL.library = "SL.ranger"
  )
  ph <- bound_probability(as.numeric(
    predict(imp_fit, newdata = Xall, onlySL = TRUE)$pred
  ))
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
  # RF exposure imputation in both outcome-nuisance comparison arms.
  meth["A"] <- "rf"

  pm <- mice::make.predictorMatrix(d)
  pm[,] <- 0
  pm["A", c("Y", "X1", "X2", "X3")] <- 1

  imp <- mice::mice(
    d, m = m, maxit = maxit, method = meth, predictorMatrix = pm,
    seed = seed, printFlag = FALSE, rfPackage = "ranger", ntree = ntree
  )

  rr <- lapply(seq_len(m), function(k) {
    z <- mice::complete(imp, k)
    z$A <- as.integer(as.character(z$A))
    run_aipw(z, nuisance_spec, trim)
  })

  q <- vapply(rr, `[[`, numeric(1), "estimate")
  se <- vapply(rr, `[[`, numeric(1), "se")
  pool_rubin(q, se^2)
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
